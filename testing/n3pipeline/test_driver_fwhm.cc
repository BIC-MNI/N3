/* Cycle 15, item 1: -fwhm and -distance are distinct at this level.
 *
 * At nu_estimate/nu_correct level -fwhm is the sharpening width
 * (nu_estimate.in:205 maps it to 'sharpen', passed on as -sharpen <width> 0.01
 * at :498); only the inner nu_estimate_np_and_em.in:1163 calls the sharpening
 * width -fwhm.  -distance is the knot spacing (:164-165).  A re-aliasing of
 * -fwhm to -distance must fail an assertion here, not merely change the
 * output -- so the test drives the compiled binary and compares its outputs:
 *
 *   -fwhm 0.15 reproduces the default run byte-for-byte;
 *   -fwhm 0.3  moves the correction (through the sharpening channel);
 *   -distance 100 moves it through the knot channel.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>

#ifndef N3_DRIVER_BIN
#error "N3_DRIVER_BIN must be the path of the built nu_correct_cxx binary"
#endif

static std::string outdir()
{
  const char *tmp = getenv("TMPDIR");
  return std::string(tmp ? tmp : "/tmp");
}

/* Run the driver on chunk.mnc with the given extra options, writing to a
 * per-pid file in TMPDIR.  Returns the output path, or "" on failure. */
static std::string run(const std::string &opts, const char *tag)
{
  char path[256], cmd[1024];
  snprintf(path, sizeof(path), "%s/n3cxx_fwhm_%d_%s.mnc",
           outdir().c_str(), (int) getpid(), tag);
  /* -V1.0 -nolegacy_rounding pins this test's numerics to the protocol they
   * were measured against (fwhm 0.15, triangular window, legacy_rounding
   * off), independent of whichever protocol the driver's own implicit
   * default currently selects (2026-08-06: that default moved to -V1.1). */
  snprintf(cmd, sizeof(cmd), "\"%s\" -V1.0 -nolegacy_rounding -shrink 1 "
           "-iterations 1 -stop 0.0 -distance 200 "
           "-mask \"%s/chunk_mask.mnc.gz\" %s "
           "\"%s/chunk.mnc.gz\" \"%s\" -clobber",
           N3_DRIVER_BIN, N3_DATA_DIR, opts.c_str(),
           N3_DATA_DIR, path);
  if(system(cmd) != 0)
    {
      printf("FAIL: driver exited non-zero for -%s run\n", tag);
      return "";
    }
  return path;
}

static void cleanup(const char *path)
{
  /* The driver's .imp for a correct output <out>.mnc is <out>.imp (final
   * extension replaced, MNI::PathUtilities::replace_ext) -- not <out>.mnc.imp
   * -- so the volume, that .imp, and any <out>* .log all go here. */
  std::string ip(path);
  size_t dot = ip.find_last_of('.');
  ip = (dot == std::string::npos ? ip : ip.substr(0, dot)) + ".imp";
  char cmd[512];
  snprintf(cmd, sizeof(cmd), "rm -f \"%s\" \"%s\" \"%s\"*.log",
           path, ip.c_str(), path);
  system(cmd);
}

int main()
{
  std::string base = run("",                         "base");
  std::string same = run("-fwhm 0.15",               "same");
  std::string fw   = run("-fwhm 0.3",                "fw");
  std::string dist = run("-distance 100 -fwhm 0.15", "dist");

  bool ok = !base.empty() && !same.empty() && !fw.empty() && !dist.empty();
  if(!ok) { n3check::failures()++; return n3check::report("driver_fwhm"); }

  VIO_Volume w_base = n3::load(base);
  VIO_Volume w_same = n3::load(same);
  VIO_Volume w_fw   = n3::load(fw);
  VIO_Volume w_dist = n3::load(dist);
  int n = n3::voxel_count(w_base);
  const double *vb = n3::values(w_base);
  const double *vs = n3::values(w_same);
  const double *vf = n3::values(w_fw);
  const double *vd = n3::values(w_dist);

  printf("  (the default sharpens at 0.15 with 200 mm knots)\n");
  double d_same = n3check::rel_rms(vb, vs, n);
  double d_fw   = n3check::rel_rms(vb, vf, n);
  double d_dist = n3check::rel_rms(vb, vd, n);
  double d_diff = n3check::rel_rms(vf, vd, n);
  printf("  rel RMS: default vs -fwhm 0.15 = %.3e,  vs -fwhm 0.3 = %.3e,\n"
         "        vs -distance 100 = %.3e,  fwhm0.3 vs dist100 = %.3e\n",
         d_same, d_fw, d_dist, d_diff);

  CHECK_TRUE("default and -fwhm 0.15 agree byte-for-byte", d_same == 0.0);

  /* A "the correction moved" check must clear the output's own quantisation:
   * chunk.mnc is 12-bit, so any corrected value can differ by one level of the
   * (range/valid_steps) quantum, and rel RMS of a single-level difference is
   * of that order.  Requiring the movement to exceed one whole quantum --
   * derived from the data, not a round figure -- proves the channel changed
   * the correction more than a least-significant-bit wobble could.  Measured
   * (2026-08-06): quantum 5.4e-04 of the in-mask mean; -fwhm 0.3 moves it
   * 18.7x that, -distance 100 moves it 2.23x (see print above). */
  VIO_Volume mask = n3::load(std::string(N3_DATA_DIR) + "/chunk_mask.mnc.gz");
  n3::Stats bs = n3::masked_stats(w_base, mask);
  double quantum = (bs.maximum - bs.minimum)
    / n3fixture::valid_steps("chunk_valid_range.txt") / bs.mean;
  printf("  (one output quantum is %.3e of the in-mask mean)\n", quantum);
  CHECK_TRUE("-fwhm 0.3  moves the correction", d_fw > quantum);
  CHECK_TRUE("-distance 100 moves the correction", d_dist > quantum);
  CHECK_TRUE("and the two channels move it differently", d_diff > quantum);
  delete_volume(mask);

  delete_volume(w_dist);
  delete_volume(w_fw);
  delete_volume(w_same);
  delete_volume(w_base);
  cleanup(base.c_str());
  cleanup(same.c_str());
  cleanup(fw.c_str());
  cleanup(dist.c_str());
  return n3check::report("driver_fwhm");
}
