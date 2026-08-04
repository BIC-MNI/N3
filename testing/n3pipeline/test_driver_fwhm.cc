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
  snprintf(cmd, sizeof(cmd), "\"%s\" -shrink 1 -iterations 1 -stop 0.0 "
           "-distance 200 -mask \"%s/chunk_mask.mnc.gz\" %s "
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
  char cmd[512];
  snprintf(cmd, sizeof(cmd), "rm -f \"%s\" \"%s.imp\" \"%s\"*.log",
           path, path, path);
  system(cmd);
}

/* rel RMS over all voxels (not strided): determinism means identical runs give
 * exactly zero, which is the strongest possible statement about an alias. */
static double all_rms(const double *a, const double *b, int n)
{
  double sd = 0.0, sb = 0.0;
  for(int i = 0; i < n; i++) { double d = a[i] - b[i]; sd += d*d; sb += b[i]*b[i]; }
  return (sb == 0.0) ? (sd == 0.0 ? 0.0 : 1.0/0.0) : sqrt(sd/sb);
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
  CHECK_TRUE("default and -fwhm 0.15 agree byte-for-byte",
             all_rms(vb, vs, n) == 0.0);
  CHECK_TRUE("-fwhm 0.3  moves the correction", all_rms(vb, vf, n) > 1e-4);
  CHECK_TRUE("-distance 100 moves the correction", all_rms(vb, vd, n) > 1e-4);
  CHECK_TRUE("and the two channels move it differently",
             all_rms(vf, vd, n) > 1e-4);

  delete_volume(w_dist);
  delete_volume(w_fw);
  delete_volume(w_same);
  delete_volume(w_base);
  cleanup(base.c_str());
  cleanup(fw.c_str());
  cleanup(dist.c_str());
  return n3check::report("driver_fwhm");
}
