/* Cycle 14: the one end-to-end comparison whose bound is justified in advance.
 *
 * Every other end-to-end number is reported, not asserted, because the stopping
 * rule quantises what follows it and the installed pipeline's file round trips
 * round the intermediates.  This cycle cuts both effects away at once: one
 * iteration at -stop 0.0 runs a fixed count, and the two MINC round trips are
 * worth half a level of the log range over the file's own quantum.
 *
 * The pipeline is the Perl top level, `nu_correct`: estimate and evaluate in
 * one pass, so the oracle (nu_correct_shrink1.f64) is the installed nu_correct
 * run on the same command line:
 *
 *   nu_correct -shrink 1 -iterations 1 -stop 0.0 -distance 200 -mask <mask>
 *              <chunk> <out> -clobber
 *
 * (regenerate_reference.sh, cycle 14).  The driver is run with and without
 * -legacy_rounding: -V1.0 reproduces the Perl's default protocol (fwhm 0.15,
 * linear interpolation) and turns legacy_rounding on; -nolegacy_rounding turns
 * that off.  Both must clear the quantum bound, which is a full-order margin
 * above the %lf rounding the two sides' round trips otherwise agree to -- that
 * is the point of cycling this end to end rather than only per block.
 *
 * The bound is measured, not assumed: 0.5 (log max - log min) / valid_steps
 * over the chunk's own masked intensities and its recorded valid_range
 * (chunk_valid_range.txt, 0..4095 -> 4095), which at the 12-bit quantum is
 * 2.683e-4.  A fixed 4095 or a fixed 16-bit figure would be wrong for any
 * other file.
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

/* The round-trip bound, in the masked-log-intensity terms of PLAN §4.  It must
 * be computed before any comparison runs and without looking at the driver's
 * output, so the test's own verdict does not depend on what is being tested. */
static double round_trip_bound(VIO_Volume volume, VIO_Volume mask)
{
  n3::Stats inside = n3::masked_stats(volume, mask);
  n3fixture::must(inside.minimum > 0.0,
                  "masked minimum is not positive (log(min) would be -inf)");
  return 0.5 * (log(inside.maximum) - log(inside.minimum))
       / n3fixture::valid_steps("chunk_valid_range.txt");
}

static std::string outdir()
{
  const char *tmp = getenv("TMPDIR");
  return std::string(tmp ? tmp : "/tmp");
}

/* Run the driver on chunk.mnc at the cycle-14 protocol with the given extra
 * options, writing to a per-pid file in TMPDIR.  Returns the output path, or
 * "" on failure. */
static std::string run(const std::string &opts, const char *tag)
{
  char cmd[1024];
  snprintf(cmd, sizeof(cmd),
           "\"%s\" %s -shrink 1 -iterations 1 -stop 0.0 -distance 200 "
           "-mask \"%s/chunk_mask.mnc.gz\" "
           "\"%s/chunk.mnc.gz\" \"%s/n3cxx_e2e_%d_%s.mnc\" -clobber",
           N3_DRIVER_BIN, opts.c_str(), N3_DATA_DIR, N3_DATA_DIR,
           outdir().c_str(), (int) getpid(), tag);
  if(system(cmd) != 0)
    {
      printf("FAIL: driver exited non-zero for the %s run\n", tag);
      return "";
    }
  std::string out = std::string(outdir()) + "/n3cxx_e2e_" +
    std::to_string((int) getpid()) + "_" + tag + ".mnc";
  return out;
}

static void cleanup(const char *path)
{
  std::string cmd = std::string("rm -f \"") + path + "\" \""
    + path + ".imp\" \"" + path + "*.log\"";
  (void) system(cmd.c_str());   /* cleanup is best-effort */
}

static void compare(const char *what, const char *path, double bound)
{
  VIO_Volume mine = n3::load(path);
  int n = n3::voxel_count(mine);
  std::vector<double> ours = n3fixture::strided(n3::values(mine), n);
  std::vector<double> oracle = n3fixture::read_f64("nu_correct_shrink1.f64");
  n3fixture::must(ours.size() == oracle.size(),
                  "nu_correct_shrink1.f64: size mismatch");
  CHECK_RMS(what, &ours[0], &oracle[0], (int) ours.size(), bound);
  delete_volume(mine);
}

int main()
{
  std::string data = std::string(N3_DATA_DIR);
  VIO_Volume chunk = n3::load(data + "/chunk.mnc.gz");
  VIO_Volume mask  = n3::load(data + "/chunk_mask.mnc.gz");

  double bound = round_trip_bound(chunk, mask);
  printf("  (a 12-bit round trip of the log volume is worth %.3e)\n", bound);

  /* -V1.0 selects the Perl's default protocol (fwhm 0.15, linear window) and,
   * by default, legacy_rounding, which reproduces its %lf rounding; the extra
   * flag turns that off so the two sides of the check bracket the driver. */
  std::string legacy_on  = run("-V1.0", "on");
  std::string legacy_off = run("-V1.0 -nolegacy_rounding", "off");

  bool ok = !legacy_on.empty() && !legacy_off.empty();
  if(!ok) { n3check::failures()++; return n3check::report("driver_endtoend"); }

  compare("end to end vs nu_correct, -legacy_rounding on",
          legacy_on.c_str(), bound);
  compare("end to end vs nu_correct, -legacy_rounding off",
          legacy_off.c_str(), bound);

  cleanup(legacy_on.c_str());
  cleanup(legacy_off.c_str());
  delete_volume(mask);
  delete_volume(chunk);
  return n3check::report("driver_endtoend");
}
