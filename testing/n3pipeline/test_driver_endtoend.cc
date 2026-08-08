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

/* Run the driver on chunk.mnc at the cycle-14 protocol with the given extra
 * options, writing to a per-pid file in TMPDIR.  Returns the output path, or
 * "" on failure. */
static std::string run(const std::string &opts, const char *tag)
{
  std::string data = N3_DATA_DIR;
  std::string out = n3fixture::temp_path(std::string("e2e_") + tag + ".mnc");
  std::string all = opts + " -shrink 1 -iterations 1 -stop 0.0 -distance 200"
    " -mask \"" + data + "/chunk_mask.mnc\"";
  if(!n3fixture::run_driver(all, data + "/chunk.mnc", out))
    {
      printf("FAIL: driver exited non-zero for the %s run\n", tag);
      return "";
    }
  return out;
}

/* The driver's output on the same stride the oracle was recorded at. */
static std::vector<double> corrected(const std::string &path)
{
  VIO_Volume v = n3::load(path);
  std::vector<double> out = n3fixture::strided(n3::values(v),
                                               n3::voxel_count(v));
  delete_volume(v);
  return out;
}

static void compare(const char *what, const std::vector<double> &ours,
                    double bound)
{
  std::vector<double> oracle = n3fixture::read_f64("nu_correct_shrink1.f64");
  n3fixture::must(ours.size() == oracle.size(),
                  "nu_correct_shrink1.f64: size mismatch");
  CHECK_RMS(what, &ours[0], &oracle[0], (int) ours.size(), bound);
}

int main()
{
  std::string data = std::string(N3_DATA_DIR);
  VIO_Volume chunk = n3::load(data + "/chunk.mnc");
  VIO_Volume mask  = n3::load(data + "/chunk_mask.mnc");

  double bound = round_trip_bound(chunk, mask);
  printf("  (a 12-bit round trip of the log volume is worth %.3e)\n", bound);

  /* -V1.0 selects the Perl's default protocol (fwhm 0.15, linear window) and,
   * by default, legacy_rounding, which reproduces its %lf rounding; the extra
   * flag turns that off so the two sides of the check bracket the driver. */
  std::string legacy_on  = run("-V1.0", "on");
  std::string legacy_off = run("-V1.0 -nolegacy_rounding", "off");

  bool ok = !legacy_on.empty() && !legacy_off.empty();
  if(!ok) { n3check::failures()++; return n3check::report("driver_endtoend"); }

  std::vector<double> on = corrected(legacy_on), off = corrected(legacy_off);

  compare("end to end vs nu_correct, -legacy_rounding on", on, bound);
  compare("end to end vs nu_correct, -legacy_rounding off", off, bound);

  /* Both comparisons above go against the same oracle under the same bound,
   * and both land at 1.84e-04.  That alone would pass unchanged if
   * -legacy_rounding did nothing, which is not hypothetical: the option had
   * no CLI counterpart until 9154eac, and EstimateOptions::legacy_rounding
   * was left at its false default by the driver, so "every driver run to date
   * has been the equivalent of -nolegacy_rounding" and no test noticed.
   *
   * So the two runs are also compared against each other, with no bound
   * fitted to the result.  Below: they must differ at all -- a property, and
   * the one an inert flag fails.  Above: the difference must stay under the
   * same round-trip quantum the cycle is already bounded by, since a rounding
   * of intermediates to six decimals cannot legitimately move the output by
   * more than the 12-bit file quantisation those intermediates pass through.
   * For scale, PLAN §4 records the lookup-position component of this rounding
   * at 5.114e-07 on chunk.mnc, measured in cycle 6; end to end at one
   * iteration the whole flag is worth 8.740e-07 here, the remainder being the
   * counts and the histogram domain, which are rounded too. */
  double rounding = n3check::rel_rms(&on[0], &off[0], (int) on.size());
  printf("  (-legacy_rounding moves the output by %.3e; cycle 6 measured the\n"
         "   lookup-position component at 5.114e-07)\n", rounding);
  CHECK_TRUE("-legacy_rounding is not inert: the two runs differ",
             rounding > 0.0);
  CHECK_RMS("and it moves less than one round trip of the intermediates",
            &on[0], &off[0], (int) on.size(), bound);

  n3fixture::cleanup(legacy_on);
  n3fixture::cleanup(legacy_off);
  delete_volume(mask);
  delete_volume(chunk);
  return n3check::report("driver_endtoend");
}
