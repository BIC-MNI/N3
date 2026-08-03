/* Cycle 11: one iteration of the estimation loop.
 *
 * The first integration cycle, and it is per-iteration rather than end to end
 * because the driver will hand over its own intermediates: -save_fields writes
 * exp(log - sharpened) before the smoothing and exp(residue) after it.  Without
 * that this cycle would have been one opaque comparison.
 *
 * Run at -shrink 1 -iterations 1 -stop 0.0, which is where the fewest MINC
 * round trips stand between the two implementations: the driver's saved
 * volumes are 16-bit like everything else it writes, and that quantisation is
 * the bound here.  -stop 0.0 never fires, so both sides run exactly one
 * iteration -- the stopping rule quantises everything downstream and has to be
 * taken out of the comparison before anything else can be read.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"
#include "../../src/N3Pipeline/NuEstimate.h"

#include <cmath>
#include <string>
#include <vector>

/* Every intermediate the driver passes between programs is a MINC file, and
 * mincmath takes the header -- and so the valid_range -- from its first input.
 * For these the lineage is the 12-bit test data, so the log-intensity volumes
 * the sharpening reads and writes carry 4096 levels across their own range.
 *
 * Half a quantum of that is the bound: the field is the exponential of a
 * difference of two such volumes, so a relative error of half a quantum of the
 * log range is what the round trips cost.  On chunk.mnc the masked intensities
 * span 100030 to 900274, a log range of 2.197, and the bound is 2.7e-04.
 *
 * Not 16 bits: the driver does write its masks -short -signed, but these
 * volumes inherit the input's own valid_range. */
static double round_trip_bound(VIO_Volume volume, VIO_Volume mask)
{
  n3::Stats inside = n3::masked_stats(volume, mask);
  return 0.5 * (log(inside.maximum) - log(inside.minimum)) / 4095.0;
}

static void compare(const char *what, VIO_Volume mine, const char *fixture,
                    double bound)
{
  std::vector<double> oracle = n3fixture::read_f64(fixture);
  std::vector<double> ours =
    n3fixture::strided(n3::values(mine), n3::voxel_count(mine));
  n3fixture::must(ours.size() == oracle.size(),
                  std::string(fixture) + ": size mismatch");
  CHECK_RMS(what, &ours[0], &oracle[0], (int) ours.size(), bound);
}

int main()
{
  std::string data = std::string(N3_DATA_DIR);

  VIO_Volume chunk = n3::load(data + "/chunk.mnc.gz");
  VIO_Volume mask  = n3::load(data + "/chunk_mask.mnc.gz");

  n3::EstimateOptions options;
  options.shrink = 1.0;
  options.distance = 200.0;
  options.lambda = 1e-7;
  options.subsample = 1;
  options.window = true;
  options.iterations.assign(1, 1);
  options.stop.assign(1, 0.0);

  n3::EstimateTrace trace;
  trace.iteration = 0;

  int iterations = 0;
  double change = 0.0;
  n3::Field *field = n3::nu_estimate(chunk, mask, options, &iterations,
                                     &change, &trace);

  /* The iteration count first.  A difference here makes every number below
   * meaningless, so it is checked before them and not after. */
  CHECK_TRUE("exactly one iteration ran", iterations == 1);

  /* field_CV, which the driver prints and which decides when the real
   * protocol stops.  Its six digits are all the log carries. */
  {
    double oracle = n3fixture::read_scalar("estimate_change.txt");
    CHECK_NEAR("the change in the field is the driver's", change, oracle,
               0.5e-8 + 1e-3 * fabs(oracle));
    printf("  (change %.8f against the driver's %.8f)\n", change, oracle);
  }

  n3fixture::must(trace.before_smoothing != NULL && trace.field != NULL,
                  "the trace was not filled in");

  double bound = round_trip_bound(chunk, mask);
  printf("  (a 12-bit round trip of the log volume is worth %.3e)\n", bound);

  compare("the residual field before smoothing is the driver's",
          trace.before_smoothing, "estimate_est0.f64", bound);
  compare("the fitted field after smoothing is the driver's",
          trace.field, "estimate_field0.f64", bound);

  /* What -legacy_rounding is for: it reproduces the six decimals the histogram
   * and the lookup table pass between programs with, leaving the volume
   * quantisation as the only divergence.  Here it changes nothing measurable,
   * which places the whole of that divergence below the round trips -- worth
   * knowing before spending any effort on it. */
  {
    n3::EstimateOptions rounded = options;
    rounded.legacy_rounding = true;
    n3::EstimateTrace rounded_trace;
    rounded_trace.iteration = 0;
    int its = 0;
    double ch = 0.0;
    n3::Field *other = n3::nu_estimate(chunk, mask, rounded, &its, &ch,
                                       &rounded_trace);

    int n = n3::voxel_count(trace.field);
    double drift = n3check::rel_rms(n3::values(rounded_trace.field),
                                    n3::values(trace.field), n);
    CHECK_TRUE("the six-decimal rounding is smaller than the round trips",
               drift < bound);
    printf("  (-legacy_rounding moves the field by %.3e)\n", drift);

    delete other;
    if(rounded_trace.before_smoothing) delete_volume(rounded_trace.before_smoothing);
    if(rounded_trace.field) delete_volume(rounded_trace.field);
  }

  /* Properties of the result, which hold whatever the driver produced.
   *
   * The field is a multiplicative correction, so it must be positive
   * everywhere the mask covers; a field that reached zero would divide the
   * output to infinity. */
  {
    VIO_Volume evaluated = n3::like(chunk);
    n3::evaluate_field(field, evaluated, mask);
    double *f = n3::values(evaluated), *m = n3::values(mask);
    int n = n3::voxel_count(evaluated);
    int nonpositive = 0;
    double lo = 0.0, hi = 0.0;
    bool first = true;
    for(int i = 0; i < n; i++)
      if(m[i] > 0.5)
        {
          if(f[i] <= 0.0) nonpositive++;
          if(first) { lo = hi = f[i]; first = false; }
          if(f[i] < lo) lo = f[i];
          if(f[i] > hi) hi = f[i];
        }
    CHECK_NEAR("the field is positive throughout the mask", nonpositive, 0.0, 0.0);
    printf("  (field spans %.4f .. %.4f in the mask)\n", lo, hi);

    /* And it is smooth: at 200 mm knot spacing over a 182 mm volume the
     * field cannot swing by much, so a field that had absorbed the anatomy
     * would show up here. */
    CHECK_TRUE("and varies by less than a factor of two", hi / lo < 2.0);

    delete_volume(evaluated);
  }

  delete field;
  if(trace.before_smoothing) delete_volume(trace.before_smoothing);
  if(trace.field) delete_volume(trace.field);
  delete_volume(mask);
  delete_volume(chunk);
  return n3check::report("estimate");
}
