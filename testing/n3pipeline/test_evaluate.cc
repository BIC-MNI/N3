/* Cycle 12: nu_evaluate.
 *
 * Driven from the driver's own .imp, so that the estimation's differences do
 * not enter: both sides start from the same field and what is compared is the
 * evaluation, the extension, the floor and the division.
 *
 * The bound is the same 12-bit round trip as cycle 11 -- nu_evaluate's field
 * volume goes to disk between evaluate_field, correct_field and mincmath, and
 * inherits the input's valid_range.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"
#include "../../src/N3Pipeline/FitField.h"
#include "../../src/N3Pipeline/MincTools.h"
#include "../../src/N3Pipeline/NuEvaluate.h"

#include <cmath>
#include <string>
#include <vector>

static double round_trip_bound(VIO_Volume volume, VIO_Volume mask)
{
  n3::Stats inside = n3::masked_stats(volume, mask);
  /* log(0) is -inf and would poison the bound; a masked minimum at zero is a
   * data or mask problem the test should refuse rather than tolerate. */
  n3fixture::must(inside.minimum > 0.0,
                  "masked minimum is not positive (log(min) would be -inf)");
  return 0.5 * (log(inside.maximum) - log(inside.minimum))
       / n3fixture::valid_steps("chunk_valid_range.txt");
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
  std::string reference = std::string(N3_REFERENCE_DIR);

  VIO_Volume chunk = n3::load(data + "/chunk.mnc");
  VIO_Volume mask  = n3::load(data + "/chunk_mask.mnc");
  int n = n3::voxel_count(chunk);

  double bound = round_trip_bound(chunk, mask);
  printf("  (a 12-bit round trip of the log volume is worth %.3e)\n", bound);

  /* The field from the driver's .imp, evaluated on the input's grid.  It is
   * zero outside the mask before the extension, which is the entire reason
   * nu_evaluate runs correct_field. */
  VIO_Volume field = n3::like(chunk);
  n3::evaluate_saved_field(reference + "/estimate.imp", field, mask);
  {
    double *f = n3::values(field), *m = n3::values(mask);
    int zero_outside = 0, outside = 0;
    for(int i = 0; i < n; i++)
      if(m[i] <= 0.5) { outside++; if(f[i] == 0.0) zero_outside++; }
    CHECK_TRUE("the field from an .imp is zero outside the mask",
               zero_outside == outside && outside > 0);
  }
  delete_volume(field);

  /* The corrected volume, from the same field.  Both sides divide by the
   * field from estimate.imp, so only the evaluation stages are compared. */
  n3::Field *loaded = n3::load_field(reference + "/estimate.imp", chunk);

  n3::EvaluateOptions options;
  VIO_Volume field_out = NULL;
  VIO_Volume corrected = n3::nu_evaluate(chunk, mask, loaded, options,
                                         &field_out);

  compare("the corrected output is the driver's", corrected,
          "evaluate_corrected.f64", bound);

  /* The floor is clamped only where the minimum falls below it; on this field
   * nothing does, so the output is exactly input / field and the field_out
   * volume the driver would have gone on to divide is the same one, extended
   * and floored.  Re-derive the division and check it matches the recorded
   * output to the same bound -- what the test above compared could in
   * principle have come out right for a wrong field. */
  {
    std::vector<double> oracle = n3fixture::read_f64("evaluate_corrected.f64");
    std::vector<double> f = n3fixture::strided(n3::values(field_out), n);
    std::vector<double> v = n3fixture::strided(n3::values(chunk), n);
    std::vector<double> recomputed(f.size());
    for(size_t i = 0; i < f.size(); i++)
      recomputed[i] = (f[i] == 0.0) ? 0.0 : v[i] / f[i];
    n3fixture::must(recomputed.size() == oracle.size(),
                    "evaluate_corrected.f64: size mismatch");
    CHECK_RMS("the division alone reproduces the driver's output",
              &recomputed[0], &oracle[0], (int) oracle.size(), bound);
  }

  /* In-mask the corrected volume's spread is narrower than the input's: the
   * field is removing multiplicative non-uniformity, not adding it. */
  {
    n3::Stats in  = n3::masked_stats(chunk, mask);
    n3::Stats out = n3::masked_stats(corrected, mask);
    double in_cv  = in.stddev / in.mean;
    double out_cv = out.stddev / out.mean;
    CHECK_TRUE("correcting removes spread inside the mask",
               out_cv < in_cv);
    printf("  (in-mask input CV %.5f -> output CV %.5f)\n", in_cv, out_cv);
  }

  delete_volume(corrected);
  delete_volume(field_out);
  delete loaded;
  delete_volume(mask);
  delete_volume(chunk);
  return n3check::report("evaluate");
}
