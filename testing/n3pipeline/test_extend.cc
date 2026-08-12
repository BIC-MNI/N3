/* Cycle 9: extending the field beyond the mask.
 *
 * The spline is zero outside its mask (cycle 8), so nu_evaluate runs
 * correct_field before dividing: Laplace's equation relaxed on the outside,
 * masked values held fixed, by successive over-relaxation on a sequence of
 * grids (SmoothField.cc).  The double run is measured against
 * field_extended.f64, correct_field's own recorded output.
 *
 * Getting there took one real defect.  correct_field passes NULL for
 * dim_names, which volume_io reads as ZYX and not as file order, and the
 * prolongation interpolates along the last array axis first -- so the axis
 * order decides what a voxel odd along two axes ends up with.  Fitting the
 * relaxation to the file order instead left 5.8% of the never-relaxed voxels
 * wrong and 1.5e-03 relative overall, which looks exactly like a precision
 * difference and is not one.
 *
 * The input is the oracle's own input, kept as a float volume.  Recomputing
 * it would not do: a relaxation amplifies whatever it is given.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"
#include "../../src/N3Pipeline/SmoothField.h"

#include <cmath>
#include <string>
#include <vector>

int main()
{
  std::string data = std::string(N3_DATA_DIR);
  std::string reference = std::string(N3_REFERENCE_DIR);

  VIO_Volume field = n3::load(reference + "/field_masked.mnc");
  VIO_Volume mask  = n3::load(data + "/chunk_mask.mnc");
  int n = n3::voxel_count(field);

  /* What the input looks like before extension, so that the comparison below
   * is known to be measuring an extension and not a copy. */
  {
    double *v = n3::values(field), *m = n3::values(mask);
    int outside = 0, nonzero_outside = 0;
    for(int i = 0; i < n; i++)
      if(m[i] <= 0.5) { outside++; if(v[i] != 0.0) nonzero_outside++; }
    CHECK_TRUE("the input is zero outside the mask", nonzero_outside == 0);
    printf("  (%d of %d voxels lie outside the mask)\n", outside, n);
  }

  std::vector<double> before(n3::values(field), n3::values(field) + n);

  n3::extend_field(field, mask);

  /* Not a bit-exact match to correct_field's own recorded output, and not by
   * rounding.  correct_field's relaxation stops when the mean absolute update
   * falls below 1e-10 -- an absolute threshold on values of order 1e5, where
   * the smallest update float can represent is about 1e-2, so it stops when
   * its updates vanish into its own precision.  The double solve here goes on
   * converging past that point. */
  {
    std::vector<double> oracle = n3fixture::read_f64("field_extended.f64");
    std::vector<double> mine = n3fixture::strided(n3::values(field), n);
    double drift = n3check::rel_rms(&mine[0], &oracle[0], (int) mine.size());
    CHECK_TRUE("in double it differs from correct_field by more than float "
               "rounding", drift > 1.2e-07);
    printf("  (the two solves sit %.3e apart)\n", drift);
  }

  /* Properties of a harmonic extension, which hold whatever the oracle says.
   *
   * Inside the mask nothing may move: those values are the data. */
  {
    double *v = n3::values(field), *m = n3::values(mask);
    double worst = 0.0;
    for(int i = 0; i < n; i++)
      if(m[i] > 0.5)
        { double d = fabs(v[i] - before[i]); if(d > worst) worst = d; }
    CHECK_NEAR("values inside the mask are untouched", worst, 0.0, 0.0);
  }

  /* Outside, the field must now be filled and must stay within the range of
   * the data it was extended from -- a solution of Laplace's equation attains
   * its extrema on the boundary. */
  {
    double *v = n3::values(field), *m = n3::values(mask);
    double lo = 0.0, hi = 0.0;
    bool first = true;
    for(int i = 0; i < n; i++)
      if(m[i] > 0.5)
        {
          if(first) { lo = hi = v[i]; first = false; }
          if(v[i] < lo) lo = v[i];
          if(v[i] > hi) hi = v[i];
        }

    int still_zero = 0, out_of_range = 0;
    for(int i = 0; i < n; i++)
      if(m[i] <= 0.5)
        {
          if(v[i] == 0.0) still_zero++;
          /* A tolerance of one part in 1e9 of the span, for the relaxation's
           * own residual. */
          if(v[i] < lo - 1e-9 * (hi - lo) || v[i] > hi + 1e-9 * (hi - lo))
            out_of_range++;
        }
    CHECK_TRUE("the outside is filled", still_zero == 0);
    CHECK_NEAR("and stays inside the masked range", out_of_range, 0.0, 0.0);
    printf("  (masked range %.1f .. %.1f)\n", lo, hi);
  }

  delete_volume(mask);
  delete_volume(field);
  return n3check::report("extend");
}
