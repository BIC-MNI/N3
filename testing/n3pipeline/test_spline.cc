/* Cycle 8: the regularized spline fit.
 *
 * Both bases and both fit loops are the legacy code; what is tested is
 * spline_smooth_volume's choices around them, and evaluating a fit on a grid
 * other than the one it was made on -- which is what makes an estimation grid
 * and a full-resolution output possible at all.
 *
 * The oracle's input is a float copy of chunk.mnc, because spline_smooth reads
 * every volume as float whatever its storage (splineSmooth.cc:101).  Both
 * sides therefore fit identical numbers and the only difference left is the
 * float the oracle's output was written in, one ulp of which is 1.2e-7
 * relative.
 *
 * Fields are compared and coefficients are not.  The normal equations are
 * nearly singular at 200 mm -- the knots are further apart than the volume is
 * wide -- so no solver determines the coefficients to better than about 1e-4
 * while the field they describe is determined to 1e-6.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"
#include "../../src/N3Pipeline/FitField.h"

#include <cmath>
#include <string>
#include <vector>

/* One ulp of float32, which is what the oracle's storage type costs. */
static const double FLOAT_ULP = 1.2e-7;

static void compare(const char *what, VIO_Volume mine, const char *fixture)
{
  std::vector<double> oracle = n3fixture::read_f64(fixture);
  std::vector<double> ours =
    n3fixture::strided(n3::values(mine), n3::voxel_count(mine));
  n3fixture::must(ours.size() == oracle.size(),
                  std::string(fixture) + ": size mismatch");
  CHECK_RMS(what, &ours[0], &oracle[0], (int) ours.size(), FLOAT_ULP);
}

int main()
{
  std::string data = std::string(N3_DATA_DIR);
  std::string reference = std::string(N3_REFERENCE_DIR);

  VIO_Volume chunk = n3::load(reference + "/chunk_float.mnc");
  VIO_Volume mask  = n3::load(data + "/chunk_mask.mnc");

  /* The two bases are given different domains, and not by accident: the driver
   * passes -full_support for one and not the other, so one covers the volume
   * and the other the mask's bounding box.
   *
   * chunk_mask reaches every face of chunk.mnc, so on it the two rules happen
   * to agree; a mask trimmed along the first axis separates them.  Asserting
   * the rule rather than a difference on one volume is the point -- the
   * difference is data-dependent and the rule is not. */
  {
    VIO_Volume trimmed = n3::like(mask);
    int sizes[VIO_N_DIMENSIONS];
    get_volume_sizes(mask, sizes);
    double *t_data = n3::values(trimmed), *m_data = n3::values(mask);
    for(int i = 0; i < sizes[0]; i++)
      for(int j = 0; j < sizes[1]; j++)
        for(int k = 0; k < sizes[2]; k++)
          {
            int index = (i * sizes[1] + j) * sizes[2] + k;
            t_data[index] = (i >= 10 && i <= 80) ? m_data[index] : 0.0;
          }

    n3::Field *b = n3::fit_field(chunk, trimmed, b_spline, 200.0, 1e-7, 1);
    n3::Field *t = n3::fit_field(chunk, trimmed, thin_plate_spline, 200.0,
                                 1e-7, 1);

    VIO_Real separations[VIO_N_DIMENSIONS];
    get_volume_separations(chunk, separations);

    CHECK_NEAR("b_spline covers the whole volume, as -full_support asks",
               b->domain(0, 1), (sizes[0] - 0.5) * separations[0], 0.0);
    CHECK_NEAR("tp_spline covers the mask's bounding box only",
               t->domain(0, 1), 80.5 * separations[0], 0.0);
    printf("  (b_spline x %.1f..%.1f, tp_spline x %.1f..%.1f)\n",
           b->domain(0, 0), b->domain(0, 1), t->domain(0, 0), t->domain(0, 1));

    delete b;
    delete t;
    delete_volume(trimmed);
  }

  /* B-spline at three knot spacings. */
  {
    const double distances[] = { 200.0, 100.0, 50.0 };
    const char *names[] = { "fit_b200.f64", "fit_b100.f64", "fit_b50.f64" };
    for(int d = 0; d < 3; d++)
      {
        n3::Field *field =
          n3::fit_field(chunk, mask, b_spline, distances[d], 1e-7, 1);
        VIO_Volume out = n3::like(chunk);
        n3::evaluate_field(field, out, mask);

        char what[160];
        sprintf(what, "b_spline at %g mm is spline_smooth's", distances[d]);
        compare(what, out, names[d]);

        delete_volume(out);
        delete field;
      }
  }

  /* Thin plate spline, on its own domain. */
  {
    n3::Field *field = n3::fit_field(chunk, mask, thin_plate_spline, 200.0,
                                     1e-7, 1);
    VIO_Volume out = n3::like(chunk);
    n3::evaluate_field(field, out, mask);
    compare("tp_spline at 200 mm is spline_smooth's", out, "fit_tp200.f64");
    delete_volume(out);
    delete field;
  }

  /* Fitted on the estimation grid, evaluated at full resolution.  This is the
   * path the pipeline takes and the one the .imp file exists for; the oracle
   * is spline_smooth -compact followed by evaluate_field -like.
   *
   * Nothing is asserted about the coefficients themselves; see the header
   * comment. */
  {
    n3::Field *field = n3::fit_field(chunk, mask, b_spline, 200.0, 1e-7, 1);
    VIO_Volume out = n3::like(chunk);
    /* evaluate_field was given no mask, so the whole grid is evaluated. */
    n3::evaluate_field(field, out, NULL);
    compare("a fit evaluated over the whole grid is evaluate_field's",
            out, "field_b200_full.f64");
    delete_volume(out);
    delete field;
  }

  /* Two properties that hold whatever the oracle says.
   *
   * A constant inside the mask must come back as that constant: the basis
   * spans constants, so any regularization weight leaves them untouched. */
  {
    VIO_Volume flat = n3::like(chunk);
    double *v = n3::values(flat);
    int n = n3::voxel_count(flat);
    for(int i = 0; i < n; i++) v[i] = 7.5;

    n3::Field *field = n3::fit_field(flat, mask, b_spline, 200.0, 1e-7, 1);
    VIO_Volume out = n3::like(chunk);
    n3::evaluate_field(field, out, mask);

    double *o = n3::values(out);
    double *m = n3::values(mask);
    double worst = 0.0;
    for(int i = 0; i < n; i++)
      if(m[i] > 0.5) { double d = fabs(o[i] - 7.5); if(d > worst) worst = d; }
    /* To the accuracy the fit has, which is not the accuracy of double: the
     * normal equations are nearly singular at 200 mm and the field they
     * describe is determined to about 1e-6 relative. */
    CHECK_NEAR("a constant is fitted to the accuracy the fit has",
               worst, 0.0, 1e-6 * 7.5);

    delete_volume(out);
    delete field;
    delete_volume(flat);
  }

  /* The spline is zero outside its mask, which is why nu_evaluate has to run
   * correct_field before dividing by the field. */
  {
    n3::Field *field = n3::fit_field(chunk, mask, b_spline, 200.0, 1e-7, 1);
    VIO_Volume out = n3::like(chunk);
    n3::evaluate_field(field, out, mask);

    double *o = n3::values(out);
    double *m = n3::values(mask);
    int n = n3::voxel_count(out);
    bool zero_outside = true;
    for(int i = 0; i < n; i++)
      if(m[i] <= 0.5 && o[i] != 0.0) zero_outside = false;
    CHECK_TRUE("the field is zero outside the mask", zero_outside);

    delete_volume(out);
    delete field;
  }

  delete_volume(mask);
  delete_volume(chunk);
  return n3check::report("spline");
}
