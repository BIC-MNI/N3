/* The regularized spline fit, as spline_smooth performs it.
 *
 * Both bases and both fit loops are the legacy code, linked unmodified
 * (Splines/{Spline,TBSpline}.cc, SplineSmooth/{splineSmooth,fieldIO}.cc).
 * What is here is spline_smooth_volume's choices: which basis, which domain,
 * and evaluating the result on a grid other than the one it was fitted on.
 *
 * The domain differs between the two bases, and not incidentally.  The driver
 * passes -full_support for b_spline and omits it for tp_spline
 * (nu_estimate_np_and_em.in:571-579), so splineSmooth.cc:124-130 gives the
 * first the whole volume's domain and the second the mask's bounding box.
 */

#ifndef N3_FITFIELD_H
#define N3_FITFIELD_H

#include <config.h>

#include <volume_io.h>

#include <EBTKS/Matrix.h>

#include "../Splines/Spline.h"
#include "../Splines/TBSpline.h"
#include "../SplineSmooth/splineSmooth.h"

namespace n3 {

/* A fitted field: the basis, the domain it was defined on, and the knot
 * spacing.  Those three are what the .imp file carries, and what is needed to
 * rebuild the basis on another grid. */
class Field
{
public:
  Field(Spline *spline, const DblMat &domain, double distance,
        enum spline_type type);
  ~Field();

  Spline *spline;
  DblMat domain;
  double distance;
  enum spline_type type;

private:
  Field(const Field &);
  Field &operator = (const Field &);
};

/* spline_smooth_volume (nu_estimate_np_and_em.in:567).  Only voxels where the
 * mask exceeds 0.5 contribute, every subsample-th one along each axis. */
Field *fit_field(VIO_Volume volume, VIO_Volume mask, enum spline_type type,
                 double distance, double lambda, int subsample);

/* Evaluate onto target, which need not be the grid the fit was made on: the
 * spline's coordinates are index times step from voxel (0,0,0)
 * (splineSmooth.cc:140), and ShrinkVolume keeps start, so a field fitted on
 * the estimation grid evaluates correctly at full resolution.  This is what
 * evaluate_field does through the .imp file.
 *
 * Outside the mask the result is zero, as smoothVolume does it; that is why
 * nu_evaluate then runs correct_field. */
void evaluate_field(Field *field, VIO_Volume target, VIO_Volume mask);

}  // namespace n3

#endif
