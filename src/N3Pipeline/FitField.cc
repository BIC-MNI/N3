#include "FitField.h"

#include <cstdio>
#include <cstdlib>

#include "../SplineSmooth/fieldIO.h"

/* Declared in splineSmooth.cc, which has no header of its own for them. */
DblMat volume_domain(VIO_Volume volume);
DblMat reduced_domain(VIO_Volume volume);
void fitSplinesToVolume(Spline *spline, VIO_Volume volume, int subsample);
void fitSplinesToVolume(Spline *spline, VIO_Volume volume, VIO_Volume mask_volume,
                        const DblMat &domain, int subsample);
void fitSplinesToVolumeLookup(TBSplineVolume *spline, VIO_Volume volume,
                              int subsample);
void fitSplinesToVolumeLookup(TBSplineVolume *spline, VIO_Volume volume,
                              VIO_Volume mask_volume, const DblMat &domain,
                              int subsample);

namespace n3 {

Field::Field(Spline *s, const DblMat &d, double dist, enum spline_type t)
  : spline(s), domain(d), distance(dist), type(t) {}

Field::~Field() { delete spline; }

Field *fit_field(VIO_Volume volume, VIO_Volume mask, enum spline_type type,
                 double distance, double lambda, int subsample)
{
  /* The driver passes -full_support with -b_spline and omits it with
   * -tp_spline (nu_estimate_np_and_em.in:571-579), which is what sends the two
   * bases down splineSmooth.cc:124-130's two branches. */
  bool full_support = (type == b_spline);
  DblMat domain = (!full_support && mask) ? reduced_domain(mask)
                                          : volume_domain(volume);

  Spline *spline;
  if(type == b_spline)
    {
      VIO_Real separations[VIO_N_DIMENSIONS];
      int sizes[VIO_N_DIMENSIONS];
      VIO_Real start[VIO_N_DIMENSIONS] = { 0.0, 0.0, 0.0 };
      get_volume_separations(volume, separations);
      get_volume_sizes(volume, sizes);

      TBSplineVolume *b = new TBSplineVolume(domain, start, separations, sizes,
                                             distance, lambda);
      if(mask) fitSplinesToVolumeLookup(b, volume, mask, domain, subsample);
      else     fitSplinesToVolumeLookup(b, volume, subsample);
      spline = b;
    }
  else
    {
      spline = createThinPlateSpline(domain, distance, lambda, FALSE);
      if(mask) fitSplinesToVolume(spline, volume, mask, domain, subsample);
      else     fitSplinesToVolume(spline, volume, subsample);
    }

  return new Field(spline, domain, distance, type);
}

void evaluate_field(Field *field, VIO_Volume target, VIO_Volume mask)
{
  double real_min, real_max;

  if(field->type == b_spline)
    {
      /* The basis has to be laid out on the target's grid: TBSplineVolume
       * caches per-grid lookup tables.  The coefficients carry across because
       * the number of them depends on the domain and the knot spacing alone,
       * and the spline's coordinates are index times step from voxel (0,0,0).
       * This is what inputCompactField does (fieldIO.cc:312-326), with the
       * coefficients coming from memory instead of from a file. */
      VIO_Real separations[VIO_N_DIMENSIONS];
      int sizes[VIO_N_DIMENSIONS];
      VIO_Real start[VIO_N_DIMENSIONS] = { 0.0, 0.0, 0.0 };
      get_volume_separations(target, separations);
      get_volume_sizes(target, sizes);

      TBSplineVolume *on_target =
        new TBSplineVolume(field->domain, start, separations, sizes,
                           field->distance, 1.0, FALSE);
      if(on_target->putCoefficients(field->spline->getCoefficients()) == FALSE)
        {
          fprintf(stderr, "n3::evaluate_field: coefficient count mismatch\n");
          exit(1);
        }

      if(mask)
        smoothVolumeLookup(on_target, target, mask, &real_min, &real_max);
      else
        smoothVolumeLookup(on_target, target, &real_min, &real_max);

      delete on_target;
    }
  else
    {
      /* A thin plate spline is evaluated at a point, so it needs no per-grid
       * basis and carries across unchanged. */
      if(mask)
        smoothVolume(field->spline, target, mask, &real_min, &real_max);
      else
        smoothVolume(field->spline, target, &real_min, &real_max);
    }
}

}  // namespace n3
