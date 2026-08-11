#include "NuEvaluate.h"

#include <cstdio>
#include <cstdlib>

#include "Buffers.h"
#include "MincTools.h"
#include "SmoothField.h"

namespace n3 {

VIO_Volume nu_evaluate(VIO_Volume input, VIO_Volume user_mask, Field *field,
                       const EvaluateOptions &options,
                       VIO_Volume *field_volume)
{
  int n = voxel_count(input);

  /* 1. A mask, if none was given: mincstats' bimodal threshold, and >= it
   * (nu_evaluate.in:49-56).  Not volume_stats' rule, which CreateMask uses. */
  VIO_Volume mask = like(input);
  {
    double *m = values(mask), *v = values(input);
    if(user_mask)
      {
        /* evaluate_field requires the mask and the -like volume to be the
         * same size, with no resampling (evaluateField.cc's compareVolumes,
         * "Mask volume and input volume must be the same size."); real
         * nu_correct crashes exactly this way if -mask is not already on the
         * input's grid. Matched here rather than silently reading past
         * user_mask's buffer, which a size mismatch would otherwise do. */
        if(voxel_count(user_mask) != n)
          {
            fprintf(stderr,
                    "n3::nu_evaluate: mask and input volume must be the "
                    "same size.\n");
            exit(1);
          }
        double *u = values(user_mask);
        for(int i = 0; i < n; i++) m[i] = (u[i] != 0.0) ? 1.0 : 0.0;
      }
    else
      {
        double threshold = bimodal_threshold_mincstats(input);
        if(options.verbose) printf("bimodal threshold: %g\n", threshold);
        for(int i = 0; i < n; i++) m[i] = (v[i] >= threshold) ? 1.0 : 0.0;
      }
  }

  /* 2. The field on the input's grid.  The spline was fitted on the estimation
   * grid; its coordinates are index times step from voxel (0,0,0), which
   * ShrinkVolume preserved. */
  VIO_Volume evaluated = like(input);
  evaluate_field(field, evaluated, mask);

  /* 3. correct_field: the spline is zero outside its mask, so the field has to
   * be extended before anything is divided by it. */
  extend_field(evaluated, mask);

  /* 4. A floor under the field, and only if something fell below it
   * (:71-76). */
  {
    double *f = values(evaluated);
    double lo = f[0];
    for(int i = 1; i < n; i++) if(f[i] < lo) lo = f[i];
    if(lo < options.field_floor)
      for(int i = 0; i < n; i++)
        if(f[i] < options.field_floor) f[i] = options.field_floor;
  }

  /* 5. output = input / field (:78).  mincmath -zero leaves a zero divisor as
   * zero rather than as an infinity. */
  VIO_Volume output = like(input);
  {
    double *o = values(output), *v = values(input), *f = values(evaluated);
    for(int i = 0; i < n; i++) o[i] = (f[i] == 0.0) ? 0.0 : v[i] / f[i];
  }

  delete_volume(mask);
  if(field_volume) *field_volume = evaluated;
  else delete_volume(evaluated);

  return output;
}

}  // namespace n3
