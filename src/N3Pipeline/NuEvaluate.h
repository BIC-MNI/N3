/* nu_evaluate, in memory.
 *
 * The Perl driver (nu_evaluate.in:46-80): build a mask if none was given,
 * evaluate the field on the input's grid, extend it beyond the mask, put a
 * floor under it, and divide.
 */

#ifndef N3_NUEVALUATE_H
#define N3_NUEVALUATE_H

#include "FitField.h"

namespace n3 {

struct EvaluateOptions
{
  double field_floor;   /* -floor, nu_evaluate.in:300 */
  bool verbose;

  EvaluateOptions() : field_floor(0.1), verbose(false) {}
};

/* The corrected volume, and optionally the field itself.
 *
 * With no mask, one is made by thresholding the input at mincstats' bimodal
 * threshold (:49-56) -- a different rule from the one CreateMask uses, and
 * both are needed.
 *
 * The caller owns both results. */
VIO_Volume nu_evaluate(VIO_Volume input, VIO_Volume user_mask, Field *field,
                       const EvaluateOptions &options,
                       VIO_Volume *field_volume = NULL);

}  // namespace n3

#endif
