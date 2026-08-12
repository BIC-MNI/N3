/* correct_field's smooth(), in double.
 *
 * The spline evaluates to zero outside its mask, so nu_evaluate extends the
 * field beyond the mask before dividing by it: Laplace's equation relaxed on
 * the outside with the masked values held fixed, by successive over-relaxation
 * on a sequence of grids.
 *
 * Transcribed rather than linked, because the original does its whole solve on
 * float arrays (correctField.cc:39-40) and this pipeline is in double.  The
 * sweep order is kept exactly -- raster order, which is what makes Gauss-Seidel
 * what it is -- so the only difference from the original is the precision.
 */

#ifndef N3_SMOOTHFIELD_H
#define N3_SMOOTHFIELD_H

#include <volume_io.h>

namespace n3 {

/* In place.  Voxels where the mask exceeds 0.5 are kept; the rest are
 * replaced by the harmonic extension of them. */
void extend_field(VIO_Volume volume, VIO_Volume mask);

}  // namespace n3

#endif
