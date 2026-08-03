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

/* The same routine with float working storage: what correct_field itself
 * runs.  A verification instrument, not a feature.
 *
 * The two do not differ by rounding alone.  The relaxation stops when the mean
 * absolute update falls below 1e-10, an absolute threshold on values of order
 * 1e5, where the smallest update float can represent is about 1e-2: the
 * original therefore stops when its updates vanish into its own precision
 * while the double solve goes on converging.  That is worth 6.2e-07 relative
 * on chunk.mnc, and it is the double solve that is closer to solving the
 * equation -- testing/n3pipeline/test_extend.cc measures the residual of each,
 * 1.7e-10 against 9.0e-02. */
void extend_field_single_precision(VIO_Volume volume, VIO_Volume mask);

}  // namespace n3

#endif
