/* Optional spatial denoising for the C++ N3 pipeline.
 *
 * A thin bridge onto the `nlm` library (minc-toolkit-v2/NLM): non-local means,
 * Coupe et al. 2008.  N3 has no spatial filter of its own -- Sharpen.cc's
 * weiner() is a 1-D histogram deconvolution, not an image filter -- and the
 * bias field estimate is sensitive to noise, so the intended use is to estimate
 * the field from a denoised copy while correcting the original.  nu_correct_cxx
 * -denoise does exactly that.
 *
 * VIO_Volume in and out for the same reasons Buffers.h gives: the geometry has
 * to be carried anyway, and every other stage of the pipeline speaks it.
 * Dimensions are in file order, as everywhere else here.
 */

#ifndef N3_DENOISE_H
#define N3_DENOISE_H

#include <volume_io.h>

namespace n3 {

struct DenoiseOptions
{
  double sigma = 0.0;        /* noise sigma; 0 => estimate it from the volume */
  double beta = 1.0;         /* smoothing strength, nlm's beta               */
  int patch_radius = 1;      /* mincnlm -v                                   */
  int search_radius = 5;     /* mincnlm -d                                   */
  bool rician = false;       /* Rician rather than Gaussian noise model      */
  int threads = 4;
  bool verbose = false;
};

/* An NLM-filtered copy of `input` on the same grid, allocated here and owned by
 * the caller.  If sigma_used is non-null it receives the sigma the filter ran
 * with, which is the estimated one when opts.sigma is 0.
 *
 * `input` must be NC_DOUBLE, as everything n3::load and n3::like produce is.
 * Throws std::runtime_error if the filter rejects the parameters, and dies with
 * "built without NLM support" if N3 was configured without the nlm library. */
VIO_Volume denoise(VIO_Volume input, const DenoiseOptions &opts,
                   double *sigma_used = NULL);

}  // namespace n3

#endif
