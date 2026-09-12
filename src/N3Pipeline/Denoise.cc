#include "Denoise.h"

#include <config.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "Buffers.h"

#ifdef N3_WITH_NLM
#include <nlm_denoise.h>
#endif

namespace n3 {

#ifndef N3_WITH_NLM

VIO_Volume denoise(VIO_Volume, const DenoiseOptions &, double *)
{
  fprintf(stderr,
          "n3: denoising requested, but N3 was built without NLM support\n");
  exit(1);
  return NULL;  /* not reached */
}

#else

VIO_Volume denoise(VIO_Volume input, const DenoiseOptions &opts,
                   double *sigma_used)
{
  int sizes[VIO_N_DIMENSIONS];
  VIO_Real steps[VIO_N_DIMENSIONS];
  get_volume_sizes(input, sizes);
  get_volume_separations(input, steps);

  /* The pipeline's volumes are in file order, where the LAST dimension varies
   * fastest (Buffers.h); nlm::denoise indexes x-fastest.  Reversing the size
   * and step triples makes the two descriptions agree on the same flat buffer,
   * so the data itself is copied straight through -- no transpose.  Getting
   * this wrong would not crash: it would silently filter with the wrong
   * anisotropy and, for a non-cubic volume, out of bounds. */
  const int nlm_sizes[3]    = { sizes[2], sizes[1], sizes[0] };
  const double nlm_steps[3] = { (double) steps[2], (double) steps[1],
                                (double) steps[0] };

  const int count = voxel_count(input);
  const double *in = values(input);

  /* nlm works in float; N3 volumes are NC_DOUBLE. */
  std::vector<float> src(count), dst(count, 0.0f);
  for(int i = 0; i < count; i++) src[i] = (float) in[i];

  nlm::denoise_params p;
  p.sigma         = opts.sigma;
  p.beta          = opts.beta;
  p.patch_radius  = opts.patch_radius;
  p.search_radius = opts.search_radius;
  /* 0 = L2/Gaussian, 2 = L2 with the Rician bias correction. */
  p.weight_method = opts.rician ? 2 : 0;
  p.threads       = opts.threads;
  p.verbose       = opts.verbose ? 1 : 0;

  const double sigma = nlm::denoise(&src[0], &dst[0], nlm_sizes, nlm_steps, p);
  if(sigma_used) *sigma_used = sigma;

  VIO_Volume out = like(input);
  double *dest = values(out);
  for(int i = 0; i < count; i++) dest[i] = dst[i];

  return out;
}

#endif  /* N3_WITH_NLM */

}  // namespace n3
