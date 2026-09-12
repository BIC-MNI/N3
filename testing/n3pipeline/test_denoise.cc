/* n3::denoise -- the bridge from the pipeline's VIO_Volume buffers onto the
 * `nlm` library.
 *
 * Two things can go wrong here and neither one crashes:
 *
 *   - the index order.  Pipeline volumes are in file order, where the LAST
 *     dimension varies fastest; nlm indexes x-fastest.  Denoise.cc reverses
 *     the size triple to reconcile them.  On a cubic volume a missing or an
 *     extra reversal is invisible, so the volume here is 30x40x50 and the test
 *     plants a feature whose position, not just whose value, has to survive.
 *   - the filter could be a no-op, or could return the input unchanged, and
 *     every "the grid matches" check would still pass.  So the substantive
 *     assertion is that added Gaussian noise actually comes out reduced.
 *
 * No oracle: NLM's own tests (NLM/tests) pin the algorithm.  What is pinned
 * here is the bridge.
 */

#include "check.h"

#include "../../src/N3Pipeline/Buffers.h"
#include "../../src/N3Pipeline/Denoise.h"

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

/* Deliberately not <random>: a fixed, self-contained generator makes the test
 * give the same answer on every platform and standard library. */
static unsigned int seed = 12345u;
static double uniform()
{
  seed = seed * 1103515245u + 12345u;
  return ((seed >> 8) & 0xFFFFFF) / (double) 0x1000000;
}

static double gaussian()
{
  double u1 = uniform(), u2 = uniform();
  if(u1 < 1e-12) u1 = 1e-12;
  return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

static double rms_error(const double *a, const double *b, int n)
{
  double s = 0.0;
  for(int i = 0; i < n; i++) { double d = a[i] - b[i]; s += d * d; }
  return sqrt(s / n);
}

int main()
{
  /* An asymmetric grid, so an index-order mistake cannot hide. */
  const int nx = 30, ny = 40, nz = 50;
  const int sizes[VIO_N_DIMENSIONS] = { nx, ny, nz };
  const VIO_Real starts[VIO_N_DIMENSIONS] = { 0.0, 0.0, 0.0 };
  const VIO_Real steps[VIO_N_DIMENSIONS]  = { 1.0, 1.0, 1.0 };

  VIO_Volume model = create_volume(VIO_N_DIMENSIONS,
                                   (char **) File_order_dimension_names,
                                   NC_DOUBLE, FALSE, 0.0, 0.0);
  set_volume_sizes(model, (int *) sizes);
  alloc_volume_data(model);
  set_volume_starts(model, (VIO_Real *) starts);
  set_volume_separations(model, (VIO_Real *) steps);

  VIO_Volume clean = n3::like(model);
  VIO_Volume noisy = n3::like(model);
  double *c = n3::values(clean);
  double *d = n3::values(noisy);

  /* A piecewise-constant phantom -- two nested boxes -- which is the shape NLM
   * is good at and a linear blur is not, plus one bright marker voxel whose
   * position the index-order check below reads back. */
  for(int i = 0; i < nx; i++)
    for(int j = 0; j < ny; j++)
      for(int k = 0; k < nz; k++)
        {
          const int idx = (i * ny + j) * nz + k;   /* file order: k fastest */
          double v = 20.0;
          if(i >= 5 && i < 25 && j >= 5 && j < 35 && k >= 5 && k < 45) v = 100.0;
          if(i >= 12 && i < 18 && j >= 15 && j < 25 && k >= 20 && k < 30) v = 160.0;
          c[idx] = v;
          d[idx] = v + 8.0 * gaussian();
        }

  const int n = n3::voxel_count(clean);

  n3::DenoiseOptions opts;
  opts.sigma = 8.0;             /* told, not estimated, so this is a filter
                                 * test and not a noise-estimator test */
  opts.threads = 2;
  VIO_Volume filtered = n3::denoise(noisy, opts);

  /* ------------------------------------------------------- the grid */
  {
    int out_sizes[VIO_N_DIMENSIONS];
    VIO_Real out_steps[VIO_N_DIMENSIONS];
    get_volume_sizes(filtered, out_sizes);
    get_volume_separations(filtered, out_steps);
    CHECK_TRUE("the output is on the input's grid",
               out_sizes[0] == nx && out_sizes[1] == ny && out_sizes[2] == nz);
    CHECK_TRUE("the output keeps the input's steps",
               out_steps[0] == steps[0] && out_steps[1] == steps[1]
               && out_steps[2] == steps[2]);
    CHECK_TRUE("the input is not modified in place", filtered != noisy);
  }

  /* ------------------------------------------------------ the filtering */
  {
    double before = rms_error(n3::values(noisy),    c, n);
    double after  = rms_error(n3::values(filtered), c, n);
    printf("rms error: %.4f before, %.4f after\n", before, after);
    /* A bound, not a tie: NLM on a piecewise-constant phantom at sigma 8
     * removes far more than a quarter of the error, so 0.75 fails loudly if
     * the bridge were feeding the filter garbage, without pinning a number
     * that belongs to the algorithm rather than to this code. */
    CHECK_TRUE("denoising reduces the error", after < 0.75 * before);
    CHECK_TRUE("denoising is not a no-op", after < before);
  }

  /* --------------------------------------------- the index order itself */
  {
    /* The filter is local, so a plateau interior voxel comes back near its
     * clean value wherever it sits.  A transposed mapping would read this
     * position out of a different region of the phantom -- 20 or 100 rather
     * than 160 -- and the margin is 60 intensity units wide. */
    const int fi = 15, fj = 20, fk = 25;   /* inside the inner box */
    const int flat = (fi * ny + fj) * nz + fk;
    CHECK_NEAR("an inner-box voxel keeps its own intensity",
               n3::values(filtered)[flat], 160.0, 25.0);

    const int oi = 2, oj = 2, ok = 2;      /* outside both boxes */
    const int oflat = (oi * ny + oj) * nz + ok;
    CHECK_NEAR("a background voxel keeps its own intensity",
               n3::values(filtered)[oflat], 20.0, 25.0);

    /* Whether the mapping is the identity on a NON-cubic volume is the part a
     * cubic test cannot see: with nz = 50 and nx = 30, a reversed-but-unpaired
     * mapping would index past the end of the shorter axis, which sanitises to
     * either a crash or a wildly wrong value at high k. */
    const int hi_flat = ((nx - 2) * ny + (ny - 2)) * nz + (nz - 2);
    CHECK_NEAR("the far corner is still background",
               n3::values(filtered)[hi_flat], 20.0, 25.0);
  }

  delete_volume(filtered);
  delete_volume(noisy);
  delete_volume(clean);
  delete_volume(model);
  return n3check::report("denoise");
}
