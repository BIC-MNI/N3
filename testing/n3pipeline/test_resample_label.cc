/* Cycle 2: the mask resample.
 *
 * CheckSampling sends a label volume through resample_labels, which is
 * trilinear thresholded at 0.5 (resample_labels.in:176-177) -- not nearest
 * neighbour.  At an integer -shrink the estimation grid falls on input voxel
 * centres and the two agree exactly, which is how the difference stays hidden;
 * factor 2.5 puts the samples between voxels and separates them.
 *
 * The oracle is the driver's own resampled mask, kept by -tmpdir/-keeptmp.
 * It is compared as labels rather than as numbers: resample_labels writes a
 * byte volume whose 1 comes back as 1 - 6e-16, and a tolerance chosen to
 * absorb that would be a tolerance on nothing.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"

#include <string>
#include <vector>

static const char *FACTORS[] = { "3", "2.5" };
static const double FACTOR_VALUE[] = { 3.0, 2.5 };
static const int N_FACTORS = 2;

static void check_geometry(const char *what, VIO_Volume volume,
                           const std::vector<double> &geom)
{
  int sizes[VIO_N_DIMENSIONS];
  VIO_Real starts[VIO_N_DIMENSIONS], steps[VIO_N_DIMENSIONS];
  get_volume_sizes(volume, sizes);
  get_volume_starts(volume, starts);
  get_volume_separations(volume, steps);

  bool ok = geom.size() == 9;
  for(int i = 0; ok && i < VIO_N_DIMENSIONS; i++)
    ok = (sizes[i] == (int) geom[i]) && (starts[i] == geom[3 + i])
         && (steps[i] == geom[6 + i]);
  CHECK_TRUE(what, ok);
}

int main()
{
  std::string data = std::string(N3_DATA_DIR);
  VIO_Volume chunk = n3::load(data + "/chunk.mnc.gz");
  VIO_Volume mask  = n3::load(data + "/chunk_mask.mnc.gz");

  /* A mask resampled onto its own grid is that mask.  True whatever the
   * interpolation, so it holds before any oracle is consulted. */
  {
    VIO_Volume same = n3::resample_label(mask, mask);
    int n = n3::voxel_count(mask);
    double *a = n3::values(same), *b = n3::values(mask);
    int differing = 0;
    for(int i = 0; i < n; i++) if((a[i] != 0.0) != (b[i] != 0.0)) differing++;
    CHECK_NEAR("resampling a mask onto its own grid changes nothing",
               differing, 0.0, 0.0);
    delete_volume(same);
  }

  for(int f = 0; f < N_FACTORS; f++)
    {
      char what[160];
      std::string name = std::string("mask") + FACTORS[f];

      VIO_Volume model = n3::shrink(chunk, FACTOR_VALUE[f]);
      VIO_Volume small = n3::resample_label(mask, model);

      sprintf(what, "shrink %s: the resampled mask is on the estimation grid",
              FACTORS[f]);
      check_geometry(what, small, n3fixture::read_text(name + ".geom"));

      std::vector<double> oracle = n3fixture::read_f64(name + ".f64");
      int n = n3::voxel_count(small);
      n3fixture::must((int) oracle.size() == n, "size mismatch against oracle");

      double *mine = n3::values(small);
      int differing = 0, inside = 0;
      for(int i = 0; i < n; i++)
        {
          bool theirs = oracle[i] >= 0.5;
          if(theirs) inside++;
          if((mine[i] != 0.0) != theirs) differing++;
        }
      sprintf(what, "shrink %s: %d of %d labels are the driver's own",
              FACTORS[f], n - differing, n);
      CHECK_NEAR(what, differing, 0.0, 0.0);

      /* Nearest neighbour is what a port reaches for here, and at an integer
       * factor it is indistinguishable.  Establish which case each factor is,
       * so that a future change to either routine cannot quietly swap them. */
      VIO_Volume nearest = n3::shrink(mask, FACTOR_VALUE[f]);
      double *nn = n3::values(nearest);
      int nn_differing = 0;
      for(int i = 0; i < n; i++)
        if((nn[i] != 0.0) != (oracle[i] >= 0.5)) nn_differing++;
      sprintf(what, "shrink %s: nearest neighbour differs on %d voxels",
              FACTORS[f], nn_differing);
      CHECK_TRUE(what, FACTOR_VALUE[f] == 3.0 ? nn_differing == 0
                                              : nn_differing > 0);
      printf("  (%d of %d voxels inside the mask)\n", inside, n);

      delete_volume(nearest);
      delete_volume(small);
      delete_volume(model);
    }

  delete_volume(mask);
  delete_volume(chunk);
  return n3check::report("resample_label");
}
