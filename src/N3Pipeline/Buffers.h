/* In-memory volume buffers for the C++ N3 pipeline.
 *
 * Every intermediate the Perl drivers pass between programs as a MINC file is
 * one of these instead: a volume_io volume of type NC_DOUBLE, allocated in
 * memory and written out only when the caller asks for a result.
 *
 * VIO_Volume rather than a plain array because the reused routines
 * (fitSplinesToVolumeLookup, smoothVolumeLookup, volume_domain,
 * outputCompactField) take one, and because the geometry has to be carried
 * anyway.  Dimensions are in file order, which is what the legacy programs use
 * (File_order_dimension_names, fieldIO.cc:101) and therefore what the spline's
 * index-by-step coordinate system is expressed in.
 */

#ifndef N3_BUFFERS_H
#define N3_BUFFERS_H

#include <string>

#include <volume_io.h>

namespace n3 {

/* Read a MINC file as NC_DOUBLE in file order.  Handles .gz, as every
 * volume_io program here does. */
VIO_Volume load(const std::string &path);

/* An NC_DOUBLE volume with the model's dimension names and direction cosines
 * but the given sizes, starts and steps, zero filled. */
VIO_Volume create_like(VIO_Volume model,
                       const int sizes[VIO_N_DIMENSIONS],
                       const VIO_Real starts[VIO_N_DIMENSIONS],
                       const VIO_Real steps[VIO_N_DIMENSIONS]);

/* The same grid as model, zero filled. */
VIO_Volume like(VIO_Volume model);

/* Flat access.  An NC_DOUBLE volume created here has no voxel-to-real scaling,
 * so the flat buffer holds real values; values() checks that rather than
 * assuming it. */
double *values(VIO_Volume volume);
int voxel_count(VIO_Volume volume);

/* ShrinkVolume, nu_estimate_np_and_em.in:955.  Keeps start, multiplies step by
 * the factor and takes ceil((n-1)/factor)+1 samples, on axes finer than
 * factor * min|step| only, sampled nearest neighbour.  Samples that fall
 * outside the input are zero, which is what mincresample fills them with. */
VIO_Volume shrink(VIO_Volume in, double factor);

/* CheckSampling's label branch (nu_estimate_np_and_em.in:861-868), which
 * drives resample_labels: trilinear onto the model's grid, then thresholded at
 * 0.5 (resample_labels.in:176-177).  Not nearest neighbour -- at an integer
 * -shrink the two coincide, which is what makes the difference easy to miss.
 * Samples outside the input contribute zero, as mincresample fills them. */
VIO_Volume resample_label(VIO_Volume in, VIO_Volume model);

/* What volume_stats prints, over the voxels where mask is non-zero
 * (volumeStats.cc:236).  A null mask takes the whole volume.  The standard
 * deviation is the population one, sqrt(E[x^2] - E[x]^2) (:276) -- the same
 * quantity the stopping rule needs, and not the default of most libraries. */
struct Stats
{
  double mean, stddev, minimum, maximum;
  int count;
};

Stats masked_stats(VIO_Volume volume, VIO_Volume mask);

/* The storage type of a file, which an output has to carry to be what
 * `mincmath -copy_header` would have written. */
nc_type storage_type(const std::string &path, VIO_BOOL *signed_flag);

/* Write, taking the header from like_path as -copy_header does. */
void save(VIO_Volume volume, const std::string &path,
          const std::string &like_path, nc_type type, VIO_BOOL signed_flag,
          const std::string &history);

}  // namespace n3

#endif
