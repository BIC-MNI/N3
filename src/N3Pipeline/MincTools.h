/* The MINC command-line tools the drivers depend on and this tree has no C++
 * for: minclookup and mincstats' bimodal threshold.
 *
 * torch_n3/minc_tools.py is a validated specification for both, measured
 * against the same binaries; the recipes here follow it.
 */

#ifndef N3_MINCTOOLS_H
#define N3_MINCTOOLS_H

#include <vector>

#include <volume_io.h>

namespace n3 {

/* minclookup -continuous: each voxel is placed against the table's own entry
 * positions, linearly interpolated between the two that bracket it, and
 * clamped outside.  In place.
 *
 * The positions are passed in rather than derived from a range because they
 * are not always the range divided evenly.  sharpen_hist writes them with six
 * decimals, so the table minclookup reads has 0.005025 where the exact
 * position is 0.005025125628; interpolating against the exact positions
 * instead moves the result by 5.1e-07 relative on chunk.mnc.  In memory the
 * positions are the bin centres and no such rounding occurs, which is one of
 * the divergences from the Perl this port exists to remove. */
void apply_lookup(VIO_Volume volume, const std::vector<double> &positions,
                  const std::vector<double> &values);

/* The same, for entries spread evenly over [lo, hi] -- what the pipeline
 * itself always has, the positions being the histogram's bin centres. */
void apply_lookup(VIO_Volume volume, const std::vector<double> &values,
                  double lo, double hi);

/* volume_stats -biModalT: an EBTKS histogram over the selected voxels' real
 * range, with as many bins as the file's voxel range has whole steps
 * (volumeStats.cc:286-299), and that class's own biModalThreshold.  Used by
 * CreateMask.
 *
 * The bin count comes from how the file was stored, so it is passed in rather
 * than derived from the buffer: on a volume held in double there is no such
 * quantity.  Legacy behaviour, reproduced and not repaired. */
double bimodal_threshold_volume_stats(VIO_Volume volume, VIO_Volume mask,
                                      int bins);

/* mincstats -biModalT: Otsu over 2000 bins, returning the winning bin's
 * centre.  A different rule from the one above, giving a different answer on
 * the same volume; nu_evaluate uses this one and CreateMask the other. */
double bimodal_threshold_mincstats(VIO_Volume volume, int bins = 2000);

}  // namespace n3

#endif
