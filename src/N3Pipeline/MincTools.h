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

}  // namespace n3

#endif
