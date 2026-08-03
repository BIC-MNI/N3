/* The masked histogram, as volume_hist builds it.
 *
 * The three estimators themselves are the legacy classes, linked unmodified
 * (VolumeHist/{DHistogram,WHistogram,GHistogram}.cc).  What is here is what
 * minchist.cc's main() does around them: choose the class, take the automatic
 * range over the selected voxels, and bin them.
 */

#ifndef N3_HISTOGRAM_H
#define N3_HISTOGRAM_H

/* config.h first, as every translation unit here does: EBTKS' SimpleArray.h
 * refuses to compile without HAVE_ISFINITE or HAVE_FINITE. */
#include <config.h>

#include <volume_io.h>

#include "../VolumeHist/DHistogram.h"

namespace n3 {

struct HistogramOptions
{
  int bins;
  bool window;          /* -parzen/-window: the linear split between the two
                         * nearest bin centres, WHistogram.h:65-93 */
  double parzen_sigma;  /* -gaussian_window: this tree's Gaussian window, off
                         * at 0.  A stated width wins over -window. */

  HistogramOptions() : bins(200), window(false), parzen_sigma(0.0) {}
};

/* The automatic range over the voxels where the mask rounds to label
 * (minchist.cc:150-197).  Seeded from the whole volume's extrema and updated
 * with an else-if, so a selected voxel sitting exactly at the volume maximum
 * moves the upper bound and not the lower one -- reproduced rather than
 * tidied. */
void auto_range(VIO_Volume volume, VIO_Volume mask, int label,
                double *min_value, double *max_value);

/* The histogram of those same voxels.  The caller owns the result. */
DHistogram *histogram(VIO_Volume volume, VIO_Volume mask, int label,
                      const HistogramOptions &options);

}  // namespace n3

#endif
