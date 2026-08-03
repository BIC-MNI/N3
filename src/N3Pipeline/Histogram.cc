#include "Histogram.h"

#include <cstdio>
#include <cstdlib>

#include "Buffers.h"

#include "../VolumeHist/WHistogram.h"
#include "../VolumeHist/GHistogram.h"

namespace n3 {

/* minchist.cc rounds the mask value and uses it as a class index; the driver
 * always asks for class 1. */
static inline int label_of(double value)
{
  return (int) (value < 0.0 ? value - 0.5 : value + 0.5);
}

void auto_range(VIO_Volume volume, VIO_Volume mask, int label,
                double *min_value, double *max_value)
{
  int sizes[VIO_N_DIMENSIONS];
  get_volume_sizes(volume, sizes);
  int n = sizes[0] * sizes[1] * sizes[2];

  double *data = values(volume);
  double *mask_data = mask ? values(mask) : NULL;

  /* Seeded from the whole volume's extrema, as minchist.cc:157-158 does.  For
   * a file volume_io reports those; here they are the data's own, which is the
   * same number. */
  double lo = data[0], hi = data[0];
  for(int i = 1; i < n; i++)
    { if(data[i] < lo) lo = data[i]; if(data[i] > hi) hi = data[i]; }

  /* class_min starts at the volume maximum and class_max at its minimum, and
   * the update is an else-if (minchist.cc:172-175), so a selected voxel
   * sitting exactly at the volume maximum moves the upper bound rather than
   * the lower.  Reproduced, not tidied. */
  double class_min = hi, class_max = lo;
  for(int i = 0; i < n; i++)
    {
      if(mask_data && label_of(mask_data[i]) != label) continue;
      double value = data[i];
      if(value < class_min) class_min = value;
      else if(value > class_max) class_max = value;
    }

  if(class_max <= class_min) class_max = class_min + 1.0;   /* :193-194 */

  *min_value = class_min;
  *max_value = class_max;
}

DHistogram *histogram(VIO_Volume volume, VIO_Volume mask, int label,
                      const HistogramOptions &options)
{
  double lo, hi;
  auto_range(volume, mask, label, &lo, &hi);

  /* new_histogram's precedence: a stated width wins over the linear split
   * (minchist.cc, new_histogram). */
  DHistogram *h;
  if(options.parzen_sigma > 0.0)
    h = (DHistogram *) new GHistogram(lo, hi, (unsigned) options.bins,
                                      options.parzen_sigma);
  else if(options.window)
    h = (DHistogram *) new WHistogram(lo, hi, (unsigned) options.bins);
  else
    h = new DHistogram(lo, hi, (unsigned) options.bins);

  int sizes[VIO_N_DIMENSIONS];
  get_volume_sizes(volume, sizes);
  int n = sizes[0] * sizes[1] * sizes[2];

  double *data = values(volume);
  double *mask_data = mask ? values(mask) : NULL;

  for(int i = 0; i < n; i++)
    {
      if(mask_data && label_of(mask_data[i]) != label) continue;
      h->add(data[i]);
    }

  return h;
}

}  // namespace n3
