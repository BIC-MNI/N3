#include "MincTools.h"

#include <cstdio>
#include <cstdlib>

#include <config.h>
#include <EBTKS/Histogram.h>

namespace n3 {

void apply_lookup(VIO_Volume volume, const std::vector<double> &positions,
                  const std::vector<double> &values)
{
  int entries = (int) values.size();
  if(entries < 1) { fprintf(stderr, "n3::apply_lookup: empty table\n"); exit(1); }
  if((int) positions.size() != entries)
    { fprintf(stderr, "n3::apply_lookup: positions and values differ\n"); exit(1); }

  int sizes[VIO_N_DIMENSIONS];
  get_volume_sizes(volume, sizes);
  int n = sizes[0] * sizes[1] * sizes[2];

  void *raw;
  GET_VOXEL_PTR_3D(raw, volume, 0, 0, 0);
  double *v = (double *) raw;

  if(entries == 1)
    {
      for(int i = 0; i < n; i++) v[i] = values[0];
      return;
    }

  double lo = positions[0], hi = positions[entries - 1];

  for(int i = 0; i < n; i++)
    {
      double value = v[i];

      if(value <= lo) { v[i] = values[0]; continue; }
      if(value >= hi) { v[i] = values[entries - 1]; continue; }

      /* The last entry not past the value.  Binary search rather than
       * arithmetic on the index, since the positions need not be evenly
       * spaced -- and where they are, this is the same answer. */
      int low = 0, high = entries - 1;
      while(high - low > 1)
        {
          int mid = (low + high) / 2;
          if(positions[mid] <= value) low = mid; else high = mid;
        }

      double slope = (values[low + 1] - values[low])
                     / (positions[low + 1] - positions[low]);
      v[i] = values[low] + slope * (value - positions[low]);
    }
}

void apply_lookup(VIO_Volume volume, const std::vector<double> &values,
                  double lo, double hi)
{
  int entries = (int) values.size();
  if(entries < 2) { apply_lookup(volume, values, values); return; }
  if(hi == lo) { fprintf(stderr, "n3::apply_lookup: empty range\n"); exit(1); }

  std::vector<double> positions(entries);
  double step = (hi - lo) / (entries - 1);
  for(int i = 0; i < entries; i++) positions[i] = lo + i * step;
  positions[entries - 1] = hi;

  apply_lookup(volume, positions, values);
}

double bimodal_threshold_volume_stats(VIO_Volume volume, VIO_Volume mask,
                                      int bins)
{
  int sizes[VIO_N_DIMENSIONS];
  get_volume_sizes(volume, sizes);
  int n = sizes[0] * sizes[1] * sizes[2];

  void *raw;
  GET_VOXEL_PTR_3D(raw, volume, 0, 0, 0);
  double *v = (double *) raw;
  double *m = NULL;
  if(mask) { GET_VOXEL_PTR_3D(raw, mask, 0, 0, 0); m = (double *) raw; }

  /* volume_stats histograms the values it selected, over their own extrema
   * (volumeStats.cc:239-259, 287). */
  bool any = false;
  double lo = 0.0, hi = 0.0;
  for(int i = 0; i < n; i++)
    {
      if(m && m[i] == 0.0) continue;
      if(!any) { lo = hi = v[i]; any = true; }
      if(v[i] < lo) lo = v[i];
      if(v[i] > hi) hi = v[i];
    }
  if(!any) return 0.0;

  /* EBTKS' own class, so the threshold is the legacy arithmetic and not a
   * second implementation of it. */
  Histogram hist(lo, hi, (unsigned) bins);
  for(int i = 0; i < n; i++)
    {
      if(m && m[i] == 0.0) continue;
      hist.add(v[i]);
    }

  return hist.biModalThreshold();
}

double bimodal_threshold_mincstats(VIO_Volume volume, int bins)
{
  int sizes[VIO_N_DIMENSIONS];
  get_volume_sizes(volume, sizes);
  int n = sizes[0] * sizes[1] * sizes[2];

  void *raw;
  GET_VOXEL_PTR_3D(raw, volume, 0, 0, 0);
  double *v = (double *) raw;

  double lo = v[0], hi = v[0];
  for(int i = 1; i < n; i++)
    { if(v[i] < lo) lo = v[i]; if(v[i] > hi) hi = v[i]; }
  if(hi <= lo) return lo;

  double width = (hi - lo) / bins;
  std::vector<double> counts(bins, 0.0), centres(bins);
  for(int i = 0; i < bins; i++) centres[i] = lo + (i + 0.5) * width;

  for(int i = 0; i < n; i++)
    {
      int index = (int) ((v[i] - lo) / width);
      if(index < 0) index = 0;
      if(index > bins - 1) index = bins - 1;
      counts[index] += 1.0;
    }

  /* Otsu: the split maximising the variance between the two groups.  The
   * answer reported is the winning bin's centre, which is what mincstats
   * prints and what nu_evaluate thresholds on. */
  double total = 0.0, total_sum = 0.0;
  for(int i = 0; i < bins; i++)
    { total += counts[i]; total_sum += counts[i] * centres[i]; }

  double weight_low = 0.0, sum_low = 0.0;
  double best = -1.0;
  int best_bin = 0;
  for(int i = 0; i < bins; i++)
    {
      weight_low += counts[i];
      sum_low += counts[i] * centres[i];
      double weight_high = total - weight_low;
      if(weight_low <= 0.0 || weight_high <= 0.0) continue;
      double difference = sum_low / weight_low - (total_sum - sum_low) / weight_high;
      double between = weight_low * weight_high * difference * difference;
      if(between > best) { best = between; best_bin = i; }
    }

  return centres[best_bin];
}

}  // namespace n3
