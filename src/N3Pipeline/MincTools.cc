#include "MincTools.h"

#include <cstdio>
#include <cstdlib>

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

}  // namespace n3
