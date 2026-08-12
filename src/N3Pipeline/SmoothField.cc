#include "SmoothField.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace n3 {

/* correctField.cc:25-205 over flat arrays, in double.  The sweep stays in
 * raster order: Gauss-Seidel reads values already updated in the same pass,
 * so reordering it would be a different iteration. */
static void relax(double *val, const char *mm, const int sizes[3],
                  const double seps[3])
{
  double SOR = 1.9;

  double fx = 1.0 / (seps[0] * seps[0]);
  double fy = 1.0 / (seps[1] * seps[1]);
  double fz = 1.0 / (seps[2] * seps[2]);

  /* The coarsest grid is at least 4 mm. */
  double min_sep = fabs(seps[0]);
  if(fabs(seps[1]) < min_sep) min_sep = fabs(seps[1]);
  if(fabs(seps[2]) < min_sep) min_sep = fabs(seps[2]);
  int inc = 2;
  while(inc < (int) (4.0 / min_sep)) inc *= 2;

  double thresh = 1.0e-10;
  int n_iters = 4 * 200;

  for(; inc >= 2; inc /= 2)
    {
      for(int iter = 0; iter < n_iters; iter++)
        {
          int count = 0;
          double res = 0.0;
          for(int i = 0; i < sizes[0]; i += inc)
            for(int j = 0; j < sizes[1]; j += inc)
              for(int k = 0; k < sizes[2]; k += inc)
                {
                  if(mm[(i * sizes[1] + j) * sizes[2] + k] != 0) continue;

                  double norm = 0.0;
                  double uxx = 0.0, uyy = 0.0, uzz = 0.0;
                  if(i > (inc - 1))
                    { uxx += fx * val[((i-inc)*sizes[1]+j)*sizes[2]+k]; norm += fx; }
                  if(i < sizes[0] - inc)
                    { uxx += fx * val[((i+inc)*sizes[1]+j)*sizes[2]+k]; norm += fx; }
                  if(j > (inc - 1))
                    { uyy += fy * val[(i*sizes[1]+j-inc)*sizes[2]+k]; norm += fy; }
                  if(j < sizes[1] - inc)
                    { uyy += fy * val[(i*sizes[1]+j+inc)*sizes[2]+k]; norm += fy; }
                  if(k > (inc - 1))
                    { uzz += fz * val[(i*sizes[1]+j)*sizes[2]+k-inc]; norm += fz; }
                  if(k < sizes[2] - inc)
                    { uzz += fz * val[(i*sizes[1]+j)*sizes[2]+k+inc]; norm += fz; }

                  double oldValue = val[(i*sizes[1]+j)*sizes[2]+k];
                  double tmpValue = (uxx + uyy + uzz) / norm;
                  double newValue = oldValue + SOR * (tmpValue - oldValue);
                  val[(i*sizes[1]+j)*sizes[2]+k] = newValue;
                  count++;
                  res += fabs(oldValue - newValue);
                }
          if(count > 0) res /= (double) count;
          if(res < thresh) break;
        }
      n_iters = n_iters / 2 + 1;

      if(inc == 1) break;

      /* Prolongation to the next finer level: along the last array axis, then
       * the middle one, then the first.  The order is load-bearing -- a voxel
       * odd along two axes keeps whatever the later pass gave it -- and after
       * the last level (inc == 2) the odd voxels are only interpolated, never
       * relaxed. */
      int i, j, k;
      for(i = 0; i < sizes[0]; i += inc)
        {
          for(j = 0; j < sizes[1]; j += inc)
            {
              for(k = inc/2; k < sizes[2] - inc/2; k += inc)
                if(mm[(i*sizes[1]+j)*sizes[2]+k] == 0)
                  val[(i*sizes[1]+j)*sizes[2]+k] =
                    0.5 * (val[(i*sizes[1]+j)*sizes[2]+k-inc/2] +
                           val[(i*sizes[1]+j)*sizes[2]+k+inc/2]);
              if(k < sizes[2])
                if(mm[(i*sizes[1]+j)*sizes[2]+k] == 0)
                  val[(i*sizes[1]+j)*sizes[2]+k] =
                    val[(i*sizes[1]+j)*sizes[2]+k-inc/2];
            }

          for(j = inc/2; j < sizes[1] - inc/2; j += inc)
            for(k = 0; k < sizes[2]; k += inc/2)
              if(mm[(i*sizes[1]+j)*sizes[2]+k] == 0)
                val[(i*sizes[1]+j)*sizes[2]+k] =
                  0.5 * (val[(i*sizes[1]+j-inc/2)*sizes[2]+k] +
                         val[(i*sizes[1]+j+inc/2)*sizes[2]+k]);
          if(j < sizes[1])
            for(k = 0; k < sizes[2]; k += inc/2)
              if(mm[(i*sizes[1]+j)*sizes[2]+k] == 0)
                val[(i*sizes[1]+j)*sizes[2]+k] =
                  val[(i*sizes[1]+j-inc/2)*sizes[2]+k];
        }

      for(i = inc/2; i < sizes[0] - inc/2; i += inc)
        for(j = 0; j < sizes[1]; j += inc/2)
          for(k = 0; k < sizes[2]; k += inc/2)
            if(mm[(i*sizes[1]+j)*sizes[2]+k] == 0)
              val[(i*sizes[1]+j)*sizes[2]+k] =
                0.5 * (val[((i+inc/2)*sizes[1]+j)*sizes[2]+k] +
                       val[((i-inc/2)*sizes[1]+j)*sizes[2]+k]);
      if(i < sizes[0])
        for(j = 0; j < sizes[1]; j += inc/2)
          for(k = 0; k < sizes[2]; k += inc/2)
            if(mm[(i*sizes[1]+j)*sizes[2]+k] == 0)
              val[(i*sizes[1]+j)*sizes[2]+k] =
                val[((i-inc/2)*sizes[1]+j)*sizes[2]+k];
    }
}

/* Where zspace, yspace and xspace sit in this volume's own dimension order.
 *
 * correct_field passes NULL for dim_names (correctField.cc:224) and volume_io
 * reads that as its default, which is ZYX -- not the file order that every
 * other N3 program here asks for by name.  The routine is not indifferent to
 * the choice: the prolongation interpolates along the last array axis first
 * and the first axis last, so a voxel odd along two axes takes a different
 * value under a different ordering.  On chunk.mnc, whose file order is
 * xspace zspace yspace, the two differ by 1.5e-03 relative at the 5.8% of
 * voxels the finest level never relaxes.  Reproduced, not repaired. */
static void zyx_permutation(VIO_Volume volume, int perm[3])
{
  VIO_STR *names = get_volume_dimension_names(volume);
  const char *wanted[3] = { MIzspace, MIyspace, MIxspace };

  for(int a = 0; a < 3; a++)
    {
      perm[a] = -1;
      for(int b = 0; b < VIO_N_DIMENSIONS; b++)
        if(strcmp(names[b], wanted[a]) == 0) perm[a] = b;
      if(perm[a] < 0)
        {
          fprintf(stderr, "n3::extend_field: volume has no %s\n", wanted[a]);
          exit(1);
        }
    }

  delete_dimension_names(volume, names);
}

static void extend(VIO_Volume volume, VIO_Volume mask)
{
  int perm[3];
  zyx_permutation(volume, perm);

  int own_sizes[VIO_N_DIMENSIONS];
  VIO_Real own_seps[VIO_N_DIMENSIONS];
  get_volume_sizes(volume, own_sizes);
  get_volume_separations(volume, own_seps);

  int sizes[3];
  double seps[3];
  for(int a = 0; a < 3; a++)
    { sizes[a] = own_sizes[perm[a]]; seps[a] = own_seps[perm[a]]; }

  int total = sizes[0] * sizes[1] * sizes[2];
  double *val = new double[total];
  char *mm = new char[total];

  int index[VIO_N_DIMENSIONS];
  for(int i = 0; i < sizes[0]; i++)
    for(int j = 0; j < sizes[1]; j++)
      for(int k = 0; k < sizes[2]; k++)
        {
          index[perm[0]] = i; index[perm[1]] = j; index[perm[2]] = k;
          int at = (i * sizes[1] + j) * sizes[2] + k;
          if(get_volume_real_value(mask, index[0], index[1], index[2], 0, 0) > 0.5)
            {
              val[at] = get_volume_real_value(volume, index[0], index[1],
                                              index[2], 0, 0);
              mm[at] = 1;
            }
          else
            {
              mm[at] = 0;
              val[at] = 0;
            }
        }

  relax(val, mm, sizes, seps);

  for(int i = 0; i < sizes[0]; i++)
    for(int j = 0; j < sizes[1]; j++)
      for(int k = 0; k < sizes[2]; k++)
        {
          index[perm[0]] = i; index[perm[1]] = j; index[perm[2]] = k;
          set_volume_real_value(volume, index[0], index[1], index[2], 0, 0,
                                val[(i * sizes[1] + j) * sizes[2] + k]);
        }

  delete [] val;
  delete [] mm;
}

void extend_field(VIO_Volume volume, VIO_Volume mask)
{
  extend(volume, mask);
}

}  // namespace n3
