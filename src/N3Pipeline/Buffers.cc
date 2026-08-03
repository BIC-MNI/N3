#include "Buffers.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace n3 {

static void fail(const char *what, const std::string &detail = "")
{
  fprintf(stderr, "n3: %s%s%s\n", what, detail.empty() ? "" : ": ", detail.c_str());
  exit(1);
}

VIO_Volume load(const std::string &path)
{
  VIO_Volume volume;

  if(input_volume((char *) path.c_str(), VIO_N_DIMENSIONS,
                  File_order_dimension_names,
                  NC_DOUBLE, FALSE, 0.0, 0.0, TRUE, &volume,
                  (minc_input_options *) NULL) != VIO_OK)
    fail("cannot read volume", path);

  return volume;
}

VIO_Volume create_like(VIO_Volume model, const int sizes[VIO_N_DIMENSIONS],
                       const VIO_Real starts[VIO_N_DIMENSIONS],
                       const VIO_Real steps[VIO_N_DIMENSIONS])
{
  VIO_Volume volume = copy_volume_definition_no_alloc(model, NC_DOUBLE, FALSE,
                                                      0.0, 0.0);
  set_volume_sizes(volume, (int *) sizes);
  alloc_volume_data(volume);
  set_volume_starts(volume, (VIO_Real *) starts);
  set_volume_separations(volume, (VIO_Real *) steps);

  int n = voxel_count(volume);
  double *v = values(volume);
  for(int i = 0; i < n; i++) v[i] = 0.0;

  return volume;
}

VIO_Volume like(VIO_Volume model)
{
  int sizes[VIO_N_DIMENSIONS];
  VIO_Real starts[VIO_N_DIMENSIONS], steps[VIO_N_DIMENSIONS];
  get_volume_sizes(model, sizes);
  get_volume_starts(model, starts);
  get_volume_separations(model, steps);
  return create_like(model, sizes, starts, steps);
}

int voxel_count(VIO_Volume volume)
{
  int sizes[VIO_N_DIMENSIONS];
  get_volume_sizes(volume, sizes);
  return sizes[0] * sizes[1] * sizes[2];
}

/* volume_io reaches its data through arrays of pointers, so a flat view is an
 * assumption about how alloc_multidim_array laid the block out.  Check it
 * rather than trust it: a volume that is not contiguous would give every bulk
 * loop in the pipeline silently wrong answers. */
double *values(VIO_Volume volume)
{
  VIO_BOOL signed_flag;
  if(get_volume_nc_data_type(volume, &signed_flag) != NC_DOUBLE)
    fail("values() on a volume that is not NC_DOUBLE");

  int sizes[VIO_N_DIMENSIONS];
  get_volume_sizes(volume, sizes);

  /* The macro yields void *, so the flat view is spelled out here rather than
   * assigned through a typed pointer. */
  void *raw;
  GET_VOXEL_PTR_3D(raw, volume, 0, 0, 0);
  double *base = (double *) raw;

  if(sizes[2] > 1)
    {
      GET_VOXEL_PTR_3D(raw, volume, 0, 0, 1);
      if((double *) raw != base + 1) fail("volume data is not contiguous along z");
    }
  if(sizes[1] > 1)
    {
      GET_VOXEL_PTR_3D(raw, volume, 0, 1, 0);
      if((double *) raw != base + sizes[2])
        fail("volume data is not contiguous along y");
    }
  if(sizes[0] > 1)
    {
      GET_VOXEL_PTR_3D(raw, volume, 1, 0, 0);
      if((double *) raw != base + sizes[1] * sizes[2])
        fail("volume data is not contiguous along x");
    }

  /* An NC_DOUBLE volume carries no voxel-to-real scaling, so the flat buffer
   * holds real values.  If that ever stopped being true every bulk loop would
   * read voxel values and call them intensities. */
  double first = get_volume_real_value(volume, 0, 0, 0, 0, 0);
  if(first != base[0]) fail("NC_DOUBLE volume has a voxel-to-real scaling");

  return base;
}

VIO_Volume shrink(VIO_Volume in, double factor)
{
  int sizes[VIO_N_DIMENSIONS], out_sizes[VIO_N_DIMENSIONS];
  VIO_Real starts[VIO_N_DIMENSIONS], steps[VIO_N_DIMENSIONS];
  VIO_Real out_steps[VIO_N_DIMENSIONS];
  double stride[VIO_N_DIMENSIONS];

  get_volume_sizes(in, sizes);
  get_volume_starts(in, starts);
  get_volume_separations(in, steps);

  /* The finest sampling sets the target; coarser axes are left alone
   * (nu_estimate_np_and_em.in:962-973). */
  double finest = fabs(steps[0]);
  for(int i = 1; i < VIO_N_DIMENSIONS; i++)
    if(fabs(steps[i]) < finest) finest = fabs(steps[i]);
  double newstep = finest * factor;

  for(int i = 0; i < VIO_N_DIMENSIONS; i++)
    {
      if(fabs(steps[i]) < newstep)
        {
          out_steps[i] = steps[i] * factor;
          out_sizes[i] = (int) ceil((sizes[i] - 1) / factor) + 1;
          stride[i] = factor;
        }
      else
        {
          out_steps[i] = steps[i];
          out_sizes[i] = sizes[i];
          stride[i] = 1.0;
        }
    }

  /* start is not passed to mincresample, so it is kept, which is what makes
   * output voxel i sit on input voxel i*factor and lets a spline fitted here
   * be evaluated on the full-resolution grid unchanged. */
  VIO_Volume out = create_like(in, out_sizes, starts, out_steps);

  double *dst = values(out);
  double *src = values(in);

  for(int i = 0; i < out_sizes[0]; i++)
    {
      int si = (int) floor(i * stride[0] + 0.5);
      for(int j = 0; j < out_sizes[1]; j++)
        {
          int sj = (int) floor(j * stride[1] + 0.5);
          for(int k = 0; k < out_sizes[2]; k++)
            {
              int sk = (int) floor(k * stride[2] + 0.5);
              double value = 0.0;
              /* ceil((n-1)/factor)+1 can ask for a sample past the last
               * voxel; mincresample fills those with zero. */
              if(si < sizes[0] && sj < sizes[1] && sk < sizes[2])
                value = src[(si * sizes[1] + sj) * sizes[2] + sk];
              dst[(i * out_sizes[1] + j) * out_sizes[2] + k] = value;
            }
        }
    }

  return out;
}

VIO_Volume resample_label(VIO_Volume in, VIO_Volume model)
{
  int in_sizes[VIO_N_DIMENSIONS], out_sizes[VIO_N_DIMENSIONS];
  get_volume_sizes(in, in_sizes);
  get_volume_sizes(model, out_sizes);

  VIO_Volume out = like(model);
  double *dst = values(out);
  double *src = values(in);

  for(int i = 0; i < out_sizes[0]; i++)
    for(int j = 0; j < out_sizes[1]; j++)
      for(int k = 0; k < out_sizes[2]; k++)
        {
          /* Through world coordinates rather than index arithmetic: the two
           * grids share an origin in this pipeline, but nothing in
           * resample_labels requires that. */
          VIO_Real x, y, z, u, v, w;
          convert_3D_voxel_to_world(model, (VIO_Real) i, (VIO_Real) j,
                                    (VIO_Real) k, &x, &y, &z);
          convert_3D_world_to_voxel(in, x, y, z, &u, &v, &w);

          int u0 = (int) floor(u), v0 = (int) floor(v), w0 = (int) floor(w);
          double du = u - u0, dv = v - v0, dw = w - w0;

          double sum = 0.0;
          for(int a = 0; a < 2; a++)
            for(int b = 0; b < 2; b++)
              for(int c = 0; c < 2; c++)
                {
                  double weight = (a ? du : 1.0 - du) * (b ? dv : 1.0 - dv)
                                  * (c ? dw : 1.0 - dw);
                  if(weight == 0.0) continue;
                  int ii = u0 + a, jj = v0 + b, kk = w0 + c;
                  /* Outside the input contributes nothing, which is the fill
                   * value mincresample uses. */
                  if(ii < 0 || jj < 0 || kk < 0) continue;
                  if(ii >= in_sizes[0] || jj >= in_sizes[1] || kk >= in_sizes[2])
                    continue;
                  sum += weight * src[(ii * in_sizes[1] + jj) * in_sizes[2] + kk];
                }

          dst[(i * out_sizes[1] + j) * out_sizes[2] + k] = sum >= 0.5 ? 1.0 : 0.0;
        }

  return out;
}

Stats masked_stats(VIO_Volume volume, VIO_Volume mask)
{
  int n = voxel_count(volume);
  double *v = values(volume);
  double *m = mask ? values(mask) : NULL;

  if(mask && voxel_count(mask) != n)
    fail("masked_stats: mask and volume are on different grids");

  Stats s;
  s.count = 0;
  s.mean = s.stddev = 0.0;
  s.minimum = s.maximum = 0.0;

  double sum = 0.0, sum2 = 0.0;
  for(int i = 0; i < n; i++)
    {
      if(m && m[i] == 0.0) continue;   /* volumeStats.cc:236 */
      if(s.count == 0) { s.minimum = s.maximum = v[i]; }
      if(v[i] < s.minimum) s.minimum = v[i];
      if(v[i] > s.maximum) s.maximum = v[i];
      sum += v[i];
      sum2 += v[i] * v[i];
      s.count++;
    }

  if(s.count > 0)
    {
      s.mean = sum / s.count;
      /* Population, as volumeStats.cc:276 computes it. */
      double variance = sum2 / s.count - s.mean * s.mean;
      s.stddev = variance > 0.0 ? sqrt(variance) : 0.0;
    }

  return s;
}

nc_type storage_type(const std::string &path, VIO_BOOL *signed_flag)
{
  VIO_Volume volume;
  volume_input_struct input_info;

  if(start_volume_input((char *) path.c_str(), VIO_N_DIMENSIONS,
                        File_order_dimension_names, NC_UNSPECIFIED, FALSE,
                        0.0, 0.0, TRUE, &volume,
                        (minc_input_options *) NULL, &input_info) != VIO_OK)
    fail("cannot read volume header", path);

  nc_type type = get_volume_nc_data_type(volume, signed_flag);

  delete_volume_input(&input_info);
  delete_volume(volume);

  return type;
}

void voxel_range(const std::string &path, double *lo, double *hi)
{
  VIO_Volume volume;
  volume_input_struct input_info;

  if(start_volume_input((char *) path.c_str(), VIO_N_DIMENSIONS,
                        File_order_dimension_names, NC_UNSPECIFIED, FALSE,
                        0.0, 0.0, TRUE, &volume,
                        (minc_input_options *) NULL, &input_info) != VIO_OK)
    fail("cannot read volume header", path);

  get_volume_voxel_range(volume, lo, hi);

  delete_volume_input(&input_info);
  delete_volume(volume);
}

void save(VIO_Volume volume, const std::string &path,
          const std::string &like_path, nc_type type, VIO_BOOL signed_flag,
          const std::string &history)
{
  /* output_modified_volume is what `mincmath -copy_header` amounts to: the
   * header comes from like_path and the data from this buffer.  Passing 0,0
   * for the range lets volume_io take it from the data. */
  if(output_modified_volume((char *) path.c_str(), type, signed_flag, 0.0, 0.0,
                            volume, (char *) like_path.c_str(),
                            (char *) history.c_str(),
                            (minc_output_options *) NULL) != VIO_OK)
    fail("cannot write volume", path);
}

}  // namespace n3
