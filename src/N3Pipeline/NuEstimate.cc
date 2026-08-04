#include "NuEstimate.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "Buffers.h"
#include "MincTools.h"
#include "Sharpen.h"
#include "Stopping.h"

namespace n3 {

EstimateOptions::EstimateOptions()
  : shrink(4.0), distance(200.0), lambda(1e-7), subsample(1), type(b_spline),
    bins(200), fwhm(0.15), noise(0.01), window(true), parzen_sigma(0.0),
    blur(false), background_threshold(1.0), bimodalT(false),
    bimodal_bins(4096), normalize_field(false), verbose(false),
    legacy_rounding(false)
{
  iterations.push_back(50);
  stop.push_back(0.001);
}

/* What %lf writes and scanf reads back.  Only used under -legacy_rounding. */
static double six_decimals(double value)
{
  char text[64];
  sprintf(text, "%lf", value);
  return atof(text);
}

Field *nu_estimate(VIO_Volume input, VIO_Volume user_mask,
                   const EstimateOptions &options,
                   int *iterations_run, double *final_change,
                   EstimateTrace *trace, const std::string *mapping_path)
{
  /* 1. The estimation grid (nu_estimate_np_and_em.in:63-64). */
  VIO_Volume work = (options.shrink != 1.0) ? shrink(input, options.shrink)
                                            : like(input);
  if(options.shrink == 1.0)
    {
      double *dst = values(work), *src = values(input);
      for(int i = 0, n = voxel_count(work); i < n; i++) dst[i] = src[i];
    }

  VIO_Volume mask_wks = NULL;
  if(user_mask)
    mask_wks = (options.shrink != 1.0) ? resample_label(user_mask, work)
                                       : resample_label(user_mask, work);

  int n = voxel_count(work);

  /* 2. log(max(v, 1)): mincmath -clamp -const2 1 1.7e308 then -log (:657-664). */
  VIO_Volume log_volume = like(work);
  {
    double *dst = values(log_volume), *src = values(work);
    for(int i = 0; i < n; i++) dst[i] = log(src[i] < 1.0 ? 1.0 : src[i]);
  }

  /* 3. The mask: input above the background threshold, intersected with the
   * user's (CreateMask, :297).  -bimodalT takes the threshold from
   * volume_stats' rule rather than from -background. */
  double background = options.background_threshold;
  if(options.bimodalT)
    background = bimodal_threshold_volume_stats(work, user_mask ? mask_wks : NULL,
                                                options.bimodal_bins);

  VIO_Volume mask = like(work);
  {
    double *dst = values(mask), *src = values(work);
    double *user = mask_wks ? values(mask_wks) : NULL;
    int inside = 0;
    for(int i = 0; i < n; i++)
      {
        bool in = src[i] > background;
        if(user && user[i] == 0.0) in = false;
        dst[i] = in ? 1.0 : 0.0;
        if(in) inside++;
      }
    if(inside == 0)
      {
        fprintf(stderr, "n3::nu_estimate: the composite mask is empty\n");
        exit(1);
      }
  }

  /* 4. apply_mask (:86) -- the log volume is zero outside. */
  {
    double *v = values(log_volume), *m = values(mask);
    for(int i = 0; i < n; i++) if(m[i] == 0.0) v[i] = 0.0;
  }

  VIO_Volume residue = like(work);      /* zero: the flat initial field */
  VIO_Volume corrected = like(work);
  VIO_Volume estimate = like(work);
  VIO_Volume working = like(work);

  HistogramOptions histogram_options;
  histogram_options.bins = options.bins;
  histogram_options.window = options.window;
  histogram_options.parzen_sigma = options.parzen_sigma;

  int total = total_iterations(options.iterations);
  int iter = 0;
  double change = 0.0;

  for(iter = 0; iter < total; iter++)
    {
      /* corrected = log - residue */
      {
        double *c = values(corrected), *l = values(log_volume),
               *r = values(residue);
        for(int i = 0; i < n; i++) c[i] = l[i] - r[i];
      }

      /* sharpen_estimate (:500): histogram, deconvolve, look up, re-mask. */
      {
        DHistogram *h = histogram(corrected, mask, 1, histogram_options);

        double min_bin = h->binCenter(0);
        double max_bin = h->binCenter(h->nBins() - 1);
        std::vector<double> counts(h->nBins());
        for(int i = 0; i < h->nBins(); i++) counts[i] = (*h)[i];
        delete h;

        if(options.legacy_rounding)
          {
            min_bin = six_decimals(min_bin);
            max_bin = six_decimals(max_bin);
            for(size_t i = 0; i < counts.size(); i++)
              counts[i] = six_decimals(counts[i]);
          }

        std::vector<double> lut = sharpen_lookup(counts, min_bin, max_bin,
                                                 options.fwhm, options.noise,
                                                 options.blur);

        double *e = values(estimate), *c = values(corrected);
        for(int i = 0; i < n; i++) e[i] = c[i];

        if(options.legacy_rounding)
          {
            /* sharpen_hist writes the entry positions with six decimals too,
             * and minclookup interpolates against those. */
            std::vector<double> positions(lut.size());
            double step = 1.0 / (lut.size() - 1);
            for(size_t i = 0; i < lut.size(); i++)
              {
                positions[i] = min_bin
                  + six_decimals(i * step) * (max_bin - min_bin);
                lut[i] = six_decimals(lut[i]);
              }
            apply_lookup(estimate, positions, lut);
          }
        else
          apply_lookup(estimate, lut, min_bin, max_bin);

        double *m = values(mask);
        for(int i = 0; i < n; i++) if(m[i] == 0.0) e[i] = 0.0;
      }

      /* working = log - estimate; the old residue is kept for the stopping
       * rule (:150-152). */
      {
        double *w = values(working), *l = values(log_volume),
               *e = values(estimate), *r = values(residue), *old = values(corrected);
        for(int i = 0; i < n; i++) { w[i] = l[i] - e[i]; old[i] = r[i]; }
      }

      if(trace && trace->iteration == iter)
        {
          /* exp(working), which is what -save_fields writes before the
           * smoothing (:154-155). */
          trace->before_smoothing = like(work);
          double *t = values(trace->before_smoothing), *w = values(working);
          for(int i = 0; i < n; i++) t[i] = exp(w[i]);
        }

      /* The new field estimate. */
      {
        Field *fit = fit_field(working, mask, options.type, options.distance,
                               options.lambda, options.subsample);
        evaluate_field(fit, residue, mask);
        delete fit;
      }

      /* field_CV (:702): the population standard deviation, inside the mask,
       * of the change in the field.  Misnamed -- it is not a coefficient of
       * variation. */
      {
        double *w = values(working), *old = values(corrected),
               *r = values(residue);
        for(int i = 0; i < n; i++) w[i] = old[i] - r[i];
        Stats s = masked_stats(working, mask);
        change = s.stddev;
      }

      if(options.verbose)
        printf("CV for change in field estimate at iteration %d: %g\n",
               iter, change);

      if(trace && trace->iteration == iter)
        {
          /* exp(residue), which is what -save_fields writes after it
           * (:193-194). */
          trace->field = like(work);
          double *t = values(trace->field), *r = values(residue);
          for(int i = 0; i < n; i++) t[i] = exp(r[i]);
        }

      if(should_stop(iter, change, options.iterations, options.stop)) break;
    }

  *iterations_run = (iter < total) ? iter + 1 : total;
  *final_change = change;

  /* 8. field = exp(residue), optionally normalised to mean 1 in the mask
   * (:199-200). */
  VIO_Volume field = like(work);
  {
    double *f = values(field), *r = values(residue);
    for(int i = 0; i < n; i++) f[i] = exp(r[i]);

    if(options.normalize_field)
      {
        Stats s = masked_stats(field, mask);
        if(s.mean != 0.0)
          for(int i = 0; i < n; i++) f[i] /= s.mean;
        else
          fprintf(stderr, "Final field volume could not be normalized.\n");
      }
  }

  /* 9. compact_spline_volume (:202): a second, fresh fit, to the exponentiated
   * field rather than to the log one the loop accumulated. */
  Field *compact = fit_field(field, mask, options.type, options.distance,
                             options.lambda, options.subsample);

  if(mapping_path)
    {
      /* compact_spline_volume (:202) writes the .imp with the estimation
       * grid's geometry as its header, which is what makes the world-domain
       * conversion (fieldIO.cc:120-133) come out right. */
      char command[256];
      snprintf(command, sizeof(command), "nu_correct_cxx");
      save_field(*mapping_path, compact, work, command);
    }

  delete_volume(field);
  delete_volume(working);
  delete_volume(estimate);
  delete_volume(corrected);
  delete_volume(residue);
  delete_volume(mask);
  delete_volume(log_volume);
  if(mask_wks) delete_volume(mask_wks);
  delete_volume(work);

  return compact;
}

}  // namespace n3
