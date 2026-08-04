/* nu_estimate_np_and_em's iteration, in memory.
 *
 * The Perl driver's loop (nu_estimate_np_and_em.in:101-205) with every
 * intermediate held as a buffer instead of written to a MINC file: log
 * transform, mask, then alternately sharpen the corrected volume and fit a
 * spline to what is left over.
 */

#ifndef N3_NUESTIMATE_H
#define N3_NUESTIMATE_H

#include <string>
#include <vector>

#include "FitField.h"
#include "Histogram.h"

namespace n3 {

struct EstimateOptions
{
  double shrink;              /* -shrink, 1 for none */
  double distance;            /* -distance, mm */
  double lambda;              /* -lambda */
  int subsample;              /* -spline_subsample */
  enum spline_type type;      /* -b_spline or -tp_spline */

  int bins;                   /* -bins */
  double fwhm, noise;         /* -sharpen */
  bool window;                /* -parzen */
  double parzen_sigma;        /* -parzen_sigma, 0 for off */
  bool blur;                  /* -blur: skip the deconvolution */

  double background_threshold;  /* -background */
  bool bimodalT;                /* -bimodalT */
  int bimodal_bins;             /* ceil(voxelMax-voxelMin+1) of the input file,
                                 * which is what volume_stats would have used */

  std::vector<int> iterations;   /* -iterations */
  std::vector<double> stop;      /* -stop */

  bool normalize_field;       /* -normalize_field */
  bool verbose;

  /* Round the histogram, its range and the lookup table to the six decimals
   * the Perl passes them between programs with.  A verification instrument,
   * not a feature: with it on, what is left against the Perl is the volume
   * quantisation. */
  bool legacy_rounding;

  EstimateOptions();
};

/* Intermediates of one iteration, for comparison against the driver's own
 * -save_fields output.  Both are exponentiated, as the driver writes them.
 * The caller owns whatever is filled in. */
struct EstimateTrace
{
  int iteration;               /* which iteration to capture, zero based */
  VIO_Volume before_smoothing; /* the driver's <base>_est<iter>.mnc */
  VIO_Volume field;            /* the driver's <base>_field<iter>.mnc */

  EstimateTrace() : iteration(0), before_smoothing(NULL), field(NULL) {}
};

/* The estimated field, as compact_spline_volume fits it: a spline through
 * exp(residue), not through the log field the loop accumulates.
 *
 * iterations_run and final_change report what the stopping rule did, which is
 * the first thing to check when two implementations disagree end to end.  The
 * caller owns the result.
 *
 * When mapping_path is set the compact fit is also written there (the .imp),
 * like Perl's compact_spline_volume.  It is written from the estimation grid,
 * whose separations are the ones Field::domain is expressed in, so giving any
 * other volume to the header would put the domain in the wrong world place --
 * the .imp's Domain block is world coordinates (fieldIO.cc:120-133). */
Field *nu_estimate(VIO_Volume input, VIO_Volume user_mask,
                   const EstimateOptions &options,
                   int *iterations_run, double *final_change,
                   EstimateTrace *trace = NULL,
                   const std::string *mapping_path = NULL);

}  // namespace n3

#endif
