/* Cycle 5: the histogram deconvolution.
 *
 * The comparison that matters is against hist_log.txt, because that is the
 * domain nu_estimate sharpens in: it takes log(max(v,1)) first
 * (nu_estimate_np_and_em.in:661-662, NuEstimate.cc:55-59), so a bin is 0.011
 * wide and -fwhm 0.15 spans 13.6 of them.  The Weiner filter is then a real
 * one, the deconvolved density never collapses to a few ulps, and the mapping
 * contracts the histogram's range by about 5%.
 *
 * The oracle prints the mapping with %lf, so half a unit in the sixth decimal
 * is what it pins.  On log intensities of order 12 the arithmetic itself is
 * good to 2.6e-11, four orders under that, so the bound is the oracle's
 * printing plus a little room for the arithmetic (PRINTED + ARITHMETIC
 * below) -- this comparison cannot fail on a rounding difference, only on a
 * real one.
 *
 * The second block runs the same code on hist_window.txt, a raw-intensity
 * histogram, where the same -fwhm 0.15 is 3.7e-05 of a 4021-unit bin.  That
 * is a genuine edge case rather than a second measurement: the kernel
 * underflows to a discrete delta and the whole transform becomes the
 * identity.  See there for why its oracle is not compared entry by entry.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Sharpen.h"

#include <cfloat>
#include <cmath>
#include <string>
#include <vector>

static const double PRINTED = 0.5e-6;

/* What the arithmetic adds on top of the oracle's printing.  The mapping
 * divides by the smoothed density, so an entry survives to
 * centre(i)*DBL_EPSILON*max(counts)/density(i), which on hist_log.txt is
 * 2.6e-11 at worst.  Without this term the bound would be exactly the printed
 * half-ulp, and an entry landing on a rounding boundary could round the other
 * way and fail on nothing.  1e-9 is forty times the arithmetic and a
 * millionth of the 1e-3 by which a broken deconvolution misses. */
static const double ARITHMETIC = 1e-9;

/* The oracle is two columns; the mapping is the second. */
static std::vector<double> second_column(const std::string &name, int expected)
{
  std::vector<double> table = n3fixture::read_text(name);
  n3fixture::must((int) table.size() == 2 * expected,
                  name + ": unexpected number of entries");
  std::vector<double> out(expected);
  for(int i = 0; i < expected; i++) out[i] = table[2 * i + 1];
  return out;
}

/* A histogram fixture is a "# domain:" line and two columns; the counts are
 * the second.  volume_hist writes both files the tests read this way. */
static std::vector<double> load_counts(const std::string &name,
                                       double *min_bin, double *max_bin)
{
  n3fixture::read_domain(name, min_bin, max_bin);
  std::vector<double> table = n3fixture::read_text(name);
  std::vector<double> counts(table.size() / 2);
  for(size_t i = 0; i < counts.size(); i++) counts[i] = table[2 * i + 1];
  return counts;
}

static double max_abs_difference(const std::vector<double> &a,
                                 const std::vector<double> &b)
{
  double worst = 0.0;
  for(size_t i = 0; i < a.size(); i++)
    { double d = fabs(a[i] - b[i]); if(d > worst) worst = d; }
  return worst;
}

int main()
{
  /* ---- the deconvolution, in the domain the pipeline uses --------------- */
  {
    double min_bin, max_bin;
    std::vector<double> counts = load_counts("hist_log.txt", &min_bin, &max_bin);
    int bins = (int) counts.size();

    std::vector<double> mine =
      n3::sharpen_lookup(counts, min_bin, max_bin, 0.15, 0.01, false);
    std::vector<double> oracle = second_column("sharp_log.txt", bins);
    CHECK_TRUE("the mapping has one entry per bin", (int) mine.size() == bins);
    CHECK_ALL("the mapping is sharpen_hist's", &mine[0], &oracle[0], bins,
              PRINTED + ARITHMETIC);

    /* How much of that margin is arithmetic rather than the oracle's own
     * printing.  Zero says the two agree in every digit the oracle carries.
     * Reported and not asserted: a value landing within 2.6e-11 of a rounding
     * boundary may round either way, which is not a defect. */
    int differs = 0;
    for(int i = 0; i < bins; i++)
      if(floor(mine[i] * 1e6 + 0.5) != floor(oracle[i] * 1e6 + 0.5)) differs++;
    printf("  (entries differing in the sixth decimal: %d of %d)\n",
           differs, bins);

    std::vector<double> blurred =
      n3::sharpen_lookup(counts, min_bin, max_bin, 0.15, 0.01, true);
    std::vector<double> oracle_blur = second_column("sharp_log_blur.txt", bins);
    CHECK_ALL("-blur skips the deconvolution, as sharpen_hist's does",
              &blurred[0], &oracle_blur[0], bins, PRINTED + ARITHMETIC);

    /* The property the method rests on, and the one check a port that had
     * dropped the deconvolution could not pass: the mapping is a conditional
     * expectation of intensity given intensity, so it must pull intensities
     * together rather than spread them.  sharpen_hist itself only ever warned
     * about this (sharpen_hist.cc:179-183, commented out). */
    double lo = mine[0], hi = mine[0];
    for(int i = 1; i < bins; i++)
      { if(mine[i] < lo) lo = mine[i]; if(mine[i] > hi) hi = mine[i]; }
    double contraction = (hi - lo) / (max_bin - min_bin);
    CHECK_TRUE("the deconvolution contracts the histogram's range",
               contraction < 1.0);
    printf("  (mapping spans %.6f of the histogram's %.6f, contraction %.4f)\n",
           hi - lo, max_bin - min_bin, contraction);

    /* -blur and -noise are inert on a delta kernel, so neither could be
     * measured before this fixture existed.  Both bounds are a hundred times
     * under the measured separation and far above any rounding. */
    double from_blur = max_abs_difference(mine, blurred);
    CHECK_TRUE("the deconvolution moves the mapping off the -blur answer",
               from_blur > 1e-3);

    std::vector<double> noisier =
      n3::sharpen_lookup(counts, min_bin, max_bin, 0.15, 0.5, false);
    double from_noise = max_abs_difference(mine, noisier);
    CHECK_TRUE("-noise reaches the Weiner filter", from_noise > 1e-3);
    printf("  (deconvolved vs -blur %.3e, noise 0.01 vs 0.5 %.3e)\n",
           from_blur, from_noise);

    bool finite = true;
    for(int i = 0; i < bins; i++) if(!std::isfinite(mine[i])) finite = false;
    CHECK_TRUE("every entry is finite", finite);

    /* An empty histogram has no conditional expectation anywhere, so the
     * whole mapping is the zero that replaces the non-finite entries
     * (sharpen_hist.cc:173-177).  A port that dropped that replacement would
     * return a lookup table full of NaN. */
    std::vector<double> empty(bins, 0.0);
    std::vector<double> degenerate =
      n3::sharpen_lookup(empty, min_bin, max_bin, 0.15, 0.01, false);
    double worst = 0.0;
    for(int i = 0; i < bins; i++)
      if(fabs(degenerate[i]) > worst) worst = fabs(degenerate[i]);
    CHECK_NEAR("an empty histogram maps everything to zero", worst, 0.0, 0.0);
  }

  /* ---- a kernel narrower than one bin -----------------------------------
   *
   * hist_window.txt is a raw-intensity histogram whose bins are 4021 units
   * wide, so -fwhm 0.15 is 3.7e-05 of a bin: gaussian() evaluates
   * scale*exp(-i*i*1.99e9) and underflows to a discrete delta, weiner()
   * becomes exactly its reciprocal, and the transform reduces to
   * centre(i)*X(i)/X(i).  The identity is what is asserted.
   *
   * sharp_window.txt is deliberately not compared entry by entry.  Three of
   * its bins carry a count of exactly zero, where the mapping is 0/0, and
   * sharpen_hist's own answers there are ratios of small ulp counts
   * (166*2^16/17, 51*2^16/5, 68*2^16/11) -- a fingerprint of one machine's
   * rounding rather than of the algorithm.  Cycle 6 still feeds the table to
   * minclookup, and test_lookup.cc still reads it.
   */
  {
    double min_bin, max_bin;
    std::vector<double> counts =
      load_counts("hist_window.txt", &min_bin, &max_bin);
    int bins = (int) counts.size();
    double slope = (max_bin - min_bin) / double(bins - 1);

    std::vector<double> mine =
      n3::sharpen_lookup(counts, min_bin, max_bin, 0.15, 0.01, false);

    double max_count = 0.0;
    for(int i = 0; i < bins; i++) if(counts[i] > max_count) max_count = counts[i];

    /* The identity is still reached through four FFTs and a division by the
     * smoothed density, so entry i survives only to
     *   centre(i) * DBL_EPSILON * max(counts) / counts(i),
     * which this fixture pushes to 2.8e-06 at its sparsest bin (count
     * 0.150275) -- five times the printed half-ulp, and the reason a constant
     * bound cannot serve here.  The bound below sits far above that and some
     * seven million times under the ~7e+03 a working deconvolution would move
     * these intensities, so it separates the two without pinning either
     * machine's rounding. */
    double worst = 0.0, predicted = 0.0;
    int compared = 0;
    for(int i = 0; i < bins; i++)
      {
        if(counts[i] <= 0.0) continue;
        double centre = min_bin + i * slope;
        double d = fabs(mine[i] - centre);
        double p = centre * DBL_EPSILON * max_count / counts[i];
        if(d > worst) worst = d;
        if(p > predicted) predicted = p;
        compared++;
      }
    CHECK_NEAR("a sub-bin kernel leaves the mapping the identity", worst, 0.0,
               1e-3);
    printf("  (%d of %d bins carry a count; roundoff there predicts %.3e)\n",
           compared, bins, predicted);

    bool finite = true;
    for(int i = 0; i < bins; i++) if(!std::isfinite(mine[i])) finite = false;
    CHECK_TRUE("every entry is finite where the density vanishes", finite);
  }

  return n3check::report("sharpen");
}
