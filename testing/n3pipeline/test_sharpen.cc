/* Cycle 5: the histogram deconvolution.
 *
 * sharpen_hist's input here is the same text the test reads, so what is
 * compared is the deconvolution alone.  The six decimals that text carries are
 * a separate divergence and are measured where the pipeline is assembled, not
 * folded into this bound.
 *
 * The oracle prints the mapping with %lf, so half a unit in the sixth decimal
 * is what it pins.  On intensities of order 1e5 that is a relative 5e-12,
 * which is a demanding comparison for a chain of four transforms and is the
 * reason this cycle is worth having separately from the pipeline.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Sharpen.h"

#include <cmath>
#include <string>
#include <vector>

static const double PRINTED = 0.5e-6;

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

int main()
{
  double min_bin, max_bin;
  n3fixture::read_domain("hist_window.txt", &min_bin, &max_bin);

  std::vector<double> table = n3fixture::read_text("hist_window.txt");
  int bins = (int) table.size() / 2;
  std::vector<double> counts(bins), centres(bins);
  for(int i = 0; i < bins; i++)
    { centres[i] = table[2 * i]; counts[i] = table[2 * i + 1]; }

  {
    std::vector<double> mine =
      n3::sharpen_lookup(counts, min_bin, max_bin, 0.15, 0.01, false);
    std::vector<double> oracle = second_column("sharp_window.txt", bins);
    CHECK_TRUE("the mapping has one entry per bin", (int) mine.size() == bins);
    CHECK_ALL("the mapping is sharpen_hist's", &mine[0], &oracle[0], bins,
              PRINTED);
  }

  {
    std::vector<double> mine =
      n3::sharpen_lookup(counts, min_bin, max_bin, 0.15, 0.01, true);
    std::vector<double> oracle = second_column("sharp_window_blur.txt", bins);
    CHECK_ALL("-blur skips the deconvolution, as sharpen_hist's does",
              &mine[0], &oracle[0], bins, PRINTED);
  }

  /* Properties that hold whatever the oracle says.  The mapping is a
   * conditional expectation of intensity given intensity, so it must be
   * finite everywhere -- sharpen_hist replaces non-finite entries with zero
   * (sharpen_hist.cc:173-177) and a port that dropped that would produce a
   * lookup table full of NaN wherever the deconvolved density vanished. */
  {
    std::vector<double> mine =
      n3::sharpen_lookup(counts, min_bin, max_bin, 0.15, 0.01, false);
    bool finite = true;
    for(int i = 0; i < bins; i++) if(!std::isfinite(mine[i])) finite = false;
    CHECK_TRUE("every entry is finite", finite);

    /* An empty histogram has no conditional expectation anywhere, so the
     * whole mapping is the zero that replaces the non-finite entries. */
    std::vector<double> empty(bins, 0.0);
    std::vector<double> degenerate =
      n3::sharpen_lookup(empty, min_bin, max_bin, 0.15, 0.01, false);
    double worst = 0.0;
    for(int i = 0; i < bins; i++)
      if(fabs(degenerate[i]) > worst) worst = fabs(degenerate[i]);
    CHECK_NEAR("an empty histogram maps everything to zero", worst, 0.0, 0.0);
  }

  /* The deconvolution sharpens: the mapping must pull intensities together
   * rather than spread them, so its range is inside the histogram's.  This is
   * the property the whole method rests on, and it is checked here rather
   * than assumed because sharpen_hist itself only ever warned about it
   * (sharpen_hist.cc:179-183, commented out). */
  {
    std::vector<double> mine =
      n3::sharpen_lookup(counts, min_bin, max_bin, 0.15, 0.01, false);
    double lo = mine[0], hi = mine[0];
    for(int i = 1; i < bins; i++)
      { if(mine[i] < lo) lo = mine[i]; if(mine[i] > hi) hi = mine[i]; }
    CHECK_TRUE("the mapping's range sits inside the histogram's",
               lo >= min_bin && hi <= max_bin);
    printf("  (mapping spans %.1f, histogram %.1f)\n", hi - lo, max_bin - min_bin);
  }

  return n3check::report("sharpen");
}
