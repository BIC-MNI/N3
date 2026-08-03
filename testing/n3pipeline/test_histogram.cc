/* Cycle 4: the masked histogram and its automatic range.
 *
 * The estimators are the legacy classes themselves, so what is under test is
 * the driver logic around them: which voxels are selected, what range they are
 * binned over, and that the option precedence matches volume_hist's.
 *
 * The oracle is volume_hist's text output, which carries six decimals of every
 * number (%lf, minchist.cc), so half a unit in the sixth is what it pins.
 * -gaussian_window exists only in this tree, so that one row is measured
 * against the locally built volume_hist and not the installed one.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"
#include "../../src/N3Pipeline/Histogram.h"

#include <string>
#include <vector>

/* Half a unit in the sixth decimal place, which is what %lf writes. */
static const double PRINTED = 0.5e-6;

static void compare(const char *label, DHistogram *mine, const char *fixture)
{
  char what[160];

  double min_bin, max_bin;
  n3fixture::read_domain(fixture, &min_bin, &max_bin);

  sprintf(what, "%s: the domain is volume_hist's", label);
  double margin = fabs(mine->binCenter(0) - min_bin);
  double upper = fabs(mine->binCenter(mine->nBins() - 1) - max_bin);
  if(upper > margin) margin = upper;
  n3check::record(what, margin <= PRINTED, margin, PRINTED);

  /* The text is two columns, bin centre then count, after the comment
   * lines. */
  std::vector<double> table = n3fixture::read_text(fixture);
  n3fixture::must((int) table.size() == 2 * mine->nBins(),
                  std::string(fixture) + ": unexpected number of entries");

  double centre_margin = 0.0, count_margin = 0.0, total = 0.0;
  for(int i = 0; i < mine->nBins(); i++)
    {
      double dc = fabs(mine->binCenter(i) - table[2 * i]);
      double dn = fabs((*mine)[i] - table[2 * i + 1]);
      if(dc > centre_margin) centre_margin = dc;
      if(dn > count_margin) count_margin = dn;
      total += (*mine)[i];
    }

  sprintf(what, "%s: %d bin centres", label, mine->nBins());
  n3check::record(what, centre_margin <= PRINTED, centre_margin, PRINTED);
  sprintf(what, "%s: %d counts", label, mine->nBins());
  n3check::record(what, count_margin <= PRINTED, count_margin, PRINTED);
  printf("  (%s holds %.1f samples)\n", label, total);
}

int main()
{
  std::string data = std::string(N3_DATA_DIR);
  VIO_Volume chunk = n3::load(data + "/chunk.mnc.gz");
  VIO_Volume mask  = n3::load(data + "/chunk_mask.mnc.gz");

  /* The range is taken over the selected voxels only, so it must be the
   * masked extrema and not the volume's.  Established before the oracle is
   * consulted, since a range that ignored the mask would still produce a
   * plausible-looking histogram. */
  {
    double lo, hi;
    n3::auto_range(chunk, mask, 1, &lo, &hi);
    n3::Stats inside = n3::masked_stats(chunk, mask);
    n3::Stats all = n3::masked_stats(chunk, NULL);
    CHECK_NEAR("the range is the masked minimum", lo, inside.minimum, 0.0);
    CHECK_NEAR("the range is the masked maximum", hi, inside.maximum, 0.0);
    CHECK_TRUE("which is not the volume's range", lo > all.minimum);
  }

  n3::HistogramOptions plain;
  n3::HistogramOptions window;   window.window = true;
  n3::HistogramOptions gauss;    gauss.parzen_sigma = 2.0;

  DHistogram *h;

  h = n3::histogram(chunk, mask, 1, plain);
  compare("plain", h, "hist_plain.txt");
  delete h;

  h = n3::histogram(chunk, mask, 1, window);
  compare("linear split", h, "hist_window.txt");
  delete h;

  h = n3::histogram(chunk, mask, 1, gauss);
  compare("gaussian sigma 2", h, "hist_gauss2.txt");
  delete h;

  /* A stated width wins over -window, which is how sharpen_volume keeps the
   * two from ever reaching volume_hist together (sharpen_volume.in:162-173). */
  {
    n3::HistogramOptions both;
    both.window = true;
    both.parzen_sigma = 2.0;
    DHistogram *a = n3::histogram(chunk, mask, 1, both);
    DHistogram *b = n3::histogram(chunk, mask, 1, gauss);
    double margin = 0.0;
    for(int i = 0; i < a->nBins(); i++)
      { double d = fabs((*a)[i] - (*b)[i]); if(d > margin) margin = d; }
    CHECK_NEAR("a stated width wins over the linear split", margin, 0.0, 0.0);
    delete a;
    delete b;
  }

  /* Every selected voxel lands somewhere, for the two windowed estimators.
   * The unwindowed one discards samples outside the first and last half bin,
   * so it is not held to this. */
  {
    n3::Stats inside = n3::masked_stats(chunk, mask);
    DHistogram *a = n3::histogram(chunk, mask, 1, window);
    DHistogram *b = n3::histogram(chunk, mask, 1, gauss);
    double sa = 0.0, sb = 0.0;
    for(int i = 0; i < a->nBins(); i++) { sa += (*a)[i]; sb += (*b)[i]; }
    CHECK_NEAR("the linear split conserves the sample count",
               sa, inside.count, 1e-6 * inside.count);
    CHECK_NEAR("the gaussian window conserves the sample count",
               sb, inside.count, 1e-6 * inside.count);
    delete a;
    delete b;
  }

  delete_volume(mask);
  delete_volume(chunk);
  return n3check::report("histogram");
}
