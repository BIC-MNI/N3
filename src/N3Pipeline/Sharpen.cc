#include "Sharpen.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <EBTKS/Matrix.h>

/* The three free functions sharpen_hist.cc does export.  Declared here rather
 * than in a header because the original has none: the file is a program. */
DblMat gaussian(double fwhm, int size);
CompMat weiner(CompMat &blur, double noise);
void non_negative(DblMat *X);

namespace n3 {

std::vector<double> sharpen_lookup(const std::vector<double> &counts,
                                   double min_bin, double max_bin,
                                   double fwhm, double noise, bool blur_flag)
{
  int bins = (int) counts.size();
  if(bins < 2 || min_bin == max_bin)
    {
      fprintf(stderr, "n3::sharpen_lookup: degenerate histogram\n");
      exit(1);
    }

  DblMat X(bins, 1, 0.0);
  for(int i = 0; i < bins; i++) X(i, 0) = counts[i];

  /* 2^(ceil(log2(bins)) + 1): at 200 bins that is 512, not 256, and the
   * histogram sits centred at offset rather than at index 0. */
  int padded_size = int(pow(2, ceil(log((double) bins) / log(2.0)) + 1) + .5);
  int offset = (padded_size - bins) / 2;

  /* The kernel width is in bin units, not intensity units. */
  double slope = (max_bin - min_bin) / double(bins - 1);
  CompMat blur = fft(gaussian(fwhm / slope, padded_size), 0, 1);
  CompMat filter = weiner(blur, noise);

  DblMat X_padded(padded_size, 1, 0.0);
  X_padded.insert(X, offset, 0);

  DblMat f;
  if(!blur_flag)
    {
      f = real(ifft(pmultEquals(asCompMat(X_padded).fft(0, 1), filter), 0, 1));
      non_negative(&f);
    }
  else
    f = X_padded;

  DblMat moment(padded_size, 1);
  for(int i = 0; i < padded_size; i++)
    moment(i, 0) = (min_bin + (i - offset) * slope) * f(i, 0);

  /* E[u | v] under the Gaussian kernel. */
  DblMat Y_padded =
    pdiv(real(ifft(pmultEquals(asCompMat(moment).fft(0, 1), blur), 0, 1)),
         real(ifft(pmultEquals(asCompMat(f).fft(0, 1), blur), 0, 1)));

  std::vector<double> Y(bins);
  for(int i = 0; i < bins; i++)
    {
      double value = Y_padded(offset + i, 0);
      /* Wherever the deconvolved density vanishes the ratio is not a number;
       * sharpen_hist.cc:173-177 replaces those with zero. */
      Y[i] = std::isfinite(value) ? value : 0.0;
    }

  return Y;
}

}  // namespace n3
