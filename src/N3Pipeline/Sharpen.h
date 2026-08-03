/* Histogram deconvolution: the body of sharpen_hist's main(), over buffers.
 *
 * This is the one numerical stage that could not be reused by linking.
 * gaussian(), weiner() and non_negative() are free functions and are called
 * here, but the sequence around them -- padding, the offset, the two forward
 * transforms, the moment vector and the Nadaraya-Watson ratio -- lives inside
 * main() (sharpen_hist.cc:98-192) and is therefore transcribed.
 */

#ifndef N3_SHARPEN_H
#define N3_SHARPEN_H

#include <config.h>

#include <vector>

namespace n3 {

/* The intensity mapping, one value per histogram bin, over bin centres
 * min_bin .. max_bin.  With blur set the deconvolution is skipped, which is
 * sharpen_hist's -blur.
 *
 * The lookup table sharpen_hist writes is this vector on a domain normalised
 * to [0, 1]; sharpen_volume then hands minclookup the same range it took the
 * histogram over, so the normalisation cancels and is not reproduced here. */
std::vector<double> sharpen_lookup(const std::vector<double> &counts,
                                   double min_bin, double max_bin,
                                   double fwhm, double noise, bool blur);

}  // namespace n3

#endif
