/* Cycle 7: the two bimodal thresholds.
 *
 * The pipeline needs both and they are not the same rule.  CreateMask uses
 * volume_stats -biModalT, an EBTKS histogram sized from the *file's* voxel
 * range; nu_evaluate uses mincstats -biModalT, 2000-bin Otsu reporting a bin
 * centre.  On chunk.mnc they differ by 33, so a port that used either for both
 * would still produce a plausible mask.  That they disagree is asserted here
 * alongside each matching its own oracle.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"
#include "../../src/N3Pipeline/MincTools.h"

#include <cmath>
#include <string>

/* volume_stats prints through cout at six significant digits; mincstats prints
 * ten.  Half a unit in the last of them is what each oracle pins. */
static double printed_precision(double value, int digits)
{
  if(value == 0.0) return 0.5 * pow(10.0, -(digits - 1));
  return 0.5 * pow(10.0, floor(log10(fabs(value))) - (digits - 1));
}

int main()
{
  std::string data = std::string(N3_DATA_DIR);
  std::string chunk_path = data + "/chunk.mnc";

  VIO_Volume chunk = n3::load(chunk_path);
  VIO_Volume mask  = n3::load(data + "/chunk_mask.mnc");

  /* volume_stats sizes its histogram from the file's valid_range, not from the
   * intensities, so the bin count is a property of how the volume was stored.
   * Check that the number read back is the one the file carries before using
   * it. */
  double voxel_lo, voxel_hi;
  n3::voxel_range(chunk_path, &voxel_lo, &voxel_hi);
  std::vector<double> recorded = n3fixture::read_text("chunk_valid_range.txt");
  n3fixture::must(recorded.size() == 2, "chunk_valid_range.txt");
  CHECK_NEAR("the voxel range is the file's valid_range (low)",
             voxel_lo, recorded[0], 0.0);
  CHECK_NEAR("the voxel range is the file's valid_range (high)",
             voxel_hi, recorded[1], 0.0);

  int bins = (int) ceil(voxel_hi - voxel_lo + 1.0);
  printf("  (%d bins, from a valid_range of %g to %g)\n", bins, voxel_lo, voxel_hi);

  double whole = n3::bimodal_threshold_volume_stats(chunk, NULL, bins);
  double oracle_whole = n3fixture::read_scalar("bimodal_volume_stats.txt");
  CHECK_NEAR("volume_stats' threshold, whole volume", whole, oracle_whole,
             printed_precision(oracle_whole, 6));

  double inside = n3::bimodal_threshold_volume_stats(chunk, mask, bins);
  double oracle_inside =
    n3fixture::read_scalar("bimodal_volume_stats_masked.txt");
  CHECK_NEAR("volume_stats' threshold, in mask", inside, oracle_inside,
             printed_precision(oracle_inside, 6));

  double otsu = n3::bimodal_threshold_mincstats(chunk);
  double oracle_otsu = n3fixture::read_scalar("bimodal_mincstats.txt");
  CHECK_NEAR("mincstats' threshold", otsu, oracle_otsu,
             printed_precision(oracle_otsu, 10));

  /* The two rules are different rules.  Recorded so that a later
   * simplification cannot quietly use one for both. */
  CHECK_TRUE("the two rules disagree on the same volume",
             fabs(whole - otsu) > printed_precision(oracle_whole, 6));
  printf("  (volume_stats %.4f, mincstats %.4f, apart by %.4f)\n",
         whole, otsu, fabs(whole - otsu));

  /* Both are thresholds on intensity, so both must lie inside the data. */
  n3::Stats all = n3::masked_stats(chunk, NULL);
  CHECK_TRUE("both thresholds lie within the volume's range",
             whole > all.minimum && whole < all.maximum
             && otsu > all.minimum && otsu < all.maximum);

  delete_volume(mask);
  delete_volume(chunk);
  return n3check::report("bimodal");
}
