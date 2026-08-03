/* Cycle 6: applying the intensity mapping.
 *
 * minclookup -continuous is how the sharpened histogram becomes a sharpened
 * volume.  It interpolates linearly between the two table entries bracketing
 * each voxel and clamps outside the range; nearest-neighbour is the wrong
 * answer and looks almost right.
 *
 * The oracle is written as double, so no storage type enters the comparison,
 * and both sides read the same lookup table text.  What is left is the
 * arithmetic of a linear interpolation, so the two agree to double rounding.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"
#include "../../src/N3Pipeline/MincTools.h"

#include <cmath>
#include <string>
#include <vector>

int main()
{
  std::string data = std::string(N3_DATA_DIR);

  double lo, hi;
  n3fixture::read_domain("hist_window.txt", &lo, &hi);

  std::vector<double> table = n3fixture::read_text("sharp_window.txt");
  int entries = (int) table.size() / 2;
  std::vector<double> lut(entries);
  for(int i = 0; i < entries; i++) lut[i] = table[2 * i + 1];

  VIO_Volume chunk = n3::load(data + "/chunk.mnc.gz");
  int n = n3::voxel_count(chunk);

  /* Properties first, on a table small enough to reason about.  A two-entry
   * table over [0, 1] is the identity on its own range and constant outside
   * it, whatever the implementation. */
  {
    int sizes[VIO_N_DIMENSIONS] = { 5, 1, 1 };
    VIO_Real starts[VIO_N_DIMENSIONS] = { 0.0, 0.0, 0.0 };
    VIO_Real steps[VIO_N_DIMENSIONS] = { 1.0, 1.0, 1.0 };
    VIO_Volume tiny = n3::create_like(chunk, sizes, starts, steps);
    double *v = n3::values(tiny);
    v[0] = -1.0; v[1] = 0.0; v[2] = 0.25; v[3] = 1.0; v[4] = 2.0;

    std::vector<double> ramp(2);
    ramp[0] = 10.0; ramp[1] = 20.0;
    n3::apply_lookup(tiny, ramp, 0.0, 1.0);

    CHECK_NEAR("below the range clamps to the first entry", v[0], 10.0, 0.0);
    CHECK_NEAR("the low end is the first entry", v[1], 10.0, 0.0);
    CHECK_NEAR("between entries it interpolates, not rounds", v[2], 12.5, 0.0);
    CHECK_NEAR("the high end is the last entry", v[3], 20.0, 0.0);
    CHECK_NEAR("above the range clamps to the last entry", v[4], 20.0, 0.0);

    delete_volume(tiny);
  }

  /* The whole volume against minclookup.
   *
   * minclookup reads the table's first column, which sharpen_hist wrote with
   * six decimals: 0.005025 where the exact position is 0.005025125628.  Given
   * those same positions the two implementations agree to double rounding.
   * Given the exact ones they do not, and the gap is measured below rather
   * than absorbed into a bound. */
  {
    std::vector<double> oracle = n3fixture::read_f64("lookup_applied.f64");
    n3fixture::must((int) oracle.size() == n, "size mismatch against oracle");

    std::vector<double> positions(entries);
    for(int i = 0; i < entries; i++)
      positions[i] = lo + table[2 * i] * (hi - lo);

    n3::apply_lookup(chunk, positions, lut);
    double *mine = n3::values(chunk);
    CHECK_RMS("the whole volume is minclookup's", mine, &oracle[0], n, 1e-14);

    double worst = 0.0, span = 0.0;
    for(int i = 0; i < n; i++)
      { double d = fabs(mine[i] - oracle[i]); if(d > worst) worst = d; }
    for(int i = 0; i < entries; i++)
      { double d = fabs(lut[i]); if(d > span) span = d; }
    printf("  (worst voxel %.3e on entries of order %.3e)\n", worst, span);

    /* What the six decimals cost.  This is the in-memory pipeline's own
     * answer, and it is deliberately not the oracle's. */
    VIO_Volume exact = n3::load(data + "/chunk.mnc.gz");
    n3::apply_lookup(exact, lut, lo, hi);
    double drift = n3check::rel_rms(n3::values(exact), &oracle[0], n);
    CHECK_TRUE("exact entry positions move the result, and by more than "
               "rounding", drift > 1e-9);
    printf("  (six decimals of entry position are worth %.3e relative)\n", drift);
    delete_volume(exact);

    /* Nearest neighbour instead of interpolation would be a plausible port and
     * is measurably wrong; establish the size of the difference so that a
     * regression to it cannot hide inside the bound above. */
    VIO_Volume again = n3::load(data + "/chunk.mnc.gz");
    double *nn = n3::values(again);
    double step = (hi - lo) / (entries - 1);
    for(int i = 0; i < n; i++)
      {
        double t = (nn[i] - lo) / step;
        int k = (int) floor(t + 0.5);
        if(k < 0) k = 0;
        if(k > entries - 1) k = entries - 1;
        nn[i] = lut[k];
      }
    double nn_rms = n3check::rel_rms(nn, &oracle[0], n);
    CHECK_TRUE("nearest neighbour is not interpolation", nn_rms > 1e-6);
    printf("  (nearest neighbour would sit %.3e away)\n", nn_rms);
    delete_volume(again);
  }

  delete_volume(chunk);
  return n3check::report("lookup");
}
