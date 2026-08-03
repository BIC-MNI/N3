/* Cycle 3: the buffers themselves, and the statistics volume_stats prints.
 *
 * Two of these checks exist because of traps rather than because of an oracle:
 *
 *   - an NC_DOUBLE volume created here must carry no voxel-to-real scaling,
 *     or every bulk loop in the pipeline reads voxel values and calls them
 *     intensities;
 *   - the standard deviation must be the population one.  volume_stats
 *     computes E[x^2] - E[x]^2 (volumeStats.cc:276), the stopping rule needs
 *     that quantity, and the sample deviation is what most libraries give by
 *     default.  Asserted against a hand-computed value rather than an oracle,
 *     since an oracle at n = 236600 cannot tell the two apart.
 *
 * The four statistics are held to six significant digits, which is what
 * volume_stats' cout prints.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"

#include <cmath>
#include <string>
#include <vector>

/* volume_stats prints through cout, whose default precision is six
 * significant digits.  What the oracle therefore pins is the value to within
 * half a unit in the sixth of them -- an absolute quantity that grows with the
 * number, not a relative one. */
static double printed_precision(double value)
{
  if(value == 0.0) return 0.5e-5;
  return 0.5 * pow(10.0, floor(log10(fabs(value))) - 5.0);
}

static void check_stats(const char *label, const n3::Stats &s, const char *tag)
{
  static const char *names[] = { "mean", "stddev", "min", "max" };
  double mine[4] = { s.mean, s.stddev, s.minimum, s.maximum };

  for(int i = 0; i < 4; i++)
    {
      char what[160];
      double oracle = n3fixture::read_scalar(std::string("stats_") + tag + "_"
                                             + names[i] + ".txt");
      sprintf(what, "%s %s is volume_stats'", label, names[i]);
      CHECK_NEAR(what, mine[i], oracle, printed_precision(oracle));
    }
}

int main()
{
  std::string data = std::string(N3_DATA_DIR);
  VIO_Volume chunk = n3::load(data + "/chunk.mnc.gz");
  VIO_Volume mask  = n3::load(data + "/chunk_mask.mnc.gz");

  /* ------------------------------------------------ the buffer itself */
  {
    VIO_Volume empty = n3::like(chunk);
    int n = n3::voxel_count(empty);
    double *v = n3::values(empty);

    bool zeroed = true;
    for(int i = 0; i < n; i++) if(v[i] != 0.0) zeroed = false;
    CHECK_TRUE("a new buffer is zero", zeroed);

    int sizes[VIO_N_DIMENSIONS], model_sizes[VIO_N_DIMENSIONS];
    get_volume_sizes(empty, sizes);
    get_volume_sizes(chunk, model_sizes);
    CHECK_TRUE("a new buffer has the model's grid",
               sizes[0] == model_sizes[0] && sizes[1] == model_sizes[1]
               && sizes[2] == model_sizes[2]);

    /* The identity that lets the pipeline touch the flat buffer at all. */
    v[7] = 12345.678;
    CHECK_NEAR("the flat buffer holds real values, not voxel values",
               get_volume_real_value(empty, 0, 0, 7, 0, 0), 12345.678, 0.0);

    delete_volume(empty);
  }

  /* ------------------------------------------- the population deviation */
  {
    int sizes[VIO_N_DIMENSIONS] = { 2, 1, 2 };
    VIO_Real starts[VIO_N_DIMENSIONS] = { 0.0, 0.0, 0.0 };
    VIO_Real steps[VIO_N_DIMENSIONS] = { 1.0, 1.0, 1.0 };
    VIO_Volume tiny = n3::create_like(chunk, sizes, starts, steps);
    double *v = n3::values(tiny);
    v[0] = 1.0; v[1] = 2.0; v[2] = 3.0; v[3] = 6.0;

    n3::Stats s = n3::masked_stats(tiny, NULL);
    /* mean 3, E[x^2] = 12.5, population sd = sqrt(3.5); the sample deviation
     * of the same four numbers is sqrt(14/3) = 2.16. */
    CHECK_NEAR("stddev is the population one", s.stddev, sqrt(3.5), 1e-12);
    CHECK_TRUE("and is not the sample one",
               fabs(s.stddev - sqrt(14.0 / 3.0)) > 0.2);
    CHECK_NEAR("mean", s.mean, 3.0, 1e-12);
    CHECK_NEAR("min", s.minimum, 1.0, 0.0);
    CHECK_NEAR("max", s.maximum, 6.0, 0.0);
    CHECK_NEAR("count", s.count, 4.0, 0.0);

    delete_volume(tiny);
  }

  /* ------------------------------------------------------ against the oracle */
  check_stats("whole volume", n3::masked_stats(chunk, NULL), "whole");
  check_stats("in mask", n3::masked_stats(chunk, mask), "masked");

  /* volume_stats' mask rule is value != 0, not value > 0.5 or a threshold
   * (volumeStats.cc:236). */
  {
    n3::Stats all = n3::masked_stats(chunk, NULL);
    n3::Stats inside = n3::masked_stats(chunk, mask);
    CHECK_TRUE("the mask selects a strict subset",
               inside.count > 0 && inside.count < all.count);
  }

  /* ---------------------------------------------------------- save/load */
  {
    VIO_BOOL signed_flag;
    nc_type type = n3::storage_type(data + "/chunk.mnc.gz", &signed_flag);
    CHECK_TRUE("chunk.mnc is 16-bit", type == NC_SHORT);

    /* Written as double and read back, a buffer is unchanged.  This is the
     * property the whole in-memory pipeline rests on. */
    n3::save(chunk, "test_buffers_double.mnc", data + "/chunk.mnc.gz",
             NC_DOUBLE, FALSE, "test_buffers");
    VIO_Volume back = n3::load("test_buffers_double.mnc");
    int n = n3::voxel_count(chunk);
    CHECK_ALL("a double round trip is exact",
              n3::values(back), n3::values(chunk), n, 0.0);
    delete_volume(back);

    /* Written in the file's own storage type it is not, and the amount is the
     * quantum of that type over the volume's range -- the cost the Perl
     * drivers pay at every step. */
    n3::save(chunk, "test_buffers_short.mnc", data + "/chunk.mnc.gz",
             type, signed_flag, "test_buffers");
    back = n3::load("test_buffers_short.mnc");
    n3::Stats s = n3::masked_stats(chunk, NULL);
    CHECK_ALL("a 16-bit round trip costs half a quantum",
              n3::values(back), n3::values(chunk), n,
              0.5 * (s.maximum - s.minimum) / 65535.0);
    delete_volume(back);
  }

  delete_volume(mask);
  delete_volume(chunk);
  return n3check::report("buffers");
}
