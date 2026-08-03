/* Cycle 1: ShrinkVolume.
 *
 * First cycle after the harness, because every later comparison happens on the
 * grid this produces: a wrong grid invalidates all of them, and would do so
 * silently.
 *
 * Two oracles per factor, both from the legacy programs (regenerate_reference.sh):
 *
 *   shrinkN.geom            the grid the Perl driver's own ShrinkVolume made
 *   shrinkN_double.f64      the same mincresample run written as double
 *   shrinkN_quantised.f64   the Perl's actual output, 12-bit like every
 *                           intermediate the drivers pass between programs
 *
 * The first two are held exactly.  The third is the divergence the in-memory
 * pipeline exists to remove, and appears here as a bound on the oracle rather
 * than on the code: the Perl's shrunk volume sits within half a 12-bit quantum
 * of the unquantised resample of the same data.
 *
 * That load() is right is not asserted separately.  It cannot be wrong and
 * still let shrink() reproduce 250k resampled voxels exactly.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"

#include <string>
#include <vector>

static const int FACTORS[] = { 2, 3, 4 };
static const int N_FACTORS = 3;

static void check_geometry(const char *what, VIO_Volume volume,
                           const std::vector<double> &geom)
{
  int sizes[VIO_N_DIMENSIONS];
  VIO_Real starts[VIO_N_DIMENSIONS], steps[VIO_N_DIMENSIONS];
  get_volume_sizes(volume, sizes);
  get_volume_starts(volume, starts);
  get_volume_separations(volume, steps);

  n3fixture::must(geom.size() == 9, std::string(what) + ": short geometry fixture");

  bool ok = true;
  for(int i = 0; i < VIO_N_DIMENSIONS; i++)
    {
      ok = ok && (sizes[i]  == (int) geom[i]);
      ok = ok && (starts[i] == geom[3 + i]);
      ok = ok && (steps[i]  == geom[6 + i]);
    }
  if(!ok)
    printf("  sizes %d %d %d (want %d %d %d)  steps %g %g %g (want %g %g %g)\n",
           sizes[0], sizes[1], sizes[2],
           (int) geom[0], (int) geom[1], (int) geom[2],
           steps[0], steps[1], steps[2], geom[6], geom[7], geom[8]);
  CHECK_TRUE(what, ok);
}

int main()
{
  std::string data = std::string(N3_DATA_DIR) + "/chunk.mnc.gz";
  VIO_Volume chunk = n3::load(data);

  check_geometry("the input grid is read as the file has it",
                 chunk, n3fixture::read_text("chunk.geom"));

  for(int f = 0; f < N_FACTORS; f++)
    {
      char name[64], what[128];
      sprintf(name, "shrink%d", FACTORS[f]);

      VIO_Volume small = n3::shrink(chunk, (double) FACTORS[f]);

      sprintf(what, "shrink %d: the grid is the driver's own", FACTORS[f]);
      check_geometry(what, small, n3fixture::read_text(std::string(name) + ".geom"));

      std::vector<double> oracle =
        n3fixture::read_f64(std::string(name) + "_double.f64");
      int n = n3::voxel_count(small);
      sprintf(what, "shrink %d: %d voxels equal mincresample -double", FACTORS[f], n);
      n3fixture::must((int) oracle.size() == n, "size mismatch against oracle");
      CHECK_ALL(what, n3::values(small), &oracle[0], n, 0.0);

      /* The oracle against itself: what a MINC round trip costs.  chunk.mnc
       * is 12-bit (valid_range 0 4095), so half a quantum over the resampled
       * volume's own range is the most the Perl's output can differ from the
       * same resample written as double. */
      std::vector<double> quantised =
        n3fixture::read_f64(std::string(name) + "_quantised.f64");
      double lo = oracle[0], hi = oracle[0];
      for(int i = 0; i < n; i++)
        { if(oracle[i] < lo) lo = oracle[i]; if(oracle[i] > hi) hi = oracle[i]; }
      sprintf(what, "shrink %d: the Perl's own output is 12-bit", FACTORS[f]);
      CHECK_ALL(what, &quantised[0], &oracle[0], n, 0.5 * (hi - lo) / 4095.0);

      delete_volume(small);
    }

  delete_volume(chunk);
  return n3check::report("shrink");
}
