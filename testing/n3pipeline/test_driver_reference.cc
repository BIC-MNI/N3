/* Cycle 15: nu_correct_cxx against nu_reference_1's own oracle-generation
 * protocol, run at the real default iteration count instead of cycle 14's
 * single -stop 0.0 pass.
 *
 * ../compare_nu_result.pl (nu_reference_1) drives nu_estimate then a bare
 * nu_evaluate, and needs two golden volumes -- brain_nu_ref.mnc and
 * brain_nu_ref_external_blas.mnc -- because that pipeline's dsysv_ solve runs
 * through separate installed binaries and MINC file round trips at every
 * stage, so a bundled-vs-external-BLAS difference compounds over 30-40
 * iterations into an RMS the tests must not blur together with the
 * deconvolution being wrong (BLAS.md; EBTKS/CMakeLists.txt).
 *
 * nu_correct_cxx cannot take that same -mask, though: nu_correct forwards
 * -mask to nu_evaluate as well as nu_estimate (nu_estimate.in:62-63), and
 * evaluate_field requires the mask and the volume it is -like to already
 * share a grid, with no resampling (evaluateField.cc's compareVolumes) --
 * confirmed by running real nu_correct with the ICBM mask on brain.mnc: it
 * crashes with "Mask volume and input volume must be the same size."
 * nu_reference_1 only gets away with that mask because compare_nu_result.pl
 * never passes it to the second, separate nu_evaluate call. nu_correct_cxx,
 * like nu_correct, offers no such split, so this cycle uses brain_mask.mnc,
 * which already shares brain.mnc's grid, and compares against nu_correct
 * itself (not the two-step script) run the same way
 * (regenerate_reference.sh, cycle 15).
 *
 * Unlike the Perl pipeline, one golden volume covers every EBTKS_BLAS_BACKEND
 * here. Measured directly, running both nu_correct and nu_correct_cxx under
 * both the bundled and lapacke backends on brain.mnc/brain_mask.mnc: the
 * driver-vs-nu_correct gap is 8.1e-4 to 1.7e-3 (in-memory double precision
 * throughout vs. the Perl's per-stage MINC round trips) and the BLAS choice
 * moves the driver's own answer by only 1.2e-5 -- two orders under that gap,
 * because nothing here is quantised to a reduced MINC storage type between
 * iterations the way the Perl's intermediates are. All four combinations of
 * {backend the driver was built with} x {backend the oracle was recorded
 * under} land at 8.1e-4 to 1.7e-3; the bound below sits at three times the
 * worst of those, comfortably clear of it without pinning either backend's
 * rounding.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"

#ifndef N3_DRIVER_BIN
#error "N3_DRIVER_BIN must be the path of the built nu_correct_cxx binary"
#endif
#ifndef N3_DATA_DIR
#error "N3_DATA_DIR must be the directory holding brain.mnc/brain_mask.mnc"
#endif

static const double BOUND = 5e-3;

int main()
{
  std::string data = N3_DATA_DIR;
  std::string out = n3fixture::temp_path("brain_nu_correct.mnc");

  bool ok = n3fixture::run_driver("-V1.0 -mask \"" + data + "/brain_mask.mnc\"",
                                  data + "/brain.mnc", out);
  if(!ok)
    {
      printf("FAIL: nu_correct_cxx exited non-zero\n");
      n3check::failures()++;
      return n3check::report("driver_reference");
    }

  VIO_Volume corrected = n3::load(out);
  std::vector<double> mine =
    n3fixture::strided(n3::values(corrected), n3::voxel_count(corrected));
  delete_volume(corrected);

  std::vector<double> oracle = n3fixture::read_f64("brain_nu_correct.f64");
  n3fixture::must(mine.size() == oracle.size(),
                  "brain_nu_correct.f64: size mismatch");

  CHECK_RMS("nu_correct_cxx matches nu_correct on brain.mnc/brain_mask.mnc",
            &mine[0], &oracle[0], (int) mine.size(), BOUND);

  n3fixture::cleanup(out);
  return n3check::report("driver_reference");
}
