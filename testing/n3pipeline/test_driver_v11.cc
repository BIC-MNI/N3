/* -V1.1, run to its own real convergence, on brain.mnc.
 *
 * -V1.1 is the driver's own default (nu_correct_cxx.cc: version = 2), not a
 * Perl-compatible protocol: nu_estimate.in defines only -V0.9 and -V1.0
 * (nu_estimate.in:197-199), so there is no installed nu_correct run and no
 * oracle volume for -V1.1 the way cycles 14/15 have for -V1.0.  Every other
 * driver test in this directory pins itself to -V1.0 -nolegacy_rounding
 * specifically so that it is not affected by "whichever protocol the driver's
 * own implicit default currently selects" (test_driver_fwhm.cc,
 * test_driver_properties.cc) -- which also means none of them ever drives the
 * protocol that actually ships as the default.  This file is the one that
 * does, asserting properties rather than a bound against a recorded answer
 * (PLAN §7 cycle 13's approach, applied to the one protocol cycle 13 does not
 * cover):
 *
 *   1. the field, read from the driver's own .imp, is strictly positive and
 *      finite in the mask;
 *   2. in-mask CV of the corrected output is below the input's;
 *   3. -stop 1e-5 ends the run before the -iterations 1000 ceiling does --
 *      the first test in this directory to let the outer loop stop on its own
 *      data-driven criterion at the real default tolerance rather than being
 *      capped by an explicit small -iterations value;
 *   4. omitting -V entirely reproduces -V1.1 byte for bit, i.e. -V1.1 really
 *      is what "no -V flag" runs.
 *
 * One driver invocation covers 1-3: -verbose on a -correct run (not only
 * -estimate_only) prints nu_estimate's own per-iteration line
 * (NuEstimate.cc:201-203), which already carries the stopping iteration.
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"
#include "../../src/N3Pipeline/FitField.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>

#ifndef N3_DRIVER_BIN
#error "N3_DRIVER_BIN must be the path of the built nu_correct_cxx binary"
#endif
#ifndef N3_DATA_DIR
#error "N3_DATA_DIR must be the directory holding brain.mnc/brain_mask.mnc"
#endif

using n3fixture::imp_of;
using n3fixture::cleanup;

/* -iterations 1000 is -V1.1's own ceiling (nu_correct_cxx.cc:390); this test
 * asserts convergence stops strictly short of it. */
static const int CEILING = 1000;

int main()
{
  const std::string data = N3_DATA_DIR;
  const std::string inp = data + "/brain.mnc";
  const std::string mask_path = data + "/brain_mask.mnc";
  const std::string v11 = "-V1.1 -mask \"" + mask_path + "\"";

  VIO_Volume input = n3::load(inp);
  VIO_Volume mask = n3::load(mask_path);
  int n = n3::voxel_count(input);
  const double *mv = n3::values(mask);

  std::string out = n3fixture::temp_path("v11.mnc");
  std::string log = n3fixture::temp_path("v11.log");
  if(!n3fixture::run_driver(v11 + " -verbose", inp, out, log))
    {
      printf("driver did not finish the -V1.1 correction on brain.mnc\n");
      n3check::failures()++;
      delete_volume(mask); delete_volume(input);
      unlink(log.c_str());
      return n3check::report("driver_v11");
    }

  VIO_Volume corr = n3::load(out.c_str());

  /* ---- 1. the field is strictly positive and finite in the mask ---- */
  {
    VIO_Volume field = n3::like(input);
    n3::evaluate_saved_field(imp_of(out), field, mask);
    const double *fv = n3::values(field);
    int nonfinite = 0, nonpositive = 0;
    for(int i = 0; i < n; i++)
      if(mv[i] > 0.0)
        {
          if(!std::isfinite(fv[i])) nonfinite++;
          else if(!(fv[i] > 0.0)) nonpositive++;
        }
    printf("  in-mask field: %d non-finite, %d non-positive values\n",
           nonfinite, nonpositive);
    CHECK_TRUE("the -V1.1 field is strictly positive in the mask",
               nonfinite == 0 && nonpositive == 0);
    delete_volume(field);
  }

  /* ---- 2. in-mask CV of the output below the input's ---- */
  {
    n3::Stats in = n3::masked_stats(input, mask);
    n3::Stats o  = n3::masked_stats(corr, mask);
    double cvin = in.stddev / in.mean, cvout = o.stddev / o.mean;
    CHECK_TRUE("in-mask CV of the -V1.1 output is below the input's",
               cvout < cvin);
    printf("  CV in %g -> CV out %g\n", cvin, cvout);
  }

  /* ---- 3. -stop 1e-5 ends the run short of the -iterations ceiling ---- */
  {
    FILE *f = fopen(log.c_str(), "r");
    int last_iter = -1;
    char line[256];
    while(f && fgets(line, sizeof(line), f))
      {
        int it; double change;
        if(sscanf(line, "CV for change in field estimate at iteration %d: %lf",
                  &it, &change) == 2)
          last_iter = it;
      }
    if(f) fclose(f);
    int total = last_iter + 1;   /* NuEstimate.cc's iter is 0-based. */
    printf("  -V1.1 stopped after %d iteration%s (ceiling %d)\n",
           total, total == 1 ? "" : "s", CEILING);
    CHECK_TRUE("the stop threshold ends the run before the iteration ceiling",
               last_iter >= 0 && total < CEILING);
  }
  unlink(log.c_str());

  /* ---- 4. no -V flag at all reproduces -V1.1 byte for bit ---- */
  {
    std::string out_default = n3fixture::temp_path("v11_default.mnc");
    bool ok = n3fixture::run_driver("-mask \"" + mask_path + "\"", inp,
                                    out_default);
    if(!ok)
      {
        printf("driver did not finish the no--V default correction on "
               "brain.mnc\n");
        n3check::failures()++;
      }
    else
      {
        VIO_Volume def = n3::load(out_default.c_str());
        double d = n3check::rel_rms(n3::values(corr), n3::values(def), n);
        printf("  rel RMS, -V1.1 vs no -V flag: %.3e\n", d);
        CHECK_TRUE("the implicit default agrees with -V1.1 byte for bit",
                   d == 0.0);
        delete_volume(def);
      }
    cleanup(out_default);
  }

  delete_volume(corr);
  delete_volume(mask);
  delete_volume(input);
  cleanup(out);
  return n3check::report("driver_v11");
}
