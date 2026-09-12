/* nu_correct_cxx -denoise, driven as a user drives it.
 *
 * The in-process test of the same feature (test_denoise.cc) pins the
 * VIO_Volume <-> nlm bridge on a synthetic phantom and nothing else.  Three
 * things about the option are only observable through the binary, and all
 * three are the kind that break silently:
 *
 *   1. off by default is INERT.  The whole safety argument for adding the
 *      option is that a run without it is what it was before; est_input
 *      aliases input and nothing is copied (nu_correct_cxx.cc:533).  Asserted
 *      as exact equality, on the corrected volume and on the field, because
 *      that is what "byte for byte what it was" means.
 *   2. the field is estimated from the DENOISED copy while the correction is
 *      applied to the ORIGINAL (nu_correct_cxx.cc:570 is handed `input`, not
 *      `est_input`).  That one argument is what a future tidy-up would
 *      "simplify", and every other test here would stay green while the
 *      output quietly became a denoised image.  Pinned by the identity
 *      corrected * field == original, which holds to 5.5e-04 today -- twice
 *      the stored output's own quantisation and no more.  Measured, not
 *      argued: passing est_input to nu_evaluate and rebuilding moves that
 *      residual to 6.4e-02, 23x this file's bound, while the control run
 *      stays at 3.1e-04.  That is the principal assertion of this file.
 *   3. -denoise_sigma / -denoise_beta / -denoise_rician cross two translation
 *      units and a library boundary to reach nlm.  A dropped assignment is
 *      invisible, so each is asserted to move the result at all -- a property,
 *      in the style of test_driver_endtoend.cc's -legacy_rounding check, with
 *      no bound fitted to the measurement.
 *
 * -denoise_threads 1 everywhere.  NLM's block aggregation is partitioned
 * across threads, so its output is reproducible at a fixed thread count and
 * not across counts: at this protocol 1 vs 2 threads moved the corrected
 * volume by up to 7.2e3 on a 9.0e5 range, and at a single iteration by far
 * more.  That is pre-existing mincnlm behaviour, not something the -denoise
 * work introduced, but it is why both this test and the fixture pin the
 * count.
 *
 * Not covered here, deliberately:
 *   - nu_estimate_cxx -denoise: built from this same source file
 *     (N3/CMakeLists.txt:327-330), so the parsing is byte-identical, and
 *     N3_DRIVER_BIN points only at nu_correct_cxx;
 *   - a bare -denoise_sigma with no value: run_driver appends the input path
 *     after the options, so the flag would swallow chunk.mnc and the run would
 *     fail for the wrong reason;
 *   - the N3_WITH_NLM stub message: unreachable from a build in which this
 *     test exists at all (the CMake block is inside IF(TARGET nlm)).
 */

#include "check.h"
#include "fixture.h"

#include "../../src/N3Pipeline/Buffers.h"
#include "../../src/N3Pipeline/FitField.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>
#include <vector>

#ifndef N3_DRIVER_BIN
#error "N3_DRIVER_BIN must be the path of the built nu_correct_cxx binary"
#endif
#ifndef N3_DATA_DIR
#error "N3_DATA_DIR must be the data directory (testing/)"
#endif

using n3fixture::imp_of;
using n3fixture::cleanup;

/* The protocol every run in this file shares.  -stop 0.0 prevents any stage
 * from stopping early, so every run here executes exactly 30 iterations and
 * nothing depends on the stopping rule; -shrink 1 keeps the resampling out of
 * it; -V1.0 pins the protocol independently of whichever default the driver
 * currently selects (test_driver_endtoend.cc).
 *
 * 30 iterations rather than the single one the other driver tests use: the
 * denoised volume is the input to EVERY iteration of the estimation loop, so
 * a one-iteration run measures the option in a regime it is not used in.  The
 * six successful runs take about 33 s together at this count. */
static std::string protocol()
{
  return std::string("-V1.0 -shrink 1 -iterations 30 -stop 0.0 -distance 200"
                     " -mask \"") + N3_DATA_DIR + "/chunk_mask.mnc\"";
}

static std::string out_path(const char *tag)
{
  return n3fixture::temp_path(std::string("dn_") + tag + ".mnc");
}

/* Run the driver on chunk.mnc with the given extra options, logging to
 * <out>.log so cleanup() removes it.  Returns the driver's success. */
static bool run(const std::string &opts, const std::string &out)
{
  return n3fixture::run_driver(opts + " " + protocol(),
                               std::string(N3_DATA_DIR) + "/chunk.mnc",
                               out, out + ".log");
}

static std::vector<double> voxels(const std::string &path)
{
  VIO_Volume v = n3::load(path);
  const double *d = n3::values(v);
  std::vector<double> out(d, d + n3::voxel_count(v));
  delete_volume(v);
  return out;
}

/* The field the run recorded in its own .imp, evaluated back onto the input
 * grid -- the field itself, not a quotient of the input and an output that is
 * input/field by construction (test_driver_properties.cc). */
static std::vector<double> field_of(const std::string &out,
                                    VIO_Volume model, VIO_Volume mask)
{
  VIO_Volume f = n3::like(model);
  n3::evaluate_saved_field(imp_of(out), f, mask);
  const double *d = n3::values(f);
  std::vector<double> v(d, d + n3::voxel_count(f));
  delete_volume(f);
  return v;
}

/* The precision the identity can hold to.  It passes through two files: the
 * corrected volume is stored in the input's type (12 bits for chunk.mnc, from
 * chunk_valid_range.txt) and the field is read back from the .imp.  Half a
 * quantum of the output's range, relative to the in-mask RMS of the input, is
 * the first of those and is computed here from the data rather than assumed --
 * a fixed figure would be wrong for a file stored at any other depth.  It is
 * the same derivation test_driver_endtoend.cc's round_trip_bound uses.
 *
 * The .imp round trip is not separately derivable, so the assertions take ten
 * times this figure.  For scale: the default run measures 1.1x it, the
 * -denoise run 1.9x (a more strongly varying field multiplies the stored
 * volume's rounding error), and correcting the wrong volume measures 225x it.
 * The factor of ten therefore sits well clear of both sides. */
static double quantisation(VIO_Volume input, VIO_Volume mask)
{
  n3::Stats whole = n3::masked_stats(input, NULL);
  n3::Stats inside = n3::masked_stats(input, mask);
  double rms = sqrt(inside.mean*inside.mean + inside.stddev*inside.stddev);
  return 0.5 * (whole.maximum - whole.minimum)
       / n3fixture::valid_steps("chunk_valid_range.txt") / rms;
}

/* rel RMS of corrected * field against the original, over the mask only:
 * outside it the field is zero by construction (evaluate_saved_field). */
static double identity(const std::vector<double> &corr,
                       const std::vector<double> &field,
                       const std::vector<double> &orig,
                       const double *mv)
{
  double sd = 0.0, sb = 0.0;
  for(size_t i = 0; i < orig.size(); i++)
    if(mv[i] > 0.0)
      { double d = corr[i]*field[i] - orig[i]; sd += d*d; sb += orig[i]*orig[i]; }
  return sb > 0.0 ? sqrt(sd/sb) : 1.0/0.0;
}

static double rel(const std::vector<double> &a, const std::vector<double> &b)
{
  return n3check::rel_rms(&a[0], &b[0], (int) a.size());
}

/* The sigma the driver reported, or -1 if it reported none. */
static double sigma_from_log(const std::string &out)
{
  FILE *f = fopen((out + ".log").c_str(), "r");
  if(!f) return -1.0;
  char line[512];
  double sigma = -1.0;
  while(fgets(line, sizeof(line), f))
    if(sscanf(line, "Denoised with sigma %lf", &sigma) == 1) break;
  fclose(f);
  return sigma;
}

/* Whole-line match, so that "Denoised with sigma 5000" is not satisfied by
 * "Denoised with sigma 50000". */
static bool log_has_line(const std::string &out, const char *want)
{
  FILE *f = fopen((out + ".log").c_str(), "r");
  if(!f) return false;
  char line[512];
  bool found = false;
  while(!found && fgets(line, sizeof(line), f))
    {
      size_t n = strlen(line);
      while(n > 0 && (line[n-1] == '\n' || line[n-1] == '\r')) line[--n] = '\0';
      found = strcmp(line, want) == 0;
    }
  fclose(f);
  return found;
}

static bool log_contains(const std::string &out, const char *want)
{
  FILE *f = fopen((out + ".log").c_str(), "r");
  if(!f) return false;
  char line[512];
  bool found = false;
  while(!found && fgets(line, sizeof(line), f)) found = strstr(line, want) != NULL;
  fclose(f);
  return found;
}

int main()
{
  const std::string data = N3_DATA_DIR;
  const std::string inp = data + "/chunk.mnc";
  const std::string mask_in = data + "/chunk_mask.mnc";

  VIO_Volume input = n3::load(inp);
  VIO_Volume mask = n3::load(mask_in);
  const double *mv = n3::values(mask);

  const std::string one_thread = "-denoise_threads 1";

  std::string p_base = out_path("base"), p_nodn = out_path("nodn"),
              p_dn = out_path("dn"), p_sig = out_path("sig"),
              p_beta = out_path("beta"), p_ric = out_path("ric"),
              p_ethreads = out_path("ethreads"), p_eopt = out_path("eopt");

  bool ok = run("", p_base)
         && run("-nodenoise", p_nodn)
         && run("-denoise -verbose " + one_thread, p_dn)
         && run("-denoise_sigma 5000 -verbose " + one_thread, p_sig)
         && run("-denoise_beta 0.5 " + one_thread, p_beta)
         && run("-denoise_rician " + one_thread, p_ric);

  if(!ok)
    {
      printf("driver exited non-zero on a run that must succeed\n");
      n3check::failures()++;
      cleanup(p_base); cleanup(p_nodn); cleanup(p_dn);
      cleanup(p_sig); cleanup(p_beta); cleanup(p_ric);
      delete_volume(mask); delete_volume(input);
      return n3check::report("driver_denoise");
    }

  std::vector<double> orig(n3::values(input),
                           n3::values(input) + n3::voxel_count(input));
  std::vector<double> base = voxels(p_base), nodn = voxels(p_nodn),
                      dn = voxels(p_dn), sig = voxels(p_sig),
                      beta = voxels(p_beta), ric = voxels(p_ric);

  /* ---- A. the recorded answer ------------------------------------------
   * Unlike every other oracle in reference/ this one is not a legacy answer:
   * the Perl nu_correct has no -denoise, so there is nothing to record it
   * against, and regenerate_reference.sh's cycle 16 says so.  It is a
   * regression lock on this tree's own output, which is worth having because
   * the run is deterministic at a pinned thread count. */
  {
    std::vector<double> ours = n3fixture::strided(&dn[0], (int) dn.size());
    std::vector<double> oracle = n3fixture::read_f64("nu_correct_denoise.f64");
    n3fixture::must(ours.size() == oracle.size(),
                    "nu_correct_denoise.f64: size mismatch");
    /* 1e-4 rather than exact equality: the comparison should survive a
     * compiler or libm change, not only this machine.  It still
     * discriminates -- the printed figure below is how far -denoise moves the
     * corrected volume at this protocol, and the bound is well under it. */
    CHECK_RMS("-denoise output matches the recorded answer",
              &ours[0], &oracle[0], (int) ours.size(), 1e-4);

    printf("  (-denoise moves the corrected volume by %.3e)\n", rel(dn, base));

    double got = sigma_from_log(p_dn), want = n3fixture::read_scalar("denoise_sigma.txt");
    printf("  (auto-estimated sigma %g, recorded %g)\n", got, want);
    CHECK_NEAR("the auto-estimated sigma matches the recorded one",
               got, want, fabs(want) * 1e-6);
  }

  /* ---- B. off by default is inert --------------------------------------
   * Exact equality, not a bound: est_input aliases input when -denoise is
   * off, so the default path executes the same code on the same buffer and
   * any difference at all is a defect, not a tolerance question. */
  {
    double d = rel(nodn, base);
    printf("  (-nodenoise vs the default run: %.3e)\n", d);
    CHECK_TRUE("-nodenoise leaves the corrected volume exactly unchanged", d == 0.0);

    std::vector<double> f_nodn = field_of(p_nodn, input, mask);
    std::vector<double> f_base = field_of(p_base, input, mask);
    CHECK_TRUE("-nodenoise leaves the estimated field exactly unchanged",
               rel(f_nodn, f_base) == 0.0);
  }

  /* ---- C. estimate from the denoised copy, correct the original -------- */
  {
    std::vector<double> f_dn = field_of(p_dn, input, mask);
    std::vector<double> f_base = field_of(p_base, input, mask);

    double moved = rel(f_dn, f_base);
    printf("  (-denoise moves the estimated field by %.3e)\n", moved);
    CHECK_TRUE("-denoise is not inert: it changes the estimated field", moved > 0.0);

    /* The centrepiece.  Both runs must satisfy corrected * field == original;
     * the -denoise one is the assertion, the default one is the control that
     * says the bound measures the correction's own round trip and not the
     * denoising.  A regression that corrected est_input instead of input
     * would leave the removed noise in the residual: the same build with
     * est_input passed to nu_evaluate measures 6.375e-02 here, 23x this
     * bound, against 5.489e-04 for the correct one -- and the control below
     * stays unmoved at 3.102e-04, so the check separates the two volumes and
     * not the two runs. */
    double q = quantisation(input, mask);
    double bound = 10.0 * q;
    double m_dn = identity(dn, f_dn, orig, mv);
    double m_base = identity(base, f_base, orig, mv);
    printf("  (a half quantum of the stored output is worth %.3e;\n"
           "   corrected * field vs original: %.3e with -denoise, %.3e without;\n"
           "   correcting the denoised copy instead measures 6.375e-02)\n",
           q, m_dn, m_base);
    n3check::record("-denoise corrects the ORIGINAL, not the denoised copy",
                    m_dn < bound, m_dn, bound);
    n3check::record("  (control) the same identity for the default run",
                    m_base < bound, m_base, bound);
  }

  /* ---- D. every option reaches nlm -------------------------------------
   * Strict inequalities with the measured margin printed, not fitted bounds:
   * what is being asserted is that the option is plumbed through at all,
   * which is what a dropped assignment fails. */
  {
    double d_sig = rel(sig, dn), d_beta = rel(beta, dn), d_ric = rel(ric, dn);
    printf("  (vs the default -denoise run: sigma %.3e, beta %.3e, rician %.3e)\n",
           d_sig, d_beta, d_ric);
    CHECK_TRUE("-denoise_sigma changes the result", d_sig > 0.0);
    /* and it is used verbatim rather than re-estimated */
    CHECK_TRUE("-denoise_sigma 5000 is used as given, not re-estimated",
               log_has_line(p_sig, "Denoised with sigma 5000"));
    CHECK_TRUE("-denoise_beta changes the result", d_beta > 0.0);
    CHECK_TRUE("-denoise_rician changes the result", d_ric > 0.0);
  }

  /* ---- E. the error paths ----------------------------------------------
   * nlm reports an unusable parameter combination by throwing; the driver
   * translates that to its own die() (nu_correct_cxx.cc:546).  Both paths
   * must fail loudly and leave nothing behind -- a half-written output would
   * be worse than no output. */
  {
    bool failed = !run("-denoise -denoise_threads 0", p_ethreads);
    bool said = log_contains(p_ethreads, "nlm::denoise: threads must be >= 1");
    bool nothing = access(p_ethreads.c_str(), F_OK) != 0;
    printf("  (-denoise_threads 0: failed %d, reported %d, no output %d)\n",
           failed, said, nothing);
    CHECK_TRUE("-denoise_threads 0 fails with nlm's message and writes nothing",
               failed && said && nothing);

    bool failed2 = !run("-denoise_bogus", p_eopt);
    bool said2 = log_contains(p_eopt, "unknown option -denoise_bogus");
    bool nothing2 = access(p_eopt.c_str(), F_OK) != 0;
    printf("  (-denoise_bogus: failed %d, reported %d, no output %d)\n",
           failed2, said2, nothing2);
    CHECK_TRUE("an unknown -denoise* option is rejected, not ignored",
               failed2 && said2 && nothing2);
  }

  cleanup(p_base); cleanup(p_nodn); cleanup(p_dn);
  cleanup(p_sig); cleanup(p_beta); cleanup(p_ric);
  cleanup(p_ethreads); cleanup(p_eopt);
  delete_volume(mask);
  delete_volume(input);
  return n3check::report("driver_denoise");
}
