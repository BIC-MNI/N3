/* nu_correct_cxx -denoise, driven as a user drives it, on a volume whose
 * non-uniformity is known because this test put it there.
 *
 * The substrate is synthetic rather than chunk.mnc.  A recorded answer taken
 * from this same tree can only say "the output has not changed"; it cannot say
 * the output is right, and a run that estimated a field of exactly 1 and
 * copied its input to its output would reproduce any recording of itself
 * perfectly.  So the input here is built in main(): a two-tissue ellipsoid
 * (anatomy, no shading), multiplied by an analytic smooth field -- a product
 * of cosines, one quarter period on each axis -- and then given additive
 * Gaussian noise.  Because the field is known voxel by voxel, the volume the
 * correction is supposed to produce is known too:
 *
 *     ideal[i] = biased[i] / injected_field(i)
 *
 * read back from the file the driver was given, so both sides carry the same
 * quantisation.  A corrected output is scored against that, and the score is
 * made scale invariant because N3 determines the field only up to a
 * multiplicative constant.  The uncorrected input scored the same way gives
 * the baseline the correction has to beat: it is the injected non-uniformity
 * itself, since biased/ideal is the injected field by construction.
 *
 * Note that the noise cancels out of that score.  Both arms correct the same
 * stored volume, so corrected/ideal is (biased/estimated)/(biased/injected) =
 * injected/estimated: the score measures the field and nothing else.
 *
 * What is asserted:
 *
 *   1. the correction recovers the injected field, with and without
 *      -denoise: the residual non-uniformity of each corrected volume is
 *      a small fraction of the injected non-uniformity.  Measured 0.46% and
 *      0.39% against an injected 7.1%, so the bound of a quarter of the
 *      injection sits about four times above either arm.
 *
 *      This is the check that the pipeline does its job at all, and it is
 *      also the one that catches the defect the option is most exposed to --
 *      the field is estimated from the DENOISED copy while the correction is
 *      applied to the ORIGINAL (nu_correct_cxx.cc:570 is handed `input`, not
 *      `est_input`), and that one argument is what a future simplification
 *      would remove.  Verified by making the change and rebuilding: the
 *      denoised arm's residual rises to 5.7%, three times the bound, while
 *      the plain arm stays at 0.46%, so the failure names the path it is in;
 *   2. the noise level the denoiser estimates for itself is the noise level
 *      this test injected: 15.15 against an injected 15;
 *   3. off by default is INERT.  The whole safety argument for adding the
 *      option is that a run without it is what it was before; est_input
 *      aliases input and nothing is copied (nu_correct_cxx.cc:533).  Asserted
 *      as exact equality, on the corrected volume and on the field, because
 *      that is what "byte for byte what it was" means;
 *   4. -denoise_sigma / -denoise_beta / -denoise_rician cross two translation
 *      units and a library boundary to reach nlm.  A dropped assignment is
 *      invisible, so each is asserted to move the result at all -- a property,
 *      in the style of test_driver_endtoend.cc's -legacy_rounding check, with
 *      no bound fitted to the measurement;
 *   5. the two error paths fail loudly and leave nothing behind.
 *
 * -denoise_threads 1 everywhere.  NLM's block aggregation is partitioned
 * across threads, so its output is reproducible at a fixed thread count and
 * not across counts.  That is pre-existing mincnlm behaviour, not something
 * the -denoise work introduced, but it is why the test pins the count.
 *
 * Not covered here, deliberately:
 *   - nu_estimate_cxx -denoise: built from this same source file
 *     (N3/CMakeLists.txt:327-330), so the parsing is byte-identical, and
 *     N3_DRIVER_BIN points only at nu_correct_cxx;
 *   - a bare -denoise_sigma with no value: run_driver appends the input path
 *     after the options, so the flag would swallow the input and the run would
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
#include <random>
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

/* ---- the phantom's constants, all in the real intensity units it is built
 * in.  The two tissues are a sharp boundary, which is structure and not
 * shading: the only smooth multiplicative term in the volume is the one
 * injected below, so a correct estimate has something to find and nothing
 * else to be confused by.  The noise is not decoration -- a noise-free or
 * single-valued volume is a degenerate histogram for the sharpen
 * deconvolution and diverges (test_driver_properties.cc) -- and here it is
 * also the thing -denoise exists to remove. */
static const double BACKGROUND = 30.0;
static const double OUTER_TISSUE = 200.0;
static const double INNER_TISSUE = 320.0;
static const double NOISE_SIGMA = 15.0;
static const double FIELD_AMPLITUDE = 0.5;

/* The injected non-uniformity: 1 + a*cos*cos*cos, a quarter period on each
 * axis, so the field is smooth and monotone across the volume rather than
 * oscillating within it.  A half period on each axis would put the extremes
 * of the product at the eight corners, all of which are outside the mask,
 * and would leave the field nearly flat over the part of the volume the
 * estimation actually sees. */
static double injected_field(int i0, int i1, int i2, const int sizes[])
{
  double u0 = sizes[0] > 1 ? (double) i0 / (sizes[0] - 1) : 0.0;
  double u1 = sizes[1] > 1 ? (double) i1 / (sizes[1] - 1) : 0.0;
  double u2 = sizes[2] > 1 ? (double) i2 / (sizes[2] - 1) : 0.0;
  return 1.0 + FIELD_AMPLITUDE * cos(0.5*M_PI*u0) * cos(0.5*M_PI*u1)
                               * cos(0.5*M_PI*u2);
}

/* The protocol every run in this file shares.  -stop 0.0 prevents any stage
 * from stopping early, so every run here executes exactly 30 iterations and
 * nothing depends on the stopping rule; -shrink 1 keeps the resampling out of
 * it; -V1.0 pins the protocol independently of whichever default the driver
 * currently selects (test_driver_endtoend.cc).
 *
 * 30 iterations rather than the single one the other driver tests use: the
 * denoised volume is the input to EVERY iteration of the estimation loop, so
 * a one-iteration run measures the option in a regime it is not used in.
 *
 * -distance 100 rather than the 200 the other driver tests use: the phantom
 * is 182 x 156 x 100 mm, so knots 200 mm apart would be a basis too coarse to
 * represent the injected field, and the test would be measuring the spline's
 * reach rather than the estimation. */
static std::string g_mask;

static std::string protocol()
{
  return std::string("-V1.0 -shrink 1 -iterations 30 -stop 0.0 -distance 100"
                     " -mask \"") + g_mask + "\"";
}

static std::string out_path(const char *tag)
{
  return n3fixture::temp_path(std::string("dn_") + tag + ".mnc");
}

/* Run the driver on the biased phantom with the given extra options, logging
 * to <out>.log so cleanup() removes it.  Returns the driver's success. */
static std::string g_input;

static bool run(const std::string &opts, const std::string &out)
{
  return n3fixture::run_driver(opts + " " + protocol(), g_input, out,
                               out + ".log");
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

/* How much non-uniformity is left in `got` relative to the volume the
 * correction was supposed to produce: the in-mask RMS departure of got/ideal
 * from its own mean, as a fraction.  Divided by that mean because N3 fixes
 * the field only up to a multiplicative constant, so an output scaled by 1.05
 * everywhere is a perfect correction, not a 5% error; the constant itself is
 * returned through `scale` and printed rather than asserted. */
static double residual(const std::vector<double> &got,
                       const std::vector<double> &ideal,
                       const double *mv, double *scale)
{
  double sum = 0.0, cnt = 0.0;
  for(size_t i = 0; i < got.size(); i++)
    if(mv[i] > 0.0) { sum += got[i]/ideal[i]; cnt += 1.0; }
  double c = cnt > 0.0 ? sum/cnt : 1.0;
  double d2 = 0.0;
  for(size_t i = 0; i < got.size(); i++)
    if(mv[i] > 0.0)
      {
        double e = got[i]/ideal[i]/c - 1.0;
        d2 += e*e;
      }
  if(scale) *scale = c;
  return cnt > 0.0 ? sqrt(d2/cnt) : 0.0;
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

/* Whole-line match, so that "Denoised with sigma 30" is not satisfied by
 * "Denoised with sigma 30.4". */
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
  /* chunk.mnc is here only as a header donor: its grid, dimension names and
   * direction cosines make the phantom a MINC file the driver will accept,
   * and its storage type makes the quantisation the same one the rest of the
   * suite works at.  None of its voxels are used. */
  const std::string like = std::string(N3_DATA_DIR) + "/chunk.mnc";

  VIO_Volume model = n3::load(like);
  int sizes[VIO_N_DIMENSIONS];
  get_volume_sizes(model, sizes);
  const int n = n3::voxel_count(model);

  /* ---- build the phantom ------------------------------------------------
   * Indexed in values()'s storage order: contiguous along sizes[2], then
   * sizes[1], then sizes[0], so the fastest counter walks the last axis and
   * each radius sits on its own (test_driver_properties.cc). */
  VIO_Volume biased_v = n3::like(model);
  VIO_Volume mask_v = n3::like(model);
  double *bd = n3::values(biased_v), *md = n3::values(mask_v);
  std::vector<double> field(n);

  std::mt19937 rng(20260912);
  std::normal_distribution<double> noise(0.0, NOISE_SIGMA);
  for(int i = 0; i < n; i++)
    {
      int i2 = i % sizes[2], r = i / sizes[2];
      int i1 = r % sizes[1], i0 = r / sizes[1];
      double d0 = (i0 - 0.5*sizes[0]) / (0.40*sizes[0]);
      double d1 = (i1 - 0.5*sizes[1]) / (0.40*sizes[1]);
      double d2 = (i2 - 0.5*sizes[2]) / (0.40*sizes[2]);
      double outer = d0*d0 + d1*d1 + d2*d2;
      double inner = outer * (0.40/0.20) * (0.40/0.20);
      double clean = outer > 1.0 ? BACKGROUND
                   : (inner <= 1.0 ? INNER_TISSUE : OUTER_TISSUE);
      double f = injected_field(i0, i1, i2, sizes);
      field[i] = f;
      /* field first, noise after it: the bias is multiplicative on the signal
       * and the noise is added by the receiver downstream of it. */
      bd[i] = clean*f + noise(rng);
      md[i] = outer <= 1.0 ? 1.0 : 0.0;
    }

  VIO_BOOL sf; nc_type type = n3::storage_type(like, &sf);
  g_input = n3fixture::temp_path("dn_phantom.mnc");
  g_mask = n3fixture::temp_path("dn_phantom_mask.mnc");
  n3fixture::save_phantom(biased_v, g_input, like, type, sf,
                          "test_driver_denoise");
  n3fixture::save_phantom(mask_v, g_mask, like, type, sf,
                          "test_driver_denoise");
  delete_volume(biased_v);
  delete_volume(mask_v);

  /* Read both back: the driver sees the quantised file, so the ground truth
   * has to be derived from the quantised file too. */
  VIO_Volume input = n3::load(g_input);
  VIO_Volume mask = n3::load(g_mask);
  const double *mv = n3::values(mask);
  const double *iv = n3::values(input);
  std::vector<double> biased(iv, iv + n), ideal(n);
  for(int i = 0; i < n; i++) ideal[i] = biased[i] / field[i];

  const std::string one_thread = "-denoise_threads 1";

  std::string p_base = out_path("base"), p_nodn = out_path("nodn"),
              p_dn = out_path("dn"), p_sig = out_path("sig"),
              p_beta = out_path("beta"), p_ric = out_path("ric"),
              p_ethreads = out_path("ethreads"), p_eopt = out_path("eopt");

  bool ok = run("", p_base)
         && run("-nodenoise", p_nodn)
         && run("-denoise -verbose " + one_thread, p_dn)
         && run("-denoise_sigma 30 -verbose " + one_thread, p_sig)
         && run("-denoise_beta 0.5 " + one_thread, p_beta)
         && run("-denoise_rician " + one_thread, p_ric);

  if(!ok)
    {
      printf("driver exited non-zero on a run that must succeed\n");
      n3check::failures()++;
      cleanup(p_base); cleanup(p_nodn); cleanup(p_dn);
      cleanup(p_sig); cleanup(p_beta); cleanup(p_ric);
      delete_volume(mask); delete_volume(input);
      unlink(g_input.c_str()); unlink(g_mask.c_str());
      return n3check::report("driver_denoise");
    }

  std::vector<double> base = voxels(p_base), nodn = voxels(p_nodn),
                      dn = voxels(p_dn), sig = voxels(p_sig),
                      beta = voxels(p_beta), ric = voxels(p_ric);

  /* ---- A. the injected field is recovered -------------------------------
   * The baseline is the input scored against the same ideal, which is the
   * injected non-uniformity itself: biased/ideal is the injected field by
   * construction, so `injected` below is exactly what a correction that did
   * nothing would score.  The bound is a quarter of it -- the assertion is
   * that most of the injected non-uniformity is gone, not that the residual
   * has a particular size, and it is stated as a fraction of the injection
   * rather than as a number fitted to this measurement. */
  {
    double c_in = 0.0, c_base = 0.0, c_dn = 0.0;
    double injected = residual(biased, ideal, mv, &c_in);
    double r_base = residual(base, ideal, mv, &c_base);
    double r_dn = residual(dn, ideal, mv, &c_dn);

    double lo = 0.0, hi = 0.0; bool first = true;
    for(int i = 0; i < n; i++)
      if(mv[i] > 0.0)
        {
          if(first) { lo = hi = field[i]; first = false; }
          else if(field[i] < lo) lo = field[i];
          else if(field[i] > hi) hi = field[i];
        }
    printf("  (injected field spans %.3f to %.3f in the mask, %.1f%% RMS"
           " non-uniformity)\n", lo, hi, 100.0*injected);
    printf("  (residual after correction: %.1f%% plain, %.1f%% denoised;"
           " output scale %.4f and %.4f)\n",
           100.0*r_base, 100.0*r_dn, c_base, c_dn);

    n3check::record("the plain correction recovers the injected field",
                    r_base < 0.25*injected, r_base, 0.25*injected);
    n3check::record("the -denoise correction recovers the injected field",
                    r_dn < 0.25*injected, r_dn, 0.25*injected);
  }

  /* ---- B. the denoiser estimates the noise that was injected ------------
   * Not a recorded number: NOISE_SIGMA is what this test added, and the
   * estimator (Coupe 2009: MAD of the finest wavelet sub-band over the
   * detected object) has to find it in the volume.  The bound is a factor of
   * two either way, which is a statement about the estimator being right to
   * within its own modelling assumptions rather than a tolerance fitted to
   * the measurement. */
  {
    double got = sigma_from_log(p_dn);
    printf("  (estimated noise sigma %g, injected %g)\n", got, NOISE_SIGMA);
    n3check::record("the estimated noise sigma is the injected one",
                    got > 0.5*NOISE_SIGMA && got < 2.0*NOISE_SIGMA,
                    got, NOISE_SIGMA);
  }

  /* ---- C. off by default is inert --------------------------------------
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

  /* ---- D. every option reaches nlm -------------------------------------
   * Strict inequalities with the measured margin printed, not fitted bounds:
   * what is being asserted is that the option is plumbed through at all,
   * which is what a dropped assignment fails. */
  {
    double d_dn = rel(dn, base), d_sig = rel(sig, dn),
           d_beta = rel(beta, dn), d_ric = rel(ric, dn);
    printf("  (-denoise vs plain: %.3e; vs the default -denoise run:"
           " sigma %.3e, beta %.3e, rician %.3e)\n",
           d_dn, d_sig, d_beta, d_ric);
    CHECK_TRUE("-denoise is not inert: it changes the corrected volume", d_dn > 0.0);
    CHECK_TRUE("-denoise_sigma changes the result", d_sig > 0.0);
    /* and it is used verbatim rather than re-estimated */
    CHECK_TRUE("-denoise_sigma 30 is used as given, not re-estimated",
               log_has_line(p_sig, "Denoised with sigma 30"));
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
  delete_volume(model);
  unlink(g_input.c_str());
  unlink(g_mask.c_str());
  return n3check::report("driver_denoise");
}
