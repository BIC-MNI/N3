/* Cycle 13: end-to-end properties of the whole pipeline, asserted by driving
 * the compiled binary as a user would.  None of these needs a recorded oracle;
 * they are properties (PLAN §7 cycle 13), so they can be asserted before any
 * end-to-end bound exists -- and they survive iteration and machine changes the
 * way a measured bound does not.
 *
 *   1. the field, read from the driver's own .imp, is strictly positive and
 *      finite everywhere in the mask (a field that reached zero would divide
 *      the output to infinity); this checks the field itself, not a quotient
 *      of the input and an output that is input/field by construction;
 *   2. the in-mask coefficient of variation of the corrected output is below
 *      that of the input, on a real volume (chunk) with real inhomogeneity;
 *   3. a volume with NO planted non-uniformity -- a two-tissue phantom plus
 *      additive noise, anatomy with no smooth multiplicative bias for N3 to
 *      remove -- yields a bounded, finite field whose variation stays below
 *      the phantom's own tissue contrast (the estimator must not absorb the
 *      structure it is meant to be invariant to).  Neither this port nor the
 *      legacy reaches a flat field here (see the assertion's comment);
 *   4. iteration counts at -stop 0: total_iterations() is the LAST staged count
 *      (Perl :1630-1631), so with -stop 0 a stage can never stop early and the
 *      driver reports exactly the requested count.
 */

#include "check.h"

#include "../../src/N3Pipeline/Buffers.h"
#include "../../src/N3Pipeline/FitField.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

#ifndef N3_DRIVER_BIN
#error "N3_DRIVER_BIN must be the path of the built nu_correct_cxx binary"
#endif
#ifndef N3_DATA_DIR
#error "N3_DATA_DIR must be the data directory (testing/)"
#endif

static std::string outdir()
{
  const char *tmp = getenv("TMPDIR");
  return std::string(tmp ? tmp : "/tmp");
}

static std::string tag_path(const char *ext)
{
  char p[256];
  snprintf(p, sizeof(p), "%s/n3cxx_prop_%d%s",
           outdir().c_str(), (int) getpid(), ext);
  return p;
}

/* The .imp the driver writes beside a correct output: the final extension is
 * replaced by .imp (nu_correct_cxx.cc::imp_path, MNI::PathUtilities::replace_ext
 * = s/\.[^\.]*$/\.imp/), so _corr.mnc -> _corr.imp, never _corr.mnc.imp. */
static std::string imp_of(const std::string &path)
{
  size_t dot = path.find_last_of('.');
  return (dot == std::string::npos ? path : path.substr(0, dot)) + ".imp";
}

/* Remove a correct run's outputs: the volume, its .imp, and any log we asked
 * the shell to write as <volume>.log. */
static void cleanup(const std::string &out_path)
{
  unlink(out_path.c_str());
  unlink(imp_of(out_path).c_str());
  unlink((out_path + ".log").c_str());
}

int main()
{
  const std::string data = N3_DATA_DIR;
  const std::string inp = data + "/chunk.mnc.gz";
  const std::string mask_in = data + "/chunk_mask.mnc.gz";

  /* ---- a correction run of the real volume: properties 1 and 2 ---- */
  VIO_Volume input = n3::load(inp);
  VIO_Volume mask = n3::load(mask_in);
  int n = n3::voxel_count(input);
  const double *vi = n3::values(input);
  const double *mv = n3::values(mask);

  std::string corr_path = tag_path("_corr.mnc");
  char cmd[1600];
  snprintf(cmd, sizeof(cmd),
           "\"%s\" -shrink 2 -iterations 20 -stop 0.001 -distance 200 "
           "-mask \"%s\" \"%s\" \"%s\" -clobber",
           N3_DRIVER_BIN, mask_in.c_str(), inp.c_str(), corr_path.c_str());
  if(system(cmd) != 0)
    {
      printf("driver did not finish the correction on chunk\n");
      n3check::failures()++;
      return n3check::report("driver_properties");
    }

  VIO_Volume corr = n3::load(corr_path.c_str());
  const double *vo = n3::values(corr);

  /* ---- 1. the field is strictly positive in the mask ---- */
  {
    /* The driver's correct run writes the .imp (imp_path), and the field is
     * what it holds, evaluated onto the input grid: a positivity test on the
     * field itself, not a quotient whose numerator/denominator the correction
     * already forces into agreement. */
    VIO_Volume field = n3::like(input);
    n3::evaluate_saved_field(imp_of(corr_path), field, mask);
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
    CHECK_TRUE("the field is strictly positive in the mask",
               nonfinite == 0 && nonpositive == 0);
    delete_volume(field);
  }

  /* ---- 2. in-mask CV of the output below the input's ---- */
  {
    n3::Stats in = n3::masked_stats(input, mask);
    n3::Stats out = n3::masked_stats(corr, mask);
    double cvin = in.stddev/in.mean, cvout = out.stddev/out.mean;
    CHECK_TRUE("in-mask CV of the output is below the input's", cvout < cvin);
    printf("  CV in %g -> CV out %g\n", cvin, cvout);
  }

  /* ---- 3. a no-non-uniformity volume -> a bounded, finite field ----
   * A two-tissue phantom (sharp boundary = structure, not shading, so the
   * anatomy contrast is a hard floor the estimation must stay below) with
   * additive noise: noise is additive anatomy, not the smooth multiplicative
   * term N3 removes, so there is no bias to estimate.  A noise-free or
   * single-valued volume is a degenerate histogram for the sharpen
   * deconvolution and diverges (the Wiener filter divides by an empty
   * spectrum); the noise keeps it estimatable. */
  {
    /* The phantom is indexed in values()'s storage order: contiguous along
     * sizes[2] (z), then sizes[1] (y), then sizes[0] (x) -- i = x*(ny*nz) +
     * y*nz + z -- so the fastest-varying counter walks z, not x, and rx/ry/rz
     * sit on their own axes.  The old `x = i % sizes[0]` marched across
     * rows and the region came out aliased stripes, not an ellipsoid. */
    int sizes[VIO_N_DIMENSIONS];
    get_volume_sizes(input, sizes);
    int nx = sizes[0], ny = sizes[1], nz = sizes[2];
    double cx = 0.5*nx, cy = 0.5*ny, cz = 0.5*nz;
    double rx = 0.36*nx, ry = 0.36*ny, rz = 0.36*nz;
    VIO_Volume ph = n3::like(input);
    double *pd = n3::values(ph);
    std::mt19937 rng(12345);
    std::normal_distribution<double> noise(0.0, 15.0);
    for(int i = 0; i < n; i++)
      {
        int z = i % nz, yz = i / nz, y = yz % ny, x = yz / ny;
        double dx = (x - cx)/rx, dy = (y - cy)/ry, dz = (z - cz)/rz;
        bool in = dx*dx + dy*dy + dz*dz <= 1.0;
        pd[i] = (in ? 250.0 : 100.0) + noise(rng);
      }
    VIO_BOOL sf; nc_type type = n3::storage_type(inp, &sf);
    std::string ph_path = tag_path("_ph.mnc");
    n3::save(ph, ph_path, inp, type, sf, "test_driver_properties");

    std::string out_path = tag_path("_ph_out.mnc");
    snprintf(cmd, sizeof(cmd),
             "\"%s\" -shrink 2 -iterations 15 -stop 0.0 -distance 100 "
             "-mask \"%s\" \"%s\" \"%s\" -clobber",
             N3_DRIVER_BIN, mask_in.c_str(), ph_path.c_str(), out_path.c_str());
    if(system(cmd) != 0)
      {
        printf("driver did not finish the correction on the phantom\n");
        n3check::failures()++;
        delete_volume(ph); cleanup(out_path); cleanup(ph_path);
        return n3check::report("driver_properties");
      }
    VIO_Volume out = n3::load(out_path.c_str());
    /* Compare against the phantom as stored on disk (reload it): the field is
     * phantom/output, and the phantom must be read at the same quantisation as
     * the driver saw it, or reloading-noise shows up as field structure. */
    VIO_Volume ph_disk = n3::load(ph_path.c_str());
    const double *od = n3::values(out), *pd_disk = n3::values(ph_disk);
    double sum = 0, sumsq = 0, cnt = 0;
    bool nan = false, inf = false;
    for(int i = 0; i < n; i++)
      if(mv[i] > 0.0)
        {
          double f = pd_disk[i]/od[i];  /* field = phantom / output */
          if(std::isnan(f)) nan = true;
          if(!std::isfinite(f)) { inf = true; continue; }
          sum += f; sumsq += f*f; cnt += 1.0;
        }
    double mean = sum/cnt;
    /* A mean-square coefficient of variation over the mask, not an extreme
     * value: a handful of quantisation-boundary voxels must not drive the
     * pass/fail the way lo/hi of a ratio would (CLAUDE.md, whole-volume
     * extremes).  The phantom is two-valued at 100/250, so a field that
     * absorbed the anatomy would sit at a CV of order the tissue spread; the
     * bound below is a quarter of that contrast, so the property is "the
     * estimator does not absorb the structure it is meant to be blind to",
     * not a tautological near-constant ideal no N3 reaches.  Measured field CV
     * is ~0.05 here (legacy ~0.0064), a port over-correction churned in cycle
     * 14, so the 0.3 derived bound leaves that gap clearly under it. */
    double rms_cv = sqrt(fabs(sumsq/cnt - mean*mean))/mean;
    double contrast = (250.0 - 100.0) / ((250.0 + 100.0) / 2.0);
    printf("  no-bias phantom field: mean %.4g  RMS CV %.4g  (tissue contrast %.4g)\n",
           mean, rms_cv, contrast);
    CHECK_TRUE("a no-nonuniformity volume yields a bounded, finite field",
               !nan && !inf && cnt > 0
               && std::isfinite(rms_cv) && rms_cv < 0.25 * contrast);
    delete_volume(ph_disk);
    delete_volume(out);
    delete_volume(ph);
    cleanup(out_path);
    cleanup(ph_path);
  }

  /* ---- 4. iteration counts at -stop 0 = the requested total ---- */
  {
    std::string outlog = tag_path("_est.log");
    std::string imp = tag_path("_est.imp");
    snprintf(cmd, sizeof(cmd),
             "\"%s\" -shrink 2 -iterations 6 -stop 0.0 -distance 100 "
             "-mask \"%s\" -estimate_only \"%s\" \"%s\" > \"%s\" 2>&1",
             N3_DRIVER_BIN, mask_in.c_str(), inp.c_str(), imp.c_str(), outlog.c_str());
    if(system(cmd) != 0)
      {
        printf("driver did not finish the estimate on chunk\n");
        n3check::failures()++;
        unlink(outlog.c_str()); unlink(imp.c_str());
        return n3check::report("driver_properties");
      }
    FILE *f = fopen(outlog.c_str(), "r");
    int got = -1;
    char line[256];
    while(f && fgets(line, sizeof(line), f))
      if(sscanf(line, "Number of iterations: %d", &got) == 1) break;
    if(f) fclose(f);
    /* total_iterations is the LAST staged count (Perl :1630-1631), not the
     * sum, and -stop 0 never stops a stage early, so an estimate requested at
     * -iterations 6 -stop 0 runs exactly 6 stages and reports 6.  This asserts
     * the driver's reading of the staged rule (cycle 10) end to end; it does
     * not claim agreement with the Perl, which the test never invokes. */
    CHECK_NEAR("-stop 0 runs all stages: iterations == the requested 6",
               (double) got, 6.0, 0.0);
    unlink(outlog.c_str()); unlink(imp.c_str());
  }

  delete_volume(corr);
  delete_volume(mask);
  delete_volume(input);
  cleanup(corr_path);
  return n3check::report("driver_properties");
}
