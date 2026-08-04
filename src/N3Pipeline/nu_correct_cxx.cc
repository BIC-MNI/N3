/* nu_correct / nu_estimate, entirely in memory.
 *
 * The top-level driver of the Perl pipeline (nu_estimate.in), with the file
 * round trips the Perl is built out of removed: nu_estimate_np_and_em's loop
 * and nu_evaluate run here in one process on double buffers (NuEstimate.cc,
 * NuEvaluate.cc), and the two communicate through a Field in memory instead of
 * through an .imp on disk.
 *
 * The .imp is still written -- as nu_estimate it is the whole output, and as
 * nu_correct it is the leftover the Perl leaves behind -- but it is a record
 * of the fit, not a hand-off between stages.
 *
 * Scope, from /app/PLAN.md: the np method with b_spline and tp_spline, and
 * the -sharpen/-parzen/-parzen_sigma/-shrink/-mask/-distance/-lambda/
 * -iterations/-stop/-normalize_field/-auto_mask/-bimodalT/-floor/-subsample/
 * -mapping_dir options.  Anything outside it -- the EM and WM branches, fir
 * smoothing, -real, -differential, -initial, -islands and their allies -- must
 * fail loudly rather than be silently ignored.
 */

#include <config.h>

#include <cstdarg>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <volume_io.h>

#include "Buffers.h"
#include "FitField.h"
#include "NuEstimate.h"
#include "NuEvaluate.h"

namespace {

struct Arguments
{
  std::string input, output;
  bool estimate_only = false;

  std::string mask;
  std::string mapping_dir;

  double distance = 200.0;      /* -distance: the field knot spacing, mm */
  double lambda = 1e-7;
  int subsample = 1;
  enum spline_type type = b_spline;

  bool sharpen = true;          /* sharpening on by default (V1.0) */
  double fwhm = 0.15, noise = 0.01; /* -sharpen/-fwhm width, and its noise */
  bool window = true;           /* -parzen, on by default */
  double parzen_sigma = 0.0;    /* 0 = off */
  int bins = 200;
  bool blur = false;

  double shrink = 4.0;

  std::vector<int> iterations;      /* {50} by default */
  std::vector<double> stop;          /* {0.001} by default */

  bool normalize_field = false;
  bool auto_mask = true;             /* nu_estimate.in adds -auto_mask always */
  bool bimodalT = false;
  double background = 1.0;
  bool have_background = false;

  double floor = 0.1;
  bool clobber = false;
  bool verbose = false;

  Arguments()
  {
    iterations.push_back(50);
    stop.push_back(0.001);
  }
};

[[noreturn]] void die(const char *fmt, ...)
{
  fprintf(stderr, "nu_correct_cxx: ");
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fprintf(stderr, "\n");
  exit(1);
}

void usage()
{
  fprintf(stderr,
    "Usage: nu_correct_cxx [options] input.mnc output.mnc\n"
    "   or: nu_estimate_cxx [options] input.mnc output.imp\n"
    "\n"
    "  -mask <file>       region to process\n"
    "  -distance <mm>     field knot spacing (default 200)\n"
    "  -lambda <x>        spline regularisation (default 1e-7)\n"
    "  -b_spline [-x]     tensor B-splines (default); optional lambda\n"
    "  -tp_spline [-x]    thin-plate splines; optional lambda\n"
    "  -subsample <n>     fit every n-th voxel each axis (also -spline_subsample)\n"
    "  -sharpen <f> <n>   histogram deconvolution (default 0.15 0.01)\n"
    "  -fwhm <x>          sharpening width; same as -sharpen <x> <noise>\n"
    "  -parzen            triangular window (default); -noparzen disables\n"
    "  -parzen_sigma <x>  Gaussian Parzen window, in bin widths\n"
    "  -bins <n>          histogram bins (default 200)\n"
    "  -nodeblur          skip the deconvolution\n"
    "  -shrink <factor>   estimation grid (default 4)\n"
    "  -iterations <n...> staged iteration counts (default 50)\n"
    "  -stop <x...>       staged stopping thresholds (default 0.001)\n"
    "  -normalize_field   scale the field to mean 1 in the mask\n"
    "  -auto_mask         automatic background mask (default); -bimodalT\n"
    "  -background <t>    background threshold (default 1)\n"
    "  -floor <x>         field floor, applied only if needed (default 0.1)\n"
    "  -mapping_dir <dir>  where to write the .imp\n"
    "  -estimate_only     write only the .imp; -correct overrides\n"
    "  -clobber           overwrite outputs\n"
    "  -verbose           be noisy\n"
    "  -V0.9              original protocol: shrink 3, stop 0.001 0.005,\n"
    "                     iterations 10 20\n");
}

/* The .imp path for a correct run: the output's directory and basename with a
 * .imp extension instead of the volume's, relocated into mapping_dir if that
 * was given (nu_estimate.in:57-58, replace_ext/replace_dir). */
std::string imp_path(const Arguments &args)
{
  /* The .imp's name is the output's basename with its volume extension
   * replaced (nu_estimate.in:57-58, replace_ext), under mapping_dir if one was
   * given and in the output's directory otherwise (replace_dir). */
  std::string out = args.output;
  size_t dot = out.find_last_of('.');
  if(dot != std::string::npos &&
     (out.compare(dot, 12, ".mnc.gz") == 0 ||
      out.compare(dot, 10, ".mnc.Z") == 0 ||
      out.compare(dot, 6, ".mnc") == 0))
    out = out.substr(0, dot);
  out += ".imp";

  if(!args.mapping_dir.empty())
    {
      size_t slash = out.find_last_of('/');
      out = args.mapping_dir + "/" +
        (slash == std::string::npos ? out : out.substr(slash + 1));
    }
  return out;
}

/* Whether a token is a number, so the multi-value options can tell their own
 * values from a following option or the positional arguments. */
bool is_number(const std::string &s)
{
  if(s.empty()) return false;
  size_t k = 0, n = s.size(), digits = 0;
  if(s[k] == '+' || s[k] == '-') k++;
  while(k < n)
    {
      if(s[k] >= '0' && s[k] <= '9') { digits++; k++; }
      else if((s[k] == '.' || s[k] == 'e' || s[k] == 'E') && k + 1 < n)
        { k++; if(s[k] == '+' || s[k] == '-') k++; }
      else break;
    }
  return digits > 0 && k == n;
}

struct Parsed
{
  Arguments args;
};

/* The options that are out of the ported scope.  They must fail, not be
 * ignored: silently dropping -em would run the wrong algorithm. */
bool out_of_scope(const std::string &tok)
{
  static const char *names[] = {
    "-em", "-expectation_maximization", "-white_matter", "-fir",
    "-real", "-differential", "-initial", "-islands", "-mean", "-tag",
    "-variance", "-probability", "-reduce", "-scale", "-workspace",
    "-save_fields", "-save_histograms", "-debug", "-mapping",
  };
  for(size_t i = 0; i < sizeof(names)/sizeof(names[0]); i++)
    if(tok.compare(names[i]) == 0)
      {
        fprintf(stderr, "nu_correct_cxx: %s is outside the ported scope\n",
                names[i]);
        exit(1);
      }
  return false;
}

}  // namespace

int main(int argc, char *argv[])
{
  /* Which program was asked for: an installed name of nu_estimate* means
   * estimate-only, as it does for the Perl (nu_estimate.in:420). */
  Arguments A;
  {
    std::string name = argv[0];
    size_t slash = name.find_last_of('/');
    if(slash != std::string::npos) name = name.substr(slash + 1);
    A.estimate_only = (name.find("nu_estimate") != std::string::npos);
  }

  std::vector<std::string> pos;
  for(int i = 1; i < argc; i++)
    {
      std::string tok = argv[i];
      if(tok.empty()) continue;

      if(tok[0] == '-' && tok.size() > 1)
        {
          if(out_of_scope(tok)) continue;

          if(tok == "-mask") { if(i+1>=argc) die("-mask needs a value"); A.mask = argv[++i]; }
          else if(tok == "-distance") { if(i+1>=argc) die("-distance needs a value"); A.distance = atof(argv[++i]); }
          else if(tok == "-lambda") { if(i+1>=argc) die("-lambda needs a value"); A.lambda = atof(argv[++i]); }
          else if(tok == "-b_spline") { A.type = b_spline; if(i+1<argc && is_number(argv[i+1])) A.lambda = atof(argv[++i]); }
          else if(tok == "-tp_spline") { A.type = thin_plate_spline; if(i+1<argc && is_number(argv[i+1])) A.lambda = atof(argv[++i]); }
          else if(tok == "-subsample" || tok == "-spline_subsample") { if(i+1>=argc) die("-subsample needs a value"); A.subsample = atoi(argv[++i]); }
          else if(tok == "-sharpen") {
            A.sharpen = true;
            if(i+1<argc && is_number(argv[i+1]))
              { A.fwhm = atof(argv[++i]);
                if(i+1<argc && is_number(argv[i+1])) A.noise = atof(argv[++i]); }
          }
          else if(tok == "-fwhm") {
            /* The sharpening width, not the knot spacing.  The top-level
             * nu_estimate.in:205 maps -fwhm to $user_options{'sharpen'} and
             * passes it on as -sharpen <width> 0.01 (:498); -distance is the
             * knots (:164-165).  Only the inner nu_estimate_np_and_em.in:1163
             * calls the sharpen width -fwhm. */
            if(i+1>=argc) die("-fwhm needs a value");
            A.fwhm = atof(argv[++i]);
          }
          else if(tok == "-parzen") { A.window = true; }
          else if(tok == "-noparzen") { A.window = false; }
          else if(tok == "-parzen_sigma") { if(i+1>=argc) die("-parzen_sigma needs a value"); A.parzen_sigma = atof(argv[++i]); }
          else if(tok == "-bins") { if(i+1>=argc) die("-bins needs a value"); A.bins = atoi(argv[++i]); }
          else if(tok == "-nodeblur" || tok == "-blur") { A.blur = true; }
          else if(tok == "-shrink") { if(i+1>=argc) die("-shrink needs a value"); A.shrink = atof(argv[++i]); }
          else if(tok == "-iterations") {
            A.iterations.clear();
            while(i+1<argc && is_number(argv[i+1]) && strchr(argv[i+1], '.') == NULL)
              A.iterations.push_back(atoi(argv[++i]));
            if(A.iterations.empty()) die("-iterations needs values");
          }
          else if(tok == "-stop") {
            A.stop.clear();
            while(i+1<argc && is_number(argv[i+1]))
              A.stop.push_back(atof(argv[++i]));
            if(A.stop.empty()) die("-stop needs values");
          }
          else if(tok == "-normalize_field") { A.normalize_field = true; }
          else if(tok == "-auto_mask") { A.auto_mask = true; }
          else if(tok == "-bimodalT") { A.bimodalT = true; }
          else if(tok == "-background") { if(i+1>=argc) die("-background needs a value"); A.background = atof(argv[++i]); A.have_background = true; }
          else if(tok == "-floor") { if(i+1>=argc) die("-floor needs a value"); A.floor = atof(argv[++i]); }
          else if(tok == "-mapping_dir") { if(i+1>=argc) die("-mapping_dir needs a value"); A.mapping_dir = argv[++i]; }
          else if(tok == "-estimate_only") { A.estimate_only = true; }
          else if(tok == "-correct") { A.estimate_only = false; }
          else if(tok == "-clobber") { A.clobber = true; }
          else if(tok == "-noclobber") { A.clobber = false; }
          else if(tok == "-verbose") { A.verbose = true; }
          else if(tok == "-quiet") { A.verbose = false; }
          else if(tok == "-V0.9") {
            A.iterations.assign({10, 20});
            A.stop.assign({0.001, 0.005});
            A.shrink = 3.0;
          }
          else if(tok == "-V1.0") { /* defaults already are 1.0 */ }
          else if(tok == "-help" || tok == "-h") { usage(); return 0; }
          else if(tok == "-version") { printf("nu_correct_cxx\n"); return 0; }
          else die("unknown option %s", tok.c_str());
        }
      else pos.push_back(tok);
    }

  if(A.iterations.size() != A.stop.size())
    die("-iterations and -stop must match in number of arguments");
  if(!A.sharpen && A.parzen_sigma > 0)
    die("-parzen_sigma requires -sharpen");
  if(A.lambda <= 0) die("smoothing parameter must be positive");
  if(A.subsample <= 0) die("subsampling factor must be positive");

  if(pos.size() != 2)
    {
      fprintf(stderr, "%s: expected an input and an output\n", argv[0]);
      usage();
      return 2;
    }
  A.input  = pos[0];
  A.output = pos[1];

  if(!A.clobber)
    {
      FILE *f = fopen(A.output.c_str(), "r");
      if(f) { fclose(f); die("output %s exists; use -clobber", A.output.c_str()); }
    }

  /* Load the input and, if given, the user mask (both at full resolution). */
  VIO_Volume input = n3::load(A.input);
  VIO_Volume user_mask = A.mask.empty() ? NULL : n3::load(A.mask);
  if(user_mask && n3::voxel_count(user_mask) != n3::voxel_count(input))
    die("mask and input have different sizes");

  /* Assemble the estimation options.  With no explicit threshold, -background
   * defaults to 1 (nu_estimate_np_and_em.in:1564); -auto_mask, which is on by
   * default, takes its threshold from volume_stats' bimodal rule when there is
   * no user mask. */
  n3::EstimateOptions e;
  e.shrink = A.shrink;
  e.distance = A.distance;
  e.lambda = A.lambda;
  e.subsample = A.subsample;
  e.type = A.type;
  e.bins = A.bins;
  e.fwhm = A.fwhm;
  e.noise = A.noise;
  e.window = A.window;
  e.parzen_sigma = A.parzen_sigma;
  e.blur = A.blur;
  e.background_threshold = A.background;
  e.bimodalT = A.bimodalT || (A.auto_mask && user_mask == NULL && !A.have_background);
  e.iterations = A.iterations;
  e.stop = A.stop;
  e.normalize_field = A.normalize_field;
  e.verbose = A.verbose;

  /* volume_stats sizes its histogram from the file's voxel range
   * (volumeStats.cc:286); in memory the file is gone, so the count has to be
   * carried here. */
  if(e.bimodalT)
    {
      double lo, hi;
      n3::voxel_range(A.input, &lo, &hi);
      e.bimodal_bins = (int) ceil(hi - lo + 1);
    }

  /* The estimate.  When writing a .imp is wanted -- every estimate-only run,
   * and every correct run for a record like the Perl leaves -- it is written
   * here, from the estimation grid, so its world domain comes out right. */
  int iterations_run = 0;
  double final_change = 0.0;
  std::string imp;
  if(A.estimate_only) imp = A.output;   /* the output is the .imp itself */
  else if(!A.mapping_dir.empty()) imp = imp_path(A); /* record like Perl leaves */

  n3::Field *field = n3::nu_estimate(input, user_mask, e, &iterations_run,
                                      &final_change, NULL,
                                      imp.empty() ? NULL : &imp);

  if(A.estimate_only)
    printf("Number of iterations: %d\nCV of field change: %g\n",
           iterations_run, final_change);

  if(!A.estimate_only)
    {
      n3::EvaluateOptions ev;
      ev.field_floor = A.floor;
      ev.verbose = A.verbose;

      VIO_Volume corrected = n3::nu_evaluate(input, user_mask, field, ev, NULL);

      VIO_BOOL signed_flag;
      nc_type type = n3::storage_type(A.input, &signed_flag);
      std::string history = "nu_correct_cxx " + std::string(argv[0]);
      n3::save(corrected, A.output, A.input, type, signed_flag, history);

      delete_volume(corrected);
    }

  delete field;
  if(user_mask) delete_volume(user_mask);
  delete_volume(input);
  return 0;
}
