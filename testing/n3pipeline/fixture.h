/* Reading the recorded oracle answers.
 *
 * No test runs an N3 program: regenerate_reference.sh drives the installed
 * binaries once and writes this directory, and the tests read it back.  The
 * same policy is in force in the PyTorch tree (tests/regenerate_reference.py),
 * for the same reason -- a test that shells out measures the machine's PATH as
 * much as the code.
 *
 * Volumes are raw little-endian float64, smaller arrays are text.  Only what
 * an assertion looks at is recorded.
 */

#ifndef N3_FIXTURE_H
#define N3_FIXTURE_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>
#include <vector>

namespace n3fixture {

inline std::string dir()
{
  const char *env = getenv("N3_REFERENCE_DIR");
  return env ? std::string(env) : std::string(N3_REFERENCE_DIR);
}

inline std::string path(const std::string &name) { return dir() + "/" + name; }

/* A missing fixture is a harness failure, not a test failure: it means
 * regenerate_reference.sh was not run.  Say so and stop rather than reporting
 * a comparison against nothing. */
inline void must(bool ok, const std::string &what)
{
  if(!ok) { fprintf(stderr, "fixture: %s\n", what.c_str()); exit(2); }
}

inline std::vector<double> read_f64(const std::string &name)
{
  std::string p = path(name);
  FILE *f = fopen(p.c_str(), "rb");
  must(f != NULL, "cannot open " + p);
  fseek(f, 0, SEEK_END);
  long bytes = ftell(f);
  fseek(f, 0, SEEK_SET);
  must(bytes % 8 == 0, p + " is not a whole number of float64");
  std::vector<double> out(bytes / 8);
  if(!out.empty()) must(fread(&out[0], 8, out.size(), f) == out.size(),
                        "short read on " + p);
  fclose(f);
  return out;
}

/* Every whitespace-separated number in a text file, comments (#) skipped.
 * volume_hist's tables, sharpen_hist's lookup tables and volume_stats' output
 * are all read this way. */
inline std::vector<double> read_text(const std::string &name)
{
  std::string p = path(name);
  FILE *f = fopen(p.c_str(), "r");
  must(f != NULL, "cannot open " + p);
  std::vector<double> out;
  char line[4096];
  while(fgets(line, sizeof(line), f))
    {
      char *hash = strchr(line, '#');
      if(hash) *hash = '\0';
      char *p = line, *end;
      while(*p)
        {
          double v = strtod(p, &end);
          if(end == p) break;
          out.push_back(v);
          p = end;
        }
    }
  fclose(f);
  return out;
}

/* The "# domain: min max" line volume_hist writes. */
inline void read_domain(const std::string &name, double *min_bin, double *max_bin)
{
  std::string p = path(name);
  FILE *f = fopen(p.c_str(), "r");
  must(f != NULL, "cannot open " + p);
  char line[4096];
  bool found = false;
  while(fgets(line, sizeof(line), f))
    if(sscanf(line, "# domain: %lf %lf", min_bin, max_bin) == 2) { found = true; break; }
  fclose(f);
  must(found, "no domain line in " + p);
}

/* Volume-sized oracles are recorded on every fourth voxel in file order
 * (regenerate_reference.sh).  A relative RMS over a systematic quarter says
 * what one over the whole volume says, without ten megabytes of float64 in the
 * source tree; the tests apply the same stride to their own answer. */
const int STRIDE = 4;

inline std::vector<double> strided(const double *values, int n)
{
  std::vector<double> out;
  for(int i = 0; i < n; i += STRIDE) out.push_back(values[i]);
  return out;
}

inline double read_scalar(const std::string &name)
{
  std::vector<double> v = read_text(name);
  must(v.size() >= 1, "no number in " + name);
  return v[0];
}

/* The number of representable levels between a MINC file's valid_range
 * endpoints (chunk.mnc records "0 4095", regenerated as chunk_valid_range.txt).
 * The MINC-round-trip bounds are derived from this and from the data rather
 * than assuming 16 bits: a 12-bit file has 4095 steps, and quoting a fixed
 * denominator is wrong for any other. */
inline double valid_steps(const std::string &name)
{
  std::vector<double> r = read_text(name);
  must(r.size() == 2, name + " must hold the valid range");
  must(r[1] > r[0], name + " has an empty range");
  return r[1] - r[0];
}

/* ---- driving the compiled binary ------------------------------------------
 *
 * Only the driver tests define N3_DRIVER_BIN; the eleven block tests include
 * this header for the readers above alone and must not pay for a declaration
 * they cannot satisfy.
 *
 * These four are shared because each was previously copied per test file and
 * one copy was wrong: test_driver_endtoend's cleanup removed <out>.mnc.imp
 * where the driver writes <out>.imp, so every run leaked two .imp files into
 * TMPDIR while the test stayed green (2026-08-07 review, item 2; predicted by
 * the 2026-08-06 review's item 22).  PLAN §9: a rule implemented twice is a
 * defect.
 */
#ifdef N3_DRIVER_BIN

inline std::string outdir()
{
  const char *tmp = getenv("TMPDIR");
  return std::string(tmp ? tmp : "/tmp");
}

/* The .imp the driver writes beside a correct output <out>.mnc is <out>.imp:
 * nu_correct_cxx.cc's imp_path() replaces the FINAL extension
 * (MNI::PathUtilities::replace_ext = s/\.[^\.]*$/\.imp/), so it is never
 * <out>.mnc.imp, and out.mnc.gz gives out.mnc.imp. */
inline std::string imp_of(const std::string &path)
{
  size_t dot = path.find_last_of('.');
  return (dot == std::string::npos ? path : path.substr(0, dot)) + ".imp";
}

/* A per-pid path in TMPDIR, so concurrent ctest jobs do not collide. */
inline std::string temp_path(const std::string &stem)
{
  return outdir() + "/n3cxx_" + std::to_string((int) getpid()) + "_" + stem;
}

/* Run the driver: <bin> <opts> "<input>" "<output>" -clobber, with stdout and
 * stderr redirected to `log` when one is given.  Returns true if it exited
 * zero.  The command is composed with std::string rather than into a fixed
 * buffer: an snprintf whose return is discarded turns a long TMPDIR into a
 * different, truncated command and the test then reports a driver failure
 * that did not occur (2026-08-07 review, item 6; 2026-08-06's item 21). */
inline bool run_driver(const std::string &opts,
                       const std::string &input,
                       const std::string &output,
                       const std::string &log = std::string())
{
  std::string cmd = std::string("\"") + N3_DRIVER_BIN + "\" " + opts
    + " \"" + input + "\" \"" + output + "\" -clobber";
  if(!log.empty()) cmd += " > \"" + log + "\" 2>&1";
  return system(cmd.c_str()) == 0;
}

/* Remove a correct run's outputs: the volume, its .imp, and a <volume>.log if
 * the caller redirected one. */
inline void cleanup(const std::string &output)
{
  unlink(output.c_str());
  unlink(imp_of(output).c_str());
  unlink((output + ".log").c_str());
}

#endif  /* N3_DRIVER_BIN */

}  // namespace n3fixture

#endif
