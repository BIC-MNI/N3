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

inline double read_scalar(const std::string &name)
{
  std::vector<double> v = read_text(name);
  must(v.size() >= 1, "no number in " + name);
  return v[0];
}

}  // namespace n3fixture

#endif
