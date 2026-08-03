/* Minimal assertion harness for the N3Pipeline tests.
 *
 * There is no gtest or catch2 in this environment and nothing may be
 * installed, so the tests are plain executables registered with CTest: each
 * check prints its measured margin against its bound, and report() returns the
 * process exit code.  Printing the margin is the point -- a green test that
 * says nothing about how much room it had is not evidence that anything is
 * still true tomorrow.
 */

#ifndef N3_CHECK_H
#define N3_CHECK_H

#include <cmath>
#include <cstdio>

namespace n3check {

inline int &failures() { static int n = 0; return n; }
inline int &total()    { static int n = 0; return n; }

inline void record(const char *what, bool ok, double margin, double bound)
{
  total()++;
  if(!ok) failures()++;
  printf("%-58s %s  margin %10.3e  bound %10.3e\n",
         what, ok ? "ok  " : "FAIL", margin, bound);
  fflush(stdout);
}

inline void record(const char *what, bool ok)
{
  total()++;
  if(!ok) failures()++;
  printf("%-58s %s\n", what, ok ? "ok  " : "FAIL");
  fflush(stdout);
}

inline void near(const char *what, double a, double b, double tol)
{
  double margin = fabs(a - b);
  record(what, margin <= tol, margin, tol);
}

/* Relative RMS difference, which is what a comparison over a whole volume
 * supports; max |a-b| over a volume is decided by a handful of mask-edge
 * voxels and does not support a fixed bound. */
inline double rel_rms(const double *a, const double *b, int n)
{
  double sd = 0.0, sb = 0.0;
  for(int i = 0; i < n; i++) { double d = a[i] - b[i]; sd += d*d; sb += b[i]*b[i]; }
  if(sb == 0.0) return sd == 0.0 ? 0.0 : 1.0 / 0.0;
  return sqrt(sd / sb);
}

inline void rms(const char *what, const double *a, const double *b, int n, double tol)
{
  double margin = rel_rms(a, b, n);
  record(what, margin <= tol, margin, tol);
}

/* Elementwise, for the comparisons whose bound is a per-value quantity (the
 * six decimals of %lf, or exactness). */
inline void all_near(const char *what, const double *a, const double *b, int n,
                     double tol)
{
  double margin = 0.0;
  for(int i = 0; i < n; i++) { double d = fabs(a[i] - b[i]); if(d > margin) margin = d; }
  record(what, margin <= tol, margin, tol);
}

inline int report(const char *name)
{
  printf("%s: %d checks, %d failed\n", name, total(), failures());
  return failures() == 0 ? 0 : 1;
}

}  // namespace n3check

#define CHECK_TRUE(what, cond)          n3check::record((what), (cond))
#define CHECK_NEAR(what, a, b, tol)     n3check::near((what), (a), (b), (tol))
#define CHECK_RMS(what, a, b, n, tol)   n3check::rms((what), (a), (b), (n), (tol))
#define CHECK_ALL(what, a, b, n, tol)   n3check::all_near((what), (a), (b), (n), (tol))

#endif
