/* Cycle 0: the harness, before anything is trusted to it.
 *
 * Run with --fail it makes one deliberate failing check; CTest registers that
 * invocation with WILL_FAIL, so "a failing check exits non-zero" stays
 * asserted rather than being demonstrated once and deleted.
 */

#include "check.h"

#include <cstring>

int main(int argc, char **argv)
{
  bool want_failure = (argc > 1 && strcmp(argv[1], "--fail") == 0);

  if(want_failure)
    {
      CHECK_NEAR("a deliberate failure is reported as one", 1.0, 2.0, 1e-9);
      return n3check::report("harness (expected to fail)");
    }

  CHECK_TRUE("a true condition passes", true);
  CHECK_NEAR("equal values pass at any bound", 1.0, 1.0, 0.0);
  CHECK_NEAR("a difference exactly at the bound passes", 1.0, 1.5, 0.5);

  /* rel_rms is the comparison every volume-sized check uses; an error in it
   * would silently loosen all of them. */
  const double a[4] = { 1.0, 2.0, 3.0, 4.0 };
  const double b[4] = { 1.0, 2.0, 3.0, 4.0 };
  const double c[4] = { 1.1, 2.0, 3.0, 4.0 };
  CHECK_NEAR("rel_rms of identical arrays is zero",
             n3check::rel_rms(a, b, 4), 0.0, 0.0);
  CHECK_NEAR("rel_rms is sqrt(sum d^2 / sum b^2)",
             n3check::rel_rms(c, b, 4), 0.1 / sqrt(30.0), 1e-15);

  /* Scale invariance: relative RMS must not depend on the units of the
   * quantity compared, which is why it and not max|a-b| is the volume-level
   * comparison. */
  double sa[4], sb[4];
  for(int i = 0; i < 4; i++) { sa[i] = 1e6 * c[i]; sb[i] = 1e6 * b[i]; }
  CHECK_NEAR("rel_rms is invariant to a common scale",
             n3check::rel_rms(sa, sb, 4), n3check::rel_rms(c, b, 4), 1e-15);

  CHECK_TRUE("no failure was recorded", n3check::failures() == 0);
  return n3check::report("harness");
}
