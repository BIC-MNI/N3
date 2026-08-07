/* Cycle 10: the staged stopping rule.
 *
 * Pure logic, so it is asserted against a table of cases rather than against
 * an oracle -- and it deserves its own cycle because the rule quantises
 * everything downstream.  Two implementations whose per-iteration change
 * differs in the fifth decimal can fall on either side of it and run different
 * numbers of iterations, which moves the output by far more than any block
 * difference does.
 *
 * The case that matters most for the rest of this work is the last one: with
 * -stop 0 the rule never fires, so the two pipelines can be made to run the
 * same number of iterations and compared.
 */

#include "check.h"

#include "../../src/N3Pipeline/Stopping.h"

#include <string>
#include <vector>

static std::vector<int> stages(int a, int b = -1)
{
  std::vector<int> v;
  v.push_back(a);
  if(b >= 0) v.push_back(b);
  return v;
}

static std::vector<double> thresholds(double a, double b = -1.0)
{
  std::vector<double> v;
  v.push_back(a);
  if(b >= 0.0) v.push_back(b);
  return v;
}

int main()
{
  /* The shipped protocol: -iterations 50 -stop 0.001. */
  {
    std::vector<int> s = stages(50);
    std::vector<double> t = thresholds(0.001);

    CHECK_TRUE("50 iterations in total", n3::total_iterations(s) == 50);
    CHECK_TRUE("a change below the threshold stops",
               n3::should_stop(0, 0.0009, s, t));
    CHECK_TRUE("a change above it does not",
               !n3::should_stop(0, 0.0011, s, t));
    CHECK_TRUE("a change exactly at it does not",
               !n3::should_stop(0, 0.001, s, t));
    CHECK_TRUE("the first threshold applies from the first iteration",
               n3::should_stop(0, 0.0001, s, t));
  }

  /* A *stricter* second threshold: -iterations 10 20 -stop 0.01 0.001.
   *
   * The first threshold applies from the start and the second only once the
   * first stage's iteration count has been reached -- but the two are tested
   * independently, so a stricter second threshold can never change the
   * outcome: anything below 0.001 is already below 0.01 and has stopped the
   * run.  Under this ordering the staged rule is the first threshold alone.
   * Recorded because it reads as though it tightens over time and does not.
   *
   * This is NOT the -V0.9 ordering, which an earlier version of this comment
   * claimed it was: V0.9 is -stop 0.001 0.005 (nu_estimate.in:438,
   * 'np:stop:0.9' => '0.001 0.005'), the looser-second case in the next
   * block, where the staged rule is not the first threshold alone. */
  {
    std::vector<int> s = stages(10, 20);
    std::vector<double> t = thresholds(0.01, 0.001);

    CHECK_TRUE("20 iterations in total", n3::total_iterations(s) == 20);
    CHECK_TRUE("the first stage's threshold applies from the start",
               n3::should_stop(0, 0.009, s, t));
    CHECK_TRUE("a change above both stops at neither",
               !n3::should_stop(15, 0.02, s, t));

    std::vector<double> alone = thresholds(0.01);
    std::vector<int> one = stages(20);
    bool same = true;
    for(int iter = 0; iter < 20; iter++)
      {
        double change = 0.05;
        for(int step = 0; step < 6; step++, change /= 10.0)
          if(n3::should_stop(iter, change, s, t)
             != n3::should_stop(iter, change, one, alone))
            same = false;
      }
    CHECK_TRUE("a stricter second threshold never changes the outcome", same);
  }

  /* A *looser* second threshold does change it: the stages are tested
   * independently, so the later one fires once its iteration count is
   * reached and not before.
   *
   * This is the -V0.9 ordering (-iterations 10 20 -stop 0.001 0.005,
   * nu_estimate.in:436-441), so on the one protocol that ships a second
   * threshold the staged rule really does relax after the first stage. */
  {
    std::vector<int> s = stages(10, 20);
    std::vector<double> t = thresholds(0.001, 0.01);
    CHECK_TRUE("a later, looser threshold fires once reached",
               n3::should_stop(12, 0.009, s, t));
    CHECK_TRUE("but not before its stage",
               !n3::should_stop(5, 0.009, s, t));
    CHECK_TRUE("and the strict one still fires at any time",
               n3::should_stop(0, 0.0009, s, t));
  }

  /* No thresholds: -iterations without -stop runs the count out. */
  {
    std::vector<int> s = stages(3);
    std::vector<double> t;
    CHECK_TRUE("without -stop the rule never fires",
               !n3::should_stop(0, 0.0, s, t) && !n3::should_stop(2, 0.0, s, t));
  }

  /* -stop 0 never fires either, however small the change gets.  This is what
   * lets the two pipelines be held to the same iteration count. */
  {
    std::vector<int> s = stages(3);
    std::vector<double> t = thresholds(0.0);
    CHECK_TRUE("-stop 0 never fires",
               !n3::should_stop(0, 0.0, s, t)
               && !n3::should_stop(1, 1e-300, s, t));
  }

  return n3check::report("stopping");
}
