#include "Stopping.h"

#include <cstdio>
#include <cstdlib>

namespace n3 {

int total_iterations(const std::vector<int> &stages)
{
  if(stages.empty())
    {
      fprintf(stderr, "n3::total_iterations: no iteration count\n");
      exit(1);
    }
  return stages[stages.size() - 1];   /* :1630-1631 */
}

bool should_stop(int iter, double change, const std::vector<int> &stages,
                 const std::vector<double> &thresholds)
{
  if(thresholds.empty()) return false;          /* :171 */

  if(change < thresholds[0]) return true;       /* :177-178 */

  /* Each later stage's threshold applies only once the previous stage's
   * iteration count has been reached (:179-184).  The stages are tested
   * independently, so a later threshold can be the looser one and still
   * fire. */
  for(size_t stage = 1; stage < stages.size(); stage++)
    {
      if(stage >= thresholds.size()) break;
      if(iter >= stages[stage - 1] && change < thresholds[stage]) return true;
    }

  return false;
}

}  // namespace n3
