/* The staged stopping rule.
 *
 * -iterations and -stop take matching lists (nu_estimate_np_and_em.in:1477-1478
 * requires them to be the same length); the last iteration count is the total
 * (:1630-1631), and each stage's threshold applies once that stage's iteration
 * count has been reached.  The shipped protocol is one stage,
 * -iterations 50 -stop 0.001; the -V0.9 protocol is two.
 *
 * Pure logic, transcribed from :171-185, so it has no oracle and needs none.
 */

#ifndef N3_STOPPING_H
#define N3_STOPPING_H

#include <vector>

namespace n3 {

/* The total number of iterations: the last stage's count. */
int total_iterations(const std::vector<int> &stages);

/* Whether iteration `iter` (zero based) is the last one, given the change in
 * the field it produced.
 *
 * The first threshold applies from the start; each later one applies only
 * once the previous stage's iteration count has been reached.  With no
 * thresholds at all the rule never fires, and with a threshold of zero it
 * never fires either -- which is what makes a fixed-iteration comparison
 * against the Perl possible. */
bool should_stop(int iter, double change, const std::vector<int> &stages,
                 const std::vector<double> &thresholds);

}  // namespace n3

#endif
