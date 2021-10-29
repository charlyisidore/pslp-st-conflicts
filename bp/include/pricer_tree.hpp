/* -*-c++-*-
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef PSLP_BP_PRICER_TREE_HPP
#define PSLP_BP_PRICER_TREE_HPP

#include "graph.hpp"
#include "pricer.hpp"
#include "symmetric_matrix.hpp"
#include <vector>

namespace pslp
{

namespace bp
{

/**
 * Tree search variable pricer for the PSLP.
 */
class PricerTree : public Pricer
{
public:
    /**
     * Variable pricer properties.
     */
    struct Properties
    {
        /**
         * Name of variable pricer.
         */
        static const char *Name;

        /**
         * Description of variable pricer.
         */
        static const char *Desc;

        /**
         * Priority of variable pricer.
         */
        static const int Priority;

        /**
         * Delay the variable pricer?
         */
        static const SCIP_Bool Delay;
    };

    /**
     * Pricer problem type definition.
     */
    using Problem = Pricer::Problem;

    /**
     * Pricer solution type definition.
     */
    using Solution = Pricer::Solution;

    /**
     * Constructor.
     */
    PricerTree(SCIP *scip);

    /**
     * Check whether the pricer finds exact solutions.
     */
    bool is_exact() const override
    {
        return false;
    }

    /**
     * Initialize the pricer.
     */
    [[nodiscard]] SCIP_RETCODE init(SCIP *scip, const Problem &prob) override;

    /**
     * Solve the pricing problem.
     */
    [[nodiscard]] SCIP_RETCODE solve(SCIP *scip,
                                     const Problem &prob,
                                     std::vector<Solution> &sols) override;

private:
    /**
     * Queue ID definition.
     */
    enum class QueueID
    {
        BFS,
        BrFS,
        CBFS,
        DFS
    };

    std::vector<std::vector<std::size_t>> _unions;
    std::vector<std::size_t> _union_of;
    SymmetricMatrix<bool, true> _diff_matrix;
    UndirectedGraph _diff_graph;
    QueueID _queue_id = QueueID::BFS;
    char *_nodeselection = nullptr;
    SCIP_Real _timelimit = 1.e20;
    SCIP_Bool _extratime = false;
    SCIP_Real _global_time_limit = 1.e20;
};

} // namespace bp

} // namespace pslp

#endif