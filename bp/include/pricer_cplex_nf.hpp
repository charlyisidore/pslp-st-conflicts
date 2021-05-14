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

#ifndef PSLP_BP_PRICER_CPLEX_NF_HPP
#define PSLP_BP_PRICER_CPLEX_NF_HPP

#include "pricer_cplex.hpp"
#include <ilcplex/ilocplex.h>
#include <vector>

namespace pslp
{

namespace bp
{

/**
 * CPLEX-based variable pricer with network flow formulation for the PSLP.
 */
class PricerCplexNf : public PricerCplex
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
    PricerCplexNf(SCIP *scip);

    /**
     * Check whether the pricer supports a given problem.
     */
    bool supports(const Problem &prob) const override
    {
        // Only adjacent conflicts
        return prob.has_adjacent_conflicts;
    }

    /**
     * Initialize the pricer.
     */
    [[nodiscard]] SCIP_RETCODE init(SCIP *scip, const Problem &prob) override;

    /**
     * Process a given solution.
     */
    [[nodiscard]] SCIP_RETCODE solution(SCIP *scip,
                                        const Problem &prob,
                                        std::size_t k,
                                        std::vector<Solution> &sols) override;

private:
    IloArray<IloNumVarArray> _y;
    IloNumVarArray _s;
    IloNumVarArray _t;
    IloNumVarArray _l;
};

} // namespace bp

} // namespace pslp

#endif