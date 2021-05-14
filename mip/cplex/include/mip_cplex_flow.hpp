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

#ifndef PSLP_MIP_MIP_CPLEX_FLOW_HPP
#define PSLP_MIP_MIP_CPLEX_FLOW_HPP

#include "mip_cplex.hpp"

namespace pslp
{

namespace mip
{

/**
 * Network flow model for the PSLP.
 */
class MipCplexFlow : public MipCplex
{
public:
    /**
     * Problem type definition.
     */
    using Problem = MipCplex::Problem;

    /**
     * Solution type definition.
     */
    using Solution = MipCplex::Solution;

    /**
     * Create the model.
     */
    void create(const Problem &prob) override;

    /**
     * Create an initial solution.
     */
    void initialize(const Problem &prob) override;

    /**
     * Check the solution.
     */
    Solution solution(const Problem &prob) override;

    /**
     * Write the solution and statistics.
     */
    void write(const Problem &prob,
               const Solution &sol,
               nlohmann::json &json) override;

protected:
    IloArray<IloNumVarArray> _x;
    IloNumVarArray _s;
    IloNumVarArray _t;
    IloNumVarArray _l;
};

} // namespace mip

} // namespace pslp

#endif