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

#ifndef PSLP_BP_PRICER_GREEDY_HPP
#define PSLP_BP_PRICER_GREEDY_HPP

#include "pricer.hpp"
#include "symmetric_matrix.hpp"
#include <vector>

namespace pslp
{

namespace bp
{

/**
 * Greedy variable pricer for the PSLP.
 */
class PricerGreedy : public Pricer
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
    PricerGreedy(SCIP *scip);

    /**
     * Check whether the pricer finds exact solutions.
     */
    bool is_exact() const override
    {
        return false;
    }

    /**
     * Solve the pricing problem.
     */
    [[nodiscard]] SCIP_RETCODE solve(SCIP *scip,
                                     const Problem &prob,
                                     std::vector<Solution> &sols) override;

private:
    std::vector<std::vector<std::size_t>> _unions;
    std::vector<std::size_t> _union_of;
    SymmetricMatrix<bool, false> _diff_matrix;
};

} // namespace bp

} // namespace pslp

#endif