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

#ifndef PSLP_BP_PRICER_MIP_LAZY_HPP
#define PSLP_BP_PRICER_MIP_LAZY_HPP

#include "pricer_mip.hpp"
#include <vector>

namespace pslp
{

namespace bp
{

/**
 * MIP-based variable pricer with lazy constraints formulation for the PSLP.
 */
class PricerMipLazy : public PricerMip
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
    PricerMipLazy(SCIP *scip);

    /**
     * Check whether the pricer supports a given problem.
     */
    bool supports(const Problem &prob) const override
    {
        // Only distant conflicts
        return !prob.has_adjacent_conflicts;
    }

    /**
     * Initialize the pricer.
     */
    [[nodiscard]] SCIP_RETCODE init(SCIP *scip,
                                    const Problem &prob) override;

    /**
     * Process a given solution.
     */
    [[nodiscard]] SCIP_RETCODE solutions(SCIP *scip,
                                         const Problem &prob,
                                         std::vector<Solution> &sols) override;

private:
    SCIP_Bool _add2cyclecons = true;
    SCIP_Bool _adddfscyclecons = false;
    int _maxcuts = 0;
};

} // namespace bp

} // namespace pslp

#endif