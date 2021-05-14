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

#ifndef PSLP_BP_PRICER_MIP_HPP
#define PSLP_BP_PRICER_MIP_HPP

#include "pricer.hpp"
#include <scip/scip.h>
#include <vector>

namespace pslp
{

namespace bp
{

/**
 * Base class for MIP-based variable pricers.
 */
class PricerMip : public Pricer
{
public:
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
    PricerMip(SCIP *scip,
              const char *name,
              const char *desc,
              int priority,
              SCIP_Bool delay);

    /**
     * Destructor.
     */
    virtual ~PricerMip();

    /**
     * Check whether the pricer finds exact solutions.
     */
    bool is_exact() const override
    {
        return true;
    }

    /**
     * Initialize the pricer.
     */
    [[nodiscard]] SCIP_RETCODE init(SCIP *scip, const Problem &prob) override;

    /**
     * Deinitialize the pricer.
     */
    [[nodiscard]] SCIP_RETCODE exit(SCIP *scip) override;

    /**
     * Solve the pricing problem.
     */
    [[nodiscard]] SCIP_RETCODE solve(SCIP *scip,
                                     const Problem &prob,
                                     std::vector<Solution> &sols) override;

    /**
     * Called when branching constraints have changed.
     */
    [[nodiscard]] virtual SCIP_RETCODE change(SCIP *scip, const Problem &prob);

    /**
     * Process a given solution.
     */
    [[nodiscard]] virtual SCIP_RETCODE solutions(SCIP *scip,
                                                 const Problem &prob,
                                                 std::vector<Solution> &sols) = 0;

protected:
    SCIP *_subscip = nullptr;
    std::vector<SCIP_VAR *> _x;
    std::vector<SCIP_CONS *> _conss_same_diff;

private:
    /**
     * Create a new model.
     */
    [[nodiscard]] SCIP_RETCODE _create(SCIP *scip, const Problem &prob);

    SCIP_Real _limits_time = 1.e20;
    SCIP_Real _limits_memory = SCIP_MEM_NOLIMIT;
    int _display_verblevel = 0;
    SCIP_Bool _misc_catchctrlc = false;
    char *_lpfilename = nullptr;
    SCIP_Bool _extratime = false;
    SCIP_Real _global_time_limit = 1.e20;
};

} // namespace bp

} // namespace pslp

#endif