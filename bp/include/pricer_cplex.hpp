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

#ifndef PSLP_BP_PRICER_CPLEX_HPP
#define PSLP_BP_PRICER_CPLEX_HPP

#include "pricer.hpp"
#include <ilcplex/ilocplex.h>
#include <vector>

namespace pslp
{

namespace bp
{

/**
 * Base class for CPLEX-based variable pricers.
 */
class PricerCplex : public Pricer
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
     * Maximum int value.
     */
    static constexpr int IntMax = 2100000000;

    /**
     * Maximum Longint value.
     */
    static constexpr SCIP_Longint LongintMax = 9223372036800000000;

    /**
     * Minimum Real value.
     */
    static constexpr SCIP_Real RealMin = -1.e75;

    /**
     * Maximum Real value.
     */
    static constexpr SCIP_Real RealMax = 1.e75;

    /**
     * Constructor.
     */
    PricerCplex(SCIP *scip,
                const char *name,
                const char *desc,
                int priority,
                SCIP_Bool delay);

    /**
     * Destructor.
     */
    virtual ~PricerCplex();

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
    [[nodiscard]] virtual SCIP_RETCODE solution(SCIP *scip,
                                                const Problem &prob,
                                                std::size_t k,
                                                std::vector<Solution> &sols) = 0;

protected:
    IloEnv _env;
    IloCplex _cplex;
    IloModel _model;
    IloObjective _objective;
    IloNumVarArray _x;
    IloRangeArray _conss_same_diff;

private:
    int _threads = 0;
    SCIP_Real _timelimit = RealMax;
    SCIP_Real _workmem = 2048;
    int _mip_strategy_file = 1;
    SCIP_Real _mip_limits_treememory = RealMax;
    SCIP_Longint _mip_limits_solutions = LongintMax;
    SCIP_Real _simplex_limits_upperobj = RealMax;
    SCIP_Bool _display_out = false;
    SCIP_Bool _display_warning = false;
    char *_lpfilename = nullptr;
    SCIP_Bool _extratime = false;
    SCIP_Real _global_time_limit = RealMax;
};

} // namespace bp

} // namespace pslp

#endif