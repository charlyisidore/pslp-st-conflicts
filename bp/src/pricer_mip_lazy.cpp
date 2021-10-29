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

#include "pricer_mip_lazy.hpp"
#include "conshdlr_cycle.hpp"
#include "graph_algorithm.hpp"
#include "logging.hpp"
#include "pprint.hpp"
#include <cassert>
#include <scip/cons_linear.h>
#include <scip/cons_varbound.h>
#include <sstream>

// Print initial 2-cycles
constexpr bool Print2Cycles = false;

// Print initial DFS cycles
constexpr bool PrintDfsCycles = false;

// Print the solution
constexpr bool PrintSolution = false;

// Print SCIP variables
constexpr bool PrintScipVariables = false;

using namespace pslp::bp;

const char *PricerMipLazy::Properties::Name = "mip_lazy";
const char *PricerMipLazy::Properties::Desc = "MIP variable pricer with lazy constraints";
const int PricerMipLazy::Properties::Priority = 0;
const SCIP_Bool PricerMipLazy::Properties::Delay = true;

// For stopping a callback execution if SCIP returns an error
struct SCIPError
{
    SCIPError(SCIP_RETCODE retcode)
        : retcode(retcode) {}
    SCIP_RETCODE retcode;
};

// Macro throwing a SCIPError in case of failure
#define SCIP_THROW(x)                     \
    {                                     \
        SCIP_RETCODE retcode;             \
        if ((retcode = (x)) != SCIP_OKAY) \
        {                                 \
            throw SCIPError(retcode);     \
        }                                 \
    }

PricerMipLazy::PricerMipLazy(SCIP *scip)
    : PricerMip(scip,
                Properties::Name,
                Properties::Desc,
                Properties::Priority,
                Properties::Delay)
{
    SCIPincludeObjConshdlr(_subscip, new ConshdlrCycle(_subscip), true);

    {
        std::ostringstream oss;
        oss << "pricers/" << Properties::Name << "/add2cyclecons";
        SCIPaddBoolParam(scip,
                         oss.str().c_str(),
                         "(pslp) add initial 2-cycle elimination constraints?",
                         &_add2cyclecons,
                         false,
                         _add2cyclecons,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << Properties::Name << "/adddfscyclecons";
        SCIPaddBoolParam(scip,
                         oss.str().c_str(),
                         "(pslp) add initial cycle elimination constraints found by DFS?",
                         &_adddfscyclecons,
                         false,
                         _adddfscyclecons,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << Properties::Name << "/maxcuts";
        SCIPaddIntParam(scip,
                        oss.str().c_str(),
                        "(pslp) maximal number of cuts added per round (0: unlimited)",
                        &_maxcuts,
                        false,
                        _maxcuts,
                        0,
                        INT_MAX,
                        nullptr,
                        nullptr);
    }
}

SCIP_RETCODE PricerMipLazy::init(SCIP *scip, const Problem &prob)
{
    if (prob.has_adjacent_conflicts)
    {
        LOG_ERROR("adjacent conflicts not supported");
        return SCIP_ERROR;
    }

    // Create the subscip
    SCIP_CALL(PricerMip::init(scip, prob));

    const auto n = prob.conflict_graph.n_vertices();

    // Initialize constraints on cycles of length 2
    if (_add2cyclecons)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            // Get out-neighbors (i -> j edges)
            for (auto j : prob.conflict_graph.adjacent_vertices(i))
            {
                // There is a 2-cycle only if there exists a j -> i edge
                // Avoid adding the constraint twice by imposing i > j
                if (i > j || !prob.conflict_matrix(j, i))
                {
                    continue;
                }

                if constexpr (Print2Cycles)
                {
                    LOG_TRACE("2-cycle: %zu %zu", i + 1, j + 1);
                }

                SCIP_CONS *cons;
                std::ostringstream oss;

                oss << "cycle{" << i + 1 << ',' << j + 1 << '}';

                SCIP_CALL(SCIPcreateConsBasicVarbound(_subscip,
                                                      &cons,
                                                      oss.str().c_str(),
                                                      _x[i],
                                                      _x[j],
                                                      1.,
                                                      -SCIPinfinity(_subscip),
                                                      1.));

                SCIP_CALL(SCIPaddCons(_subscip, cons));
                SCIP_CALL(SCIPreleaseCons(_subscip, &cons));
            }
        }
    }

    // Run DFS to add more cycles
    if (_adddfscyclecons)
    {
        // This callback creates a constraint for every cycle found
        auto create_constraint = [this, &prob](auto first, auto last) {
            std::size_t cycle_length = std::distance(first, last);

            assert(cycle_length >= 2);

            // Avoid creating redundant constraints
            if (cycle_length >= prob.capacity ||
                (_add2cyclecons && cycle_length == 2))
            {
                return;
            }

            assert(!_add2cyclecons || cycle_length > 2);

            if constexpr (PrintDfsCycles)
            {
                std::ostringstream oss;
                for (auto it = first; it != last; ++it)
                {
                    oss << ' ' << *it + 1;
                }
                LOG_TRACE("dfs-cycle:%s", oss.str().c_str());
            }

            SCIP_CONS *cons;

            SCIP_THROW(SCIPcreateConsBasicLinear(_subscip,
                                                 &cons,
                                                 "cycle",
                                                 0,
                                                 nullptr,
                                                 nullptr,
                                                 -SCIPinfinity(_subscip),
                                                 cycle_length - 1));

            for (auto i = first; i != last; ++i)
            {
                SCIP_THROW(SCIPaddCoefLinear(_subscip, cons, _x[*i], 1.));
            }

            SCIP_THROW(SCIPaddCons(_subscip, cons));
            SCIP_THROW(SCIPreleaseCons(_subscip, &cons));
        };

        // Find cycles
        try
        {
            CycleFinder cycle_finder(create_constraint);
            depth_first_visit(prob.conflict_graph, cycle_finder);
        }
        catch (const SCIPError &e)
        {
            return e.retcode;
        }
    }

    // Lazy constraints
    {
        SCIP_CONS *cons;

        SCIP_CALL(ConshdlrCycle::create_cons(_subscip,
                                             &cons,
                                             "cycle",
                                             _x,
                                             prob.conflict_matrix,
                                             _maxcuts));

        SCIP_CALL(SCIPaddCons(_subscip, cons));
        SCIP_CALL(SCIPreleaseCons(_subscip, &cons));
    }

    return SCIP_OKAY;
}

SCIP_RETCODE PricerMipLazy::solutions(SCIP *scip,
                                      const Problem &prob,
                                      std::vector<Solution> &sols)
{
    const auto n = prob.conflict_graph.n_vertices();

    SCIP_SOL **scip_sols = SCIPgetSols(_subscip);
    int n_sols = SCIPgetNSols(_subscip);

    for (int k = 0; k < n_sols; ++k)
    {
        SCIP_SOL *sol = scip_sols[k];
        double obj = prob.offset;
        std::vector<std::size_t> items;

        // We recompute the objective because of high numerical imprecision
        for (std::size_t i = 0; i < n; ++i)
        {
            if (SCIPgetSolVal(_subscip, sol, _x[i]) > 0.5)
            {
                obj += prob.coefs[i];
                items.push_back(i);
            }
        }

#ifndef NDEBUG
        double scip_obj = SCIPgetSolOrigObj(_subscip, sol) + prob.offset;
        if (SCIPisPositive(scip, obj) != SCIPisPositive(scip, scip_obj))
        {
            LOG_WARN("numerical imprecision (SCIP obj: %g [%s], recomputed: %g [%s])",
                     scip_obj,
                     SCIPisPositive(scip, scip_obj) ? "positive" : "non positive",
                     obj,
                     SCIPisPositive(scip, obj) ? "positive" : "non positive");
        }
#endif

        // Can this solution improve the master problem's objective?
        if (!SCIPisPositive(scip, obj))
        {
            // SCIP solutions are sorted by objective value
            break;
        }

        LOG_DEBUG("found MIP pricing sol (obj: %g)", obj);

        if constexpr (PrintScipVariables)
        {
            SCIP_CALL(SCIPprintSol(_subscip, sol, nullptr, false));
        }

        Solution s(std::size(items));

        // Build a subgraph from selected items
        auto subgraph = Pricer::subgraph(prob.conflict_matrix, items);

        // Find a topological sorting
        try
        {
            // Elements are in reversed order
            std::size_t pos = std::size(s);
            TopologicalSorter topo_sorter([&](auto v) {
                assert(pos > 0);
                s[--pos] = items[v];
            });
            depth_first_visit(subgraph, topo_sorter);
        }
        catch (const NotADag &)
        {
            LOG_ERROR("found cycle");
            return SCIP_ERROR;
        }

        if constexpr (PrintSolution)
        {
            LOG_TRACE("-> sol: %s",
                      pprint_n(std::size(s),
                               [&s](auto i) { return s[i] + 1; })
                          .c_str());
        }

        sols.push_back(s);
    }

    return SCIP_OKAY;
}