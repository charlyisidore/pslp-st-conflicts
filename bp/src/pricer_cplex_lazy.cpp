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

#include "pricer_cplex_lazy.hpp"
#include "cons_callback_cycle.hpp"
#include "graph_algorithm.hpp"
#include "logging.hpp"
#include "pprint.hpp"
#include <cassert>
#include <iostream>
#include <sstream>

// Print initial 2-cycles
constexpr bool Print2Cycles = false;

// Print initial DFS cycles
constexpr bool PrintDfsCycles = false;

// Print the solution
constexpr bool PrintSolution = false;

// Print CPLEX variables
constexpr bool PrintCplexVariables = false;

using namespace pslp::bp;

const char *PricerCplexLazy::Properties::Name = "cplex_lazy";
const char *PricerCplexLazy::Properties::Desc = "CPLEX variable pricer with lazy constraints formulation";
const int PricerCplexLazy::Properties::Priority = 0;
const SCIP_Bool PricerCplexLazy::Properties::Delay = true;

PricerCplexLazy::PricerCplexLazy(SCIP *scip)
    : PricerCplex(scip,
                  Properties::Name,
                  Properties::Desc,
                  Properties::Priority,
                  Properties::Delay)
{
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

SCIP_RETCODE PricerCplexLazy::init(SCIP *scip, const Problem &prob)
{
    if (prob.has_adjacent_conflicts)
    {
        LOG_ERROR("adjacent conflicts not supported");
        return SCIP_ERROR;
    }

    SCIP_CALL(PricerCplex::init(scip, prob));

    // Removes one warning
    _cplex.setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);

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

                std::ostringstream oss;
                oss << "cycle{" << i + 1 << ',' << j + 1 << '}';

                IloRange cons(_env, -IloInfinity, 1., oss.str().c_str());
                cons.setLinearCoef(_x[i], 1.);
                cons.setLinearCoef(_x[j], 1.);
                _model.add(cons);
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

            IloRange cons(_env, -IloInfinity, cycle_length - 1, "cycle");
            for (auto i = first; i != last; ++i)
            {
                cons.setLinearCoef(_x[*i], 1.);
            }
            _model.add(cons);
        };

        // Find cycles
        try
        {
            CycleFinder cycle_finder(create_constraint);
            depth_first_visit(prob.conflict_graph, cycle_finder);
        }
        catch (...)
        {
            return SCIP_ERROR;
        }
    }

    // Adjacency matrix of the conflict graph
    IloArray<IloBoolArray> adjacency_matrix(_env, n);
    for (std::size_t i = 0; i < n; ++i)
    {
        adjacency_matrix[i] = IloBoolArray(_env, n);
        for (std::size_t j = 0; j < n; ++j)
        {
            adjacency_matrix[i][j] = prob.conflict_matrix(i, j);
        }
    }

    _cplex.use(ConsCallbackCycle::create(_env,
                                         _x,
                                         adjacency_matrix,
                                         _maxcuts));

    return SCIP_OKAY;
}

SCIP_RETCODE PricerCplexLazy::solution(SCIP *scip,
                                       const Problem &prob,
                                       std::size_t k,
                                       std::vector<Solution> &sols)
{
    const auto n = prob.conflict_graph.n_vertices();
    double obj = _cplex.getObjValue(k);

    // Can this solution improve the master problem's objective?
    if (!SCIPisPositive(scip, obj))
    {
        return SCIP_OKAY;
    }

    LOG_DEBUG("found CPLEX pricing sol (obj: %g)", obj);

    if constexpr (PrintCplexVariables)
    {
        for (IloNum i = 0; i < _x.getSize(); ++i)
        {
            double value = _cplex.getValue(_x[i], k);
            if (!SCIPisEQ(scip, value, 0.))
            {
                LOG_DEBUG("%s  %g", _x[i].getName(), value);
            }
        }
    }

    std::vector<std::size_t> items;

    // Get all selected items, sort them later
    for (std::size_t i = 0; i < n; ++i)
    {
        if (_cplex.getValue(_x[i], k) > 0.5)
        {
            items.push_back(i);
        }
    }

    Solution s(std::size(items));

    // Build a subgraph from selected items
    auto subgraph = Pricer::subgraph(prob.conflict_matrix, items);

    // Find a topological sorting
    try
    {
        std::size_t pos = std::size(items);
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

    return SCIP_OKAY;
}