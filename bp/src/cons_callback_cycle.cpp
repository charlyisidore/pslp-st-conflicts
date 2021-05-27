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

#include "cons_callback_cycle.hpp"
#include "graph.hpp"
#include "graph_algorithm.hpp"
#include <cassert>
#include <iostream>

// Print cuts
constexpr bool PrintCuts = false;

using namespace pslp::bp;

ConsCallbackCycle::ConsCallbackCycle(IloEnv env,
                                     const IloNumVarArray &x,
                                     const IloArray<IloBoolArray> &adjacency_matrix,
                                     std::size_t max_cuts)
    : IloCplex::LazyConstraintCallbackI(env),
      _x(x),
      _adjacency_matrix(adjacency_matrix),
      _max_cuts(max_cuts)
{
}

void ConsCallbackCycle::main()
{
    IloEnv env = getEnv();

    // Get the current solution
    std::vector<std::size_t> sol;

    for (IloInt i = 0; i < _x.getSize(); ++i)
    {
        const auto value = getValue(_x[i]);
        if (value > 0.5)
        {
            assert(value >= 0.9999);
            sol.push_back(i);
        }
        else
        {
            assert(value <= 0.0001);
        }
    }

    std::size_t n = std::size(sol);

    // This constraint cannot be violated if zero or one item is selected
    if (n < 2)
    {
        return;
    }

    // Build a subgraph from selected items
    DirectedGraph subgraph(n);

    for (std::size_t u = 0; u < n; ++u)
    {
        const auto i = sol[u];

        for (std::size_t v = 0; v < n; ++v)
        {
            if (u == v)
            {
                continue;
            }

            const auto j = sol[v];

            if (_adjacency_matrix[i][j])
            {
                subgraph.add_edge(u, v);
            }
        }
    }

    // Count the number of cuts for this round
    std::size_t n_cuts = 0;

    // Exception thrown when the maximum number of cuts is reached
    struct MaxCuts
    {
    };

    // This callback is executed for every cycle detected
    auto create_cut = [&](auto first, auto last) {
        auto cycle_length = std::distance(first, last);

#ifndef NDEBUG
        // Check whether it is a valid cycle
        assert(cycle_length >= 2);
        for (auto i = first; i != last; ++i)
        {
            assert((i + 1 == last) || _adjacency_matrix[sol[*i]][sol[*(i + 1)]]);
            assert((i + 1 != last) || _adjacency_matrix[sol[*i]][sol[*first]]);
        }
#endif

        if constexpr (PrintCuts)
        {
            std::cout << "cut:";
        }

        IloRange cons(env, -IloInfinity, IloNum(cycle_length - 1));
        for (auto i = first; i != last; ++i)
        {
            if constexpr (PrintCuts)
            {
                std::cout << ' ' << sol[*i] + 1;
            }

            cons.setLinearCoef(_x[sol[*i]], 1.);
        }
        add(cons).end();

        if constexpr (PrintCuts)
        {
            std::cout << std::endl;
        }

        if (_max_cuts > 0 && ++n_cuts >= _max_cuts)
        {
            throw MaxCuts();
        }
    };

    // Find cycles
    try
    {
        CycleFinder cycle_finder(create_cut);
        subgraph.depth_first_visit(cycle_finder);
    }
    catch (const MaxCuts &)
    {
        // Maximum number of cuts reached
    }
}

IloCplex::CallbackI *ConsCallbackCycle::duplicateCallback() const
{
    return (new (getEnv()) ConsCallbackCycle(*this));
}

IloCplex::Callback ConsCallbackCycle::create(IloEnv env,
                                             const IloNumVarArray &x,
                                             const IloArray<IloBoolArray> &adjacency_matrix,
                                             std::size_t max_cuts)
{
    return IloCplex::Callback(
        new (env) ConsCallbackCycle(env,
                                    x,
                                    adjacency_matrix,
                                    max_cuts));
}