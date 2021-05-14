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

#include "cycle_constraint_callback.hpp"
#include "graph.hpp"
#include "graph_algorithm.hpp"
#include <cassert>
#include <iostream>

// Print cuts
constexpr bool PrintCuts = false;

using namespace pslp::bp;

CycleConstraintCallback::CycleConstraintCallback(IloEnv env,
                                                 const IloNumVarArray &x,
                                                 const IloArray<IloBoolArray> &adjacency_matrix,
                                                 IloInt max_cuts)
    : IloCplex::LazyConstraintCallbackI(env),
      _x(x),
      _adjacency_matrix(adjacency_matrix),
      _max_cuts(max_cuts)
{
}

void CycleConstraintCallback::main()
{
    IloEnv env = getEnv();

    // Get the current solution
    IloIntArray s(env);

    for (IloInt i = 0; i < _x.getSize(); ++i)
    {
        const auto value = getValue(_x[i]);
        if (value > 0.5)
        {
            assert(value >= 0.9999);
            s.add(i);
        }
        else
        {
            assert(value <= 0.0001);
        }
    }

    IloInt n = s.getSize();

    // This constraint cannot be violated if zero or one item is selected
    if (n < 2)
    {
        return;
    }

    // Build a subgraph from selected items
    DirectedGraph subgraph(n);

    for (IloInt u = 0; u < n; ++u)
    {
        IloInt i = s[u];

        for (IloInt v = 0; v < n; ++v)
        {
            if (u == v)
            {
                continue;
            }

            IloInt j = s[v];

            if (_adjacency_matrix[i][j])
            {
                subgraph.add_edge(u, v);
            }
        }
    }

    // Count the number of cuts for this round
    IloInt n_cuts = 0;

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
            assert((i + 1 == last) || _adjacency_matrix[s[*i]][s[*(i + 1)]]);
            assert((i + 1 != last) || _adjacency_matrix[s[*i]][s[*first]]);
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
                std::cout << ' ' << s[*i] + 1;
            }

            cons.setLinearCoef(_x[s[*i]], 1.);
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

    s.end();
}

IloCplex::CallbackI *CycleConstraintCallback::duplicateCallback() const
{
    return (new (getEnv()) CycleConstraintCallback(*this));
}

IloCplex::Callback CycleConstraintCallback::create(IloEnv env,
                                                   const IloNumVarArray &x,
                                                   const IloArray<IloBoolArray> &adjacency_matrix,
                                                   IloInt max_cuts)
{
    return IloCplex::Callback(
        new (env) CycleConstraintCallback(env,
                                          x,
                                          adjacency_matrix,
                                          max_cuts));
}