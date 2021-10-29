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

#include "mip_cplex_bp.hpp"
#include "graph_algorithm.hpp"
#include "logging.hpp"
#include <sstream>

constexpr double Epsilon = 1e-6;

using namespace pslp::mip;

void MipCplexBp::create(const Problem &prob)
{
    const auto n_stacks = std::size(_heuristic_sol);
    const auto n_initial_items = std::size(prob.initial_positions);

    _x = IloArray<IloNumVarArray>(_env, prob.n_items);
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        _x[i] = IloNumVarArray(_env, n_stacks, 0., 1., IloNumVar::Bool);
        for (std::size_t k = 0; k < n_stacks; ++k)
        {
            std::ostringstream oss;
            oss << "x{" << i + 1 << ',' << k + 1 << '}';
            _x[i][k].setName(oss.str().c_str());
        }
        _model.add(_x[i]);
    }

    _y = IloNumVarArray(_env, n_stacks, 0., 1., IloNumVar::Bool);
    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        std::ostringstream oss;
        oss << "y{" << k + 1 << '}';
        _y[k].setName(oss.str().c_str());
    }
    _model.add(_y);

    // Minimize the number of used stacks
    {
        _objective = IloMinimize(_env, IloSum(_y), "St");
        _model.add(_objective);
    }

    // Fixed items
    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        bool used = false;

        for (std::size_t i = 0; i < n_initial_items; ++i)
        {
            IloNum val = 0.;
            if (prob.initial_positions[i] == k)
            {
                val = 1.;
                used = true;
            }
            _x[i][k].setBounds(val, val);
        }

        if (used)
        {
            _y[k].setBounds(1., 1.);
        }
    }

    // Capacity constraints
    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        std::ostringstream oss;
        oss << "capacity{" << k + 1 << '}';

        IloRange cons(_env, -IloInfinity, 0., oss.str().c_str());
        for (std::size_t i = 0; i < prob.n_items; ++i)
        {
            cons.setLinearCoef(_x[i][k], 1.);
        }
        cons.setLinearCoef(_y[k], -IloNum(prob.capacity));
        _model.add(cons);
    }

    // Putting an item in stack k means using k
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        for (std::size_t k = 0; k < n_stacks; ++k)
        {
            std::ostringstream oss;
            oss << "use{" << i + 1 << ',' << k + 1 << '}';

            IloRange cons(_env, -IloInfinity, 0., oss.str().c_str());
            cons.setLinearCoef(_x[i][k], 1.);
            cons.setLinearCoef(_y[k], -1.);
            _model.add(cons);
        }
    }

    // Assignment constraints
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        std::ostringstream oss;
        oss << "assign{" << i + 1 << '}';

        IloRange cons(_env, 1., 1., oss.str().c_str());
        for (std::size_t k = 0; k < n_stacks; ++k)
        {
            cons.setLinearCoef(_x[i][k], 1.);
        }
        _model.add(cons);
    }

    // Conflicts in empty stacks
    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        if (prob.initial_layout.height(k) > 0)
        {
            continue;
        }

        for (std::size_t i = n_initial_items; i < prob.n_items; ++i)
        {
            std::ostringstream oss;
            oss << "conflict{" << i + 1 << ',' << k + 1 << '}';

            IloRange cons(_env,
                          -IloInfinity,
                          IloNum(prob.capacity),
                          oss.str().c_str());
            for (std::size_t j = n_initial_items; j < prob.n_items; ++j)
            {
                if (i != j &&
                    prob.conflict_matrix(i, j) &&
                    prob.conflict_matrix(j, i))
                {
                    cons.setLinearCoef(_x[j][k], 1.);
                }
            }
            cons.setLinearCoef(_x[i][k], IloNum(prob.capacity));
            _model.add(cons);
        }
    }

    // Conflicts in stacks used by initial items
    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        const auto h = prob.initial_layout.height(k);
        if (h == 0)
        {
            continue;
        }

        for (std::size_t i = 0; i < prob.n_items; ++i)
        {
            // If item i is initial, we require it is topmost
            if (i < n_initial_items &&
                prob.initial_layout.level_of(i) < h - 1)
            {
                continue;
            }

            std::ostringstream oss;
            oss << "conflict{" << i + 1 << ',' << k + 1 << '}';

            IloRange cons(_env,
                          -IloInfinity,
                          IloNum(prob.capacity - h),
                          oss.str().c_str());
            for (std::size_t j = n_initial_items; j < prob.n_items; ++j)
            {
                if (i != j &&
                    prob.conflict_matrix(i, j) &&
                    prob.conflict_matrix(j, i))
                {
                    cons.setLinearCoef(_x[j][k], 1.);
                }
            }
            cons.setLinearCoef(_x[i][k], IloNum(prob.capacity - h));
            _model.add(cons);
        }
    }
}

void MipCplexBp::initialize(const Problem &prob)
{
    const auto n_stacks = std::size(_heuristic_sol);
    IloNumVarArray vars(_env);
    IloNumArray vals(_env);

    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        for (std::size_t i = 0; i < prob.n_items; ++i)
        {
            bool found = std::find(std::begin(_heuristic_sol[k]),
                                   std::end(_heuristic_sol[k]),
                                   i) != std::end(_heuristic_sol[k]);
            vars.add(_x[i][k]);
            vals.add(found ? 1. : 0.);
        }

        vars.add(_y[k]);
        vals.add(std::empty(_heuristic_sol[k]) ? 0. : 1.);
    }

    _cplex.addMIPStart(vars, vals);

    vars.end();
    vals.end();
}

MipCplexBp::Solution MipCplexBp::solution(const Problem &prob)
{
    const auto n_stacks = std::size(_heuristic_sol);
    double obj = _cplex.getValue(_objective);
    std::vector<std::size_t> stack_of(prob.n_items, prob.n_items);
    std::vector<std::size_t> level_of(prob.n_items, prob.capacity);
    std::vector<std::vector<std::size_t>> layout(n_stacks);

    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        for (std::size_t i = 0; i < prob.n_items; ++i)
        {
            if (_cplex.getValue(_x[i][k]) < 0.5)
            {
                continue;
            }

            if (stack_of[i] != prob.n_items)
            {
                LOG_ERROR("item %zu already assigned", i + 1);
            }

            stack_of[i] = k;
            layout[k].push_back(i);
        }
    }

    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        const auto s = layout[k];
        const auto h = std::size(s);
        DirectedGraph subgraph(h);

        for (std::size_t u = 0; u < h; ++u)
        {
            const auto i = s[u];
            for (std::size_t v = 0; v < h; ++v)
            {
                if (u == v)
                {
                    continue;
                }
                const auto j = s[v];
                if (prob.conflict_matrix(i, j))
                {
                    subgraph.add_edge(u, v);
                }
            }
        }

        // Find a topological sorting
        try
        {
            std::size_t pos = std::size(s);
            TopologicalSorter topo_sorter([&](auto v) {
                assert(pos > 0);
                layout[k][--pos] = s[v];
            });
            depth_first_visit(subgraph, topo_sorter);
        }
        catch (const NotADag &)
        {
            LOG_ERROR("could not find a feasible order");
            layout[k] = s;
        }

        for (std::size_t l = 0; l < std::size(layout[k]); ++l)
        {
            level_of[layout[k][l]] = l;
        }
    }

    return {obj, stack_of, level_of};
}

void MipCplexBp::write(const Problem &prob,
                       const Solution &sol,
                       nlohmann::json &json)
{
    const auto n_stacks = std::size(_heuristic_sol);

    MipCplex::write(prob, sol, json);

    if (!has_solution())
    {
        return;
    }

    json["vars"] = nlohmann::json::object();
    json["x"] = nlohmann::json::array();
    json["y"] = nlohmann::json::array();

    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        for (std::size_t k = 0; k < n_stacks; ++k)
        {
            auto value = _cplex.getValue(_x[i][k]);
            if (value > Epsilon)
            {
                const char *name = _x[i][k].getName();
                json["vars"][name] = value;
                json["x"].push_back({i + 1, k + 1});
            }
        }
    }

    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        auto value = _cplex.getValue(_y[k]);
        if (value > Epsilon)
        {
            const char *name = _y[k].getName();
            json["vars"][name] = value;
            json["y"].push_back(k + 1);
        }
    }
}