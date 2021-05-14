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

#include "mip_cplex_3ind.hpp"
#include "logging.hpp"
#include <cassert>
#include <sstream>

constexpr double Epsilon = 1e-6;

using namespace pslp::mip;

void MipCplex3ind::create(const Problem &prob)
{
    const auto n_stacks = std::size(_heuristic_sol);

    _x = IloArray<IloArray<IloNumVarArray>>(_env, prob.n_items);
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        _x[i] = IloArray<IloNumVarArray>(_env, n_stacks);
        for (std::size_t k = 0; k < n_stacks; ++k)
        {
            _x[i][k] = IloNumVarArray(_env,
                                      prob.capacity,
                                      0.,
                                      1.,
                                      IloNumVar::Bool);

            for (std::size_t l = 0; l < prob.capacity; ++l)
            {
                std::ostringstream oss;
                oss << "x{" << i + 1 << ',' << k + 1 << ',' << l + 1 << '}';
                _x[i][k][l].setName(oss.str().c_str());
            }

            _model.add(_x[i][k]);
        }
    }

    // Minimize the number of used stacks
    {
        _objective = IloMinimize(_env, 0., "St");
        for (std::size_t i = 0; i < prob.n_items; ++i)
        {
            for (std::size_t k = 0; k < n_stacks; ++k)
            {
                _objective.setLinearCoef(_x[i][k][0], 1.);
            }
        }
        _model.add(_objective);
    }

    // Fixed items
    for (std::size_t i = 0; i < std::size(prob.initial_positions); ++i)
    {
        assert(prob.initial_layout.contains(i));

        for (std::size_t k = 0; k < n_stacks; ++k)
        {
            for (std::size_t l = 0; l < prob.capacity; ++l)
            {
                IloNum val = ((prob.initial_layout.stack_of(i) == k &&
                               prob.initial_layout.level_of(i) == l)
                                  ? 1.
                                  : 0.);
                _x[i][k][l].setBounds(val, val);
            }
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
            for (std::size_t l = 0; l < prob.capacity; ++l)
            {
                cons.setLinearCoef(_x[i][k][l], 1.);
            }
        }
        _model.add(cons);
    }

    // At most one item per slot
    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        for (std::size_t l = 0; l < prob.capacity; ++l)
        {
            std::ostringstream oss;
            oss << "slot{" << k + 1 << ',' << l + 1 << '}';

            IloRange cons(_env, -IloInfinity, 1., oss.str().c_str());
            for (std::size_t i = 0; i < prob.n_items; ++i)
            {
                cons.setLinearCoef(_x[i][k][l], 1.);
            }
            _model.add(cons);
        }
    }

    // Conflict constraints
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        for (std::size_t k = 0; k < n_stacks; ++k)
        {
            for (std::size_t l = 1; l < prob.capacity; ++l)
            {
                std::ostringstream oss;
                oss << "conflict{" << i + 1 << ',' << k + 1 << ',' << l + 1 << '}';

                IloRange cons(_env, 0., IloInfinity, oss.str().c_str());
                for (std::size_t j = 0; j < prob.n_items; ++j)
                {
                    if (i == j || prob.conflict_matrix(i, j))
                    {
                        continue;
                    }
                    cons.setLinearCoef(_x[j][k][l - 1], 1.);
                }
                cons.setLinearCoef(_x[i][k][l], -1.);
                _model.add(cons);
            }
        }
    }
}

void MipCplex3ind::initialize(const Problem &prob)
{
    const auto n_stacks = std::size(_heuristic_sol);
    IloNumVarArray vars(_env);
    IloNumArray vals(_env);

    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        const auto h = std::size(_heuristic_sol[k]);

        for (std::size_t l = 0; l < h; ++l)
        {
            const auto j = _heuristic_sol[k][l];
            for (std::size_t i = 0; i < prob.n_items; ++i)
            {
                vars.add(_x[i][k][l]);
                vals.add((i == j) ? 1. : 0.);
            }
        }

        for (std::size_t l = h; l < prob.capacity; ++l)
        {
            for (std::size_t i = 0; i < prob.n_items; ++i)
            {
                vars.add(_x[i][k][l]);
                vals.add(0.);
            }
        }
    }

    _cplex.addMIPStart(vars, vals);

    vars.end();
    vals.end();
}

MipCplex3ind::Solution MipCplex3ind::solution(const Problem &prob)
{
    const auto n_stacks = std::size(_heuristic_sol);
    double obj = _cplex.getValue(_objective);
    std::vector<std::size_t> stack_of(prob.n_items, prob.n_items);
    std::vector<std::size_t> level_of(prob.n_items, prob.capacity);

    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        for (std::size_t k = 0; k < n_stacks; ++k)
        {
            for (std::size_t l = 0; l < prob.capacity; ++l)
            {
                if (_cplex.getValue(_x[i][k][l]) < 0.5)
                {
                    continue;
                }

                if (stack_of[i] != prob.n_items ||
                    level_of[i] != prob.capacity)
                {
                    LOG_ERROR("item %zu already assigned", i + 1);
                }

                stack_of[i] = k;
                level_of[i] = l;
            }
        }
    }

    return {obj, stack_of, level_of};
}

void MipCplex3ind::write(const Problem &prob,
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

    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        for (std::size_t k = 0; k < n_stacks; ++k)
        {
            for (std::size_t l = 0; l < prob.capacity; ++l)
            {
                auto value = _cplex.getValue(_x[i][k][l]);
                if (value > Epsilon)
                {
                    const char *name = _x[i][k][l].getName();
                    json["vars"][name] = value;
                    json["x"].push_back({i + 1, k + 1, l + 1});
                }
            }
        }
    }
}