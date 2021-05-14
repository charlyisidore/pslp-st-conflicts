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

#include "mip_cplex_flow.hpp"
#include "logging.hpp"
#include <cassert>
#include <sstream>

constexpr double Epsilon = 1e-6;

using namespace pslp::mip;

void MipCplexFlow::create(const Problem &prob)
{
    const auto n_initial_items = std::size(prob.initial_positions);

    _x = IloArray<IloNumVarArray>(_env, prob.n_items);
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        _x[i] = IloNumVarArray(_env, prob.n_items, 0., 1., IloNumVar::Bool);
        for (std::size_t j = 0; j < prob.n_items; ++j)
        {
            if (i == j)
            {
                continue;
            }

            std::ostringstream oss;
            oss << "x{" << i + 1 << ',' << j + 1 << '}';
            _x[i][j].setName(oss.str().c_str());
            _model.add(_x[i][j]);
        }
    }

    _s = IloNumVarArray(_env, prob.n_items, 0., 1., IloNumVar::Bool);
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        std::ostringstream oss;
        oss << "s{" << i + 1 << '}';
        _s[i].setName(oss.str().c_str());
    }
    _model.add(_s);

    _t = IloNumVarArray(_env, prob.n_items, 0., 1., IloNumVar::Bool);
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        std::ostringstream oss;
        oss << "t{" << i + 1 << '}';
        _t[i].setName(oss.str().c_str());
    }
    _model.add(_t);

    _l = IloNumVarArray(_env,
                        prob.n_items,
                        0.,
                        IloNum(prob.capacity) - 1.,
                        IloNumVar::Float);
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        std::ostringstream oss;
        oss << "l{" << i + 1 << '}';
        _l[i].setName(oss.str().c_str());
    }
    _model.add(_l);

    // Minimize the number of used stacks
    {
        _objective = IloMinimize(_env, IloSum(_s), "St");
        _model.add(_objective);
    }

    // Fixed items
    for (std::size_t i = 0; i < n_initial_items; ++i)
    {
        assert(prob.initial_layout.contains(i));

        const auto k_i = prob.initial_layout.stack_of(i);
        const auto l_i = prob.initial_layout.level_of(i);

        for (std::size_t j = 0; j < n_initial_items; ++j)
        {
            if (i == j)
            {
                continue;
            }

            IloNum val = ((k_i == prob.initial_layout.stack_of(j) &&
                           l_i == prob.initial_layout.level_of(j) + 1)
                              ? 1.
                              : 0.);
            _x[i][j].setBounds(val, val);
        }
    }

    // In flow
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        std::ostringstream oss;
        oss << "in{" << i + 1 << '}';

        IloRange cons(_env, 1., 1., oss.str().c_str());
        for (std::size_t j = 0; j < prob.n_items; ++j)
        {
            if (i == j || prob.conflict_matrix(j, i))
            {
                continue;
            }
            cons.setLinearCoef(_x[j][i], 1.);
        }
        cons.setLinearCoef(_s[i], 1.);
        _model.add(cons);
    }

    // Out flow
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        std::ostringstream oss;
        oss << "out{" << i + 1 << '}';

        IloRange cons(_env, 1., 1., oss.str().c_str());
        for (std::size_t j = 0; j < prob.n_items; ++j)
        {
            if (i == j || prob.conflict_matrix(i, j))
            {
                continue;
            }
            cons.setLinearCoef(_x[i][j], 1.);
        }
        cons.setLinearCoef(_t[i], 1.);
        _model.add(cons);
    }

    // MTZ constraints
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        for (std::size_t j = 0; j < prob.n_items; ++j)
        {
            if (i == j || prob.conflict_matrix(i, j))
            {
                continue;
            }

            std::ostringstream oss;
            oss << "mtz{" << i + 1 << ',' << j + 1 << '}';

            IloRange cons(_env,
                          -IloInfinity,
                          IloNum(prob.capacity) - 1.,
                          oss.str().c_str());
            cons.setLinearCoef(_l[i], 1.);
            cons.setLinearCoef(_l[j], -1.);
            cons.setLinearCoef(_x[i][j], IloNum(prob.capacity));
            _model.add(cons);
        }
    }
}

void MipCplexFlow::initialize(const Problem &prob)
{
    const auto n_stacks = std::size(_heuristic_sol);
    IloNumVarArray vars(_env);
    IloNumArray vals(_env);

    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        const auto h = std::size(_heuristic_sol[k]);

        for (std::size_t l = 0; l < h; ++l)
        {
            const auto i = _heuristic_sol[k][l];
            if (l > 0)
            {
                const auto i2 = _heuristic_sol[k][l - 1];
                for (std::size_t j = 0; j < prob.n_items; ++j)
                {
                    if (i == j)
                    {
                        continue;
                    }
                    vars.add(_x[i][j]);
                    vals.add((j == i2) ? 1. : 0.);
                }
            }
            if (l == 0)
            {
                vars.add(_t[i]);
                vals.add(1.);
            }
            else if (l == h - 1)
            {
                vars.add(_s[i]);
                vals.add(1.);
            }
            vars.add(_l[i]);
            vals.add(IloNum(h - l - 1));
        }
    }

    _cplex.addMIPStart(vars, vals);

    vars.end();
    vals.end();
}

MipCplexFlow::Solution MipCplexFlow::solution(const Problem &prob)
{
    double obj = _cplex.getValue(_objective);
    const auto n_initial_items = std::size(prob.initial_positions);
    std::vector<std::size_t> stack_of(prob.n_items, prob.n_items);
    std::vector<std::size_t> level_of(prob.n_items, prob.capacity);
    std::vector<bool> used_stacks(prob.n_items, false);

    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        if (_cplex.getValue(_t[i]) < 0.5)
        {
            continue;
        }

        if (stack_of[i] != prob.n_items ||
            level_of[i] != prob.capacity)
        {
            LOG_ERROR("item %zu already assigned", i + 1);
        }

        if (i < n_initial_items)
        {
            const auto k = prob.initial_positions[i];
            if (used_stacks[k])
            {
                LOG_ERROR("stack %zu already used", k + 1);
                return {};
            }
            stack_of[i] = k;
            used_stacks[k] = true;
        }
        else
        {
            std::size_t k = 0;
            while (used_stacks[k])
            {
                ++k;
                assert(k < std::size(used_stacks));
            }
            stack_of[i] = k;
            used_stacks[k] = true;
        }
        level_of[i] = 0;

        std::size_t j = i;

        while (_cplex.getValue(_s[j]) < 0.5)
        {
            bool found = false;

            for (std::size_t k = 0; k < prob.n_items && !found; ++k)
            {
                if (j == k || _cplex.getValue(_x[k][j]) < 0.5)
                {
                    continue;
                }

                if (stack_of[k] != prob.n_items ||
                    level_of[k] != prob.capacity)
                {
                    LOG_ERROR("item %zu already assigned", k + 1);
                    return {};
                }

                assert(_cplex.getValue(_l[j]) > _cplex.getValue(_l[k]));

                stack_of[k] = stack_of[j];
                level_of[k] = level_of[j] + 1;
                j = k;
                found = true;
            }

            if (!found)
            {
                LOG_ERROR("cannot build solution");
                return {};
            }
        }
    }

    LOG_DEBUG("obj: %g", obj);
    LOG_DEBUG("used stacks: %zu",
              std::count(std::begin(used_stacks), std::end(used_stacks), true));

    return {obj, stack_of, level_of};
}

void MipCplexFlow::write(const Problem &prob,
                         const Solution &sol,
                         nlohmann::json &json)
{
    MipCplex::write(prob, sol, json);

    if (!has_solution())
    {
        return;
    }

    json["vars"] = nlohmann::json::object();
    json["x"] = nlohmann::json::array();
    json["s"] = nlohmann::json::array();
    json["t"] = nlohmann::json::array();
    json["l"] = nlohmann::json::array();

    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        for (std::size_t j = 0; j < prob.n_items; ++j)
        {
            if (i == j)
            {
                continue;
            }

            auto value = _cplex.getValue(_x[i][j]);
            if (value > Epsilon)
            {
                const char *name = _x[i][j].getName();
                json["vars"][name] = value;
                json["x"].push_back({i + 1, j + 1});
            }
        }
    }

    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        auto value = _cplex.getValue(_s[i]);
        if (value > Epsilon)
        {
            const char *name = _s[i].getName();
            json["vars"][name] = value;
            json["s"].push_back(i + 1);
        }
    }

    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        auto value = _cplex.getValue(_t[i]);
        if (value > Epsilon)
        {
            const char *name = _t[i].getName();
            json["vars"][name] = value;
            json["t"].push_back(i + 1);
        }
    }

    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        auto value = _cplex.getValue(_l[i]);
        if (value > Epsilon)
        {
            const char *name = _l[i].getName();
            json["vars"][name] = value;
        }
        json["l"].push_back(value);
    }
}