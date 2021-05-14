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

#include "pricer_cplex_pa.hpp"
#include "logging.hpp"
#include "pprint.hpp"
#include <cassert>
#include <sstream>

// Print the solution
constexpr bool PrintSolution = false;

// Print CPLEX variables
constexpr bool PrintCplexVariables = false;

using namespace pslp::bp;

const char *PricerCplexPa::Properties::Name = "cplex_pa";
const char *PricerCplexPa::Properties::Desc = "CPLEX variable pricer with position assignment formulation";
const int PricerCplexPa::Properties::Priority = 0;
const SCIP_Bool PricerCplexPa::Properties::Delay = true;

PricerCplexPa::PricerCplexPa(SCIP *scip)
    : PricerCplex(scip,
                  Properties::Name,
                  Properties::Desc,
                  Properties::Priority,
                  Properties::Delay)
{
}

SCIP_RETCODE PricerCplexPa::init(SCIP *scip, const Problem &prob)
{
    SCIP_CALL(PricerCplex::init(scip, prob));

    const auto n = prob.conflict_graph.n_vertices();
    const auto b = prob.capacity;

    _y = IloArray<IloNumVarArray>(_env, n);
    for (std::size_t i = 0; i < n; ++i)
    {
        _y[i] = IloNumVarArray(_env, b, 0, 1, IloNumVar::Bool);
        for (std::size_t l = 0; l < b; ++l)
        {
            std::ostringstream oss;
            oss << "y{" << i + 1 << ',' << l + 1 << '}';
            _y[i][l].setName(oss.str().c_str());
        }
        _model.add(_y[i]);
    }

    // Selection constraints: a_i = sum_l y_il
    for (std::size_t i = 0; i < n; ++i)
    {
        std::ostringstream oss;
        oss << "select{" << i + 1 << '}';

        IloRange cons(_env, 0., 0., oss.str().c_str());
        cons.setLinearCoef(_x[i], 1.);
        for (std::size_t l = 0; l < b; ++l)
        {
            cons.setLinearCoef(_y[i][l], -1.);
        }
        _model.add(cons);
    }

    // Assignment constraints: sum_i y_il <= 1
    for (std::size_t l = 0; l < b; ++l)
    {
        std::ostringstream oss;
        oss << "assign{" << l + 1 << '}';

        IloRange cons(_env, -IloInfinity, 1., oss.str().c_str());
        for (std::size_t i = 0; i < n; ++i)
        {
            cons.setLinearCoef(_y[i][l], 1.);
        }
        _model.add(cons);
    }

    // Gravity constraints: sum_i y_il >= sum_i y_il+1
    for (std::size_t l = 0; l < b - 1; ++l)
    {
        std::ostringstream oss;
        oss << "gravity{" << l + 1 << '}';

        IloRange cons(_env, 0., IloInfinity, oss.str().c_str());
        for (std::size_t i = 0; i < n; ++i)
        {
            cons.setLinearCoef(_y[i][l], 1.);
            cons.setLinearCoef(_y[i][l + 1], -1.);
        }
        _model.add(cons);
    }

    if (prob.has_adjacent_conflicts)
    {
        // Adjacent conflicts: y_il + y_jl' <= 1
        for (std::size_t i = 0; i < n; ++i)
        {
            for (auto j : prob.conflict_graph.adjacent_vertices(i))
            {
                assert(i != j);

                for (std::size_t l = 0; l < b - 1; ++l)
                {
                    std::ostringstream oss;
                    oss << "adjc{"
                        << i + 1 << ','
                        << j + 1 << ','
                        << l + 1 << '}';

                    IloRange cons(_env, -IloInfinity, 1., oss.str().c_str());
                    cons.setLinearCoef(_y[i][l + 1], 1.);
                    cons.setLinearCoef(_y[j][l], 1.);
                    _model.add(cons);
                }
            }
        }
    }
    else
    {
        // Distant conflicts: y_il + y_jl' <= 1
        for (std::size_t i = 0; i < n; ++i)
        {
            for (auto j : prob.conflict_graph.adjacent_vertices(i))
            {
                assert(i != j);

                for (std::size_t l_i = 0; l_i < b; ++l_i)
                {
                    for (std::size_t l_j = 0; l_j < l_i; ++l_j)
                    {
                        std::ostringstream oss;
                        oss << "distc{"
                            << i + 1 << ','
                            << j + 1 << ','
                            << l_i + 1 << ','
                            << l_j + 1 << '}';

                        IloRange cons(_env, -IloInfinity, 1., oss.str().c_str());
                        cons.setLinearCoef(_y[i][l_i], 1.);
                        cons.setLinearCoef(_y[j][l_j], 1.);
                        _model.add(cons);
                    }
                }
            }
        }
    }

    return SCIP_OKAY;
}

SCIP_RETCODE PricerCplexPa::solution(SCIP *scip,
                                     const Problem &prob,
                                     std::size_t k,
                                     std::vector<Solution> &sols)
{
    const auto n = prob.conflict_graph.n_vertices();
    const auto b = prob.capacity;

    SCIP_Real obj = _cplex.getObjValue(k);

    // Can this solution improve the master problem's objective?
    if (!SCIPisPositive(scip, obj))
    {
        return SCIP_OKAY;
    }

    LOG_DEBUG("found CPLEX PA pricing sol (obj: %g)", obj);

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
        for (IloNum i = 0; i < _y.getSize(); ++i)
        {
            for (IloNum l = 0; l < _y[i].getSize(); ++l)
            {
                double value = _cplex.getValue(_y[i][l], k);
                if (!SCIPisEQ(scip, value, 0.))
                {
                    LOG_DEBUG("%s  %g", _y[i][l].getName(), value);
                }
            }
        }
    }

    Solution sol;

    // Get all selected items in the right order, level by level
    for (std::size_t l = 0; l < b; ++l)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            if (_cplex.getValue(_y[i][l], k) > 0.5)
            {
#ifndef NDEBUG
                assert(_cplex.getValue(_x[i], k) > 0.5);
                for (std::size_t j = i + 1; j < n; ++j)
                {
                    assert(_cplex.getValue(_y[j][l], k) < 0.5);
                }
#endif

                sol.push_back(i);
                break;
            }
        }

        // Because of gravity constraints, we should not find any item above
        if (std::size(sol) != l + 1)
        {
#ifndef NDEBUG
            assert(std::size(sol) == l);

            // Ensure that no item is put above
            for (std::size_t l_i = l + 1; l_i < b; ++l_i)
            {
                for (std::size_t i = 0; i < n; ++i)
                {
                    assert(_cplex.getValue(_y[i][l_i], k) < 0.5);
                }
            }
#endif

            break;
        }
    }

    if constexpr (PrintSolution)
    {
        LOG_TRACE("-> sol: %s",
                  pprint_n(std::size(sol),
                           [&sol](auto i) { return sol[i] + 1; })
                      .c_str());
    }

    sols.push_back(sol);

    return SCIP_OKAY;
}