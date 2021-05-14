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

#include "pricer_cplex_lo.hpp"
#include "logging.hpp"
#include "pprint.hpp"
#include <cassert>
#include <sstream>

// Print the solution
constexpr bool PrintSolution = false;

// Print CPLEX variables
constexpr bool PrintCplexVariables = false;

using namespace pslp::bp;

const char *PricerCplexLo::Properties::Name = "cplex_lo";
const char *PricerCplexLo::Properties::Desc = "CPLEX variable pricer with linear ordering formulation";
const int PricerCplexLo::Properties::Priority = 0;
const SCIP_Bool PricerCplexLo::Properties::Delay = true;

PricerCplexLo::PricerCplexLo(SCIP *scip)
    : PricerCplex(scip,
                  Properties::Name,
                  Properties::Desc,
                  Properties::Priority,
                  Properties::Delay)
{
}

SCIP_RETCODE PricerCplexLo::init(SCIP *scip, const Problem &prob)
{
    if (prob.has_adjacent_conflicts)
    {
        LOG_ERROR("adjacent conflicts not supported");
        return SCIP_ERROR;
    }

    SCIP_CALL(PricerCplex::init(scip, prob));

    const auto n = prob.conflict_graph.n_vertices();

    _y = IloArray<IloNumVarArray>(_env, n);
    for (std::size_t i = 0; i < n; ++i)
    {
        _y[i] = IloNumVarArray(_env, n, 0, 1, IloNumVar::Bool);
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j || prob.conflict_matrix(i, j))
            {
                continue;
            }

            std::ostringstream oss;
            oss << "y{" << i + 1 << ',' << j + 1 << '}';
            _y[i][j].setName(oss.str().c_str());
            _model.add(_y[i][j]);
        }
    }

    // Either above or below: y_ij + y_ji <= 1
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            // Skip if y_ij = 0 or y_ji = 0
            if (prob.conflict_matrix(i, j) || prob.conflict_matrix(j, i))
            {
                continue;
            }

            std::ostringstream oss;
            oss << "order{" << i + 1 << ',' << j + 1 << '}';

            IloRange cons(_env, -IloInfinity, 1., oss.str().c_str());
            cons.setLinearCoef(_y[i][j], 1.);
            cons.setLinearCoef(_y[j][i], 1.);
            _model.add(cons);
        }
    }

    // Transitivity: y_ij + y_jk + y_ki <= 2
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            // If y_ij = 0, the constraint is satisfied
            if (i == j || prob.conflict_matrix(i, j))
            {
                continue;
            }

            for (std::size_t k = 0; k < n; ++k)
            {
                // If y_jk = 0 or y_ki = 0, the constraint is satisfied
                if (i == k || j == k ||
                    prob.conflict_matrix(j, k) ||
                    prob.conflict_matrix(k, i))
                {
                    continue;
                }

                std::ostringstream oss;
                oss << "trans{" << i + 1 << ',' << j + 1 << ',' << k + 1 << '}';

                IloRange cons(_env, -IloInfinity, 2., oss.str().c_str());
                cons.setLinearCoef(_y[i][j], 1.);
                cons.setLinearCoef(_y[j][k], 1.);
                cons.setLinearCoef(_y[k][i], 1.);
                _model.add(cons);
            }
        }
    }

    // Selection: x_i + x_j - y_ij - y_ji <= 1
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            std::ostringstream oss;
            oss << "select{" << i + 1 << ',' << j + 1 << '}';

            IloRange cons(_env, -IloInfinity, 1., oss.str().c_str());
            cons.setLinearCoef(_x[i], 1.);
            cons.setLinearCoef(_x[j], 1.);
            if (!prob.conflict_matrix(i, j))
            {
                cons.setLinearCoef(_y[i][j], -1.);
            }
            if (!prob.conflict_matrix(j, i))
            {
                cons.setLinearCoef(_y[j][i], -1.);
            }
            _model.add(cons);
        }
    }

    return SCIP_OKAY;
}

SCIP_RETCODE PricerCplexLo::solution(SCIP *scip,
                                     const Problem &prob,
                                     std::size_t k,
                                     std::vector<Solution> &sols)
{
    const auto n = prob.conflict_graph.n_vertices();

    SCIP_Real obj = _cplex.getObjValue(k);

    // Can this solution improve the master problem's objective?
    if (!SCIPisPositive(scip, obj))
    {
        return SCIP_OKAY;
    }

    LOG_DEBUG("found CPLEX LO pricing sol (obj: %g)", obj);

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
            for (IloNum j = 0; j < _y[i].getSize(); ++j)
            {
                if (i == j || prob.conflict_matrix(i, j))
                {
                    continue;
                }
                double value = _cplex.getValue(_y[i][j], k);
                if (!SCIPisEQ(scip, value, 0.))
                {
                    LOG_DEBUG("%s  %g", _y[i][j].getName(), value);
                }
            }
        }
    }

    Solution sol;

    // Get all selected items, sort them later
    for (std::size_t i = 0; i < n; ++i)
    {
        if (_cplex.getValue(_x[i], k) > 0.5)
        {
            sol.push_back(i);
        }
    }

#ifndef NDEBUG
    // Check that some constraints were satisfied
    for (std::size_t i = 0; i < n; ++i)
    {
        bool x_i = _cplex.getValue(_x[i], k) > 0.5;

        if (!x_i)
        {
            continue;
        }

        for (std::size_t j = 0; j < n; ++j)
        {
            bool x_j = _cplex.getValue(_x[j], k) > 0.5;

            if (i == j || !x_j)
            {
                continue;
            }

            bool y_ij = !prob.conflict_matrix(i, j) &&
                        (_cplex.getValue(_y[i][j], k) > 0.5);
            bool y_ji = !prob.conflict_matrix(j, i) &&
                        (_cplex.getValue(_y[j][i], k) > 0.5);

            assert(y_ij || y_ji);
        }
    }
#endif

    // Sort items according to y_ij
    std::sort(std::begin(sol),
              std::end(sol),
              [this, &prob, k](auto i, auto j) {
                  return !prob.conflict_matrix(j, i) &&
                         (_cplex.getValue(_y[j][i], k) > 0.5);
              });

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