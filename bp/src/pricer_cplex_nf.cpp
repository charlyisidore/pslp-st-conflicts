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

#include "pricer_cplex_nf.hpp"
#include "logging.hpp"
#include "pprint.hpp"
#include <cassert>
#include <sstream>

// Use the valid inequality?
constexpr bool UseValidInequality = false;

// Print the solution
constexpr bool PrintSolution = false;

// Print CPLEX variables
constexpr bool PrintCplexVariables = false;

using namespace pslp::bp;

const char *PricerCplexNf::Properties::Name = "cplex_nf";
const char *PricerCplexNf::Properties::Desc = "CPLEX variable pricer with network flow formulation";
const int PricerCplexNf::Properties::Priority = 0;
const SCIP_Bool PricerCplexNf::Properties::Delay = true;

PricerCplexNf::PricerCplexNf(SCIP *scip)
    : PricerCplex(scip,
                  Properties::Name,
                  Properties::Desc,
                  Properties::Priority,
                  Properties::Delay)
{
}

SCIP_RETCODE PricerCplexNf::init(SCIP *scip, const Problem &prob)
{
    if (!prob.has_adjacent_conflicts)
    {
        LOG_ERROR("distant conflicts not supported");
        return SCIP_ERROR;
    }

    SCIP_CALL(PricerCplex::init(scip, prob));

    const auto n = prob.conflict_graph.n_vertices();
    const auto b = prob.capacity;

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

    _s = IloNumVarArray(_env, n, 0, 1, IloNumVar::Bool);
    for (std::size_t i = 0; i < n; ++i)
    {
        std::ostringstream oss;
        oss << "s{" << i + 1 << '}';
        _s[i].setName(oss.str().c_str());
    }
    _model.add(_s);

    _t = IloNumVarArray(_env, n, 0, 1, IloNumVar::Bool);
    for (std::size_t i = 0; i < n; ++i)
    {
        std::ostringstream oss;
        oss << "t{" << i + 1 << '}';
        _t[i].setName(oss.str().c_str());
    }
    _model.add(_t);

    _l = IloNumVarArray(_env, n, 0, b - 1, IloNumVar::Float);
    for (std::size_t i = 0; i < n; ++i)
    {
        std::ostringstream oss;
        oss << "l{" << i + 1 << '}';
        _l[i].setName(oss.str().c_str());
    }
    _model.add(_l);

    // Source constraint: sum s_i = 1
    {
        IloRange cons(_env, 1., 1., "source");
        for (std::size_t i = 0; i < n; ++i)
        {
            cons.setLinearCoef(_s[i], 1.);
        }
        _model.add(cons);
    }

    // Sink constraint: sum t_i = 1
    {
        IloRange cons(_env, 1., 1., "sink");
        for (std::size_t i = 0; i < n; ++i)
        {
            cons.setLinearCoef(_t[i], 1.);
        }
        _model.add(cons);
    }

    // In flow constraints: x_i = s_i + sum_j y_ji
    for (std::size_t i = 0; i < n; ++i)
    {
        std::ostringstream oss;
        oss << "in{" << i + 1 << '}';

        IloRange cons(_env, 0., 0., oss.str().c_str());
        cons.setLinearCoef(_x[i], 1.);
        cons.setLinearCoef(_s[i], -1.);

        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j || prob.conflict_matrix(j, i))
            {
                continue;
            }

            cons.setLinearCoef(_y[j][i], -1.);
        }
        _model.add(cons);
    }

    // Out flow constraints:
    // x_i = t_i + sum_j y_ij
    for (std::size_t i = 0; i < n; ++i)
    {
        std::ostringstream oss;
        oss << "out{" << i + 1 << '}';

        IloRange cons(_env, 0., 0., oss.str().c_str());
        cons.setLinearCoef(_x[i], 1.);
        cons.setLinearCoef(_t[i], -1.);

        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j || prob.conflict_matrix(i, j))
            {
                continue;
            }

            cons.setLinearCoef(_y[i][j], -1.);
        }
        _model.add(cons);
    }

    // Level constraint: l_i + b (1 - y_ij) >= l_j + 1
    //                => l_i - l_j - b y_ij >= 1 - b
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j || prob.conflict_matrix(i, j))
            {
                continue;
            }

            std::ostringstream oss;
            oss << "level{" << i + 1 << ',' << j + 1 << '}';

            IloRange cons(_env, 1. - IloNum(b), IloInfinity, oss.str().c_str());
            cons.setLinearCoef(_l[i], 1.);
            cons.setLinearCoef(_l[j], -1.);
            cons.setLinearCoef(_y[i][j], -IloNum(b));
            _model.add(cons);
        }
    }

    // Valid inequality: sum a_i = sum y_ij + 1
    if constexpr (UseValidInequality)
    {
        IloRange cons(_env, 1., 1., "path");

        for (std::size_t i = 0; i < n; ++i)
        {
            cons.setLinearCoef(_x[i], 1.);

            for (std::size_t j = 0; j < n; ++j)
            {
                if (i == j || prob.conflict_matrix(i, j))
                {
                    continue;
                }

                cons.setLinearCoef(_y[i][j], -1.);
            }
        }
        _model.add(cons);
    }

    return SCIP_OKAY;
}

SCIP_RETCODE PricerCplexNf::solution(SCIP *scip,
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

    LOG_DEBUG("found CPLEX NF pricing sol (obj: %g)", obj);

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
        for (IloNum i = 0; i < _s.getSize(); ++i)
        {
            double value = _cplex.getValue(_s[i], k);
            if (!SCIPisEQ(scip, value, 0.))
            {
                LOG_DEBUG("%s  %g", _s[i].getName(), value);
            }
        }
        for (IloNum i = 0; i < _t.getSize(); ++i)
        {
            double value = _cplex.getValue(_t[i], k);
            if (!SCIPisEQ(scip, value, 0.))
            {
                LOG_DEBUG("%s  %g", _t[i].getName(), value);
            }
        }
        for (IloNum i = 0; i < _l.getSize(); ++i)
        {
            double value = _cplex.getValue(_l[i], k);
            if (!SCIPisEQ(scip, value, 0.))
            {
                LOG_DEBUG("%s  %g", _l[i].getName(), value);
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

    // Check that some constraints were satisfied
#ifndef NDEBUG
    for (std::size_t i = 0; i < n; ++i)
    {
        bool x_i = _cplex.getValue(_x[i], k) > 0.5;
        bool s_i = _cplex.getValue(_s[i], k) > 0.5;
        bool t_i = _cplex.getValue(_t[i], k) > 0.5;
        double l_i = _cplex.getValue(_l[i], k);

        assert(!s_i || x_i);
        assert(!t_i || x_i);

        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j || prob.conflict_matrix(i, j))
            {
                continue;
            }

            bool x_j = _cplex.getValue(_x[j], k) > 0.5;
            double l_j = _cplex.getValue(_l[j], k);
            bool y_ij = _cplex.getValue(_y[i][j], k) > 0.5;

            assert(!y_ij || (x_i && x_j && l_i > l_j));
        }
    }
#endif

    // Sort items according to l_i
    std::sort(std::begin(sol),
              std::end(sol),
              [this, k](auto i, auto j) {
                  double l_i = _cplex.getValue(_l[i], k);
                  double l_j = _cplex.getValue(_l[j], k);
                  return l_i < l_j;
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