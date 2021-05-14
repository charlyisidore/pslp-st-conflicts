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

#include "pricer_cplex.hpp"
#include "logging.hpp"
#include <cassert>
#include <chrono>
#include <sstream>
#include <string>

// Print CPLEX log
constexpr bool PrintCplexLog = false;

// Print CPLEX warnings
constexpr bool PrintCplexWarnings = false;

using namespace pslp::bp;

PricerCplex::PricerCplex(SCIP *scip,
                         const char *name,
                         const char *desc,
                         int priority,
                         SCIP_Bool delay)
    : Pricer(scip, name, desc, priority, delay)
{
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/threads";
        SCIPaddIntParam(scip,
                        oss.str().c_str(),
                        "(pslp) global thread count (0: auto)",
                        &_threads,
                        false,
                        _threads,
                        0,
                        IntMax,
                        nullptr,
                        nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/timelimit";
        SCIPaddRealParam(scip,
                         oss.str().c_str(),
                         "(pslp) optimizer time limit in seconds",
                         &_timelimit,
                         false,
                         _timelimit,
                         0.,
                         RealMax,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/workmem";
        SCIPaddRealParam(scip,
                         oss.str().c_str(),
                         "(pslp) memory available for working storage",
                         &_workmem,
                         false,
                         _workmem,
                         0.,
                         RealMax,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/mip/strategy/file";
        SCIPaddIntParam(scip,
                        oss.str().c_str(),
                        "(pslp) node storage file switch (0: no node file, 1: in memory and compressed, 2: on disk, 3: on disk and compressed)",
                        &_mip_strategy_file,
                        false,
                        _mip_strategy_file,
                        0,
                        3,
                        nullptr,
                        nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/mip/limits/treememory";
        SCIPaddRealParam(scip,
                         oss.str().c_str(),
                         "(pslp) tree memory limit",
                         &_mip_limits_treememory,
                         false,
                         _mip_limits_treememory,
                         0.,
                         RealMax,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/mip/limits/solutions";
        SCIPaddLongintParam(scip,
                            oss.str().c_str(),
                            "(pslp) MIP integer solution limit",
                            &_mip_limits_solutions,
                            false,
                            _mip_limits_solutions,
                            1,
                            LongintMax,
                            nullptr,
                            nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/simplex/limits/upperobj";
        SCIPaddRealParam(scip,
                         oss.str().c_str(),
                         "(pslp) upper objective value limit",
                         &_simplex_limits_upperobj,
                         false,
                         _simplex_limits_upperobj,
                         RealMin,
                         RealMax,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/display/out";
        SCIPaddBoolParam(scip,
                         oss.str().c_str(),
                         "(pslp) display pricer output log?",
                         &_display_out,
                         false,
                         _display_out,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/display/warning";
        SCIPaddBoolParam(scip,
                         oss.str().c_str(),
                         "(pslp) display pricer warnings?",
                         &_display_warning,
                         false,
                         _display_warning,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/lpfilename";
        SCIPaddStringParam(scip,
                           oss.str().c_str(),
                           "(pslp) file name to export lp files",
                           &_lpfilename,
                           false,
                           "",
                           nullptr,
                           nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/extratime";
        SCIPaddBoolParam(scip,
                         oss.str().c_str(),
                         "(pslp) can the pricer exceed the global time limit?",
                         &_extratime,
                         false,
                         _extratime,
                         nullptr,
                         nullptr);
    }
}

PricerCplex::~PricerCplex()
{
    _env.end();
}

SCIP_RETCODE PricerCplex::init(SCIP *scip, const Problem &prob)
{
    SCIP_CALL(Pricer::init(scip, prob));

    const auto n = prob.conflict_graph.n_vertices();
    const auto b = prob.capacity;

    _env.end();
    _env = IloEnv();

    _cplex = IloCplex(_env);
    _model = IloModel(_env);

    if (!_display_out)
    {
        _cplex.setOut(_env.getNullStream());
    }

    if (!_display_warning)
    {
        _cplex.setWarning(_env.getNullStream());
    }

    // Get the global time limit of SCIP
    SCIP_CALL(SCIPgetRealParam(scip, "limits/time", &_global_time_limit));

    // Settings
    _cplex.setParam(IloCplex::Param::Threads, _threads);
    _cplex.setParam(IloCplex::Param::TimeLimit, _timelimit);
    _cplex.setParam(IloCplex::Param::WorkMem, _workmem);
    _cplex.setParam(IloCplex::Param::MIP::Strategy::File, _mip_strategy_file);
    _cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, _mip_limits_treememory);
    _cplex.setParam(IloCplex::Param::MIP::Limits::Solutions, _mip_limits_solutions);
    _cplex.setParam(IloCplex::Param::Simplex::Limits::UpperObj, _simplex_limits_upperobj);

    _x = IloNumVarArray(_env, n, 0, 1, IloNumVar::Bool);
    for (std::size_t i = 0; i < n; ++i)
    {
        std::ostringstream oss;
        oss << "x{" << i + 1 << '}';
        _x[i].setName(oss.str().c_str());
    }
    _model.add(_x);

    // Minimize reduced cost
    _objective = IloMaximize(_env);
    _model.add(_objective);

    // Capacity constraint
    if (b < n)
    {
        IloRange cons(_env, -IloInfinity, b, "capacity");
        for (std::size_t i = 0; i < n; ++i)
        {
            cons.setLinearCoef(_x[i], 1.);
        }
        _model.add(cons);
    }

    _conss_same_diff = IloRangeArray(_env);

    return SCIP_OKAY;
}

SCIP_RETCODE PricerCplex::exit(SCIP *scip)
{
    _env.end();
    _env = IloEnv();

    SCIP_CALL(Pricer::exit(scip));

    return SCIP_OKAY;
}

SCIP_RETCODE PricerCplex::solve(SCIP *scip,
                                const Problem &prob,
                                std::vector<Solution> &sols)
{
#ifndef NDEBUG
    using Clock = std::chrono::steady_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    using Duration = std::chrono::duration<double>;
#endif

    const auto n = prob.conflict_graph.n_vertices();

    try
    {
        if (prob.has_changed)
        {
            SCIP_CALL(change(scip, prob));
        }

        // Update coefs
        for (std::size_t i = 0; i < n; ++i)
        {
            _objective.setLinearCoef(_x[i], prob.coefs[i]);
        }
        _objective.setConstant(prob.offset);

        _cplex.extract(_model);

        if (_lpfilename && *_lpfilename)
        {
            std::string filename = make_filename(scip, _lpfilename);
            LOG_DEBUG("write pricing lp in <%s>", filename.c_str());
            _cplex.exportModel(filename.c_str());
        }

        // If _extratime = FALSE and SCIP time limit < inf,
        // align the time limit of CPLEX on SCIP
        if (!_extratime && _global_time_limit < 1.e20)
        {
            // _global_time_limit is the value of the "limits/time" parameter
            SCIP_Real elapsed_time = SCIPgetSolvingTime(scip);
            SCIP_Real remaining_time = _global_time_limit - elapsed_time;

            if (remaining_time <= 0)
            {
                LOG_TRACE("no time remaining (elapsed: %g / %g s)",
                          elapsed_time,
                          _global_time_limit);
                return SCIP_OKAY;
            }

            if (remaining_time < _timelimit)
            {
                LOG_TRACE("set pricing time limit to %g s (elapsed: %g / %g s)",
                          remaining_time,
                          elapsed_time,
                          _global_time_limit);
                _cplex.setParam(IloCplex::Param::TimeLimit, remaining_time);
            }
        }

#ifndef NDEBUG
        TimePoint start_time = Clock::now();
#endif

        _cplex.solve();

#ifndef NDEBUG
        TimePoint finish_time = Clock::now();
        double solving_time = Duration(finish_time - start_time).count();

        LOG_DEBUG("CPLEX pricing solving time: %g s", solving_time);
#endif

        std::size_t n_sols = _cplex.getSolnPoolNsolns();

        for (std::size_t k = 0; k < n_sols; ++k)
        {
            SCIP_CALL(solution(scip, prob, k, sols));
        }
    }
    catch (const IloException &e)
    {
        LOG_ERROR("%s", e.getMessage());
        return SCIP_ERROR;
    }

    return SCIP_OKAY;
}

SCIP_RETCODE PricerCplex::change([[maybe_unused]] SCIP *scip,
                                 const Problem &prob)
{
    _model.remove(_conss_same_diff);
    _conss_same_diff.clear();

    // Constraints of type same
    for (auto [i, j] : prob.same)
    {
        std::ostringstream oss;
        oss << "same{" << i + 1 << ',' << j + 1 << '}';

        IloRange cons(_env, 0., 0., oss.str().c_str());
        cons.setLinearCoef(_x[i], 1.);
        cons.setLinearCoef(_x[j], -1.);

        _conss_same_diff.add(cons);
    }

    // Constraints of type diff
    for (auto [i, j] : prob.diff)
    {
        std::ostringstream oss;
        oss << "diff{" << i + 1 << ',' << j + 1 << '}';

        IloRange cons(_env, -IloInfinity, 1., oss.str().c_str());
        cons.setLinearCoef(_x[i], 1.);
        cons.setLinearCoef(_x[j], 1.);
        _conss_same_diff.add(cons);
    }

    _model.add(_conss_same_diff);

    return SCIP_OKAY;
}
