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

#include "pricer_mip.hpp"
#include "graph_algorithm.hpp"
#include "logging.hpp"
#include <cassert>
#include <scip/scipdefplugins.h>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace pslp::bp;

PricerMip::PricerMip(SCIP *scip,
                     const char *name,
                     const char *desc,
                     int priority,
                     SCIP_Bool delay)
    : Pricer(scip, name, desc, priority, delay)
{
    SCIPcreate(&_subscip);
    SCIPincludeDefaultPlugins(_subscip);

    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/limits/time";
        SCIPaddRealParam(scip,
                         oss.str().c_str(),
                         "(pslp) maximal time in seconds to run",
                         &_limits_time,
                         false,
                         _limits_time,
                         0.,
                         1.e20,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/limits/memory";
        SCIPaddRealParam(scip,
                         oss.str().c_str(),
                         "(pslp) maximal memory usage in MB",
                         &_limits_memory,
                         false,
                         _limits_memory,
                         0.,
                         SCIP_MEM_NOLIMIT,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/display/verblevel";
        SCIPaddIntParam(scip,
                        oss.str().c_str(),
                        "(pslp) verbosity level of output",
                        &_display_verblevel,
                        false,
                        _display_verblevel,
                        0,
                        5,
                        nullptr,
                        nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/misc/catchctrlc";
        SCIPaddBoolParam(scip,
                         oss.str().c_str(),
                         "(pslp) should the CTRL-C interrupt be caught by the pricer?",
                         &_misc_catchctrlc,
                         false,
                         _misc_catchctrlc,
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
}

PricerMip::~PricerMip()
{
    SCIPfree(&_subscip);
}

SCIP_RETCODE PricerMip::init(SCIP *scip, const Problem &prob)
{
    SCIP_CALL(Pricer::init(scip, prob));

    const auto n = prob.conflict_graph.n_vertices();
    const auto b = prob.capacity;

    // Get the time limit of the master problem
    SCIP_CALL(SCIPgetRealParam(scip, "limits/time", &_global_time_limit));

    SCIP_CALL(SCIPsetIntParam(_subscip, "display/verblevel", _display_verblevel));
    SCIP_CALL(SCIPsetBoolParam(_subscip, "misc/catchctrlc", _misc_catchctrlc));
    SCIP_CALL(SCIPsetRealParam(_subscip, "limits/time", _limits_time));
    SCIP_CALL(SCIPsetRealParam(_subscip, "limits/memory", _limits_memory));

    SCIP_CALL(SCIPcreateProbBasic(_subscip, "mip"));
    SCIP_CALL(SCIPsetObjsense(_subscip, SCIP_OBJSENSE_MAXIMIZE));

    _x.resize(n);

    for (std::size_t i = 0; i < n; ++i)
    {
        SCIP_VAR *var;
        std::ostringstream oss;

        oss << "x{" << i + 1 << '}';

        SCIP_CALL(SCIPcreateVarBasic(_subscip,
                                     &var,
                                     oss.str().c_str(),
                                     0.,
                                     1.,
                                     0.,
                                     SCIP_VARTYPE_BINARY));

        _x[i] = var;

        SCIP_CALL(SCIPaddVar(_subscip, var));
        SCIP_CALL(SCIPreleaseVar(_subscip, &var));
    }

    // Capacity constraint
    if (b < n)
    {
        SCIP_CONS *cons;

        SCIP_CALL(SCIPcreateConsBasicKnapsack(_subscip,
                                              &cons,
                                              "capacity",
                                              0,
                                              nullptr,
                                              nullptr,
                                              b));

        for (std::size_t i = 0; i < n; ++i)
        {
            SCIP_CALL(SCIPaddCoefKnapsack(_subscip,
                                          cons,
                                          _x[i],
                                          1));
        }

        SCIP_CALL(SCIPaddCons(_subscip, cons));
        SCIP_CALL(SCIPreleaseCons(_subscip, &cons));
    }

    return SCIP_OKAY;
}

SCIP_RETCODE PricerMip::exit(SCIP *scip)
{
    SCIP_CALL(SCIPfreeProb(_subscip));

    _conss_same_diff.clear();

    SCIP_CALL(Pricer::exit(scip));

    return SCIP_OKAY;
}

SCIP_RETCODE PricerMip::solve(SCIP *scip,
                              const Problem &prob,
                              std::vector<Solution> &sols)
{
    const auto n = prob.conflict_graph.n_vertices();

    if (prob.has_changed)
    {
        SCIP_CALL(change(scip, prob));
    }

    // Update coefs
    for (std::size_t i = 0; i < n; ++i)
    {
        SCIP_CALL(SCIPchgVarObj(_subscip, _x[i], prob.coefs[i]));
    }

    if (_lpfilename && *_lpfilename)
    {
        std::string filename = make_filename(scip, _lpfilename);
        LOG_DEBUG("write pricing lp in <%s>", filename.c_str());
        SCIP_CALL(SCIPwriteOrigProblem(_subscip,
                                       filename.c_str(),
                                       nullptr,
                                       false));
    }

    // If _extratime = FALSE and master time limit < inf,
    // align the time limit of the pricer on the master
    if (!_extratime && _global_time_limit < 1.e20)
    {
        // _global_time_limit is the value of the global "limits/time" parameter
        SCIP_Real elapsed_time = SCIPgetSolvingTime(scip);
        SCIP_Real remaining_time = _global_time_limit - elapsed_time;

        if (remaining_time <= 0)
        {
            LOG_TRACE("no time remaining (elapsed: %g / %g s)",
                      elapsed_time,
                      _global_time_limit);
            return SCIP_OKAY;
        }

        if (remaining_time < _limits_time)
        {
            LOG_TRACE("set pricing time limit to %g s (elapsed: %g / %g s)",
                      remaining_time,
                      elapsed_time,
                      _global_time_limit);
            SCIP_CALL(SCIPsetRealParam(_subscip,
                                       "limits/time",
                                       remaining_time));
        }
    }

    SCIP_CALL(SCIPsolve(_subscip));

#ifndef NDEBUG
    LOG_DEBUG("MIP pricing solving time: %g s", SCIPgetSolvingTime(_subscip));
#endif

    SCIP_CALL(solutions(scip, prob, sols));

    SCIP_CALL(SCIPfreeTransform(_subscip));

    return SCIP_OKAY;
}

SCIP_RETCODE PricerMip::change([[maybe_unused]] SCIP *scip, const Problem &prob)
{
    for (auto cons : _conss_same_diff)
    {
        SCIP_CALL(SCIPdelCons(_subscip, cons));
    }

    _conss_same_diff.clear();

    // Constraints of type same
    for (auto [i, j] : prob.same)
    {
        SCIP_CONS *cons;
        std::ostringstream oss;

        oss << "same{" << i + 1 << ',' << j + 1 << '}';

        SCIP_CALL(SCIPcreateConsBasicVarbound(_subscip,
                                              &cons,
                                              oss.str().c_str(),
                                              _x[i],
                                              _x[j],
                                              -1.,
                                              0.,
                                              0.));

        _conss_same_diff.push_back(cons);

        SCIP_CALL(SCIPaddCons(_subscip, cons));
        SCIP_CALL(SCIPreleaseCons(_subscip, &cons));
    }

    // Constraints of type diff
    for (auto [i, j] : prob.diff)
    {
        SCIP_CONS *cons;
        std::ostringstream oss;

        oss << "diff{" << i + 1 << ',' << j + 1 << '}';

        SCIP_CALL(SCIPcreateConsBasicVarbound(_subscip,
                                              &cons,
                                              oss.str().c_str(),
                                              _x[i],
                                              _x[j],
                                              1.,
                                              -SCIPinfinity(_subscip),
                                              1.));

        _conss_same_diff.push_back(cons);

        SCIP_CALL(SCIPaddCons(_subscip, cons));
        SCIP_CALL(SCIPreleaseCons(_subscip, &cons));
    }

    return SCIP_OKAY;
}