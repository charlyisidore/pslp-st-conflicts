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

#include "conshdlr_cycle.hpp"
#include "graph_algorithm.hpp"
#include "logging.hpp"
#include "pprint.hpp"
#include "pricer.hpp"
#include <algorithm>
#include <scip/scip_cut.h>
#include <scip/scip_lp.h>
#include <tuple>
#include <vector>

// Print events
constexpr bool PrintEvents = false;

// Print cuts
constexpr bool PrintCuts = false;

using namespace pslp::bp;

const char *ConshdlrCycle::Properties::Name = "cycle";
const char *ConshdlrCycle::Properties::Desc = "cycle elimination constraint handler for MIP pricer";
const int ConshdlrCycle::Properties::Sepapriority = 0;
const int ConshdlrCycle::Properties::Enfopriority = -2000000;
const int ConshdlrCycle::Properties::Checkpriority = -2000000;
const int ConshdlrCycle::Properties::Sepafreq = -1;
const int ConshdlrCycle::Properties::Propfreq = -1;
const int ConshdlrCycle::Properties::Eagerfreq = 1;
const int ConshdlrCycle::Properties::Maxprerounds = 0;
const SCIP_Bool ConshdlrCycle::Properties::Delaysepa = false;
const SCIP_Bool ConshdlrCycle::Properties::Delayprop = false;
const SCIP_Bool ConshdlrCycle::Properties::Needscons = true;
const SCIP_PROPTIMING ConshdlrCycle::Properties::Proptiming = SCIP_PROPTIMING_BEFORELP;
const SCIP_PRESOLTIMING ConshdlrCycle::Properties::Presoltiming = SCIP_PRESOLTIMING_FAST;

struct SCIP_ConsData
{
    SCIP_ConsData(const std::vector<SCIP_VAR *> &vars,
                  const SquareMatrix<bool> &adjacency_matrix,
                  std::size_t max_cuts)
        : vars(vars),
          adjacency_matrix(adjacency_matrix),
          max_cuts(max_cuts)
    {
    }

    std::vector<SCIP_VAR *> vars;
    SquareMatrix<bool> adjacency_matrix;
    std::size_t max_cuts = 0;
};

// For stopping a callback execution if SCIP returns an error
struct SCIPError
{
    SCIPError(SCIP_RETCODE retcode)
        : retcode(retcode) {}
    SCIP_RETCODE retcode;
};

// Macro throwing a SCIPError in case of failure
#define SCIP_THROW(x)                     \
    {                                     \
        SCIP_RETCODE retcode;             \
        if ((retcode = (x)) != SCIP_OKAY) \
        {                                 \
            throw SCIPError(retcode);     \
        }                                 \
    }

[[nodiscard]] static SCIP_RETCODE enfo(SCIP *scip,
                                       SCIP_CONSHDLR *conshdlr,
                                       SCIP_CONS **conss,
                                       int nconss,
                                       SCIP_SOL *sol,
                                       SCIP_RESULT *result)
{
    assert(result != nullptr);

    *result = SCIP_FEASIBLE;

    // Exception thrown when a cutoff is detected
    struct Cutoff
    {
    };

    // Exception thrown when the maximum number of cuts is reached
    struct MaxCuts
    {
    };

    for (int c = 0; c < nconss; ++c)
    {
        assert(conss != nullptr);
        assert(conss[c] != nullptr);

        SCIP_CONSDATA *consdata = SCIPconsGetData(conss[c]);

        assert(consdata != nullptr);

        const auto n_vars = std::size(consdata->vars);
        std::vector<std::size_t> s;

        // Get the selected items
        for (std::size_t i = 0; i < n_vars; ++i)
        {
            auto value = SCIPgetSolVal(scip, sol, consdata->vars[i]);
            if (value > 0.5)
            {
                assert(SCIPisFeasEQ(scip, value, 1.));
                s.push_back(i);
            }
            else
            {
                assert(SCIPisFeasEQ(scip, value, 0.));
            }
        }

        // This constraint cannot be violated if zero or one item is selected
        if (std::size(s) < 2)
        {
            continue;
        }

        // Count the number of cuts added
        std::size_t n_cuts = 0;

        // This callback is executed for every cycle detected
        // Creates a constraint to avoid selecting all items in [first,last)
        auto create_cut = [&](auto first, auto last) {
            auto cycle_length = std::distance(first, last);

#ifndef NDEBUG
            // Check whether it is a valid cycle
            assert(cycle_length >= 2);
            for (auto i = first; i != last; ++i)
            {
                assert((i + 1 == last) ||
                       consdata->adjacency_matrix(s[*i], s[*(i + 1)]));
                assert((i + 1 != last) ||
                       consdata->adjacency_matrix(s[*i], s[*first]));
            }
#endif

            SCIP_ROW *row;

            SCIP_THROW(SCIPcreateEmptyRowConshdlr(scip,
                                                  &row,
                                                  conshdlr,
                                                  "c",
                                                  -SCIPinfinity(scip),
                                                  SCIP_Real(cycle_length - 1),
                                                  false,
                                                  false,
                                                  true));

            SCIP_THROW(SCIPcacheRowExtensions(scip, row));

            for (auto i = first; i != last; ++i)
            {
                auto var = consdata->vars[s[*i]];
                SCIP_THROW(SCIPaddVarToRow(scip, row, var, 1.));
            }

            SCIP_THROW(SCIPflushRowExtensions(scip, row));

            if constexpr (PrintCuts)
            {
                LOG_TRACE("enfo cut:");
                SCIP_THROW(SCIPprintRow(scip, row, nullptr));
            }

            SCIP_Bool infeasible;
            SCIP_THROW(SCIPaddRow(scip, row, false, &infeasible));
            SCIP_THROW(SCIPreleaseRow(scip, &row));

            if (infeasible)
            {
                throw Cutoff();
            }

            if (consdata->max_cuts > 0 && ++n_cuts >= consdata->max_cuts)
            {
                throw MaxCuts();
            }
        };

        // Build a subgraph from selected items
        auto subgraph = Pricer::subgraph(consdata->adjacency_matrix, s);

        // Find cycles
        try
        {
            CycleFinder cycle_finder(create_cut);
            depth_first_visit(subgraph, cycle_finder);
        }
        catch (const MaxCuts &)
        {
            // Reached the maximum number of cycles allowed per round
        }
        catch (const Cutoff &)
        {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
        }
        catch (const SCIPError &e)
        {
            return e.retcode;
        }

        if (n_cuts > 0)
        {
            *result = SCIP_SEPARATED;
        }
    }

    return SCIP_OKAY;
}

[[nodiscard]] static SCIP_RETCODE check(SCIP *scip,
                                        SCIP_CONS **conss,
                                        int nconss,
                                        SCIP_SOL *sol,
                                        SCIP_Bool printreason,
                                        SCIP_RESULT *result)
{
    assert(result != nullptr);

    *result = SCIP_FEASIBLE;

    // Exception thrown when a cycle is detected
    struct HasCycle
    {
    };

    for (int c = 0; c < nconss; ++c)
    {
        assert(conss != nullptr);
        assert(conss[c] != nullptr);

        SCIP_CONSDATA *consdata = SCIPconsGetData(conss[c]);

        assert(consdata != nullptr);

        const auto n_vars = std::size(consdata->vars);
        std::vector<std::size_t> s;

        // Get the selected items
        for (std::size_t i = 0; i < n_vars; ++i)
        {
            auto value = SCIPgetSolVal(scip, sol, consdata->vars[i]);
            if (value > 0.5)
            {
                assert(SCIPisFeasEQ(scip, value, 1.));
                s.push_back(i);
            }
            else
            {
                assert(SCIPisFeasEQ(scip, value, 0.));
            }
        }

        // This constraint cannot be violated if zero or one item is selected
        if (std::size(s) < 2)
        {
            continue;
        }

        // Build a subgraph from selected items
        auto subgraph = Pricer::subgraph(consdata->adjacency_matrix, s);

        // We just want to know if there exists a cycle
        struct
        {
            void back_edge([[maybe_unused]] Edge e,
                           [[maybe_unused]] const decltype(subgraph) &g)
            {
                throw HasCycle();
            }
        } cycle_detector;

        // Find a single cycle
        try
        {
            depth_first_visit(subgraph, cycle_detector);
        }
        catch (const HasCycle &)
        {
            *result = SCIP_INFEASIBLE;

            if (printreason)
            {
                SCIP_CALL(SCIPprintCons(scip, conss[c], nullptr));
                SCIPinfoMessage(scip, nullptr, "violation: cycle");
            }
        }
        catch (const SCIPError &e)
        {
            return e.retcode;
        }
    }

    return SCIP_OKAY;
}

ConshdlrCycle::ConshdlrCycle(SCIP *scip)
    : scip::ObjConshdlr(scip,
                        Properties::Name,
                        Properties::Desc,
                        Properties::Sepapriority,
                        Properties::Enfopriority,
                        Properties::Checkpriority,
                        Properties::Sepafreq,
                        Properties::Propfreq,
                        Properties::Eagerfreq,
                        Properties::Maxprerounds,
                        Properties::Delaysepa,
                        Properties::Delayprop,
                        Properties::Needscons,
                        Properties::Proptiming,
                        Properties::Presoltiming)
{
}

SCIP_RETCODE ConshdlrCycle::create_cons(SCIP *scip,
                                        SCIP_CONS **cons,
                                        const char *name,
                                        const std::vector<SCIP_VAR *> &vars,
                                        const SquareMatrix<bool> &adjacency_matrix,
                                        std::size_t max_cuts)
{
    SCIP_CONSHDLR *conshdlr = SCIPfindConshdlr(scip, Properties::Name);

    assert(conshdlr != nullptr);

    if constexpr (PrintEvents)
    {
        LOG_DEBUG("create constraint <%s>", name);
    }

    SCIP_CONSDATA *consdata = new SCIP_ConsData(vars,
                                                adjacency_matrix,
                                                max_cuts);

    SCIP_CALL(SCIPcreateCons(scip,
                             cons,
                             name,
                             conshdlr,
                             consdata,
                             false,   // initial
                             true,    // separate
                             true,    // enforce
                             true,    // check
                             true,    // propagate
                             false,   // local
                             false,   // modifiable
                             false,   // dynamic
                             true,    // removable
                             false)); // stickingatnode

    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrCycle::scip_delete([[maybe_unused]] SCIP *scip,
                                        [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                        [[maybe_unused]] SCIP_CONS *cons,
                                        SCIP_CONSDATA **consdata)
{
    assert(consdata != nullptr);
    assert(*consdata != nullptr);

    if constexpr (PrintEvents)
    {
        LOG_TRACE("delete constraint <%s>", SCIPconsGetName(cons));
    }

    delete *consdata;
    *consdata = nullptr;

    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrCycle::scip_trans(SCIP *scip,
                                       [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                       SCIP_CONS *sourcecons,
                                       SCIP_CONS **targetcons)
{
    assert(conshdlr != nullptr);
    assert(sourcecons != nullptr);
    assert(targetcons != nullptr);

    if constexpr (PrintEvents)
    {
        LOG_TRACE("transform constraint <%s>", SCIPconsGetName(sourcecons));
    }

    SCIP_CONSDATA *sourcedata = SCIPconsGetData(sourcecons);

    assert(sourcedata != nullptr);

    SCIP_CONSDATA *targetdata = new SCIP_ConsData(*sourcedata);

    SCIP_CALL(SCIPcreateCons(scip,
                             targetcons,
                             SCIPconsGetName(sourcecons),
                             conshdlr,
                             targetdata,
                             SCIPconsIsInitial(sourcecons),
                             SCIPconsIsSeparated(sourcecons),
                             SCIPconsIsEnforced(sourcecons),
                             SCIPconsIsChecked(sourcecons),
                             SCIPconsIsPropagated(sourcecons),
                             SCIPconsIsLocal(sourcecons),
                             SCIPconsIsModifiable(sourcecons),
                             SCIPconsIsDynamic(sourcecons),
                             SCIPconsIsRemovable(sourcecons),
                             SCIPconsIsStickingAtNode(sourcecons)));

    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrCycle::scip_enfolp(SCIP *scip,
                                        SCIP_CONSHDLR *conshdlr,
                                        SCIP_CONS **conss,
                                        int nconss,
                                        [[maybe_unused]] int nusefulconss,
                                        [[maybe_unused]] SCIP_Bool solinfeasible,
                                        SCIP_RESULT *result)
{
    SCIP_CALL(enfo(scip,
                   conshdlr,
                   conss,
                   nconss,
                   nullptr,
                   result));
    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrCycle::scip_enfops(SCIP *scip,
                                        SCIP_CONSHDLR *conshdlr,
                                        SCIP_CONS **conss,
                                        int nconss,
                                        [[maybe_unused]] int nusefulconss,
                                        [[maybe_unused]] SCIP_Bool solinfeasible,
                                        [[maybe_unused]] SCIP_Bool objinfeasible,
                                        SCIP_RESULT *result)
{
    SCIP_CALL(enfo(scip,
                   conshdlr,
                   conss,
                   nconss,
                   nullptr,
                   result));
    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrCycle::scip_check(SCIP *scip,
                                       [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                       SCIP_CONS **conss,
                                       int nconss,
                                       SCIP_SOL *sol,
                                       [[maybe_unused]] SCIP_Bool checkintegrality,
                                       [[maybe_unused]] SCIP_Bool checklprows,
                                       SCIP_Bool printreason,
                                       [[maybe_unused]] SCIP_Bool completely,
                                       SCIP_RESULT *result)
{
    SCIP_CALL(check(scip,
                    conss,
                    nconss,
                    sol,
                    printreason,
                    result));
    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrCycle::scip_prop([[maybe_unused]] SCIP *scip,
                                      [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                      [[maybe_unused]] SCIP_CONS **conss,
                                      [[maybe_unused]] int nconss,
                                      [[maybe_unused]] int nusefulconss,
                                      [[maybe_unused]] int nmarkedconss,
                                      [[maybe_unused]] SCIP_PROPTIMING proptiming,
                                      SCIP_RESULT *result)
{
    assert(result != nullptr);
    *result = SCIP_DIDNOTRUN;
    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrCycle::scip_lock(SCIP *scip,
                                      [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                      SCIP_CONS *cons,
                                      SCIP_LOCKTYPE locktype,
                                      int nlockspos,
                                      int nlocksneg)
{
    SCIP_CONSDATA *consdata = SCIPconsGetData(cons);

    assert(consdata != nullptr);

    for (auto var : consdata->vars)
    {
        SCIP_CALL(SCIPaddVarLocksType(scip,
                                      var,
                                      locktype,
                                      nlocksneg,
                                      nlockspos));
    }

    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrCycle::scip_print(SCIP *scip,
                                       [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                       [[maybe_unused]] SCIP_CONS *cons,
                                       FILE *file)
{
    SCIPinfoMessage(scip, file, "cycle elimination constraint\n");

    return SCIP_OKAY;
}