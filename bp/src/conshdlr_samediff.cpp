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

#include "conshdlr_samediff.hpp"
#include "logging.hpp"
#include "pprint.hpp"
#include "probdata.hpp"
#include "vardata.hpp"
#include <algorithm>
#include <cassert>

// Print all constraint handler events
constexpr bool PrintEvents = false;

using namespace pslp::bp;

const char *ConshdlrSamediff::Properties::Name = "samediff";
const char *ConshdlrSamediff::Properties::Desc = "same/diff constraint handler for the PSLP";
const int ConshdlrSamediff::Properties::Sepapriority = 0;
const int ConshdlrSamediff::Properties::Enfopriority = 0;
const int ConshdlrSamediff::Properties::Checkpriority = 2000000;
const int ConshdlrSamediff::Properties::Sepafreq = -1;
const int ConshdlrSamediff::Properties::Propfreq = 1;
const int ConshdlrSamediff::Properties::Eagerfreq = 1;
const int ConshdlrSamediff::Properties::Maxprerounds = 0;
const SCIP_Bool ConshdlrSamediff::Properties::Delaysepa = false;
const SCIP_Bool ConshdlrSamediff::Properties::Delayprop = false;
const SCIP_Bool ConshdlrSamediff::Properties::Needscons = true;
const SCIP_PROPTIMING ConshdlrSamediff::Properties::Proptiming = SCIP_PROPTIMING_BEFORELP;
const SCIP_PRESOLTIMING ConshdlrSamediff::Properties::Presoltiming = SCIP_PROPTIMING_ALWAYS;

struct SCIP_ConsData
{
    pslp::bp::ConshdlrSamediff::Type type;
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t n_propagated_vars = 0;
    std::size_t n_propagations = 0;
    bool is_propagated = false;
    SCIP_NODE *node = nullptr;
};

ConshdlrSamediff::ConshdlrSamediff(SCIP *scip)
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

SCIP_RETCODE ConshdlrSamediff::create_cons(SCIP *scip,
                                           SCIP_CONS **cons,
                                           const char *name,
                                           Type type,
                                           std::size_t i,
                                           std::size_t j,
                                           SCIP_NODE *node)
{
    SCIP_CONSHDLR *conshdlr = SCIPfindConshdlr(scip, Properties::Name);

    assert(conshdlr != nullptr);

    if constexpr (PrintEvents)
    {
        LOG_DEBUG("create constraint <%s>", name);
    }

    assert(i != j);

    SCIP_CONSDATA *consdata = new SCIP_ConsData;

    consdata->type = type;
    consdata->i = i;
    consdata->j = j;
    consdata->node = node;

    bool local = node != nullptr;

    SCIP_CALL(SCIPcreateCons(scip,
                             cons,
                             name,
                             conshdlr,
                             consdata,
                             false,  // initial
                             false,  // separate
                             false,  // enforce
                             false,  // check
                             true,   // propagate
                             local,  // local
                             false,  // modifiable
                             false,  // dynamic
                             false,  // removable
                             true)); // stickingatnode

    return SCIP_OKAY;
}

ConshdlrSamediff::Consdata ConshdlrSamediff::data(SCIP_CONS *cons)
{
    SCIP_CONSDATA *consdata = SCIPconsGetData(cons);
    assert(consdata != nullptr);
    return {consdata->type, consdata->i, consdata->j};
}

bool ConshdlrSamediff::check(SCIP *scip, SCIP_CONS *cons, bool is_propagated)
{
    SCIP_CONSDATA *consdata = SCIPconsGetData(cons);

    assert(consdata != nullptr);

    const auto &vars = Probdata::get(scip).vars();
    const auto n_vars = is_propagated
                            ? std::size(vars)
                            : consdata->n_propagated_vars;

    assert(n_vars <= std::size(vars));

    for (std::size_t v = 0; v < n_vars; ++v)
    {
        auto var = vars[v];

        // Already fixed?
        if (SCIPvarGetUbLocal(var) < 0.5)
        {
            continue;
        }

        const auto &vardata = Vardata::get(scip, var);
        const auto &s = vardata.items();

        bool has_i = std::find(std::begin(s),
                               std::end(s),
                               consdata->i) != std::end(s);

        bool has_j = std::find(std::begin(s),
                               std::end(s),
                               consdata->j) != std::end(s);

        if ((consdata->type == Type::Same && has_i != has_j) ||
            (consdata->type == Type::Diff && has_i && has_j))
        {
            LOG_ERROR("variable <%s> violates constraint <%s>",
                      SCIPvarGetName(var),
                      SCIPconsGetName(cons));
            LOG_ERROR("-> items: %s",
                      pprint_n(std::size(s),
                               [&s](auto i) { return s[i] + 1; })
                          .c_str());
            return false;
        }
    }
    return true;
}

SCIP_RETCODE ConshdlrSamediff::scip_delete([[maybe_unused]] SCIP *scip,
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

SCIP_RETCODE ConshdlrSamediff::scip_trans(SCIP *scip,
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

SCIP_RETCODE ConshdlrSamediff::scip_enfolp([[maybe_unused]] SCIP *scip,
                                           [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                           [[maybe_unused]] SCIP_CONS **conss,
                                           [[maybe_unused]] int nconss,
                                           [[maybe_unused]] int nusefulconss,
                                           [[maybe_unused]] SCIP_Bool solinfeasible,
                                           SCIP_RESULT *result)
{
    *result = SCIP_FEASIBLE;
    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrSamediff::scip_enfops([[maybe_unused]] SCIP *scip,
                                           [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                           [[maybe_unused]] SCIP_CONS **conss,
                                           [[maybe_unused]] int nconss,
                                           [[maybe_unused]] int nusefulconss,
                                           [[maybe_unused]] SCIP_Bool solinfeasible,
                                           [[maybe_unused]] SCIP_Bool objinfeasible,
                                           SCIP_RESULT *result)
{
    *result = SCIP_FEASIBLE;
    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrSamediff::scip_check([[maybe_unused]] SCIP *scip,
                                          [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                          [[maybe_unused]] SCIP_CONS **conss,
                                          [[maybe_unused]] int nconss,
                                          [[maybe_unused]] SCIP_SOL *sol,
                                          [[maybe_unused]] SCIP_Bool checkintegrality,
                                          [[maybe_unused]] SCIP_Bool checklprows,
                                          [[maybe_unused]] SCIP_Bool printreason,
                                          [[maybe_unused]] SCIP_Bool completely,
                                          SCIP_RESULT *result)
{
    *result = SCIP_FEASIBLE;
    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrSamediff::scip_prop(SCIP *scip,
                                         [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                         [[maybe_unused]] SCIP_CONS **conss,
                                         int nconss,
                                         [[maybe_unused]] int nusefulconss,
                                         [[maybe_unused]] int nmarkedconss,
                                         [[maybe_unused]] SCIP_PROPTIMING proptiming,
                                         SCIP_RESULT *result)
{
    assert(result != nullptr);

    const auto &vars = Probdata::get(scip).vars();
    const auto n_vars = std::size(vars);

    *result = SCIP_DIDNOTFIND;

    for (int c = 0; c < nconss; ++c)
    {
        SCIP_CONS *cons = conss[c];
        SCIP_CONSDATA *consdata = SCIPconsGetData(cons);

        assert(consdata != nullptr);
        assert(check(scip, cons, false));

#ifndef NDEBUG
        for (int c2 = 0; c2 < c; ++c2)
        {
            SCIP_CONSDATA *consdata2 = SCIPconsGetData(conss[c2]);
            assert(consdata->type != consdata2->type ||
                   consdata->i != consdata2->i ||
                   consdata->j != consdata2->j);
            assert(consdata->type != consdata2->type ||
                   consdata->i != consdata2->j ||
                   consdata->j != consdata2->i);
        }
#endif

        if (!consdata->is_propagated)
        {
            if constexpr (PrintEvents)
            {
                LOG_TRACE("propagate constraint <%s>", SCIPconsGetName(cons));

                // Fix variables
                LOG_TRACE("check variables %zu to %zu",
                          consdata->n_propagated_vars,
                          n_vars);
            }

            std::size_t n_fixed_vars = 0;
            bool cutoff = false;

            for (std::size_t v = consdata->n_propagated_vars;
                 v < n_vars && !cutoff;
                 ++v)
            {
                SCIP_VAR *var = vars[v];

                if (SCIPvarGetUbLocal(var) < 0.5)
                {
                    continue;
                }

                const auto &vardata = Vardata::get(scip, var);
                const auto &s = vardata.items();

                bool has_i = std::find(std::begin(s),
                                       std::end(s),
                                       consdata->i) != std::end(s);

                bool has_j = std::find(std::begin(s),
                                       std::end(s),
                                       consdata->j) != std::end(s);

                if ((consdata->type == Type::Same && has_i != has_j) ||
                    (consdata->type == Type::Diff && has_i && has_j))
                {
                    SCIP_Bool infeasible;
                    SCIP_Bool fixed;

                    if constexpr (PrintEvents)
                    {
                        LOG_TRACE("fix variable <%s> = 0", SCIPvarGetName(var));
                    }

                    SCIP_CALL(SCIPfixVar(scip, var, 0., &infeasible, &fixed));

                    if (infeasible)
                    {
                        assert(SCIPvarGetLbLocal(var) > 0.5);

                        if constexpr (PrintEvents)
                        {
                            LOG_TRACE("-> cutoff");
                        }

                        cutoff = true;
                    }
                    else
                    {
                        assert(fixed);
                        ++n_fixed_vars;
                    }
                }
            }

            if constexpr (PrintEvents)
            {
                LOG_TRACE("fixed %zu variables locally", n_fixed_vars);
            }

            if (cutoff)
            {
                *result = SCIP_CUTOFF;
            }
            else if (n_fixed_vars > 0)
            {
                *result = SCIP_REDUCEDDOM;
            }

            ++consdata->n_propagations;

            if (*result != SCIP_CUTOFF)
            {
                consdata->is_propagated = true;
                consdata->n_propagated_vars = n_vars;
            }
            else
            {
                break;
            }
        }

        assert(check(scip, cons, true));
    }

    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrSamediff::scip_lock([[maybe_unused]] SCIP *scip,
                                         [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                         [[maybe_unused]] SCIP_CONS *cons,
                                         [[maybe_unused]] SCIP_LOCKTYPE locktype,
                                         [[maybe_unused]] int nlockspos,
                                         [[maybe_unused]] int nlocksneg)
{
    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrSamediff::scip_active(SCIP *scip,
                                           [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                           SCIP_CONS *cons)
{
    SCIP_CONSDATA *consdata = SCIPconsGetData(cons);

    assert(consdata != nullptr);

    const auto n_vars = std::size(Probdata::get(scip).vars());

    assert(consdata->n_propagated_vars <= n_vars);

    if constexpr (PrintEvents)
    {
        if (consdata->node == nullptr)
        {
            LOG_TRACE("activate constraint <%s>", SCIPconsGetName(cons));
        }
        else
        {
            LOG_TRACE("activate constraint <%s> at node %lld, depth %d",
                      SCIPconsGetName(cons),
                      SCIPnodeGetNumber(consdata->node),
                      SCIPnodeGetDepth(consdata->node));
        }
    }

    if (consdata->n_propagated_vars < n_vars)
    {
        if constexpr (PrintEvents)
        {
            LOG_TRACE("-> mark constraint to be repropagated");
        }

        consdata->is_propagated = false;

        assert(consdata->node != nullptr);

        SCIP_CALL(SCIPrepropagateNode(scip, consdata->node));
    }

    return SCIP_OKAY;
}

SCIP_RETCODE ConshdlrSamediff::scip_deactive([[maybe_unused]] SCIP *scip,
                                             [[maybe_unused]] SCIP_CONSHDLR *conshdlr,
                                             SCIP_CONS *cons)
{
    SCIP_CONSDATA *consdata = SCIPconsGetData(cons);

    assert(consdata != nullptr);

    const auto n_vars = std::size(Probdata::get(scip).vars());

    assert(consdata->is_propagated || SCIPgetNChildren(scip) == 0);

    if constexpr (PrintEvents)
    {
        if (consdata->node == nullptr)
        {
            LOG_TRACE("deactivate constraint <%s>", SCIPconsGetName(cons));
        }
        else
        {
            LOG_TRACE("deactivate constraint <%s> at node %lld, depth %d",
                      SCIPconsGetName(cons),
                      SCIPnodeGetNumber(consdata->node),
                      SCIPnodeGetDepth(consdata->node));
        }
    }

    consdata->n_propagated_vars = n_vars;

    return SCIP_OKAY;
}