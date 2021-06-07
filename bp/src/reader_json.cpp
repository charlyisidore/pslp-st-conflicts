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

#include "reader_json.hpp"
#include "logging.hpp"
#include "probdata.hpp"
#include "vardata.hpp"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <vector>

using namespace pslp::bp;

const char *ReaderJson::Properties::Name = "jsonreader";
const char *ReaderJson::Properties::Desc = "file reader for the PSLP";
const char *ReaderJson::Properties::Extension = "json";

static const char *to_string(SCIP_STATUS status)
{
    switch (status)
    {
        case SCIP_STATUS_UNKNOWN:
            return "UNKNOWN";
        case SCIP_STATUS_USERINTERRUPT:
            return "USERINTERRUPT";
        case SCIP_STATUS_NODELIMIT:
            return "NODELIMIT";
        case SCIP_STATUS_TOTALNODELIMIT:
            return "TOTALNODELIMIT";
        case SCIP_STATUS_STALLNODELIMIT:
            return "STALLNODELIMIT";
        case SCIP_STATUS_TIMELIMIT:
            return "TIMELIMIT";
        case SCIP_STATUS_MEMLIMIT:
            return "MEMLIMIT";
        case SCIP_STATUS_GAPLIMIT:
            return "GAPLIMIT";
        case SCIP_STATUS_SOLLIMIT:
            return "SOLLIMIT";
        case SCIP_STATUS_BESTSOLLIMIT:
            return "BESTSOLLIMIT";
        case SCIP_STATUS_RESTARTLIMIT:
            return "RESTARTLIMIT";
        case SCIP_STATUS_OPTIMAL:
            return "OPTIMAL";
        case SCIP_STATUS_INFEASIBLE:
            return "INFEASIBLE";
        case SCIP_STATUS_UNBOUNDED:
            return "UNBOUNDED";
        case SCIP_STATUS_INFORUNBD:
            return "INFORUNBD";
        case SCIP_STATUS_TERMINATE:
            return "TERMINATE";
        default:
            return "?";
    }
}

static const char *to_string(SCIP_SOLTYPE type)
{
    switch (type)
    {
        case SCIP_SOLTYPE_UNKNOWN:
            return "UNKNOWN";
        case SCIP_SOLTYPE_HEUR:
            return "HEUR";
        case SCIP_SOLTYPE_RELAX:
            return "RELAX";
        case SCIP_SOLTYPE_LPRELAX:
            return "LPRELAX";
        case SCIP_SOLTYPE_STRONGBRANCH:
            return "STRONGBRANCH";
        case SCIP_SOLTYPE_PSEUDO:
            return "PSEUDO";
        default:
            return "?";
    }
}

ReaderJson::ReaderJson(SCIP *scip)
    : scip::ObjReader(scip,
                      Properties::Name,
                      Properties::Desc,
                      Properties::Extension)
{
    {
        std::ostringstream oss;
        oss << "reading/" << Properties::Name << "/useadjacentconflicts";
        SCIPaddBoolParam(scip,
                         oss.str().c_str(),
                         "(pslp) use adjacent conflicts when unspecified?",
                         &_useadjacentconflicts,
                         false,
                         _useadjacentconflicts,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "reading/" << Properties::Name << "/transformconflicts";
        SCIPaddBoolParam(scip,
                         oss.str().c_str(),
                         "(pslp) transform conflicts to remove equivalent items?",
                         &_transformconflicts,
                         false,
                         _transformconflicts,
                         nullptr,
                         nullptr);
    }
}

SCIP_RETCODE ReaderJson::scip_read(SCIP *scip,
                                   [[maybe_unused]] SCIP_READER *reader,
                                   const char *filename,
                                   SCIP_RESULT *result)
{
    assert(result != nullptr);

    std::ifstream file(filename);
    Probdata::Type type = Probdata::Type::MinNStacks;
    std::size_t n_items = 0;
    std::size_t n_stacks = 0;
    std::size_t capacity = 0;
    bool has_adjacent_conflicts = _useadjacentconflicts;
    nlohmann::json json;

    if (!file.is_open())
    {
        LOG_ERROR("cannot open file <%s>", filename);
        return SCIP_READERROR;
    }

    file >> json;

    auto m_type = json.find("type");
    auto m_n_items = json.find("n_items");
    auto m_n_stacks = json.find("n_stacks");
    auto m_capacity = json.find("capacity");
    auto m_arrivals = json.find("arrivals");
    auto m_departures = json.find("departures");
    auto m_blocking_graph = json.find("blocking_graph");
    auto m_blocking_matrix = json.find("blocking_matrix");
    auto m_weights = json.find("weights");
    auto m_conflict_graph = json.find("conflict_graph");
    auto m_conflict_matrix = json.find("conflict_matrix");
    auto m_initial_positions = json.find("initial_positions");

    // "n_items" is mandatory
    if (m_n_items == json.end())
    {
        LOG_ERROR("missing value: n_items");
        return SCIP_READERROR;
    }

    // minimize St by default
    if (m_type != json.end())
    {
        if (m_type->get<std::string>() == "St")
        {
            type = Probdata::Type::MinNStacks;
        }
        else if (m_type->get<std::string>() == "BI")
        {
            type = Probdata::Type::MinNBlockingItems;
        }
        else
        {
            LOG_ERROR("wrong value: type");
            return SCIP_READERROR;
        }
    }

    n_items = m_n_items->get<std::size_t>();

    // n_stacks = n_items by default
    if (m_n_stacks == json.end())
    {
        n_stacks = n_items;
    }
    else
    {
        n_stacks = m_n_stacks->get<std::size_t>();
    }

    // capacity = n_items by default
    if (m_capacity == json.end())
    {
        capacity = n_items;
    }
    else
    {
        capacity = m_capacity->get<std::size_t>();
    }

    // Quick bound checking

    if (m_arrivals != json.end() && m_arrivals->size() != n_items)
    {
        LOG_ERROR("wrong size: arrivals");
        return SCIP_READERROR;
    }

    if (m_weights != json.end() && m_weights->size() != n_items)
    {
        LOG_ERROR("wrong size: weights");
        return SCIP_READERROR;
    }

    if (m_departures != json.end() && m_departures->size() != n_items)
    {
        LOG_ERROR("wrong size: departures");
        return SCIP_READERROR;
    }

    if (m_initial_positions != json.end() && m_initial_positions->size() > n_items)
    {
        LOG_ERROR("wrong size: initial_positions");
        return SCIP_READERROR;
    }

    if (m_blocking_matrix != json.end() && m_blocking_matrix->size() != n_items)
    {
        LOG_ERROR("wrong size: blocking_matrix");
        return SCIP_READERROR;
    }

    if (m_conflict_matrix != json.end() && m_conflict_matrix->size() != n_items)
    {
        LOG_ERROR("wrong size: conflict_matrix");
        return SCIP_READERROR;
    }

    // Create the empty SCIP problem
    SCIP_CALL(Probdata::create_prob(scip,
                                    filename,
                                    type,
                                    n_items,
                                    n_stacks,
                                    capacity,
                                    has_adjacent_conflicts));

    auto &probdata = Probdata::get(scip);

    // Set arrival times
    if (m_arrivals != json.end())
    {
        std::vector<std::size_t> arrivals(n_items);
        for (std::size_t i = 0; i < n_items; ++i)
        {
            arrivals[i] = (*m_arrivals)[i].get<std::size_t>();
        }
        probdata.set_arrivals(arrivals);
    }

    // Set weights
    if (m_weights != json.end())
    {
        std::vector<double> weights(n_items);
        for (std::size_t i = 0; i < n_items; ++i)
        {
            weights[i] = (*m_weights)[i].get<double>();
        }
        probdata.set_weights(weights);
    }
    // Set conflict graph
    else if (m_conflict_graph != json.end())
    {
        std::vector<std::vector<bool>> conflict_matrix(n_items,
                                                       std::vector<bool>(n_items,
                                                                         false));

        for (auto edge : *m_conflict_graph)
        {
            if (edge.size() != 2)
            {
                LOG_ERROR("malformed edge in conflict graph");
                return SCIP_READERROR;
            }
            std::size_t i = edge[0];
            std::size_t j = edge[1];
            if (i <= 0 || i > n_items || j <= 0 || j > n_items)
            {
                LOG_ERROR("invalid edge (%zu, %zu) in conflict graph", i, j);
                return SCIP_READERROR;
            }
            conflict_matrix[i - 1][j - 1] = true;
        }
        probdata.set_conflict_matrix(conflict_matrix);
    }
    // Set conflict matrix
    else if (m_conflict_matrix != json.end())
    {
        std::vector<std::vector<bool>> conflict_matrix(n_items,
                                                       std::vector<bool>(n_items));

        for (std::size_t i = 0; i < n_items; ++i)
        {
            for (std::size_t j = 0; j < n_items; ++j)
            {
                conflict_matrix[i][j] = (*m_conflict_matrix)[i][j].get<int>() != 0;
            }
        }
        probdata.set_conflict_matrix(conflict_matrix);
    }

    // Set departures
    if (m_departures != json.end())
    {
        std::vector<double> departures(n_items);
        for (std::size_t i = 0; i < n_items; ++i)
        {
            departures[i] = (*m_departures)[i].get<double>();
        }
        probdata.set_departures(departures);
    }
    // Set blocking graph
    else if (m_blocking_graph != json.end())
    {
        std::vector<std::vector<bool>> blocking_matrix(n_items,
                                                       std::vector<bool>(n_items,
                                                                         false));

        for (auto edge : *m_blocking_graph)
        {
            if (edge.size() != 2)
            {
                LOG_ERROR("malformed edge in blocking graph");
                return SCIP_READERROR;
            }
            std::size_t i = edge[0];
            std::size_t j = edge[1];
            if (i <= 0 || i > n_items || j <= 0 || j > n_items)
            {
                LOG_ERROR("invalid edge (%zu, %zu) in blocking graph", i, j);
                return SCIP_READERROR;
            }
            blocking_matrix[i - 1][j - 1] = true;
        }
        probdata.set_blocking_matrix(blocking_matrix);
    }
    // Set blocking matrix
    else if (m_blocking_matrix != json.end())
    {
        std::vector<std::vector<bool>> blocking_matrix(n_items,
                                                       std::vector<bool>(n_items));

        for (std::size_t i = 0; i < n_items; ++i)
        {
            for (std::size_t j = 0; j < n_items; ++j)
            {
                blocking_matrix[i][j] = (*m_blocking_matrix)[i][j].get<int>() != 0;
            }
        }
        probdata.set_blocking_matrix(blocking_matrix);
    }

    // Set position of initial items
    if (m_initial_positions != json.end())
    {
        std::size_t n_fix = m_initial_positions->size();
        std::vector<std::size_t> initial_positions(n_fix);
        for (std::size_t i = 0; i < n_fix; ++i)
        {
            initial_positions[i] = (*m_initial_positions)[i].get<std::size_t>() - 1;
        }
        probdata.set_initial_positions(initial_positions);
    }

    // Transform conflicts if the stacking matrix is transitive
    if (_transformconflicts &&
        probdata.type() == Probdata::Type::MinNStacks &&
        probdata.is_stacking_transitive())
    {
        SCIP_CALL(probdata.transform_conflicts());
        LOG_INFO("conflicts transformed");
    }

    SCIP_CALL(probdata.initialize(scip));

    *result = SCIP_SUCCESS;

    return SCIP_OKAY;
}

SCIP_RETCODE ReaderJson::scip_write(SCIP *scip,
                                    [[maybe_unused]] SCIP_READER *reader,
                                    FILE *file,
                                    [[maybe_unused]] const char *name,
                                    [[maybe_unused]] SCIP_PROBDATA *probdata,
                                    [[maybe_unused]] SCIP_Bool transformed,
                                    [[maybe_unused]] SCIP_OBJSENSE objsense,
                                    [[maybe_unused]] SCIP_Real objscale,
                                    [[maybe_unused]] SCIP_Real objoffset,
                                    [[maybe_unused]] SCIP_VAR **vars,
                                    [[maybe_unused]] int nvars,
                                    [[maybe_unused]] int nbinvars,
                                    [[maybe_unused]] int nintvars,
                                    [[maybe_unused]] int nimplvars,
                                    [[maybe_unused]] int ncontvars,
                                    [[maybe_unused]] SCIP_VAR **fixedvars,
                                    [[maybe_unused]] int nfixedvars,
                                    [[maybe_unused]] int startnvars,
                                    [[maybe_unused]] SCIP_CONS **conss,
                                    [[maybe_unused]] int nconss,
                                    [[maybe_unused]] int maxnconss,
                                    [[maybe_unused]] int startnconss,
                                    [[maybe_unused]] SCIP_Bool genericnames,
                                    SCIP_RESULT *result)
{
    nlohmann::json json;

    json["solver"] = "scip";

    json["prob_name"] = SCIPgetProbName(scip);

    SCIP_HEUR **heurs = SCIPgetHeurs(scip);
    int n_heurs = SCIPgetNHeurs(scip);

    json["heurs"] = nlohmann::json::array();

    for (int t = 0; t < n_heurs; ++t)
    {
        SCIP_HEUR *heur = heurs[t];

        json["heurs"].push_back({{"name", SCIPheurGetName(heur)},
                                 {"priority", SCIPheurGetPriority(heur)},
                                 {"dispchar", SCIPheurGetDispchar(heur)},
                                 {"timingmask", SCIPheurGetTimingmask(heur)},
                                 {"uses_subscip", SCIPheurUsesSubscip(heur)},
                                 {"priority", SCIPheurGetPriority(heur)},
                                 {"freq", SCIPheurGetFreq(heur)},
                                 {"freqofs", SCIPheurGetFreqofs(heur)},
                                 {"maxdepth", SCIPheurGetMaxdepth(heur)},
                                 {"n_calls", SCIPheurGetNCalls(heur)},
                                 {"n_sols_found", SCIPheurGetNSolsFound(heur)},
                                 {"n_best_sols_found", SCIPheurGetNBestSolsFound(heur)},
                                 {"setup_time", SCIPheurGetSetupTime(heur)},
                                 {"time", SCIPheurGetTime(heur)}});
    }

    const auto &probvars = Probdata::get(scip).vars();
    SCIP_SOL **sols = SCIPgetSols(scip);
    int n_sols = SCIPgetNSols(scip);

    json["sols"] = nlohmann::json::array();

    for (int t = 0; t < n_sols; ++t)
    {
        SCIP_SOL *sol = sols[t];

        assert(sol != nullptr);

        SCIP_Bool feasible;
        SCIP_CALL(SCIPcheckSolOrig(scip, sol, &feasible, false, false));

        if (!feasible)
        {
            LOG_WARN("solution in master problem infeasible");
        }

        double obj = SCIPgetSolOrigObj(scip, sol);
        nlohmann::json vars = nlohmann::json::array();

        for (std::size_t v = 0; v < std::size(probvars); ++v)
        {
            if (SCIPgetSolVal(scip, sol, probvars[v]) > 0.5)
            {
                vars.push_back(v + 1);
            }
        }

        SCIP_HEUR *heur = SCIPsolGetHeur(sol);
        int heur_id = -1;

        if (heur != nullptr)
        {
            heur_id = std::distance(heurs, std::find(heurs,
                                                     heurs + n_heurs,
                                                     heur));

            if (heur_id >= n_heurs)
            {
                heur_id = -1;
            }
        }

        json["sols"].push_back({{"obj", obj},
                                {"vars", vars},
                                {"heur", heur_id + 1},
                                {"time", SCIPgetSolTime(scip, sol)},
                                {"runnum", SCIPgetSolRunnum(scip, sol)},
                                {"nodenum", SCIPgetSolNodenum(scip, sol)},
                                {"depth", SCIPsolGetDepth(sol)},
                                {"type_n", SCIPsolGetType(sol)},
                                {"type", to_string(SCIPsolGetType(sol))}});
    }

    SCIP_PRICER **pricers = SCIPgetPricers(scip);
    int n_pricers = SCIPgetNPricers(scip);
    int n_active_pricers = SCIPgetNActivePricers(scip);

    json["pricers"] = nlohmann::json::array();

    for (int t = 0; t < n_active_pricers; ++t)
    {
        SCIP_PRICER *pricer = pricers[t];

        json["pricers"].push_back({{"name", SCIPpricerGetName(pricer)},
                                   {"priority", SCIPpricerGetPriority(pricer)},
                                   {"n_calls", SCIPpricerGetNCalls(pricer)},
                                   {"n_vars_found", SCIPpricerGetNVarsFound(pricer)},
                                   {"setup_time", SCIPpricerGetSetupTime(pricer)},
                                   {"time", SCIPpricerGetTime(pricer)}});
    }

    json["vars"] = nlohmann::json::array();

    for (auto var : probvars)
    {
        const auto vardata = Vardata::get(scip, var);
        const auto &s = vardata.items();
        SCIP_PRICER *pricer = vardata.pricer();
        int pricer_id = -1;
        double obj = SCIPvarGetObj(var);
        nlohmann::json items = nlohmann::json::array();

        for (auto i : s)
        {
            items.push_back(i + 1);
        }

        if (pricer != nullptr)
        {
            pricer_id = std::distance(pricers, std::find(pricers,
                                                         pricers + n_pricers,
                                                         pricer));

            if (pricer_id >= n_active_pricers)
            {
                LOG_WARN("var generated by inactive pricer?");
            }

            if (pricer_id >= n_pricers)
            {
                pricer_id = -1;
            }
        }

        json["vars"].push_back({{"obj", obj},
                                {"items", items},
                                {"name", SCIPvarGetName(var)},
                                {"pricer", pricer_id + 1}});
    }

    SCIP_BRANCHRULE **branchrules = SCIPgetBranchrules(scip);
    int n_branchrules = SCIPgetNBranchrules(scip);

    json["branchrules"] = nlohmann::json::array();

    for (int t = 0; t < n_branchrules; ++t)
    {
        SCIP_BRANCHRULE *branchrule = branchrules[t];

        json["branchrules"].push_back({{"name", SCIPbranchruleGetName(branchrule)},
                                       {"priority", SCIPbranchruleGetPriority(branchrule)},
                                       {"maxdepth", SCIPbranchruleGetMaxdepth(branchrule)},
                                       {"maxbounddist", SCIPbranchruleGetMaxbounddist(branchrule)},
                                       {"setup_time", SCIPbranchruleGetSetupTime(branchrule)},
                                       {"time", SCIPbranchruleGetTime(branchrule)},
                                       {"n_lp_calls", SCIPbranchruleGetNLPCalls(branchrule)},
                                       {"n_extern_calls", SCIPbranchruleGetNExternCalls(branchrule)},
                                       {"n_pseudo_calls", SCIPbranchruleGetNPseudoCalls(branchrule)},
                                       {"n_cutoffs", SCIPbranchruleGetNCutoffs(branchrule)},
                                       {"n_cuts_found", SCIPbranchruleGetNCutsFound(branchrule)},
                                       {"n_conss_found", SCIPbranchruleGetNConssFound(branchrule)},
                                       {"n_domreds_found", SCIPbranchruleGetNDomredsFound(branchrule)},
                                       {"n_children", SCIPbranchruleGetNChildren(branchrule)}});
    }

    SCIP_CONSHDLR **conshdlrs = SCIPgetConshdlrs(scip);
    int n_conshdlrs = SCIPgetNConshdlrs(scip);

    json["conshdlrs"] = nlohmann::json::array();

    for (int t = 0; t < n_conshdlrs; ++t)
    {
        SCIP_CONSHDLR *conshdlr = conshdlrs[t];
        auto needs_cons = SCIPconshdlrNeedsCons(conshdlr);
        auto max_n_active_conss = SCIPconshdlrGetMaxNActiveConss(conshdlr);

        if (needs_cons && max_n_active_conss == 0)
        {
            continue;
        }

        auto total_time = SCIPconshdlrGetSepaTime(conshdlr) +
                          SCIPconshdlrGetPropTime(conshdlr) +
                          SCIPconshdlrGetStrongBranchPropTime(conshdlr) +
                          SCIPconshdlrGetEnfoLPTime(conshdlr) +
                          SCIPconshdlrGetEnfoPSTime(conshdlr) +
                          SCIPconshdlrGetEnfoRelaxTime(conshdlr) +
                          SCIPconshdlrGetCheckTime(conshdlr) +
                          SCIPconshdlrGetRespropTime(conshdlr) +
                          SCIPconshdlrGetSetupTime(conshdlr);

        json["conshdlrs"].push_back({{"name", SCIPconshdlrGetName(conshdlr)},
                                     {"sepa_priority", SCIPconshdlrGetSepaPriority(conshdlr)},
                                     {"enfo_priority", SCIPconshdlrGetEnfoPriority(conshdlr)},
                                     {"check_priority", SCIPconshdlrGetCheckPriority(conshdlr)},
                                     {"sepa_freq", SCIPconshdlrGetSepaFreq(conshdlr)},
                                     {"prop_freq", SCIPconshdlrGetPropFreq(conshdlr)},
                                     {"eager_freq", SCIPconshdlrGetEagerFreq(conshdlr)},
                                     {"needs_cons", needs_cons},
                                     {"does_presolve", SCIPconshdlrDoesPresolve(conshdlr)},
                                     {"delay_sepa", SCIPconshdlrIsSeparationDelayed(conshdlr)},
                                     {"delay_prop", SCIPconshdlrIsPropagationDelayed(conshdlr)},
                                     {"prop_timing", SCIPconshdlrGetPropTiming(conshdlr)},
                                     {"presol_timing", SCIPconshdlrGetPresolTiming(conshdlr)},
                                     {"n_conss", SCIPconshdlrGetNConss(conshdlr)},
                                     {"total_time", total_time},
                                     {"setup_time", SCIPconshdlrGetSetupTime(conshdlr)},
                                     {"presol_time", SCIPconshdlrGetPresolTime(conshdlr)},
                                     {"sepa_time", SCIPconshdlrGetSepaTime(conshdlr)},
                                     {"enfo_lp_time", SCIPconshdlrGetEnfoLPTime(conshdlr)},
                                     {"enfo_ps_time", SCIPconshdlrGetEnfoPSTime(conshdlr)},
                                     {"enfo_relax_time", SCIPconshdlrGetEnfoRelaxTime(conshdlr)},
                                     {"prop_time", SCIPconshdlrGetPropTime(conshdlr)},
                                     {"strong_branch_prop_time", SCIPconshdlrGetStrongBranchPropTime(conshdlr)},
                                     {"check_time", SCIPconshdlrGetCheckTime(conshdlr)},
                                     {"resprop_time", SCIPconshdlrGetRespropTime(conshdlr)},
                                     {"n_sepa_calls", SCIPconshdlrGetNSepaCalls(conshdlr)},
                                     {"n_enfo_lp_calls", SCIPconshdlrGetNEnfoLPCalls(conshdlr)},
                                     {"n_enfo_ps_calls", SCIPconshdlrGetNEnfoPSCalls(conshdlr)},
                                     {"n_enfo_relax_calls", SCIPconshdlrGetNEnfoRelaxCalls(conshdlr)},
                                     {"n_prop_calls", SCIPconshdlrGetNPropCalls(conshdlr)},
                                     {"n_check_calls", SCIPconshdlrGetNCheckCalls(conshdlr)},
                                     {"n_resprop_calls", SCIPconshdlrGetNRespropCalls(conshdlr)},
                                     {"n_cutoffs", SCIPconshdlrGetNCutoffs(conshdlr)},
                                     {"n_cuts_found", SCIPconshdlrGetNCutsFound(conshdlr)},
                                     {"n_cuts_applied", SCIPconshdlrGetNCutsApplied(conshdlr)},
                                     {"n_conss_found", SCIPconshdlrGetNConssFound(conshdlr)},
                                     {"n_domreds_found", SCIPconshdlrGetNDomredsFound(conshdlr)},
                                     {"n_children", SCIPconshdlrGetNChildren(conshdlr)},
                                     {"max_n_active_conss", max_n_active_conss},
                                     {"start_n_active_conss", SCIPconshdlrGetStartNActiveConss(conshdlr)},
                                     {"n_fixed_vars", SCIPconshdlrGetNFixedVars(conshdlr)},
                                     {"n_aggr_vars", SCIPconshdlrGetNAggrVars(conshdlr)},
                                     {"n_chg_var_types", SCIPconshdlrGetNChgVarTypes(conshdlr)},
                                     {"n_chg_bds", SCIPconshdlrGetNChgBds(conshdlr)},
                                     {"n_add_holes", SCIPconshdlrGetNAddHoles(conshdlr)},
                                     {"n_del_conss", SCIPconshdlrGetNDelConss(conshdlr)},
                                     {"n_add_conss", SCIPconshdlrGetNAddConss(conshdlr)},
                                     {"n_upgd_conss", SCIPconshdlrGetNUpgdConss(conshdlr)},
                                     {"n_chg_coefs", SCIPconshdlrGetNChgCoefs(conshdlr)},
                                     {"n_chg_sides", SCIPconshdlrGetNChgSides(conshdlr)},
                                     {"n_presol_calls", SCIPconshdlrGetNPresolCalls(conshdlr)}});
    }

    json["stats"] = {
        {"status_n", SCIPgetStatus(scip)},
        {"status", to_string(SCIPgetStatus(scip))},
        {"pressed_ctrl_c", SCIPpressedCtrlC(scip)},
        {"n_vars", SCIPgetNVars(scip)},
        {"n_bin_vars", SCIPgetNBinVars(scip)},
        {"n_int_vars", SCIPgetNIntVars(scip)},
        {"n_impl_vars", SCIPgetNImplVars(scip)},
        {"n_cont_vars", SCIPgetNContVars(scip)},
        {"n_obj_vars", SCIPgetNObjVars(scip)},
        {"n_fixed_vars", SCIPgetNFixedVars(scip)},
        {"n_orig_vars", SCIPgetNOrigVars(scip)},
        {"n_orig_bin_vars", SCIPgetNOrigBinVars(scip)},
        {"n_orig_int_vars", SCIPgetNOrigIntVars(scip)},
        {"n_orig_impl_vars", SCIPgetNOrigImplVars(scip)},
        {"n_orig_cont_vars", SCIPgetNOrigContVars(scip)},
        {"n_total_vars", SCIPgetNTotalVars(scip)},
        {"n_upgr_conss", SCIPgetNUpgrConss(scip)},
        {"n_conss", SCIPgetNConss(scip)},
        {"n_orig_conss", SCIPgetNOrigConss(scip)},
        {"n_check_conss", SCIPgetNCheckConss(scip)},
        {"total_time", SCIPgetTotalTime(scip)},
        {"solving_time", SCIPgetSolvingTime(scip)},
        {"reading_time", SCIPgetReadingTime(scip)},
        {"presolving_time", SCIPgetPresolvingTime(scip)},
        {"first_lp_time", SCIPgetFirstLPTime(scip)},
        {"n_runs", SCIPgetNRuns(scip)},
        {"n_reopt_runs", SCIPgetNReoptRuns(scip)},
        {"n_nodes", SCIPgetNNodes(scip)},
        {"n_total_nodes", SCIPgetNTotalNodes(scip)},
        {"n_feasible_leaves", SCIPgetNFeasibleLeaves(scip)},
        {"n_infeasible_leaves", SCIPgetNInfeasibleLeaves(scip)},
        {"n_objlim_leaves", SCIPgetNObjlimLeaves(scip)},
        {"n_delayed_cutoffs", SCIPgetNDelayedCutoffs(scip)},
        {"n_lps", SCIPgetNLPs(scip)},
        {"n_lp_iterations", SCIPgetNLPIterations(scip)},
        {"n_nzs", SCIPgetNNZs(scip)},
        {"n_root_lp_iterations", SCIPgetNRootLPIterations(scip)},
        {"n_root_first_lp_iterations", SCIPgetNRootFirstLPIterations(scip)},
        {"n_primal_lps", SCIPgetNPrimalLPs(scip)},
        {"n_primal_lp_iterations", SCIPgetNPrimalLPIterations(scip)},
        {"n_dual_lps", SCIPgetNDualLPs(scip)},
        {"n_dual_lp_iterations", SCIPgetNDualLPIterations(scip)},
        {"n_barrier_lps", SCIPgetNBarrierLPs(scip)},
        {"n_barrier_lp_iterations", SCIPgetNBarrierLPIterations(scip)},
        {"n_resolve_lps", SCIPgetNResolveLPs(scip)},
        {"n_resolve_lp_iterations", SCIPgetNResolveLPIterations(scip)},
        {"n_primal_resolve_lps", SCIPgetNPrimalResolveLPs(scip)},
        {"n_primal_resolve_lp_iterations", SCIPgetNPrimalResolveLPIterations(scip)},
        {"n_dual_resolve_lps", SCIPgetNDualResolveLPs(scip)},
        {"n_dual_resolve_lp_iterations", SCIPgetNDualResolveLPIterations(scip)},
        {"n_node_lps", SCIPgetNNodeLPs(scip)},
        {"n_node_lp_iterations", SCIPgetNNodeLPIterations(scip)},
        {"n_node_init_lps", SCIPgetNNodeInitLPs(scip)},
        {"n_node_init_lp_iterations", SCIPgetNNodeInitLPIterations(scip)},
        {"n_diving_lps", SCIPgetNDivingLPs(scip)},
        {"n_diving_lp_iterations", SCIPgetNDivingLPIterations(scip)},
        {"n_strongbranchs", SCIPgetNStrongbranchs(scip)},
        {"n_strongbranch_lp_iterations", SCIPgetNStrongbranchLPIterations(scip)},
        {"n_root_strongbranchs", SCIPgetNRootStrongbranchs(scip)},
        {"n_root_strongbranch_lp_iterations", SCIPgetNRootStrongbranchLPIterations(scip)},
        {"n_pricevars", SCIPgetNPricevars(scip)},
        {"n_pricevars_found", SCIPgetNPricevarsFound(scip)},
        {"n_pricevars_applied", SCIPgetNPricevarsApplied(scip)},
        {"n_cuts_found", SCIPgetNCutsFound(scip)},
        {"n_cuts_applied", SCIPgetNCutsApplied(scip)},
        {"n_conflict_conss_found", SCIPgetNConflictConssFound(scip)},
        {"n_conflict_conss_applied", SCIPgetNConflictConssApplied(scip)},
        {"max_depth", SCIPgetMaxDepth(scip)},
        {"max_total_depth", SCIPgetMaxTotalDepth(scip)},
        {"n_backtracks", SCIPgetNBacktracks(scip)},
        {"avg_dualbound", SCIPgetAvgDualbound(scip)},
        {"avg_lowerbound", SCIPgetAvgLowerbound(scip)},
        {"dualbound", SCIPgetDualbound(scip)},
        {"lowerbound", SCIPgetLowerbound(scip)},
        {"dualbound_root", SCIPgetDualboundRoot(scip)},
        {"lowerbound_root", SCIPgetLowerboundRoot(scip)},
        {"first_lp_dualbound_root", SCIPgetFirstLPDualboundRoot(scip)},
        {"first_lp_lowerbound_root", SCIPgetFirstLPLowerboundRoot(scip)},
        {"first_primalbound", SCIPgetFirstPrimalBound(scip)},
        {"primalbound", SCIPgetPrimalbound(scip)},
        {"upperbound", SCIPgetUpperbound(scip)},
        {"cutoffbound", SCIPgetCutoffbound(scip)},
        {"gap", SCIPgetGap(scip)},
        {"trans_gap", SCIPgetTransGap(scip)},
        {"n_sols_found", SCIPgetNSolsFound(scip)},
        {"n_lim_sols_found", SCIPgetNLimSolsFound(scip)},
        {"n_best_sols_found", SCIPgetNBestSolsFound(scip)},
        {"deterministic_time", SCIPgetDeterministicTime(scip)},
        {"mem_used", SCIPgetMemUsed(scip)},
        {"mem_total", SCIPgetMemTotal(scip)},
        {"mem_extern_estim", SCIPgetMemExternEstim(scip)}};

    SCIP_SOL *best_sol = SCIPgetBestSol(scip);
    if (best_sol != nullptr)
    {
        json["stats"]["best_sol_orig_obj"] = SCIPgetSolOrigObj(scip, best_sol);
    }

    std::ostringstream oss;

    oss << json;

    SCIPinfoMessage(scip, file, oss.str().c_str());

    *result = SCIP_SUCCESS;

    return SCIP_OKAY;
}