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

#include "branchrule_samediff.hpp"
#include "conshdlr_samediff.hpp"
#include "heur_init.hpp"
#include "pricer.hpp"
#include "reader_json.hpp"
#include <cstdlib>
#include <objscip/objscip.h>
#include <scip/branch_allfullstrong.h>
#include <scip/branch_fullstrong.h>
#include <scip/branch_inference.h>
#include <scip/branch_leastinf.h>
#include <scip/branch_mostinf.h>
#include <scip/branch_pscost.h>
#include <scip/branch_random.h>
#include <scip/branch_relpscost.h>
#include <scip/cons_integral.h>
#include <scip/cons_knapsack.h>
#include <scip/cons_setppc.h>
#include <scip/dialog_default.h>
#include <scip/disp_default.h>
#include <scip/nodesel_bfs.h>
#include <scip/nodesel_dfs.h>
#include <scip/nodesel_estimate.h>
#include <scip/nodesel_hybridestim.h>
#include <scip/nodesel_restartdfs.h>
#include <scip/reader_lp.h>
#include <scip/scip.h>
#include <scip/scipshell.h>
#include <scip/table_default.h>
#include <sstream>
#include <tuple>

using pslp::bp::BranchruleSamediff;
using pslp::bp::ConshdlrSamediff;
using pslp::bp::HeurInit;
using pslp::bp::Pricer;
using pslp::bp::ReaderJson;

// Pricers
#include "pricer_cplex_lazy.hpp"
#include "pricer_cplex_lo.hpp"
#include "pricer_cplex_nf.hpp"
#include "pricer_cplex_pa.hpp"
#include "pricer_greedy.hpp"
#include "pricer_mip_lazy.hpp"
#include "pricer_multi.hpp"

static void register_pricers()
{
    Pricer::register_pricer<pslp::bp::PricerGreedy>(true);
    Pricer::register_pricer<pslp::bp::PricerCplexLo>(false);
    Pricer::register_pricer<pslp::bp::PricerCplexNf>(false);
    Pricer::register_pricer<pslp::bp::PricerCplexPa>(true);
    Pricer::register_pricer<pslp::bp::PricerCplexLazy>(false);
    Pricer::register_pricer<pslp::bp::PricerMipLazy>(false);

    using PricerMulti = pslp::bp::PricerMulti<pslp::bp::PricerGreedy,
                                              pslp::bp::PricerCplexLo,
                                              pslp::bp::PricerCplexNf,
                                              pslp::bp::PricerCplexPa,
                                              pslp::bp::PricerCplexLazy,
                                              pslp::bp::PricerMipLazy>;

    Pricer::register_pricer<PricerMulti>(false);
}

[[nodiscard]] static SCIP_RETCODE include_plugins(SCIP *scip)
{
    SCIP_CALL(SCIPincludeBranchruleAllfullstrong(scip));
    SCIP_CALL(SCIPincludeBranchruleFullstrong(scip));
    SCIP_CALL(SCIPincludeBranchruleInference(scip));
    SCIP_CALL(SCIPincludeBranchruleMostinf(scip));
    SCIP_CALL(SCIPincludeBranchruleLeastinf(scip));
    SCIP_CALL(SCIPincludeBranchrulePscost(scip));
    SCIP_CALL(SCIPincludeBranchruleRandom(scip));
    SCIP_CALL(SCIPincludeBranchruleRelpscost(scip));

    SCIP_CALL(SCIPincludeConshdlrIntegral(scip));
    SCIP_CALL(SCIPincludeConshdlrSetppc(scip));
    SCIP_CALL(SCIPincludeConshdlrKnapsack(scip));

    SCIP_CALL(SCIPincludeNodeselBfs(scip));
    SCIP_CALL(SCIPincludeNodeselDfs(scip));
    SCIP_CALL(SCIPincludeNodeselEstimate(scip));
    SCIP_CALL(SCIPincludeNodeselHybridestim(scip));
    SCIP_CALL(SCIPincludeNodeselRestartdfs(scip));

    SCIP_CALL(SCIPincludeReaderLp(scip));

    SCIP_CALL(SCIPincludeTableDefault(scip));

    SCIP_CALL(SCIPincludeObjBranchrule(scip, new BranchruleSamediff(scip), true));
    SCIP_CALL(SCIPincludeObjConshdlr(scip, new ConshdlrSamediff(scip), true));
    SCIP_CALL(SCIPincludeObjHeur(scip, new HeurInit(scip), true));
    SCIP_CALL(SCIPincludeObjReader(scip, new ReaderJson(scip), true));

    SCIP_CALL(SCIPincludeDispDefault(scip));
    SCIP_CALL(SCIPincludeDialogDefault(scip));

    const auto &pricers = Pricer::list();

    for (const auto &info : pricers)
    {
        SCIP_CALL(SCIPincludeObjPricer(scip, info.create(scip), true));
    }

    return SCIP_OKAY;
}

[[nodiscard]] static SCIP_RETCODE add_params(SCIP *scip)
{
    const auto &pricers = Pricer::list();

    for (const auto &info : pricers)
    {
        std::ostringstream oss_name;
        std::ostringstream oss_desc;

        oss_name << "pricers/" << info.name << "/activate";
        oss_desc << "(pslp) should pricer <" << info.name << "> be activated?";

        SCIP_CALL(SCIPaddBoolParam(scip,
                                   oss_name.str().c_str(),
                                   oss_desc.str().c_str(),
                                   nullptr,
                                   false,
                                   info.activate,
                                   nullptr,
                                   nullptr));
    }

    return SCIP_OKAY;
}

[[nodiscard]] static SCIP_RETCODE run(int argc, char **argv)
{
    SCIP *scip;

    register_pricers();

    SCIP_CALL(SCIPcreate(&scip));

#ifndef NDEBUG
    SCIPenableDebugSol(scip);
#endif

    SCIP_CALL(include_plugins(scip));
    SCIP_CALL(add_params(scip));

    SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrestarts", 0));
    SCIP_CALL(SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, true));

    SCIP_CALL(SCIPprocessShellArguments(scip, argc, argv, nullptr));

    SCIP_CALL(SCIPfree(&scip));

#ifndef NDEBUG
    BMScheckEmptyMemory();
#endif

    return SCIP_OKAY;
}

int main(int argc, char **argv)
{
    SCIP_RETCODE retcode = run(argc, argv);

    if (retcode != SCIP_OKAY)
    {
        SCIPprintError(retcode);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}