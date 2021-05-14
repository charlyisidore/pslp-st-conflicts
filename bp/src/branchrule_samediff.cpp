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
#include "logging.hpp"
#include "probdata.hpp"
#include "symmetric_matrix.hpp"
#include "vardata.hpp"
#include <algorithm>
#include <cassert>
#include <sstream>
#include <vector>

using namespace pslp::bp;

const char *BranchruleSamediff::Properties::Name = "samediff";
const char *BranchruleSamediff::Properties::Desc = "same/diff branching rule for the PSLP";
const int BranchruleSamediff::Properties::Priority = 50000;
const int BranchruleSamediff::Properties::Maxdepth = -1;
const SCIP_Real BranchruleSamediff::Properties::Maxbounddist = 1.;

BranchruleSamediff::BranchruleSamediff(SCIP *scip)
    : scip::ObjBranchrule(scip,
                          Properties::Name,
                          Properties::Desc,
                          Properties::Priority,
                          Properties::Maxdepth,
                          Properties::Maxbounddist)
{
}

SCIP_RETCODE BranchruleSamediff::scip_execlp(SCIP *scip,
                                             [[maybe_unused]] SCIP_BRANCHRULE *branchrule,
                                             [[maybe_unused]] SCIP_Bool allowaddcons,
                                             SCIP_RESULT *result)
{
    assert(result != nullptr);

    const auto &probdata = Probdata::get(scip);
    auto n = probdata.n_items();

    SCIP_VAR **lpcands;
    SCIP_Real *lpcandsfrac;
    int nlpcands;

    SCIP_CALL(SCIPgetLPBranchCands(scip,
                                   &lpcands,
                                   nullptr,
                                   &lpcandsfrac,
                                   nullptr,
                                   &nlpcands,
                                   nullptr));

    assert(lpcands != nullptr);
    assert(lpcandsfrac != nullptr);
    assert(nlpcands > 0);

    // Compute weights for each pair of items
    SymmetricMatrix<double, true> weights(n, 0.);

    for (int v = 0; v < nlpcands; ++v)
    {
        SCIP_VAR *var = lpcands[v];
        SCIP_Real frac = lpcandsfrac[v];

        assert(var != nullptr);

        const auto &vardata = Vardata::get(scip, var);
        const auto &s = vardata.items();

        for (auto i = std::begin(s); i != std::end(s); ++i)
        {
            // w_{ii} = sum_s a_i^s x_s
            weights(*i, *i) += frac;

            for (auto j = std::begin(s); j != i; ++j)
            {
                // w_{ij} = sum_s a_i^s a_j^s x_s
                weights(*i, *j) += frac;

                assert(SCIPisGE(scip, weights(*i, *j), 0.));
            }
        }
    }

    // Select the pair with the weight closest to 0.5
    double best_w = 0.;
    std::size_t best_i = n;
    std::size_t best_j = n;

    for (std::size_t i = 0; i < n; ++i)
    {
        const auto w_ii = weights(i, i);

        for (std::size_t j = 0; j < i; ++j)
        {
            const auto w_ij = weights(i, j);
            const auto w_jj = weights(j, j);
            const auto w = std::min(w_ij, 1. - w_ij);

            if (best_w < w)
            {
                // w_{ii} = w_{ij} means that no variable contains only one of i and j
                // We need a variable that contain i or j but not both
                if (SCIPisEQ(scip, w_ij, w_ii) &&
                    SCIPisEQ(scip, w_ij, w_jj))
                {
                    continue;
                }

                best_w = w;
                best_i = i;
                best_j = j;
            }
        }
    }

    assert(best_w > 0.);
    assert(best_i < n);
    assert(best_j < n);

    // When conflicts are distant, if i and j have conflict in both directions,
    // they cannot be put together at all. When conflicts are adjacent, it is
    // possible to put i and j together if there is at least an intermediate item
    assert(probdata.has_adjacent_conflicts() ||
           !probdata.conflict_matrix()(best_i, best_j) ||
           !probdata.conflict_matrix()(best_j, best_i));

    LOG_DEBUG("branch on <%zu,%zu> at node %lld, depth %d",
              best_i + 1,
              best_j + 1,
              SCIPgetNNodes(scip),
              SCIPgetDepth(scip));

    SCIP_NODE *child_same;
    SCIP_NODE *child_diff;
    SCIP_CONS *cons_same;
    SCIP_CONS *cons_diff;

    SCIP_CALL(SCIPcreateChild(scip,
                              &child_same,
                              0.,
                              SCIPgetLocalTransEstimate(scip)));

    SCIP_CALL(SCIPcreateChild(scip,
                              &child_diff,
                              0.,
                              SCIPgetLocalTransEstimate(scip)));

    {
        std::ostringstream oss;
        oss << "same{" << best_i + 1 << ',' << best_j + 1 << '}';
        SCIP_CALL(ConshdlrSamediff::create_cons(scip,
                                                &cons_same,
                                                oss.str().c_str(),
                                                ConshdlrSamediff::Type::Same,
                                                best_i,
                                                best_j,
                                                child_same));
    }

    {
        std::ostringstream oss;
        oss << "diff{" << best_i + 1 << ',' << best_j + 1 << '}';
        SCIP_CALL(ConshdlrSamediff::create_cons(scip,
                                                &cons_diff,
                                                oss.str().c_str(),
                                                ConshdlrSamediff::Type::Diff,
                                                best_i,
                                                best_j,
                                                child_diff));
    }

    SCIP_CALL(SCIPaddConsNode(scip, child_same, cons_same, nullptr));
    SCIP_CALL(SCIPaddConsNode(scip, child_diff, cons_diff, nullptr));

    SCIP_CALL(SCIPreleaseCons(scip, &cons_same));
    SCIP_CALL(SCIPreleaseCons(scip, &cons_diff));

    *result = SCIP_BRANCHED;

    return SCIP_OKAY;
}