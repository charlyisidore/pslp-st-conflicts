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

#include "heur_init.hpp"
#include "layout.hpp"
#include "logging.hpp"
#include "probdata.hpp"
#include <algorithm>
#include <cassert>
#include <numeric>
#include <vector>

using namespace pslp::bp;

// Print the layout
constexpr bool PrintLayout = false;

const char *HeurInit::Properties::Name = "pslp_init";
const char *HeurInit::Properties::Desc = "initial primal heuristic for the PSLP";
const char HeurInit::Properties::Dispchar = SCIP_HEURDISPCHAR_TRIVIAL;
const int HeurInit::Properties::Priority = 1;
const int HeurInit::Properties::Freq = 1;
const int HeurInit::Properties::Freqofs = 0;
const int HeurInit::Properties::Maxdepth = 0;
const SCIP_HEURTIMING HeurInit::Properties::Timing = SCIP_HEURTIMING_BEFORENODE;
const SCIP_Bool HeurInit::Properties::Usesubscip = false;

HeurInit::HeurInit(SCIP *scip)
    : scip::ObjHeur(scip,
                    Properties::Name,
                    Properties::Desc,
                    Properties::Dispchar,
                    Properties::Priority,
                    Properties::Freq,
                    Properties::Freqofs,
                    Properties::Maxdepth,
                    Properties::Timing,
                    Properties::Usesubscip)
{
}

SCIP_RETCODE HeurInit::scip_exec(SCIP *scip,
                                 SCIP_HEUR *heur,
                                 [[maybe_unused]] SCIP_HEURTIMING heurtiming,
                                 [[maybe_unused]] SCIP_Bool nodeinfeasible,
                                 SCIP_RESULT *result)
{
    assert(result != nullptr);

    auto &probdata = Probdata::get(scip);
    const auto n = probdata.n_items();
    const auto b = probdata.capacity();
    const auto n_initial = probdata.n_initial_items();
    const auto n_incoming = n - n_initial;
    const auto &conflict_graph = probdata.conflict_graph();
    const auto &conflict_matrix = probdata.conflict_matrix();
    std::vector<std::size_t> incoming_items(n_incoming);
    std::vector<std::vector<std::size_t>> layout;

    // Place the initial items
    for (std::size_t i = 0; i < n_initial; ++i)
    {
        const auto k = probdata.initial_position(i);
        if (k >= std::size(layout))
        {
            layout.resize(k + 1);
        }
        layout[k].push_back(i);
    }

    // Sort incoming items
    std::iota(std::begin(incoming_items), std::end(incoming_items), n_initial);
    std::sort(std::begin(incoming_items),
              std::end(incoming_items),
              [&conflict_graph](auto i, auto j) {
                  return conflict_graph.degree(i) > conflict_graph.degree(j);
              });

    // Insert each incoming item in the first available stack
    for (auto i : incoming_items)
    {
        bool placed = false;

        // Try stacks one by one until the item is placed
        for (std::size_t k = 0; !placed; ++k)
        {
            // Add empty stacks if necessary
            if (k >= std::size(layout))
            {
                layout.resize(k + 1);
            }

            const auto h = std::size(layout[k]);

            // Skip full stacks
            if (h >= b)
            {
                assert(h == b);
                continue;
            }

            // Check for conflicts
            if (probdata.has_adjacent_conflicts())
            {
                // Find the lowest level not causing conflicts
                for (std::size_t l = 0; l <= h; ++l)
                {
                    if (l > 0)
                    {
                        // Item just below
                        const auto j = layout[k][l - 1];
                        if (conflict_matrix(i, j))
                        {
                            continue;
                        }
                    }

                    if (l < h)
                    {
                        // Item just above
                        const auto j = layout[k][l];
                        if (conflict_matrix(j, i))
                        {
                            continue;
                        }
                    }

                    layout[k].insert(std::begin(layout[k]) + l, i);
                    placed = true;
                    break;
                }
            }
            else
            {
                // Find the min and max levels s.t. no conflict happens (O(h))
                std::size_t l_min = 0;
                std::size_t l_max = h;

                for (std::size_t l = 0; l < h && (l_min == 0 || l_max == h); ++l)
                {
                    const auto j_min = layout[k][h - l - 1];
                    const auto j_max = layout[k][l];

                    // Can item i be put below j_min?
                    if (l_min == 0 && conflict_matrix(j_min, i))
                    {
                        // Item i must be stricly above j_min (i.e. > h-l-1)
                        l_min = h - l;
                    }

                    // Can item i be put above j_max?
                    if (l_max == h && conflict_matrix(i, j_max))
                    {
                        // Item i must be strictly below j_max (i.e. <= l)
                        l_max = l;
                    }
                }

#ifndef NDEBUG
                {
                    // This O(h^2) algorithm should return the same result
                    std::size_t l_i = h + 1;
                    for (std::size_t l = 0; l <= h; ++l)
                    {
                        // Check all items below and above
                        bool is_violating = false;
                        for (std::size_t l_j = 0; l_j < h; ++l_j)
                        {
                            const auto j = layout[k][l_j];
                            if ((l_j < l && conflict_matrix(i, j)) ||
                                (l_j >= l && conflict_matrix(j, i)))
                            {
                                is_violating = true;
                                break;
                            }
                        }
                        if (is_violating)
                        {
                            continue;
                        }
                        l_i = l;
                        break;
                    }

                    if (l_i < h + 1)
                    {
                        assert(l_min <= l_max && l_min == l_i);
                    }
                    else
                    {
                        assert(l_min > l_max);
                    }
                }
#endif

                // No feasible insertion level found?
                if (l_min > l_max)
                {
                    continue;
                }

                layout[k].insert(std::begin(layout[k]) + l_min, i);
                placed = true;
            }
        }
    }

#ifndef NDEBUG
    // Check solution
    for (std::size_t k = 0; k < std::size(layout); ++k)
    {
        const auto h = std::size(layout[k]);

        if (probdata.has_adjacent_conflicts())
        {
            for (std::size_t l = 0; l < h - 1; ++l)
            {
                const auto i = layout[k][l + 1];
                const auto j = layout[k][l];

                if (conflict_matrix(i, j))
                {
                    LOG_ERROR("infeasible initial heuristic solution: %zu on top of %zu in stack %zu",
                              i + 1,
                              j + 1,
                              k + 1);
                }
            }
        }
    }
#endif

    if constexpr (PrintLayout)
    {
        const auto m = std::size(layout);
        Layout lyt(n, m, b);

        for (std::size_t k = 0; k < m; ++k)
        {
            for (auto i : layout[k])
            {
                lyt.push(k, i);
            }
        }

        lyt.print();
    }

    LOG_INFO("initial heuristic solution obj: %zu", std::size(layout));

    // Create a new SCIP solution
    SCIP_SOL *heur_sol;
    SCIP_CALL(SCIPcreateSol(scip, &heur_sol, heur));

    for (const auto &s : layout)
    {
        SCIP_CALL(probdata.add_var(scip, s, nullptr));

        SCIP_VAR *var = probdata.vars().back();
        SCIP_CALL(SCIPsetSolVal(scip, heur_sol, var, 1.));
    }

    SCIP_Bool stored;
    SCIP_CALL(SCIPaddSolFree(scip, &heur_sol, &stored));

    assert(stored);

    *result = SCIP_FOUNDSOL;

    return SCIP_OKAY;
}
