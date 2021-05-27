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

#include "pricer_greedy.hpp"
#include "logging.hpp"
#include "pprint.hpp"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <numeric>

// Use the "push front" strategy?
constexpr bool UsePushFront = false;

// Print the solution
constexpr bool PrintSolution = false;

using namespace pslp::bp;

const char *PricerGreedy::Properties::Name = "greedy";
const char *PricerGreedy::Properties::Desc = "greedy variable pricer";
const int PricerGreedy::Properties::Priority = 1000000;
const SCIP_Bool PricerGreedy::Properties::Delay = true;

PricerGreedy::PricerGreedy(SCIP *scip)
    : Pricer(scip,
             Properties::Name,
             Properties::Desc,
             Properties::Priority,
             Properties::Delay)
{
}

SCIP_RETCODE PricerGreedy::solve(SCIP *scip,
                                 const Problem &prob,
                                 std::vector<Solution> &sols)
{
#ifndef NDEBUG
    using Clock = std::chrono::steady_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    using Duration = std::chrono::duration<double>;

    TimePoint start_time;
    TimePoint finish_time;

    start_time = Clock::now();
#endif

    const auto n_items = prob.conflict_graph.n_vertices();

    // Recompute unions when same/diff constraints have changed
    if (prob.has_changed)
    {
        _unions.resize(n_items);
        _union_of.resize(n_items);

        for (std::size_t i = 0; i < n_items; ++i)
        {
            _unions[i] = {i};
            _union_of[i] = i;
        }

        for (auto [i, j] : prob.same)
        {
            const auto u = _union_of[i];
            const auto v = _union_of[j];

            if (u == v)
            {
                continue;
            }

            for (auto k : _unions[v])
            {
                // Sort items by out-degree in each union
                auto it = std::upper_bound(std::begin(_unions[u]),
                                           std::end(_unions[u]),
                                           k,
                                           [&prob](auto i, auto j) {
                                               return prob.conflict_graph.degree(i) >
                                                      prob.conflict_graph.degree(j);
                                           });

                _unions[u].insert(it, k);
                _union_of[k] = u;
            }

            _unions.erase(std::begin(_unions) + v);

            for (std::size_t k = 0; k < n_items; ++k)
            {
                assert(_union_of[k] != v);

                if (_union_of[k] > v)
                {
                    --_union_of[k];
                }
            }
        }

        // Express diff constraints as a matrix for quick access

        _diff_matrix.initialize(std::size(_unions), false);

        for (auto [i, j] : prob.diff)
        {
            const auto u = _union_of[i];
            const auto v = _union_of[j];

            assert(u != v);

            _diff_matrix(u, v) = true;
        }
    }

    assert(!std::empty(_unions) && !std::empty(_union_of));

    // Update coefs for each union
    const auto n_unions = std::size(_unions);
    std::vector<double> coefs(n_unions, 0.);

    for (std::size_t u = 0; u < n_unions; ++u)
    {
        assert(!std::empty(_unions[u]));

        for (auto i : _unions[u])
        {
            coefs[u] += prob.coefs[i];
        }
    }

    // Sort unions by coef
    std::vector<std::size_t> order(n_unions);
    std::iota(std::begin(order), std::end(order), 0);
    std::sort(std::begin(order),
              std::end(order),
              [this, &coefs](auto u, auto v) {
                  return coefs[u] / double(std::size(_unions[u])) >
                         coefs[v] / double(std::size(_unions[v]));
              });

    // Remove non-positive unions
    {
        auto it = std::find_if_not(std::begin(order),
                                   std::end(order),
                                   [&](auto u) {
                                       return SCIPisPositive(scip, coefs[u]);
                                   });
        order.erase(it, std::end(order));
    }

    // We store the selected unions for checking diff constraints
    std::vector<std::size_t> selected_unions;

    // We remember which union has been pushed front
    std::vector<bool> is_pushed_front;

    Solution sol;
    double obj = 0.;

    if constexpr (UsePushFront)
    {
        is_pushed_front.resize(std::size(order), false);
    }

    while (true)
    {
        // Try to select unions one by one
        for (auto u : order)
        {
            // We removed non-positve unions
            assert(SCIPisPositive(scip, coefs[u]));

            // Capacity is satisfied?
            if (std::size(sol) + std::size(_unions[u]) > prob.capacity)
            {
                continue;
            }

            bool is_infeasible = false;

            // The union violates diff constraints?
            for (auto v : selected_unions)
            {
                if (_diff_matrix(u, v))
                {
                    is_infeasible = true;
                    break;
                }
            }

            if (is_infeasible)
            {
                continue;
            }

            // Try to insert items of the union in the column
            auto s = sol;

            for (auto i : _unions[u])
            {
                const auto h = std::size(s);

                // Insert the item in the lowest available level
                if (prob.has_adjacent_conflicts)
                {
                    bool is_inserted = false;

                    for (std::size_t l = 0; l <= h; ++l)
                    {
                        // Check only adjacent items
                        if (l > 0)
                        {
                            const auto j = s[l - 1];
                            if (prob.conflict_matrix(i, j))
                            {
                                continue;
                            }
                        }

                        if (l < h)
                        {
                            const auto j = s[l];
                            if (prob.conflict_matrix(j, i))
                            {
                                continue;
                            }
                        }

                        s.insert(std::begin(s) + l, i);
                        is_inserted = true;
                        break;
                    }

                    // No feasible insertion level found
                    if (!is_inserted)
                    {
                        is_infeasible = true;
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
                        const auto j_min = s[h - l - 1];
                        const auto j_max = s[l];

                        // Can item i be put below j_min?
                        if (l_min == 0 && prob.conflict_matrix(j_min, i))
                        {
                            // Item i must be stricly above j_min (i.e. > h-l-1)
                            l_min = h - l;
                        }

                        // Can item i be put above j_max?
                        if (l_max == h && prob.conflict_matrix(i, j_max))
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
                                const auto j = s[l_j];
                                if ((l_j < l && prob.conflict_matrix(i, j)) ||
                                    (l_j >= l && prob.conflict_matrix(j, i)))
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
                        is_infeasible = true;
                        break;
                    }
                    else
                    {
                        // Choose the minimum level
                        s.insert(std::begin(s) + l_min, i);
                    }
                }
            }

            if (is_infeasible)
            {
                continue;
            }

            // Add the union and its items to the solution
            selected_unions.push_back(u);
            sol = s;

            assert(std::size(sol) <= prob.capacity);
        }

        obj = Pricer::obj_value(prob, sol);

        if constexpr (UsePushFront)
        {
            // We search for a column having a strictly positive objective value
            if (SCIPisPositive(scip, obj))
            {
                break;
            }
            else
            {
                // Search for an union that is not selected and has not been moved
                // to the front yet
                auto it = std::find_if(std::begin(order),
                                       std::end(order),
                                       [&](auto u) {
                                           return !is_pushed_front[u] &&
                                                  std::find(std::begin(selected_unions),
                                                            std::end(selected_unions),
                                                            u) == std::end(selected_unions);
                                       });

                // All non-selected unions have been pushed front?
                if (it == std::end(order))
                {
                    break;
                }
                else
                {
                    // Move the union to the front
                    auto u = *it;
                    order.erase(it);
                    order.insert(std::begin(order), u);
                    is_pushed_front[u] = true;
                    selected_unions.clear();
                    sol.clear();
                }
            }
        }
        else
        {
            break;
        }
    }

#ifndef NDEBUG
    finish_time = Clock::now();

    double time = Duration(finish_time - start_time).count();

    LOG_DEBUG("greedy pricing solving time: %g s", time);
#endif

    if (SCIPisPositive(scip, obj))
    {
        if constexpr (PrintSolution)
        {
            LOG_TRACE("-> sol: %s",
                      pprint_n(std::size(sol),
                               [&sol](auto i) {
                                   return sol[i] + 1;
                               })
                          .c_str());
        }

        sols.push_back(sol);
    }

    return SCIP_OKAY;
}