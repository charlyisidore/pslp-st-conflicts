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

#include "probdata.hpp"
#include "conshdlr_samediff.hpp"
#include "logging.hpp"
#include "pprint.hpp"
#include "pricer.hpp"
#include "vardata.hpp"
#include <algorithm>
#include <fstream>
#include <scip/cons_knapsack.h>
#include <scip/cons_setppc.h>
#include <sstream>

// Print the problem
constexpr bool PrintProblem = false;

// Print the solution
constexpr bool PrintSolution = false;

// Write final linear program in lp file
constexpr bool WriteLp = false;

// Check all the variables
constexpr bool CheckVars = true;

using namespace pslp::bp;

Probdata::Probdata(Type type,
                   std::size_t n,
                   std::size_t m,
                   std::size_t b,
                   bool has_adj_conflicts)
    : _type(type),
      _n_stacks(m),
      _capacity(b),
      _has_adjacent_conflicts(has_adj_conflicts),
      _conflict_matrix(n, false),
      _blocking_matrix(n, false),
      _conflict_graph(n),
      _blocking_graph(n),
      _initial_layout(n, m, b),
      _conss_assign(n, nullptr)
{
}

void Probdata::set_arrivals(const std::vector<std::size_t> &arrivals)
{
    [[maybe_unused]] const auto n = n_items();

    // Items must be sorted by arrival time
    assert(std::size(arrivals) == n);
    assert(std::is_sorted(std::begin(arrivals), std::end(arrivals)));
    assert(std::all_of(std::begin(arrivals),
                       std::end(arrivals),
                       [n](auto i) { return i <= n; }));

    _arrivals = arrivals;
}

void Probdata::set_weights(const std::vector<double> &weights)
{
    const auto n = n_items();

    assert(std::size(weights) == n);

    _weights = weights;

    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            _conflict_matrix(i, j) = (_weights[i] > _weights[j]);
        }
    }
}

void Probdata::set_departures(const std::vector<double> &departures)
{
    const auto n = n_items();

    assert(std::size(departures) == n);

    _departures = departures;

    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            _blocking_matrix(i, j) = (_departures[i] > _departures[j]);
        }
    }
}

void Probdata::set_initial_positions(const std::vector<std::size_t> &init_pos)
{
    [[maybe_unused]] const auto n = n_items();
    [[maybe_unused]] const auto m = n_stacks();

    assert(std::size(init_pos) < n);
    assert(std::all_of(std::begin(init_pos),
                       std::end(init_pos),
                       [m](auto k) { return k < m; }));

    _initial_positions = init_pos;
}

void Probdata::set_conflict_matrix(const std::vector<std::vector<bool>> &matrix)
{
    const auto n = n_items();

    assert(std::size(matrix) == n);
    assert(std::all_of(std::begin(matrix),
                       std::end(matrix),
                       [n](const auto &v) { return std::size(v) == n; }));

    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            _conflict_matrix(i, j) = matrix[i][j];
        }
    }
}

void Probdata::set_blocking_matrix(const std::vector<std::vector<bool>> &matrix)
{
    const auto n = n_items();

    assert(std::size(matrix) == n);
    assert(std::all_of(std::begin(matrix),
                       std::end(matrix),
                       [n](const auto &v) { return std::size(v) == n; }));

    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            _blocking_matrix(i, j) = matrix[i][j];
        }
    }
}

SCIP_Real Probdata::dual_assign(SCIP *scip, std::size_t i, bool farkas) const
{
    return farkas
               ? SCIPgetDualfarkasSetppc(scip, _conss_assign[i])
               : SCIPgetDualsolSetppc(scip, _conss_assign[i]);
}

SCIP_Real Probdata::dual_space(SCIP *scip, bool farkas) const
{
    // This constraint exists only in the PSLP minimizing BI
    assert(_type == Type::MinNBlockingItems);

    return farkas
               ? SCIPgetDualfarkasKnapsack(scip, _cons_space)
               : SCIPgetDualsolKnapsack(scip, _cons_space);
}

bool Probdata::has_var(SCIP *scip, const std::vector<std::size_t> &s) const
{
    auto ordered_s = s;

    std::sort(std::begin(ordered_s), std::end(ordered_s));

    for (auto var : _vars)
    {
        const auto &vardata = Vardata::get(scip, var);
        auto ordered_var = vardata.items();

        std::sort(std::begin(ordered_var), std::end(ordered_var));

        if (ordered_s == ordered_var)
        {
            return true;
        }
    }

    return false;
}

bool Probdata::check_var(const std::vector<std::size_t> &s) const
{
    const auto n = n_items();
    const auto m = n_stacks();

    // Items must have a valid index between 0 and n-1
    if (!std::all_of(std::begin(s), std::end(s), [n](auto i) { return i < n; }))
    {
        LOG_ERROR("invalid index");
        return false;
    }

    // Items must be unique
    auto sorted = s;
    std::sort(std::begin(sorted), std::end(sorted));
    if (std::adjacent_find(std::begin(sorted),
                           std::end(sorted)) != std::end(sorted))
    {
        LOG_ERROR("duplicate items");
        return false;
    }

    // Items must satisfy arrival order and conflict constraints
    for (auto i = std::begin(s); i != std::end(s); ++i)
    {
        for (auto j = std::begin(s); j != i; ++j)
        {
            if (_arrivals[*i] < _arrivals[*j])
            {
                LOG_ERROR("violation of arrival order");
                return false;
            }

            if (!_has_adjacent_conflicts && _conflict_matrix(*i, *j))
            {
                LOG_ERROR("violation of distant conflict constraints");
                return false;
            }
        }

        if (_has_adjacent_conflicts &&
            i + 1 != std::end(s) &&
            _conflict_matrix(*(i + 1), *i))
        {
            LOG_ERROR("violation of adjacent conflict constraints");
            return false;
        }
    }

    // The number of items cannot exceed the capacity
    if (std::size(s) > _capacity)
    {
        LOG_ERROR("violation of capacity");
        return false;
    }

    // Initial items must be in either the same column or different columns
    for (auto i : s)
    {
        assert((_arrivals[i] == 0) == (_initial_layout.contains(i)));

        // Skip incoming items
        if (_arrivals[i] > 0)
        {
            continue;
        }

        auto k_i = _initial_layout.stack_of(i);

        for (std::size_t k = 0; k < m; ++k)
        {
            for (auto j : _initial_layout.stack(k))
            {
                bool found = std::find(std::begin(s), std::end(s), j) != std::end(s);

                if (k == k_i)
                {
                    // All initial items in k_i must be in s
                    if (!found)
                    {
                        LOG_ERROR("required initial item %zu not found (stack %zu)",
                                  j + 1,
                                  k + 1);
                        return false;
                    }
                }
                else
                {
                    // All initial items not in k_i must be not in s
                    if (found)
                    {
                        LOG_ERROR("forbidden initial item %zu found (stack %zu)",
                                  j + 1,
                                  k + 1);
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

double Probdata::get_var_obj(const std::vector<std::size_t> &s) const
{
    double obj = 0.;

    switch (_type)
    {
        case Type::MinNStacks:
            obj = 1.;
            break;

        case Type::MinNBlockingItems:
            std::size_t n_blocking = 0;
            for (auto i = std::begin(s); i != std::end(s); ++i)
            {
                for (auto j = std::begin(s); j != i; ++j)
                {
                    if (_blocking_matrix(*i, *j))
                    {
                        ++n_blocking;
                        break;
                    }
                }
            }
            obj = n_blocking;
    }

    return obj;
}

SCIP_RETCODE Probdata::add_var(SCIP *scip,
                               const std::vector<std::size_t> &s,
                               SCIP_PRICER *pricer)
{
    SCIP_VAR *var;
    std::ostringstream oss;

    // Name of the column
    oss << "x{" << std::size(_vars) + 1 << "}";

    // Get the obj coefficient of the column
    double obj = get_var_obj(s);

    if constexpr (LOG_ACTIVE_LEVEL <= LOG_LEVEL_DEBUG)
    {
        if (pricer)
        {
            LOG_DEBUG("create var <%s> (pricer: %s, obj: %g)",
                      oss.str().c_str(),
                      SCIPpricerGetName(pricer),
                      obj);
        }
        else
        {
            LOG_DEBUG("create var <%s> (obj: %g)", oss.str().c_str(), obj);
        }

        LOG_TRACE("-> items: %s",
                  pprint_n(std::size(s),
                           [&s](auto i) { return s[i] + 1; })
                      .c_str());

#ifndef NDEBUG
        if (has_var(scip, s))
        {
            LOG_DEBUG("-> redundant");
        }
#endif
    }

    assert(check_var(s));

    // Create the SCIP variable
    SCIP_CALL(SCIPcreateObjVar(scip,
                               &var,
                               oss.str().c_str(),
                               0.,
                               1.,
                               obj,
                               SCIP_VARTYPE_BINARY,
                               true, // initial
                               true, // removable
                               new Vardata(s, pricer),
                               true));

    // Avoid adding "x <= 1" constraints
    SCIP_CALL(SCIPchgVarUbLazy(scip, var, 1.));

    // Add the variable to the assignment constraints
    for (auto i : s)
    {
        SCIP_CALL(SCIPaddCoefSetppc(scip, _conss_assign[i], var));
    }

    assert(_cons_space != nullptr || _type == Type::MinNStacks);
    assert(_cons_space == nullptr || _type == Type::MinNBlockingItems);

    // This constraints exists only for the PSLP minimizing BI
    if (_cons_space)
    {
        SCIP_CALL(SCIPaddCoefKnapsack(scip, _cons_space, var, 1));
    }

    if (pricer != nullptr)
    {
        SCIP_CALL(SCIPaddPricedVar(scip, var, 1.));
    }
    else
    {
        SCIP_CALL(SCIPaddVar(scip, var));
    }

    _vars.push_back(var);

    return SCIP_OKAY;
}

SCIP_RETCODE Probdata::prepare()
{
    const auto n = n_items();
    const auto n_initial = n_initial_items();

    // Build the initial layout
    _initial_layout.clear();

    for (std::size_t i = 0; i < n_initial; ++i)
    {
        auto k = _initial_positions[i];
        _initial_layout.push(k, i);
    }

    assert(n_initial == _initial_layout.n_stored_items());

    // Initialize arrival times if not specified
    if (std::empty(_arrivals))
    {
        _arrivals.resize(n, 0);

        for (std::size_t i = n_initial; i < n; ++i)
        {
            _arrivals[i] = 1 + i - n_initial;
        }
    }

    // Arrival of initial items is 0, arrival of incoming items is > 0
    assert(std::all_of(std::begin(_arrivals),
                       std::begin(_arrivals) + n_initial,
                       [](auto t) { return t == 0; }));
    assert(std::all_of(std::begin(_arrivals) + n_initial,
                       std::end(_arrivals),
                       [](auto t) { return t > 0; }));

    // Create additional conflict edges
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            // Items are sorted by arrival time (i > j)
            // Item cannot be stacked above items arriving strictly later
            if (_arrivals[i] > _arrivals[j])
            {
                _conflict_matrix(j, i) = true;
            }
            else
            {
                assert(_arrivals[i] == _arrivals[j]);
            }

            // For each pair of initial item (i, j), i > j
            if (_arrivals[i] == 0 && _arrivals[j] == 0)
            {
                // Same stack
                if (_initial_positions[i] == _initial_positions[j])
                {
                    // Items from the same stack have a predefined order
                    // => i is above j, so j cannot be stacked above i
                    _conflict_matrix(j, i) = true;

                    if (_conflict_matrix(i, j))
                    {
                        const auto k = _initial_positions[i];

                        if (_has_adjacent_conflicts &&
                            _initial_layout.level_of(i) == _initial_layout.level_of(j) + 1)
                        {
                            _conflict_matrix(i, j) = false;

                            LOG_WARN("ignore adjacent conflict between initial items %zu and %zu in stack %zu",
                                     i + 1,
                                     j + 1,
                                     k + 1);
                        }
                        else if (!_has_adjacent_conflicts)
                        {
                            assert(_initial_layout.level_of(i) > _initial_layout.level_of(j));

                            _conflict_matrix(i, j) = false;

                            LOG_WARN("ignore distant conflict between initial items %zu and %zu in stack %zu",
                                     i + 1,
                                     j + 1,
                                     k + 1);
                        }
                    }
                }
                // Different stacks
                else
                {
                    // Items from distinct stacks cannot be put together at all
                    _conflict_matrix(j, i) = true;
                    _conflict_matrix(i, j) = true;
                }
            }
        }
    }

    // Remove redundant blocking constraints (dominated by conflict constraints)
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i != j && _conflict_matrix(i, j) && _blocking_matrix(i, j))
            {
                _blocking_matrix(i, j) = false;
            }
        }
    }

    _conflict_graph.clear_edges();
    _blocking_graph.clear_edges();

    // Build conflict constraint graph
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i != j && _conflict_matrix(i, j))
            {
                _conflict_graph.add_edge(i, j);
            }
        }
    }

    // Build blocking constraint graph
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i != j && _blocking_matrix(i, j))
            {
                _blocking_graph.add_edge(i, j);
            }
        }
    }

    // Build the conflict graph complement
    _conflict_graph_complement = _conflict_graph.complement();

#ifndef NDEBUG
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j)
            {
                continue;
            }

            assert(_conflict_matrix(i, j) == _conflict_graph.has_edge(i, j));
            assert(_blocking_matrix(i, j) == _blocking_graph.has_edge(i, j));
            assert(_conflict_matrix(i, j) != _conflict_graph_complement.has_edge(i, j));
        }
    }
#endif

    return SCIP_OKAY;
}

SCIP_RETCODE Probdata::initialize(SCIP *scip)
{
    // Some structures need to be initialized
    SCIP_CALL(prepare());

    LOG_INFO("problem has %s conflicts", _has_adjacent_conflicts
                                             ? "adjacent"
                                             : "distant");

    // The objective value is an integer
    SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));
    SCIP_CALL(SCIPsetObjIntegral(scip));

    // Assignment constraints
    for (std::size_t i = 0; i < std::size(_conss_assign); ++i)
    {
        SCIP_CONS *cons;
        std::ostringstream oss;

        oss << "assign{" << i + 1 << "}";

        SCIP_CALL(SCIPcreateConsBasicSetpart(scip,
                                             &cons,
                                             oss.str().c_str(),
                                             0,
                                             nullptr));

        SCIP_CALL(SCIPsetConsModifiable(scip, cons, true));
        SCIP_CALL(SCIPaddCons(scip, cons));

        _conss_assign[i] = cons;
    }

    // This constraint exists only for the PSLP minimizing BI
    if (_type == Type::MinNBlockingItems)
    {
        SCIP_CONS *cons;

        SCIP_CALL(SCIPcreateConsBasicKnapsack(scip,
                                              &cons,
                                              "space",
                                              0,
                                              nullptr,
                                              nullptr,
                                              _n_stacks));

        SCIP_CALL(SCIPsetConsModifiable(scip, cons, true));
        SCIP_CALL(SCIPaddCons(scip, cons));

        _cons_space = cons;
    }

    const auto n_initial = n_initial_items();

    // Add constraints for initial items
    for (std::size_t i = 0; i < n_initial; ++i)
    {
        assert(_arrivals[i] == 0);

        for (std::size_t j = 0; j < i; ++j)
        {
            if (i == j)
            {
                continue;
            }

            SCIP_CONS *cons;
            std::ostringstream oss;
            ConshdlrSamediff::Type type;

            if (_initial_positions[i] == _initial_positions[j])
            {
                oss << "same{" << i + 1 << ',' << j + 1 << '}';
                type = ConshdlrSamediff::Type::Same;
            }
            else
            {
                oss << "diff{" << i + 1 << ',' << j + 1 << '}';
                type = ConshdlrSamediff::Type::Diff;
            }

            SCIP_CALL(ConshdlrSamediff::create_cons(scip,
                                                    &cons,
                                                    oss.str().c_str(),
                                                    type,
                                                    i,
                                                    j,
                                                    nullptr));

            SCIP_CALL(SCIPaddCons(scip, cons));
            SCIP_CALL(SCIPreleaseCons(scip, &cons));
        }
    }

    // Activate pricers
    const auto &pricers = Pricer::list();

    for (const auto &info : pricers)
    {
        // Request for the pricers/<name>/activate parameter
        std::ostringstream oss;
        oss << "pricers/" << info.name << "/activate";

        SCIP_Bool activate;
        SCIP_CALL(SCIPgetBoolParam(scip, oss.str().c_str(), &activate));

        if (activate)
        {
            auto *pricer = SCIPfindPricer(scip, info.name);

            if (pricer)
            {
                LOG_INFO("activate pricer <%s>", info.name);
                SCIP_CALL(SCIPactivatePricer(scip, pricer));
            }
            else
            {
                LOG_ERROR("cannot find pricer <%s>", info.name);
                return SCIP_ERROR;
            }
        }
    }

    return SCIP_OKAY;
}

SCIP_RETCODE Probdata::release(SCIP *scip)
{
    for (auto &cons : _conss_assign)
    {
        SCIP_CALL(SCIPreleaseCons(scip, &cons));
    }

    if (_cons_space)
    {
        SCIP_CALL(SCIPreleaseCons(scip, &_cons_space));
    }

    for (auto &var : _vars)
    {
        SCIP_CALL(SCIPreleaseVar(scip, &var));
    }

    return SCIP_OKAY;
}

SCIP_RETCODE Probdata::build_sol(SCIP *scip,
                                 SCIP_SOL *sol,
                                 std::vector<std::size_t> &stack_of,
                                 std::vector<std::size_t> &level_of) const
{
    assert(sol != nullptr);

    const auto n = n_items();
    const auto b = capacity();
    const auto n_initial = n_initial_items();
    std::vector<bool> used_stacks(n, false);

    stack_of.clear();
    level_of.clear();

    stack_of.resize(n, n);
    level_of.resize(n, b);

    // First build columns having initial items
    for (auto var : _vars)
    {
        if (SCIPgetSolVal(scip, sol, var) < 0.5)
        {
            continue;
        }

        const auto &s = Vardata::get(scip, var).items();

        if (std::empty(s))
        {
            LOG_ERROR("no empty stack should be selected");
            return SCIP_ERROR;
        }

        // Only columns having initial items
        if (s[0] >= n_initial)
        {
            continue;
        }

        // Stack number is predefined by initial items
        std::size_t k = _initial_positions[s[0]];

        if (used_stacks[k])
        {
            LOG_ERROR("stack %zu already used", k + 1);
            return SCIP_ERROR;
        }

        used_stacks[k] = true;

        for (std::size_t l = 0; l < std::size(s); ++l)
        {
            const auto i = s[l];

            if (i >= n)
            {
                LOG_ERROR("invalid item index");
                return SCIP_ERROR;
            }

            if (stack_of[i] != n)
            {
                LOG_ERROR("item i %zu already assigned", i + 1);
                return SCIP_ERROR;
            }

            stack_of[i] = k;
            level_of[i] = l;
        }
    }

    // Build the rest of the columns
    for (auto var : _vars)
    {
        if (SCIPgetSolVal(scip, sol, var) < 0.5)
        {
            continue;
        }

        const auto &s = Vardata::get(scip, var).items();

        // Only columns having non-initial items
        if (s[0] < n_initial)
        {
            continue;
        }

        // Pick any available stack
        std::size_t k = 0;
        while (used_stacks[k])
        {
            ++k;
            if (k >= std::size(used_stacks))
            {
                LOG_ERROR("cannot find any available stack");
                return SCIP_ERROR;
            }
        }

        used_stacks[k] = true;

        for (std::size_t l = 0; l < std::size(s); ++l)
        {
            const auto i = s[l];

            if (i >= n)
            {
                LOG_ERROR("invalid item index");
                return SCIP_ERROR;
            }

            if (stack_of[i] != n)
            {
                LOG_ERROR("item i %zu already assigned", i + 1);
                return SCIP_ERROR;
            }

            stack_of[i] = k;
            level_of[i] = l;
        }
    }

    return SCIP_OKAY;
}

SCIP_RETCODE Probdata::check_sol(SCIP *scip,
                                 SCIP_Real obj,
                                 const std::vector<std::size_t> &stack_of,
                                 const std::vector<std::size_t> &level_of) const
{
    const auto n = n_items();
    const auto m = n_stacks();
    const auto b = capacity();
    const auto n_initial = n_initial_items();

    if (std::size(stack_of) != n || std::size(level_of) != n)
    {
        LOG_ERROR("invalid solution size");
        return SCIP_ERROR;
    }

    std::vector<bool> used_stacks(n, false);

    for (auto k : stack_of)
    {
        used_stacks[k] = true;
    }

    std::size_t n_used_stacks = std::count(std::begin(used_stacks),
                                           std::end(used_stacks),
                                           true);

    // Check whether the number of used stacks exceeds the number of stacks
    if (n_used_stacks > m)
    {
        switch (_type)
        {
            case Type::MinNStacks:
                LOG_WARN("too many stacks used (%zu used / %zu available)",
                         n_used_stacks,
                         m);
                break;

            case Type::MinNBlockingItems:
                LOG_ERROR("too many stacks used (%zu used / %zu available)",
                          n_used_stacks,
                          m);
                return SCIP_ERROR;
        }
    }

    // Check whether each item is assigned to a stack
    for (std::size_t i = 0; i < n; ++i)
    {
        if (stack_of[i] >= n)
        {
            LOG_ERROR("item %zu not assigned", i + 1);
            return SCIP_ERROR;
        }
    }

    // Check initial items stacks
    for (std::size_t i = 0; i < n_initial; ++i)
    {
        const auto k = _initial_positions[i];
        const auto k2 = stack_of[i];
        if (k != k2)
        {
            LOG_ERROR("item %zu must be in stack %zu, not %zu",
                      i + 1,
                      k + 1,
                      k2 + 1);
            return SCIP_ERROR;
        }
    }

    // Check initial items levels
    for (std::size_t i = 0; i < n_initial; ++i)
    {
        const auto l = _initial_layout.level_of(i);
        const auto l2 = level_of[i];
        if (l != l2)
        {
            LOG_ERROR("item %zu must be at level %zu, not %zu",
                      i + 1,
                      l + 1,
                      l2 + 1);
            return SCIP_ERROR;
        }
    }

    // Check initial items
    for (std::size_t i = 0; i < n_initial; ++i)
    {
        for (std::size_t j = 0; j < n_initial; ++j)
        {
            if (i != j &&
                _initial_positions[i] == _initial_positions[j] &&
                stack_of[i] != stack_of[j])
            {
                LOG_ERROR("item %zu and %zu must be together", i + 1, j + 1);
                return SCIP_ERROR;
            }
        }
    }

    // Check arrival order
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            // Conflicts are adjacent
            if (stack_of[i] == stack_of[j] &&
                ((level_of[i] > level_of[j] &&
                  _arrivals[i] < _arrivals[j]) ||
                 (level_of[j] > level_of[i] &&
                  _arrivals[j] < _arrivals[i])))
            {
                LOG_ERROR("invalid order between items %zu and %zu",
                          i + 1,
                          j + 1);
                return SCIP_ERROR;
            }
        }
    }

    // Check conflict constraint violation
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j || !_conflict_matrix(i, j) || stack_of[i] != stack_of[j])
            {
                continue;
            }

            if (_has_adjacent_conflicts && level_of[i] == level_of[j] + 1)
            {
                LOG_ERROR("item %zu cannot be put on top of %zu", i + 1, j + 1);
                return SCIP_ERROR;
            }
            else if (!_has_adjacent_conflicts && level_of[i] > level_of[j])
            {
                LOG_ERROR("item %zu cannot be put above %zu", i + 1, j + 1);
                return SCIP_ERROR;
            }
        }
    }

    // Check capacity constraint
    for (std::size_t k = 0; k < n_stacks(); ++k)
    {
        std::size_t h = 0;

        for (auto l : stack_of)
        {
            if (k == l)
            {
                ++h;
            }
        }

        if (h > b)
        {
            LOG_ERROR("stack %zu of height %zu violates capacity %zu",
                      k + 1,
                      h,
                      b);
            return SCIP_ERROR;
        }
    }

    // Check objective value
    switch (_type)
    {
        case Type::MinNStacks:
            if (!SCIPisEQ(scip, SCIP_Real(n_used_stacks), obj))
            {
                LOG_ERROR("objectives values are different: %zu != %g",
                          n_used_stacks,
                          obj);
                return SCIP_ERROR;
            }
            break;

        case Type::MinNBlockingItems:
        {
            std::vector<bool> is_blocking(n, false);

            for (std::size_t i = 0; i < n; ++i)
            {
                for (std::size_t j = 0; j < n; ++j)
                {
                    if (i == j)
                    {
                        continue;
                    }

                    if (stack_of[i] == stack_of[j] &&
                        level_of[i] > level_of[j] &&
                        _blocking_matrix(i, j))
                    {
                        is_blocking[i] = true;
                    }
                }
            }

            auto n_blocking = std::count(std::begin(is_blocking),
                                         std::end(is_blocking),
                                         true);

            if (!SCIPisEQ(scip, SCIP_Real(n_blocking), obj))
            {
                LOG_ERROR("objectives values are different: %zu != %g",
                          n_blocking,
                          obj);
                return SCIP_ERROR;
            }
        }
        break;
    }

    return SCIP_OKAY;
}

SCIP_RETCODE Probdata::print_sol(SCIP *scip,
                                 SCIP_SOL *sol,
                                 std::ostream &os) const
{
    assert(sol != nullptr);

    const auto n = n_items();
    const auto m = n_stacks();
    const auto b = capacity();
    Layout layout(n, n, b);
    std::size_t k = 0;

    for (auto var : _vars)
    {
        if (SCIPgetSolVal(scip, sol, var) > 0.5)
        {
            const auto &vardata = Vardata::get(scip, var);
            const auto &s = vardata.items();

            if constexpr (LOG_ACTIVE_LEVEL <= LOG_LEVEL_TRACE)
            {
                std::ostringstream oss;
                for (auto i : s)
                {
                    oss << ' ' << i + 1;
                }

                LOG_DEBUG("var %s:%s", SCIPvarGetName(var), oss.str().c_str());
            }

            for (auto i : s)
            {
                assert(!layout.contains(i));

                layout.push(k, i);
            }

            ++k;
        }
    }

    LOG_DEBUG("final layout (indices):");
    layout.print(
        [this](auto i) {
            return i + 1;
        },
        os);

    if (has_departures())
    {
        LOG_DEBUG("final layout (departures):");
        layout.print(
            [this](auto i) {
                return departure(i);
            },
            os);
    }

    if (has_weights())
    {
        LOG_DEBUG("final layout (weights):");
        layout.print(
            [this](auto i) {
                return weight(i);
            },
            os);
    }

    const auto n_used_stacks = layout.n_used_stacks();

    if (n_used_stacks > m)
    {
        LOG_WARN("used more stacks (%zu) than available (%zu)",
                 n_used_stacks,
                 m);
    }

    return SCIP_OKAY;
}

SCIP_RETCODE Probdata::create_prob(SCIP *scip,
                                   const char *name,
                                   Type type,
                                   std::size_t n_items,
                                   std::size_t n_stacks,
                                   std::size_t capacity,
                                   bool has_adj_conflicts)
{
    auto *probdata = new Probdata(type,
                                  n_items,
                                  n_stacks,
                                  capacity,
                                  has_adj_conflicts);

    SCIP_CALL(SCIPcreateObjProb(scip, name, probdata, true));

    return SCIP_OKAY;
}

Probdata &Probdata::get(SCIP *scip)
{
    auto *probdata = static_cast<Probdata *>(SCIPgetObjProbData(scip));
    assert(probdata != nullptr);
    return *probdata;
}

SCIP_RETCODE Probdata::scip_delorig(SCIP *scip)
{
    return release(scip);
}

SCIP_RETCODE Probdata::scip_trans(SCIP *scip,
                                  scip::ObjProbData **objprobdata,
                                  SCIP_Bool *deleteobject)
{
    assert(objprobdata != nullptr);
    assert(deleteobject != nullptr);

    auto *t_probdata = new Probdata(*this);

    SCIP_CALL(SCIPtransformVars(scip,
                                std::size(t_probdata->_vars),
                                std::data(t_probdata->_vars),
                                std::data(t_probdata->_vars)));

    SCIP_CALL(SCIPtransformConss(scip,
                                 std::size(t_probdata->_conss_assign),
                                 std::data(t_probdata->_conss_assign),
                                 std::data(t_probdata->_conss_assign)));

    assert(t_probdata->_cons_space != nullptr ||
           t_probdata->_type == Type::MinNStacks);
    assert(t_probdata->_cons_space == nullptr ||
           t_probdata->_type == Type::MinNBlockingItems);

    if (t_probdata->_cons_space)
    {
        SCIP_CALL(SCIPtransformCons(scip,
                                    t_probdata->_cons_space,
                                    &t_probdata->_cons_space));
    }

    *objprobdata = t_probdata;
    *deleteobject = true;

    return SCIP_OKAY;
}

SCIP_RETCODE Probdata::scip_deltrans(SCIP *scip)
{
    return release(scip);
}

SCIP_RETCODE Probdata::scip_initsol([[maybe_unused]] SCIP *scip)
{
    if constexpr (PrintProblem && LOG_ACTIVE_LEVEL <= LOG_LEVEL_TRACE)
    {
        const auto n = n_items();
        const auto m = n_stacks();
        const auto b = capacity();

        LOG_TRACE("n_items: %zu", n);
        LOG_TRACE("n_stacks: %zu", m);
        LOG_TRACE("capacity: %zu", b);

        if (_initial_layout.n_stored_items() > 0)
        {
            LOG_TRACE("initial layout (indices):");
            _initial_layout.print([](auto i) { return i + 1; });

            if (has_departures())
            {
                LOG_TRACE("initial layout (departures):");
                _initial_layout.print([this](auto i) { return departure(i); });
            }

            if (has_weights())
            {
                LOG_TRACE("initial layout (weights):");
                _initial_layout.print([this](auto i) { return weight(i); });
            }
        }
        else
        {
            LOG_TRACE("empty initial layout");
        }

        if (_type == Type::MinNBlockingItems)
        {
            LOG_TRACE("blocking graph:");
            _blocking_graph.print();
        }

        LOG_TRACE("conflict graph:");
        _conflict_graph.print();

        LOG_TRACE("conflicts are %s", has_adjacent_conflicts()
                                          ? "adjacent"
                                          : "distant");
    }

    return SCIP_OKAY;
}

SCIP_RETCODE Probdata::scip_exitsol(SCIP *scip,
                                    [[maybe_unused]] SCIP_Bool restart)
{
    if constexpr (PrintSolution && LOG_ACTIVE_LEVEL <= LOG_LEVEL_DEBUG)
    {
        SCIP_SOL *sol = SCIPgetBestSol(scip);
        if (sol)
        {
            LOG_DEBUG("solution:");
            SCIP_CALL(print_sol(scip, sol));
        }
        else
        {
            LOG_DEBUG("no solution");
        }
    }

    if constexpr (WriteLp)
    {
        std::ostringstream oss;

        oss << SCIPgetProbName(scip) << ".lp";

        LOG_INFO("write transformed master problem in <%s>", oss.str().c_str());

        SCIP_CALL(SCIPwriteTransProblem(scip,
                                        oss.str().c_str(),
                                        nullptr,
                                        false));
    }

    SCIP_SOL **sols = SCIPgetSols(scip);
    std::size_t n_sols = SCIPgetNSols(scip);

    for (std::size_t t = 0; t < n_sols; ++t)
    {
        LOG_INFO("checking solution %zu / %zu", t + 1, n_sols);

        SCIP_SOL *sol = sols[t];
        SCIP_Real obj = SCIPgetSolOrigObj(scip, sol);
        std::vector<std::size_t> stack_of;
        std::vector<std::size_t> level_of;

        SCIP_CALL(build_sol(scip, sol, stack_of, level_of));

        SCIP_RETCODE retcode = check_sol(scip, obj, stack_of, level_of);
        if (retcode != SCIP_OKAY)
        {
            SCIP_CALL(print_sol(scip, sol));
            return retcode;
        }
    }

    LOG_INFO("solutions checked");

    // Detect if pricers generate invalid variables
    if constexpr (CheckVars)
    {
        for (auto var : _vars)
        {
            const auto &s = Vardata::get(scip, var).items();
            if (!check_var(s))
            {
                LOG_ERROR("invalid variable");
                return SCIP_ERROR;
            }
        }

        LOG_INFO("variables checked");
    }

    return SCIP_OKAY;
}
