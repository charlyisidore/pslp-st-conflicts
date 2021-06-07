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

#include "mip_cplex.hpp"
#include "logging.hpp"
#include "symmetric_matrix.hpp"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <ilconcert/ilothread.h>
#include <numeric>

constexpr int MIPDisplay = 3;
constexpr int MIPInterval = 0;
constexpr int Threads = 0;
constexpr double TimeLimit = 3600.;
constexpr double WorkMem = 4096.;
constexpr int MIPStrategyFile = 3;
constexpr double MIPLimitsTreeMemory = 1.e75;
constexpr bool UseMipStart = true;
constexpr bool UseIncumbentCallback = false;

using namespace pslp::mip;

using Clock = std::chrono::steady_clock;
using TimePoint = std::chrono::time_point<Clock>;
using Duration = std::chrono::duration<double>;

static const char *to_string(IloAlgorithm::Status status)
{
    switch (status)
    {
        case IloAlgorithm::Unknown:
            return "Unknown";
        case IloAlgorithm::Feasible:
            return "Feasible";
        case IloAlgorithm::Optimal:
            return "Optimal";
        case IloAlgorithm::Infeasible:
            return "Infeasible";
        case IloAlgorithm::Unbounded:
            return "Unbounded";
        case IloAlgorithm::InfeasibleOrUnbounded:
            return "InfeasibleOrUnbounded";
        case IloAlgorithm::Error:
            return "Error";
        default:
            return "?";
    }
}

static const char *to_string(IloCplex::CplexStatus status)
{
    switch (status)
    {
        case IloCplex::Unknown:
            return "Unknown";
        case IloCplex::Optimal:
            return "Optimal";
        case IloCplex::Unbounded:
            return "Unbounded";
        case IloCplex::Infeasible:
            return "Infeasible";
        case IloCplex::InfOrUnbd:
            return "InfOrUnbd";
        case IloCplex::OptimalInfeas:
            return "OptimalInfeas";
        case IloCplex::NumBest:
            return "NumBest";
        case IloCplex::FeasibleRelaxedSum:
            return "FeasibleRelaxedSum";
        case IloCplex::OptimalRelaxedSum:
            return "OptimalRelaxedSum";
        case IloCplex::FeasibleRelaxedInf:
            return "FeasibleRelaxedInf";
        case IloCplex::OptimalRelaxedInf:
            return "OptimalRelaxedInf";
        case IloCplex::FeasibleRelaxedQuad:
            return "FeasibleRelaxedQuad";
        case IloCplex::OptimalRelaxedQuad:
            return "OptimalRelaxedQuad";
        case IloCplex::AbortRelaxed:
            return "AbortRelaxed";
        case IloCplex::AbortObjLim:
            return "AbortObjLim";
        case IloCplex::AbortPrimObjLim:
            return "AbortPrimObjLim";
        case IloCplex::AbortDualObjLim:
            return "AbortDualObjLim";
        case IloCplex::AbortItLim:
            return "AbortItLim";
        case IloCplex::AbortTimeLim:
            return "AbortTimeLim";
        case IloCplex::AbortDetTimeLim:
            return "AbortDetTimeLim";
        case IloCplex::AbortUser:
            return "AbortUser";
        case IloCplex::OptimalFaceUnbounded:
            return "OptimalFaceUnbounded";
        case IloCplex::OptimalTol:
            return "OptimalTol";
        case IloCplex::SolLim:
            return "SolLim";
        case IloCplex::PopulateSolLim:
            return "PopulateSolLim";
        case IloCplex::NodeLimFeas:
            return "NodeLimFeas";
        case IloCplex::NodeLimInfeas:
            return "NodeLimInfeas";
        case IloCplex::FailFeas:
            return "FailFeas";
        case IloCplex::FailInfeas:
            return "FailInfeas";
        case IloCplex::MemLimFeas:
            return "MemLimFeas";
        case IloCplex::MemLimInfeas:
            return "MemLimInfeas";
        case IloCplex::FailFeasNoTree:
            return "FailFeasNoTree";
        case IloCplex::FailInfeasNoTree:
            return "FailInfeasNoTree";
        case IloCplex::ConflictFeasible:
            return "ConflictFeasible";
        case IloCplex::ConflictMinimal:
            return "ConflictMinimal";
        case IloCplex::ConflictAbortContradiction:
            return "ConflictAbortContradiction";
        case IloCplex::ConflictAbortTimeLim:
            return "ConflictAbortTimeLim";
        case IloCplex::ConflictAbortDetTimeLim:
            return "ConflictAbortDetTimeLim";
        case IloCplex::ConflictAbortItLim:
            return "ConflictAbortItLim";
        case IloCplex::ConflictAbortNodeLim:
            return "ConflictAbortNodeLim";
        case IloCplex::ConflictAbortObjLim:
            return "ConflictAbortObjLim";
        case IloCplex::ConflictAbortMemLim:
            return "ConflictAbortMemLim";
        case IloCplex::ConflictAbortUser:
            return "ConflictAbortUser";
        case IloCplex::Feasible:
            return "Feasible";
        case IloCplex::OptimalPopulated:
            return "OptimalPopulated";
        case IloCplex::OptimalPopulatedTol:
            return "OptimalPopulatedTol";
        case IloCplex::RelaxationUnbounded:
            return "RelaxationUnbounded";
        case IloCplex::FirstOrder:
            return "FirstOrder";
        case IloCplex::MultiObjOptimal:
            return "MultiObjOptimal";
        case IloCplex::MultiObjNonOptimal:
            return "MultiObjNonOptimal";
        case IloCplex::MultiObjInfeasible:
            return "MultiObjInfeasible";
        case IloCplex::MultiObjUnbounded:
            return "MultiObjUnbounded";
        case IloCplex::MultiObjInfOrUnbd:
            return "MultiObjInfOrUnbd";
        case IloCplex::MultiObjStopped:
            return "MultiObjStopped";
        default:
            return "?";
    }
}

static const char *to_string(IloCplex::IncumbentCallbackI::SolutionSource source)
{
    switch (source)
    {
        case IloCplex::IncumbentCallbackI::NodeSolution:
            return "NodeSolution";
        case IloCplex::IncumbentCallbackI::HeuristicSolution:
            return "HeuristicSolution";
        case IloCplex::IncumbentCallbackI::UserSolution:
            return "UserSolution";
        case IloCplex::IncumbentCallbackI::MIPStartSolution:
            return "MIPStartSolution";
        default:
            return "?";
    }
}

/**
 * Callback for every incumbent solution.
 */
class IncumbentCallback : public IloCplex::IncumbentCallbackI
{
public:
    using IloCplex::MIPInfoCallbackI::getNcuts;
    struct Incumbent
    {
        TimePoint time_point;
        IloNum obj_value;
        IloNum mip_relative_gap;
        IloInt64 n_nodes;
        IloInt64 n_remaining_nodes;
        IloCplex::IncumbentCallbackI::SolutionSource solution_source;
    };

    IncumbentCallback(IloEnv env,
                      std::vector<Incumbent> &incumbents,
                      IloFastMutex &mutex)
        : IloCplex::IncumbentCallbackI(env),
          _incumbents(incumbents),
          _mutex(mutex) {}

    void main() override
    {
        _mutex.lock();
        auto time_point = Clock::now();
        auto obj_value = getObjValue();
        auto mip_relative_gap = getMIPRelativeGap();
        auto n_nodes = getNnodes64();
        auto n_remaining_nodes = getNremainingNodes64();
        auto solution_source = getSolutionSource();
        _incumbents.push_back({time_point,
                               obj_value,
                               mip_relative_gap,
                               n_nodes,
                               n_remaining_nodes,
                               solution_source});
        _mutex.unlock();
    }

    IloCplex::CallbackI *duplicateCallback() const
    {
        return (new (getEnv()) IncumbentCallback(*this));
    }

private:
    std::vector<Incumbent> &_incumbents;
    IloFastMutex &_mutex;
};

MipCplex::~MipCplex()
{
    _env.end();
}

MipCplex::Problem MipCplex::read(const nlohmann::json &json)
{
    Problem prob;
    std::size_t n = 0;

    auto m_n_items = json.find("n_items");
    auto m_n_stacks = json.find("n_stacks");
    auto m_capacity = json.find("capacity");
    auto m_arrivals = json.find("arrivals");
    auto m_weights = json.find("weights");
    auto m_conflict_graph = json.find("conflict_graph");
    auto m_conflict_matrix = json.find("conflict_matrix");
    auto m_initial_positions = json.find("initial_positions");

    // "n_items" is mandatory
    if (m_n_items == json.end())
    {
        LOG_ERROR("missing value: n_items");
        std::abort();
    }

    n = m_n_items->get<std::size_t>();

    prob.n_items = n;

    // n_stacks = n_items by default
    if (m_n_stacks == json.end())
    {
        prob.n_stacks = n;
    }
    else
    {
        prob.n_stacks = m_n_stacks->get<std::size_t>();
    }

    // capacity = n_items by default
    if (m_capacity == json.end())
    {
        prob.capacity = n;
    }
    else
    {
        prob.capacity = m_capacity->get<std::size_t>();
    }

    // Quick bound checking

    if (m_arrivals != json.end() && m_arrivals->size() != n)
    {
        LOG_ERROR("wrong size: arrivals");
        std::abort();
    }

    if (m_weights != json.end() && m_weights->size() != n)
    {
        LOG_ERROR("wrong size: weights");
        std::abort();
    }

    if (m_initial_positions != json.end() && m_initial_positions->size() > n)
    {
        LOG_ERROR("wrong size: initial_positions");
        std::abort();
    }

    if (m_conflict_matrix != json.end() && m_conflict_matrix->size() != n)
    {
        LOG_ERROR("wrong size: conflict_matrix");
        std::abort();
    }

    prob.conflict_matrix.initialize(n, false);

    // Set arrival times
    if (m_arrivals != json.end())
    {
        prob.arrivals.resize(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            prob.arrivals[i] = (*m_arrivals)[i].get<std::size_t>();
        }
    }

    // Set weights
    if (m_weights != json.end())
    {
        prob.weights.resize(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            prob.weights[i] = (*m_weights)[i].get<double>();
        }

        for (std::size_t i = 0; i < n; ++i)
        {
            for (std::size_t j = 0; j < n; ++j)
            {
                prob.conflict_matrix(i, j) = (prob.weights[i] > prob.weights[j]);
            }
        }
    }
    // Set conflict graph
    else if (m_conflict_graph != json.end())
    {
        for (auto edge : *m_conflict_graph)
        {
            if (edge.size() != 2)
            {
                LOG_ERROR("malformed edge in conflict graph");
                std::abort();
            }
            std::size_t i = edge[0];
            std::size_t j = edge[1];
            if (i <= 0 || i > n || j <= 0 || j > n)
            {
                LOG_ERROR("invalid edge (%zu, %zu) in conflict graph", i, j);
                std::abort();
            }
            prob.conflict_matrix(i - 1, j - 1) = true;
        }
    }
    // Set conflict matrix
    else if (m_conflict_matrix != json.end())
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            for (std::size_t j = 0; j < n; ++j)
            {
                prob.conflict_matrix(i, j) = (*m_conflict_matrix)[i][j].get<int>() != 0;
            }
        }
    }

    // Set position of initial items
    if (m_initial_positions != json.end())
    {
        std::size_t n_fix = m_initial_positions->size();
        prob.initial_positions.resize(n_fix);
        for (std::size_t i = 0; i < n_fix; ++i)
        {
            prob.initial_positions[i] = (*m_initial_positions)[i].get<std::size_t>() - 1;
        }
    }

    // Preprocess
    const auto n_initial = std::size(prob.initial_positions);

    // Build the initial layout
    prob.initial_layout.initialize(n, n, prob.capacity);

    for (std::size_t i = 0; i < n_initial; ++i)
    {
        auto k = prob.initial_positions[i];
        prob.initial_layout.push(k, i);
    }

    assert(n_initial == prob.initial_layout.n_stored_items());

    // Initialize arrival times if not specified
    if (std::empty(prob.arrivals))
    {
        prob.arrivals.resize(n, 0);

        for (std::size_t i = n_initial; i < n; ++i)
        {
            prob.arrivals[i] = 1 + i - n_initial;
        }
    }

    // Arrival of initial items is 0, arrival of incoming items is > 0
    assert(std::all_of(std::begin(prob.arrivals),
                       std::begin(prob.arrivals) + n_initial,
                       [](auto t) { return t == 0; }));
    assert(std::all_of(std::begin(prob.arrivals) + n_initial,
                       std::end(prob.arrivals),
                       [](auto t) { return t > 0; }));

    // Create additional conflict edges
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            // Items are sorted by arrival time (i > j)
            // Item cannot be stacked above items arriving strictly later
            if (prob.arrivals[i] > prob.arrivals[j])
            {
                prob.conflict_matrix(j, i) = true;
            }
            else
            {
                assert(prob.arrivals[i] == prob.arrivals[j]);
            }

            // For each pair of initial item (i, j), i > j
            if (prob.arrivals[i] == 0 && prob.arrivals[j] == 0)
            {
                // Same stack
                if (prob.initial_positions[i] == prob.initial_positions[j])
                {
                    // Items from the same stack have a predefined order
                    // => i is above j, so j cannot be stacked above i
                    prob.conflict_matrix(j, i) = true;

                    if (prob.conflict_matrix(i, j))
                    {
                        const auto k = prob.initial_positions[i];

                        // Conflicts are adjacent
                        if (prob.initial_layout.level_of(i) == prob.initial_layout.level_of(j) + 1)
                        {
                            prob.conflict_matrix(i, j) = false;

                            LOG_WARN("ignore adjacent conflict between initial items %zu and %zu in stack %zu",
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
                    prob.conflict_matrix(j, i) = true;
                    prob.conflict_matrix(i, j) = true;
                }
            }
        }
    }

    // Build conflict constraint graph
    prob.conflict_graph.initialize(n);

    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i != j && prob.conflict_matrix(i, j))
            {
                prob.conflict_graph.add_edge(i, j);
            }
        }
    }

#ifndef NDEBUG
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j)
            {
                continue;
            }

            assert(prob.conflict_matrix(i, j) == prob.conflict_graph.has_edge(i, j));
        }
    }
#endif

    return prob;
}

bool MipCplex::is_stacking_transitive(const Problem &prob) const
{
    const auto n = prob.conflict_matrix.n_rows();

    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (i == j)
            {
                continue;
            }

            if (prob.conflict_matrix(i, j))
            {
                continue;
            }

            for (std::size_t k = 0; k < n; ++k)
            {
                if (i == k || j == k)
                {
                    continue;
                }

                // ... && !prob.conflict_matrix(i, j)
                if (!prob.conflict_matrix(j, k) && prob.conflict_matrix(i, k))
                {
                    return false;
                }
            }
        }
    }

    return true;
}

void MipCplex::transform_conflicts(Problem &prob)
{
    const auto n = prob.conflict_matrix.n_rows();
    const auto n_initial = std::size(prob.initial_positions);
    bool has_max_set = true;

    assert(prob.type == Type::MinNStacks);
    assert(is_stacking_transitive(prob));

    while (has_max_set)
    {
        has_max_set = false;
        for (std::size_t k = n_initial; k < n; ++k)
        {
            std::vector<std::size_t> max_set = {k};
            for (std::size_t i = n_initial; i < n; ++i)
            {
                if (i == k)
                {
                    continue;
                }
                bool include = true;
                for (auto j : max_set)
                {
                    if (prob.conflict_matrix(i, j) || prob.conflict_matrix(j, i))
                    {
                        include = false;
                        break;
                    }
                }
                if (include)
                {
                    max_set.push_back(i);
                }
            }
            if (std::size(max_set) >= 2)
            {
                has_max_set = true;
                std::sort(std::begin(max_set), std::end(max_set));
                for (std::size_t u = 0; u < std::size(max_set); ++u)
                {
                    const auto i = max_set[u];
                    for (std::size_t v = 0; v < u; ++v)
                    {
                        const auto j = max_set[v];
                        prob.conflict_matrix(i, j) = true;
                    }
                }
            }
        }
    }
}

void MipCplex::heuristic(const Problem &prob)
{
    const auto n_initial = std::size(prob.initial_positions);
    const auto n_incoming = prob.n_items - n_initial;
    std::vector<std::size_t> incoming_items(n_incoming);

    _heuristic_sol.clear();

    // Place the initial items
    for (std::size_t i = 0; i < n_initial; ++i)
    {
        const auto k = prob.initial_positions[i];
        if (k >= std::size(_heuristic_sol))
        {
            _heuristic_sol.resize(k + 1);
        }
        _heuristic_sol[k].push_back(i);
    }

    // Sort incoming items
    std::iota(std::begin(incoming_items), std::end(incoming_items), n_initial);
    std::sort(std::begin(incoming_items),
              std::end(incoming_items),
              [&prob](auto i, auto j) {
                  return prob.conflict_graph.degree(i) >
                         prob.conflict_graph.degree(j);
              });

    // Insert each incoming item in the first available stack
    for (auto i : incoming_items)
    {
        bool placed = false;

        // Try stacks one by one until the item is placed
        for (std::size_t k = 0; !placed; ++k)
        {
            // Add empty stacks if necessary
            if (k >= std::size(_heuristic_sol))
            {
                _heuristic_sol.resize(k + 1);
            }

            const auto h = std::size(_heuristic_sol[k]);

            // Skip full stacks
            if (h >= prob.capacity)
            {
                assert(h == prob.capacity);
                continue;
            }

            // Conflicts are adjacent
            // Find the lowest level not causing conflicts
            for (std::size_t l = 0; l <= h; ++l)
            {
                if (l > 0)
                {
                    // Item just below
                    const auto j = _heuristic_sol[k][l - 1];
                    if (prob.conflict_matrix(i, j))
                    {
                        continue;
                    }
                }

                if (l < h)
                {
                    // Item just above
                    const auto j = _heuristic_sol[k][l];
                    if (prob.conflict_matrix(j, i))
                    {
                        continue;
                    }
                }

                _heuristic_sol[k].insert(std::begin(_heuristic_sol[k]) + l, i);
                placed = true;
                break;
            }
        }
    }

#ifndef NDEBUG
    // Conflicts are adjacent
    for (std::size_t k = 0; k < std::size(_heuristic_sol); ++k)
    {
        const auto h = std::size(_heuristic_sol[k]);
        for (std::size_t l = 0; l < h - 1; ++l)
        {
            const auto i = _heuristic_sol[k][l + 1];
            const auto j = _heuristic_sol[k][l];
            assert(!prob.conflict_matrix(i, j));
        }
    }
#endif
}

void MipCplex::build(const Problem &prob)
{
    _env.end();
    _env = IloEnv();
    _cplex = IloCplex(_env);
    _model = IloModel(_env);

    _cplex.setParam(IloCplex::Param::MIP::Display, MIPDisplay);
    _cplex.setParam(IloCplex::Param::MIP::Interval, MIPInterval);
    _cplex.setParam(IloCplex::Param::Threads, Threads);
    _cplex.setParam(IloCplex::Param::TimeLimit, TimeLimit);
    _cplex.setParam(IloCplex::Param::WorkMem, WorkMem);
    _cplex.setParam(IloCplex::Param::MIP::Strategy::File, MIPStrategyFile);
    _cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, MIPLimitsTreeMemory);

    TimePoint start_time = Clock::now();
    heuristic(prob);
    create(prob);
    _cplex.extract(_model);
    if constexpr (UseMipStart)
    {
        initialize(prob);
    }
    TimePoint finish_time = Clock::now();

    _reading_time = Duration(finish_time - start_time).count();
    _reading_memory = _env.getMemoryUsage();
    _reading_total_memory = _env.getTotalMemoryUsage();
}

void MipCplex::solve()
{
    constexpr double MB = 1024. * 1024.;
    std::vector<IncumbentCallback::Incumbent> incumbents;
    IloFastMutex mutex;

    if constexpr (UseIncumbentCallback)
    {
        _cplex.use(IloCplex::Callback(new (_env) IncumbentCallback(_env,
                                                                   incumbents,
                                                                   mutex)));
    }

    TimePoint start_time = Clock::now();
    _has_solution = _cplex.solve();
    TimePoint finish_time = Clock::now();

    _solving_time = Duration(finish_time - start_time).count();
    _solving_memory = _env.getMemoryUsage();
    _solving_total_memory = _env.getTotalMemoryUsage();

    if constexpr (UseIncumbentCallback)
    {
        _incumbents.clear();
        for (const auto &incumbent : incumbents)
        {
            double time = Duration(incumbent.time_point - start_time).count();
            _incumbents.push_back({time,
                                   incumbent.obj_value,
                                   incumbent.mip_relative_gap,
                                   incumbent.n_nodes,
                                   incumbent.n_remaining_nodes,
                                   incumbent.solution_source});
        }
    }

    LOG_INFO("status: %s", to_string(_cplex.getCplexStatus()));
    LOG_INFO("objective: %g", _cplex.getBestObjValue());
    LOG_INFO("reading time: %g s", _reading_time);
    LOG_INFO("solving time: %g s", _solving_time);
    LOG_INFO("total time: %g s", _reading_time + _solving_time);
    LOG_INFO("reading memory: %g / %g MB",
             _reading_memory / MB,
             _reading_total_memory / MB);
    LOG_INFO("solving memory: %g / %g MB",
             _solving_memory / MB,
             _solving_total_memory / MB);
}

bool MipCplex::has_solution() const
{
    return _has_solution;
}

bool MipCplex::check(const Problem &prob, const Solution &sol)
{
    const auto n_stacks = std::size(_heuristic_sol);

    // Check size of vectors
    if (std::size(sol.stack_of) != prob.n_items)
    {
        LOG_ERROR("stack_of has invalid size");
        return false;
    }

    if (std::size(sol.level_of) != prob.n_items)
    {
        LOG_ERROR("level_of has invalid size");
        return false;
    }

    // Check whether each item is assigned to a stack
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        if (sol.stack_of[i] >= prob.n_items)
        {
            LOG_ERROR("item %zu not assigned to a stack", i + 1);
            return false;
        }
    }

    // Check whether each item is assigned to a level
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        if (sol.level_of[i] >= prob.capacity)
        {
            LOG_ERROR("item %zu not assigned to a level", i + 1);
            return false;
        }
    }

    // Check initial items
    for (std::size_t i = 0; i < std::size(prob.initial_positions); ++i)
    {
        const auto k = prob.initial_positions[i];
        const auto k2 = sol.stack_of[i];
        if (k != k2)
        {
            LOG_ERROR("item %zu must be in stack %zu, not %zu",
                      i + 1,
                      k + 1,
                      k2 + 1);
            return false;
        }
    }

    // Check initial items
    for (std::size_t i = 0; i < std::size(prob.initial_positions); ++i)
    {
        for (std::size_t j = 0; j < std::size(prob.initial_positions); ++j)
        {
            if (i != j &&
                prob.initial_positions[i] == prob.initial_positions[j] &&
                sol.stack_of[i] != sol.stack_of[j])
            {
                LOG_ERROR("item %zu and %zu must be together", i + 1, j + 1);
                return false;
            }
        }
    }

    // Check arrival order
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            // Conflicts are adjacent
            if (sol.stack_of[i] == sol.stack_of[j] &&
                ((sol.level_of[i] > sol.level_of[j] &&
                  prob.arrivals[i] < prob.arrivals[j]) ||
                 (sol.level_of[j] > sol.level_of[i] &&
                  prob.arrivals[j] < prob.arrivals[i])))
            {
                LOG_ERROR("invalid order between items %zu and %zu",
                          i + 1,
                          j + 1);
                return false;
            }
        }
    }

    // Check conflict constraint violation
    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            // Conflicts are adjacent
            if (sol.stack_of[i] == sol.stack_of[j] &&
                ((sol.level_of[i] == sol.level_of[j] + 1 &&
                  prob.conflict_matrix(i, j)) ||
                 (sol.level_of[j] == sol.level_of[i] + 1 &&
                  prob.conflict_matrix(j, i))))
            {
                LOG_ERROR("conflict between items %zu and %zu",
                          i + 1,
                          j + 1);
                return false;
            }
        }
    }

    // Check capacity constraint
    for (std::size_t k = 0; k < n_stacks; ++k)
    {
        std::size_t count = 0;

        for (auto l : sol.stack_of)
        {
            if (k == l)
            {
                ++count;
            }
        }

        if (count > prob.capacity)
        {
            LOG_ERROR("stack %zu violates capacity", k + 1);
            return false;
        }
    }

    // Check objective value
    std::vector<bool> is_used(n_stacks, false);

    for (std::size_t i = 0; i < prob.n_items; ++i)
    {
        is_used[sol.stack_of[i]] = true;
    }

    auto n_used_stacks = std::count(std::begin(is_used),
                                    std::end(is_used),
                                    true);

    if (std::abs(double(n_used_stacks) - sol.obj) > 1e-6)
    {
        LOG_ERROR("objective values are different: %zu != %g",
                  n_used_stacks,
                  sol.obj);
        return false;
    }

    LOG_INFO("solution checked");
    return true;
}

void MipCplex::write_lp(const char *filename)
{
    _cplex.exportModel(filename);
}

void MipCplex::write([[maybe_unused]] const Problem &prob,
                     const Solution &sol,
                     nlohmann::json &json)
{
    json["solver"] = "cplex";

    if (!std::empty(sol.stack_of))
    {
        nlohmann::json stack_of = nlohmann::json::array();
        nlohmann::json level_of = nlohmann::json::array();

        for (auto k : sol.stack_of)
        {
            stack_of.push_back(k + 1);
        }

        for (auto l : sol.level_of)
        {
            level_of.push_back(l + 1);
        }

        json["sol"] = {
            {"obj", sol.obj},
            {"stack_of", stack_of},
            {"level_of", level_of}};
    }

    json["stats"] = {
        {"reading_time", _reading_time},
        {"solving_time", _solving_time},
        {"reading_memory", _reading_memory},
        {"solving_memory", _solving_memory},
        {"reading_total_memory", _reading_total_memory},
        {"solving_total_memory", _solving_total_memory},
        {"total_time", _reading_time + _solving_time},
        {"best_obj_value", _cplex.getBestObjValue()},
        {"cplex_status", to_string(_cplex.getCplexStatus())},
        {"cplex_status_n", _cplex.getCplexStatus()},
        {"cplex_sub_status", to_string(_cplex.getCplexSubStatus())},
        {"cplex_sub_status_n", _cplex.getCplexSubStatus()},
        {"cutoff", _cplex.getCutoff()},
        {"n_barrier_iterations", _cplex.getNbarrierIterations64()},
        {"n_bin_vars", _cplex.getNbinVars()},
        {"n_cols", _cplex.getNcols()},
        {"n_cross_d_exch", _cplex.getNcrossDExch64()},
        {"n_cross_d_push", _cplex.getNcrossDPush64()},
        {"n_cross_p_exch", _cplex.getNcrossPExch64()},
        {"n_cross_p_push", _cplex.getNcrossPPush64()},
        {"n_dual_superbasics", _cplex.getNdualSuperbasics()},
        {"n_filters", _cplex.getNfilters()},
        {"n_indicators", _cplex.getNindicators()},
        {"n_int_vars", _cplex.getNintVars()},
        {"n_iterations", _cplex.getNiterations64()},
        {"n_lcs", _cplex.getNLCs()},
        {"n_mip_starts", _cplex.getNMIPStarts()},
        {"n_nodes", _cplex.getNnodes64()},
        {"n_nodes_left", _cplex.getNnodesLeft64()},
        {"n_nzs", _cplex.getNNZs64()},
        {"n_phase_one_iterations", _cplex.getNphaseOneIterations64()},
        {"n_primal_superbasics", _cplex.getNprimalSuperbasics()},
        {"n_qcs", _cplex.getNQCs()},
        {"n_rows", _cplex.getNrows()},
        {"n_semi_cont_vars", _cplex.getNsemiContVars()},
        {"n_semi_int_vars", _cplex.getNsemiIntVars()},
        {"n_sifting_iterations", _cplex.getNsiftingIterations64()},
        {"n_sifting_phase_one_iterations", _cplex.getNsiftingPhaseOneIterations64()},
        {"n_soss", _cplex.getNSOSs()},
        {"n_ucs", _cplex.getNUCs()},
        {"num_cores", _cplex.getNumCores()},
        {"status", to_string(_cplex.getStatus())},
        {"status_n", _cplex.getStatus()},
        {"version_number", _cplex.getVersionNumber()},
        {"is_dual_feasible", _cplex.isDualFeasible()},
        {"is_mip", _cplex.isMIP()},
        {"is_primal_feasible", _cplex.isPrimalFeasible()},
        {"is_qc", _cplex.isQC()},
        {"is_qo", _cplex.isQO()}};

    if (_has_solution)
    {
        json["stats"].update({
            {"obj_value", _cplex.getObjValue()},
            {"mip_relative_gap", _cplex.getMIPRelativeGap()},
        });
    }

    if constexpr (UseIncumbentCallback)
    {
        json["incumbents"] = nlohmann::json::array();

        for (const auto &incumbent : _incumbents)
        {
            json["incumbents"].push_back({{"time", incumbent.time},
                                          {"obj", incumbent.obj},
                                          {"gap", incumbent.gap},
                                          {"n_nodes", incumbent.n_nodes},
                                          {"n_remaining_nodes", incumbent.n_remaining_nodes},
                                          {"source", to_string(incumbent.source)},
                                          {"source_n", incumbent.source}});
        }
    }
}
