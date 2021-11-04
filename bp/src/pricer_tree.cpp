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

#include "pricer_tree.hpp"
#include "graph_algorithm.hpp"
#include "logging.hpp"
#include "pprint.hpp"
#include "tree_search.hpp"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <memory>
#include <numeric>
#include <set>
#include <vector>

// Print the search tree
constexpr bool PrintTree = false;

// Print the solution
constexpr bool PrintSolution = false;

using namespace pslp::bp;

const char *PricerTree::Properties::Name = "tree";
const char *PricerTree::Properties::Desc = "tree search variable pricer";
const int PricerTree::Properties::Priority = 1000;
const SCIP_Bool PricerTree::Properties::Delay = true;

/**
 * Node type definition for distant conflicts.
 */
struct NodeDistant : tree_search::Node
{
    /**
     * Partial objective value.
     */
    double obj = 0.;

    /**
     * Partial volume.
     */
    std::size_t vol = 0;

    /**
     * Upper bound (relaxation value).
     */
    double upper_bound = 0.;

    /**
     * Partial solution.
     */
    std::vector<std::size_t> sol;

    /**
     * Ordered list of candidates for branching.
     */
    std::vector<std::size_t> cands;
};

/**
 * Determine if node1 is more promising than node2.
 */
bool operator<(const NodeDistant &node1, const NodeDistant &node2)
{
    return node1.upper_bound > node2.upper_bound;
}

struct AlgorithmDistant : tree_search::Algorithm<NodeDistant>
{
    /**
     * Node type definition.
     */
    using Node = NodeDistant;

    /**
     * Context type definition.
     */
    using Context = tree_search::Context<Node>;

    /**
     * Constructor.
     */
    AlgorithmDistant(SCIP *scip,
                     const Pricer::Problem &prob,
                     const std::vector<std::vector<std::size_t>> &unions,
                     const std::vector<std::size_t> &union_of,
                     const SymmetricMatrix<bool, true> &diff_matrix)
        : _scip(scip),
          _prob(prob),
          _unions(unions),
          _union_of(union_of),
          _diff_matrix(diff_matrix)
    {
        // Compute profits for each union
        _profits.resize(std::size(_unions), 0.);

        for (std::size_t u = 0; u < std::size(_unions); ++u)
        {
            for (auto i : unions[u])
            {
                _profits[u] += prob.coefs[i];
            }
        }
    }

    /**
     * Create the root node.
     */
    void start(Context &ctx) override
    {
        Node root;

        root.obj = _prob.offset;
        root.vol = 0;
        root.upper_bound = _prob.offset;
        root.cands.resize(std::size(_unions));

        // Sort candidates by profit / volume
        std::iota(std::begin(root.cands), std::end(root.cands), 0);
        std::sort(std::begin(root.cands),
                  std::end(root.cands),
                  [this](auto u, auto v) {
                      return _profits[u] / double(std::size(_unions[u])) >
                             _profits[v] / double(std::size(_unions[v]));
                  });

        // Remove non-positive candidates (they cannot improve the objective)
        {
            auto it = std::find_if_not(std::begin(root.cands),
                                       std::end(root.cands),
                                       [this](auto u) {
                                           return SCIPisPositive(_scip, _profits[u]);
                                       });
            root.cands.erase(it, std::end(root.cands));
        }

        // Remove candidates that exceed the capacity or having internal conflicts
        for (auto u = std::begin(root.cands); u != std::end(root.cands);)
        {
            if (_diff_matrix(*u, *u) ||
                std::size(_unions[*u]) > _prob.capacity)
            {
                u = root.cands.erase(u);
            }
            else
            {
                ++u;
            }
        }

        ctx.add_root(std::move(root));
    }

    /**
     * Called when the search is finished.
     */
    void finish(Context &ctx) override
    {
        _total_time = ctx.elapsed_time();
    }

    /**
     * Estimate a node before adding it to the queue.
     */
    bool estimate([[maybe_unused]] Context &ctx, Node &node) override
    {
        // Compute upper bound (relaxation)
        const auto n_unions = std::size(_unions);
        double obj = node.obj;
        std::size_t vol = node.vol;
        std::size_t critical_index = n_unions;

        // Candidates are sorted by ratio profit/volume
        for (auto u : node.cands)
        {
            if (vol + std::size(_unions[u]) <= _prob.capacity)
            {
                obj += _profits[u];
                vol += std::size(_unions[u]);
            }
            else
            {
                critical_index = u;
                break;
            }
        }

        // Add the fractional item for the upper bound
        if (critical_index < std::size(_unions))
        {
            obj += double(_prob.capacity - vol) *
                   double(_profits[critical_index]) /
                   double(std::size(_unions[critical_index]));
        }

        node.upper_bound = obj;

        // If the upper bound is not improving, reject the node
        if (!SCIPisPositive(_scip, node.upper_bound) ||
            SCIPisLE(_scip, node.upper_bound, _best_obj))
        {
            return false;
        }

        return true;
    }

    /**
     * Evaluate a node before branching.
     */
    bool evaluate([[maybe_unused]] Context &ctx, Node &node) override
    {
        assert(SCIPisPositive(_scip, node.upper_bound));

        std::vector<std::size_t> items;

        for (auto u : node.sol)
        {
            for (auto i : _unions[u])
            {
                items.push_back(i);
            }
        }

        // Get the subgraph built with the items
        auto subgraph = Pricer::subgraph(_prob.conflict_matrix, items);
        std::vector<std::size_t> s(std::size(items));

        // Detect feasibility by finding a topological sorting
        try
        {
            std::size_t pos = std::size(items);
            TopologicalSorter topo_sorter([&](auto v) {
                assert(pos > 0);
                s[--pos] = items[v];
            });
            depth_first_visit(subgraph, topo_sorter);
        }
        catch (const NotADag &)
        {
            // Infeasible
            return false;
        }

        // Is the node terminal (leaf)?
        if (std::empty(node.cands) ||                    // No candidate for branching
            node.vol >= _prob.capacity ||                // Capacity reached
            SCIPisGE(_scip, node.obj, node.upper_bound)) // LB >= UB
        {
            assert(check_solution(s, node.obj));

            // Check whether the node is improving and at least as good as the best
            if (SCIPisPositive(_scip, node.obj) &&
                SCIPisGE(_scip, node.obj, _best_obj))
            {
#ifndef NDEBUG
                // Check whether the solution is redundant
                {
                    auto sol = s;
                    std::sort(std::begin(sol), std::end(sol));
                    assert(std::adjacent_find(std::begin(sol), std::end(sol)) == std::end(sol));
                    for (const auto &s2 : _best_sols)
                    {
                        auto sol2 = s2;
                        std::sort(std::begin(sol2), std::end(sol2));
                        assert(sol != sol2);
                    }
                }
#endif

                if (_best_obj < node.obj)
                {
                    _best_obj = node.obj;
                }

                _best_sols.push_back(s);
            }

            return false;
        }

        return true;
    }

    /**
     * Create child nodes.
     */
    void branch(Context &ctx, Node &node) override
    {
        assert(!std::empty(node.cands));

        const auto u = node.cands.front();

        assert(u < std::size(_unions));
        assert(std::find(std::begin(node.sol), std::end(node.sol), u) == std::end(node.sol));

        node.cands.erase(std::begin(node.cands));

        // Include candidate
        if (node.vol + std::size(_unions[u]) <= _prob.capacity)
        {
            auto child = node;
            child.obj += _profits[u];
            child.vol += std::size(_unions[u]);
            child.sol.push_back(u);
            for (auto v = std::begin(child.cands); v != std::end(child.cands);)
            {
                if (_diff_matrix(u, *v) ||
                    node.vol + std::size(_unions[*v]) > _prob.capacity)
                {
                    v = child.cands.erase(v);
                }
                else
                {
                    ++v;
                }
            }
            ctx.add_child(std::move(child), node);
        }

        // Exclude candidate
        {
            auto child = node;
            ctx.add_child(std::move(child), node);
        }
    }

    /**
     * Check whether the algorithm should finish.
     */
    bool termination(ContextT &ctx) override
    {
        return _time_limit < 1.e20 && ctx.elapsed_time() >= _time_limit;
    }

    /**
     * Print a node.
     */
    void print(const Node &node, std::ostream &os = std::cout) const override
    {
        tree_search::Algorithm<Node>::print(node, os);
        os << std::endl;
        os << "obj:" << node.obj << " vol:" << node.vol << std::endl;
        os << "sol:";
        for (auto u : node.sol)
        {
            os << ' ' << u + 1;
        }
        os << std::endl;
        os << "cands:";
        for (auto u : node.cands)
        {
            os << ' ' << u + 1;
        }
        os << std::endl;
    }

    /**
     * Return the best solutions.
     */
    const auto &best_sols() const
    {
        return _best_sols;
    }

    /**
     * Return the total time.
     */
    double total_time() const
    {
        return _total_time;
    }

    /**
     * Set the time limit in seconds.
     */
    void set_time_limit(double time_limit)
    {
        _time_limit = time_limit;
    }

    /**
     * Check whether a solution satisfies all the constraints.
     */
    bool check_solution(const std::vector<std::size_t> &sol, double obj) const
    {
        const auto h = std::size(sol);

        // Check objective value
        double obj2 = _prob.offset;
        for (auto i : sol)
        {
            obj2 += _prob.coefs[i];
        }

        if (!SCIPisEQ(_scip, obj2, obj))
        {
            LOG_ERROR("wrong objective %g (expected %g)",
                      obj2,
                      obj);
            return false;
        }

        // Check capacity
        if (h > _prob.capacity)
        {
            LOG_ERROR("violated capacity %zu > %zu", h, _prob.capacity);
            return false;
        }

        // Check if conflict constraints are satisfied
        for (std::size_t l_i = 1; l_i < h; ++l_i)
        {
            const auto i = sol[l_i];
            for (std::size_t l_j = 0; l_j < l_i + 1; ++l_j)
            {
                const auto j = sol[l_j];
                if (_prob.conflict_matrix(i, j))
                {
                    LOG_ERROR("item %zu cannot be above %zu", i + 1, j + 1);
                    return false;
                }
            }
        }

        // Check if the same constraints are satisfied
        for (auto [i, j] : _prob.same)
        {
            bool has_i = std::find(std::begin(sol),
                                   std::end(sol),
                                   i) != std::end(sol);
            bool has_j = std::find(std::begin(sol),
                                   std::end(sol),
                                   j) != std::end(sol);

            if (has_i != has_j)
            {
                LOG_ERROR("items %zu and %zu must be together", i + 1, j + 1);
                return false;
            }
        }

        // Check if the diff constraints are satisfied
        for (auto [i, j] : _prob.diff)
        {
            bool has_i = std::find(std::begin(sol),
                                   std::end(sol),
                                   i) != std::end(sol);
            bool has_j = std::find(std::begin(sol),
                                   std::end(sol),
                                   j) != std::end(sol);

            if (has_i && has_j)
            {
                LOG_ERROR("items %zu and %zu must be separated", i + 1, j + 1);
                return false;
            }
        }

        return true;
    }

    SCIP *_scip;
    const Pricer::Problem &_prob;
    const std::vector<std::vector<std::size_t>> &_unions;
    const std::vector<std::size_t> &_union_of;
    const SymmetricMatrix<bool, true> &_diff_matrix;
    std::vector<double> _profits;
    std::vector<std::vector<std::size_t>> _best_sols;
    double _best_obj = 0.;
    double _total_time = 0.;
    double _time_limit = 1.e20;
};

PricerTree::PricerTree(SCIP *scip)
    : Pricer(scip,
             Properties::Name,
             Properties::Desc,
             Properties::Priority,
             Properties::Delay)
{
    {
        std::ostringstream oss;
        oss << "pricers/" << Properties::Name << "/nodeselection";
        SCIPaddStringParam(scip,
                           oss.str().c_str(),
                           "(pslp) node selection rule (bfs, brfs, cbfs, dfs)",
                           &_nodeselection,
                           false,
                           "dfs",
                           nullptr,
                           nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << Properties::Name << "/timelimit";
        SCIPaddRealParam(scip,
                         oss.str().c_str(),
                         "(pslp) maximal time in seconds to run",
                         &_timelimit,
                         false,
                         _timelimit,
                         0.,
                         1.e20,
                         nullptr,
                         nullptr);
    }
    {
        std::ostringstream oss;
        oss << "pricers/" << Properties::Name << "/extratime";
        SCIPaddBoolParam(scip,
                         oss.str().c_str(),
                         "(pslp) can the pricer exceed the global time limit?",
                         &_extratime,
                         false,
                         _extratime,
                         nullptr,
                         nullptr);
    }
}

SCIP_RETCODE PricerTree::init(SCIP *scip, const Problem &prob)
{
    SCIP_CALL(Pricer::init(scip, prob));

    // Get the time limit of the master problem
    SCIP_CALL(SCIPgetRealParam(scip, "limits/time", &_global_time_limit));

    // Get the node selection rule
    std::string nodeselection = _nodeselection;
    if (nodeselection == "bfs")
    {
        LOG_INFO("node selection rule: best first search");
        _queue_id = QueueID::BFS;
    }
    else if (nodeselection == "brfs")
    {
        LOG_INFO("node selection rule: breadth first search");
        _queue_id = QueueID::BrFS;
    }
    else if (nodeselection == "cbfs")
    {
        LOG_INFO("node selection rule: cyclic best first search");
        _queue_id = QueueID::CBFS;
    }
    else if (nodeselection == "dfs")
    {
        LOG_INFO("node selection rule: depth first search");
        _queue_id = QueueID::DFS;
    }
    else
    {
        LOG_ERROR("unknown node selection rule <%s>", _nodeselection);
        return SCIP_ERROR;
    }

    return SCIP_OKAY;
}

SCIP_RETCODE PricerTree::solve(SCIP *scip,
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

        if (!prob.has_adjacent_conflicts)
        {
            // If conflicts are distant, when i -> j and j -> i,
            // it is equivalent to a diff constraint
            for (auto i : prob.conflict_graph.vertices())
            {
                for (auto j : prob.conflict_graph.adjacent_vertices(i))
                {
                    const auto u = _union_of[i];
                    const auto v = _union_of[j];

                    if (!_diff_matrix(u, v) && prob.conflict_matrix(j, i))
                    {
                        _diff_matrix(u, v) = true;
                    }
                }
            }
        }
    }

    assert(!std::empty(_unions) && !std::empty(_union_of));

    if (prob.has_adjacent_conflicts)
    {
        throw std::runtime_error("adjacent conflicts not supported");
    }
    else ///////////////////////////////////////////////////////////////////////
    {
        AlgorithmDistant algorithm(scip, prob, _unions, _union_of, _diff_matrix);

        // Get the node selection rule
        std::unique_ptr<tree_search::Queue<NodeDistant>> queue;
        switch (_queue_id)
        {
            case QueueID::BFS:
            {
                using BFS = tree_search::BasicQueue<NodeDistant, tree_search::BFSCompare<>>;
                queue = std::make_unique<BFS>();
            }
            break;
            case QueueID::BrFS:
            {
                using BrFS = tree_search::BasicQueue<NodeDistant, tree_search::BrFSCompare<>>;
                queue = std::make_unique<BrFS>();
            }
            break;
            case QueueID::CBFS:
            {
                using CBFS = tree_search::CyclicBestFirstSearchQueue<NodeDistant>;
                queue = std::make_unique<CBFS>();
            }
            break;
            case QueueID::DFS:
            {
                using DFS = tree_search::BasicQueue<NodeDistant, tree_search::DFSCompare<>>;
                queue = std::make_unique<DFS>();
            }
            break;
        }

        // Pricer time limit
        algorithm.set_time_limit(_timelimit);

        // If _extratime = FALSE and master time limit < inf,
        // align the pricer time limit on the master
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

            if (remaining_time < _timelimit)
            {
                LOG_TRACE("set pricing time limit to %g s (elapsed: %g / %g s)",
                          remaining_time,
                          elapsed_time,
                          _global_time_limit);
                algorithm.set_time_limit(remaining_time);
            }
        }

        if constexpr (PrintTree)
        {
            tree_search::run(algorithm,
                             *queue,
                             tree_search::print_visitor<NodeDistant>(2));
        }
        else
        {
            tree_search::run(algorithm, *queue);
        }

#ifndef NDEBUG
        finish_time = Clock::now();

        double time = Duration(finish_time - start_time).count();

        LOG_DEBUG("tree pricing solving time: %g s", time);
#endif

        for (auto &sol : algorithm.best_sols())
        {
#ifndef NDEBUG
            double obj = Pricer::obj_value(prob, sol);
            LOG_DEBUG("found tree pricing sol (obj: %g)", obj);
            assert(SCIPisPositive(scip, obj));
#endif

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
    }

    return SCIP_OKAY;
}