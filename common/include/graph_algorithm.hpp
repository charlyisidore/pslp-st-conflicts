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

#ifndef GRAPH_ALGORITHM_HPP
#define GRAPH_ALGORITHM_HPP

#include "graph.hpp"
#include <deque>
#include <limits>
#include <optional>
#include <stdexcept>

/**
 * Inform that the graph is not a directed acyclic graph.
 */
struct NotADag : std::runtime_error
{
    NotADag()
        : std::runtime_error("not a directed acyclic graph") {}
};

/**
 * Wrapper for Breadth-First Search visitors.
 */
template <typename VisitorT>
class _BFSVisitor
{
public:
    enum class Mark
    {
        Unvisited = 0,
        Visiting = 1,
        Visited = 2
    };

    _BFSVisitor(VisitorT &visitor)
        : _visitor(visitor) {}

    template <typename GraphT>
    void discover_vertex(Vertex v, const GraphT &g)
    {
        _discover_vertex(_visitor, v, g);
    }

    template <typename GraphT>
    void examine_vertex(Vertex v, const GraphT &g)
    {
        _examine_vertex(_visitor, v, g);
    }

    template <typename GraphT>
    void examine_edge(Edge e, const GraphT &g)
    {
        _examine_edge(_visitor, e, g);
    }

    template <typename GraphT>
    void tree_edge(Edge e, const GraphT &g)
    {
        _tree_edge(_visitor, e, g);
    }

    template <typename GraphT>
    void non_tree_edge(Edge e, const GraphT &g)
    {
        _non_tree_edge(_visitor, e, g);
    }

    template <typename GraphT>
    void visiting_target(Edge e, const GraphT &g)
    {
        _visiting_target(_visitor, e, g);
    }

    template <typename GraphT>
    void visited_target(Edge e, const GraphT &g)
    {
        _visited_target(_visitor, e, g);
    }

    template <typename GraphT>
    void finish_vertex(Vertex v, const GraphT &g)
    {
        _finish_vertex(_visitor, v, g);
    }

    template <typename GraphT>
    bool terminate(Vertex v, const GraphT &g)
    {
        return _terminate(_visitor, v, g);
    }

private:
    template <typename T, typename GraphT>
    auto _discover_vertex(T &vis, Vertex v, const GraphT &g)
        -> decltype(vis.discover_vertex(v, g), void())
    {
        vis.discover_vertex(v, g);
    }

    template <typename T, typename GraphT>
    auto _examine_vertex(T &vis, Vertex v, const GraphT &g)
        -> decltype(vis.examine_vertex(v, g), void())
    {
        vis.examine_vertex(v, g);
    }

    template <typename T, typename GraphT>
    auto _examine_edge(T &vis, Edge e, const GraphT &g)
        -> decltype(vis.examine_edge(e, g), void())
    {
        vis.examine_edge(e, g);
    }

    template <typename T, typename GraphT>
    auto _tree_edge(T &vis, Edge e, const GraphT &g)
        -> decltype(vis.tree_edge(e, g), void())
    {
        vis.tree_edge(e, g);
    }

    template <typename T, typename GraphT>
    auto _non_tree_edge(T &vis, Edge e, const GraphT &g)
        -> decltype(vis.non_tree_edge(e, g), void())
    {
        vis.non_tree_edge(e, g);
    }

    template <typename T, typename GraphT>
    auto _visiting_target(T &vis, Edge e, const GraphT &g)
        -> decltype(vis.visiting_target(e, g), void())
    {
        vis.visiting_target(e, g);
    }

    template <typename T, typename GraphT>
    auto _visited_target(T &vis, Edge e, const GraphT &g)
        -> decltype(vis.visited_target(e, g), void())
    {
        vis.visited_target(e, g);
    }

    template <typename T, typename GraphT>
    auto _finish_vertex(T &vis, Vertex v, const GraphT &g)
        -> decltype(vis.finish_vertex(v, g), void())
    {
        vis.finish_vertex(v, g);
    }

    template <typename T, typename GraphT>
    auto _terminate(T &vis, Vertex v, const GraphT &g)
        -> decltype(vis.terminate(v, g), bool())
    {
        return vis.terminate(v, g);
    }

    void _discover_vertex(...) {}
    void _examine_vertex(...) {}
    void _examine_edge(...) {}
    void _tree_edge(...) {}
    void _non_tree_edge(...) {}
    void _visiting_target(...) {}
    void _visited_target(...) {}
    void _finish_vertex(...) {}
    bool _terminate(...)
    {
        return false;
    }

    VisitorT &_visitor;
};

/**
 * Wrapper for Depth-First Search visitors.
 */
template <typename VisitorT>
class _DFSVisitor
{
public:
    enum class Mark
    {
        Unvisited = 0,
        Visiting = 1,
        Visited = 2
    };

    _DFSVisitor(VisitorT &visitor)
        : _visitor(visitor) {}

    template <typename GraphT>
    void start_vertex(Vertex v, const GraphT &g)
    {
        _start_vertex(_visitor, v, g);
    }

    template <typename GraphT>
    void discover_vertex(Vertex v, const GraphT &g)
    {
        _discover_vertex(_visitor, v, g);
    }

    template <typename GraphT>
    void examine_edge(Edge e, const GraphT &g)
    {
        _examine_edge(_visitor, e, g);
    }

    template <typename GraphT>
    void tree_edge(Edge e, const GraphT &g)
    {
        _tree_edge(_visitor, e, g);
    }

    template <typename GraphT>
    void back_edge(Edge e, const GraphT &g)
    {
        _back_edge(_visitor, e, g);
    }

    template <typename GraphT>
    void forward_or_cross_edge(Edge e, const GraphT &g)
    {
        _forward_or_cross_edge(_visitor, e, g);
    }

    template <typename GraphT>
    void finish_edge(Edge e, const GraphT &g)
    {
        _finish_edge(_visitor, e, g);
    }

    template <typename GraphT>
    void finish_vertex(Vertex v, const GraphT &g)
    {
        _finish_vertex(_visitor, v, g);
    }

    template <typename GraphT>
    bool terminate(Vertex v, const GraphT &g)
    {
        return _terminate(_visitor, v, g);
    }

private:
    template <typename T, typename GraphT>
    auto _start_vertex(T &vis, Vertex v, const GraphT &g)
        -> decltype(vis.start_vertex(v, g), void())
    {
        vis.start_vertex(v, g);
    }

    template <typename T, typename GraphT>
    auto _discover_vertex(T &vis, Vertex v, const GraphT &g)
        -> decltype(vis.discover_vertex(v, g), void())
    {
        vis.discover_vertex(v, g);
    }

    template <typename T, typename GraphT>
    auto _examine_edge(T &vis, Edge e, const GraphT &g)
        -> decltype(vis.examine_edge(e, g), void())
    {
        vis.examine_edge(e, g);
    }

    template <typename T, typename GraphT>
    auto _tree_edge(T &vis, Edge e, const GraphT &g)
        -> decltype(vis.tree_edge(e, g), void())
    {
        vis.tree_edge(e, g);
    }

    template <typename T, typename GraphT>
    auto _back_edge(T &vis, Edge e, const GraphT &g)
        -> decltype(vis.back_edge(e, g), void())
    {
        vis.back_edge(e, g);
    }

    template <typename T, typename GraphT>
    auto _forward_or_cross_edge(T &vis, Edge e, const GraphT &g)
        -> decltype(vis.forward_or_cross_edge(e, g), void())
    {
        vis.forward_or_cross_edge(e, g);
    }

    template <typename T, typename GraphT>
    auto _finish_edge(T &vis, Edge e, const GraphT &g)
        -> decltype(vis.finish_edge(e, g), void())
    {
        vis.finish_edge(e, g);
    }

    template <typename T, typename GraphT>
    auto _finish_vertex(T &vis, Vertex v, const GraphT &g)
        -> decltype(vis.finish_vertex(v, g), void())
    {
        vis.finish_vertex(v, g);
    }

    template <typename T, typename GraphT>
    auto _terminate(T &vis, Vertex v, const GraphT &g)
        -> decltype(vis.terminate(v, g), bool())
    {
        return vis.terminate(v, g);
    }

    void _start_vertex(...) {}
    void _discover_vertex(...) {}
    void _examine_edge(...) {}
    void _tree_edge(...) {}
    void _back_edge(...) {}
    void _forward_or_cross_edge(...) {}
    void _finish_edge(...) {}
    void _finish_vertex(...) {}
    bool _terminate(...)
    {
        return false;
    }

    VisitorT &_visitor;
};

/**
 * Explore a graph using a breadth-first search.
 */
template <typename GraphT, typename VisitorT>
void breadth_first_visit(const GraphT &g, VisitorT &visitor)
{
    using Mark = typename _DFSVisitor<VisitorT>::Mark;

    _BFSVisitor<VisitorT> vis(visitor);
    std::vector<Mark> marks(g.n_vertices(), Mark::Unvisited);

    for (auto v : g.vertices())
    {
        if (marks[v] != Mark::Unvisited)
        {
            continue;
        }

        _breadth_first_visit(g, vis, v, marks);
    }
}

/**
 * Explore a graph using a breadth-first search.
 */
template <typename GraphT, typename VisitorT>
void breadth_first_visit(const GraphT &g, VisitorT &visitor, Vertex v)
{
    using Mark = typename _DFSVisitor<VisitorT>::Mark;

    _BFSVisitor<VisitorT> vis(visitor);
    std::vector<Mark> marks(g.n_vertices(), Mark::Unvisited);

    _breadth_first_visit(g, vis, v, marks);
}

/**
 * Explore a graph using a breadth-first search.
 */
template <typename GraphT, typename VisitorT>
void _breadth_first_visit(const GraphT &g,
                          _BFSVisitor<VisitorT> &vis,
                          Vertex v,
                          std::vector<typename _BFSVisitor<VisitorT>::Mark> &marks)
{
    using Mark = typename _DFSVisitor<VisitorT>::Mark;

    std::deque<Vertex> queue;

    marks[v] = Mark::Visiting;
    vis.discover_vertex(v, g);
    queue.push_back(v);

    while (!std::empty(queue))
    {
        Vertex u = queue.front();
        queue.pop_front();
        vis.examine_vertex(u, g);

        for (auto a : g.adjacent_vertices(u))
        {
            Edge e = {u, a};
            vis.examine_edge(e, g);

            if (marks[a] == Mark::Unvisited)
            {
                vis.tree_edge(e, g);
                marks[a] = Mark::Visiting;
                vis.discover_vertex(a, g);
                queue.push_back(a);
            }
            else
            {
                vis.non_tree_edge(e, g);
                if (marks[a] == Mark::Visiting)
                {
                    vis.visiting_target(e, g);
                }
                else
                {
                    assert(marks[a] == Mark::Visited);
                    vis.visited_target(e, g);
                }
            }
        }

        marks[u] = Mark::Visited;
        vis.finish_vertex(u, g);
    }
}

/**
 * Explore a graph using a depth-first search.
 */
template <typename GraphT, typename VisitorT>
void depth_first_visit(const GraphT &g, VisitorT &visitor)
{
    using Mark = typename _DFSVisitor<VisitorT>::Mark;

    _DFSVisitor<VisitorT> vis(visitor);
    std::vector<Mark> marks(g.n_vertices(), Mark::Unvisited);

    for (auto v : g.vertices())
    {
        if (marks[v] != Mark::Unvisited)
        {
            continue;
        }

        _depth_first_visit(g, vis, v, marks);
    }
}

/**
 * Explore a using a depth-first search from a single vertex.
 */
template <typename GraphT, typename VisitorT>
void depth_first_visit(const GraphT &g, VisitorT &visitor, Vertex v)
{
    using Mark = typename _DFSVisitor<VisitorT>::Mark;

    _DFSVisitor<VisitorT> vis(visitor);
    std::vector<Mark> marks(g.n_vertices(), Mark::Unvisited);

    _depth_first_visit(g, vis, v, marks);
}

/**
 * Explore a graph using a depth-first search.
 */
template <typename GraphT, typename VisitorT>
void _depth_first_visit(const GraphT &g,
                        _DFSVisitor<VisitorT> &vis,
                        Vertex v,
                        std::vector<typename _DFSVisitor<VisitorT>::Mark> &marks)
{
    using Mark = typename _DFSVisitor<VisitorT>::Mark;
    using OutEdgeIterator = typename GraphT::OutEdgeIterator;
    using State = std::tuple<Vertex,
                             std::optional<Edge>,
                             OutEdgeIterator,
                             OutEdgeIterator>;

    std::vector<State> stack;
    std::optional<Edge> e_source;
    OutEdgeIterator e;
    OutEdgeIterator e_end;

    vis.start_vertex(v, g);
    marks[v] = Mark::Visiting;
    vis.discover_vertex(v, g);

    std::tie(e, e_end) = g.out_edges(v);

    if (vis.terminate(v, g))
    {
        stack.push_back({v, {}, e_end, e_end});
    }
    else
    {
        stack.push_back({v, {}, e, e_end});
    }

    while (!std::empty(stack))
    {
        std::tie(v, e_source, e, e_end) = std::move(stack.back());
        stack.pop_back();

        if (e_source)
        {
            vis.finish_edge(*e_source, g);
        }

        while (e != e_end)
        {
            auto u = (*e).second;

            vis.examine_edge(*e, g);

            if (marks[u] == Mark::Unvisited)
            {
                vis.tree_edge(*e, g);

                e_source = *e;
                stack.push_back({v, e_source, ++e, e_end});

                v = u;
                marks[v] = Mark::Visiting;
                vis.discover_vertex(u, g);

                std::tie(e, e_end) = g.out_edges(v);

                if (vis.terminate(v, g))
                {
                    e = e_end;
                }
            }
            else
            {
                if (marks[u] == Mark::Visiting)
                {
                    vis.back_edge(*e, g);
                }
                else
                {
                    assert(marks[u] == Mark::Visited);
                    vis.forward_or_cross_edge(*e, g);
                }

                vis.finish_edge(*e, g);
                ++e;
            }
        }

        marks[v] = Mark::Visited;
        vis.finish_vertex(v, g);
    }
}

/**
 * Prints events of a Breadth-First Search.
 */
class BFSPrinter
{
public:
    BFSPrinter(std::ostream &os = std::cout)
        : _os(os) {}

    template <typename GraphT>
    void discover_vertex(Vertex v, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "discover vertex " << v + 1 << std::endl;
    }

    template <typename GraphT>
    void examine_vertex(Vertex v, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth++, ' ')
            << "examine vertex " << v + 1 << std::endl;
    }

    template <typename GraphT>
    void examine_edge(Edge e, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "examine edge " << e.first + 1
            << " -> " << e.second + 1 << std::endl;
    }

    template <typename GraphT>
    void tree_edge(Edge e, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "tree edge " << e.first + 1
            << " -> " << e.second + 1 << std::endl;
    }

    template <typename GraphT>
    void non_tree_edge(Edge e, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "non tree edge " << e.first + 1
            << " -> " << e.second + 1 << std::endl;
    }

    template <typename GraphT>
    void visiting_target(Edge e, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "visiting target " << e.first + 1
            << " -> " << e.second + 1 << std::endl;
    }

    template <typename GraphT>
    void visited_target(Edge e, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "visited target " << e.first + 1
            << " -> " << e.second + 1 << std::endl;
    }

    template <typename GraphT>
    void finish_vertex(Vertex v, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * --_depth, ' ')
            << "finish vertex " << v + 1 << std::endl;
    }

private:
    std::ostream &_os;
    std::size_t _depth = 0;
};

/**
 * Prints events of a Depth-First Search.
 */
class DFSPrinter
{
public:
    DFSPrinter(std::ostream &os = std::cout)
        : _os(os) {}

    template <typename GraphT>
    void start_vertex(Vertex v, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "start vertex " << v + 1 << std::endl;
    }

    template <typename GraphT>
    void discover_vertex(Vertex v, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth++, ' ')
            << "discover vertex " << v + 1 << std::endl;
    }

    template <typename GraphT>
    void examine_edge(Edge e, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "examine edge " << e.first + 1
            << " -> " << e.second + 1 << std::endl;
    }

    template <typename GraphT>
    void tree_edge(Edge e, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "tree edge " << e.first + 1
            << " -> " << e.second + 1 << std::endl;
    }

    template <typename GraphT>
    void back_edge(Edge e, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "back edge " << e.first + 1
            << " -> " << e.second + 1 << std::endl;
    }

    template <typename GraphT>
    void forward_or_cross_edge(Edge e, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "forward or cross edge " << e.first + 1
            << " -> " << e.second + 1 << std::endl;
    }

    template <typename GraphT>
    void finish_edge(Edge e, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * _depth, ' ')
            << "finish edge " << e.first + 1
            << " -> " << e.second + 1 << std::endl;
    }

    template <typename GraphT>
    void finish_vertex(Vertex v, [[maybe_unused]] const GraphT &g)
    {
        _os << std::string(2 * --_depth, ' ')
            << "finish vertex " << v + 1 << std::endl;
    }

private:
    std::ostream &_os;
    std::size_t _depth = 0;
};

/**
 * Find the shortest paths using Dijkstra.
 */
template <typename WeightF, typename DistanceF, typename PutF>
class DijkstraVisitor
{
public:
    DijkstraVisitor(WeightF weight,
                    DistanceF distance,
                    PutF put)
        : _weight(weight),
          _distance(distance),
          _put(put) {}

    template <typename GraphT>
    void tree_edge(Edge e, [[maybe_unused]] const GraphT &g)
    {
        relax(e);
    }

    template <typename GraphT>
    void visiting_target(Edge e, [[maybe_unused]] const GraphT &g)
    {
        relax(e);
    }

    void relax(Edge e)
    {
        Vertex u = e.first;
        Vertex v = e.second;
        const auto d = _distance(u) + _weight(e);
        if (d < _distance(v))
        {
            _put(v, d, u);
        }
    }

private:
    WeightF _weight;
    DistanceF _distance;
    PutF _put;
};

/**
 * Find the shortest paths using Dijkstra.
 */
template <typename GraphT, typename WeightF, typename DistanceF, typename PutF>
inline void dijkstra_shortest_paths(const GraphT &g,
                                    Vertex v,
                                    WeightF weight,
                                    DistanceF distance,
                                    PutF put)
{
    using DistanceT = decltype(distance(0));
    DijkstraVisitor vis(weight, distance, put);
    DistanceT inf = std::numeric_limits<DistanceT>::max();
    for (auto a : g.vertices())
    {
        put(a, (a == v) ? 0 : inf, a);
    }
    g.breadth_first_visit(vis, v);
}

/**
 * Find the shortest paths using Dijkstra.
 */
template <typename GraphT, typename WeightT, typename DistanceT>
inline void dijkstra_shortest_paths(const GraphT &g,
                                    Vertex v,
                                    const std::vector<std::vector<WeightT>> &weight,
                                    std::vector<DistanceT> &distance,
                                    std::vector<Vertex> &predecessor)
{
    const auto weight_f = [&weight](Edge e) {
        assert(e.first < std::size(weight));
        assert(e.second < std::size(weight[e.first]));
        return weight[e.first][e.second];
    };

    const auto distance_f = [&distance](Vertex v) {
        assert(v < std::size(distance));
        return distance[v];
    };

    auto put_f = [&predecessor, &distance](Vertex v, DistanceT d, Vertex u) {
        assert(v < std::size(predecessor));
        assert(v < std::size(distance));
        predecessor[v] = u;
        distance[v] = d;
    };

    dijkstra_shortest_paths(g, v, weight_f, distance_f, put_f);
}

/**
 * Topological sort using Depth-First Search.
 */
template <typename CallbackT>
struct TopologicalSorter
{
    TopologicalSorter(CallbackT callback)
        : callback(callback) {}

    template <typename GraphT>
    void back_edge([[maybe_unused]] Edge e,
                   [[maybe_unused]] const GraphT &g)
    {
        throw NotADag();
    }

    template <typename GraphT>
    void finish_vertex(Vertex v, [[maybe_unused]] const GraphT &g)
    {
        callback(v);
    }

    CallbackT callback;
};

/**
 * Find cycles using Depth-First Search.
 */
template <typename CallbackT>
struct CycleFinder
{
    CycleFinder(CallbackT callback)
        : callback(callback) {}

    template <typename GraphT>
    void discover_vertex(Vertex v, [[maybe_unused]] const GraphT &g)
    {
        assert(std::find(std::begin(stack),
                         std::end(stack),
                         v) == std::end(stack));
        stack.push_back(v);
    }

    template <typename GraphT>
    void back_edge(Edge e, [[maybe_unused]] const GraphT &g)
    {
        callback(std::find(std::begin(stack), std::end(stack), e.second),
                 std::end(stack));
    }

    template <typename GraphT>
    void finish_vertex([[maybe_unused]] Vertex v,
                       [[maybe_unused]] const GraphT &g)
    {
        assert(v == stack.back());
        stack.pop_back();
    }

    std::vector<Vertex> stack;
    CallbackT callback;
};

#endif