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

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

/**
 * Vertex type definition.
 */
using Vertex = std::size_t;

/**
 * Edge type definition.
 */
using Edge = std::pair<Vertex, Vertex>;

/**
 * Graph based on an adjacency list structure.
 */
template <bool DirectedV>
class AdjacencyList
{
private:
    class _VertexIterator;
    class _EdgeIterator;
    class _OutEdgeIterator;
    class _AdjacencyIterator;
    template <typename IteratorT>
    struct Range;
    template <typename VisitorT>
    struct _DFSVisitor;

public:
    /**
     * Vertex iterator type definition.
     */
    using VertexIterator = _VertexIterator;

    /**
     * Edge iterator type definition.
     */
    using EdgeIterator = _EdgeIterator;

    /**
     * Outgoing edge iterator type definition.
     */
    using OutEdgeIterator = _OutEdgeIterator;

    /**
     * Adjacency iterator type definition.
     */
    using AdjacencyIterator = std::vector<Vertex>::const_iterator;

    /**
     * Default constructor.
     */
    AdjacencyList() = default;

    /**
     * Create a graph with n vertices.
     */
    AdjacencyList(std::size_t n)
        : _adj_list(n) {}

    /**
     * Copy constructor.
     */
    AdjacencyList(const AdjacencyList &g) = default;

    /**
     * Move constructor.
     */
    AdjacencyList(AdjacencyList &&g) = default;

    /**
     * Copy assignment operator.
     */
    AdjacencyList &operator=(const AdjacencyList &g) = default;

    /**
     * Move assignment operator.
     */
    AdjacencyList &operator=(AdjacencyList &&g) = default;

    /**
     * Check whether the graph is directed.
     */
    bool is_directed() const
    {
        return DirectedV;
    }

    /**
     * Get the number of vertices.
     */
    std::size_t n_vertices() const
    {
        return std::size(_adj_list);
    }

    /**
     * Get the number of edges.
     */
    std::size_t n_edges() const
    {
        std::size_t m = 0;
        for (const auto &adj : _adj_list)
        {
            m += std::size(adj);
        }
        if constexpr (DirectedV)
        {
            return m;
        }
        else
        {
            assert(m % 2 == 0);
            return m;
        }
    }

    /**
     * Get the number of edges connecting a vertex.
     */
    std::size_t degree(Vertex v) const
    {
        assert(v < std::size(_adj_list));
        return std::size(_adj_list[v]);
    }

    /**
     * Check whether an edge between u and v exists.
     */
    bool has_edge(Vertex u, Vertex v) const
    {
        assert(u < std::size(_adj_list));
        assert(v < std::size(_adj_list));
        assert(u != v);
        return std::binary_search(std::begin(_adj_list[u]),
                                  std::end(_adj_list[u]),
                                  v);
    }

    /**
     * Get the range of all vertices.
     */
    Range<VertexIterator> vertices() const
    {
        return {VertexIterator(0), VertexIterator(std::size(_adj_list))};
    }

    /**
     * Returns an iterator range to access all edges.
     */
    Range<EdgeIterator> edges() const
    {
        return {EdgeIterator(_adj_list, Vertex(0)),
                EdgeIterator(_adj_list, Vertex(std::size(_adj_list)))};
    }

    /**
     * Returns an iterator range to access outgoing edges of a vertex.
     */
    Range<OutEdgeIterator> out_edges(Vertex v) const
    {
        assert(v < std::size(_adj_list));
        return {OutEdgeIterator(v, std::begin(_adj_list[v])),
                OutEdgeIterator(v, std::end(_adj_list[v]))};
    }

    /**
     * Returns an iterator range to access vertices adjacent to vertex v.
     */
    Range<AdjacencyIterator> adjacent_vertices(Vertex v) const
    {
        assert(v < std::size(_adj_list));
        return {std::begin(_adj_list[v]), std::end(_adj_list[v])};
    }

    /**
     * Initialize an empty graph with n vertices.
     */
    void initialize(std::size_t n)
    {
        _adj_list.clear();
        _adj_list.resize(n);
    }

    /**
     * Add a vertex.
     */
    Vertex add_vertex()
    {
        _adj_list.emplace_back();
        return Vertex(std::size(_adj_list) - 1);
    }

    /**
     * Add new vertices.
     */
    void add_vertices(std::size_t n)
    {
        _adj_list.resize(std::size(_adj_list) + n);
    }

    /**
     * Add an edge between vertices u and v.
     */
    bool add_edge(Vertex u, Vertex v)
    {
        assert(u < std::size(_adj_list));
        assert(v < std::size(_adj_list));
        assert(u != v);

        auto it_v = std::lower_bound(std::begin(_adj_list[u]),
                                     std::end(_adj_list[u]),
                                     v);

        // Does the edge exist?
        if (it_v != std::end(_adj_list[u]) && *it_v == v)
        {
            if constexpr (!DirectedV)
            {
                // Undirected graphs: if u -> v exists, v -> u exists
                assert(std::find(std::begin(_adj_list[v]),
                                 std::end(_adj_list[v]),
                                 u) != std::end(_adj_list[v]));
            }
            return false;
        }

        assert(it_v == std::end(_adj_list[u]) || *it_v > v);

        _adj_list[u].insert(it_v, v);

        // The adjacency list of u must remain sorted
        assert(std::is_sorted(std::begin(_adj_list[u]),
                              std::end(_adj_list[u])));
        assert(std::adjacent_find(std::begin(_adj_list[u]),
                                  std::end(_adj_list[u])) ==
               std::end(_adj_list[u]));

        if constexpr (!DirectedV)
        {
            // Undirected graphs: v -> u must be created as well
            auto it_u = std::lower_bound(std::begin(_adj_list[v]),
                                         std::end(_adj_list[v]),
                                         u);
            assert(it_u == std::end(_adj_list[v]) || *it_u > u);

            _adj_list[v].insert(it_u, u);

            // The adjacency list of v must remain sorted
            assert(std::is_sorted(std::begin(_adj_list[v]),
                                  std::end(_adj_list[v])));
            assert(std::adjacent_find(std::begin(_adj_list[v]),
                                      std::end(_adj_list[v])) ==
                   std::end(_adj_list[v]));
        }

        return true;
    }

    /**
     * Remove a given vertex.
     */
    void remove_vertex(Vertex v)
    {
        assert(v < std::size(_adj_list));

        _adj_list.erase(std::begin(_adj_list) + v);

        for (auto &adj : _adj_list)
        {
            auto u = std::lower_bound(std::begin(adj), std::end(adj), v);
            while (u != std::end(adj))
            {
                if (*u > v)
                {
                    --(*u);
                    ++u;
                }
                else
                {
                    assert(*u == v);
                    u = adj.erase(u);
                }
            }
        }
    }

    /**
     * Remove the edge between u and v.
     */
    bool remove_edge(Vertex u, Vertex v)
    {
        assert(u < std::size(_adj_list));
        assert(v < std::size(_adj_list));
        assert(u != v);
        auto a = std::lower_bound(std::begin(_adj_list[u]),
                                  std::end(_adj_list[u]),
                                  v);
        bool remove = (a != std::end(_adj_list[u]) && *a == v);
        if (remove)
        {
            _adj_list[u].erase(a);
            if constexpr (!DirectedV)
            {
                auto b = std::lower_bound(std::begin(_adj_list[v]),
                                          std::end(_adj_list[v]),
                                          u);
                assert(b != std::end(_adj_list[v]) && *b == u);
                _adj_list[v].erase(b);
            }
        }
        return remove;
    }

    /**
     * Remove all the vertices and edges.
     */
    void clear()
    {
        _adj_list.clear();
    }

    /**
     * Remove all the edges.
     */
    void clear_edges()
    {
        for (auto &adj : _adj_list)
        {
            adj.clear();
        }
    }

    /**
     * Explore the graph using a depth-first search.
     */
    template <typename VisitorT>
    void depth_first_visit(VisitorT &&visitor) const
    {
        depth_first_visit(visitor);
    }

    /**
     * Explore the graph using a depth-first search.
     */
    template <typename VisitorT>
    void depth_first_visit(VisitorT &visitor) const
    {
        enum Mark
        {
            Unvisited = 0,
            Visiting = 1,
            Visited = 2
        };

        using State = std::tuple<Vertex,
                                 std::optional<Edge>,
                                 OutEdgeIterator,
                                 OutEdgeIterator>;

        _DFSVisitor<VisitorT> vis(visitor);
        std::vector<Mark> marks(n_vertices(), Unvisited);
        std::vector<State> stack;
        std::optional<Edge> e_source;
        OutEdgeIterator e;
        OutEdgeIterator e_end;

        for (auto v : vertices())
        {
            if (marks[v] != Unvisited)
            {
                continue;
            }

            vis.start_vertex(v, *this);

            marks[v] = Visiting;
            vis.discover_vertex(v, *this);

            std::tie(e, e_end) = out_edges(v);

            if (vis.terminate(v, *this))
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
                    vis.finish_edge(*e_source, *this);
                }

                while (e != e_end)
                {
                    auto u = (*e).second;

                    vis.examine_edge(*e, *this);

                    if (marks[u] == Unvisited)
                    {
                        vis.tree_edge(*e, *this);

                        e_source = *e;
                        stack.push_back({v, e_source, ++e, e_end});

                        v = u;
                        marks[v] = Visiting;
                        vis.discover_vertex(u, *this);

                        std::tie(e, e_end) = out_edges(v);

                        if (vis.terminate(v, *this))
                        {
                            e = e_end;
                        }
                    }
                    else
                    {
                        if (marks[u] == Visiting)
                        {
                            vis.back_edge(*e, *this);
                        }
                        else // marks[u] == Visited
                        {
                            vis.forward_or_cross_edge(*e, *this);
                        }

                        vis.finish_edge(*e, *this);
                        ++e;
                    }
                }

                marks[v] = Visited;
                vis.finish_vertex(v, *this);
            }
        }
    }

    /**
     * Contract vertices u and v into a single vertex (removes v).
     */
    Vertex contract_vertices(Vertex u, Vertex v)
    {
        assert(u < std::size(_adj_list));
        assert(v < std::size(_adj_list));
        assert(u != v);
        for (auto [a, b] : edges())
        {
            if (a == v)
            {
                add_edge(u, b);
            }
            else if (b == v)
            {
                add_edge(a, u);
            }
        }
        remove_vertex(v);
        return Vertex(u < v ? u : u - 1);
    }

    /**
     * Get a subgraph of the graph.
     */
    template <typename IteratorT>
    AdjacencyList subgraph(const std::vector<Vertex> &s) const
    {
        std::size_t n = std::size(s);
        AdjacencyList g(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            Vertex u = s[i];
            for (std::size_t j = 0; j < n; ++j)
            {
                if (i == j)
                {
                    continue;
                }
                Vertex v = s[j];
                if (has_edge(u, v))
                {
                    g.add_edge(i, j);
                }
            }
        }
        return g;
    }

    /**
     * Get the complement of the graph.
     */
    AdjacencyList complement() const
    {
        const auto n = std::size(_adj_list);
        AdjacencyList g(n);
        for (Vertex u = 0; u < n; ++u)
        {
            Vertex v = 0;
            g._adj_list[u].clear();
            g._adj_list[u].reserve(n - std::size(_adj_list[u]));
            for (auto w : _adj_list[u])
            {
                while (v < w)
                {
                    if (u != v)
                    {
                        g._adj_list[u].push_back(v);
                    }
                    ++v;
                }
                assert(v == w);
                ++v;
            }
            while (v < n)
            {
                if (u != v)
                {
                    g._adj_list[u].push_back(v);
                }
                ++v;
            }
            assert(std::size(g._adj_list[u]) + std::size(_adj_list[u]) + 1 == n);
        }
        return g;
    }

    /**
     * Print the graph.
     */
    void print(std::ostream &os = std::cout) const
    {
        for (Vertex v = 0; v < std::size(_adj_list); ++v)
        {
            os << v + 1;
            if constexpr (DirectedV)
            {
                os << " ->";
            }
            else
            {
                os << " --";
            }
            for (auto u : _adj_list[v])
            {
                os << ' ' << u + 1;
            }
            os << std::endl;
        }
    }

private:
    /**
     * Vertex iterator implementation.
     */
    class _VertexIterator
    {
    public:
        _VertexIterator() = default;

        explicit _VertexIterator(Vertex v)
            : _v(v) {}

        _VertexIterator(const _VertexIterator &other) = default;

        _VertexIterator(_VertexIterator &&other) = default;

        _VertexIterator &operator=(const _VertexIterator &other) = default;

        _VertexIterator &operator=(_VertexIterator &&other) = default;

        _VertexIterator &operator++()
        {
            ++_v;
            return *this;
        }

        _VertexIterator operator++(int)
        {
            _VertexIterator it = *this;
            ++(*this);
            return it;
        }

        bool operator==(const _VertexIterator &other) const
        {
            return _v == other._v;
        }

        bool operator!=(const _VertexIterator &other) const
        {
            return !(*this == other);
        }

        Vertex operator*() const
        {
            return _v;
        }

    private:
        Vertex _v = 0;
    };

    /**
     * Edge iterator implementation.
     */
    class _EdgeIterator
    {
    private:
        using _AdjList = std::vector<std::vector<Vertex>>;
        using _AdjListReference = std::reference_wrapper<const _AdjList>;

    public:
        _EdgeIterator() = default;

        explicit _EdgeIterator(const _AdjList &adj_list, Vertex v)
            : _adj_list(adj_list), _v(v), _a(0)
        {
            ++(*this);
        }

        _EdgeIterator(const _EdgeIterator &other) = default;

        _EdgeIterator(_EdgeIterator &&other) = default;

        _EdgeIterator &operator=(const _EdgeIterator &other) = default;

        _EdgeIterator &operator=(_EdgeIterator &&other) = default;

        _EdgeIterator &operator++()
        {
            while (_v < std::size(_adj_list.get()))
            {
                while (++_a <= std::size(_adj_list.get()[_v]))
                {
                    if constexpr (DirectedV)
                    {
                        return *this;
                    }
                    else
                    {
                        if (_adj_list.get()[_v][_a - 1] >= _v)
                        {
                            return *this;
                        }
                    }
                }
                ++_v;
                _a = 0;
            }
            return *this;
        }

        _EdgeIterator operator++(int)
        {
            _EdgeIterator it = *this;
            ++(*this);
            return it;
        }

        bool operator==(const _EdgeIterator &other) const
        {
            return _v == other._v && _a == other._a;
        }

        bool operator!=(const _EdgeIterator &other) const
        {
            return !(*this == other);
        }

        Edge operator*() const
        {
            return {_v, _adj_list.get()[_v][_a - 1]};
        }

    private:
        _AdjListReference _adj_list;
        Vertex _v = 0;
        Vertex _a = 0;
    };

    /**
     * Outgoing edge iterator implementation.
     */
    class _OutEdgeIterator
    {
    private:
        using _AdjListIterator = std::vector<Vertex>::const_iterator;

    public:
        _OutEdgeIterator() = default;

        explicit _OutEdgeIterator(Vertex u, _AdjListIterator v)
            : _u(u), _v(v) {}

        _OutEdgeIterator(const _OutEdgeIterator &other) = default;

        _OutEdgeIterator(_OutEdgeIterator &&other) = default;

        _OutEdgeIterator &operator=(const _OutEdgeIterator &other) = default;

        _OutEdgeIterator &operator=(_OutEdgeIterator &&other) = default;

        _OutEdgeIterator &operator++()
        {
            ++_v;
            return *this;
        }

        _OutEdgeIterator operator++(int)
        {
            _OutEdgeIterator it = *this;
            ++(*this);
            return it;
        }

        bool operator==(const _OutEdgeIterator &other) const
        {
            return _u == other._u && _v == other._v;
        }

        bool operator!=(const _OutEdgeIterator &other) const
        {
            return !(*this == other);
        }

        Edge operator*() const
        {
            return {_u, *_v};
        }

    private:
        Vertex _u = 0;
        _AdjListIterator _v;
    };

    /**
     * Iterator range implementation.
     */
    template <typename IteratorT>
    struct Range
    {
        using Iterator = IteratorT;
        using ReverseIterator = std::reverse_iterator<IteratorT>;

        Iterator begin()
        {
            return first;
        }

        Iterator end()
        {
            return last;
        }

        ReverseIterator rbegin()
        {
            return std::make_reverse_iterator(last);
        }

        ReverseIterator rend()
        {
            return std::make_reverse_iterator(begin);
        }

        std::size_t size() const
        {
            return std::distance(first, last);
        }

        template <typename T>
        operator std::tuple<T, T>()
        {
            return {first, last};
        }

        Iterator first;
        Iterator last;
    };

    /**
     * Wrapper for Depth-First Search visitors.
     */
    template <typename VisitorT>
    class _DFSVisitor
    {
    public:
        using Graph = AdjacencyList;

        _DFSVisitor(VisitorT &visitor)
            : _visitor(visitor) {}

        void start_vertex(Vertex v, const Graph &g)
        {
            _start_vertex(_visitor, v, g);
        }

        void discover_vertex(Vertex v, const Graph &g)
        {
            _discover_vertex(_visitor, v, g);
        }

        void examine_edge(Edge e, const Graph &g)
        {
            _examine_edge(_visitor, e, g);
        }

        void tree_edge(Edge e, const Graph &g)
        {
            _tree_edge(_visitor, e, g);
        }

        void back_edge(Edge e, const Graph &g)
        {
            _back_edge(_visitor, e, g);
        }

        void forward_or_cross_edge(Edge e, const Graph &g)
        {
            _forward_or_cross_edge(_visitor, e, g);
        }

        void finish_edge(Edge e, const Graph &g)
        {
            _finish_edge(_visitor, e, g);
        }

        void finish_vertex(Vertex v, const Graph &g)
        {
            _finish_vertex(_visitor, v, g);
        }

        bool terminate(Vertex v, const Graph &g)
        {
            return _terminate(_visitor, v, g);
        }

    private:
        template <typename T>
        auto _start_vertex(T &vis, Vertex v, const Graph &g)
            -> decltype(vis.start_vertex(v, g), void())
        {
            vis.start_vertex(v, g);
        }

        template <typename T>
        auto _discover_vertex(T &vis, Vertex v, const Graph &g)
            -> decltype(vis.discover_vertex(v, g), void())
        {
            vis.discover_vertex(v, g);
        }

        template <typename T>
        auto _examine_edge(T &vis, Edge e, const Graph &g)
            -> decltype(vis.examine_edge(e, g), void())
        {
            vis.examine_edge(e, g);
        }

        template <typename T>
        auto _tree_edge(T &vis, Edge e, const Graph &g)
            -> decltype(vis.tree_edge(e, g), void())
        {
            vis.tree_edge(e, g);
        }

        template <typename T>
        auto _back_edge(T &vis, Edge e, const Graph &g)
            -> decltype(vis.back_edge(e, g), void())
        {
            vis.back_edge(e, g);
        }

        template <typename T>
        auto _forward_or_cross_edge(T &vis, Edge e, const Graph &g)
            -> decltype(vis.forward_or_cross_edge(e, g), void())
        {
            vis.forward_or_cross_edge(e, g);
        }

        template <typename T>
        auto _finish_edge(T &vis, Edge e, const Graph &g)
            -> decltype(vis.finish_edge(e, g), void())
        {
            vis.finish_edge(e, g);
        }

        template <typename T>
        auto _finish_vertex(T &vis, Vertex v, const Graph &g)
            -> decltype(vis.finish_vertex(v, g), void())
        {
            vis.finish_vertex(v, g);
        }

        template <typename T>
        auto _terminate(T &vis, Vertex v, const Graph &g)
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

    std::vector<std::vector<Vertex>> _adj_list;
};

/**
 * Undirected graph type definition.
 */
using UndirectedGraph = AdjacencyList<false>;

/**
 * Directed graph type definition.
 */
using DirectedGraph = AdjacencyList<true>;

#endif