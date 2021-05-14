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