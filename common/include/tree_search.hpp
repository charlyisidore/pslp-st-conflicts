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

#ifndef TREE_SEARCH_HPP
#define TREE_SEARCH_HPP

#include <cassert>
#include <chrono>
#include <iostream>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

namespace tree_search
{

// Forward declarations
template <typename NodeT>
class Context;

/**
 * Base class for tree search nodes.
 */
class Node
{
public:
    /**
     * Destructor.
     */
    virtual ~Node() {}

    /**
     * Return the index of the node.
     */
    std::size_t index() const
    {
        return _index;
    }

    /**
     * Return the depth of the node.
     */
    std::size_t depth() const
    {
        return _depth;
    }

    /**
     * Return the parent of the node.
     */
    std::size_t parent() const
    {
        return _parent;
    }

    /**
     * Set the index of the node.
     */
    void set_index(std::size_t index)
    {
        _index = index;
    }

    /**
     * Set the depth of the node.
     */
    void set_depth(std::size_t depth)
    {
        _depth = depth;
    }

    /**
     * Set the parent of the node.
     */
    void set_parent(std::size_t parent)
    {
        _parent = parent;
    }

private:
    std::size_t _index = 0;
    std::size_t _depth = 0;
    std::size_t _parent = 0;
};

/**
 * Base class for tree search queues.
 */
template <typename NodeT>
class Queue
{
public:
    static_assert(std::is_base_of_v<Node, NodeT>);

    /**
     * Destructor.
     */
    virtual ~Queue() {}

    /**
     * Check whether the queue is empty.
     */
    virtual bool empty() const = 0;

    /**
     * Add a node to the queue.
     */
    virtual void add(NodeT &&) = 0;

    /**
     * Extract the highest priority node from the queue.
     */
    virtual NodeT extract() = 0;
};

/**
 * Base class for tree search algorithms.
 */
template <typename NodeT>
class Algorithm
{
public:
    static_assert(std::is_base_of_v<Node, NodeT>);

    /**
     * Context type definition.
     */
    using ContextT = Context<NodeT>;

    /**
     * Destructor.
     */
    virtual ~Algorithm() {}

    /**
     * Create root nodes.
     */
    virtual void start(ContextT &) = 0;

    /**
     * Create the final solution.
     */
    virtual void finish(ContextT &) = 0;

    /**
     * Estimate a node before adding it to the queue.
     */
    virtual bool estimate(ContextT &, NodeT &) = 0;

    /**
     * Evaluate a node before branching.
     */
    virtual bool evaluate(ContextT &, NodeT &) = 0;

    /**
     * Create child nodes.
     */
    virtual void branch(ContextT &, NodeT &) = 0;

    /**
     * Check whether the algorithm should finish.
     */
    virtual bool termination(ContextT &)
    {
        return false;
    }

    /**
     * Print a node.
     */
    virtual void print(const NodeT &node, std::ostream &os) const
    {
        os << '[' << node.index() << ']'
           << " depth:" << node.depth()
           << " parent:" << node.parent();
    }
};

/**
 * Base class for visitors.
 */
template <typename NodeT>
class Visitor
{
public:
    static_assert(std::is_base_of_v<Node, NodeT>);

    /**
     * Context type definition.
     */
    using ContextT = Context<NodeT>;

    /**
     * Destructor.
     */
    virtual ~Visitor() {}

    /**
     * Called when the tree search is started.
     */
    virtual void start(const ContextT &) {}

    /**
     * Called when the tree search is going to finish.
     */
    virtual void finish(const ContextT &) {}

    /**
     * Called when a node is going to be estimated.
     */
    virtual void node_estimate(const ContextT &, const NodeT &) {}

    /**
     * Called when a node is going to be evaluated.
     */
    virtual void node_evaluate(const ContextT &, const NodeT &) {}

    /**
     * Called when a node is going to be added to the queue.
     */
    virtual void node_enqueue(const ContextT &, const NodeT &) {}

    /**
     * Called when a node is going to be branched.
     */
    virtual void node_branch(const ContextT &, const NodeT &) {}

    /**
     * Called when a node is going to be pruned.
     */
    virtual void node_prune(const ContextT &, const NodeT &) {}
};

/**
 * Context type definition.
 */
template <typename NodeT>
struct Context
{
    static_assert(std::is_base_of_v<Node, NodeT>);

    /**
     * Algorithm type definition.
     */
    using AlgorithmT = Algorithm<NodeT>;

    /**
     * Queue type definition.
     */
    using QueueT = Queue<NodeT>;

    /**
     * Visitor type definition.
     */
    using VisitorT = Visitor<NodeT>;

    /**
     * Constructor.
     */
    Context(AlgorithmT &algorithm, QueueT &queue, VisitorT &visitor)
        : _algorithm(algorithm), _queue(queue), _visitor(visitor) {}

    /**
     * Estimate and add a root node to the queue.
     */
    void add_root(NodeT &&node)
    {
        node.set_index(++_n_total_nodes);
        node.set_depth(0);
        node.set_parent(0);
        _estimate(std::move(node));
    }

    /**
     * Estimate and add a child node to the queue.
     */
    void add_child(NodeT &&node, const NodeT &parent)
    {
        node.set_index(++_n_total_nodes);
        node.set_depth(parent.depth() + 1);
        node.set_parent(parent.index());
        _estimate(std::move(node));
    }

    /**
     * Print a node.
     */
    void print(const NodeT &node, std::ostream &os = std::cout) const
    {
        _algorithm.print(node, os);
    }

    /**
     * Return the total number of nodes.
     */
    std::size_t n_total_nodes() const
    {
        return _n_total_nodes;
    }

    /**
     * Return the elapsed time in seconds from the start.
     */
    double elapsed_time() const
    {
        auto now = std::chrono::steady_clock::now();
        return std::chrono::duration<double>(now - _start_time_point).count();
    }

    /**
     * Run the algorithm.
     */
    void run()
    {
        _start_time_point = std::chrono::steady_clock::now();
        _algorithm.start(*this);
        _visitor.start(*this);
        while (!_queue.empty() && !_algorithm.termination(*this))
        {
            auto node = std::move(_queue.extract());
            _evaluate(std::move(node));
        }
        _visitor.finish(*this);
        _algorithm.finish(*this);
    }

    /**
     * Access to the algorithm.
     */
    const AlgorithmT &algorithm() const
    {
        return _algorithm;
    }

    /**
     * Access to the algorithm.
     */
    AlgorithmT &algorithm()
    {
        return _algorithm;
    }

    /**
     * Access to the queue.
     */
    const QueueT &queue() const
    {
        return _queue;
    }

    /**
     * Access to the queue.
     */
    QueueT &queue()
    {
        return _queue;
    }

    /**
     * Access to the visitor.
     */
    const VisitorT &visitor() const
    {
        return _visitor;
    }

    /**
     * Access to the visitor.
     */
    VisitorT &visitor()
    {
        return _visitor;
    }

private:
    /**
     * Estimate and add a node to the queue.
     */
    void _estimate(NodeT &&node)
    {
        _visitor.node_estimate(*this, node);

        if (_algorithm.estimate(*this, node))
        {
            _visitor.node_enqueue(*this, node);
            _queue.add(std::move(node));
        }
        else
        {
            _visitor.node_prune(*this, node);
        }
    }

    /**
     * Evaluate a node and create child nodes.
     */
    void _evaluate(NodeT &&node)
    {
        _visitor.node_evaluate(*this, node);

        if (_algorithm.evaluate(*this, node))
        {
            _visitor.node_branch(*this, node);
            _algorithm.branch(*this, node);
        }
        else
        {
            _visitor.node_prune(*this, node);
        }
    }

    AlgorithmT &_algorithm;
    QueueT &_queue;
    VisitorT &_visitor;
    std::size_t _n_total_nodes = 0;
    std::chrono::time_point<std::chrono::steady_clock> _start_time_point;
};

/**
 * Compare values using given tolerance.
 */
template <typename ValueT = double>
class Tolerance
{
public:
    /**
     * Constructor.
     */
    Tolerance(ValueT eps = _default())
        : _eps(eps) {}

    /**
     * Get the tolerance value.
     */
    ValueT get() const
    {
        return _eps;
    }

    /**
     * Set the tolerance value.
     */
    void set(ValueT eps)
    {
        _eps = eps;
    }

    /**
     * Check whether x and y are equal.
     */
    bool eq(ValueT x, ValueT y) const
    {
        return std::abs(x - y) <= _eps;
    }

    /**
     * Check whether x is less than y.
     */
    bool lt(ValueT x, ValueT y) const
    {
        return x - y < -_eps;
    }

    /**
     * Check whether x is less than or equal to y.
     */
    bool le(ValueT x, ValueT y) const
    {
        return x - y <= _eps;
    }

    /**
     * Check whether x is greater than y.
     */
    bool gt(ValueT x, ValueT y) const
    {
        return x - y > _eps;
    }

    /**
     * Check whether x is greater than or equal to y.
     */
    bool ge(ValueT x, ValueT y) const
    {
        return x - y >= -_eps;
    }

private:
    /**
     * Return the default tolerance value.
     */
    static ValueT _default()
    {
        if constexpr (std::is_floating_point_v<ValueT>)
        {
            return 1e-6;
        }
        else
        {
            return 0;
        }
    }

    ValueT _eps;
};

/**
 * Basic queue based on std::set and a comparator.
 */
template <typename NodeT, typename CompareT>
class BasicQueue : public Queue<NodeT>
{
public:
    /**
     * Constructor.
     */
    BasicQueue(const CompareT &comp = CompareT())
        : _set(comp) {}

    /**
     * Check whether the queue is empty.
     */
    bool empty() const override
    {
        return std::empty(_set);
    }

    /**
     * Add a node to the queue.
     */
    void add(NodeT &&node) override
    {
        _set.insert(std::move(node));
    }

    /**
     * Extract the highest priority node from the queue.
     */
    NodeT extract() override
    {
        assert(!std::empty(_set));
        auto it = std::begin(_set);
        assert(it != std::end(_set));
        return std::move(_set.extract(it).value());
    }

protected:
    std::set<NodeT, CompareT> _set;
};

/**
 * Best First Search comparator.
 */
template <typename CompareT = std::less<>>
class BFSCompare
{
public:
    BFSCompare(const CompareT &comp = CompareT())
        : _comp(comp) {}

    template <typename NodeT>
    bool operator()(const NodeT &node1, const NodeT &node2) const
    {
        return _comp(node1, node2) ||
               (!_comp(node2, node1) &&
                node1.index() < node2.index());
    }

private:
    CompareT _comp;
};

/**
 * Depth First Search comparator.
 */
template <typename CompareT = std::less<>>
class DFSCompare
{
public:
    DFSCompare(const CompareT &comp = CompareT())
        : _comp(comp) {}

    template <typename NodeT>
    bool operator()(const NodeT &node1, const NodeT &node2) const
    {
        return node1.depth() > node2.depth() ||
               (node1.depth() == node2.depth() &&
                (_comp(node1, node2) ||
                 (!_comp(node2, node1) &&
                  node1.index() < node2.index())));
    }

private:
    CompareT _comp;
};

/**
 * Breadth First Search comparator.
 */
template <typename CompareT = std::less<>>
class BrFSCompare
{
public:
    BrFSCompare(const CompareT &comp = CompareT())
        : _comp(comp) {}

    template <typename NodeT>
    bool operator()(const NodeT &node1, const NodeT &node2) const
    {
        return node1.depth() < node2.depth() ||
               (node1.depth() == node2.depth() &&
                (_comp(node1, node2) ||
                 (!_comp(node2, node1) &&
                  node1.index() < node2.index())));
    }

private:
    CompareT _comp;
};

/**
 * Cyclic Best First Search queue.
 */
template <typename NodeT, typename CompareT = std::less<>>
class CyclicBestFirstSearchQueue : public Queue<NodeT>
{
public:
    /**
     * Constructor.
     */
    CyclicBestFirstSearchQueue(const CompareT &comp = CompareT())
        : _comp(comp) {}

    /**
     * Check whether the queue is empty.
     */
    bool empty() const override
    {
#ifndef NDEBUG
        bool is_empty = true;
        for (const auto &set : _sets)
        {
            is_empty = is_empty && std::empty(set);
        }
        assert(is_empty == (_size == 0));
#endif
        return _size == 0;
    }

    /**
     * Add a node to the queue.
     */
    void add(NodeT &&node) override
    {
        const auto d = node.depth();
        if (std::size(_sets) <= d)
        {
            _sets.resize(d + 1,
                         std::set<NodeT, BFSCompare<CompareT>>(BFSCompare<CompareT>(_comp)));
        }
        auto [it, inserted] = _sets[d].insert(std::move(node));
        if (inserted)
        {
            ++_size;
        }
#ifndef NDEBUG
        std::size_t new_size = 0;
        for (const auto &set : _sets)
        {
            new_size += std::size(set);
        }
        assert(new_size == _size);
#endif
    }

    /**
     * Extract the highest priority node from the queue.
     */
    NodeT extract() override
    {
        assert(_size > 0);
        if (_depth >= std::size(_sets))
        {
            _depth = 0;
        }
        while (std::empty(_sets[_depth]))
        {
            _depth = (_depth + 1) % std::size(_sets);
        }
        auto it = std::begin(_sets[_depth]);
        assert(it != std::end(_sets[_depth]));
        --_size;
        return std::move(_sets[_depth++].extract(it).value());
    }

protected:
    CompareT _comp;
    std::vector<std::set<NodeT, BFSCompare<CompareT>>> _sets;
    std::size_t _depth = 0;
    std::size_t _size = 0;
};

/**
 * Beam Search queue.
 */
template <typename NodeT, typename CompareT = std::less<>>
class BeamSearchQueue : public BasicQueue<NodeT, BrFSCompare<CompareT>>
{
public:
    using Parent = BasicQueue<NodeT, BrFSCompare<CompareT>>;
    using Parent::_set;

    BeamSearchQueue(std::size_t width, const CompareT &comp = CompareT())
        : Parent(BrFSCompare<CompareT>(comp)),
          _width(width) {}

    void add(NodeT &&node) override
    {
        _set.insert(std::move(node));
        if (std::size(_set) > _width)
        {
            _set.erase(std::prev(std::end(_set)));
        }
    }

protected:
    std::size_t _width;
};

/**
 * A visitor that prints events.
 */
template <typename NodeT>
class PrintVisitor : public Visitor<NodeT>
{
    class IndentStreambuf;

public:
    static_assert(std::is_base_of_v<Node, NodeT>);

    /**
     * Context type definition.
     */
    using ContextT = Context<NodeT>;

    /**
     * Constructor.
     */
    PrintVisitor(std::size_t tab_width = 1,
                 char c = ' ',
                 std::ostream &os = std::cout)
        : _os(os), _tw(tab_width), _c(c) {}

    /**
     * Called when the tree search is started.
     */
    virtual void start(const ContextT &)
    {
        _os << "start tree search" << std::endl;
    }

    /**
     * Called when the tree search is going to finish.
     */
    virtual void finish(const ContextT &)
    {
        _os << "end tree search" << std::endl;
    }

    /**
     * Called when a node is going to be estimated.
     */
    virtual void node_estimate(const ContextT &ctx, const NodeT &node)
    {
        _print("estimate", ctx, node);
    }

    /**
     * Called when a node is going to be evaluated.
     */
    virtual void node_evaluate(const ContextT &ctx, const NodeT &node)
    {
        _print("evaluate", ctx, node);
    }

    /**
     * Called when a node is going to be added to the queue.
     */
    virtual void node_enqueue(const ContextT &ctx, const NodeT &node)
    {
        _print("enqueue", ctx, node);
    }

    /**
     * Called when a node is going to be branched.
     */
    virtual void node_branch(const ContextT &ctx, const NodeT &node)
    {
        _print("branch", ctx, node);
    }

    /**
     * Called when a node is going to be pruned.
     */
    virtual void node_prune(const ContextT &ctx, const NodeT &node)
    {
        _print("prune", ctx, node);
    }

private:
    /**
     * Print an event.
     */
    virtual void _print(const char *name,
                        const ContextT &ctx,
                        const NodeT &node)
    {
        _os << std::string(_tw * node.depth(), _c) << name << std::endl;

        IndentStreambuf buf(_os, _tw * (node.depth() + 1), _c);
        ctx.print(node, _os);
        _os << std::endl;
    }

    /**
     * Make output indented.
     */
    class IndentStreambuf : public std::streambuf
    {
    public:
        explicit IndentStreambuf(std::ostream &os,
                                 std::size_t width,
                                 char ch = ' ')
            : _os(&os), _buf(_os->rdbuf()), _blank(width, ch)
        {
            _os->rdbuf(this);
        }

        explicit IndentStreambuf(std::streambuf *buf,
                                 std::size_t width,
                                 char ch = ' ')
            : _buf(buf), _blank(width, ch) {}

        ~IndentStreambuf()
        {
            if (_os)
            {
                _os->rdbuf(_buf);
            }
        }

    protected:
        int overflow(int ch) override
        {
            if (_begin && ch != '\n')
            {
                _buf->sputn(std::data(_blank), std::size(_blank));
            }
            _begin = (ch == '\n');
            return _buf->sputc(ch);
        }

    private:
        std::ostream *_os = nullptr;
        std::streambuf *_buf = nullptr;
        std::string _blank;
        bool _begin = true;
    };

    std::ostream &_os;
    std::size_t _tw = 1;
    char _c;
};

/**
 * Create a Best First Search queue.
 */
template <typename NodeT, typename CompareT = std::less<>>
inline BasicQueue<NodeT, BFSCompare<CompareT>> bfs_queue(const CompareT &comp = CompareT())
{
    return BasicQueue<NodeT, BFSCompare<CompareT>>(BFSCompare<CompareT>(comp));
}

/**
 * Create a Depth First Search queue.
 */
template <typename NodeT, typename CompareT = std::less<>>
inline BasicQueue<NodeT, DFSCompare<CompareT>> dfs_queue(const CompareT &comp = CompareT())
{
    return BasicQueue<NodeT, DFSCompare<CompareT>>(DFSCompare<CompareT>(comp));
}

/**
 * Create a Breadth First Search queue.
 */
template <typename NodeT, typename CompareT = std::less<>>
inline BasicQueue<NodeT, BrFSCompare<CompareT>> brfs_queue(const CompareT &comp = CompareT())
{
    return BasicQueue<NodeT, BrFSCompare<CompareT>>(BrFSCompare<CompareT>(comp));
}

/**
 * Create a Cyclic Best First Search queue.
 */
template <typename NodeT, typename CompareT = std::less<>>
inline CyclicBestFirstSearchQueue<NodeT, CompareT> cbfs_queue(const CompareT &comp = CompareT())
{
    return CyclicBestFirstSearchQueue<NodeT, CompareT>(comp);
}

/**
 * Create a Beam Search queue.
 */
template <typename NodeT, typename CompareT = std::less<>>
inline BeamSearchQueue<NodeT, CompareT> beam_search_queue(std::size_t width,
                                                          const CompareT &comp = CompareT())
{
    return BeamSearchQueue<NodeT, CompareT>(width, comp);
}

/**
 * Create a null visitor.
 */
template <typename NodeT>
inline Visitor<NodeT> null_visitor()
{
    return Visitor<NodeT>();
}

/**
 * Create a print visitor.
 */
template <typename NodeT>
inline PrintVisitor<NodeT> print_visitor(std::size_t tab_width = 1,
                                         char c = ' ')
{
    return PrintVisitor<NodeT>(tab_width, c);
}

/**
 * Run a tree search algorithm.
 */
template <typename NodeT>
inline void run(Algorithm<NodeT> &algo,
                Queue<NodeT> &queue,
                Visitor<NodeT> &visitor)
{
    static_assert(std::is_base_of_v<Node, NodeT>);

    Context<NodeT> ctx(algo, queue, visitor);
    ctx.run();
}

/**
 * Run a tree search algorithm.
 */
template <typename NodeT,
          typename VisitorT = Visitor<NodeT>>
inline void run(Algorithm<NodeT> &algo,
                Queue<NodeT> &queue,
                VisitorT visitor = VisitorT())
{
    static_assert(std::is_base_of_v<Node, NodeT>);
    static_assert(std::is_base_of_v<Visitor<NodeT>, VisitorT>);

    Context<NodeT> ctx(algo, queue, visitor);
    ctx.run();
}

/**
 * Run a tree search algorithm.
 */
template <typename NodeT,
          typename QueueT = BasicQueue<NodeT, BFSCompare<>>,
          typename VisitorT = Visitor<NodeT>>
inline void run(Algorithm<NodeT> &algo,
                QueueT queue = QueueT(),
                VisitorT visitor = VisitorT())
{
    static_assert(std::is_base_of_v<Node, NodeT>);
    static_assert(std::is_base_of_v<Queue<NodeT>, QueueT>);
    static_assert(std::is_base_of_v<Visitor<NodeT>, VisitorT>);

    Context<NodeT> ctx(algo, queue, visitor);
    ctx.run();
}

} // namespace tree_search

#endif