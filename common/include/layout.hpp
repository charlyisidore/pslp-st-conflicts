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

#ifndef LAYOUT_HPP
#define LAYOUT_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

/**
 * Structure to manage layouts of items.
 */
class Layout
{
private:
    template <typename IteratorT>
    struct Range;

public:
    /**
     * Default constructor.
     */
    Layout() = default;

    /**
     * Creates an empty layout.
     */
    Layout(std::size_t n_items, std::size_t n_stacks)
        : _layout(n_stacks),
          _stack_of(n_items, n_stacks),
          _level_of(n_items, n_items),
          _capacity(n_items) {}

    /**
     * Creates an empty layout.
     */
    Layout(std::size_t n_items, std::size_t n_stacks, std::size_t capacity)
        : _layout(n_stacks),
          _stack_of(n_items, n_stacks),
          _level_of(n_items, n_items),
          _capacity(capacity) {}

    /**
     * Copy constructor.
     */
    Layout(const Layout &layout) = default;

    /**
     * Move constructor.
     */
    Layout(Layout &&layout) = default;

    /**
     * Copy assignment operator.
     */
    Layout &operator=(const Layout &layout) = default;

    /**
     * Move assignment operator.
     */
    Layout &operator=(Layout &&layout) = default;

    /**
     * Creates an empty layout.
     */
    void initialize(std::size_t n_items, std::size_t n_stacks)
    {
        _layout.clear();
        _layout.resize(n_stacks);
        _stack_of.clear();
        _stack_of.resize(n_items, n_stacks);
        _level_of.clear();
        _level_of.resize(n_items, n_items);
        _capacity = n_items;
    }

    /**
     * Creates an empty layout.
     */
    void initialize(std::size_t n_items,
                    std::size_t n_stacks,
                    std::size_t capacity)
    {
        initialize(n_items, n_stacks);
        _capacity = capacity;
    }

    /**
     * Removes all the items.
     */
    void clear()
    {
        initialize(n_items(), n_stacks(), capacity());
    }

    /**
     * Returns the number of items.
     */
    std::size_t n_items() const
    {
        return std::size(_stack_of);
    }

    /**
     * Returns the number of stacks.
     */
    std::size_t n_stacks() const
    {
        return std::size(_layout);
    }

    /**
     * Returns the maximum number of items in a stack.
     */
    std::size_t capacity() const
    {
        return _capacity;
    }

    /**
     * Returns the height of a stack.
     */
    std::size_t height(std::size_t k) const
    {
        assert(k < n_stacks());

        return std::size(_layout[k]);
    }

    /**
     * Returns the number of items in the storage.
     */
    std::size_t n_stored_items() const
    {
        std::size_t n = 0;
        for (const auto &s : _layout)
        {
            n += std::size(s);
        }
        return n;
    }

    /**
     * Returns the number of used stacks.
     */
    std::size_t n_used_stacks() const
    {
        return std::count_if(std::begin(_layout),
                             std::end(_layout),
                             [](const auto &s) {
                                 return !std::empty(s);
                             });
    }

    /**
     * Returns the topmost item of a stack.
     */
    std::size_t top(std::size_t k) const
    {
        assert(k < n_stacks());
        assert(!std::empty(_layout[k]));

        return _layout[k].back();
    }

    /**
     * Checks whether an item is placed.
     */
    bool contains(std::size_t i) const
    {
        assert(i < n_items());

        return _stack_of[i] < n_stacks();
    }

    /**
     * Returns the item at given stack and level.
     */
    std::size_t at(std::size_t k, std::size_t l) const
    {
        assert(k < n_stacks());
        assert(l < height(k));

        return _layout[k][l];
    }

    /**
     * Returns the stack of an item.
     */
    std::size_t stack_of(std::size_t i) const
    {
        assert(i < n_items());
        assert(_stack_of[i] < n_stacks());

        return _stack_of[i];
    }

    /**
     * Returns the level of an item.
     */
    std::size_t level_of(std::size_t i) const
    {
        assert(i < n_items());
        assert(_level_of[i] < n_items());

        return _level_of[i];
    }

    /**
     * Returns a range of iterators to explore a stack from bottom to top.
     */
    auto stack(std::size_t k) const
    {
        using range = Range<std::vector<std::size_t>::const_iterator>;

        assert(k < n_stacks());

        return range{std::begin(_layout[k]), std::end(_layout[k])};
    }

    /**
     * Returns a range of iterators to explore a stack from top to bottom.
     */
    auto rstack(std::size_t k) const
    {
        using range = Range<std::vector<std::size_t>::const_reverse_iterator>;

        assert(k < n_stacks());

        return range{std::rbegin(_layout[k]), std::rend(_layout[k])};
    }

    /**
     * Puts an item on top of a stack.
     */
    void push(std::size_t k, std::size_t i)
    {
        assert(i < n_items());
        assert(k < n_stacks());
        assert(_stack_of[i] == n_stacks());
        assert(_level_of[i] == n_items());
        assert(height(k) < capacity());

        _stack_of[i] = k;
        _level_of[i] = std::size(_layout[k]);
        _layout[k].push_back(i);
    }

    /**
     * Removes the topmost item from a stack.
     */
    std::size_t pop(std::size_t k)
    {
        assert(k < n_stacks());
        assert(!std::empty(_layout[k]));

        auto i = _layout[k].back();

        assert(_stack_of[i] < n_stacks());
        assert(_level_of[i] < n_items());

        _stack_of[i] = n_stacks();
        _level_of[i] = n_items();
        _layout[k].pop_back();
        return i;
    }

    /**
     * Pretty prints a layout.
     */
    void print(std::ostream &os = std::cout) const
    {
        print([](auto i) { return i + 1; }, os);
    }

    /**
     * Pretty prints a layout, prints data returned by a function.
     */
    template <typename F>
    void print(F f, std::ostream &os = std::cout) const
    {
        for (std::size_t k = 0; k < std::size(_layout); ++k)
        {
            os << k + 1 << " >";
            for (auto i : _layout[k])
            {
                os << ' ' << f(i);
            }
            os << std::endl;
        }
    }

    /**
     * Prints the contents of the structure.
     */
    void dump(std::ostream &os = std::cout) const
    {
        os << "n_items: " << n_items() << std::endl
           << "n_stacks: " << n_stacks() << std::endl
           << "stack_of:";
        for (auto k : _stack_of)
        {
            os << ' ' << k;
        }
        os << std::endl
           << "level_of:";
        for (auto l : _level_of)
        {
            os << ' ' << l;
        }
        os << std::endl
           << "layout:" << std::endl;
        for (const auto &s : _layout)
        {
            os << '-';
            for (auto i : s)
            {
                os << ' ' << i;
            }
            os << std::endl;
        }
    }

private:
    /**
     * Iterator range implementation.
     */
    template <typename IteratorT>
    struct Range
    {
        using Iterator = IteratorT;

        Iterator begin()
        {
            return first;
        }

        Iterator end()
        {
            return last;
        }

        Iterator first;
        Iterator last;
    };

    std::vector<std::vector<std::size_t>> _layout;
    std::vector<std::size_t> _stack_of;
    std::vector<std::size_t> _level_of;
    std::size_t _capacity;
};

#endif