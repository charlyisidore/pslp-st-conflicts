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

#ifndef SYMMETRIC_MATRIX_HPP
#define SYMMETRIC_MATRIX_HPP

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

/**
 * A matrix such that $M(i, j) = M(j, i)$. 
 */
template <typename ValueT = double, bool IncludeDiag = true>
class SymmetricMatrix
{
public:
    /**
     * Value type definition.
     */
    using Value = ValueT;

    /**
     * Const reference type definition.
     */
    using ConstReference = typename std::vector<Value>::const_reference;

    /**
     * Reference type definition.
     */
    using Reference = typename std::vector<Value>::reference;

    /**
     * Default constructor.
     */
    SymmetricMatrix() = default;

    /**
     * Constructor.
     */
    SymmetricMatrix(std::size_t n, const Value &v = Value())
        : _data(index(0, n), v), _n(n)
    {
    }

    /**
     * Copy constructor.
     */
    SymmetricMatrix(const SymmetricMatrix &g) = default;

    /**
     * Move constructor.
     */
    SymmetricMatrix(SymmetricMatrix &&g) = default;

    /**
     * Copy assignment operator.
     */
    SymmetricMatrix &operator=(const SymmetricMatrix &g) = default;

    /**
     * Move assignment operator.
     */
    SymmetricMatrix &operator=(SymmetricMatrix &&g) = default;

    /**
     * Initializes a n*n matrix with default values.
     */
    void initialize(std::size_t n, const Value &v = Value())
    {
        _data.clear();
        _data.resize(index(0, n), v);
        _n = n;
    }

    /**
     * Returns the number of rows in the matrix.
     */
    std::size_t n_rows() const
    {
        return _n;
    }

    /**
     * Returns the number of columns in the matrix.
     */
    std::size_t n_cols() const
    {
        return _n;
    }

    /**
     * Returns the number of elements in the matrix.
     */
    std::size_t n_elements() const
    {
        return std::size(_data);
    }

    /**
     * Returns all the elements as a vector.
     */
    const auto &elements() const
    {
        return _data;
    }

    /**
     * Returns the element at coordinates (i, j).
     */
    ConstReference operator()(std::size_t i, std::size_t j) const
    {
        assert(i < _n);
        assert(j < _n);
        assert(index(i, j) < std::size(_data));
        assert(index(i, j) == index(j, i));

        return _data[index(i, j)];
    }

    /**
     * Returns the element at coordinates (i, j).
     */
    Reference operator()(std::size_t i, std::size_t j)
    {
        assert(i < _n);
        assert(j < _n);
        assert(index(i, j) < std::size(_data));
        assert(index(i, j) == index(j, i));

        return _data[index(i, j)];
    }

    /**
     * Returns the element at coordinates (i, j).
     */
    ConstReference at(std::size_t i, std::size_t j) const
    {
        assert(i < _n);
        assert(j < _n);
        assert(index(i, j) == index(j, i));

        return _data.at(index(i, j));
    }

    /**
     * Returns the element at coordinates (i, j).
     */
    Reference at(std::size_t i, std::size_t j)
    {
        assert(i < _n);
        assert(j < _n);
        assert(index(i, j) == index(j, i));

        return _data.at(index(i, j));
    }

    /**
     * Delete a row and column.
     */
    void remove(std::size_t i)
    {
        assert(i < _n);

        std::vector<std::size_t> to_remove;

        for (std::size_t j = 0; j < _n; ++j)
        {
            if (i != j || IncludeDiag)
            {
                to_remove.push_back(index(i, j));
            }
        }

        std::sort(std::rbegin(to_remove), std::rend(to_remove));

        assert(std::adjacent_find(std::begin(to_remove),
                                  std::end(to_remove)) == std::end(to_remove));

        for (auto k : to_remove)
        {
            _data.erase(std::begin(_data) + k);
        }

        --_n;

        assert(std::size(_data) == index(0, _n));
    }

    /**
     * Print the matrix.
     */
    void print(std::ostream &os = std::cout) const
    {
        for (std::size_t i = 0; i < _n; ++i)
        {
            for (std::size_t j = 0; IncludeDiag ? (j <= i) : (j < i); ++j)
            {
                if (j > 0)
                {
                    os << ' ';
                }
                os << _data[index(i, j)];
            }
            os << std::endl;
        }
    }

    /**
     * Returns the index of coordinates (i, j).
     */
    std::size_t index(std::size_t i, std::size_t j) const
    {
        if constexpr (IncludeDiag)
        {
            return i < j
                       ? j * (j + 1) / 2 + i
                       : i * (i + 1) / 2 + j;
        }
        else
        {
            assert(i != j);
            return i < j
                       ? j * (j - 1) / 2 + i
                       : i * (i - 1) / 2 + j;
        }
    }

private:
    std::vector<Value> _data;
    std::size_t _n = 0;
};

#endif