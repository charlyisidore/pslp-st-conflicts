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

#ifndef SQUARE_MATRIX_HPP
#define SQUARE_MATRIX_HPP

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

/**
 * A square matrix. 
 */
template <typename ValueT = double>
class SquareMatrix
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
    SquareMatrix() = default;

    /**
     * Constructor.
     */
    SquareMatrix(std::size_t n, const Value &v = Value())
        : _data(n, std::vector<Value>(n, v))
    {
    }

    /**
     * Copy constructor.
     */
    SquareMatrix(const SquareMatrix &g) = default;

    /**
     * Move constructor.
     */
    SquareMatrix(SquareMatrix &&g) = default;

    /**
     * Copy assignment operator.
     */
    SquareMatrix &operator=(const SquareMatrix &g) = default;

    /**
     * Move assignment operator.
     */
    SquareMatrix &operator=(SquareMatrix &&g) = default;

    /**
     * Initializes a n*n matrix with default values.
     */
    void initialize(std::size_t n, const Value &v = Value())
    {
        _data.clear();
        _data.resize(n, std::vector<Value>(n, v));
    }

    /**
     * Returns the number of rows in the matrix.
     */
    std::size_t n_rows() const
    {
        return std::size(_data);
    }

    /**
     * Returns the number of columns in the matrix.
     */
    std::size_t n_cols() const
    {
        return std::size(_data);
    }

    /**
     * Returns the element at coordinates (i, j).
     */
    ConstReference operator()(std::size_t i, std::size_t j) const
    {
        assert(i < std::size(_data));
        assert(j < std::size(_data));

        return _data[i][j];
    }

    /**
     * Returns the element at coordinates (i, j).
     */
    Reference operator()(std::size_t i, std::size_t j)
    {
        assert(i < std::size(_data));
        assert(j < std::size(_data));

        return _data[i][j];
    }

    /**
     * Returns the element at coordinates (i, j).
     */
    ConstReference at(std::size_t i, std::size_t j) const
    {
        assert(i < std::size(_data));
        assert(j < std::size(_data));

        return _data.at(i).at(j);
    }

    /**
     * Returns the element at coordinates (i, j).
     */
    Reference at(std::size_t i, std::size_t j)
    {
        assert(i < std::size(_data));
        assert(j < std::size(_data));

        return _data.at(i).at(j);
    }

    /**
     * Delete a row and column.
     */
    void remove(std::size_t i)
    {
        assert(i < std::size(_data));

        _data.erase(std::begin(_data) + i);

        for (auto &row : _data)
        {
            row.erase(std::begin(row) + i);
        }
    }

    /**
     * Print the matrix.
     */
    void print(std::ostream &os = std::cout) const
    {
        for (std::size_t i = 0; i < std::size(_data); ++i)
        {
            for (std::size_t j = 0; j < std::size(_data[i]); ++j)
            {
                if (j > 0)
                {
                    os << ' ';
                }
                os << _data[i][j];
            }
            os << std::endl;
        }
    }

private:
    std::vector<std::vector<Value>> _data;
};

#endif