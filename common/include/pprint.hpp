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

#ifndef PPRINT_HPP
#define PPRINT_HPP

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

/**
 * Pretty-prints vectors.
 */
template <typename T>
inline std::string pprint(const std::vector<T> &vec)
{
    std::ostringstream os;
    std::size_t w = 0;

    for (const auto &v : vec)
    {
        std::ostringstream oss;
        oss << v;
        if (w < std::size(oss.str()))
        {
            w = std::size(oss.str());
        }
    }

    os << '[';
    for (auto it = std::begin(vec); it != std::end(vec); ++it)
    {
        if (it != std::begin(vec))
        {
            os << ' ';
        }
        os << std::setw(w) << *it;
    }
    os << ']';

    return os.str();
}

/**
 * Pretty-prints matrices.
 */
template <typename T>
inline std::string pprint(const std::vector<std::vector<T>> &vec)
{
    std::ostringstream os;
    std::size_t w = 0;

    for (const auto &row : vec)
    {
        for (const auto &v : row)
        {
            std::ostringstream oss;
            oss << v;
            if (w < std::size(oss.str()))
            {
                w = std::size(oss.str());
            }
        }
    }

    for (const auto &row : vec)
    {
        os << '[';
        for (auto it = std::begin(row); it != std::end(row); ++it)
        {
            if (it != std::begin(row))
            {
                os << ' ';
            }
            os << std::setw(w) << *it;
        }
        os << ']' << std::endl;
    }

    return os.str();
}

/**
 * Pretty-prints vectors generated by a function.
 */
template <typename FuncT>
inline std::string pprint_n(std::size_t n, FuncT fn)
{
    std::ostringstream os;
    std::size_t w = 0;

    for (std::size_t i = 0; i < n; ++i)
    {
        std::ostringstream oss;
        oss << fn(i);
        if (w < std::size(oss.str()))
        {
            w = std::size(oss.str());
        }
    }

    os << '[';
    for (std::size_t i = 0; i < n; ++i)
    {
        if (i > 0)
        {
            os << ' ';
        }
        os << std::setw(w) << fn(i);
    }
    os << ']';

    return os.str();
}

/**
 * Pretty-prints matrices generated by a function.
 */
template <typename FuncT>
inline std::string pprint_n(std::size_t m, std::size_t n, FuncT fn)
{
    std::ostringstream os;
    std::size_t w = 0;

    for (std::size_t i = 0; i < m; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            std::ostringstream oss;
            oss << fn(i, j);
            if (w < std::size(oss.str()))
            {
                w = std::size(oss.str());
            }
        }
    }

    for (std::size_t i = 0; i < m; ++i)
    {
        os << '[';
        for (std::size_t j = 0; j < n; ++j)
        {
            if (j > 0)
            {
                os << ' ';
            }
            os << std::setw(w) << fn(i, j);
        }
        os << ']' << std::endl;
    }

    return os.str();
}

/**
 * Prints a table from multiple vectors.
 */
template <typename... ColsT>
static std::string pprint_table(const std::vector<std::string> &names,
                                ColsT... cols)
{
    std::ostringstream os;
    std::string sep = "  ";
    std::vector<std::size_t> w(sizeof...(cols), 0);

    {
        std::size_t j = 0;

        auto get_width = [&](const auto &col) {
            for (const auto &item : col)
            {
                std::ostringstream oss;
                oss << item;
                w[j] = std::max({w[j],
                                 std::size(names[j]),
                                 std::size(oss.str())});
            }
            ++j;
        };

        (get_width(cols), ...);
    }

    std::size_t total_width = std::accumulate(std::begin(w),
                                              std::end(w),
                                              std::size(sep) * (std::size(w) - 1));

    os << std::string(total_width, '-') << std::endl;

    std::size_t m = 0;

    ((m = std::max(m, std::size(cols))), ...);

    for (std::size_t j = 0; j < sizeof...(cols); ++j)
    {
        if (j > 0)
        {
            os << sep;
        }
        os << std::setw(w[j]) << names[j];
    }
    os << std::endl;

    os << std::string(total_width, '-') << std::endl;

    for (std::size_t i = 0; i < m; ++i)
    {
        std::size_t j = 0;

        auto print_row = [&](const auto &col) {
            if (j > 0)
            {
                os << sep;
            }
            os << std::setw(w[j]) << col[i];
            ++j;
        };

        (print_row(cols), ...);
        os << std::endl;
    }

    os << std::string(total_width, '-') << std::endl;

    return os.str();
}

#endif