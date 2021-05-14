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

#ifndef PSLP_MIP_MIP_CPLEX_HPP
#define PSLP_MIP_MIP_CPLEX_HPP

#include "graph.hpp"
#include "layout.hpp"
#include "square_matrix.hpp"
#include <ilcplex/ilocplex.h>
#include <nlohmann/json.hpp>
#include <utility>
#include <vector>

namespace pslp
{

namespace mip
{

/**
 * Base class for CPLEX-based MIP models.
 */
class MipCplex
{
public:
    /**
     * Problem type definition.
     */
    struct Problem
    {
        std::size_t n_items = 0;                    // Number of items
        std::size_t n_stacks = 0;                   // Number of stacks
        std::size_t capacity = 0;                   // Capacity per stack
        std::vector<std::size_t> arrivals;          // Arrival time of items
        std::vector<double> weights;                // Weight of items, if available
        std::vector<std::size_t> initial_positions; // Stack number of initial items
        SquareMatrix<bool> conflict_matrix;         // Conflict constraints
        DirectedGraph conflict_graph;               // Conflict graph
        Layout initial_layout;                      // The initial layout
    };

    /**
     * Solution type definition.
     */
    struct Solution
    {
        double obj = 0.;
        std::vector<std::size_t> stack_of;
        std::vector<std::size_t> level_of;
    };

    /**
     * Constructor.
     */
    MipCplex() = default;

    /**
     * Destructor.
     */
    virtual ~MipCplex();

    /**
     * Read the problem.
     */
    Problem read(const nlohmann::json &json);

    /**
     * Find an initial heuristic solution, also determines an upper bound.
     */
    void heuristic(const Problem &prob);

    /**
     * Build the model.
     */
    void build(const Problem &prob);

    /**
     * Solve the model.
     */
    void solve();

    /**
     * Check whether the solver found a solution.
     */
    bool has_solution() const;

    /**
     * Check whether a solution is valid.
     */
    bool check(const Problem &prob, const Solution &sol);

    /**
     * Write the LP in a file.
     */
    void write_lp(const char *filename);

    /**
     * Create the model.
     */
    virtual void create(const Problem &prob) = 0;

    /**
     * Create an initial solution.
     */
    virtual void initialize(const Problem &prob) = 0;

    /**
     * Return the solution.
     */
    virtual Solution solution(const Problem &prob) = 0;

    /**
     * Write the solution and statistics.
     */
    virtual void write(const Problem &prob,
                       const Solution &sol,
                       nlohmann::json &json);

protected:
    struct Incumbent
    {
        double time;
        IloNum obj;
        IloNum gap;
        IloInt64 n_nodes;
        IloInt64 n_remaining_nodes;
        IloCplex::IncumbentCallbackI::SolutionSource source;
    };

    IloEnv _env;
    IloCplex _cplex;
    IloModel _model;
    IloObjective _objective;
    IloBool _has_solution = false;
    double _reading_time = 0.;
    IloInt _reading_memory = 0;
    IloInt _reading_total_memory = 0;
    double _solving_time = 0.;
    IloInt _solving_memory = 0;
    IloInt _solving_total_memory = 0;
    std::vector<Incumbent> _incumbents;
    std::vector<std::vector<std::size_t>> _heuristic_sol;
};

} // namespace mip

} // namespace pslp

#endif
