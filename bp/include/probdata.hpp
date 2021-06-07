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

#ifndef PSLP_BP_PROBDATA_HPP
#define PSLP_BP_PROBDATA_HPP

#include "graph.hpp"
#include "layout.hpp"
#include "square_matrix.hpp"
#include <cassert>
#include <iostream>
#include <objscip/objprobdata.h>
#include <tuple>
#include <vector>

namespace pslp
{

namespace bp
{

/**
 * Problem data for the PSLP.
 */
class Probdata : public scip::ObjProbData
{
public:
    /**
     * Problem type definition.
     */
    enum class Type
    {
        MinNStacks,
        MinNBlockingItems
    };

    /**
     * Default constructor.
     */
    Probdata(Type type,
             std::size_t n,
             std::size_t m,
             std::size_t b,
             bool has_adj_conflicts);

    /**
     * Returns the type of problem.
     */
    auto type() const
    {
        return _type;
    }

    /**
     * Returns the number of items.
     */
    auto n_items() const
    {
        return std::size(_conss_assign);
    }

    /**
     * Returns the number of stacks.
     */
    auto n_stacks() const
    {
        return _n_stacks;
    }

    /**
     * Returns the capacity of stacks.
     */
    auto capacity() const
    {
        return _capacity;
    }

    /**
     * Checks whether the conflicts are adjacent or not.
     */
    auto has_adjacent_conflicts() const
    {
        return _has_adjacent_conflicts;
    }

    /**
     * Returns the arrival time of an item.
     */
    auto arrival(std::size_t i) const
    {
        assert(i < n_items());

        return _arrivals[i];
    }

    /**
     * Returns the weight of an item.
     */
    auto weight(std::size_t i) const
    {
        assert(i < std::size(_weights));

        return _weights[i];
    }

    /**
     * Returns the departure date of an item.
     */
    auto departure(std::size_t i) const
    {
        assert(i < std::size(_departures));

        return _departures[i];
    }

    /**
     * Returns the stack number of an initial item.
     */
    auto initial_position(std::size_t i) const
    {
        assert(i < std::size(_initial_positions));

        return _initial_positions[i];
    }

    /**
     * Checks whether the instance has weights.
     */
    auto has_weights() const
    {
        return !std::empty(_weights);
    }

    /**
     * Checks whether the instance has departure dates.
     */
    auto has_departures() const
    {
        return !std::empty(_departures);
    }

    /**
     * Returns the number of initial items.
     */
    auto n_initial_items() const
    {
        return std::size(_initial_positions);
    }

    /**
     * Returns the matrix of conflict constraints.
     */
    const auto &conflict_matrix() const
    {
        return _conflict_matrix;
    }

    /**
     * Returns the matrix of blocking constraints.
     */
    const auto &blocking_matrix() const
    {
        return _blocking_matrix;
    }

    /**
     * Returns the graph of conflict constraints.
     */
    const auto &conflict_graph() const
    {
        return _conflict_graph;
    }

    /**
     * Returns the graph complement of conflict constraints.
     */
    const auto &conflict_graph_complement() const
    {
        return _conflict_graph_complement;
    }

    /**
     * Returns the graph of blocking constraints.
     */
    const auto &blocking_graph() const
    {
        return _blocking_graph;
    }

    /**
     * Returns the initial layout.
     */
    const auto &initial_layout() const
    {
        return _initial_layout;
    }

    /**
     * Initialize arrival times of items.
     */
    void set_arrivals(const std::vector<std::size_t> &arrivals);

    /**
     * Initialize weights of items.
     */
    void set_weights(const std::vector<double> &weights);

    /**
     * Initialize departure times of items.
     */
    void set_departures(const std::vector<double> &departures);

    /**
     * Initialize initial positions.
     */
    void set_initial_positions(const std::vector<std::size_t> &init_pos);

    /**
     * Initialize the blocking matrix.
     */
    void set_blocking_matrix(const std::vector<std::vector<bool>> &matrix);

    /**
     * Initialize the conflict matrix.
     */
    void set_conflict_matrix(const std::vector<std::vector<bool>> &matrix);

    /**
     * Checks whether the stacking constraints are transitive.
     */
    bool is_stacking_transitive() const;

    /**
     * Transforms the stacking constraints to remove equivalent items.
     */
    [[nodiscard]] SCIP_RETCODE transform_conflicts();

    /**
     * Returns the dual value of an assignment constraint.
     */
    SCIP_Real dual_assign(SCIP *scip, std::size_t i, bool farkas) const;

    /**
     * Returns the dual value of the space constraint.
     */
    SCIP_Real dual_space(SCIP *scip, bool farkas) const;

    /**
     * Returns the SCIP variables of the master problem.
     */
    const auto &vars() const
    {
        return _vars;
    }

    /**
     * Checks if a variable exists in the master problem.
     */
    bool has_var(SCIP *scip, const std::vector<std::size_t> &s) const;

    /**
     * Checks whether a variable satisfies the constraints.
     */
    bool check_var(const std::vector<std::size_t> &s) const;

    /**
     * Get the objective coefficient of a variable.
     */
    double get_var_obj(const std::vector<std::size_t> &s) const;

    /**
     * Adds a variable to the master problem.
     */
    [[nodiscard]] SCIP_RETCODE add_var(SCIP *scip,
                                       const std::vector<std::size_t> &s,
                                       SCIP_PRICER *pricer = nullptr);

    /**
     * Prepare the problem data.
     */
    [[nodiscard]] SCIP_RETCODE prepare();

    /**
     * Initializes the model.
     */
    [[nodiscard]] SCIP_RETCODE initialize(SCIP *scip);

    /**
     * Releases the model.
     */
    [[nodiscard]] SCIP_RETCODE release(SCIP *scip);

    /**
     * Builds a solution.
     */
    [[nodiscard]] SCIP_RETCODE build_sol(SCIP *scip,
                                         SCIP_SOL *sol,
                                         std::vector<std::size_t> &stack_of,
                                         std::vector<std::size_t> &level_of) const;

    /**
     * Checks whether a solution is valid.
     */
    [[nodiscard]] SCIP_RETCODE check_sol(SCIP *scip,
                                         SCIP_Real obj,
                                         const std::vector<std::size_t> &stack_of,
                                         const std::vector<std::size_t> &level_of) const;

    /**
     * Prints a solution.
     */
    [[nodiscard]] SCIP_RETCODE print_sol(SCIP *scip,
                                         SCIP_SOL *sol,
                                         std::ostream &os = std::cout) const;

    /**
     * Creates the problem.
     */
    [[nodiscard]] static SCIP_RETCODE create_prob(SCIP *scip,
                                                  const char *name,
                                                  Type type,
                                                  std::size_t n_items,
                                                  std::size_t n_stacks,
                                                  std::size_t capacity,
                                                  bool has_adj_conflicts);

    /**
     * Returns the problem data.
     */
    static Probdata &get(SCIP *scip);

    /**
     * Destructor of user problem data to free original user data (called when
     * original problem is freed).
     */
    SCIP_RETCODE scip_delorig(SCIP *scip) override;

    /**
     * Creates user data of transformed problem by transforming the original
     * user problem data (called after problem was transformed).
     */
    SCIP_RETCODE scip_trans(SCIP *scip,
                            scip::ObjProbData **objprobdata,
                            SCIP_Bool *deleteobject) override;

    /**
     * Destructor of user problem data to free transformed user data (called
     * when transformed problem is freed).
     */
    SCIP_RETCODE scip_deltrans(SCIP *scip) override;

    /**
     * Solving process initialization method of transformed data (called before
     * the branch and bound process begins).
     */
    SCIP_RETCODE scip_initsol(SCIP *scip) override;

    /**
     * Solving process deinitialization method of transformed data (called
     * before the branch and bound data is freed).
     */
    SCIP_RETCODE scip_exitsol(SCIP *scip, SCIP_Bool restart) override;

private:
    Type _type;                                  // Type of problem
    std::size_t _n_stacks = 0;                   // Number of stacks
    std::size_t _capacity = 0;                   // Capacity per stack
    bool _has_adjacent_conflicts = false;        // Are conflicts adjacent only?
    std::vector<std::size_t> _arrivals;          // Arrival time of items
    std::vector<double> _weights;                // Weight of items, if available
    std::vector<double> _departures;             // Departure time of items, if available
    std::vector<std::size_t> _initial_positions; // Stack number of initial items
    SquareMatrix<bool> _conflict_matrix;         // Conflict constraints
    SquareMatrix<bool> _blocking_matrix;         // Blocking constraints
    DirectedGraph _conflict_graph;               // Conflict graph
    DirectedGraph _blocking_graph;               // Blocking graph
    DirectedGraph _conflict_graph_complement;    // Conflict graph complement
    Layout _initial_layout;                      // The initial layout
    std::vector<SCIP_CONS *> _conss_assign;      // Assignment constraints
    SCIP_CONS *_cons_space = nullptr;            // Space constraint
    std::vector<SCIP_VAR *> _vars;               // Variables of the master problem
};

} // namespace bp

} // namespace pslp

#endif