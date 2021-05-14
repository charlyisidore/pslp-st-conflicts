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

#ifndef PSLP_BP_PRICER_HPP
#define PSLP_BP_PRICER_HPP

#include "graph.hpp"
#include "probdata.hpp"
#include "square_matrix.hpp"
#include <objscip/objpricer.h>
#include <string>
#include <vector>

namespace pslp
{

namespace bp
{

/**
 * Base class for variable pricers.
 */
class Pricer : public scip::ObjPricer
{
public:
    /**
     * Pricer info type definition.
     */
    struct Info
    {
        /**
         * Name of the pricer.
         */
        const char *name = nullptr;

        /**
         * Factory for creating pricer objects.
         */
        scip::ObjPricer *(*create)(SCIP *);

        /**
         * Determine whether the pricer should be activated by default.
         */
        bool activate = false;
    };

    /**
     * Pricer problem data type definition.
     */
    struct Problem
    {
        /**
         * Problem type definition.
         */
        using Type = Probdata::Type;

        /**
         * Problem type.
         */
        Type type;

        /**
         * Stack capacity.
         */
        std::size_t capacity = 0;

        /**
         * Determine whether conflicts are adjacent only.
         */
        bool has_adjacent_conflicts = false;

        /**
         * Conflict graph.
         */
        const DirectedGraph &conflict_graph;

        /**
         * Conflict matrix.
         */
        const SquareMatrix<bool> &conflict_matrix;

        /**
         * Blocking graph.
         */
        const DirectedGraph &blocking_graph;

        /**
         * Blocking matrix.
         */
        const SquareMatrix<bool> &blocking_matrix;

        /**
         * Constraints of type "same".
         */
        const std::vector<std::pair<std::size_t, std::size_t>> &same;

        /**
         * Constraints of type "diff".
         */
        const std::vector<std::pair<std::size_t, std::size_t>> &diff;

        /**
         * Coefficients for items.
         */
        const std::vector<double> &coefs;

        /**
         * Objective offset.
         */
        double offset = 0.;

        /**
         * Determine whether the pricing is in farkas mode.
         */
        bool farkas = false;

        /**
         * Inform that the constraints have changed.
         */
        bool has_changed = true;
    };

    /**
     * Solution type definition.
     */
    using Solution = std::vector<std::size_t>;

    /**
     * Constructor.
     */
    Pricer(SCIP *scip,
           const char *name,
           const char *desc,
           int priority,
           SCIP_Bool delay);

    /**
     * Destructor.
     */
    virtual ~Pricer() {}

    /**
     * Return the name of the pricer.
     */
    const char *name() const
    {
        return _name.c_str();
    }

    /**
     * Check whether the pricer finds exact solutions.
     */
    virtual bool is_exact() const = 0;

    /**
     * Check whether the pricer supports a given problem.
     */
    virtual bool supports([[maybe_unused]] const Problem &prob) const
    {
        return true;
    }

    /**
     * Initialize the pricer.
     */
    [[nodiscard]] virtual SCIP_RETCODE init(SCIP *scip, const Problem &prob);

    /**
     * Deinitialize the pricer.
     */
    [[nodiscard]] virtual SCIP_RETCODE exit(SCIP *scip);

    /**
     * Solve the pricing problem.
     */
    [[nodiscard]] virtual SCIP_RETCODE solve(SCIP *scip,
                                             const Problem &prob,
                                             std::vector<Solution> &sols) = 0;

    /**
     * Pricing method.
     */
    [[nodiscard]] SCIP_RETCODE pricing(SCIP *scip,
                                       SCIP_PRICER *pricer,
                                       bool farkas,
                                       SCIP_RESULT *result);

    /**
     * Replace '{}' by the call number in a filename.
     */
    std::string make_filename(SCIP *scip, std::string filename);

    /**
     * Returns the objective value of a solution.
     */
    static double obj_value(const Problem &prob, const Solution &sol);

    /**
     * Returns a conflict subgraph built from a subset of items.
     */
    static DirectedGraph subgraph(const SquareMatrix<bool> &matrix,
                                  const std::vector<std::size_t> &s);

    /**
     * Returns a conflict subgraph complement built from a subset of items.
     */
    static DirectedGraph subgraph_complement(const SquareMatrix<bool> &matrix,
                                             const std::vector<std::size_t> &s);

    /**
     * Outputs the pricing problem data in a file.
     */
    static void write_prob(const Problem &prob, const char *filename);

    /**
     * Returns the list of registered pricers.
     */
    static std::vector<Info> &list();

    /**
     * Registers a pricer.
     */
    template <typename PricerT>
    static void register_pricer(bool activate)
    {
        Info info;
        info.name = PricerT::Properties::Name;
        info.create = _create<PricerT>;
        info.activate = activate;
        list().push_back(info);
    }

    /**
     * Initialization method of variable pricer (called after problem was
     * transformed).
     */
    SCIP_RETCODE scip_init(SCIP *scip, SCIP_PRICER *pricer) override;

    /**
     * Deinitialization method of variable pricer (called before transformed
     * problem is freed).
     */
    SCIP_RETCODE scip_exit(SCIP *scip, SCIP_PRICER *pricer) override;

    /**
     * Reduced cost pricing method of variable pricer for feasible LPs.
     */
    SCIP_RETCODE scip_redcost(SCIP *scip,
                              SCIP_PRICER *pricer,
                              SCIP_Real *lowerbound,
                              SCIP_Bool *stopearly,
                              SCIP_RESULT *result) override;

    /**
     * Farkas pricing method of variable pricer for infeasible LPs.
     */
    SCIP_RETCODE scip_farkas(SCIP *scip,
                             SCIP_PRICER *pricer,
                             SCIP_RESULT *result) override;

private:
    /**
     * Creates a pricer object.
     */
    template <typename PricerT>
    static scip::ObjPricer *_create(SCIP *scip)
    {
        return new PricerT(scip);
    }

    SCIP_NODE *_node = nullptr;
    std::vector<std::pair<std::size_t, std::size_t>> _same;
    std::vector<std::pair<std::size_t, std::size_t>> _diff;
    std::string _name;
    char *_probfilename = nullptr;
};

} // namespace bp

} // namespace pslp

#endif