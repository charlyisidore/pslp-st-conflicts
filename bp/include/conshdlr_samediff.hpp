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

#ifndef PSLP_BP_CONSHDLR_SAMEDIFF_HPP
#define PSLP_BP_CONSHDLR_SAMEDIFF_HPP

#include <objscip/objconshdlr.h>

namespace pslp
{

namespace bp
{

/**
 * Same/diff constraint handler for the PSLP.
 */
class ConshdlrSamediff : public scip::ObjConshdlr
{
public:
    /**
     * Constraint type definition.
     */
    enum class Type
    {
        Same,
        Diff
    };

    /**
     * Constraint data type definition.
     */
    struct Consdata
    {
        Type type;
        std::size_t i;
        std::size_t j;
    };

    /**
     * Constraint handler properties.
     */
    struct Properties
    {
        /**
         * Name of constraint handler.
         */
        static const char *Name;

        /**
         * Description of constraint handler.
         */
        static const char *Desc;

        /**
         * Priority of constraint handler for separation.
         */
        static const int Sepapriority;

        /**
         * Priority of constraint handler for constraint enforcing.
         */
        static const int Enfopriority;

        /**
         * Priority of constraint handler for checking infeasibility (and
         * propagation).
         */
        static const int Checkpriority;

        /**
         * Frequency for separating cuts.
         */
        static const int Sepafreq;

        /**
         * Frequency for propagating domains.
         */
        static const int Propfreq;

        /**
         * Frequency for using all instead of only the useful constraints in
         * separation, propagation and enforcement.
         */
        static const int Eagerfreq;

        /**
         * Maximal number of presolving rounds the constraint handler participates
         * in.
         */
        static const int Maxprerounds;

        /**
         * Should separation method be delayed, if other separators found cuts?
         */
        static const SCIP_Bool Delaysepa;

        /**
         * Should propagation method be delayed, if other propagators found
         * reductions?
         */
        static const SCIP_Bool Delayprop;

        /**
         * Should the constraint handler be skipped, if no constraints are
         * available?
         */
        static const SCIP_Bool Needscons;

        /**
         * Positions in the node solving loop where propagation method of constraint
         * handlers should be executed.
         */
        static const SCIP_PROPTIMING Proptiming;

        /**
         * Timing mask of the constraint handler's presolving method.
         */
        static const SCIP_PRESOLTIMING Presoltiming;
    };

    /**
     * Default constructor.
     */
    ConshdlrSamediff(SCIP *scip);

    /**
     * Creates and captures a constraint.
     */
    [[nodiscard]] static SCIP_RETCODE create_cons(SCIP *scip,
                                                  SCIP_CONS **cons,
                                                  const char *name,
                                                  Type type,
                                                  std::size_t i,
                                                  std::size_t j,
                                                  SCIP_NODE *node);

    /**
     * Returns the data of a constraint.
     */
    static Consdata data(SCIP_CONS *cons);

    /**
     * Check whether a constraint is satisfied.
     */
    bool check(SCIP *scip, SCIP_CONS *cons, bool is_propagated);

    /**
     * Frees specific constraint data.
     */
    SCIP_RETCODE scip_delete(SCIP *scip,
                             SCIP_CONSHDLR *conshdlr,
                             SCIP_CONS *cons,
                             SCIP_CONSDATA **consdata) override;

    /**
     * Transforms constraint data into data belonging to the transformed
     * problem.
     */
    SCIP_RETCODE scip_trans(SCIP *scip,
                            SCIP_CONSHDLR *conshdlr,
                            SCIP_CONS *sourcecons,
                            SCIP_CONS **targetcons) override;

    /**
     * Constraint enforcing method of constraint handler for LP solutions.
     */
    SCIP_RETCODE scip_enfolp(SCIP *scip,
                             SCIP_CONSHDLR *conshdlr,
                             SCIP_CONS **conss,
                             int nconss,
                             int nusefulconss,
                             SCIP_Bool solinfeasible,
                             SCIP_RESULT *result) override;

    /**
     * Constraint enforcing method of constraint handler for pseudo solutions.
     */
    SCIP_RETCODE scip_enfops(SCIP *scip,
                             SCIP_CONSHDLR *conshdlr,
                             SCIP_CONS **conss,
                             int nconss,
                             int nusefulconss,
                             SCIP_Bool solinfeasible,
                             SCIP_Bool objinfeasible,
                             SCIP_RESULT *result) override;

    /**
     * Feasibility check method of constraint handler for primal solutions.
     */
    SCIP_RETCODE scip_check(SCIP *scip,
                            SCIP_CONSHDLR *conshdlr,
                            SCIP_CONS **conss,
                            int nconss,
                            SCIP_SOL *sol,
                            SCIP_Bool checkintegrality,
                            SCIP_Bool checklprows,
                            SCIP_Bool printreason,
                            SCIP_Bool completely,
                            SCIP_RESULT *result) override;

    /**
     * Domain propagation method of constraint handler.
     */
    SCIP_RETCODE scip_prop(SCIP *scip,
                           SCIP_CONSHDLR *conshdlr,
                           SCIP_CONS **conss,
                           int nconss,
                           int nusefulconss,
                           int nmarkedconss,
                           SCIP_PROPTIMING proptiming,
                           SCIP_RESULT *result) override;

    /**
     * Variable rounding lock method of constraint handler.
     */
    SCIP_RETCODE scip_lock(SCIP *scip,
                           SCIP_CONSHDLR *conshdlr,
                           SCIP_CONS *cons,
                           SCIP_LOCKTYPE locktype,
                           int nlockspos,
                           int nlocksneg) override;

    /**
     * Constraint activation notification method of constraint handler.
     */
    SCIP_RETCODE scip_active(SCIP *scip,
                             SCIP_CONSHDLR *conshdlr,
                             SCIP_CONS *cons) override;

    /**
     * Constraint deactivation notification method of constraint handler.
     */
    SCIP_RETCODE scip_deactive(SCIP *scip,
                               SCIP_CONSHDLR *conshdlr,
                               SCIP_CONS *cons) override;
};

} // namespace bp

} // namespace pslp

#endif