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

#ifndef PSLP_BP_HEUR_INIT_HPP
#define PSLP_BP_HEUR_INIT_HPP

#include <objscip/objheur.h>

namespace pslp
{

namespace bp
{

/**
 * Initial primal heuristic for the PSLP.
 */
class HeurInit : public scip::ObjHeur
{
public:
    /**
     * Primal heuristic properties.
     */
    struct Properties
    {
        /**
         * Name of primal heuristic.
         */
        static const char *Name;

        /**
         * Description of primal heuristic.
         */
        static const char *Desc;

        /**
         * Display character of primal heuristic.
         */
        static const char Dispchar;

        /**
         * Priority of primal heuristic.
         */
        static const int Priority;

        /**
         * Frequency for calling primal heuristic.
         */
        static const int Freq;

        /**
         * Frequency offset for calling primal heuristic.
         */
        static const int Freqofs;

        /**
         * Maximal depth level to call heuristic at (-1: no limit).
         */
        static const int Maxdepth;

        /**
         * Positions in the node solving loop where heuristic should be executed.
         */
        static const SCIP_HEURTIMING Timing;

        /**
         * Does the heuristic use a secondary SCIP instance?
         */
        static const SCIP_Bool Usesubscip;
    };

    /**
     * Default constructor.
     */
    HeurInit(SCIP *scip);

    /**
     * Execution method of primal heuristic.
     */
    SCIP_RETCODE scip_exec(SCIP *scip,
                           SCIP_HEUR *heur,
                           SCIP_HEURTIMING heurtiming,
                           SCIP_Bool nodeinfeasible,
                           SCIP_RESULT *result) override;
};

} // namespace bp

} // namespace pslp

#endif