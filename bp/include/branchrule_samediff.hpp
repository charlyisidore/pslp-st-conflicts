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

#ifndef PSLP_BP_BRANCHRULE_SAMEDIFF_HPP
#define PSLP_BP_BRANCHRULE_SAMEDIFF_HPP

#include <objscip/objbranchrule.h>

namespace pslp
{

namespace bp
{

/**
 * Same/diff branching rule for the PSLP.
 */
class BranchruleSamediff : public scip::ObjBranchrule
{
public:
    /**
     * Branching rule properties.
     */
    struct Properties
    {
        /**
         * Name of branching rule.
         */
        static const char *Name;

        /**
         * Description of branching rule.
         */
        static const char *Desc;

        /**
         * Priority of branching rule.
         */
        static const int Priority;

        /**
         * Maximal level depth.
         */
        static const int Maxdepth;

        /**
         * Maximal relative distance between dual bounds.
         */
        static const SCIP_Real Maxbounddist;
    };

    /**
     * Default constructor.
     */
    BranchruleSamediff(SCIP *scip);

    /**
     * Branching execution method for fractional LP solutions.
     */
    SCIP_RETCODE scip_execlp(SCIP *scip,
                             SCIP_BRANCHRULE *branchrule,
                             SCIP_Bool allowaddcons,
                             SCIP_RESULT *result) override;
};

} // namespace bp

} // namespace pslp

#endif