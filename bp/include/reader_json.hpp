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

#ifndef PSLP_BP_READER_JSON_HPP
#define PSLP_BP_READER_JSON_HPP

#include <objscip/objreader.h>

namespace pslp
{

namespace bp
{

/**
 * Reader for PSLP instances in JSON format.
 */
class ReaderJson : public scip::ObjReader
{
public:
    /**
     * File reader properties.
     */
    struct Properties
    {
        /**
         * Name of file reader.
         */
        static const char *Name;

        /**
         * Description of file reader.
         */
        static const char *Desc;

        /**
         * File extension that reader processes.
         */
        static const char *Extension;
    };

    /**
     * Default constructor.
     */
    ReaderJson(SCIP *scip);

    /**
     * Problem reading method of reader.
     */
    SCIP_RETCODE scip_read(SCIP *scip,
                           SCIP_READER *reader,
                           const char *filename,
                           SCIP_RESULT *result) override;

    /**
     * Problem writing method of reader.
     */
    SCIP_RETCODE scip_write(SCIP *scip,
                            SCIP_READER *reader,
                            FILE *file,
                            const char *name,
                            SCIP_PROBDATA *probdata,
                            SCIP_Bool transformed,
                            SCIP_OBJSENSE objsense,
                            SCIP_Real objscale,
                            SCIP_Real objoffset,
                            SCIP_VAR **vars,
                            int nvars,
                            int nbinvars,
                            int nintvars,
                            int nimplvars,
                            int ncontvars,
                            SCIP_VAR **fixedvars,
                            int nfixedvars,
                            int startnvars,
                            SCIP_CONS **conss,
                            int nconss,
                            int maxnconss,
                            int startnconss,
                            SCIP_Bool genericnames,
                            SCIP_RESULT *result) override;

private:
    SCIP_Bool _useadjacentconflicts = true;
    SCIP_Bool _transformconflicts = false;
};

} // namespace bp

} // namespace pslp

#endif