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

#ifndef PSLP_BP_VARDATA_HPP
#define PSLP_BP_VARDATA_HPP

#include <cassert>
#include <objscip/objvardata.h>
#include <vector>

namespace pslp
{

namespace bp
{

/**
 * Variable data for the PSLP.
 */
class Vardata : public scip::ObjVardata
{
public:
    /**
     * Constructor.
     */
    Vardata(const std::vector<std::size_t> &s, SCIP_PRICER *pricer)
        : _items(s),
          _pricer(pricer) {}

    /**
     * Returns the set of items.
     */
    const auto &items() const
    {
        return _items;
    }

    /**
     * Returns the pricer that generates the variable.
     */
    SCIP_PRICER *pricer() const
    {
        return _pricer;
    }

    /**
     * Returns the data associated to a variable.
     */
    static Vardata &get(SCIP *scip, SCIP_VAR *var)
    {
        auto *vardata = static_cast<Vardata *>(SCIPgetObjVardata(scip, var));
        assert(vardata != nullptr);
        return *vardata;
    }

private:
    std::vector<std::size_t> _items;
    SCIP_PRICER *_pricer = nullptr;
};

} // namespace bp

} // namespace pslp

#endif