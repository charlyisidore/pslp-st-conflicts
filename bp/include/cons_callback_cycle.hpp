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

#ifndef PSLP_BP_CONS_CALLBACK_CYCLE_HPP
#define PSLP_BP_CONS_CALLBACK_CYCLE_HPP

#include <ilcplex/ilocplex.h>

namespace pslp
{

namespace bp
{

/**
 * CPLEX cycle elimination constraint callback.
 */
class ConsCallbackCycle : public IloCplex::LazyConstraintCallbackI
{
public:
    /**
     * Constructor.
     */
    ConsCallbackCycle(IloEnv env,
                      const IloNumVarArray &x,
                      const IloArray<IloBoolArray> &adjacency_matrix,
                      std::size_t max_cuts);

    /**
     * Callback method.
     */
    void main() override;

    /**
     * Create a copy of the invoking callback object.
     */
    IloCplex::CallbackI *duplicateCallback() const override;

    /**
     * Create a callback object.
     */
    static IloCplex::Callback create(IloEnv env,
                                     const IloNumVarArray &x,
                                     const IloArray<IloBoolArray> &adjacency_matrix,
                                     std::size_t max_cuts);

private:
    IloNumVarArray _x;
    IloArray<IloBoolArray> _adjacency_matrix;
    std::size_t _max_cuts = 0;
};

} // namespace bp

} // namespace pslp

#endif