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

#ifndef PSLP_BP_PRICER_MULTI_HPP
#define PSLP_BP_PRICER_MULTI_HPP

#include "logging.hpp"
#include "pprint.hpp"
#include "pricer.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <scip/scip.h>
#include <sstream>
#include <vector>

namespace pslp
{

namespace bp
{

/**
 * Executes multiple pricers for comparison.
 */
template <typename MainPricerT, typename... PricersT>
class PricerMulti : public Pricer
{
public:
    /**
     * Variable pricer properties.
     */
    struct Properties
    {
        /**
         * Name of variable pricer.
         */
        static constexpr const char *Name = "multi";

        /**
         * Description of variable pricer.
         */
        static constexpr const char *Desc = "executes multiple pricers for comparison";

        /**
         * Priority of variable pricer.
         */
        static constexpr const int Priority = 0;

        /**
         * Delay the variable pricer?
         */
        static constexpr const SCIP_Bool Delay = true;
    };

    /**
     * Pricer problem type definition.
     */
    using Problem = Pricer::Problem;

    /**
     * Pricer solution type definition.
     */
    using Solution = Pricer::Solution;

    /**
     * Clock type definition.
     */
    using Clock = std::chrono::steady_clock;

    /**
     * Time point type definition.
     */
    using TimePoint = std::chrono::time_point<Clock>;

    /**
     * Duration type definition.
     */
    using Duration = std::chrono::duration<double>;

    /**
     * Constructor.
     */
    PricerMulti(SCIP *scip)
        : Pricer(scip,
                 Properties::Name,
                 Properties::Desc,
                 Properties::Priority,
                 Properties::Delay)
    {
        SCIPcreate(&_subscip);

        _pricers.push_back(std::make_unique<MainPricerT>(_subscip));
        _names.push_back(MainPricerT::Properties::Name);
        _priority.push_back(MainPricerT::Properties::Priority);

        (_pricers.push_back(std::make_unique<PricersT>(_subscip)), ...);
        (_names.push_back(PricersT::Properties::Name), ...);
        (_priority.push_back(PricersT::Properties::Priority), ...);

        const auto n_pricers = std::size(_pricers);

        _use.resize(n_pricers, true);

        for (std::size_t i = 0; i < n_pricers; ++i)
        {
            {
                std::ostringstream oss_name;
                std::ostringstream oss_desc;
                oss_name << "pricers/" << Properties::Name << '/' << _names[i]
                         << "/use";
                oss_desc << "(pslp) use <" << _names[i] << "> pricer?";
                SCIPaddBoolParam(scip,
                                 oss_name.str().c_str(),
                                 oss_desc.str().c_str(),
                                 &_use[i],
                                 false,
                                 _use[i],
                                 nullptr,
                                 nullptr);
            }
            {
                std::ostringstream oss_name;
                std::ostringstream oss_desc;
                oss_name << "pricers/" << Properties::Name << '/' << _names[i]
                         << "/priority";
                oss_desc << "(pslp) priority of <" << _names[i] << "> pricer";
                SCIPaddIntParam(scip,
                                oss_name.str().c_str(),
                                oss_desc.str().c_str(),
                                &_priority[i],
                                false,
                                _priority[i],
                                INT_MIN / 4,
                                INT_MAX / 4,
                                nullptr,
                                nullptr);
            }
        }

        {
            std::ostringstream oss;
            oss << "pricers/" << Properties::Name << "/addvars";
            SCIPaddBoolParam(scip,
                             oss.str().c_str(),
                             "(pslp) should the pricer add variables?",
                             &_addvars,
                             false,
                             _addvars,
                             nullptr,
                             nullptr);
        }
    }

    /**
     * Destructor.
     */
    ~PricerMulti()
    {
        SCIPfree(&_subscip);
    }

    /**
     * Check whether the pricer finds exact solutions.
     */
    bool is_exact() const override
    {
        const auto n_pricers = std::size(_pricers);

        for (std::size_t i = 0; i < n_pricers; ++i)
        {
            if (_use[i] && _pricers[i]->is_exact())
            {
                return true;
            }
        }

        return false;
    }

    /**
     * Check whether the pricer supports a given problem.
     */
    bool supports(const Problem &prob) const override
    {
        const auto n_pricers = std::size(_pricers);

        for (std::size_t i = 0; i < n_pricers; ++i)
        {
            if (_use[i] && _pricers[i]->supports(prob))
            {
                return true;
            }
        }

        return false;
    }

    /**
     * Initialize the pricers.
     */
    [[nodiscard]] SCIP_RETCODE init([[maybe_unused]] SCIP *scip,
                                    const Problem &prob) override
    {
        const auto n_pricers = std::size(_pricers);

        _order.clear();
        for (std::size_t i = 0; i < n_pricers; ++i)
        {
            if (_use[i] && _pricers[i]->supports(prob))
            {
                _order.push_back(i);
            }
        }

        std::stable_sort(std::begin(_order),
                         std::end(_order),
                         [this](auto i, auto j) {
                             return _priority[i] > _priority[j];
                         });

        const auto n_active_pricers = std::size(_order);
        std::vector<int> priorities(n_active_pricers);

        _active_names.resize(n_active_pricers);

        for (std::size_t k = 0; k < n_active_pricers; ++k)
        {
            const auto i = _order[k];

            _active_names[k] = _names[i];
            priorities[k] = _priority[i];
        }

        std::cout << "====== PricerMulti init ======" << std::endl;
        std::cout << pprint_table({"pricer", "priority"},
                                  _active_names,
                                  priorities);

        _n_total_sols.resize(n_active_pricers, 0);
        _init_times.resize(n_active_pricers, 0.);
        _solve_times.resize(n_active_pricers, 0.);
        _exit_times.resize(n_active_pricers, 0.);

        for (std::size_t k = 0; k < n_active_pricers; ++k)
        {
            const auto i = _order[k];
            TimePoint start_point;
            TimePoint finish_point;

            start_point = Clock::now();
            SCIP_CALL(_pricers[i]->init(_subscip, prob));
            finish_point = Clock::now();

            _init_times[k] += Duration(finish_point - start_point).count();
        }

        return SCIP_OKAY;
    }

    /**
     * Deinitialize the pricers.
     */
    [[nodiscard]] SCIP_RETCODE exit([[maybe_unused]] SCIP *scip) override
    {
        const auto n_active_pricers = std::size(_order);
        std::vector<double> total_times(n_active_pricers, 0.);

        for (std::size_t k = 0; k < n_active_pricers; ++k)
        {
            const auto i = _order[k];
            TimePoint start_point;
            TimePoint finish_point;

            start_point = Clock::now();
            SCIP_CALL(_pricers[i]->exit(_subscip));
            finish_point = Clock::now();

            _exit_times[k] += Duration(finish_point - start_point).count();
            total_times[k] = _init_times[k] + _solve_times[k] + _exit_times[k];
        }

        std::cout << "====== PricerMulti exit ======" << std::endl;
        std::cout << pprint_table({"pricer",
                                   "n_sols",
                                   "init",
                                   "solve",
                                   "exit",
                                   "total"},
                                  _active_names,
                                  _n_total_sols,
                                  _init_times,
                                  _solve_times,
                                  _exit_times,
                                  total_times);

        return SCIP_OKAY;
    }

    /**
     * Solve the pricing problem.
     */
    [[nodiscard]] SCIP_RETCODE solve([[maybe_unused]] SCIP *scip,
                                     const Problem &prob,
                                     std::vector<Solution> &sols) override
    {

        const auto n_active_pricers = std::size(_order);
        std::vector<double> objs(n_active_pricers, 0.);
        std::vector<double> times(n_active_pricers, 0.);
        std::vector<std::size_t> n_sols(n_active_pricers, 0);

        for (std::size_t k = 0; k < n_active_pricers; ++k)
        {
            const auto i = _order[k];
            std::vector<Solution> subsols;
            TimePoint start_point;
            TimePoint finish_point;

            start_point = Clock::now();
            SCIP_CALL(_pricers[i]->solve(_subscip, prob, subsols));
            finish_point = Clock::now();

            for (const auto &sol : subsols)
            {
                double obj = obj_value(prob, sol);
                if (objs[k] < obj)
                {
                    objs[k] = obj;
                }
            }

            if (_addvars && std::empty(sols))
            {
                sols = subsols;
            }

            times[k] = Duration(finish_point - start_point).count();
            _solve_times[k] += times[k];

            n_sols[k] = std::size(subsols);
            _n_total_sols[k] += n_sols[k];
        }

        std::cout << pprint_table({"pricer", "n_sols", "obj", "time"},
                                  _active_names, n_sols, objs, times);

        double best_obj = 0.;

        for (std::size_t k = 0; k < n_active_pricers; ++k)
        {
            if (best_obj < objs[k])
            {
                best_obj = objs[k];
            }
        }

        for (std::size_t k = 0; k < n_active_pricers; ++k)
        {
            const auto i = _order[k];

            if (_pricers[i]->is_exact() &&
                !SCIPisEQ(_subscip, objs[k], best_obj))
            {
                LOG_ERROR("multiple exact values");
                return SCIP_ERROR;
            }
            else if (!_pricers[i]->is_exact() &&
                     SCIPisGT(_subscip, objs[k], best_obj))
            {
                LOG_ERROR("heuristic value better than exact value");
                return SCIP_ERROR;
            }
        }

        return SCIP_OKAY;
    }

private:
    SCIP *_subscip = nullptr;
    std::vector<std::unique_ptr<Pricer>> _pricers;
    std::vector<const char *> _names;
    std::vector<std::size_t> _order;
    std::vector<const char *> _active_names;
    std::vector<SCIP_Bool> _use;
    std::vector<int> _priority;
    SCIP_Bool _addvars = false;
    std::vector<std::size_t> _n_total_sols;
    std::vector<double> _init_times;
    std::vector<double> _solve_times;
    std::vector<double> _exit_times;
};

} // namespace bp

} // namespace pslp

#endif