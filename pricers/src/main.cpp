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

#include "logging.hpp"
#include "pprint.hpp"
#include "pricer.hpp"
#include "pricer_cplex_lazy.hpp"
#include "pricer_cplex_lo.hpp"
#include "pricer_cplex_nf.hpp"
#include "pricer_cplex_pa.hpp"
#include "pricer_greedy.hpp"
#include "pricer_mip_lazy.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using pslp::bp::Pricer;
using pslp::bp::PricerCplexLazy;
using pslp::bp::PricerCplexLo;
using pslp::bp::PricerCplexNf;
using pslp::bp::PricerCplexPa;
using pslp::bp::PricerGreedy;
using pslp::bp::PricerMipLazy;

template <typename... PricersT>
[[nodiscard]] static SCIP_RETCODE benchmark(const Pricer::Problem &prob,
                                            std::size_t n_runs)
{
    using Clock = std::chrono::steady_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    using Duration = std::chrono::duration<double>;

    SCIP *scip;

    SCIP_CALL(SCIPcreate(&scip));

    std::vector<std::unique_ptr<Pricer>> pricers;
    std::vector<std::size_t> active_pricers;
    std::vector<std::string> names;

    (pricers.push_back(std::make_unique<PricersT>(scip)), ...);

    const std::size_t n_pricers = std::size(pricers);

    for (std::size_t i = 0; i < n_pricers; ++i)
    {
        LOG_DEBUG("prepare pricer <%s>", pricers[i]->name());

        if (!pricers[i]->supports(prob))
        {
            LOG_INFO("pricer <%s> disabled (problem not supported)",
                     pricers[i]->name());
            continue;
        }

        active_pricers.push_back(i);
        names.push_back(pricers[i]->name());
    }

    const std::size_t n_active_pricers = std::size(active_pricers);
    std::vector<double> init_times(n_active_pricers, 0.);
    std::vector<double> solve_times(n_active_pricers, 0.);
    std::vector<double> exit_times(n_active_pricers, 0.);
    std::vector<double> total_times(n_active_pricers, 0.);
    std::vector<double> objs(n_active_pricers, 0.);

    LOG_INFO("start benchmark");

    for (std::size_t t = 0; t < n_runs; ++t)
    {
        std::vector<double> init_times_run(n_active_pricers, 0.);
        std::vector<double> solve_times_run(n_active_pricers, 0.);
        std::vector<double> exit_times_run(n_active_pricers, 0.);
        std::vector<double> total_times_run(n_active_pricers, 0.);
        std::vector<double> objs_run(n_active_pricers, 0.);

        for (std::size_t i = 0; i < n_active_pricers; ++i)
        {
            const auto k = active_pricers[i];
            std::vector<Pricer::Solution> sols;
            TimePoint start_point;
            TimePoint init_point;
            TimePoint exit_point;
            TimePoint finish_point;

            assert(pricers[k]->supports(prob));

            start_point = Clock::now();

            SCIP_CALL(pricers[k]->init(scip, prob));

            init_point = Clock::now();

            SCIP_CALL(pricers[k]->solve(scip, prob, sols));

            exit_point = Clock::now();

            SCIP_CALL(pricers[k]->exit(scip));

            finish_point = Clock::now();

            init_times_run[i] = Duration(init_point - start_point).count();
            solve_times_run[i] = Duration(exit_point - init_point).count();
            exit_times_run[i] = Duration(finish_point - exit_point).count();
            total_times_run[i] = Duration(finish_point - start_point).count();

            init_times[i] += init_times_run[i];
            solve_times[i] += solve_times_run[i];
            exit_times[i] += exit_times_run[i];
            total_times[i] += total_times_run[i];

            for (const auto &sol : sols)
            {
                double obj = Pricer::obj_value(prob, sol);
                if (objs_run[i] < obj)
                {
                    objs_run[i] = obj;
                }
            }

            if (objs[i] < objs_run[i])
            {
                objs[i] = objs_run[i];
            }

            LOG_INFO("pricer: <%s>", names[i].c_str());
            LOG_INFO("- obj: %g", objs_run[i]);
            LOG_INFO("- init: %g", init_times_run[i]);
            LOG_INFO("- solve: %g", solve_times_run[i]);
            LOG_INFO("- exit: %g", exit_times_run[i]);
            LOG_INFO("- total: %g", total_times_run[i]);
        }

        double best_obj = 0.;

        for (std::size_t i = 0; i < n_active_pricers; ++i)
        {
            if (best_obj < objs_run[i])
            {
                best_obj = objs_run[i];
            }
        }

        for (std::size_t i = 0; i < n_active_pricers; ++i)
        {
            const auto k = active_pricers[i];
            if (pricers[k]->is_exact() &&
                !SCIPisEQ(scip, objs_run[i], best_obj))
            {
                LOG_ERROR("multiple exact values");
                return SCIP_ERROR;
            }
            else if (!pricers[k]->is_exact() &&
                     SCIPisGT(scip, objs_run[i], best_obj))
            {
                LOG_ERROR("heuristic value better than exact value");
                return SCIP_ERROR;
            }
        }

        if (n_runs > 1)
        {
            std::cout << "====== run: " << t + 1 << " ======" << std::endl;
            std::cout << pprint_table({"pricer",
                                       "obj",
                                       "init",
                                       "solve",
                                       "exit",
                                       "total"},
                                      names,
                                      objs_run,
                                      init_times_run,
                                      solve_times_run,
                                      exit_times_run,
                                      total_times_run);
        }
    }

    SCIP_CALL(SCIPfree(&scip));

    BMScheckEmptyMemory();

    std::cout << "====== summary ======" << std::endl;
    std::cout << pprint_table({"pricer", "obj", "init", "solve", "exit", "total"},
                              names,
                              objs,
                              init_times,
                              solve_times,
                              exit_times,
                              total_times);

    return SCIP_OKAY;
}

[[nodiscard]] static SCIP_RETCODE run(int argc, char **argv)
{
    if (argc <= 1)
    {
        std::cout << "Usage: " << argv[0] << " FILE" << std::endl;
        return SCIP_OKAY;
    }

    const char *filename = argv[1];
    std::size_t n_runs = 1;

    if (argc > 2)
    {
        std::istringstream(argv[2]) >> n_runs;
    }

    // Open file

    LOG_INFO("read pricing problem in <%s>", filename);

    std::ifstream is(filename);

    if (!is.is_open())
    {
        SCIPerrorMessage("cannot open <%s>\n", filename);
        return SCIP_ERROR;
    }

    // Read instance

    std::string typestr;
    Pricer::Problem::Type type;
    std::size_t n_items;
    std::size_t capacity;
    bool has_adjacent_conflicts;
    bool farkas;
    std::vector<double> coefs;
    double offset;
    std::size_t n_same;
    std::vector<std::pair<std::size_t, std::size_t>> same;
    std::size_t n_diff;
    std::vector<std::pair<std::size_t, std::size_t>> diff;
    std::size_t n_conflict_edges;
    DirectedGraph conflict_graph;
    SquareMatrix<bool> conflict_matrix;
    std::size_t n_blocking_edges;
    DirectedGraph blocking_graph;
    SquareMatrix<bool> blocking_matrix;

    is >> typestr;

    if (typestr == "St")
    {
        type = Pricer::Problem::Type::MinNStacks;
    }
    else if (typestr == "BI")
    {
        type = Pricer::Problem::Type::MinNBlockingItems;
    }
    else
    {
        std::cerr << "problem type '"
                  << typestr << "' not supported"
                  << std::endl;
        return SCIP_ERROR;
    }

    is >> n_items >> capacity >> has_adjacent_conflicts >> farkas;

    coefs.resize(n_items);

    for (std::size_t i = 0; i < n_items; ++i)
    {
        is >> coefs[i];
    }

    is >> offset;

    is >> n_same;

    same.reserve(n_same);

    for (std::size_t e = 0; e < n_same; ++e)
    {
        std::size_t i, j;
        is >> i >> j;
        --i;
        --j;
        same.push_back({i, j});
    }

    is >> n_diff;

    diff.reserve(n_diff);

    for (std::size_t e = 0; e < n_diff; ++e)
    {
        std::size_t i, j;
        is >> i >> j;
        --i;
        --j;
        diff.push_back({i, j});
    }

    is >> n_conflict_edges;

    conflict_graph.initialize(n_items);
    conflict_matrix.initialize(n_items);

    for (std::size_t e = 0; e < n_conflict_edges; ++e)
    {
        std::size_t i, j;
        is >> i >> j;
        --i;
        --j;
        conflict_graph.add_edge(i, j);
        conflict_matrix(i, j) = true;
    }

    is >> n_blocking_edges;

    blocking_graph.initialize(n_items);
    blocking_matrix.initialize(n_items);

    for (std::size_t e = 0; e < n_blocking_edges; ++e)
    {
        std::size_t i, j;
        is >> i >> j;
        --i;
        --j;
        blocking_graph.add_edge(i, j);
        blocking_matrix(i, j) = true;
    }

    Pricer::Problem prob{type,
                         capacity,
                         has_adjacent_conflicts,
                         conflict_graph,
                         conflict_matrix,
                         blocking_graph,
                         blocking_matrix,
                         same,
                         diff,
                         coefs,
                         offset,
                         farkas,
                         true};

    // Benchmark

    SCIP_CALL((benchmark<
               PricerGreedy,
               PricerCplexLo,
               PricerCplexNf,
               PricerCplexPa,
               PricerCplexLazy,
               PricerMipLazy>(prob, n_runs)));

    return SCIP_OKAY;
}

int main(int argc, char **argv)
{
    SCIP_RETCODE retcode = run(argc, argv);

    if (retcode != SCIP_OKAY)
    {
        SCIPprintError(retcode);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}