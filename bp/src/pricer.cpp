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

#include "pricer.hpp"
#include "conshdlr_samediff.hpp"
#include "logging.hpp"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <numeric>
#include <scip/scip.h>
#include <set>

using namespace pslp::bp;

std::vector<Pricer::Info> &Pricer::list()
{
    static std::vector<Info> list;
    return list;
}

Pricer::Pricer(SCIP *scip,
               const char *name,
               const char *desc,
               int priority,
               SCIP_Bool delay)
    : scip::ObjPricer(scip, name, desc, priority, delay),
      _name(name)
{
    {
        std::ostringstream oss;
        oss << "pricers/" << name << "/probfilename";
        SCIPaddStringParam(scip,
                           oss.str().c_str(),
                           "(pslp) file name to export pricer instances",
                           &_probfilename,
                           false,
                           "",
                           nullptr,
                           nullptr);
    }
}

SCIP_RETCODE Pricer::init([[maybe_unused]] SCIP *scip,
                          [[maybe_unused]] const Problem &prob)
{
    return SCIP_OKAY;
}

SCIP_RETCODE Pricer::exit([[maybe_unused]] SCIP *scip)
{
    return SCIP_OKAY;
}

SCIP_RETCODE Pricer::pricing(SCIP *scip,
                             SCIP_PRICER *pricer,
                             bool farkas,
                             SCIP_RESULT *result)
{
    assert(result != nullptr);

    // Get the data to create an instance of the pricing problem
    auto &probdata = Probdata::get(scip);
    const auto n_items = probdata.n_items();
    std::vector<double> coefs(n_items, 0.);
    double offset = 0.;

    switch (probdata.type())
    {
        case Probdata::Type::MinNStacks:
            offset = farkas ? 0. : -1.;
            break;

        case Probdata::Type::MinNBlockingItems:
            offset = probdata.dual_space(scip, farkas);
            break;
    }

    for (std::size_t i = 0; i < n_items; ++i)
    {
        coefs[i] = probdata.dual_assign(scip, i, farkas);
    }

    SCIP_NODE *current = SCIPgetCurrentNode(scip);
    bool has_changed = false;

    if (_node == nullptr || _node != current)
    {
        _node = current;
        has_changed = true;

        SCIP_CONSHDLR *conshdlr = SCIPfindConshdlr(scip,
                                                   ConshdlrSamediff::Properties::Name);
        SCIP_CONS **conss = SCIPconshdlrGetConss(conshdlr);
        int n_active_conss = SCIPconshdlrGetNActiveConss(conshdlr);

        _same.clear();
        _diff.clear();
        _same.reserve(n_active_conss);
        _diff.reserve(n_active_conss);

        for (int c = 0; c < n_active_conss; ++c)
        {
            SCIP_CONS *cons = conss[c];

            assert(SCIPconsIsActive(cons));

            const auto [type, i, j] = ConshdlrSamediff::data(cons);

            if (type == ConshdlrSamediff::Type::Same)
            {
                _same.push_back({i, j});
            }
            else // type == ConshdlrSamediff::Type::Diff
            {
                _diff.push_back({i, j});
            }
        }

#ifndef NDEBUG
        // The constraints after n_active_conss should be inactive
        int n_conss = SCIPconshdlrGetNConss(conshdlr);
        for (int c = n_active_conss; c < n_conss; ++c)
        {
            assert(!SCIPconsIsActive(conss[c]));
        }
#endif
    }

    // Create an instance of the pricing problem
    Problem prob{probdata.type(),
                 probdata.capacity(),
                 probdata.has_adjacent_conflicts(),
                 probdata.conflict_graph(),
                 probdata.conflict_matrix(),
                 probdata.blocking_graph(),
                 probdata.blocking_matrix(),
                 _same,
                 _diff,
                 coefs,
                 offset,
                 farkas,
                 has_changed};

    if (_probfilename && *_probfilename)
    {
        std::string filename = make_filename(scip, _probfilename);
        LOG_DEBUG("write pricing problem in <%s>", filename.c_str());
        write_prob(prob, filename.c_str());
    }

    std::vector<Solution> sols;

    // This is a virtual method
    SCIP_CALL(solve(scip, prob, sols));

    // Add each of the solutions to the master problem
    for (const auto &sol : sols)
    {
#ifndef NDEBUG
        // Check if the solution is improving
        double obj = obj_value(prob, sol);
        assert(SCIPisPositive(scip, obj));

        // Check if the same constraints are satisfied
        for (auto [i, j] : _same)
        {
            bool has_i = std::find(std::begin(sol),
                                   std::end(sol),
                                   i) != std::end(sol);
            bool has_j = std::find(std::begin(sol),
                                   std::end(sol),
                                   j) != std::end(sol);

            assert(has_i == has_j);
        }

        // Check if the diff constraints are satisfied
        for (auto [i, j] : _diff)
        {
            bool has_i = std::find(std::begin(sol),
                                   std::end(sol),
                                   i) != std::end(sol);
            bool has_j = std::find(std::begin(sol),
                                   std::end(sol),
                                   j) != std::end(sol);

            assert(!has_i || !has_j);
        }
#endif

        SCIP_CALL(probdata.add_var(scip, sol, pricer));
    }

    *result = SCIP_SUCCESS;

    return SCIP_OKAY;
}

std::string Pricer::make_filename(SCIP *scip, std::string filename)
{
    std::string placeholder = "{}";
    std::size_t pos = filename.find(placeholder);

    if (pos != std::string::npos)
    {
        std::ostringstream oss;
        std::size_t n = 1;

        if (SCIPgetStage(scip) == SCIP_STAGE_SOLVING)
        {
            auto *pricer = SCIPfindPricer(scip, _name.c_str());
            if (pricer)
            {
                n = 1 + SCIPpricerGetNCalls(pricer);
            }
        }

        oss << n;

        filename.replace(pos, std::size(placeholder), oss.str());
    }

    return filename;
}

double Pricer::obj_value(const Problem &prob, const Solution &sol)
{
    double obj = prob.offset;

    for (auto i : sol)
    {
        obj += prob.coefs[i];
    }

    if (prob.type == Pricer::Problem::Type::MinNBlockingItems)
    {
        std::size_t n_blocking = 0;

        for (auto i = std::begin(sol); i != std::end(sol); ++i)
        {
            for (auto j = std::begin(sol); j != i; ++j)
            {
                if (prob.blocking_matrix(*i, *j))
                {
                    ++n_blocking;
                    break;
                }
            }
        }

        obj -= double(n_blocking);
    }

    return obj;
}

DirectedGraph Pricer::subgraph(const SquareMatrix<bool> &matrix,
                               const std::vector<std::size_t> &s)
{
    const auto n = std::size(s);
    DirectedGraph graph(n);

    for (std::size_t u = 0; u < n; ++u)
    {
        const auto i = s[u];

        for (std::size_t v = 0; v < n; ++v)
        {
            if (u == v)
            {
                continue;
            }

            const auto j = s[v];

            if (matrix(i, j))
            {
                graph.add_edge(u, v);
            }
        }
    }

    return graph;
}

DirectedGraph Pricer::subgraph_complement(const SquareMatrix<bool> &matrix,
                                          const std::vector<std::size_t> &s)
{
    const auto n = std::size(s);
    DirectedGraph graph(n);

    for (std::size_t u = 0; u < n; ++u)
    {
        const auto i = s[u];

        for (std::size_t v = 0; v < n; ++v)
        {
            if (u == v)
            {
                continue;
            }

            const auto j = s[v];

            if (!matrix(i, j))
            {
                graph.add_edge(u, v);
            }
        }
    }

    return graph;
}

void Pricer::write_prob(const Problem &prob, const char *filename)
{
    const auto n_items = prob.conflict_graph.n_vertices();
    std::ofstream file(filename);

    if (!file.is_open())
    {
        LOG_WARN("cannot open file %s", filename);
        return;
    }

    // Problem type
    switch (prob.type)
    {
        case Problem::Type::MinNStacks:
            file << "St";
            break;

        case Problem::Type::MinNBlockingItems:
            file << "BI";
            break;
    }

    // Number of items, capacity, adjacent conflicts, farkas
    file << ' '
         << n_items << ' '
         << prob.capacity << ' '
         << (prob.has_adjacent_conflicts ? 1 : 0) << ' '
         << (prob.farkas ? 1 : 0)
         << std::endl
         << std::endl;

    // Item coefs and offset
    for (std::size_t i = 0; i < n_items; ++i)
    {
        file << prob.coefs[i] << ' ';
    }

    file << prob.offset << std::endl
         << std::endl;

    // Constraints of type same
    file << std::size(prob.same) << std::endl;

    for (auto [i, j] : prob.same)
    {
        file << i + 1 << ' ' << j + 1 << std::endl;
    }

    file << std::endl;

    // Constraints of type diff
    file << std::size(prob.diff) << std::endl;

    for (auto [i, j] : prob.diff)
    {
        file << i + 1 << ' ' << j + 1 << std::endl;
    }

    file << std::endl;

    // Conflict graph
    file << prob.conflict_graph.n_edges() << std::endl;

    for (auto [i, j] : prob.conflict_graph.edges())
    {
        file << i + 1 << ' ' << j + 1 << std::endl;
    }

    file << std::endl;

    // Blocking graph
    file << prob.blocking_graph.n_edges() << std::endl;

    for (auto [i, j] : prob.blocking_graph.edges())
    {
        file << i + 1 << ' ' << j + 1 << std::endl;
    }

    file << std::endl;
}

SCIP_RETCODE Pricer::scip_init(SCIP *scip, [[maybe_unused]] SCIP_PRICER *pricer)
{
    auto &probdata = Probdata::get(scip);
    std::vector<double> coefs;

    // Create an initial instance of the pricing problem
    Problem prob{probdata.type(),
                 probdata.capacity(),
                 probdata.has_adjacent_conflicts(),
                 probdata.conflict_graph(),
                 probdata.conflict_matrix(),
                 probdata.blocking_graph(),
                 probdata.blocking_matrix(),
                 _same,
                 _diff,
                 coefs};

    SCIP_CALL(init(scip, prob));

    return SCIP_OKAY;
}

SCIP_RETCODE Pricer::scip_exit(SCIP *scip, [[maybe_unused]] SCIP_PRICER *pricer)
{
    assert(pricer != nullptr);

    SCIP_CALL(exit(scip));

    LOG_DEBUG("pricer <%s> called %d times, found %d variables",
              SCIPpricerGetName(pricer),
              SCIPpricerGetNCalls(pricer),
              SCIPpricerGetNVarsFound(pricer));

    return SCIP_OKAY;
}

SCIP_RETCODE Pricer::scip_redcost(SCIP *scip,
                                  SCIP_PRICER *pricer,
                                  [[maybe_unused]] SCIP_Real *lowerbound,
                                  [[maybe_unused]] SCIP_Bool *stopearly,
                                  SCIP_RESULT *result)
{
    return pricing(scip, pricer, false, result);
}

SCIP_RETCODE Pricer::scip_farkas(SCIP *scip,
                                 SCIP_PRICER *pricer,
                                 SCIP_RESULT *result)
{
    return pricing(scip, pricer, true, result);
}