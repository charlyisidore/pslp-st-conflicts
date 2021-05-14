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
#include "mip_cplex.hpp"
#include "mip_cplex_3ind.hpp"
#include "mip_cplex_bp.hpp"
#include "mip_cplex_flow.hpp"
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>

static const char USAGE[] =
    R"(CPLEX MIP models for the PSLP.

Usage:
  pslp_mip_cplex -m <model> -i <input> [options]

Options:
  -m --model <file>   Model name (3ind, bp, flow).
  -i --input <file>   Input filename.
  -o --output <file>  Output filename.
  -l --lp <file>      Export the LP to <file>.
  -h --help           Show this screen.
)";

int main(int argc, char **argv)
{
    using pslp::mip::MipCplex;
    using pslp::mip::MipCplex3ind;
    using pslp::mip::MipCplexBp;
    using pslp::mip::MipCplexFlow;

    std::string model;
    std::string input;
    std::string output;
    std::string output_lp;

    // Parse command line options

    {
        int help_flag = 0;

        static struct ::option long_options[] = {
            {"help", no_argument, &help_flag, 1},
            {"model", required_argument, 0, 'm'},
            {"input", required_argument, 0, 'i'},
            {"output", required_argument, 0, 'o'},
            {"lp", required_argument, 0, 'l'},
            {0, 0, 0, 0}};

        while (true)
        {
            int option_index = 0;

            int c = getopt_long(argc,
                                argv,
                                "hm:i:o:l:",
                                long_options,
                                &option_index);

            if (c == -1)
            {
                break;
            }

            switch (c)
            {
                case 0:
                    break;

                case 'h':
                    help_flag = 1;
                    break;

                case 'm':
                    model = optarg;
                    break;

                case 'i':
                    input = optarg;
                    break;

                case 'o':
                    output = optarg;
                    break;

                case 'l':
                    output_lp = optarg;
                    break;

                case '?':
                    // getopt_long already printed an error message
                    std::abort();

                default:
                    std::abort();
            }
        }

        if (help_flag || std::empty(model) || std::empty(input))
        {
            std::cout << USAGE << std::endl;
            return 0;
        }
    }

    std::unique_ptr<MipCplex> mip;

    if (model == "3ind")
    {
        mip = std::make_unique<MipCplex3ind>();
    }
    else if (model == "bp")
    {
        mip = std::make_unique<MipCplexBp>();
    }
    else if (model == "flow")
    {
        mip = std::make_unique<MipCplexFlow>();
    }
    else
    {
        LOG_ERROR("unknown model <%s>", model.c_str());
        return EXIT_FAILURE;
    }

    LOG_INFO("model: %s", model.c_str());

    nlohmann::json in;
    nlohmann::json out = nlohmann::json::object();

    {
        std::ifstream file(input);
        if (!file.is_open())
        {
            LOG_ERROR("cannot open file <%s>", input.c_str());
            return EXIT_FAILURE;
        }
        file >> in;
    }

    auto prob = mip->read(in);

    mip->build(prob);

    if (!std::empty(output_lp))
    {
        mip->write_lp(output_lp.c_str());
    }

    mip->solve();

    MipCplex::Solution sol;

    if (mip->has_solution())
    {
        sol = mip->solution(prob);

        if (!mip->check(prob, sol))
        {
            LOG_ERROR("solution check failed");
        }
    }

    mip->write(prob, sol, out);

    if (std::empty(output))
    {
        std::cout << out << std::endl;
    }
    else
    {
        std::ofstream file(output);
        if (!file.is_open())
        {
            LOG_ERROR("cannot open file <%s>", output.c_str());
            return EXIT_FAILURE;
        }
        file << out;
    }

    return EXIT_SUCCESS;
}
