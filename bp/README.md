# Branch-and-Price for the Parallel Stack Loading Problem with Conflicts

Implementation of a *Branch-and-Price* approach for the *Parallel Stack Loading Problem* (PSLP) in C++ using [SCIP](https://scipopt.org/).

## Install

Required libraries:

- [SCIP optimization suite](https://www.scipopt.org/).
- [JSON for Modern C++](https://json.nlohmann.me/).
- [CPLEX](https://www.ibm.com/docs/en/icos).

## Usage

```
./pslp_bp [-l <logfile>] [-q] [-s <settings>] [-r <randseed>] [-f <problem>] [-b <batchfile>] [-c "command"]
  -v, --version : print version and build options
  -l <logfile>  : copy output into log file
  -q            : suppress screen messages
  -s <settings> : load parameter settings (.set) file
  -f <problem>  : load and solve problem file
  -o <primref> <dualref> : pass primal and dual objective reference values for validation at the end of the solve
  -b <batchfile>: load and execute dialog command batch file (can be used multiple times)
  -r <randseed> : nonnegative integer to be used as random seed. Has priority over random seed specified through parameter settings (.set) file
  -c "command"  : execute single line of dialog commands (can be used multiple times)
```

See a [list of parameters](PARAMETERS.md) specific to the PSLP solver.

## License

This software is licensed under the GPLv3.
Please see the `LICENSE` file for further information.