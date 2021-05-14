# CPLEX MIP models for the PSLP with Conflicts

Implementation of MIP models for the *Parallel Stack Loading Problem* (PSLP) with conflicts in C++ using [CPLEX](https://www.ibm.com/docs/en/icos).

The MIP models are described in:

[Tobias Oelschlägel, & Sigrid Knust (2021). Solution approaches for storage loading problems with stacking constraints. Computers & Operations Research, 127, 105142.](https://doi.org/10.1016/j.cor.2020.105142)

## Install

Required libraries:

- [CPLEX](https://www.ibm.com/docs/en/icos).
- [JSON for Modern C++](https://json.nlohmann.me/).
- [Getopt](https://www.gnu.org/software/libc/manual/html_node/Getopt.html)

## Usage

```
Usage:
  pslp_mip_cplex -m <model> -i <input> [options]

Options:
  -m --model <file>   Model name (3ind, bp, flow).
  -i --input <file>   Input filename.
  -o --output <file>  Output filename.
  -l --lp <file>      Export the LP to <file>.
  -h --help           Show this screen.
```

## License

This software is licensed under the GPLv3.
Please see the `LICENSE` file for further information.