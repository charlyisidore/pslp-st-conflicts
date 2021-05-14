# Algorithms for the Parallel Stack Loading Problem with Conflicts

This repository includes:

- A [branch-and-price algorithm](bp) in `bp`.
- [MIP models](mip/cplex) in `mip/cplex`.
- [Standalone variable pricers](pricers) (for benchmarking) in `pricers`.
- [Common files](common) for all projects in `common`.

## Install

```bash
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## License

This software is licensed under the GPLv3.
Please see the `LICENSE` file for further information.