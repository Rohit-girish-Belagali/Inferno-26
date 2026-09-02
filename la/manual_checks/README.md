# Manual checks

Not part of `ctest` — these need the downloaded Netlib set
(`bash bench/download_netlib.sh` first) and are for validating `la/`
against real data by hand, not for CI.

## verify_real_bases.cpp

Solves each given Netlib instance with the Phase 1.1 dense simplex,
reconstructs the *actual* optimal basis matrix that solve landed on, and
checks the Phase 1.2 sparse Markowitz LU factors and solves it correctly
too (residual against a fixed test RHS, compared to the basis's real
columns — same "recompute from raw data" spirit as `checker/`).

```bash
cmake --build build -j
clang++ -std=c++20 -O2 -I. la/manual_checks/verify_real_bases.cpp \
  build/CMakeFiles/inferno_core.dir/io/mps_reader.cpp.o \
  build/CMakeFiles/inferno_core.dir/simplex/dense_simplex.cpp.o \
  build/CMakeFiles/inferno_core.dir/la/markowitz_lu.cpp.o \
  build/CMakeFiles/inferno_core.dir/la/lu_factors.cpp.o \
  build/CMakeFiles/inferno_core.dir/la/lu_solve.cpp.o \
  build/CMakeFiles/inferno_core.dir/la/basis_factorization.cpp.o \
  -o /tmp/verify_real_bases

/tmp/verify_real_bases                        # a small hand-picked sample
/tmp/verify_real_bases bench/netlib/mps/*.mps  # every decoded instance
```

Last full-set result: see the commit that introduced this file's history,
or just rerun it — it takes as long as the dense simplex takes to solve
whatever's passed in (which for the full set is the same ~20 minutes as
`bench/run_netlib.py`, since it's driven by the same throwaway solver).
