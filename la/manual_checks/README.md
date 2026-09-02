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

**No per-instance timeout.** Unlike `bench/run_netlib.py` (60s subprocess
timeout per instance), this script calls `SolveDense` directly with no
wrapper — on a large pathological instance the dense simplex's own
iteration cap can still mean an astronomically expensive number of O(m^3)
iterations. Running it against the *full* set once left it stuck for over
6 hours on one instance before being killed; running it against a small
hand-picked sample (as the default arg-less invocation does) is fine and
is what's actually been verified. If you want full-set coverage, prefer
timeout-wrapping each call yourself, or wait for Phase 2.1's revised
simplex (`simplex/revised_simplex.*`) to replace the dense one here too.
