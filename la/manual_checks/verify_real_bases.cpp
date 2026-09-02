// Manual validation, not part of ctest (needs the downloaded Netlib set):
// solve a handful of real instances with the Phase 1.1 dense simplex,
// reconstruct the ACTUAL optimal basis matrix each one produced, and
// confirm the new Phase 1.2 sparse Markowitz LU factors it and solves it
// correctly — connecting the two phases on real data, not just synthetic
// matrices.

#include <cmath>
#include <cstdio>
#include <vector>

#include "core/sparse.hpp"
#include "io/mps_reader.hpp"
#include "la/basis_factorization.hpp"
#include "la/lu_solve.hpp"
#include "la/markowitz_lu.hpp"
#include "simplex/dense_simplex.hpp"

using namespace inferno;

int main(int argc, char** argv) {
  std::vector<std::string> files;
  for (int i = 1; i < argc; ++i) files.push_back(argv[i]);
  if (files.empty()) {
    files = {"bench/netlib/mps/afiro.mps", "bench/netlib/mps/adlittle.mps",
              "bench/netlib/mps/blend.mps", "bench/netlib/mps/scagr7.mps",
              "bench/netlib/mps/share2b.mps", "bench/netlib/mps/recipe.mps",
              "bench/netlib/mps/israel.mps", "bench/netlib/mps/kb2.mps"};
  }

  int pass = 0, fail = 0;
  for (const auto& path : files) {
    core::LpProblem p;
    try {
      p = io::ReadMps(path);
    } catch (const std::exception& e) {
      printf("%-30s SKIP (read failed: %s)\n", path.c_str(), e.what());
      continue;
    }
    core::Solution sol = simplex::SolveDense(p);
    if (sol.status != core::SolveStatus::kOptimal) {
      printf("%-30s SKIP (dense simplex status=%d, not optimal)\n", path.c_str(), (int)sol.status);
      continue;
    }
    if ((int)sol.basis.size() != p.num_rows) {
      printf("%-30s FAIL (basis size %zu != num_rows %d)\n", path.c_str(), sol.basis.size(), p.num_rows);
      ++fail;
      continue;
    }

    // Reconstruct the m x m basis matrix from the real solve's basis.
    core::CscBuilder builder(p.num_rows, p.num_rows);
    for (int slot = 0; slot < p.num_rows; ++slot) {
      int var = sol.basis[slot];
      if (var < p.num_cols) {
        for (int k = p.a.col_ptr[var]; k < p.a.col_ptr[var + 1]; ++k) {
          builder.AddEntry(slot, p.a.row_idx[k], p.a.values[k]);
        }
      } else {
        int row = var - p.num_cols;
        builder.AddEntry(slot, row, -1.0);  // matches dense_simplex's M = [A | -I]
      }
    }
    core::CscMatrix basis_matrix = std::move(builder).Build();

    la::MarkowitzOptions opts;
    la::LuFactors lu;
    bool ok = la::FactorizeMarkowitz(basis_matrix, opts, core::DefaultTolerances().pivot, lu);
    if (!ok) {
      printf("%-30s FAIL (real optimal basis, %d x %d, did not factor)\n", path.c_str(),
             p.num_rows, p.num_rows);
      ++fail;
      continue;
    }

    std::vector<double> rhs(p.num_rows);
    for (int i = 0; i < p.num_rows; ++i) rhs[i] = 1.0 + 0.01 * i;
    std::vector<std::pair<int, double>> rhs_sparse;
    for (int i = 0; i < p.num_rows; ++i) rhs_sparse.emplace_back(i, rhs[i]);
    std::vector<double> x = la::Ftran(lu, rhs_sparse);

    std::vector<double> bx(p.num_rows, 0.0);
    basis_matrix.MultiplyAdd(x, bx);
    double res = 0.0;
    for (int i = 0; i < p.num_rows; ++i) res = std::max(res, std::abs(bx[i] - rhs[i]));

    bool good = res < 1e-9;
    printf("%-30s %s (m=%d, FTRAN residual=%.3e)\n", path.c_str(), good ? "PASS" : "FAIL",
           p.num_rows, res);
    good ? ++pass : ++fail;
  }

  printf("\n%d/%d real optimal bases factored and solved correctly\n", pass, pass + fail);
  return fail == 0 ? 0 : 1;
}
