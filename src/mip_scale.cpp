// MILP scale characterisation. Ramps binary count under a controlled time
// limit and reports what actually happened, including where the search
// stops proving optimality. A scale test that hides its unproved runs is
// worthless; the point is to find the limit, not to avoid it.
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "checker/mip_checker.hpp"
#include "core/mip_problem.hpp"
#include "core/sparse.hpp"
#include "mip/branch_and_bound.hpp"

using namespace inferno;

namespace {
struct Rng {
  uint64_t s;
  explicit Rng(unsigned seed) : s(seed * 6364136223846793005ULL + 1442695040888963407ULL) {}
  uint64_t next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return s >> 11; }
  double uniform() { return static_cast<double>(next() % 1000000) / 1000000.0; }
  int below(int n) { return static_cast<int>(next() % static_cast<uint64_t>(n)); }
};

double Now() {
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// Multi-knapsack: n binaries, m capacity rows, each item in a few rows.
// A standard hard-ish MILP family, and one where every instance is
// feasible by construction (all-zero is feasible) so a run measures search
// effort rather than feasibility hunting.
core::MipProblem MakeMultiKnapsack(int n_bin, int m_rows, unsigned seed) {
  Rng rng(seed);
  core::MipProblem p;
  p.lp.name = "MULTIKNAP";
  p.lp.num_cols = n_bin;
  p.lp.col_lo.assign(n_bin, 0.0);
  p.lp.col_hi.assign(n_bin, 1.0);
  p.lp.obj.resize(n_bin);
  p.lp.col_names.assign(n_bin, "x");
  p.is_integer.assign(n_bin, 1);

  std::vector<std::vector<std::pair<int, double>>> cols(n_bin);
  std::vector<double> row_total(m_rows, 0.0);
  for (int j = 0; j < n_bin; ++j) {
    p.lp.obj[j] = -(10.0 + 90.0 * rng.uniform());
    int touches = 2 + rng.below(3);
    for (int k = 0; k < touches; ++k) {
      int r = rng.below(m_rows);
      double w = 1.0 + 20.0 * rng.uniform();
      cols[j].emplace_back(r, w);
      row_total[r] += w;
    }
  }
  for (int i = 0; i < m_rows; ++i) {
    p.lp.row_lo.push_back(-core::kInfinity);
    p.lp.row_hi.push_back(0.4 * row_total[i]);  // binding but satisfiable
    p.lp.row_names.push_back("cap");
  }
  p.lp.num_rows = m_rows;
  core::CscBuilder b(m_rows, n_bin);
  for (int c = 0; c < n_bin; ++c) {
    for (const auto& [r, v] : cols[c]) b.AddEntry(c, r, v);
  }
  p.lp.a = std::move(b).Build();
  return p;
}
}  // namespace

int main(int argc, char** argv) {
  setvbuf(stdout, nullptr, _IONBF, 0);
  double limit = (argc > 1) ? std::stod(argv[1]) : 30.0;

  printf("\nMILP SCALE CHARACTERISATION  (time limit %.0fs per instance)\n", limit);
  printf("Every incumbent is checked by the independent MILP verifier. A run\n");
  printf("that does not prove optimality reports its true remaining gap.\n\n");
  printf("%8s %7s %7s %8s %8s %9s %14s %14s %9s %s\n", "binaries", "rows", "nnz", "intvars",
         "nodes", "time", "objective", "bound", "gap", "status");

  for (int n : {100, 500, 1000, 5000, 10000}) {
    int rows = std::max(10, n / 10);
    auto p = MakeMultiKnapsack(n, rows, 4242u);
    mip::BranchAndBoundOptions o;
    o.time_limit_seconds = limit;
    o.node_limit = 1000000;
    // Bound the work in any single LP so the wall-clock limit is real. A
    // node that hits this is recorded as unresolved, which weakens the
    // reported bound rather than corrupting it.
    o.lp_iteration_limit = 60000;

    double t0 = Now();
    auto s = mip::SolveMip(p, o);
    double dt = Now() - t0;

    const char* status = "-";
    if (s.status == core::SolveStatus::kOptimal) {
      status = s.proved_optimal ? "PROVED" : "incumbent (gap open)";
    } else if (s.status == core::SolveStatus::kInfeasible) {
      status = "infeasible";
    } else {
      status = "root LP exceeded its iteration budget";
    }

    // Print "none" rather than a number when there is no incumbent; a
    // numeric objective and gap there would be actively misleading.
    char objbuf[32], gapbuf[32];
    if (std::isfinite(s.objective_value)) snprintf(objbuf, sizeof objbuf, "%14.6g", s.objective_value);
    else snprintf(objbuf, sizeof objbuf, "%14s", "none");
    if (std::isfinite(s.gap)) snprintf(gapbuf, sizeof gapbuf, "%8.2e", s.gap);
    else snprintf(gapbuf, sizeof gapbuf, "%8s", "n/a");
    char bndbuf[32];
    if (std::isfinite(s.best_bound)) snprintf(bndbuf, sizeof bndbuf, "%14.6g", s.best_bound);
    else snprintf(bndbuf, sizeof bndbuf, "%14s", "none");
    printf("%8d %7d %7d %8d %8d %8.2fs %s %s %s %s\n", n, rows, p.lp.a.nnz(), n,
           s.nodes_explored, dt, objbuf, bndbuf, gapbuf, status);

    if (s.status == core::SolveStatus::kOptimal) {
      auto chk = checker::VerifyMipSolution(p, s);
      printf("%8s   verifier: %s\n", "", chk.message.c_str());
      if (!chk.passed) {
        printf("%8s   *** VERIFICATION FAILED ***\n", "");
        return 1;
      }
    }
  }
  printf("\nRead this honestly: the largest size where optimality is PROVED is the\n");
  printf("real capability. Sizes that return an incumbent with an open gap are\n");
  printf("useful answers but are not proofs, and are labelled as such.\n");
  return 0;
}
