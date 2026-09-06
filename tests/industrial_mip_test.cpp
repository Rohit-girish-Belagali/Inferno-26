// Industrial MILP cases. Each is verified independently by
// checker/mip_checker, and each additionally asserts a STRUCTURAL property
// that a wrong answer would violate — an independently verifiable claim
// even where the exact optimum is not hand-computed.
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "checker/mip_checker.hpp"
#include "models/crude_blending.hpp"
#include "models/industrial_mip.hpp"
#include "mip/branch_and_bound.hpp"
#include "simplex/revised_simplex.hpp"

namespace {
int g_fail = 0;
void Expect(bool c, const std::string& m) {
  if (!c) { std::cerr << "FAIL: " << m << "\n"; ++g_fail; }
  else { std::cout << "ok: " << m << "\n"; }
}
using namespace inferno;

// Solves, verifies independently, and prints the search statistics the
// brief asks to record for every MILP.
core::MipSolution RunAndVerify(const std::string& name, const core::MipProblem& p,
                               double time_limit = 20.0) {
  core::Solution relax = simplex::SolveRevised(p.lp);
  mip::BranchAndBoundOptions o;
  o.time_limit_seconds = time_limit;
  auto s = mip::SolveMip(p, o);

  printf("  %-22s rootLP=%-13.6g obj=%-13.6g bound=%-13.6g gap=%.2e nodes=%-6d %s\n",
         name.c_str(),
         relax.status == core::SolveStatus::kOptimal ? relax.objective_value : NAN,
         s.objective_value, s.best_bound, s.gap, s.nodes_explored,
         s.proved_optimal ? "PROVED" : "not proved");

  Expect(s.status == core::SolveStatus::kOptimal, name + ": solved");
  if (s.status == core::SolveStatus::kOptimal) {
    auto chk = checker::VerifyMipSolution(p, s);
    Expect(chk.passed, name + ": independent verifier — " + chk.message);
  }
  return s;
}
}  // namespace

void TestFacility() {
  auto m = models::ExampleFacilityModel();
  auto p = models::BuildFacilityMip(m);
  auto s = RunAndVerify("facility location", p);
  if (s.status != core::SolveStatus::kOptimal) return;
  const int np = 3, nc = 3;
  // Structural: a plant may ship only if it is open, and never above
  // capacity. Recomputed from the raw solution.
  for (int pl = 0; pl < np; ++pl) {
    double shipped = 0.0;
    for (int c = 0; c < nc; ++c) shipped += s.x[np + pl * nc + c];
    bool open = s.x[pl] > 0.5;
    Expect(open || shipped < 1e-6, "facility: closed plants ship nothing");
    Expect(shipped <= m.capacity[pl] + 1e-6, "facility: capacity respected");
  }
  for (int c = 0; c < nc; ++c) {
    double got = 0.0;
    for (int pl = 0; pl < np; ++pl) got += s.x[np + pl * nc + c];
    Expect(got >= m.demand[c] - 1e-6, "facility: demand met");
  }
}

void TestWorkforce() {
  auto m = models::ExampleWorkforceModel();
  auto p = models::BuildWorkforceMip(m);
  auto s = RunAndVerify("workforce scheduling", p);
  if (s.status != core::SolveStatus::kOptimal) return;
  for (size_t sh = 0; sh < m.shift_names.size(); ++sh) {
    double cover = 0.0;
    for (size_t k = 0; k < m.pattern_names.size(); ++k) cover += m.covers[sh][k] * s.x[k];
    Expect(cover >= m.required[sh] - 1e-6, "workforce: " + m.shift_names[sh] + " is staffed");
  }
  for (double v : s.x) {
    Expect(std::abs(v - std::round(v)) < 1e-6, "workforce: whole workers only");
  }
}

void TestLotSizing() {
  auto m = models::ExampleLotSizingModel();
  auto p = models::BuildLotSizingMip(m);
  auto s = RunAndVerify("lot sizing w/ setup", p);
  if (s.status != core::SolveStatus::kOptimal) return;
  const int T = static_cast<int>(m.demand.size());
  for (int t = 0; t < T; ++t) {
    double produced = s.x[t], setup = s.x[T + t];
    Expect(produced < 1e-6 || setup > 0.5,
           "lot sizing: production in t" + std::to_string(t) + " implies its setup is paid");
  }
}

void TestUnitCommitment() {
  auto m = models::ExampleUnitCommitmentModel();
  auto p = models::BuildUnitCommitmentMip(m);
  auto s = RunAndVerify("unit commitment", p, 30.0);
  if (s.status != core::SolveStatus::kOptimal) return;
  const int ng = 3, nt = 4;
  for (int t = 0; t < nt; ++t) {
    double gen = 0.0;
    for (int g = 0; g < ng; ++g) gen += s.x[g * nt + t];
    Expect(std::abs(gen - m.load[t]) < 1e-6, "UC: load met exactly in t" + std::to_string(t));
  }
  // The property the LP relaxation could NOT express: a committed unit
  // must run at or above its minimum stable output, and an uncommitted one
  // must be at zero.
  for (int g = 0; g < ng; ++g) {
    for (int t = 0; t < nt; ++t) {
      double pgen = s.x[g * nt + t], on = s.x[ng * nt + g * nt + t];
      if (on > 0.5) {
        Expect(pgen >= m.pmin[g] - 1e-6, "UC: committed unit respects minimum output");
      } else {
        Expect(pgen < 1e-6, "UC: uncommitted unit generates nothing");
      }
    }
  }
}

void TestRefineryMip() {
  auto m = models::ExampleRefineryMipModel();
  auto p = models::BuildCrudeBlendingMip(m);
  auto s = RunAndVerify("refinery MILP", p, 30.0);
  if (s.status != core::SolveStatus::kOptimal) return;

  const int nc = static_cast<int>(m.blending.crudes.size());
  const int np = static_cast<int>(m.blending.products.size());
  const int base = nc * np;
  // Activation link, and — the part that matters — every quality spec from
  // the LP still holds, recomputed from raw volumes.
  for (int i = 0; i < nc; ++i) {
    double used = 0.0;
    for (int j = 0; j < np; ++j) used += s.x[i * np + j];
    bool active = s.x[base + i] > 0.5;
    Expect(active || used < 1e-6, "refinery MILP: unactivated crude supplies nothing");
  }
  for (int j = 0; j < np; ++j) {
    double vol = 0.0;
    std::vector<double> prop(m.blending.property_names.size(), 0.0);
    for (int i = 0; i < nc; ++i) {
      double v = s.x[i * np + j];
      vol += v;
      for (size_t k = 0; k < prop.size(); ++k) prop[k] += v * m.blending.crudes[i].properties[k];
    }
    Expect(vol >= m.blending.products[j].demand - 1e-6,
           "refinery MILP: " + m.blending.products[j].name + " demand met");
    for (size_t k = 0; k < prop.size(); ++k) {
      if (vol <= 0) continue;
      double avg = prop[k] / vol;
      double hi = m.blending.products[j].max_spec[k];
      double lo = m.blending.products[j].min_spec[k];
      if (std::isfinite(hi)) Expect(avg <= hi + 1e-6, "refinery MILP: max spec preserved");
      if (std::isfinite(lo)) Expect(avg >= lo - 1e-6, "refinery MILP: min spec preserved");
    }
  }
}

int main() {
  printf("\nIndustrial MILP cases — search statistics\n");
  TestFacility();
  TestWorkforce();
  TestLotSizing();
  TestUnitCommitment();
  TestRefineryMip();
  if (g_fail) { std::cerr << g_fail << " test(s) failed\n"; return 1; }
  std::cout << "all industrial MILP tests passed\n";
  return 0;
}
