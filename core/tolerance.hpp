#pragma once

namespace inferno::core {

// Every numerical threshold in the solver lives here. Phase 2 standing rule:
// no scattered magic-number tolerances in simplex/presolve/LU code — read
// and write through this struct so a single change can be swept across the
// whole solver and justified in one place.
struct TolerancePolicy {
  double feasibility = 1e-9;      // primal/dual bound violation tolerance
  double optimality = 1e-9;       // reduced-cost / dual feasibility tolerance
  double pivot = 1e-9;            // minimum acceptable pivot magnitude
  double markowitz_threshold = 0.1;  // LU threshold pivoting stability factor
  double growth_refactor = 1e10;  // basis update growth bound before refactorizing
  double checker_residual = 1e-6; // independent solution checker pass/fail bound
  double zero = 1e-12;            // treat magnitudes below this as exact zero
};

inline const TolerancePolicy& DefaultTolerances() {
  static const TolerancePolicy kDefault{};
  return kDefault;
}

}  // namespace inferno::core
