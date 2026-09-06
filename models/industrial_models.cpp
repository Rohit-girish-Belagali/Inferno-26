#include "models/industrial_models.hpp"

#include <cmath>
#include <cstdint>
#include <vector>

#include "core/sparse.hpp"

namespace inferno::models {

using core::kInfinity;

namespace {
// Deterministic generator, so a stress instance is identical run to run
// and a regression is reproducible rather than a story about one bad
// afternoon. Plain 64-bit LCG; the statistical quality does not matter
// here, only repeatability.
struct Rng {
  uint64_t s;
  explicit Rng(unsigned seed) : s(seed * 6364136223846793005ULL + 1442695040888963407ULL) {}
  uint64_t next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return s >> 11; }
  double uniform() { return static_cast<double>(next() % 1000000) / 1000000.0; }
  int below(int n) { return static_cast<int>(next() % static_cast<uint64_t>(n)); }
};

core::LpProblem Finish(core::LpProblem lp,
                        std::vector<std::vector<std::pair<int, double>>>& cols, int rows) {
  lp.num_rows = rows;
  core::CscBuilder b(rows, lp.num_cols);
  for (int c = 0; c < lp.num_cols; ++c) {
    for (const auto& [r, v] : cols[c]) {
      if (v != 0.0) b.AddEntry(c, r, v);
    }
  }
  lp.a = std::move(b).Build();
  return lp;
}
}  // namespace

core::LpProblem BuildProductionLp(const ProductionModel& m) {
  const int np = static_cast<int>(m.product_names.size());
  const int nr = static_cast<int>(m.resource_names.size());
  core::LpProblem lp;
  lp.name = "PRODUCTION";
  lp.num_cols = np;
  lp.col_lo.assign(np, 0.0);
  lp.col_hi.assign(np, kInfinity);
  lp.obj.resize(np);
  lp.col_names = m.product_names;
  for (int j = 0; j < np; ++j) lp.obj[j] = -m.profit[j];  // minimise negated profit

  std::vector<std::vector<std::pair<int, double>>> cols(np);
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < np; ++j) cols[j].emplace_back(i, m.usage[i][j]);
    lp.row_lo.push_back(-kInfinity);
    lp.row_hi.push_back(m.available[i]);
    lp.row_names.push_back(m.resource_names[i]);
  }
  return Finish(std::move(lp), cols, nr);
}

ProductionModel ExampleProductionModel() {
  // The classic two-product mix: A earns 40 and needs 2 labour + 1
  // material; B earns 30 and needs 1 labour + 2 material; 100 labour and
  // 80 material available. Optimum is A=40, B=20 for a profit of 2200,
  // with both resources exactly exhausted — verifiable by hand, which is
  // the point of keeping it in the battery.
  ProductionModel m;
  m.product_names = {"ProductA", "ProductB"};
  m.profit = {40.0, 30.0};
  m.resource_names = {"Labour_hours", "Material_kg"};
  m.available = {100.0, 80.0};
  m.usage = {{2.0, 1.0}, {1.0, 2.0}};
  return m;
}

core::LpProblem BuildDietLp(const DietModel& m) {
  const int nf = static_cast<int>(m.food_names.size());
  const int nn = static_cast<int>(m.nutrient_names.size());
  core::LpProblem lp;
  lp.name = "DIET";
  lp.num_cols = nf;
  lp.col_lo.assign(nf, 0.0);
  lp.col_hi.assign(nf, kInfinity);
  lp.obj = m.cost;
  lp.col_names = m.food_names;

  std::vector<std::vector<std::pair<int, double>>> cols(nf);
  for (int i = 0; i < nn; ++i) {
    for (int j = 0; j < nf; ++j) cols[j].emplace_back(i, m.content[i][j]);
    lp.row_lo.push_back(m.min_required[i]);
    lp.row_hi.push_back(m.max_allowed[i]);
    lp.row_names.push_back(m.nutrient_names[i]);
  }
  return Finish(std::move(lp), cols, nn);
}

DietModel ExampleDietModel() {
  DietModel m;
  m.food_names = {"Rice", "Dal", "Milk", "Egg", "Spinach", "Oil"};
  m.cost = {40.0, 110.0, 55.0, 7.0, 30.0, 140.0};  // rupees per unit
  m.nutrient_names = {"Calories", "Protein_g", "Carbs_g", "VitaminA_mcg"};
  m.min_required = {2200.0, 55.0, 250.0, 700.0};
  m.max_allowed = {2900.0, kInfinity, kInfinity, kInfinity};
  //                Rice   Dal    Milk   Egg    Spinach Oil
  m.content = {
      {1300.0, 1200.0, 600.0,  70.0,   25.0,  880.0},   // calories
      {  24.0,   72.0,  32.0,   6.2,    3.0,    0.0},   // protein
      { 280.0,  200.0,  48.0,   0.6,    4.0,    0.0},   // carbohydrate
      {   0.0,   10.0, 320.0, 160.0, 9400.0,    0.0},   // vitamin A
  };
  return m;
}

core::LpProblem BuildDegenerateLp(int blocks) {
  // Transportation-shaped: `blocks` sources each supplying `blocks` sinks,
  // with every cost identical. Identical costs mean an enormous number of
  // distinct bases achieve the same objective, which is precisely the
  // degeneracy the PS asks to be stressed.
  const int n = blocks;
  const int ncols = n * n;
  core::LpProblem lp;
  lp.name = "DEGENERATE";
  lp.num_cols = ncols;
  lp.col_lo.assign(ncols, 0.0);
  lp.col_hi.assign(ncols, kInfinity);
  lp.obj.assign(ncols, 1.0);  // every route costs the same
  lp.col_names.assign(ncols, std::string());

  std::vector<std::vector<std::pair<int, double>>> cols(ncols);
  int row = 0;
  for (int i = 0; i < n; ++i) {  // supply, all equal
    for (int j = 0; j < n; ++j) cols[i * n + j].emplace_back(row, 1.0);
    lp.row_lo.push_back(-kInfinity);
    lp.row_hi.push_back(1.0);
    lp.row_names.push_back("S" + std::to_string(i));
    ++row;
  }
  for (int j = 0; j < n; ++j) {  // demand, all equal
    for (int i = 0; i < n; ++i) cols[i * n + j].emplace_back(row, 1.0);
    lp.row_lo.push_back(1.0);
    lp.row_hi.push_back(kInfinity);
    lp.row_names.push_back("D" + std::to_string(j));
    ++row;
  }
  return Finish(std::move(lp), cols, row);
}

core::LpProblem BuildIllConditionedLp(int n, int decades) {
  // Column j is scaled by 10^(exponent), sweeping the full range from
  // -decades to +decades. The optimum is unaffected by that scaling in
  // exact arithmetic, so any error in the reported objective is purely
  // numerical -- which makes this a clean instrument rather than just a
  // hard instance.
  core::LpProblem lp;
  lp.name = "ILLCOND";
  lp.num_cols = n;
  lp.col_lo.assign(n, 0.0);
  lp.col_hi.assign(n, kInfinity);
  lp.obj.resize(n);
  lp.col_names.assign(n, std::string());

  std::vector<double> scale(n);
  for (int j = 0; j < n; ++j) {
    double t = (n == 1) ? 0.0 : (2.0 * j / (n - 1) - 1.0);  // -1 .. +1
    scale[j] = std::pow(10.0, t * decades);
    lp.obj[j] = -scale[j];  // maximise sum of scaled variables
  }

  std::vector<std::vector<std::pair<int, double>>> cols(n);
  int row = 0;
  for (int i = 0; i < n; ++i) {
    cols[i].emplace_back(row, scale[i]);
    if (i + 1 < n) cols[i + 1].emplace_back(row, scale[i + 1]);
    lp.row_lo.push_back(-kInfinity);
    lp.row_hi.push_back(scale[i]);
    lp.row_names.push_back("C" + std::to_string(i));
    ++row;
  }
  return Finish(std::move(lp), cols, row);
}

core::LpProblem BuildLargeSparseLp(int rows, int cols_n, int nnz_per_col, unsigned seed) {
  Rng rng(seed);
  core::LpProblem lp;
  lp.name = "LARGE_SPARSE";
  lp.num_cols = cols_n;
  lp.col_lo.assign(cols_n, 0.0);
  lp.col_hi.assign(cols_n, 10.0);
  lp.obj.resize(cols_n);
  lp.col_names.assign(cols_n, std::string());

  std::vector<std::vector<std::pair<int, double>>> cols(cols_n);
  std::vector<double> row_load(rows, 0.0);
  for (int j = 0; j < cols_n; ++j) {
    lp.obj[j] = -(1.0 + rng.uniform());
    for (int k = 0; k < nnz_per_col; ++k) {
      int r = rng.below(rows);
      double v = 0.5 + rng.uniform();
      cols[j].emplace_back(r, v);
      row_load[r] += v;
    }
  }
  for (int i = 0; i < rows; ++i) {
    // Right-hand side set from the row's own coefficient mass, so the
    // instance is always feasible and always has a bounded optimum --
    // otherwise a scale test measures how fast the solver detects
    // unboundedness, which is not the thing being tested.
    lp.row_lo.push_back(-kInfinity);
    lp.row_hi.push_back(0.25 * row_load[i] + 1.0);
    lp.row_names.push_back(std::string());
  }
  return Finish(std::move(lp), cols, rows);
}

}  // namespace inferno::models
