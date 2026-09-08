"""Predefined models, so the application is fully usable with zero API
calls and no internet at all.

These exist for two reasons. The obvious one is the demo: a venue's wifi
is not something to bet a presentation on, and every one of these solves
against the real Inferno engine with no language model anywhere in the
path. The less obvious one is that they are a fixed reference: when an
AI-generated model behaves oddly, these are the known-good formulations to
compare it against.

Each is written in exactly the schema the validator accepts, so they take
the identical path an AI model takes -- validate, compile, solve, verify.
Nothing here gets a shortcut the AI does not get.
"""


def _var(name, kind, lo, hi, desc):
    return {"name": name, "type": kind, "lower_bound": lo, "upper_bound": hi,
            "description": desc}


def _term(v, c):
    return {"variable": v, "coefficient": c}


def _con(name, terms, lo, hi, desc):
    return {"name": name, "terms": terms, "lower_bound": lo, "upper_bound": hi,
            "description": desc}


def production_planning():
    """Level 1 LP. Two products competing for three shared machines."""
    return {
        "name": "production_planning",
        "objective_sense": "maximize",
        "variables": [
            _var("chairs", "continuous", 0, None, "Chairs produced per week"),
            _var("tables", "continuous", 0, None, "Tables produced per week"),
        ],
        "objective": [_term("chairs", 45.0), _term("tables", 80.0)],
        "constraints": [
            _con("cutting", [_term("chairs", 2.0), _term("tables", 4.0)], None, 480.0,
                 "Cutting shop hours available per week"),
            _con("assembly", [_term("chairs", 3.0), _term("tables", 5.0)], None, 600.0,
                 "Assembly hours available per week"),
            _con("finishing", [_term("chairs", 1.0), _term("tables", 3.0)], None, 300.0,
                 "Finishing hours available per week"),
        ],
    }


def transportation():
    """Level 1 LP. Classical Hitchcock transportation: two plants, three
    warehouses, minimize freight."""
    costs = {("p1", "w1"): 4.0, ("p1", "w2"): 6.0, ("p1", "w3"): 9.0,
             ("p2", "w1"): 5.0, ("p2", "w2"): 3.0, ("p2", "w3"): 8.0}
    supply = {"p1": 120.0, "p2": 140.0}
    demand = {"w1": 100.0, "w2": 90.0, "w3": 60.0}
    variables, objective = [], []
    for (p, w), c in costs.items():
        name = f"ship_{p}_{w}"
        variables.append(_var(name, "continuous", 0, None, f"Units shipped {p} to {w}"))
        objective.append(_term(name, c))
    constraints = []
    for p, cap in supply.items():
        constraints.append(_con(
            f"supply_{p}", [_term(f"ship_{p}_{w}", 1.0) for w in demand], None, cap,
            f"Plant {p} cannot ship more than it makes"))
    for w, d in demand.items():
        constraints.append(_con(
            f"demand_{w}", [_term(f"ship_{p}_{w}", 1.0) for p in supply], d, None,
            f"Warehouse {w} demand must be met"))
    return {"name": "transportation", "objective_sense": "minimize",
            "variables": variables, "objective": objective, "constraints": constraints}


def diet():
    """Level 1 LP. Stigler-style diet: cheapest basket meeting nutrition
    floors without exceeding sodium."""
    foods = [
        ("oats", 0.28, 389, 16.9, 6.0),
        ("milk", 0.55, 61, 3.2, 44.0),
        ("eggs", 1.10, 155, 12.6, 124.0),
        ("beans", 0.45, 127, 8.7, 2.0),
        ("spinach", 0.62, 23, 2.9, 79.0),
    ]
    variables = [_var(n, "continuous", 0, 6.0, f"Servings of {n} per day")
                 for n, _, _, _, _ in foods]
    objective = [_term(n, c) for n, c, _, _, _ in foods]
    constraints = [
        _con("energy", [_term(n, kcal) for n, _, kcal, _, _ in foods], 2000.0, None,
             "Daily energy floor in kcal"),
        _con("protein", [_term(n, pro) for n, _, _, pro, _ in foods], 60.0, None,
             "Daily protein floor in grams"),
        _con("sodium", [_term(n, na) for n, _, _, _, na in foods], None, 900.0,
             "Daily sodium ceiling in milligrams"),
    ]
    return {"name": "diet", "objective_sense": "minimize", "variables": variables,
            "objective": objective, "constraints": constraints}


def plant_selection():
    """Level 2 MILP. Fixed opening costs make this genuinely discrete: the
    linking rows are what force a plant to be opened before it can produce."""
    plants = [("A", 50000.0, 500.0, 20.0), ("B", 70000.0, 800.0, 15.0),
              ("C", 40000.0, 400.0, 26.0)]
    demand = 1000.0
    variables, objective, constraints = [], [], []
    for name, fixed, cap, unit in plants:
        variables.append(_var(f"open_{name}", "binary", 0, 1, f"Open plant {name}"))
        variables.append(_var(f"prod_{name}", "continuous", 0, cap,
                              f"Units produced at plant {name}"))
        objective.append(_term(f"open_{name}", fixed))
        objective.append(_term(f"prod_{name}", unit))
        constraints.append(_con(
            f"link_{name}", [_term(f"prod_{name}", 1.0), _term(f"open_{name}", -cap)],
            None, 0.0, f"Plant {name} can only produce if it is opened"))
    constraints.insert(0, _con(
        "demand", [_term(f"prod_{n}", 1.0) for n, _, _, _ in plants], demand, None,
        "Total production must meet demand"))
    return {"name": "plant_selection", "objective_sense": "minimize",
            "variables": variables, "objective": objective, "constraints": constraints}


def workforce_scheduling():
    """Level 2 MILP. Cyclic shift cover: each shift covers two consecutive
    periods, so the rows wrap around the day."""
    periods = ["00-06", "06-12", "12-18", "18-24"]
    need = [6, 18, 22, 12]
    cost = [1400.0, 1000.0, 1050.0, 1300.0]
    variables = [_var(f"shift_{i}", "integer", 0, 40, f"Workers starting at {periods[i]}")
                 for i in range(4)]
    objective = [_term(f"shift_{i}", cost[i]) for i in range(4)]
    constraints = []
    for i in range(4):
        prev = (i - 1) % 4
        constraints.append(_con(
            f"cover_{i}", [_term(f"shift_{i}", 1.0), _term(f"shift_{prev}", 1.0)],
            float(need[i]), None, f"Coverage required during {periods[i]}"))
    return {"name": "workforce_scheduling", "objective_sense": "minimize",
            "variables": variables, "objective": objective, "constraints": constraints}


def production_with_setup():
    """Level 3 MILP. Setup costs plus a shared capacity — the weak-relaxation
    shape that makes branch-and-bound actually work for its bound."""
    items = [("P1", 1200.0, 38.0, 260.0), ("P2", 1800.0, 25.0, 300.0),
             ("P3", 900.0, 44.0, 180.0)]
    variables, objective, constraints = [], [], []
    for name, setup, margin, cap in items:
        variables.append(_var(f"run_{name}", "binary", 0, 1, f"Set up a run of {name}"))
        variables.append(_var(f"qty_{name}", "continuous", 0, cap, f"Units of {name} made"))
        objective.append(_term(f"run_{name}", -setup))
        objective.append(_term(f"qty_{name}", margin))
        constraints.append(_con(
            f"setup_{name}", [_term(f"qty_{name}", 1.0), _term(f"run_{name}", -cap)],
            None, 0.0, f"{name} can only be made after its setup is paid"))
    constraints.append(_con(
        "machine_hours",
        [_term("qty_P1", 1.4), _term("qty_P2", 1.0), _term("qty_P3", 1.9)],
        None, 640.0, "Shared machine hours available this month"))
    return {"name": "production_with_setup", "objective_sense": "maximize",
            "variables": variables, "objective": objective, "constraints": constraints}


def refinery_blending():
    """PS 26119's own problem: blend crudes into a product slate under
    sulfur and yield specifications. Volume-weighted quality is cleared to
    linear form, which is what keeps it an LP-shaped row set; the binaries
    decide which crudes are actually brought online."""
    crudes = [
        # name, cost/bbl, sulfur %, naphtha yield, distillate yield, availability
        ("arab_light", 5400.0, 1.80, 0.24, 0.36, 40000.0),
        ("arab_heavy", 4900.0, 2.90, 0.18, 0.31, 35000.0),
        ("bombay_high", 6100.0, 0.15, 0.31, 0.39, 22000.0),
        ("murban", 5800.0, 0.78, 0.28, 0.38, 28000.0),
    ]
    throughput = 90000.0
    variables, objective, constraints = [], [], []
    for name, cost, _s, _n, _d, avail in crudes:
        variables.append(_var(f"use_{name}", "binary", 0, 1, f"Bring {name} online"))
        variables.append(_var(f"vol_{name}", "continuous", 0, avail,
                              f"Barrels per day of {name}"))
        objective.append(_term(f"vol_{name}", cost))
        constraints.append(_con(
            f"online_{name}", [_term(f"vol_{name}", 1.0), _term(f"use_{name}", -avail)],
            None, 0.0, f"{name} can only be run if its train is brought online"))
        constraints.append(_con(
            f"minrun_{name}", [_term(f"vol_{name}", 1.0), _term(f"use_{name}", -8000.0)],
            0.0, None, f"An online {name} train must run at least its minimum rate"))
    constraints.insert(0, _con(
        "throughput", [_term(f"vol_{n}", 1.0) for n, _, _, _, _, _ in crudes],
        throughput, throughput, "Crude unit throughput, barrels per day"))
    # Sulfur cap, cleared to linear form: sum(s_i v_i) <= spec * sum(v_i).
    # With throughput fixed, the right-hand side is a constant.
    constraints.append(_con(
        "sulfur_spec", [_term(f"vol_{n}", s) for n, _, s, _, _, _ in crudes],
        None, 1.55 * throughput, "Blended sulfur must stay within the 1.55% specification"))
    constraints.append(_con(
        "naphtha_floor", [_term(f"vol_{n}", y) for n, _, _, y, _, _ in crudes],
        0.25 * throughput, None, "Naphtha yield must cover downstream reformer feed"))
    constraints.append(_con(
        "distillate_floor", [_term(f"vol_{n}", y) for n, _, _, _, y, _ in crudes],
        0.355 * throughput, None, "Distillate yield must cover the diesel commitment"))
    return {"name": "refinery_blending", "objective_sense": "minimize",
            "variables": variables, "objective": objective, "constraints": constraints}


def power_dispatch():
    """Level 3 MILP. Unit commitment over four periods with ramp coupling —
    the discrete on/off decision is what separates this from economic
    dispatch, which is the continuous sub-problem."""
    units = [("coal", 30000.0, 2.6, 120.0, 420.0), ("gas", 12000.0, 4.9, 0.0, 300.0),
             ("hydro", 4000.0, 1.1, 0.0, 160.0)]
    load = [420.0, 610.0, 730.0, 540.0]
    variables, objective, constraints = [], [], []
    for name, startup, rate, minout, maxout in units:
        variables.append(_var(f"on_{name}", "binary", 0, 1, f"Commit {name} for the day"))
        objective.append(_term(f"on_{name}", startup))
        for t in range(4):
            variables.append(_var(f"gen_{name}_{t}", "continuous", 0, maxout,
                                  f"MW from {name} in period {t}"))
            objective.append(_term(f"gen_{name}_{t}", rate * 100.0))
            constraints.append(_con(
                f"cap_{name}_{t}",
                [_term(f"gen_{name}_{t}", 1.0), _term(f"on_{name}", -maxout)],
                None, 0.0, f"{name} output in period {t} requires commitment"))
            if minout > 0:
                constraints.append(_con(
                    f"min_{name}_{t}",
                    [_term(f"gen_{name}_{t}", 1.0), _term(f"on_{name}", -minout)],
                    0.0, None, f"{name} has a technical minimum when committed"))
    for t in range(4):
        constraints.append(_con(
            f"load_{t}", [_term(f"gen_{n}_{t}", 1.0) for n, _, _, _, _ in units],
            load[t], None, f"Demand of {load[t]:g} MW in period {t} must be served"))
    for t in range(1, 4):
        constraints.append(_con(
            f"ramp_coal_{t}",
            [_term(f"gen_coal_{t}", 1.0), _term(f"gen_coal_{t-1}", -1.0)],
            -180.0, 180.0, f"Coal ramp limit between periods {t-1} and {t}"))
    return {"name": "power_dispatch", "objective_sense": "minimize",
            "variables": variables, "objective": objective, "constraints": constraints}


CATALOG = [
    {"id": "production_planning", "title": "Production planning",
     "level": "Level 1 · LP", "prompt": "Maximize weekly profit from chairs and tables "
     "given cutting, assembly and finishing hours.", "build": production_planning},
    {"id": "transportation", "title": "Transportation",
     "level": "Level 1 · LP", "prompt": "Ship from two plants to three warehouses at "
     "minimum freight cost.", "build": transportation},
    {"id": "diet", "title": "Diet",
     "level": "Level 1 · LP", "prompt": "Find the cheapest daily food basket that meets "
     "energy and protein floors without exceeding sodium.", "build": diet},
    {"id": "plant_selection", "title": "Plant selection",
     "level": "Level 2 · MILP", "prompt": "Choose which plants to open when each has a "
     "fixed opening cost and a capacity.", "build": plant_selection},
    {"id": "workforce_scheduling", "title": "Workforce scheduling",
     "level": "Level 2 · MILP", "prompt": "Cover a 24-hour roster with overlapping "
     "shifts at minimum wage cost.", "build": workforce_scheduling},
    {"id": "refinery_blending", "title": "Refinery crude blending",
     "level": "Level 2 · MILP · PS 26119", "prompt": "Blend crudes to meet throughput, "
     "sulfur and yield specifications at minimum cost.", "build": refinery_blending},
    {"id": "production_with_setup", "title": "Production with setup costs",
     "level": "Level 3 · MILP", "prompt": "Maximize margin when each product run "
     "carries a fixed setup cost and machine hours are shared.",
     "build": production_with_setup},
    {"id": "power_dispatch", "title": "Power unit commitment",
     "level": "Level 3 · MILP", "prompt": "Commit generating units and dispatch them "
     "across four periods under ramp limits.", "build": power_dispatch},
]


def catalog():
    return [{k: v for k, v in item.items() if k != "build"} for item in CATALOG]


def build(example_id):
    for item in CATALOG:
        if item["id"] == example_id:
            return item["build"]()
    return None
