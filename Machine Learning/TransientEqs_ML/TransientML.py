import os
import re
import numpy as np
import pandas as pd
from pysr import PySRRegressor

print("=" * 40)
print(" START TRANSIENT SYMBOLIC REGRESSION ")
print("=" * 40)

# ============================================================================
# 1. LOAD DATA
# ============================================================================
filename = "IGA_Transient_Database.csv"

if not os.path.exists(filename):
    raise FileNotFoundError(f"CSV file '{filename}' was not found.")

df = pd.read_csv(filename)

required_cols = [
    "L", "R", "kx", "ky", "rho", "c",
    "T_initial", "T_hole", "T_outer",
    "x", "y", "t", "T"
]

missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns in CSV: {missing}")

df = df.replace([np.inf, -np.inf], np.nan).dropna().copy()
print(f"Loaded rows: {len(df)}")

# ============================================================================
# 2. FEATURE ENGINEERING
# ============================================================================
L = df["L"].astype(float)
R = df["R"].astype(float)

cx = L / 2.0
cy = L / 2.0

dx = df["x"].astype(float) - cx
dy = df["y"].astype(float) - cy
r_dist = np.sqrt(dx**2 + dy**2)

eps = 1e-12
r_safe = np.maximum(r_dist, eps)

# Dimensionless geometric features
df["rhat"] = r_dist / R
df["cost"] = dx / r_safe
df["sint"] = dy / r_safe
df["RL"]   = R / L

# Material / transient features
df["anis"] = np.sqrt(df["kx"] / df["ky"])
df["alpha_eff"] = np.sqrt(df["kx"] * df["ky"]) / (df["rho"] * df["c"])
df["Fo"] = df["alpha_eff"] * df["t"] / (df["L"] ** 2)

# Thermal feature
den = df["T_initial"] - df["T_outer"]
df = df.loc[np.abs(den) > 1e-12].copy()
df["beta_bc"] = (df["T_hole"] - df["T_outer"]) / (df["T_initial"] - df["T_outer"])

# Dimensionless target
df["u"] = (df["T"] - df["T_outer"]) / (df["T_initial"] - df["T_outer"])

feature_cols = ["rhat", "cost", "sint", "anis", "RL", "Fo", "beta_bc"]

mask_finite = np.isfinite(df[feature_cols + ["u"]]).all(axis=1)
df = df.loc[mask_finite].copy()

print(f"Usable rows after cleaning: {len(df)}")

if len(df) == 0:
    raise ValueError("No valid rows remain after preprocessing.")

# ============================================================================
# 3. BALANCE SIMULATIONS WITH WEIGHTS
#    Each simulation has many spatial samples; this makes each simulation
#    contribute approximately equally to the loss.
# ============================================================================
sim_cols = ["L", "R", "kx", "ky", "rho", "c", "T_initial", "T_hole", "T_outer"]

# Robust grouping for float columns
sim_key = (
    df[sim_cols]
    .round(10)
    .astype(str)
    .agg("|".join, axis=1)
)

sim_counts = sim_key.map(sim_key.value_counts())
weights = 1.0 / sim_counts.to_numpy(dtype=float)
weights = weights / np.mean(weights)

# ============================================================================
# 4. DATA FOR PySR
# ============================================================================
X = df[feature_cols].copy()
y = df["u"].to_numpy(dtype=np.float32)

# Variable-specific complexities:
# Penalize rhat more, so the search is less likely to collapse to rhat alone.
# Lower complexity means "cheaper" to use.
complexity_list = [4.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0]

# These are the feature groups we want in the final equation.
# Requiring all 7 can be too strict; this set ensures the final formula
# depends on space, geometry, material, transient response, and thermal BC level.
must_have = ["rhat", "anis", "RL", "Fo", "beta_bc"]

# ============================================================================
# 5. CONFIGURE PySR
# ============================================================================
model = PySRRegressor(
    niterations=1200,
    model_selection="accuracy",

    # Keep operators physically interpretable.
    binary_operators=["+", "-", "*", "/"],
    unary_operators=["square", "sqrt"],

    # Constrain operator usage to discourage pathological expressions.
    constraints={
        "/": (-1, 1),   # denominator must stay simple
        "square": 5,
        "sqrt": 5,
    },

    nested_constraints={
        "square": {"square": 1, "sqrt": 1},
        "sqrt": {"sqrt": 1, "square": 1},
    },

    complexity_of_constants=5.0,
    parsimony=0.0005,
    maxsize=35,
    warmup_maxsize_by=0.4,

    batching=True,
    batch_size=512,
    precision=32,
    turbo=True,

    verbosity=1,
)

# ============================================================================
# 6. FIT
# ============================================================================
print("\nStarting symbolic discovery...")
model.fit(
    X,
    y,
    weights=weights,
    complexity_of_variables=complexity_list,
)

# Save all equations
eq_df = model.equations_.copy()
eq_df.to_csv("Transient_PySR_Equations.csv", index=False)
print("\nSaved PySR equation table to: Transient_PySR_Equations.csv")

# ============================================================================
# 7. PICK AN EQUATION THAT USES THE REQUIRED PHYSICAL FEATURES
# ============================================================================
def count_used_features(eq_string, features):
    used = []
    for feat in features:
        if re.search(rf"\b{re.escape(feat)}\b", str(eq_string)):
            used.append(feat)
    return used

eq_df["used_features"] = eq_df["equation"].apply(lambda s: count_used_features(s, feature_cols))
eq_df["n_used"] = eq_df["used_features"].apply(len)

eq_df["uses_required"] = eq_df["equation"].apply(
    lambda s: all(re.search(rf"\b{re.escape(f)}\b", str(s)) for f in must_have)
)

full_candidates = eq_df[eq_df["uses_required"]].copy()

if len(full_candidates) > 0:
    chosen_idx = full_candidates.sort_values(["loss", "complexity"], ascending=[True, True]).index[0]
    selection_reason = "best loss among equations using all required physical features"
else:
    # fallback: maximize number of required features used, then minimize loss
    eq_df["n_required_used"] = eq_df["equation"].apply(
        lambda s: sum(bool(re.search(rf"\b{re.escape(f)}\b", str(s))) for f in must_have)
    )
    chosen_idx = eq_df.sort_values(
        ["n_required_used", "loss", "complexity"],
        ascending=[False, True, True]
    ).index[0]
    selection_reason = "no equation used all required features; chose the one using the most required features"

best_row = eq_df.loc[chosen_idx]
best_eq = model.sympy(index=int(chosen_idx))

print("\nTop discovered equations:")
try:
    print(eq_df[["complexity", "loss", "equation"]].tail(10))
except Exception:
    print(eq_df)

print("\n" + "=" * 80)
print(" SELECTED TRANSIENT ANALYTICAL EQUATION ")
print("=" * 80)
print(f"Selection reason: {selection_reason}")
print(f"Used features: {best_row['used_features']}")
print(f"u = {best_eq}")

# ============================================================================
# 8. SAFE MATLAB EXPORT
# ============================================================================
def sympy_to_matlab_elementwise(expr):
    expr_str = str(expr)

    protected = {
        "rhat": "__RHAT__",
        "cost": "__COST__",
        "sint": "__SINT__",
        "anis": "__ANIS__",
        "RL": "__RL__",
        "Fo": "__FO__",
        "beta_bc": "__BETA_BC__",
    }

    for var, token in protected.items():
        expr_str = re.sub(rf"\b{var}\b", token, expr_str)

    expr_str = expr_str.replace("**", "@@POW@@")
    expr_str = expr_str.replace("*", ".*")
    expr_str = expr_str.replace("/", "./")
    expr_str = expr_str.replace("@@POW@@", ".^")

    matlab_map = {
        "__RHAT__": "(sqrt((x - L./2).^2 + (y - L./2).^2) ./ R)",
        "__COST__": "((x - L./2) ./ max(sqrt((x - L./2).^2 + (y - L./2).^2), 1e-12))",
        "__SINT__": "((y - L./2) ./ max(sqrt((x - L./2).^2 + (y - L./2).^2), 1e-12))",
        "__ANIS__": "sqrt(kx ./ ky)",
        "__RL__": "(R ./ L)",
        "__FO__": "((sqrt(kx .* ky) ./ (rho_m .* c)) .* t ./ (L.^2))",
        "__BETA_BC__": "((T_hole - T_outer) ./ (T_initial - T_outer))",
    }

    for token, repl in matlab_map.items():
        expr_str = expr_str.replace(token, repl)

    return expr_str

eq_mat = sympy_to_matlab_elementwise(best_eq)

matlab_function = f"""
function T = Predict_Temperature_Transient_AI(x, y, t, L, R, kx, ky, rho_m, c, T_initial, T_hole, T_outer)
    % AI-discovered transient temperature predictor
    %
    % Inputs:
    %   x, y       : coordinate arrays
    %   t          : time where to analyze (e.g. 5 sec)
    %   L, R       : plate size and hole radius
    %   kx, ky     : conductivities
    %   rho_m, c   : density and specific heat
    %   T_initial  : initial temperature
    %   T_hole     : hole boundary temperature
    %   T_outer    : outer boundary temperature
    %
    % Output:
    %   T          : predicted temperature field

    % Selected equation:
    % u = {str(best_eq)}
    u_shape = {eq_mat};

    % Recover dimensional temperature
    T = T_outer + (T_initial - T_outer) .* u_shape;

    % Physical clamping
    Tmin = min([T_initial, T_hole, T_outer]);
    Tmax = max([T_initial, T_hole, T_outer]);
    T = max(Tmin, min(Tmax, T));

    % Geometry mask
    cx = L ./ 2;
    cy = L ./ 2;
    dist = sqrt((x - cx).^2 + (y - cy).^2);

    is_valid = (x >= 0) & (x <= L) & ...
               (y >= 0) & (y <= L) & ...
               (dist >= R) & ...
               (t >= 0);

    T(~is_valid) = NaN;
end
"""

print("\n" + "=" * 80)
print(" MATLAB FUNCTION ")
print("=" * 80)
print(matlab_function)

with open("Predict_Temperature_Transient_AI.m", "w", encoding="utf-8") as f:
    f.write(matlab_function)

print("\nSaved MATLAB function to: Predict_Temperature_Transient_AI.m")