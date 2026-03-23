import pandas as pd
import numpy  as np
from   pysr   import PySRRegressor
import os

print("="*30)
print(" START SYMBOLIC REGRESSION ")
print("="*30)

# 1. LOAD DATA
filename = 'IGA_DATABASE.csv'
if not os.path.exists(filename):
    print(f"Error: '{filename}' not found. Please ensure it is in the same folder.")
    exit()

df = pd.read_csv(filename)
df_train = df.dropna().copy()

# 2. FEATURE ENGINEERING
L = df_train['L']
R = df_train['R']
cx, cy = L / 2, L / 2
dx, dy = df_train['x'] - cx, df_train['y'] - cy
r_dist = np.sqrt(dx**2 + dy**2)

# Define variables for the AI
df_train['rho']   = r_dist / R                     # Normalized Radial Distance
df_train['theta'] = np.arctan2(dy, dx)             # Angular coordinate
df_train['RL']    = R / L                          # Geometry ratio
df_train['kxky']  = np.sqrt(df_train['kx'] / df_train['ky']) # Anisotropy ratio

# Normalize Target 
df_train['T_norm'] = (df_train['T'] - df_train['To']) / (df_train['Th'] - df_train['To'])

X = df_train[['rho', 'theta', 'kxky', 'RL']]
y = df_train['T_norm']

# Check if kx ky is constant 
if df_train['kxky'].std() < 1e-9:
    print("CRITICAL WARNING: 'kxky' (kx/ky) is constant in your CSV. PySR cannot learn a constant.")

# 3. CONFIGURE PySR 
model = PySRRegressor(
    niterations = 500,
    binary_operators = ["+", "-", "*", "/"],
    unary_operators = [
        "inv(x) = 1/x", 
        "square", 
        "exp",    
        "cos",
        "sin",
        "sqrt"
    ],
    # --- PHYSICAL CONSTRAINTS ---
    # Prevents "cos(cos(exp(x)))"
    nested_constraints = {
        "cos": {"cos": 0, "exp": 0, "square": 1},
        "sin": {"sin": 0, "exp": 0, "square": 1},
        "exp": {"exp": 0, "cos": 1, "square": 1},
        "sqrt": {"sqrt": 0, "exp": 0}
    },
    
    # --- VARIABLE INCENTIVES ---
    complexity_of_constants = 10.0,
    complexity_of_variables = 0.01,
    
    # General simplicity 
    parsimony = 0.001,
    maxsize = 50, 
    
    extra_sympy_mappings = {"inv": lambda x: 1/x},
    loss = "loss(prediction, target) = (prediction - target)^2",
    model_selection = "best", 
)

# 4. RUN DISCOVERY
print("Starting Discovery...")
model.fit(X, y)

# 5. GENERATE FINAL MATLAB CODE
print("\n" + "="*80)
print(" BEST ANALYTICAL EQUATION ")
print("="*80)

best_eq = model.sympy()
eq_str = str(best_eq)

print(f" Equation: T_norm = {eq_str}")

# Map internal AI variables back to MATLAB geometry
eq_str = eq_str.replace("rho", "((sqrt((x-L./2).^2 + (y-L./2).^2)) ./ R)")
eq_str = eq_str.replace("theta", "atan2(y-L./2, x-L./2)")
eq_str = eq_str.replace("RL", "(R./L)")
eq_str = eq_str.replace("kxky", "sqrt(kx./ky)") 

# Python to MATLAB Syntax conversions
eq_str = eq_str.replace("**", ".^") 
eq_str = eq_str.replace("*", ".*") 
eq_str = eq_str.replace("/", "./") 

matlab_function = f"""
function T = Predict_Temperature_AI(x, y, L, R, kx, ky, Th, To)
    % Analytical Physics Equation (Discovered by AI)
    % Parameters:
    % x, y  : Coordinate Matrices
    % L     : Plate Dimension
    % R     : Hole Radius
    % kx, ky: Thermal Conductivities
    % Th, To: Temperature Boundaries

    % 1. Evaluate Dimensionless Shape Function
    % Discovered Equation: {str(best_eq)}
    f_shape = {eq_str};

    % 2. Dimensional Scaling
    T = To + (Th - To) .* f_shape;
    
    % 3. Physical Clamping
    % Ensures T is never outside [To, Th] due to numerical noise
    T = max(min(To, Th), min(max(To, Th), T));
    
    % 4. Geometric Masking (Invalid regions = NaN)
    cx = L./2; cy = L./2;
    dist = sqrt((x-cx).^2 + (y-cy).^2);
    is_valid = (x >= 0) & (x <= L) & (y >= 0) & (y <= L) & (dist >= R);
    
    T(~is_valid) = NaN;
end
"""
print(matlab_function)