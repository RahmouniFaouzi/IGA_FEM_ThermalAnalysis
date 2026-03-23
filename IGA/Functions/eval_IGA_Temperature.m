function T_val = eval_IGA_Temperature(Pt, T_vec, Surf)
% EVALIGA_TEMPERATURE Evaluates Temperature at physical point Pt(x,y)
%   Inputs:
%      Pt    - [x, y] coordinates
%      T_vec - Solved Temperature Vector (Global DOFs)
%      Surf  - NURBS Geometry Structure

% 1. EXTRACT DATA FROM SURF
% -------------------------
uKnot = Surf.KntVect{1};
vKnot = Surf.KntVect{2};
p = Surf.Order(1);
q = Surf.Order(2);
weights = Surf.Weights;

% Calculate Number of Control Points in X (needed for indexing)
noPtsX = Surf.NCtrlPts(1);
noPtsY = Surf.NCtrlPts(1);

% 2. GET COORDINATES
% ------------------
u = Pt(1);
v = Pt(2);

% 3. FIND KNOT SPANS 
% ----------------------------------
uSpan = FindSpan(noPtsX, p, u, uKnot);
vSpan = FindSpan(noPtsY, q, v, vKnot);

% 4. EVALUATE BASIS FUNCTIONS
% ---------------------------
Nu = BasisFuns(uSpan, u, p, uKnot); % 
Nv = BasisFuns(vSpan, v, q, vKnot);

% 5. COMPUTE TEMPERATURE (Rational Sum)
% -------------------------------------
T_val = 0;
w_sum = 0;

for j = 0:q
    for i = 0:p
        % Convert 2D tensor index to 1D Global Index
        ix = uSpan - p + i;
        iy = vSpan - q + j;
        
        % Global Index
        global_idx = (iy-1) * noPtsX + ix;
        
        % Calculate Rational Basis
        val = Nu(i+1) * Nv(j+1) * weights(global_idx);
        
        T_val = T_val + val * T_vec(global_idx);
        w_sum = w_sum + val;
    end
end

% Normalize by weight sum (NURBS definition)
if w_sum ~= 0
    T_val = T_val / w_sum;
end

end