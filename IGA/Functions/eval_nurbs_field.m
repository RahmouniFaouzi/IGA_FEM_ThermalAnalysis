function [val, phys] = eval_nurbs_field(u, v, patch, pattern, T_global)
% EVAL_NURBS_FIELD Evaluates the Solution and Geometry at (u,v)
%
%   This function performs the mapping from the parametric space (u,v)
%   to the physical space (x,y) and interpolates the temperature field.
%
%   Inputs:
%       u, v      : Parametric coordinates [0,1]
%       patch     : The NURBS patch structure
%       pattern   : The Global Node ID matrix for this patch
%       T_global  : The full solution vector (Temperature)

    % 1. FIND ACTIVE ELEMENT
    % NURBS basis functions are only non-zero over a small span of knots.
    % We first identify which "element" (knot span) the point (u,v) is inside.
    % 
    spanU = FindSpan_(length(patch.uKnot)-patch.p-2, patch.p, u, patch.uKnot);
    spanV = FindSpan_(length(patch.vKnot)-patch.q-2, patch.q, v, patch.vKnot);
    
    % Calculate the starting indices of the Control Points affecting this span
    % For a degree p, (p+1) control points influence the span.
    idx_u_start = spanU - patch.p; 
    idx_v_start = spanV - patch.q;
    
    % 2. GATHER LOCAL DATA
    % We extract only the relevant Control Points, Weights, and Temperature 
    % values that influence this specific point.
    weights_local = []; 
    T_local = []; 
    CPs_local = [];
    
    count = 1;
    % Loop over the supporting control points (v direction outer, u inner)
    for jj = 0:patch.q
        for ii = 0:patch.p
             % Calculate Linear Index in the patch's control point list
             % Note: MATLAB matrices are column-major (down rows first)
             % Row index = v_start + jj
             % Col index = u_start + ii
             lin_idx = (idx_v_start + jj) + (idx_u_start + ii - 1)*patch.noPtsY;
             
             % Extract Weight (4th component of Control Point)
             w = patch.controlPts(lin_idx, 4); 
             weights_local(count) = w;
             
             % Extract Global Node ID to look up Solution T
             glob_idx = pattern(idx_v_start + jj, idx_u_start + ii);
             T_local(count) = T_global(glob_idx);
             
             % Extract Physical Coordinates (Project from Homogeneous)
             % x = wx / w, y = wy / w
             cp = patch.controlPts(lin_idx, :); 
             CPs_local = [CPs_local; cp(1)/w, cp(2)/w];
             
             count = count + 1;
        end
    end
    
    % 3. COMPUTE BASIS FUNCTIONS
    % Calculate the Rational Basis Functions (R) at (u,v) using local weights
    [R, ~, ~] = NURBS2DBasisDers([u;v], patch.p, patch.q, patch.uKnot, patch.vKnot, weights_local');
    
    % 4. INTERPOLATE FIELD AND GEOMETRY
    % The value at (u,v) is the weighted sum of the basis functions and the nodal values.
    
    % Temperature: T(u,v) = sum( R_i * T_i )
    val = dot(R, T_local); 
    
    % Physical Position: X(u,v) = sum( R_i * P_i )
    phys = R' * CPs_local;
end