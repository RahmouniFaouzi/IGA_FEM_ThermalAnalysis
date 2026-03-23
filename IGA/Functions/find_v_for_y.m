function v_sol = find_v_for_y(v_guess, u_fixed, y_target, patch)
% Newton-Raphson to find v such that Physical_Y(u_fixed, v) = y_target
% --------------------------------------------------------------------
tol = 1e-6;
max_iter = 20;
v = v_guess;

for iter = 1:max_iter
    % 1. Evaluate Position and Derivatives at current (u,v)
    spanU = FindSpan_(length(patch.uKnot)-patch.p-2, patch.p, u_fixed, patch.uKnot);
    spanV = FindSpan_(length(patch.vKnot)-patch.q-2, patch.q, v, patch.vKnot);
    
    % Extract Local Data
    idx_u_start = spanU - patch.p; 
    idx_v_start = spanV - patch.q;
    weights_local = []; 
    CPs_local = [];
    count = 1;
    for jj = 0:patch.q
        for ii = 0:patch.p
            lin_idx = (idx_v_start + jj) + (idx_u_start + ii - 1)*patch.noPtsY;
            cp = patch.controlPts(lin_idx, :);
            weights_local(count) = cp(4);
            CPs_local = [CPs_local; cp(1)/cp(4), cp(2)/cp(4)];
            count = count + 1;
        end
    end
    
    % Basis & Derivatives
    [R, ~, dRdeta] = NURBS2DBasisDers([u_fixed; v], patch.p, patch.q, patch.uKnot, patch.vKnot, weights_local');
    
    % Physical Y and dY/dv
    y_curr = dot(R, CPs_local(:,2));
    dy_dv  = dot(dRdeta, CPs_local(:,2));
    
    % Residual
    res = y_curr - y_target;
    
    if abs(res) < tol
        v_sol = v; return;
    end
    
    % Update
    v = v - res / dy_dv;
    
    % Clamp to domain
    if v < 0, v = 0; end
    if v > 1, v = 1; end
end
v_sol = v;
end

