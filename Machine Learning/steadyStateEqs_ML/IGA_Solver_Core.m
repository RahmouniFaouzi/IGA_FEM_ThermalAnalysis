function [patches, global_node_patterns, T_sol] = IGA_Solver_Core(L, R, kx, ky, T_hole, T_outer, refine_level)
    % IGA_SOLVER_CORE
    % Returns the solved geometry (patches) and temperature vector (T_sol).
    
    % Material
    D_mat = [kx 0; 0 ky];
    
    % --- 1. GEOMETRY GENERATION ---
    fac = sind(45);
    a   = fac * R;
    controlPts = zeros(4,3,3);
    
    % Inner Arc
    controlPts(1:2,1,1) = [a; a];
    controlPts(1:2,2,1) = [0; 2*a];
    controlPts(1:2,3,1) = [-a; a];
    
    % Outer Square
    controlPts(1:2,1,3) = [L/2; L/2];
    controlPts(1:2,2,3) = [0; L/2];
    controlPts(1:2,3,3) = [-L/2; L/2];
    
    % Mid-points
    controlPts(1:3,:,2) = (controlPts(1:3,:,1) + controlPts(1:3,:,3)) / 2;
    
    % Shift to Center
    for r = 1:3, for c = 1:3, controlPts(1:2,r,c) = controlPts(1:2,r,c) + [L/2; L/2]; end, end
    
    % Weights
    controlPts(4,:,:) = 1;
    controlPts(4,2,1) = fac;
    controlPts(1:2,2,1) = controlPts(1:2,2,1) * fac;
    
    uKnot = [0 0 0 1 1 1]; vKnot = [0 0 0 1 1 1];
    solid1 = nrbmak(controlPts, {uKnot, vKnot});
    
    % Rotate (4 Patches)
    T_cen = [1 0 0 -L/2; 0 1 0 -L/2; 0 0 1 0; 0 0 0 1];
    R_90  = [0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
    T_back= [1 0 0 L/2; 0 1 0 L/2; 0 0 1 0; 0 0 0 1];
    M_trans = T_back * R_90 * T_cen;
    
    solid2 = nrbtform(solid1, M_trans);
    solid3 = nrbtform(solid2, M_trans);
    solid4 = nrbtform(solid3, M_trans);
    
    % Refinement
    if refine_level > 0
        solid1 = hRefineNURBS(solid1, refine_level);
        solid2 = hRefineNURBS(solid2, refine_level);
        solid3 = hRefineNURBS(solid3, refine_level);
        solid4 = hRefineNURBS(solid4, refine_level);
    end
    
    % Convert
    patches(1) = convert2DNurbsToPatch(solid1);
    patches(2) = convert2DNurbsToPatch(solid2);
    patches(3) = convert2DNurbsToPatch(solid3);
    patches(4) = convert2DNurbsToPatch(solid4);
    
    % Need to attach the raw struct for evaluation functions later
    patches(1).nurbs_struct = solid1; 
    patches(2).nurbs_struct = solid2;
    patches(3).nurbs_struct = solid3; 
    patches(4).nurbs_struct = solid4;
    
    % --- 2. CONNECTIVITY ---
    all_CPs = [];
    for p=1:4
        cp = patches(p).controlPts;
        cp_geom = [cp(:,1)./cp(:,4), cp(:,2)./cp(:,4)];
        all_CPs = [all_CPs; cp_geom];
    end
    [~, ~, ic] = unique(round(all_CPs, 5), 'rows');
    ndof = max(ic);
    
    global_node_patterns = cell(4,1);
    curr = 0;
    for p=1:4
        n_p = size(patches(p).controlPts,1);
        indices_local = (1:n_p)' + curr;
        global_indices = ic(indices_local);
        global_node_patterns{p} = reshape(global_indices, patches(p).noPtsY, patches(p).noPtsX);
        curr = curr + n_p;
    end
    
    % --- 3. ASSEMBLY ---
    K = sparse(ndof, ndof);
    F = zeros(ndof, 1);
    [W, Q_gauss] = quadrature(3, 'GAUSS', 2);
    
    for p = 1:4
        patch = patches(p);
        for e = 1:size(patch.element, 1)
            local_sctr = patch.element(e, :);
            glob_sctr  = global_node_patterns{p}(local_sctr);
            pts        = patch.controlPts(local_sctr, :);
            
            idx_u = patch.index(e,1); idx_v = patch.index(e,2);
            xiE = patch.elRangeU(idx_u, :); etaE = patch.elRangeV(idx_v, :);
            
            for gp = 1:length(W)
                wt = W(gp);
                Xi  = parent2ParametricSpace(xiE, Q_gauss(gp,1));
                Eta = parent2ParametricSpace(etaE, Q_gauss(gp,2));
                J2  = jacobianPaPaMapping(xiE,etaE);
                
                % Call your library function
                [~, dRdxi, dRdeta] = NURBS2DBasisDers([Xi; Eta], patch.p, patch.q, patch.uKnot, patch.vKnot, pts(:,4));
                
                CP_phys = [pts(:,1)./pts(:,4), pts(:,2)./pts(:,4)];
                Jac = CP_phys' * [dRdxi, dRdeta];
                detJ = det(Jac);
                dRdx = [dRdxi, dRdeta] / Jac; 
                B = dRdx';
                
                K(glob_sctr, glob_sctr) = K(glob_sctr, glob_sctr) + B' * D_mat * B * detJ * J2 * wt;
            end
        end
    end
    
    % --- 4. BOUNDARY CONDITIONS ---
    n_hole = []; n_outer = [];
    for p = 1:4
        hole_nodes_p  = global_node_patterns{p}(1, :);
        outer_nodes_p = global_node_patterns{p}(end, :);
        n_hole  = [n_hole; hole_nodes_p(:)];
        n_outer = [n_outer; outer_nodes_p(:)];
    end
    n_hole = unique(n_hole);
    n_outer = unique(n_outer);
    
    fixed = unique([n_hole; n_outer]);
    free  = setdiff(1:ndof, fixed);
    
    T_sol = zeros(ndof, 1);
    T_sol(n_hole)  = T_hole;
    T_sol(n_outer) = T_outer;
    
    T_sol(free) = K(free, free) \ (F(free) - K(free, fixed)*T_sol(fixed));
end