function [T_hist, nodes, ndof, time_vec] = IGA_Transient_Solver(params, refineCount, dt, n_steps)
    % Unpack Parameters
    L = params.L;
    R = params.R;
    kx = params.kx;
    ky = params.ky;
    rho = params.rho;
    c = params.c;
    Th = params.Th;
    To = params.To;
    Ti = params.Ti;
    
    % 1. MESH GENERATION 
    p = 2; q = 2;
    patches = holedPlateGeometry(L, R, p, q, refineCount);
    [nodes, global_node_patterns, ndof] = mergePatches(patches);
    
    % 2. ASSEMBLY
    K = sparse(ndof, ndof);
    M = sparse(ndof, ndof);
    D = [kx 0; 0 ky];
    
    [Q_gauss, W] = HB_lgwt(p+1, -1, 1);
    
    for pat = 1:4
        P = patches(pat);
        g_map = global_node_patterns{pat};
        
        for e = 1:size(P.element, 1)
            local_sctr = P.element(e, :);
            glob_sctr = g_map(local_sctr);
            
            wts = P.controlPts(local_sctr, 4);
            CP_phys = bsxfun(@rdivide, P.controlPts(local_sctr, 1:2), wts);
            
            idx_u = P.index(e, 1);
            idx_v = P.index(e, 2);
            Xi = P.spans_u(idx_u, :);
            Eta = P.spans_v(idx_v, :);
            
            Ke = zeros(length(glob_sctr));
            Me = zeros(length(glob_sctr));
            
            for j = 1:length(Q_gauss)
                for i = 1:length(Q_gauss)
                    xi  = parent2ParametricSpace(Xi, Q_gauss(i,1));
                    eta = parent2ParametricSpace(Eta, Q_gauss(j,1));
                    J2  = jacobianPaPaMapping(Xi, Eta);
                    
                    [R_n, dR_dxi, dR_deta] = NURBS2DBasisDers([xi; eta], p, q, P.uKnot, P.vKnot, wts);
                    
                    Jac = [dR_dxi'; dR_deta'] * CP_phys;
                    detJ = det(Jac);
                    dRdx = (Jac \ [dR_dxi'; dR_deta'])';
                    B = dRdx';
                    
                    Ke = Ke + B' * D * B * detJ * W(i) * W(j) * J2;
                    Me = Me + R_n * R_n' * rho * c * detJ * W(i) * W(j) * J2;
                end
            end
            K(glob_sctr, glob_sctr) = K(glob_sctr, glob_sctr) + Ke;
            M(glob_sctr, glob_sctr) = M(glob_sctr, glob_sctr) + Me;
        end
    end
    
    % 3. BOUNDARY CONDITIONS
    cx = L/2; cy = L/2;
    dist = sqrt((nodes(:,1)-cx).^2 + (nodes(:,2)-cy).^2);
    tol = 1e-4;
    
    n_hole  = find(abs(dist - R) < 0.01);
    n_outer = find(nodes(:,1)<tol | nodes(:,1)>L-tol | nodes(:,2)<tol | nodes(:,2)>L-tol);
    
    fixed = unique([n_hole; n_outer]);
    free  = setdiff(1:ndof, fixed);
    
    % Initial Condition
    T = ones(ndof, 1) * Ti;
    T(n_hole)  = Th;
    T(n_outer) = To;
    
    % 4. TIME INTEGRATION
    LHS = M/dt + K;
    M_eff = M/dt;
    
    % Pre-factorize
    [L_f, U_f, P_f, Q_f] = lu(LHS(free, free));
    
    % Storage
    T_hist = zeros(ndof, n_steps);
    time_vec = zeros(n_steps, 1);
    
    for step = 1:n_steps
        RHS_vec = M_eff * T;
        RHS_eff = RHS_vec(free) - LHS(free, fixed) * T(fixed);
        T(free) = Q_f * (U_f \ (L_f \ (P_f * RHS_eff)));
        
        T_hist(:, step) = T;
        time_vec(step) = step * dt;
    end
end