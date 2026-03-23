clc, clear, close all
addpath(genpath('../Functions'));

%% ========================================================================
%  PROJECT: HB-IGA THERMAL ANALYSIS OF ORTHOTROPIC PLATE
%  ------------------------------------------------------------------------
%  Subject: Steady-State Heat Conduction (Convergence Study)
%  Method : Isogeometric Analysis (IGA) with Multi-Patch NURBS
%  Feature: Graded Knot Refinement for Local Concentration
%  Authors: Pr. Khennane Amar & Dr. Rahmouni Faouzi
%  Date   : 2026
%  ------------------------------------------------------------------------
%  Description:
%  This script solves the 2D steady heat equation on a plate with a 
%  central hole using orthotropic material properties. It employs a 
%  multi-patch approach with graded mesh refinement to capture high 
%  gradients near the singularity (hole).
% =========================================================================

%% 1. SIMULATION PARAMETERS & CONSTANTS
%  ========================================================================
% -- Geometry --
L = 50;              % Plate Half-Width/Height
R = 4;               % Radius of the central hole

% -- Material Properties (Orthotropic) --
kx = 1003;           % Thermal Conductivity in X [W/m.K]
ky = 171;            % Thermal Conductivity in Y [W/m.K]
D_mat = [kx 0; 
         0 ky];      % Conductivity Matrix

% -- Boundary Conditions (Dirichlet) --
T_hole  = 773;       % Temperature at Inner Boundary (Hole) [K]
T_outer = 273;       % Temperature at Outer Boundary [K]

% -- Mesh Refinement Strategy --
% Level 1: Coarse Mesh (8x8)
% Level 2: Concentrated Mesh (16x16)
% Level 3: Fine Concentrated Mesh (32x32)
refinement_levels = [2]; 
plotLevel = 1;

% -- Results Initialization --
res_elems = zeros(length(refinement_levels), 1); % Element count
res_T     = zeros(length(refinement_levels), 4); % Temp at probes A,B,C,D
res_err   = zeros(length(refinement_levels), 1); % Flux Validation Error

fprintf('------------------------------------------------------------------------\n');
fprintf('  STARTING IGA SIMULATION (Graded Knot Refinement & Multi-Patch NURBS)\n');
fprintf('------------------------------------------------------------------------\n');

%% 2. MAIN CONVERGENCE LOOP
%  ========================================================================
for step = 1:length(refinement_levels)
    level = refinement_levels(step);
    
    % ---------------------------------------------------------------------
    % 2.1 GEOMETRY GENERATION
    % ---------------------------------------------------------------------
    fac = sind(45); 
    a   = fac * R;
    cp  = zeros(4,3,3);
    
    % Row 1: Inner Circular Arc (Hole)
    cp(1:2,1,1) = [a; a];   
    cp(1:2,2,1) = [0; 2*a];   
    cp(1:2,3,1) = [-a; a];
    % Row 3: Outer Square Boundary
    cp(1:2,1,3) = [L/2; L/2]; 
    cp(1:2,2,3) = [0; L/2]; 
    cp(1:2,3,3) = [-L/2; L/2];
    % Row 2: Midpoints (Linear Interpolation)
    cp(1:3,:,2) = (cp(1:3,:,1) + cp(1:3,:,3)) / 2;
    
    % Shift Geometry to Center & Apply Weights
    for r = 1:3 
        for c = 1:3 
            cp(1:2,r,c) = cp(1:2,r,c) + [L/2; L/2];
        end
    end
    cp(4,:,:) = 1; 
    cp(4,2,1) = fac; 
    cp(1:2,2,1) = cp(1:2,2,1) * fac; % Homogeneous pre-multiplication
    
    % Base Knot Vectors (Quadratic)
    U = [0 0 0 1 1 1]; 
    V = [0 0 0 1 1 1]; 
    p = 2; q = 2;
    
    % ---------------------------------------------------------------------
    % 2.2 MESH REFINEMENT (Graded Strategy)
    % ---------------------------------------------------------------------
    for k = 1:(level+1)
        [cp, U, V] = refine_uniform_dir(cp, U, V, p, q, 1); 
    end
    
    % V-Direction: Graded Refinement 
    if level > 0
        n_rad = 6 * 2^(level-1);
        t = linspace(0, 1, n_rad+1);
        knots = t.^3;           % Power law distribution
        knots = knots(2:end-1); % Exclude 0 and 1
        
        for k = 1:length(knots)
            [cp, U, V] = insert_knot_surface(cp, U, V, p, q, knots(k), 2); 
        end
    end
    
    % ---------------------------------------------------------------------
    % 2.3 TOPOLOGY ASSEMBLY (4 Patches)
    % ---------------------------------------------------------------------
    patches = cell(4,1);
    
    % Define Rotation Matrices
    T_cen = [1 0 0 -L/2; 0 1 0 -L/2; 0 0 1 0; 0 0 0 1];
    R_90  = [0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
    T_back= [1 0 0 L/2; 0 1 0 L/2; 0 0 1 0; 0 0 0 1];
    M_rot = T_back * R_90 * T_cen;
    
    % Generate Rotated Patches
    patches{1}.cp = cp;                                         % North
    patches{2}.cp = apply_transform_4D(patches{1}.cp, M_rot);   % West
    patches{3}.cp = apply_transform_4D(patches{2}.cp, M_rot);   % South
    patches{4}.cp = apply_transform_4D(patches{3}.cp, M_rot);   % East
    
    % Process Connectivity & Flattening
    all_CPs = []; 
    patch_offsets = [0];
    
    for k = 1:4
        patches{k}.U = U; 
        patches{k}.V = V; 
        patches{k}.p = p; 
        patches{k}.q = q;
        
        [dim, nu, nv] = size(patches{k}.cp);
        patches{k}.nu = nu; 
        patches{k}.nv = nv;
        
        % Build Element Connectivity Tables
        [patches{k}.elRangeU, patches{k}.elConnU] = Connectivity(p, U);
        [patches{k}.elRangeV, patches{k}.elConnV] = Connectivity(q, V);
        
        % Flatten Control Points [nu*nv, 4]
        flat_cp = zeros(nu*nv, 4); 
        cnt = 1;
        for j = 1:nv
            for i = 1:nu
                flat_cp(cnt,:) = patches{k}.cp(:,i,j)'; 
                cnt = cnt+1; 
            end
        end
        patches{k}.cp_list = flat_cp;
        
        % Extract Physical Coords for Merging
        phys    = flat_cp(:,1:2) ./ repmat(flat_cp(:,4), 1, 2);
        all_CPs = [all_CPs; phys];
        patch_offsets(end+1) = patch_offsets(end) + nu*nv;
    end
    
    % Global Node Merging (Remove Duplicates)
    [unique_nodes, ~, ic] = unique(round(all_CPs, 5), 'rows');
    ndof = size(unique_nodes, 1);
    
    % Map Local IDs to Global IDs
    for k = 1:4
        range = (patch_offsets(k)+1) : patch_offsets(k+1);
        patches{k}.g_map = reshape(ic(range), patches{k}.nu, patches{k}.nv); 
    end
    
    % Report Element Count
    num_elems = 4 * size(patches{1}.elRangeU,1) * size(patches{1}.elRangeV,1);
    res_elems(step) = num_elems;
    fprintf('  > Level %d: %d Elements (Assembling...)\n', level, num_elems);
    
    % ---------------------------------------------------------------------
    % 2.4 STIFFNESS MATRIX ASSEMBLY
    % ---------------------------------------------------------------------
    est_nz = num_elems * 100;
    I_K = zeros(est_nz, 1); 
    J_K = zeros(est_nz, 1); 
    V_K = zeros(est_nz, 1); 
    nz_k = 0;
    
    [W, Q_gauss] = quadrature(3, 'GAUSS', 2);
    
    for k = 1:4
        patch = patches{k}; 
        nu = patch.nu; 
        
        for v_el = 1:size(patch.elRangeV,1)
            for u_el = 1:size(patch.elRangeU,1)
                % Element Geometry
                xiE  = patch.elRangeU(u_el, :); 
                etaE = patch.elRangeV(v_el, :);
                
                % Element Connectivity
                j_idx = patch.elConnV(v_el,:); 
                i_idx = patch.elConnU(u_el,:); 
                glob_sctr = patch.g_map(i_idx, j_idx); 
                glob_sctr = glob_sctr(:); 
                
                [JJ, II] = meshgrid(j_idx, i_idx);
                sctr_lin = sub2ind([nu, patch.nv], II(:), JJ(:));
                pts = patch.cp_list(sctr_lin, :);
                
                % Gauss Integration Loop
                for gp = 1:length(W)
                    wt  = W(gp);
                    Xi  = ((xiE(2)-xiE(1))*Q_gauss(gp,1) + sum(xiE))/2;
                    Eta = ((etaE(2)-etaE(1))*Q_gauss(gp,2) + sum(etaE))/2;
                    J2  = (xiE(2)-xiE(1))*(etaE(2)-etaE(1))/4;
                    
                    % NURBS Basis & Derivatives
                    [~, dRdxi, dRdeta] = NURBS2DBasisDers([Xi; Eta], p, q, U, V, pts(:,4));
                    
                    % Jacobian & Physical Derivatives
                    CP_phys = [pts(:,1)./pts(:,4), pts(:,2)./pts(:,4)];
                    Jac  = CP_phys' * [dRdxi, dRdeta]; 
                    detJ = det(Jac);
                    dRdx = [dRdxi, dRdeta] / Jac;
                    B    = dRdx';
                    
                    % Local Stiffness
                    Ke = B' * D_mat * B * detJ * J2 * wt;
                    
                    % Sparse Assembly
                    n_b = length(glob_sctr); 
                    [c_g, r_g] = meshgrid(glob_sctr, glob_sctr);
                    idx = nz_k+1 : nz_k+n_b^2;
                    I_K(idx) = r_g(:); 
                    J_K(idx) = c_g(:); 
                    V_K(idx) = Ke(:); 
                    nz_k = nz_k + n_b^2;
                end
            end
        end
    end
    K = sparse(I_K(1:nz_k), J_K(1:nz_k), V_K(1:nz_k), ndof, ndof);
    F = zeros(ndof, 1);
    
    % ---------------------------------------------------------------------
    % 2.5 BOUNDARY CONDITIONS & SOLUTION
    % ---------------------------------------------------------------------
    n_hole = []; n_outer = [];
    
    % Topological Detection (V=0 is Hole)
    for k = 1:4
        n_hole = [n_hole; patches{k}.g_map(:, 1)];
        n_outer = [n_outer; patches{k}.g_map(:, end)];
    end
    
    % Physical Detection (Outer Box Corners)
    tol = 1e-4;
    for i = 1:ndof
        pt = unique_nodes(i,:);
        if abs(pt(1)-L) < tol || abs(pt(2)-L) < tol
            n_outer = [n_outer; i]; 
        end
    end
    
    fixed = unique([n_hole; n_outer]); 
    free  = setdiff(1:ndof, fixed);
    
    T = zeros(ndof, 1); 
    T(n_hole)  = T_hole; 
    T(n_outer) = T_outer;
    
    % Solve System
    T(free) = K(free, free) \ (F(free) - K(free, fixed)*T(fixed));
    
    % -- Validation Error (Flux Balance) --
    Q_flux = K * T;
    % Sum flux at boundaries (should ideally sum to 0)
    res_err(step) = abs(sum(Q_flux(unique(n_hole))) + sum(Q_flux(unique(n_outer))));
    
    % ---------------------------------------------------------------------
    % 2.6 EXACT EVALUATION AT PROBES (Newton-Raphson)
    % ---------------------------------------------------------------------
    % Load target physical points generated by the FEM code
    pts        = load('Exact_Nodes.mat'); 
    pts_eval_y = pts.final_pts(:,2)';
    row_res = zeros(1,4); 
    patch = patches{1}; u_fix = 0.5;
    
    for pi = 1:4
        y_target = pts_eval_y(pi);
        v_sol = (y_target - (L/2+R))/(L/2-R); % Initial guess
        
        for iter = 1:10
            spanU = FindSpan_(length(U)-p-2, p, u_fix, U); 
            spanV = FindSpan_(length(V)-q-2, q, v_sol, V);
            idx_u = (spanU-p):spanU; 
            idx_v = (spanV-q):spanV;
            
            glob_sctr = patch.g_map(idx_u, idx_v); 
            glob_sctr = glob_sctr(:);
            
            [JJ, II] = meshgrid(idx_v, idx_u);
            lin_idx  = sub2ind([patch.nu, patch.nv], II(:), JJ(:));
            local_w  = patch.cp_list(lin_idx, 4); 
            local_cp = patch.cp_list(lin_idx, 1:2);
            local_cp(:,1) = local_cp(:,1)./ local_w;
            local_cp(:,2) = local_cp(:,2)./ local_w;
            
            [R_b, ~, dRdv] = NURBS2DBasisDers([u_fix; v_sol], p, q, U, V, local_w);
            phys = R_b' * local_cp;
            
            f = phys(2) - y_target; 
            if abs(f) < 1e-7, break; end
            
            dy_dv = dRdv' * local_cp; 
            v_sol = v_sol - f/dy_dv(2);
        end
        
        [R_b, ~]    = NURBS2DBasisDers([u_fix; v_sol], p, q, U, V, local_w);
        row_res(pi) = dot(R_b, T(glob_sctr));
    end
    res_T(step, :) = row_res;
    
    % Plot section
    if (step == plotLevel)
        plotGRADEDIGA(patches, T, U, V, p, q)
    end
end

%% 3. RESULTS OUTPUT
%  =================
fprintf('\n===================================================================\n');
fprintf(' IGA CONVERGENCE RESULTS (NURBS - Orthotropic Plate)\n');
fprintf('===================================================================\n');
fprintf('%-10s %-10s %-10s %-10s %-10s %-15s\n', 'Elements', 'A', 'B', 'C', 'D', 'Error (W)');
fprintf('-------------------------------------------------------------------\n');
for i = 1:length(res_elems)
    fprintf('%-10d %-10.2f %-10.2f %-10.2f %-10.2f %-15.2e\n', ...
        res_elems(i), res_T(i,1), res_T(i,2), res_T(i,3), res_T(i,4), res_err(i));
end
fprintf('-------------------------------------------------------------------\n');