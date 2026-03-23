clc, clear, close all
addpath(genpath('../Functions'));
% ====================================================================
%  NURBS: TRANSIENT ORTHOTROPIC HOLED PLATE
%  Author : Dr. Rahmouni Faouzi  | Email : rahmounifaouzi01@gmail.com
%           Pr. Khennane Amar    | Email : a.khennane@adfa.edu.au
% ====================================================================

% ------------------------
% 1. SIMULATION PARAMETERS
% ------------------------
L   = 0.5;      % Plate Width/Height
R   = 0.04;     % Hole Radius
kx  = 1003;     % W/m.K
ky  = 171;      % W/m.K
rho = 2000;     % Density
c   = 900;      % Specific Heat

% 2. TIME INTEGRATION PARAMETERS
t_final = 30;        % Final time (seconds)
dt      = 0.05;      % Time step size
n_steps = ceil(t_final/dt);
snap_times = [5, 10, 15, 20, 25, 30];

% 3. BOUNDARY & INITIAL CONDITIONS
T_initial = 453;     % Initial temp (Hot plate)
T_hole    = 273;     % Fixed Temperature at the hole
T_outer   = 273;     % Fixed Temperature at outer edges

% Define Evaluation Points (A, B, C, D)
pts       = load('Exact_Nodes.mat');
ideal_pts = pts.final_pts/100;

fprintf('Running Transient Study...\n');

% -------------------------------------------------------------------------
% 4. IGA MESH GENERATION
% -------------------------------------------------------------------------
refineCount = 6;

fprintf('  Mesh Step (Refinement Level %d)...\n', refineCount);
p = 2; q = 2;  % (Don't Adjust - Degree Fixed)
patches = holedPlateGeometry(L, R, p, q, refineCount);
[node, global_node_patterns, ndof] = mergePatches(patches);

res_elems = 0;
for i = 1:4
    res_elems = res_elems + size(patches(i).element,1);
end

% -------------------------------------------------------------------------
% 5. ASSEMBLY
% -------------------------------------------------------------------------
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
        
        % P.controlPts stores [wx, wy, 0, w]. We need [x, y] = [wx/w, wy/w].
        wts = P.controlPts(local_sctr, 4);
        CP_element_phys = bsxfun(@rdivide, P.controlPts(local_sctr, 1:2), wts);
        
        idx_u = P.index(e, 1);
        idx_v = P.index(e, 2);
        Xi = P.spans_u(idx_u, :);
        Eta = P.spans_v(idx_v, :);
        
        Ke = zeros(length(glob_sctr));
        Me = zeros(length(glob_sctr));
        
        for j = 1:length(Q_gauss)
            for i = 1:length(Q_gauss)
                % Map Gauss Point to Parametric Space
                xi  = parent2ParametricSpace(Xi, Q_gauss(i,1));
                eta = parent2ParametricSpace(Eta, Q_gauss(j,1));
                J2  = jacobianPaPaMapping(Xi,Eta);
                
                % Evaluate NURBS Basis & Derivatives
                [R_n, dR_dxi, dR_deta] = NURBS2DBasisDers([xi; eta], p, q, P.uKnot, P.vKnot, wts);
                
                % Jacobian (Parametric -> Physical) uses Physical Points
                Jac = [dR_dxi'; dR_deta'] * CP_element_phys;
                detJ = det(Jac);
                
                % Physical Derivatives [dR/dx] = inv(J) * [dR/dxi]
                dRdx = (Jac \ [dR_dxi'; dR_deta'])';
                B = dRdx';
                
                % local Stiff & Mass arrays
                Ke = Ke + B' * D * B * detJ * W(i) * W(j) * J2;
                Me = Me + R_n * R_n' * rho * c * detJ * W(i) * W(j) * J2;
            end
        end
        K(glob_sctr, glob_sctr) = K(glob_sctr, glob_sctr) + Ke;
        M(glob_sctr, glob_sctr) = M(glob_sctr, glob_sctr) + Me;
    end
end

% -------------------------------------------------------------------------
% 6. BOUNDARY CONDITIONS
% -------------------------------------------------------------------------
cx = L/2;
cy = L/2;
rad_dist = sqrt((node(:,1)-cx).^2 + (node(:,2)-cy).^2);
tol = 1e-4;
n_hole  = find(abs(rad_dist - R) < 0.01);
n_outer = find(node(:,1)<tol | node(:,1)>L-tol | node(:,2)<tol | node(:,2)>L-tol);

fixed = unique([n_hole; n_outer]);
free  = setdiff(1:ndof, fixed);

T = ones(ndof, 1) * T_initial;
T(n_hole)  = T_hole;
T(n_outer) = T_outer;

% -------------------------------------------------------------------------
% 7. TIME INTEGRATION
% -------------------------------------------------------------------------
LHS = M/dt + K;
M_eff = M/dt;

% Pre-factorize
[L_f, U_f, P_f, Q_f] = lu(LHS(free, free));

fprintf('  Solving & Plotting...\n');
timer_start = tic;

% Setup Figure
figure('Name', 'Transient Evolution', 'Color', 'w', 'Position', [100 100 900 600]);
colormap jet;

curr_time = 0;
solver_residual = 0;

% Visualization Snapshots
plot_idx = 1;

for step = 1:n_steps
    curr_time = curr_time + dt;
    
    RHS_vec = M_eff * T;
    RHS_eff = RHS_vec(free) - LHS(free, fixed) * T(fixed);
    T(free) = Q_f * (U_f \ (L_f \ (P_f * RHS_eff)));
    
    % Check residual at final step
    if step == n_steps
        solver_residual = norm(LHS(free, free) * T(free) - RHS_eff);
    end
    
    % Plotting
    if any(abs(curr_time - snap_times) < dt/2)
        subplot(2, 3, plot_idx);
        PlotTemp(patches, global_node_patterns, T, L, R);
        title(sprintf('t = %.0f s', curr_time));
        caxis([273, 453]);
        axis equal;
        axis off;
        axis([0 L 0 L]);
        plot_idx = plot_idx + 1;
        drawnow;
    end
end

time_taken = toc(timer_start);
fprintf('  Done in %.2fs\n', time_taken);

% -------------------------------------------------------------------------
% 8. EVALUATE EXACT POINTS
% -------------------------------------------------------------------------
res_T = zeros(1, 4);
for i = 1:4
    pt = ideal_pts(i, :);
    found = false;
    val = 273; % Default
    
    for k = 1:4
        [xi, eta, success] = invMapNewton(patches(k), pt);
        if success
            val = evalTempTrans(patches(k), global_node_patterns{k}, T, xi, eta, p, q);
            found = true;
            break;
        end
    end
    res_T(i) = val;
end

% -------------------------------------------------------------------------
% 9. OUTPUT TABLE
% -------------------------------------------------------------------------
fprintf('\n=========================================================================\n');
fprintf(' TRANSIENT RESULTS (t = %.2f s)\n', t_final);
fprintf('=========================================================================\n');
fprintf('%-10s %-10s %-10s %-10s %-10s %-15s\n', 'Elements', 'A', 'B', 'C', 'D', 'Valid Err');
fprintf('-------------------------------------------------------------------------\n');
fprintf('%-10d %-10.2f %-10.2f %-10.2f %-10.2f %-15.2e\n', ...
    res_elems, res_T(1), res_T(2), res_T(3), res_T(4), solver_residual);
fprintf('-------------------------------------------------------------------------\n');