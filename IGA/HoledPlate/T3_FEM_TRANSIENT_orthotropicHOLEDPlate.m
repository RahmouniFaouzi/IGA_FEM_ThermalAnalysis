clc, clear, close all
addpath(genpath('../Functions'));

%% ======================================
%  FEM: TRANSIENT ORTHOTROPIC HOLED PLATE 
% =======================================

% 1. GEOMETRY & MATERIAL
L = 0.5;         % Plate Width/Height (m)
R = 0.04;        % Hole Radius (m)
kx = 1003;       % W/m.K
ky = 171;        % W/m.K

rho = 2000;      % Density (kg/m^3)
c   = 900;       % Specific Heat (J/kg.K)

% 2. BOUNDARY & INITIAL CONDITIONS
T_initial = 453; % Initial temp (Hot plate)
T_hole    = 273; % Sudden cooling at hole
T_outer   = 273; % Sudden cooling at edges

% 3. TIME INTEGRATION PARAMETERS
t_final = 30;    % Final time (seconds)
dt      = 0.01;  % Time step size
n_steps = ceil(t_final/dt);

% Snapshots to plot (Seconds)
snap_times = [5, 10, 15, 20, 25, 30];
snap_indices = round(snap_times / dt); 

% 4. TARGET ELEMENT COUNT 
target_elems = [500]; 

res_elems = zeros(length(target_elems), 1);
res_T     = zeros(length(target_elems), 4); 
res_err   = zeros(length(target_elems), 1);

fprintf('Running Transient Study...\n');

%% ========================
%  LOOP OVER MESH DENSITIES
% =========================
for i = 1:length(target_elems)
    
    target = target_elems(i);
    fprintf('  Mesh Step...\n ');
    
    scale  = sqrt(target / (R*20)); 
    n_circ = round(20 * scale);
    n_rad  = round(5 * scale);
    
    % Symmetry adjustments
    n_circ = ceil(n_circ/R)*R;
    if n_circ < 8, n_circ = 8; end 
    if n_rad < 2, n_rad = 2; end
    
    bias = 1.0;

    % 1. GENERATE MESH
    [TR, node, elem, ndof] = triMeshGen(n_circ, n_rad, bias, R, L);
    num_el = size(elem, 1);
    res_elems(i) = num_el;
    
    % 2. ASSEMBLY (K and M)
    K = sparse(ndof, ndof);
    M = sparse(ndof, ndof);
    D = [kx 0; 0 ky]; 

    for e = 1:num_el
        n = elem(e,:);
        x_e = node(n,1); y_e = node(n,2);
        
        Ae = 0.5 * det([1 x_e(1) y_e(1); 1 x_e(2) y_e(2); 1 x_e(3) y_e(3)]);
        
        % Stiffness K
        b_c = [y_e(2)-y_e(3); y_e(3)-y_e(1); y_e(1)-y_e(2)];
        c_c = [x_e(3)-x_e(2); x_e(1)-x_e(3); x_e(2)-x_e(1)];
        B = [b_c'; c_c'] / (2*Ae);
        K(n,n) = K(n,n) + B' * D * B * Ae; 
        
        % Mass M 
        Me = (rho * c * Ae / 12) * [2 1 1; 1 2 1; 1 1 2];
        M(n,n) = M(n,n) + Me;
    end

    % 3. BOUNDARY CONDITIONS & ORPHAN FIX
    dist_c = sqrt((node(:,1)-L/2).^2 + (node(:,2)-L/2).^2);
    n_hole = find(abs(dist_c - R) < 0.01); 
    n_outer = find(node(:,1)<1e-4 | node(:,1)>L-1e-4 | node(:,2)<1e-4 | node(:,2)>L-1e-4);
    fixed_bc = unique([n_hole; n_outer]);
    
    % Detect Nodes
    used_nodes = unique(elem(:));
    all_nodes  = (1:ndof)';
    orphans    = setdiff(all_nodes, used_nodes);
    
    fixed = unique([fixed_bc; orphans]);
    free  = setdiff(1:ndof, fixed);
    
    % 4. INITIAL CONDITION
    T = ones(ndof, 1) * T_initial;
    T(fixed_bc) = T_hole; 
    T(orphans)  = 0; 

    % 5. TIME INTEGRATION
    LHS = M/dt + K;
    LHS_free = LHS(free, free);
    LHS_bc   = LHS(free, fixed);
    M_eff    = M/dt;
    
    [L_f, U_f, P_f, Q_f] = lu(LHS_free);
    
    % SETUP FIGURE
    figure('Name', 'Transient Evolution', 'Color', 'w', 'Position', [50 50 1200 600]);
    fprintf(' Solving & Plotting...\n');
    
    % 6. TIME LOOP
    timer_start = tic;
    plot_idx = 1; 
    
    % Exact Hole Boundary for Plotting
    theta = linspace(0, 2*pi, 200);
    xc = L/2 + R*cos(theta); yc = L/2 + R*sin(theta);
    
    solver_residual = 0;

    for t_step = 1:n_steps
        RHS_vec = M_eff * T;
        BC_contribution = LHS_bc * T(fixed);
        RHS_eff = RHS_vec(free) - BC_contribution;
        
        d1 = P_f * RHS_eff;
        d2 = L_f \ d1;
        d3 = U_f \ d2;
        T(free) = Q_f * d3;
        
        % Calculate Residual at final step
        if t_step == n_steps
             solver_residual = norm(LHS_free * T(free) - RHS_eff);
        end

        % CHECK IF SNAPSHOT TIME
        if ismember(t_step, snap_indices)
            
            subplot(2,3, plot_idx);
            
            trisurf(elem, node(:,1), node(:,2), T, 'EdgeColor', 'none', 'FaceColor', 'interp');
            hold on;
            fill(xc, yc, 'w', 'EdgeColor', 'k', 'LineWidth', 1); 
            
            view(2); axis equal; axis([0 L 0 L]);
            colormap jet; 
            caxis([273 453]); 
            
            % Add colorbar to every plot
            c = colorbar; c.FontSize = 8;
            
            current_time = t_step * dt;
            title(['t = ' num2str(current_time) ' s']);
            xlabel('X'); ylabel('Y');
            
            drawnow; 
            plot_idx = plot_idx + 1;
        end
    end
    fprintf('  Done in %.2fs\n', toc(timer_start));

    % 7. EXTRACT FINAL RESULTS (t=30s)
    pts       = load('Exact_Nodes.mat');
    ideal_pts = pts.final_pts/100;
    for pt = 1:4
        dists = sqrt((node(:,1)-ideal_pts(pt,1)).^2 + (node(:,2)-ideal_pts(pt,2)).^2);
        [~, min_idx] = min(dists);
        res_T(i, pt) = T(min_idx);
    end
    
    % 8. SAVE ERROR
    res_err(i) = solver_residual; 
end

%% ============================================================
%  RESULTS TABLE (t = 30s)
% ============================================================
fprintf('\n=========================================================================\n');
fprintf(' TRANSIENT RESULTS (t = %.2f s)\n', t_final);
fprintf('=========================================================================\n');
fprintf('%-10s %-10s %-10s %-10s %-10s %-15s\n', 'Elements', 'A', 'B', 'C', 'D', 'Valid Err');
fprintf('-------------------------------------------------------------------------\n');
for i = 1:length(target_elems)
    fprintf('%-10d %-10.2f %-10.2f %-10.2f %-10.2f %-15.2e\n', ...
        res_elems(i), res_T(i,1), res_T(i,2), res_T(i,3), res_T(i,4), res_err(i));
end
fprintf('-------------------------------------------------------------------------\n');