clc; clear; close all;
addpath(genpath('../Functions'));

%% ============================================================
%  FEM: ORTHOTROPIC SQUARE PLATE WITH CENTRAL HOLE
%  (Convergence Study)
% ============================================================

% 1. GEOMETRY & MATERIAL
% ----------------------
L = 50;          % Plate Width/Height
R = 4;           % Hole Radius
kx = 1003;       % W/m.K
ky = 171;        % W/m.K

% 2. BOUNDARY CONDITIONS
% ----------------------
T_hole  = 773;   % Inner Circle Temp (K)
T_outer = 273;   % Outer Edge Temp (K)

% 3. TARGET ELEMENT COUNTS
% ------------------------
% We adjust mesh parameters (n_circ, n_rad) to approximate these counts
target_elems = [50000];

% Storage for Results
res_elems = zeros(length(target_elems), 1);
res_T     = zeros(length(target_elems), 4); % Columns for A, B, C, D
res_err   = zeros(length(target_elems), 1);
res_valid = cell(length(target_elems), 1);

fprintf('Running Convergence Study...\n');

%% ========================
%  LOOP OVER MESH DENSITIES
% =========================
for i = 1:length(target_elems)
    
    target = target_elems(i);
    
    % Approximate n_circ and n_rad to get target elements
    scale  = sqrt(target / (R*20)); 
    n_circ = round(20 * scale);
    n_rad  = round(5 * scale);
    
    % Ensure n_circ is multiple of 8 for symmetry
    n_circ = ceil(n_circ/R)*R;
    if n_circ < R
        n_circ = R; 
    end
    if n_rad < R/4 
        n_rad = R/4; 
    end
    
    bias = 1.5;

    % 1. GENERATE MESH
    [TR, node, elem, ndof] = triMeshGen(n_circ, n_rad, bias, R, L);
    num_el = size(elem, 1);
    res_elems(i) = num_el;
    
    % 2. ASSEMBLY
    K = sparse(ndof, ndof);
    F = zeros(ndof, 1);
    D = [kx 0; 0 ky];

    for e = 1:num_el
        n = elem(e,:);
        x_e = node(n,1); y_e = node(n,2);
        Ae = 0.5 * det([1 x_e(1) y_e(1); 1 x_e(2) y_e(2); 1 x_e(3) y_e(3)]);
        b_c = [y_e(2)-y_e(3); y_e(3)-y_e(1); y_e(1)-y_e(2)];
        c_c = [x_e(3)-x_e(2); x_e(1)-x_e(3); x_e(2)-x_e(1)];
        B = [b_c'; c_c'] / (2*Ae);
        K(n,n) = K(n,n) + B' * D * B * Ae; 
    end

    % 3. BOUNDARY CONDITIONS
    dist_c = sqrt((node(:,1)-L/2).^2 + (node(:,2)-L/2).^2);
    n_hole = find(abs(dist_c - R) < 0.1);
    n_outer = find(node(:,1)<1e-4 | node(:,1)>L-1e-4 | node(:,2)<1e-4 | node(:,2)>L-1e-4);

    fixed = unique([n_hole; n_outer]);
    free  = setdiff(1:ndof, fixed);

    T = zeros(ndof, 1);
    T(n_hole)  = T_hole;
    T(n_outer) = T_outer;

    % 4. SOLVE
    T(free) = K(free,free) \ (F(free) - K(free,fixed)*T(fixed));

    % 5. EXACT POINT EVALUATION
    % -----------------------------------------
    ideal_pts = [L/2, L/2 + R;      % A (Hole Edge)
                 L/2, 5*L/8 + R;    % B
                 L/2, 3*L/4 + R;    % C
                 L/2, 7*L/8 + R];   % D
             
    % Find EXACT existing nodes closest to these ideal spots
    pts_eval = zeros(4, 2);
    
    for pt = 1:4
        % Distance to all nodes
        dists = sqrt((node(:,1)-ideal_pts(pt,1)).^2 + (node(:,2)-ideal_pts(pt,2)).^2);
        [min_dist, min_idx] = min(dists);
        
        % UPDATE pts_eval to be the REAL node coordinates
        pts_eval(pt, :) = node(min_idx, :);
    end
    % Store the found coordinates for printing later
    res_coords(i, :, :) = pts_eval;
    
    % =========================================================
    % 6. EXACT SEARCH (Now Guaranteed to Find the Node)
    % =========================================================
    tol = 1e-8; % Very strict tolerance
    
    for pt_idx = 1:4
        XP = pts_eval(pt_idx, :);
        
        % Strictly search for the node with these EXACT coordinates
        idXP = find( abs(node(:,1)-XP(1)) < tol & ...
                     abs(node(:,2)-XP(2)) < tol );

        if isempty(idXP)
            res_T(i, pt_idx) = NaN; % Should never happen now
        elseif length(idXP) > 1
            res_T(i, pt_idx) = T(idXP(1));
        else
            res_T(i, pt_idx) = T(idXP(1));
        end
    end
    
    % 6. VALIDATION f(Energy Balance)
    R_flux = K * T;
    Q_in   = sum(R_flux(n_hole));
    Q_out  = sum(R_flux(n_outer));
    
    Balance_Error = abs(Q_in + Q_out); 
    
    res_err(i) = Balance_Error;
    if Balance_Error < 1e-6
        res_valid{i} = 'Pass';
    else
        res_valid{i} = 'Fail';
    end
    if (i == 1)
        plot_results(pts_eval, elem, node, L, T)
    end
end

%% ============================================================
%  OUTPUT TABLES
% ============================================================
% 1. COORDINATES TABLE
fprintf('\n======================================================================================\n');
fprintf(' COORDINATES OF NODES FOUND (x, y)\n');
fprintf('======================================================================================\n');
fprintf('%-10s %-18s %-18s %-18s %-18s\n', 'Elements', 'A (x,y)', 'B (x,y)', 'C (x,y)', 'D (x,y)');
fprintf('--------------------------------------------------------------------------------------\n');
for i = 1:length(target_elems)
    % Format strings for (x,y)
    strA = sprintf('(%.2f, %.2f)', res_coords(i,1,1), res_coords(i,1,2));
    strB = sprintf('(%.2f, %.2f)', res_coords(i,2,1), res_coords(i,2,2));
    strC = sprintf('(%.2f, %.2f)', res_coords(i,3,1), res_coords(i,3,2));
    strD = sprintf('(%.2f, %.2f)', res_coords(i,4,1), res_coords(i,4,2));
    
    fprintf('%-10d %-18s %-18s %-18s %-18s\n', res_elems(i), strA, strB, strC, strD);
    if i == length(target_elems)
        final_pts = squeeze(res_coords(i, :, :)); 
        save('Exact_Nodes.mat', 'final_pts');
    end
end
fprintf('--------------------------------------------------------------------------------------\n');

fprintf('\n=============================================================\n');
fprintf(' RESULTS TABLE\n');
fprintf('=============================================================\n');
fprintf('%-20s %-10s %-10s %-10s %-10s\n', 'Number of Elements', 'A', 'B', 'C', 'D');
fprintf('-------------------------------------------------------------\n');
for i = 1:length(target_elems)
    fprintf('%-20d %-10.2f %-10.2f %-10.2f %-10.2f\n', ...
        res_elems(i), res_T(i,1), res_T(i,2), res_T(i,3), res_T(i,4));
end
fprintf('-------------------------------------------------------------\n');

fprintf('\n=============================================================\n');
fprintf(' VALIDATION TABLE\n');
fprintf('=============================================================\n');
fprintf('%-20s %-15s %-10s\n', 'Number of Elements', 'Error (W)', 'Valid');
fprintf('-------------------------------------------------------------\n');
for i = 1:length(target_elems)
    fprintf('%-20d %-15.2e %-10s\n', ...
        res_elems(i), res_err(i), res_valid{i});
end
fprintf('-------------------------------------------------------------\n');