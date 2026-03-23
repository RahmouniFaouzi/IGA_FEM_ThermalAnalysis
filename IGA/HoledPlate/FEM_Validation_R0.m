clc, clear, close all
addpath(genpath('../Functions'));

%% ==================================================================
%  FEM: VALIDATION OF ORTHOTROPIC SQUARE PLATE 
%  * Validation is performed by setting R -> 0 (approx 0), reducing 
%    the geometry to a solid square plate.
%  * Results are compared to the analytical series solution by 
%    Bruch & Zyvoloski (1974).
% ====================================================================

% 1. GEOMETRY & MATERIAL
% ----------------------
L  = 0.5;        % Plate Width/Height
R  = 0.001;      % Tiny hole radius to approximate solid plate (R -> 0)
kx = 1003;       % W/m.K
ky = 171;        % W/m.K

% 2. BOUNDARY CONDITIONS
% ----------------------
% Note: For the Bruch & Zyvoloski validation case (solid plate):
% - Top/Sides usually fixed at T0
% - Bottom/Left/Right usually fixed at Tg 
% 
% Bruch & Zyvoloski benchmark for orthotropic squares:
% T(x, L) = T_hole (773)   TOP
% T(x, 0) = T_outer (273)  BOTTOM
% T(0, y) = T_outer (273)  LEFT
% T(L, y) = T_outer (273)  RIGHT
T_hot  = 773;   % Boundary at y = L (Top)
T_cold = 273;   % Boundaries at x=0, x=L, y=0

% 3. TARGET ELEMENT COUNTS
% ------------------------
target_elems = [100];

% Storage
res_elems  = zeros(length(target_elems), 1);
res_T_fem  = zeros(length(target_elems), 4); 
res_err_L2 = zeros(length(target_elems), 1);

% Define Evaluation Points (A, B, C, D) for Validation
pts_eval_ideal = [L/2, L/2+R;      % A (Center)
                  L/2, 5*L/8+R;    % B
                  L/2, 3*L/4+R;    % C
                  L/2, 7*L/8+R];   % D

fprintf('Running Validation Study (R -> 0)...\n');

%% ========================
%  LOOP OVER MESH DENSITIES
% =========================
for i = 1:length(target_elems)
    
    target = target_elems(i);
    
    % Adjust parameters for solid-like mesh
    scale  = sqrt(target/100);
    n_circ = max(8, round(8 * scale)); % Minimal circle nodes for tiny hole
    n_rad  = max(4, round(10 * scale));
    bias   = 1.5;

    % 1. GENERATE MESH
    [TR, node, elem, ndof] = triMeshGen(n_circ, n_rad, bias, R, L);
    num_el       = size(elem, 1);
    res_elems(i) = num_el;
    
    % 2. ASSEMBLY
    K = sparse(ndof, ndof);
    F = zeros(ndof, 1);
    D_mat = [kx 0; 0 ky];

    for e = 1:num_el
        n = elem(e,:);
        x_e = node(n,1); y_e = node(n,2);
        Ae  = 0.5 * det([1 x_e(1) y_e(1); 1 x_e(2) y_e(2); 1 x_e(3) y_e(3)]);
        b_c = [y_e(2)-y_e(3); y_e(3)-y_e(1); y_e(1)-y_e(2)];
        c_c = [x_e(3)-x_e(2); x_e(1)-x_e(3); x_e(2)-x_e(1)];
        B = [b_c'; c_c'] / (2*Ae);
        K(n,n) = K(n,n) + B' * D_mat * B * Ae; 
    end

    % 3. BOUNDARY CONDITIONS (Matching Analytical Case)
    tol = 1e-4;
    n_bottom = find(abs(node(:,2)) < tol);
    n_left   = find(abs(node(:,1)) < tol);
    n_right  = find(abs(node(:,1) - L) < tol);
    n_top    = find(abs(node(:,2) - L) < tol);
    fixed_cold = unique([n_bottom; n_left; n_right]);
    fixed_hot  = n_top;
    all_fixed = [fixed_cold; fixed_hot];
    free      = setdiff(1:ndof, all_fixed);

    T = zeros(ndof, 1);
    T(fixed_cold) = T_cold;
    T(fixed_hot)  = T_hot;

    % 4. SOLVE
    T(free) = K(free,free) \ (F(free) - K(free,all_fixed)*T(all_fixed));

    % 5. EXTRACT RESULTS AT POINTS
    % Find nearest nodes
    for pt = 1:4
        dists = sqrt((node(:,1)-pts_eval_ideal(pt,1)).^2 + (node(:,2)-pts_eval_ideal(pt,2)).^2);
        [~, min_idx] = min(dists);
        res_T_fem(i, pt) = T(min_idx);
    end
    
    if (i == 1)
        figure(1);
        hold on;
        triplot(elem, node(:,1), node(:,2), 'Color', [0.5 0.5 0.5]);
        plot(pts_eval_ideal(:,1), pts_eval_ideal(:,2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
        labels = {'A', 'B', 'C', 'D'};
        % Label Points
        for k = 1:4
            text(pts_eval_ideal(k,1)+1, pts_eval_ideal(k,2), labels{k}, 'FontWeight','bold', 'FontSize',12);
        end
        axis equal; 
        axis([0 L 0 L]);
        title('Mesh & Evaluation Points');
        xlabel('X'); 
        ylabel('Y');
    end
end

%% ============================================================
%  ANALYTICAL SOLUTION (Bruch & Zyvoloski, 1974)
% ============================================================
% Equation for Orthotropic Square Plate:
% T(x,y) = T_cold + (T_hot - T_cold) * Sum( An * sin(n*pi*x/L) * sinh(beta*y)/sinh(beta*L) )
fprintf('\nCalculating Analytical Solution...\n');

T_anal = zeros(4, 1);
epsilon = sqrt(kx/ky); % Orthotropic Scaling Factor
for k = 1:4
    x = pts_eval_ideal(k, 1);
    y = pts_eval_ideal(k, 2);
    
    sum_val = 0;
    for n = 1:2:500
        % Calculate beta_n
        beta_n = (n * pi / L) * sqrt(ky/kx); 
        term_y = (n * pi * y / L) * sqrt(kx/ky);
        term_L = (n * pi * L / L) * sqrt(kx/ky);
        % Fourier Coefficient 
        An = 4 / (n * pi);
        % Ratio of sinh terms
        ratio = exp(term_y - term_L) * (1 - exp(-2*term_y)) / (1 - exp(-2*term_L));
        term = An * sin(n * pi * x / L) * ratio;
        sum_val = sum_val + term;
    end
    T_anal(k) = T_cold + (T_hot - T_cold) * sum_val;
end

%% ============================================================
%  COMPARISON TABLES
% ============================================================

fprintf('\n========================================================================\n');
fprintf(' VALIDATION RESULTS (Analytical vs FEM)\n');
fprintf(' Reference: Bruch & Zyvoloski (Solid Plate Limit)\n');
fprintf('========================================================================\n');
fprintf('%-10s | %-15s | %-15s | %-15s | %-15s\n', 'Method', 'Point A', 'Point B', 'Point C', 'Point D');
fprintf('------------------------------------------------------------------------\n');

% Print Analytical Row
fprintf('%-10s | %-15.4f | %-15.4f | %-15.4f | %-15.4f\n', ...
    'ANALYTICAL', T_anal(1), T_anal(2), T_anal(3), T_anal(4));
fprintf('------------------------------------------------------------------------\n');

% Print FEM Rows
for i = 1:length(target_elems)
    fprintf('FEM (%4d) | %-15.4f | %-15.4f | %-15.4f | %-15.4f\n', ...
        res_elems(i), res_T_fem(i,1), res_T_fem(i,2), res_T_fem(i,3), res_T_fem(i,4));
end
fprintf('------------------------------------------------------------------------\n');

% Calculate % Error for the finest mesh
err_final = abs(res_T_fem(end,:) - T_anal') ./ T_anal' * 100;
fprintf('\nFinal Mesh Error (%%): A=%.2f%%, B=%.2f%%, C=%.2f%%, D=%.2f%%\n', ...
    err_final(1), err_final(2), err_final(3), err_final(4));