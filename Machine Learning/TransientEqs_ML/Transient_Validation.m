clc, clear, close all
addpath(genpath('../../IGA/Functions'));

% ====================================================================
%  TRANSIENT VALIDATION AGAINST EXACT NODES
% ====================================================================

% 1. PARAMETERS
params.L   = 0.5;
params.R   = 0.04;
params.kx  = 1003;
params.ky  = 171;
params.rho = 2000;
params.c   = 900;
params.Th  = 273;  % Boundary Cold
params.To  = 273;  % Boundary Cold
params.Ti  = 453;  % Initial Hot

% Time Settings
t_final = 30;
dt      = 0.1;
n_steps = ceil(t_final/dt);

% Mesh Settings
refineCount = 6; % High resolution for accurate validation

fprintf('------------------------------------------------------------\n');
fprintf(' Running Transient Solver (T_final = %.1f s)...\n', t_final);
fprintf('------------------------------------------------------------\n');

% 2. CALL SOLVER FUNCTION
[T_history, nodes, ndof, time_vec] = IGA_Transient_Solver(params, refineCount, dt, n_steps);

% 3. PLOT SNAPSHOTS
snap_times = [5, 10, 15, 20, 25, 30]; 
figure('Name', 'Transient Validation', 'Color', 'w', 'Position', [100 100 1200 350]);

for i = 1:length(snap_times)
    t_target = snap_times(i);
    step = round(t_target/dt);
    if step > n_steps, step = n_steps; end
    if step < 1, step = 1; end
    
    T_snap = T_history(:, step);
    
    subplot(1, 6, i);
    scatter(nodes(:,1), nodes(:,2), 10, T_snap, 'filled');
    title(sprintf('t = %.0f s', time_vec(step)));
    axis equal; axis tight; axis off;
    colormap jet; 
    caxis([params.To, params.Ti]); 
end
colorbar;

% ====================================================================
% 4. VALIDATION AGAINST EXACT POINTS (4 Points)
% ====================================================================
try
    pts_data = load('Exact_Nodes.mat');
    exact_pts = pts_data.final_pts / 100; 
    fprintf('\nLoaded "Exact_Nodes.mat" successfully.\n');
catch
    warning('Could not load Exact_Nodes.mat. Using dummy diagonal points for test.');
end

% Create Interpolator for T at t_final
T_final = T_history(:, end);
F_interp = scatteredInterpolant(nodes(:,1), nodes(:,2), T_final, 'natural');

% Evaluate at the 4 Exact Points
T_computed = F_interp(exact_pts(:,1), exact_pts(:,2));

% 5. PRINT COMPARISON TABLE
fprintf('\n============================================================\n');
fprintf(' VALIDATION RESULTS AT t = %.2f s\n', t_final);
fprintf('============================================================\n');
fprintf('%-6s %-10s %-10s %-15s\n', 'Point', 'X [m]', 'Y [m]', 'T_IGA [K]');
fprintf('------------------------------------------------------------\n');

labels = {'A', 'B', 'C', 'D'};
for i = 1:4
    fprintf('%-6s %-10.4f %-10.4f %-15.4f\n', ...
        labels{i}, exact_pts(i,1), exact_pts(i,2), T_computed(i));
end
fprintf('------------------------------------------------------------\n');
fprintf('Mesh Refinement: %d | Elements: ~%d | DOFs: %d\n', ...
    refineCount, size(nodes,1), ndof);
fprintf('============================================================\n');