clc, clear, close all
addpath(genpath('../../IGA/Functions'));

fprintf('STARTING FULL GRID DATABASE GENERATION ...\n');

% =========================================================================
% CONFIGURATION
% =========================================================================
num_simulations = 1;
refine_db       = 6;
output_file     = 'IGA_Grid_With_NaNs.csv';

% 1. Grid Density for Material
u_vec = 0.5;
v_vec = [0, 0.274, 0.579, 0.8930];

for sim = 1:num_simulations
    fprintf('Simulation %d / %d ... \n', sim, num_simulations);
    
    % --- Randomize Inputs ---
    L  = 50;
    R  = 4;
    kx = 1003;
    ky = 171;
    Th = 773;
    To = 273;
    
    % --- Solve Physics ---
    [patches, global_patterns, T_sol] = IGA_Solver_Core(L, R, kx, ky, Th, To, refine_db);
    
    sim_block = [];
    
    % MATERIAL POINTS
    p_idx = 1;
    u     = u_vec;
    for v = v_vec
        % Evaluate Exact Physics
        [val, phys] = eval_nurbs_field(u, v, patches(p_idx), global_patterns{p_idx}, T_sol);
        fprintf('%.6f %.6f %.6f\n', phys(1), phys(2), val);
    end
    fprintf('Done\n');
end