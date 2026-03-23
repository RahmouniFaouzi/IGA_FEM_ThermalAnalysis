clear, clc, close all
addpath(genpath('../../IGA/Functions'));

fprintf('STARTING TRANSIENT DATABASE GENERATION...\n');

% =========================================================================
% CONFIGURATION
% =========================================================================
num_simulations = 100;  % Number of different plate geometries/materials
refine_db       = 6;    % Medium refinement for speed 
output_file     = 'IGA_Transient_Database2.csv';
points_per_sim  = 2000;  % Random points to harvest per simulation

% Time Settings for Database
t_final_db = 40;
dt_db      = 0.1;        % Sampling time step
n_steps_db = ceil(t_final_db/dt_db);

Full_Data = [];

for sim = 1:num_simulations
    fprintf('Simulation %d / %d ... \n', sim, num_simulations);
    
    % --- 1. Randomize Parameters ---
    prm.L   = 40 + rand() * 60;      
    prm.R   = 2  + rand() * 8;       
    if prm.R > (prm.L/2)*0.5, prm.R = (prm.L/2)*0.4; end 
    
    prm.kx  = 1000 + rand() * 1000;    
    prm.ky  = 100  + rand() * 400;    
    prm.rho = 1500 + rand() * 1000;    % Random Density
    prm.c   = 700  + rand() * 200;     % Random Specific Heat
    
    % Temperatures
    prm.Ti = 600 + rand() * 200; % Initial Hot Temp
    prm.Th = 200 + rand() * 100; % Cold Boundary
    prm.To = 200 + rand() * 100; % Cold Boundary
    
    % --- 2. Run Transient Solver ---
    [T_hist, nodes, ndof, t_vec] = IGA_Transient_Solver(prm, refine_db, dt_db, n_steps_db);
    
    % --- 3. Harvest Random Points (Space-Time) ---
    n_nodes = size(nodes, 1);
    n_times = length(t_vec);
    
    % Randomly select indices for this batch
    rand_node_indices = randi(n_nodes, points_per_sim, 1);
    rand_time_indices = randi(n_times, points_per_sim, 1);
    
    sim_block = zeros(points_per_sim, 13);
    
    for k = 1:points_per_sim
        nid = rand_node_indices(k);
        tid = rand_time_indices(k);
        
        x_val = nodes(nid, 1);
        y_val = nodes(nid, 2);
        time  = t_vec(tid);
        Temp  = T_hist(nid, tid);
        
        % Store [L, R, kx, ky, rho, c, T_initial, T_hole, T_outer, x, y, t, T]
        sim_block(k, :) = [prm.L, prm.R, prm.kx, prm.ky, prm.rho, prm.c, ...
                           prm.Ti, prm.Th, prm.To, x_val, y_val, time, Temp];
    end
    
    Full_Data = [Full_Data; sim_block];
    fprintf('  Harvested %d points. Total Database: %d rows.\n', points_per_sim, size(Full_Data,1));
end

% --- 4. Save to CSV ---
header = {'L', 'R', 'kx', 'ky', 'rho', 'c', 'T_initial', 'T_hole', 'T_outer', 'x', 'y', 't', 'T'};
T_tab = array2table(Full_Data, 'VariableNames', header);
writetable(T_tab, output_file);
fprintf('------------------------------------------------------------\n');
fprintf('Transient Database Saved: %s\n', output_file);
fprintf('------------------------------------------------------------\n');