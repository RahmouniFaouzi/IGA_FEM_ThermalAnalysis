clc; clear; close all;
addpath(genpath('../../IGA/Functions'));

fprintf('STARTING FULL GRID DATABASE GENERATION (Valid + NaN Regions)...\n');

% =========================================================================
% CONFIGURATION
% =========================================================================
num_simulations = 100;   
refine_db       = 6;    % High refinement for accuracy
output_file     = 'IGA_DATABASE.csv';

% 1. Grid Density for Material
grid_pts = 15; 
u_vec = linspace(0, 1, grid_pts);
v_vec = linspace(0, 1, grid_pts);

% 2. Density for NaN Regions (Hole/Outside)
nan_density = 0; % Number of random NaN points to add per simulation

Full_Data = [];

for sim = 1:num_simulations
    fprintf('Simulation %d / %d ... ', sim, num_simulations);
    
    % --- Randomize Inputs ---
    L  = 40 + rand() * 60;      
    R  = 2  + rand() * 8;       
    if R > (L/2)*0.6, R = (L/2)*0.5; end 
    
    kx = 700 + rand() * 900;    
    ky = 50  + rand() * 400;    
    Th = 500 + rand() * 500;    
    To = 150 + rand() * 100;    
    
    % --- Solve ---
    [patches, global_patterns, T_sol] = IGA_Solver_Core(L, R, kx, ky, Th, To, refine_db);
    
    sim_block = [];
    
    % =====================================================================
    % PART A: MATERIAL POINTS (Include Edges)
    % =====================================================================
    % Iterate parametric space to get exact mesh points
    for p_idx = 1:4
        for u = u_vec
            for v = v_vec
                
                % Evaluate Exact Physics
                [val, phys] = eval_nurbs_field(u, v, patches(p_idx), global_patterns{p_idx}, T_sol);
                
                % Append Valid Point
                new_row = [L, Th, To, kx, ky, R, phys(1), phys(2), val];
                sim_block = [sim_block; new_row];
            end
        end
    end
    
    % =====================================================================
    % PART B: HOLE POINTS (NaN)
    % =====================================================================
    % Generate points strictly inside radius R (distance < R)
    count = 0;
    while count < nan_density
        % Random point in bounding box of hole
        cx = L/2; cy = L/2;
        r_rand = rand() * (R * 0.95); % 0.95 to stay strictly inside (avoid edge)
        theta  = rand() * 2 * pi;
        
        x_h = cx + r_rand * cos(theta);
        y_h = cy + r_rand * sin(theta);
        
        % Append NaN Point
        new_row = [L, Th, To, kx, ky, R, x_h, y_h, NaN];
        sim_block = [sim_block; new_row];
        count = count + 1;
    end
    
    % =====================================================================
    % PART C: OUTSIDE POINTS (NaN)
    % =====================================================================
    % Generate points slightly outside the [0, L] box
    count = 0;
    while count < nan_density
        % Pick a random side (Left, Right, Bottom, Top)
        side = randi(4);
        margin = L * 0.1; % 10% margin outside
        
        if side == 1 % Left (x < 0)
            x_o = -rand() * margin;
            y_o = rand() * L;
        elseif side == 2 % Right (x > L)
            x_o = L + rand() * margin;
            y_o = rand() * L;
        elseif side == 3 % Bottom (y < 0)
            x_o = rand() * L;
            y_o = -rand() * margin;
        else % Top (y > L)
            x_o = rand() * L;
            y_o = L + rand() * margin;
        end
        
        % Append NaN Point
        new_row = [L, Th, To, kx, ky, R, x_o, y_o, NaN];
        sim_block = [sim_block; new_row];
        count = count + 1;
    end
    
    Full_Data = [Full_Data; sim_block];
    fprintf('Done. (%d points)\n', size(sim_block, 1));
end

% Save
header = {'L', 'Th', 'To', 'kx', 'ky', 'R', 'x', 'y', 'T'};
T_tab = array2table(Full_Data, 'VariableNames', header);
writetable(T_tab, output_file);
fprintf('Database Saved: %s\n', output_file);