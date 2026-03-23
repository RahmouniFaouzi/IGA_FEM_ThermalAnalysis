function plotHolesPlate(patches, global_node_patterns, T)
%PLOTHOLESPLATE Optimized visualization of IGA results
%   Generates two separate figures:
%   1. Temperature Contour Map (Surface)
%   2. Mesh Topology (Knot Lines)

fprintf('Generating IGA plots... \n');

% --- CONFIGURATION ---
surf_res = 50;  % Resolution for temperature surface
line_res = 30;  % Resolution for mesh lines
cmap     = 'jet';

%% --- FIGURE 1: TEMPERATURE FIELD ---
figure('Color', 'w', 'Name', 'IGA Temperature Results', 'NumberTitle', 'off');
hold on; axis equal; axis tight; box on;
title('\textbf{Temperature Field Distribution}', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('X', 'Interpreter', 'latex');
ylabel('Y', 'Interpreter', 'latex');
colormap(cmap);

% Create colorbar
c = colorbar;
c.Label.String = 'Temperature (K)';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.FontSize = 11;

for p = 1:length(patches)
    patch = patches(p);
    
    % 1. Determine parametric bounds for THIS patch
    u_min = min(patch.uKnot); u_max = max(patch.uKnot);
    v_min = min(patch.vKnot); v_max = max(patch.vKnot);
    
    % 2. Generate Parametric Grid
    uu = linspace(u_min, u_max, surf_res);
    vv = linspace(v_min, v_max, surf_res);
    [UU, VV] = meshgrid(uu, vv);
    
    % 3. Initialize buffers
    XX = zeros(size(UU));
    YY = zeros(size(UU));
    TT = zeros(size(UU));
    
    % 4. Evaluate Field
    for k = 1:numel(UU)
        [val, phys] = eval_nurbs_field(UU(k), VV(k), patch, global_node_patterns{p}, T);
        XX(k) = phys(1);
        YY(k) = phys(2);
        TT(k) = val;
    end
    
    % 5. Plot Surface (EdgeColor=none to see pure temperature)
    surf(XX, YY, zeros(size(XX)), TT, 'EdgeColor', 'none', 'FaceColor', 'interp');
end
view(2); % Top-down view

%% --- FIGURE 2: MESH TOPOLOGY ---
figure('Color', 'w', 'Name', 'IGA Mesh Structure', 'NumberTitle', 'off');
hold on; axis equal; axis tight; box on;
title('\textbf{IGA Mesh Topology (Knots)}', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('X', 'Interpreter', 'latex');
ylabel('Y', 'Interpreter', 'latex');

for p = 1:length(patches)
    patch = patches(p);
    
    % CRITICAL FIX: Re-calculate bounds for the current patch (was missing in original)
    u_min = min(patch.uKnot); u_max = max(patch.uKnot);
    v_min = min(patch.vKnot); v_max = max(patch.vKnot);
    
    u_knots = unique(patch.uKnot);
    v_knots = unique(patch.vKnot);
    
    % 1. Draw Vertical Lines (Constant U, Varying V)
    vv_line = linspace(v_min, v_max, line_res);
    for uk = u_knots
        xk = zeros(size(vv_line)); yk = zeros(size(vv_line));
        for k = 1:length(vv_line)
            [~, phys] = eval_nurbs_field(uk, vv_line(k), patch, global_node_patterns{p}, T);
            xk(k) = phys(1); yk(k) = phys(2);
        end
        plot(xk, yk, 'k-', 'LineWidth', 0.8);
    end
    
    % 2. Draw Horizontal Lines (Constant V, Varying U)
    uu_line = linspace(u_min, u_max, line_res);
    for vk = v_knots
        xk = zeros(size(uu_line)); yk = zeros(size(uu_line));
        for k = 1:length(uu_line)
            [~, phys] = eval_nurbs_field(uu_line(k), vk, patch, global_node_patterns{p}, T);
            xk(k) = phys(1); yk(k) = phys(2);
        end
        plot(xk, yk, 'k-', 'LineWidth', 0.8);
    end
end
end