clc, clear, close all

%% -----------------------------
% Analytical References
%% -----------------------------
ref_A = 287.17; 
ref_B = 309.66; 
ref_C = 367.32; 
ref_D = 507.70;

%% -----------------------------
% FEM (T3 elements)
%% -----------------------------
fem_elem = [128, 200, 288, 392, 512, 648, 800, 968, 1152, 3200];
fem_A = [289.27, 288.53, 288.12, 287.87, 287.71, 287.60, 287.52, 287.46, 287.41, 287.26];
fem_B = [313.39, NaN, NaN, NaN, 310.65, NaN, NaN, NaN, 310.11, 309.82];
fem_C = [372.19, NaN, 369.59, NaN, 368.62, NaN, 368.16, NaN, 367.90, 367.53];
fem_D = [508.31, NaN, NaN, NaN, 507.57, NaN, NaN, NaN, 507.60, 507.66];

%% -----------------------------
% NURBS IGA
%% -----------------------------
nurbs_elem = [64, 100, 144, 196, 256, 324, 400, 484, 576, 676, 784, 900, 1600];
n_A = [287.81, 287.60, 287.46, 287.38, 287.33, 287.30, 287.27, 287.26, 287.24, 287.23, 287.22, 287.22, 287.20];
n_B = [311.43, 310.72, 310.39, 310.21, 310.08, 309.99, 309.93, 309.88, 309.85, 309.82, 309.80, 309.78, 309.73];
n_C = [371.43, 369.84, 369.15, 368.64, 368.33, 368.11, 367.96, 367.85, 367.77, 367.70, 367.65, 367.60, 367.48];
n_D = [514.36, 512.13, 510.92, 510.04, 509.45, 509.09, 508.84, 508.64, 508.49, 508.37, 508.28, 508.21, 507.98];

%% -----------------------------
% HB-IGA
%% -----------------------------
hb_elem = [94, 136, 186, 244, 376, 456, 544, 640, 744, 964, 1090, 1224, 2176];
hb_A = [286.95, 287.06, 287.09, 287.11, 287.13, 287.14, 287.14, 287.15, 287.15, 287.15, 287.16, 287.16, 287.16];
hb_B = [309.35, 309.36, 309.46, 309.52, 309.55, 309.57, 309.59, 309.60, 309.61, 309.62, 309.63, 309.63, 309.64];
hb_C = [365.84, 366.57, 366.81, 366.94, 367.03, 367.09, 367.13, 367.17, 367.19, 367.21, 367.22, 367.24, 367.27];
hb_D = [505.90, 506.51, 506.77, 506.98, 507.17, 507.28, 507.37, 507.42, 507.47, 507.50, 507.53, 507.55, 507.62];

%% -----------------------------
% data
%% -----------------------------
points_label = {'A', 'B', 'C', 'D'};
points_ref = [ref_A, ref_B, ref_C, ref_D];
points_fem = {fem_A, fem_B, fem_C, fem_D};
points_nurbs = {n_A, n_B, n_C, n_D};
points_hb = {hb_A, hb_B, hb_C, hb_D};

%% -----------------------------
% Plot Loop
%% -----------------------------
for i = 1:4
    
    figure('Color','w','Units','centimeters','Position',[5 5 13 9]);
    hold on;
    
    % colors
    c_ref   = [0 0 0];
    c_fem   = [0 0.4470 0.7410];
    c_nurbs = [0.4660 0.6740 0.1880];
    c_hb    = [0.8500 0.3250 0.0980];
    
    xmin = min([fem_elem nurbs_elem hb_elem]);
    xmax = max([fem_elem nurbs_elem hb_elem]);
    
    % Analytical reference line
    plot([xmin xmax], ...
         [points_ref(i) points_ref(i)], ...
         '--','Color',c_ref,'LineWidth',1.3);
    
    % FEM
    current_fem = points_fem{i};
    idx = ~isnan(current_fem);
    plot(fem_elem(idx), current_fem(idx), ...
        '-s','Color',c_fem,'LineWidth',1.3,...
        'MarkerSize',5,'MarkerFaceColor','w');
    
    % NURBS
    plot(nurbs_elem, points_nurbs{i}, ...
        '-<','Color',c_nurbs,'LineWidth',1.3,...
        'MarkerSize',5,'MarkerFaceColor','w');
    
    % HB
    plot(hb_elem, points_hb{i}, ...
        '-o','Color',c_hb,'LineWidth',1.3,...
        'MarkerSize',5,'MarkerFaceColor','w');
    
    % Axes styling
    set(gca,...
        'FontName','Times New Roman',...
        'FontSize',12,...
        'LineWidth',1.3,...
        'TickDir','in');
    
    xlabel('Number of Elements','FontSize',13,'FontWeight','bold');
    ylabel(['Temperature at Point ',points_label{i},' (K)'],...
        'FontSize',13,'FontWeight','bold');
    
    legend({'Analytic','T3 FEM','NURBS IGA','HB-Spline'},...
        'Location','best',...
        'FontSize',14,...
        'Box','on');
    
    % y-limits
    all_data = [points_fem{i}, points_nurbs{i}, points_hb{i}];
    all_data = all_data(~isnan(all_data));
    ylim([min(all_data)-0.4, max(all_data)+0.4]);
   
    grid off
    hold off;
    box on
    
end
