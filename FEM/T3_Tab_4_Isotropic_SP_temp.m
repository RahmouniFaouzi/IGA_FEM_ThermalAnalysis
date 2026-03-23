 clear, clc, close all

addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'Functions')))
%% ============================================================
% MATERIAL PROPERTIES (ISOTROPIC PLATE)
%% ============================================================
k  = 150;        % Thermal conductivity (W/m.K)
kx = k;          % Conductivity in x-direction
ky = k;          % Conductivity in y-direction

%% ============================================================
% BOUNDARY TEMPERATURES (DIRICHLET BCs)
%% ============================================================
T0 = 273;        % Cold temperature (K) on left, right, bottom edges
Tg = 773;        % Hot temperature (K) on top edge

%% ============================================================
% GEOMETRY: ASPECT RATIOS (a/b)
%% ============================================================
AR = [1 1.5 2 2.5 3];

%% ============================================================
% MESH PARAMETERS
%% ============================================================
nx = 24;         % Number of elements in x-direction
ny = 24;         % Number of elements in y-direction

tol = 1e-10;     % Numerical tolerance for boundary detection
opts = struct('showMesh',1,'showNodes',0,'showABCD',1,'showBoundary',0);

%% ============================================================
% STORAGE FOR TABLE 4 RESULTS
% Each row corresponds to a location (A,B,C,D)
% Each column corresponds to an aspect ratio
%% ============================================================
TA = zeros(1,length(AR));   % Temperatures at point A
TB = zeros(1,length(AR));   % Temperatures at point B
TC = zeros(1,length(AR));   % Temperatures at point C
TD = zeros(1,length(AR));   % Temperatures at point D

%% ============================================================
% LOOP OVER ALL ASPECT RATIOS
%% ============================================================
for r = 1:length(AR)

    b = 0.1;            % Plate height
    a = AR(r)*b;        % Plate width for current aspect ratio

    %% ========================================================
    % MESH GENERATION (T3 TRIANGULAR ELEMENTS)
    %% ========================================================
    [node, elem] = T3_Mesh(a,b,nx,ny);
    % node : [x y] coordinates of each node
    % elem : connectivity matrix (each row = 1 triangular element)
    
    A = [a/2  b/2];
    B = [a/2  5*b/8];
    C = [a/2  3*b/4];
    D = [a/2  7*b/8];

    % Plot Mesh
    if (AR(r) == 1)
        T3_plotMesh(node, elem, a, b, A, B, C, D, opts)
    end

    ndof = size(node,1);   % Number of degrees of freedom (1 DOF per node)

    %% ========================================================
    % GLOBAL STIFFNESS MATRIX ASSEMBLY
    %% ========================================================
    K = zeros(ndof);       % Global conductivity (stiffness) matrix
    F = zeros(ndof,1);     % Global heat source vector (zero here)

    for e = 1:size(elem,1)
        n = elem(e,:);
        % n contains the 3 node numbers of element e
        Ke = T3_Stiffness(node(n,1),node(n,2),diag([kx,ky]));
        
        % Compute element conductivity matrix using FEM formulation
        K(n,n) = K(n,n) + Ke;
    end

    %% ========================================================
    % APPLY DIRICHLET BOUNDARY CONDITIONS
    %% ========================================================
    T = zeros(ndof,1);     % Initialize temperature vector

    % Identify cold boundary nodes:
    % x = 0   (left edge)
    % x = a   (right edge)
    % y = 0   (bottom edge)
    cold = find( abs(node(:,1))<tol | ...
        abs(node(:,1)-a)<tol | ...
        abs(node(:,2))<tol );

    % Identify hot boundary nodes:
    % y = b   (top edge)
    hot  = find( abs(node(:,2)-b)<tol );

    fixed = unique([cold; hot]);     % Nodes with known temperature
    free  = setdiff(1:ndof,fixed);   % Unknown temperature nodes

    % Assign prescribed temperatures
    T(cold) = T0;
    T(hot)  = Tg;

    %% ========================================================
    % SOLVE FEM SYSTEM
    %% ========================================================
    % Reduced system:
    % K_ff * T_f = - K_fc * T_c
    Ff = F(free) - K(free,fixed)*T(fixed);

    % Solve for unknown nodal temperatures
    T(free) = K(free,free)\Ff;

    % Interpolate temperature at these locations
    TA(r) = femNodeTemperature(node, A, T, tol);
    TB(r) = femNodeTemperature(node, B, T, tol);
    TC(r) = femNodeTemperature(node, C, T, tol);
    TD(r) = femNodeTemperature(node, D, T, tol);
    STE(1:ndof, r) = T;   % store nodal temperatures
    Snode{r} = node;      % store node coordinates
end

%% ============================================================
% PRINT TABLE 4 (SAME FORMAT AS PAPER)
%% ============================================================
fprintf('\nTable 4. Nodal temperatures (K) Isotropic plate\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Location ');
for r = 1:length(AR)
    fprintf('   AR=%-4.1f',AR(r));
end
fprintf('\n-------------------------------------------------------------\n');

fprintf('A        '); fprintf('%10.2f',TA); fprintf('\n');
fprintf('B        '); fprintf('%10.2f',TB); fprintf('\n');
fprintf('C        '); fprintf('%10.2f',TC); fprintf('\n');
fprintf('D        '); fprintf('%10.2f',TD); fprintf('\n');

%% ========================================================
% SMOOTH 2D TEMPERATURE CONTOUR
%% ========================================================

% Create a smooth triangular surface plot
for ARi = 1:length(AR)
    figure(ARi+1); hold on;
    node = Snode{ARi};
    trisurf(elem, node(:,1), node(:,2), STE(:,ARi), ...
        'FaceColor', 'interp', ...   
        'EdgeColor', 'none');     

    colormap('jet');      
    colorbar;
    xlabel('X-Coordinate');
    ylabel('Y-Coordinate');
    title(sprintf('Temperature field (AR = %.2f)', AR(ARi)));
    view(2);             
    axis equal tight;   
    set(gca, 'FontSize', 16);
    hold off
end