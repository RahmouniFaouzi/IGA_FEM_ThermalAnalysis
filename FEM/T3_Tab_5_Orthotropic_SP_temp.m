clc, clear, close all
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'Functions')))

%% ============================================================
% MATERIAL PROPERTIES (ORTHOTROPIC PLATE)
%% ============================================================
kx = 1003;       % Thermal conductivity in x-direction (W/m.K)
ky = 171;        % Thermal conductivity in y-direction (W/m.K)

%% ============================================================
% BOUNDARY TEMPERATURES (SPECIFIED TEMPERATURE BC)
%% ============================================================
T0 = 273;        % Cold edges temperature (K)
Tg = 773;        % Hot edge temperature (K)

%% ============================================================
% ASPECT RATIOS (a/b)
%% ============================================================
AR = [1 1.5 2 2.5 3];

%% ============================================================
% MESH PARAMETERS
%% ============================================================
nx = 24;
ny = 24;
tol = 1e-8;
opts = struct('showMesh',1,'showNodes',0,'showABCD',1,'showBoundary',0);
%% ============================================================
% STORAGE FOR TABLE 5
%% ============================================================
TA = zeros(1,length(AR));
TB = zeros(1,length(AR));
TC = zeros(1,length(AR));
TD = zeros(1,length(AR));

%% ============================================================
% LOOP OVER ASPECT RATIOS
%% ============================================================
for r = 1:length(AR)

    b = 0.1;            % Plate height
    a = AR(r)*b;        % Plate width

    A = [a/2  b/2];
    B = [a/2  5*b/8];
    C = [a/2  3*b/4];
    D = [a/2  7*b/8];
    
    %% Mesh generation
    [node, elem] = T3_Mesh(a,b,nx,ny);
    ndof = size(node,1);

    % Plot Mesh
    if (AR(r) == 1)
        T3_plotMesh(node, elem, a, b, A, B, C, D, opts)
    end
    
    %% Global stiffness matrix
    K = zeros(ndof);
    F = zeros(ndof,1);

    for e = 1:size(elem,1)
        n = elem(e,:);
        Ke = T3_Stiffness(node(n,1),node(n,2),diag([kx,ky]));
        K(n,n) = K(n,n) + Ke;
    end

    %% Boundary conditions
    T = zeros(ndof,1);

    cold = find( abs(node(:,1))<tol | ...
                 abs(node(:,1)-a)<tol | ...
                 abs(node(:,2))<tol );

    hot  = find( abs(node(:,2)-b)<tol );

    fixed = unique([cold; hot]);
    free  = setdiff(1:ndof,fixed);

    T(cold) = T0;
    T(hot)  = Tg;

    %% Solve system
    Ff = F(free) - K(free,fixed)*T(fixed);
    T(free) = K(free,free)\Ff;

    %% Locations
    TA(r) = femNodeTemperature(node, A, T, tol);
    TB(r) = femNodeTemperature(node, B, T, tol);
    TC(r) = femNodeTemperature(node, C, T, tol);
    TD(r) = femNodeTemperature(node, D, T, tol);
    STE(1:ndof, r) = T;   % store nodal temperatures
    Snode{r} = node;      % store node coordinates
end

%% ============================================================
% PRINT TABLE 5 (PAPER FORMAT)
%% ============================================================
fprintf('\nTable 5. Nodal temperatures (K) Orthotropic plate\n');
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
        'FaceColor', 'interp', ...   % smooth color interpolation
        'EdgeColor', 'none');         % hide mesh edges

    colormap('jet');      % jet colormap (rainbow)
    colorbar;
    xlabel('X-Coordinate');
    ylabel('Y-Coordinate');
    title(sprintf('Temperature field (AR = %.2f)', AR(ARi)));
    view(2);             % top-down view for 2D contour effect
    axis equal tight;    % keep aspect ratio
    set(gca, 'FontSize', 16);
    hold off
end