clc, clear, close all
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'Functions')))

%% ============================================================
% MATERIAL PROPERTIES
%% ============================================================
kx  = 1003;
ky  = 171;

%% ============================================================
% MIXED BC PARAMETERS
%% ============================================================
Tg   = 453;        % Top edge temperature (K)
Tinf = 298;        % Ambient temperature (K)
h    = 50;         % Convection coefficient

%% ============================================================
% GEOMETRY
%% ============================================================
AR = [1 1.5 2 2.5 3];

%% ============================================================
% MESH
%% ============================================================
nx = 24; ny = nx;
tol = 1e-12;
opts = struct('showMesh',1,'showNodes',0,'showABCD',1,'showBoundary',0);
%% ============================================================
% STORAGE
%% ============================================================
TA = zeros(1,length(AR));
TB = zeros(1,length(AR));
TC = zeros(1,length(AR));
TD = zeros(1,length(AR));

%% ============================================================
% LOOP OVER ASPECT RATIOS
%% ============================================================
for r = 1:length(AR)

    b = 0.1;
    a = AR(r)*b;

    [node, elem] = T3_Mesh(a,b,nx,ny);
    ndof = size(node,1);
    
    A = [a/2 b/2];
    B = [a/2 5*b/8];
    C = [a/2 3*b/4];
    D = [a/2 7*b/8];

    if (AR(r) == 1)
        T3_plotMesh(node, elem, a, b, A, B, C, D, opts)
    end

    %% ========================================================
    % ASSEMBLE K = KT + HT
    %% ========================================================
    K = zeros(ndof);
    F = zeros(ndof,1);

    for e = 1:size(elem,1)

        n = elem(e,:);
        Ke = T3_Stiffness(node(n,1),node(n,2), diag([kx,ky]));
        K(n,n) = K(n,n) + Ke;

        % Convection right edge
        xe = node(n,1);
        onRight = abs(xe-a) < tol;

        if sum(onRight)==2
            id = find(onRight);
            j = id(1); k2 = id(2);
            i = setdiff(1:3,id);

            L = norm(node(n(j),:) - node(n(k2),:));

            % HT term
            Hloc = (h*L/6)*[0 0 0; 0 2 1; 0 1 2];
            perm = [i j k2];
            Hloc = Hloc(perm,perm);
            K(n,n) = K(n,n) + Hloc;

            % r_infinity
            Rloc = (h*Tinf*L/2)*[0;1;1];
            Rloc = Rloc(perm);
            F(n) = F(n) + Rloc;
        end
    end

    %% ========================================================
    % BC top edge
    %% ========================================================
    T = zeros(ndof,1);
    hot = find(abs(node(:,2)-b)<tol);
    T(hot) = Tg;

    free = setdiff(1:ndof,hot);
    Ff = F(free) - K(free,hot)*T(hot);
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
% OUTPUT
%% ============================================================
fprintf('\nTable 7 - Isotropic - Mixed BC\n');
fprintf('----------------------------------------\n');
fprintf('A '); fprintf('%10.2f',TA); fprintf('\n');
fprintf('B '); fprintf('%10.2f',TB); fprintf('\n');
fprintf('C '); fprintf('%10.2f',TC); fprintf('\n');
fprintf('D '); fprintf('%10.2f',TD); fprintf('\n');

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

    colormap('jet');      
    colorbar;
    xlabel('X-Coordinate');
    ylabel('Y-Coordinate');
    title(sprintf('Temperature field (AR = %.2f)', AR(ARi)));
    view(2);            
    axis equal tight;    
    hold off
end