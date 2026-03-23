clc, clear, close all
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'Functions')))

% -------------------------------------------------------
%% TABLE 1-2: INCORRECT RESULTS FROM REFERENCE 1-2025
%  The correct results are from the original 1994 reference.
%  Conductivity kx was changed between the references.
%---------------------------------------------------------

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
b = 0.1;      % Plate height mm
a = b;        % Plate width mm

%% ============================================================
% MESH PARAMETERS
%% ============================================================
NX = [8, 10, 12, 14, 16, 18, 20, 22, 24, 40];
NY = NX;
tol = 1e-10;

%% ============================================================
% STORAGE FOR TABLE 2-3
%% ============================================================
storeElments = zeros(1,length(NX));
TA = zeros(1,length(NX));
TB = zeros(1,length(NX));
TC = zeros(1,length(NX));
TD = zeros(1,length(NX));

%% ============================================================
% LOOP OVER ASPECT RATIOS
%% ============================================================
for r = 1:length(NX)

    nx = NX(r);
    ny = NY(r);

    %% Mesh generation
    [node, elem] = T3_Mesh(a,b,nx,ny);
    ndof = size(node,1);
    storeElments(r) = size(elem,1);

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

    %% Locations along centerline
    A = [a/2  b/2];
    B = [a/2  5*b/8];
    C = [a/2  3*b/4];
    D = [a/2  7*b/8];
    
    TA(r) = femNodeTemperature(node, A, T, tol);
    TB(r) = femNodeTemperature(node, B, T, tol);
    TC(r) = femNodeTemperature(node, C, T, tol);
    TD(r) = femNodeTemperature(node, D, T, tol);
end

fprintf('\nTable 5. Nodal temperatures (K)\n');
fprintf('-------------------------------------------\n');

% Header
fprintf('%-10s %8s %8s %8s %8s\n', 'Number of Elements', 'A', 'B', 'C', 'D');

% Data rows
for i = 1:length(NX)
    fprintf('      %d             %8.2f %8.2f %8.2f %8.2f\n', ...
        storeElments(i), ...
        TA(i), TB(i), TC(i), TD(i));
end

figure;
trisurf(elem, node(:,1), node(:,2), T, 'EdgeColor','none');
view(2)
shading interp
colormap(jet)
colorbar
title('Temperature field')
xlabel('x')
ylabel('y')

% PLOT 2: MESH TOPOLOGY 
figure('Name', 'Finite Element Mesh', 'Color', 'w');
triplot(elem, node(:,1), node(:,2), 'k');
axis equal;
title(['Finite Element Mesh (Triangular Elements: ', num2str(size(elem,1)), ')']);
xlabel('X');
ylabel('Y');
axis off;
hold on;

% Points A, B, C, D
pts = [A; B; C; D]; 
labels = {'A','B','C','D'};
% Plot points
plot(pts(:,1), pts(:,2), 'bo', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
% Add labels
for i = 1:4
    text(pts(i,1), pts(i,2), ['  ' labels{i}], ...
        'FontWeight', 'bold', 'Color', 'k', 'FontSize', 12);
end
hold off;

%% Analytic Solution Bruch & Zyvoloski (1974)
fprintf('---------------------------------------------\n')
fprintf('Steady-State Analytical Results (Bruch & Zyvoloski Equation):\n');

% Orthotropic ratio
eps = sqrt(kx/ky);

for i = 1:4
    x = pts(i,1);
    y = pts(i,2);
    SumVal = 0;
    
    % n must be odd for the Fourier coefficients to be non-zero
    for n = 1:2:200 
        % beta_n as defined in the paper for orthotropic media
        beta_n = (n * pi / a) * eps;
        
        % Standard Fourier coefficient for a square domain
        An = 4 / (n * pi);
        
        % Using the exponential identity to ensure numerical stability
        ratio = (exp(beta_n*(y-b)) - exp(-beta_n*(y+b))) / (1 - exp(-2*beta_n*b));
        
        term = An * sin(n * pi * x / a) * ratio;
        SumVal = SumVal + term;
    end
    
    T_ss = T0 + (Tg - T0) * SumVal;
    fprintf('Point %s [%.3f, %.3f]: %.4f K\n', labels{i}, x, y, T_ss);
end