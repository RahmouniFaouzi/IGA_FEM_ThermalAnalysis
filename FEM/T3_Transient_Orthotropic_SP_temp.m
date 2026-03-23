% % ============================================================
% % T3 TRANSIENT HEAT CONDUCTION IN ORTHOTROPIC PLATE
% % ============================================================
% 
% % MATERIAL PROPERTIES (ORTHOTROPIC PLATE)
% kx = 1003;       % Thermal conductivity in x-direction (W/m.K)
% ky = 171;        % Thermal conductivity in y-direction (W/m.K)
%
% rho = 2000;      % density (kg/m^3)
% c   = 900;       % specific heat (J/kgK)
% 
% ---------------- Geometry ----------------------------------
% a = 0.1;         % plate length (m)
% b = 0.1;         % plate width  (m)
% 
% ---------------- Boundary conditions -----------------------
% Ti   = 273;      % All edges (K)
% T0   = 323;      % at t = 0 
%
% A = [a/2 b/2];
% B = [a/2 5*b/8];
% C = [a/2 3*b/4];
% D = [a/2 7*b/8];
% 
%           ___273k___
%           |        |
%     273k  |        | 273k
%           |__273k__|
%           T0 = 323k at t = 0

clc, clear, close all;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'Functions')))

% MATERIAL PROPERTIES
kx = 1003;       % W/m.K
ky = 171;        % W/m.K
rho = 2000;      % kg/m^3
c   = 900;       % J/kg.K

% GEOMETRY
a = 0.1;         % plate length (m)
b = a;           % plate width  (m)

% MESH PARAMETERS
nx = 24;         % elements in x
ny = nx;         % elements in y

% TIME PARAMETERS
t_target = 3.2;  % target time to evaluate
dt = 0.01;       % small time step for stability

% INITIAL & BOUNDARY CONDITIONS
T0 = 323;        % initial temperature 
Ti = 273;        % fixed edge temperature

% POINTS TO MONITOR
A = [a/2 b/2];
B = [a/2 5*b/8];
C = [a/2 3*b/4];
D = [a/2 7*b/8];

%% ============================================================
% MESH GENERATION
%% ============================================================
[node, elem] = T3_Mesh(a,b,nx,ny);
ndof = size(node,1);

opts = struct('showMesh',1,'showNodes',0,'showABCD',1,'showBoundary',0);
T3_plotMesh(node, elem, a, b, A, B, C, D, opts)

%% ============================================================
% ASSEMBLE GLOBAL MATRICES
%% ============================================================
K = zeros(ndof);
M = zeros(ndof);

for e = 1:size(elem,1)
    n = elem(e,:);
    Ke = T3_Stiffness(node(n,1), node(n,2), diag([kx,ky]));
    Me = T3_Mass(node(n,1), node(n,2), rho, c);
    
    K(n,n) = K(n,n) + Ke;
    M(n,n) = M(n,n) + Me;
end

%% ============================================================
% APPLY BOUNDARY CONDITIONS
%% ============================================================
T     = T0*ones(ndof,1); % initial temperature
cold  = find(abs(node(:,1))<1e-10 | abs(node(:,1)-a)<1e-10 | abs(node(:,2))<1e-10 | abs(node(:,2)-b)<1e-10);
fixed = unique(cold);
free  = setdiff(1:ndof,fixed);
T(fixed) = Ti;

%% ============================================================
% TIME INTEGRATION (BACKWARD EULER)
%% ============================================================
t_current = 0;
while t_current < t_target
    % Implicit system
    Kt = M/dt + K;
    R  = M/dt * T;   % no external heat flux
    
    % Solve for next time step
    T_next = T;
    T_next(free) = Kt(free,free) \ (R(free) - Kt(free,fixed)*T_next(fixed));
    
    % Update
    T = T_next;
    t_current = t_current + dt;
end

%% ============================================================
% INTERPOLATE TEMPERATURES AT POINTS
%% ============================================================
TA = femNodeTemperature(node, A, T, 1e-12);
TB = femNodeTemperature(node, B, T, 1e-12);
TC = femNodeTemperature(node, C, T, 1e-12);
TD = femNodeTemperature(node, D, T, 1e-12);

fprintf('T3 FEM Transient temperatures at t = %.2f s\n', t_target);
fprintf('A: %.2f K\n', TA);
fprintf('B: %.2f K\n', TB);
fprintf('C: %.2f K\n', TC);
fprintf('D: %.2f K\n', TD);

figure(2);
trisurf(elem, node(:,1), node(:,2), T, 'EdgeColor','none');
view(2)
shading interp
colormap(jet)
colorbar
% title('Temperature field at t = 3.20 sec')
% xlabel('x (m)')
% ylabel('y (m)')
set(gca, 'FontSize', 16);

%% Analytic Solution
fprintf('=========================================\n')
% Physical Constants (Bruch & Zyvoloski 1974)
x = A(1); 
y_coords = [A(2), B(2), C(2), D(2)]; 
labels = {'A','B','C','D'};

fprintf('Analytical transient temperature Results at t = 3.20 s:\n');
for i = 1:4
    y = y_coords(i);
    SumVal = 0;
    
    for m = 1:2:200 % High for convergence
        for n = 1:2:200
            alpha_mn = (1/(rho*c)) * (kx*(m*pi/a)^2 + ky*(n*pi/b)^2);
            Bn = 16 / (m * n * pi^2);
            spatial = sin(m*pi*x/a) * sin(n*pi*y/b);
            temp = exp(-alpha_mn * t_target);
            SumVal = SumVal + Bn * spatial * temp;
        end
    end
    T = Ti + (T0 - Ti) * SumVal;
    fprintf('%s: %.2f K\n', labels{i}, T);
end