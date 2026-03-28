clc, clear, close all
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'Functions')))

%% ============================================================
% MATERIAL PROPERTIES
%% ============================================================
kx   = 1003;     % W/mK
ky   = 171;      % W/mK
Tg   = 453;      % Top temperature
Tinf = 298;      % Ambient temperature
h    = 50;       % Convection coefficient

%% ============================================================
% GEOMETRY
%% ============================================================
AR = [1 1.5 2 2.5 3];
b  = 0.1;

%% ============================================================
% MESH
%% ============================================================
nx = 16; ny = 16;
tol = 1e-10;

%% ============================================================
% STORAGE
%% ============================================================
TA = zeros(1,length(AR));
TB = zeros(1,length(AR));
TC = zeros(1,length(AR));
TD = zeros(1,length(AR));

%% ============================================================
% MAIN LOOP
%% ============================================================
for r = 1:length(AR)

    a = AR(r)*b;
    [node,elem] = Q4_Mesh(a,b,nx,ny);
    ndof = size(node,1);

    K = zeros(ndof);
    F = zeros(ndof,1);

    %% ========== ASSEMBLY ==========
    for e = 1:size(elem,1)
        n = elem(e,:);
        xe = node(n,1);
        ye = node(n,2);

        % Element stiffness
        Ke = Q4_Stiffness(xe,ye,kx,ky);
        K(n,n) = K(n,n) + Ke;

        % --------- CONVECTION ON RIGHT EDGE ----------
        idx = find(abs(xe - a) < tol);  % only nodes on x=a
        if length(idx) == 2
            i1 = idx(1); i2 = idx(2);
            L  = abs(ye(i2) - ye(i1));

            Hloc = zeros(4);
            Hloc([i1 i2],[i1 i2]) = (h*L/6)*[2 1; 1 2];

            Rloc = zeros(4,1);
            Rloc([i1 i2]) = (h*Tinf*L/2)*[1;1];

            K(n,n) = K(n,n) + Hloc;
            F(n)   = F(n) + Rloc;
        end
    end

    %% ========== APPLY DIRICHLET BC ==========
    T = zeros(ndof,1);
    top = find(abs(node(:,2)-b) < tol);
    T(top) = Tg;

    free = setdiff(1:ndof,top);
    Ff = F(free) - K(free,top)*T(top);
    T(free) = K(free,free)\Ff;

    %% ========== SAMPLE POINTS ==========
    A = [a/2 b/2];
    B = [a/2 5*b/8];
    C = [a/2 3*b/4];
    D = [a/2 7*b/8];

    TA(r) = femNodeTemperature(node,A,T,tol);
    TB(r) = femNodeTemperature(node,B,T,tol);
    TC(r) = femNodeTemperature(node,C,T,tol);
    TD(r) = femNodeTemperature(node,D,T,tol);
end

%% ============================================================
% OUTPUT
%% ============================================================
fprintf('\nTable 7 Q4 Isotropic MIXED BC\n');
fprintf('A '); fprintf('%10.2f',TA); fprintf('\n');
fprintf('B '); fprintf('%10.2f',TB); fprintf('\n');
fprintf('C '); fprintf('%10.2f',TC); fprintf('\n');
fprintf('D '); fprintf('%10.2f',TD); fprintf('\n');
