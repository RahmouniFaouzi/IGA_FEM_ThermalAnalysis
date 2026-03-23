clc, clear, close all
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'Functions')))

% =====================================
% MATERIAL PROPERTIES (SiC/6061 Al MMC)
% =====================================
vf = 0.5;
km = 171;
kf = 16;

kL = kf * vf + km * (1-vf);
kT = km + ( (vf * (kf-km) * km)/ (0.5*(1-vf) * (kf-km) + km) );

% ========
% LAMINATE
% ========
plyAngles = [0, 30, 30, 0];   % degrees
nply = length(plyAngles);

% ========================
% SPECIFIED TEMPERATURE BC
% ========================
T0 = 273;     % S2, S3, S4
Tg = 773;     % S1 (top)

% ========
% GEOMETRY
% ========
AR = [1 1.5 2 2.5 3];

% =====
% MESH
% =====
nx  = 24; 
ny  = 24;    
tol = 1e-10; 

% =======
% STORAGE
% =======
TA = zeros(size(AR));
TB = zeros(size(AR));
TC = zeros(size(AR));
TD = zeros(size(AR));

% =======================
% LOOP OVER ASPECT RATIOS
% =======================
for r = 1:length(AR)

    b = 0.1;
    a = AR(r)*b;

    [node, elem] = T3_Mesh(a,b,nx,ny);
    ndof = size(node,1);

    K = zeros(ndof);
    F = zeros(ndof,1);

    %% ---------------- ASSEMBLY ----------------
    for e = 1:size(elem,1)

        n = elem(e,:);
        x = node(n,1);
        y = node(n,2);

        Ke = zeros(3);

        % Ply averaging 
        for p = 1:nply
            theta = plyAngles(p)*pi/180;
            D  = condRotate(kL,kT,theta);
            Ke = Ke + T3_Stiffness(x, y, D);
        end

        Ke = Ke / nply;
        K(n,n) = K(n,n) + Ke;
    end

    %% ---------------- DIRICHLET BC ----------------
    T = zeros(ndof,1);

    top    = find(abs(node(:,2)-b)<tol);
    bottom = find(abs(node(:,2))<tol);
    left   = find(abs(node(:,1))<tol);
    right  = find(abs(node(:,1)-a)<tol);

    cold = unique([bottom; left; right]);

    T(top)  = Tg;
    T(cold) = T0;

    fixed = unique([top; cold]);
    free  = setdiff(1:ndof,fixed);

    Ff = F(free) - K(free,fixed)*T(fixed);
    T(free) = K(free,free)\Ff;

    %% ---------------- TABLE POINTS ----------------
    A = [a/2 b/2];
    B = [a/2 5*b/8];
    C = [a/2 3*b/4];
    D = [a/2 7*b/8];

    TA(r) = femNodeTemperature(node,A,T,tol);
    TB(r) = femNodeTemperature(node,B,T,tol);
    TC(r) = femNodeTemperature(node,C,T,tol);
    TD(r) = femNodeTemperature(node,D,T,tol);

    fprintf('AR = %.1f | A = %.2f  B = %.2f  C = %.2f  D = %.2f\n', ...
        AR(r), TA(r), TB(r), TC(r), TD(r));
end

% PLOT
figure; hold on; box on
plot(AR,TA,'ks-','LineWidth',1.8)
plot(AR,TB,'ro-','LineWidth',1.8)
plot(AR,TC,'b^-','LineWidth',1.8)
plot(AR,TD,'gv-','LineWidth',1.8)
legend('Point A (k)','Point B (k)','Point C (k)','Point D (k)','Location','northwest')
set(gca, 'FontSize', 16);
xlabel('Aspect Ratios')
ylabel('Nodal Temperatures (K)')
grid off