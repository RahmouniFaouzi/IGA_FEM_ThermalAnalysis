clc, clear, close all
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'Functions')))

%% ============================================================
%  T3 FEM Laminated plate Mixed BC
%  [0/30]s [0/45]s [0/60]s [0/90]s
%% ============================================================

%% -------- Fiber / Matrix properties -------------------------
kf = 16;      % Fiber conductivity (W/m.K)
km = 171;     % Matrix conductivity (W/m.K)
vf = 0.5;

% On-axis conductivities
k1 = kf*vf + km*(1-vf);
k2 = km + (vf*(kf-km)*km)/(0.5*(1-vf)*(kf-km)+km);

%% -------- Ply orientations ---------------------------------
angles = [30 45 60 90];      % degrees
nOri   = length(angles);

%% -------- Boundary conditions -------------------------------
Tg   = 453;
Tinf = 298;
h    = 50;

%% -------- Geometry ------------------------------------------
AR = [1 1.5 2 2.5 3];

%% -------- Mesh ----------------------------------------------
nx = 24; ny = 24;
tol = 1e-12;

%% -------- Storage -------------------------------------------
TA = zeros(nOri,length(AR));
TB = zeros(nOri,length(AR));
TC = zeros(nOri,length(AR));
TD = zeros(nOri,length(AR));

%% ============================================================
%  LOOP OVER ORIENTATIONS
%% ============================================================
for o = 1:nOri

    theta = angles(o)*pi/180;

    % Effective laminate conductivities ([0/theta]s)
    kx = 0.5*( k1 + (k1*cos(theta)^2 + k2*sin(theta)^2) );
    ky = 0.5*( k2 + (k1*sin(theta)^2 + k2*cos(theta)^2) );

    %% ========================================================
    %  LOOP OVER ASPECT RATIOS
    %% ========================================================
    for r = 1:length(AR)

        b = 0.1;
        a = AR(r)*b;

        [node, elem] = T3_Mesh(a,b,nx,ny);
        ndof = size(node,1);

        K = zeros(ndof);
        F = zeros(ndof,1);

        %% -------- Assembly ----------------------------------
        for e = 1:size(elem,1)

            n = elem(e,:);
            Ke = T3_Stiffness(node(n,1),node(n,2), diag([kx,ky]));
            K(n,n) = K(n,n) + Ke;

            % Convection on right edge (x = a)
            xe = node(n,1);
            onRight = abs(xe-a) < tol;

            if sum(onRight)==2
                id = find(onRight);
                j = id(1); k2 = id(2);
                i = setdiff(1:3,id);

                L = norm(node(n(j),:) - node(n(k2),:));

                Hloc = (h*L/6)*[0 0 0;
                                0 2 1;
                                0 1 2];

                Rloc = (h*Tinf*L/2)*[0;1;1];

                perm = [i j k2];
                K(n,n) = K(n,n) + Hloc(perm,perm);
                F(n)   = F(n)   + Rloc(perm);
            end
        end

        %% -------- Essential BC (top edge) -------------------
        T = zeros(ndof,1);
        hot = find(abs(node(:,2)-b)<tol);
        T(hot) = Tg;

        free = setdiff(1:ndof,hot);
        T(free) = K(free,free)\(F(free) - K(free,hot)*T(hot));

        %% -------- Locations A B C D -------------------------
        A = [a/2 b/2];
        B = [a/2 5*b/8];
        C = [a/2 3*b/4];
        D = [a/2 7*b/8];

        TA(o,r) = femNodeTemperature(node,A,T,tol);
        TB(o,r) = femNodeTemperature(node,B,T,tol);
        TC(o,r) = femNodeTemperature(node,C,T,tol);
        TD(o,r) = femNodeTemperature(node,D,T,tol);
    end
end

%% ============================================================
%  TABLE OUTPUT
%% ============================================================
for o = 1:nOri
    fprintf('\n[0/%d]s Laminate Mixed BC\n',angles(o));
    fprintf('------------------------------------------\n');
    fprintf('AR     A        B        C        D\n');
    fprintf('------------------------------------------\n');
    for r=1:length(AR)
        fprintf('%-4.2f %8.2f %8.2f %8.2f %8.2f\n', ...
            AR(r),TA(o,r),TB(o,r),TC(o,r),TD(o,r));
    end
end

%% ============================================================
%  FIGURES (paper-style)
%% ============================================================
for o = 1:nOri
    figure; hold on; grid on; box on
    plot(AR,TA(o,:),'ks-','LineWidth',1.5)
    plot(AR,TB(o,:),'ro-','LineWidth',1.5)
    plot(AR,TC(o,:),'b^-','LineWidth',1.5)
    plot(AR,TD(o,:),'gv-','LineWidth',1.5)
    xlabel('Aspect Ratio')
    grid off
    ylabel('Nodal Temperature (K)')
    title(sprintf('[0/%d]s laminate Mixed BC (T3 FEM)',angles(o)))
    legend('Location A','Location B','Location C','Location D','Location','northwest')
end