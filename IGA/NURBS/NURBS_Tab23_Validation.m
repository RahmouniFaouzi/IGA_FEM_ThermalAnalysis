clc, clear, close all
addpath(genpath('../Functions'));

% ====================================================================
%  NURBS based-IGA Temperature Field – Orthotropic Heat Conduction
%  Author : Pr. Khennane Amar    | Email : a.khennane@adfa.edu.au
%           Dr. Rahmouni Faouzi  | Email : rahmounifaouzi01@gmail.com
% ====================================================================

% =========================================
% ORTHOTROPIC HEAT CONDUCTION – IGA (NURBS)
% =========================================
% MATERIAL
kx = 1003;
ky = 171;
% BC
T0 = 273;
Tg = 773;
% Plate Dim
a = 100; 
b = a;

% ==============
% NURBS GEOMETRY
% ==============
CtrlPts = zeros(4,2,2);
CtrlPts(1:3,1,1) = [0; 0; 0];
CtrlPts(1:3,2,1) = [a; 0; 0];
CtrlPts(1:3,1,2) = [0; b; 0];
CtrlPts(1:3,2,2) = [a; b; 0];
CtrlPts(4,:,:)   = 1;

uKnoti = [0 0 a a];
vKnoti = [0 0 b b];

% ==========
% REFINEMENT
% ==========
p = 3; q = 3;
nelX = [8 10 12 14 16 18 20 22 24 26 28 30 40];
nelY = nelX;

storeElements = zeros(numel(nelX),1);
TA = zeros(numel(nelX),1);
TB = TA; TC = TA; TD = TA; dt = TA;

fprintf('--- Analysis starts ---\n');

for ref = 1:numel(nelX)
    t_start = tic;
    Surf = KRefine({uKnoti,vKnoti},CtrlPts,...
        [nelX(ref) nelY(ref)],[p q],[p-1 q-1]);

    % ================
    % EVALUATION DATA
    % ================
    [uKnot,vKnot,noCtrPts,noElems,~,~,~, ...
        elRangeU,~,elRangeV,~,element,index,...
        noGPs,controlPts,weights,p,q] = Evaluation_Nurbs_Par(Surf);

    Surf.ControlPoints = controlPts;
    Surf.Weights       = weights;
    storeElements(ref) = size(element,1);
    
    % Plot Mesh
    if ref == 1
        figure(1)
        hold on
        axis equal
        daspect([1, 1, 1])
        PlotGeo(Surf);
        PlotKnts(Surf);
        PlotCtrlPts(Surf);
        PlotCtrlNet(Surf);
    end 

    % ==================
    % INITIALIZE SYSTEM
    % ==================
    K = sparse(noCtrPts,noCtrPts);
    D = diag([kx ky]);
    [W,Q] = quadrature(noGPs,'GAUSS',2);

    % ========
    % ASSEMBLY
    % ========
    for e = 1:noElems

        sctr = element(e,:);
        pts  = controlPts(sctr,1:2);

        idu  = index(e,1);
        idv  = index(e,2);
        xiE  = elRangeU(idu,:);
        etaE = elRangeV(idv,:);

        Ke = zeros(length(sctr));

        for gp = 1:length(W)

            xi  = parent2ParametricSpace(xiE,Q(gp,1));
            eta = parent2ParametricSpace(etaE,Q(gp,2));

            [dRdxi,dRdeta] = NURBS2Dders([xi;eta],p,q,...
                uKnot,vKnot,weights');

            J1 = [dRdxi; dRdeta] * pts;       
            detJ1 = det(J1);

            dRdx = [dRdxi' dRdeta'] / J1;     % physical derivatives
            B = dRdx';

            J2 = jacobianPaPaMapping(xiE,etaE);

            Ke = Ke + B' * D * B * detJ1 * J2 * W(gp);
        end

        K(sctr,sctr) = K(sctr,sctr) + Ke;
    end

    % =============
    % DIRICHLET BCs
    % =============
    tol = 1e-10;
    x = controlPts(:,1); y = controlPts(:,2);

    cold = find(abs(x)<tol | abs(x-a)<tol | abs(y)<tol);
    hot  = find(abs(y-b)<tol);

    fixed = unique([cold; hot]);
    free  = setdiff(1:noCtrPts,fixed);

    T = zeros(noCtrPts,1);
    T(cold) = T0;
    T(hot)  = Tg;
    T(free) = K(free,free) \ (-K(free,fixed)*T(fixed));

    % ================
    % POINT EVALUATION
    % ================
    Ap = [a/2 b/2];
    Bp = [a/2 5*b/8];
    Cp = [a/2 3*b/4];
    Dp = [a/2 7*b/8];
    
    TA(ref) = eval_IGA_Temperature(Ap,T,Surf);
    TB(ref) = eval_IGA_Temperature(Bp,T,Surf);
    TC(ref) = eval_IGA_Temperature(Cp,T,Surf);
    TD(ref) = eval_IGA_Temperature(Dp,T,Surf);
    
    % 2. DISPLAY PROGRESS
    dt(ref) = toc(t_start);
    fprintf('Mesh %d * %d done in %.4f sec\n', nelX(ref), nelY(ref), dt(ref));
end

fprintf('--- Analysis done in %.4f sec ---\n', sum(dt));

%% DISPLAY RESULTS
% ================
fprintf(' NURBS-IGA Orthotropic Results (K)\n');
fprintf('-------------------------------------------------------------\n');
fprintf('%-10s %8s %8s %8s %8s\n', 'Number of Elements', 'A', 'B', 'C', 'D');
for i = 1:numel(nelX)
    fprintf('      %d          %8.2f %8.2f %8.2f %8.2f\n', ...
        storeElements(i), TA(i), TB(i), TC(i), TD(i));
end
fprintf('-------------------------------------------------------------\n');
fprintf('%-10s         %8.2f %8.2f %8.2f %8.2f\n\n', 'Analytical', 287.17, 309.66, 367.32, 507.70);

% Plot the final result
plotTempField(Surf, T, a, b);

figure(1);
hold on
% Points A, B, C, D
pts = [Ap; Bp; Cp; Dp]; 
labels = {'A','B','C','D'};
axis off
% Plot points
plot(pts(:,1), pts(:,2), 'bo', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
% Add labels
for i = 1:4
    text(pts(i,1), pts(i,2), ['  ' labels{i}], ...
        'FontWeight', 'bold', 'Color', 'k', 'FontSize', 12);
end
hold off;
