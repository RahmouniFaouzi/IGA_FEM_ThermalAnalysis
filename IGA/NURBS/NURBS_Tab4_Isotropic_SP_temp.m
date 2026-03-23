clc, clear, close all
addpath(genpath('../Functions')); 

%% ============================================================
% PARAMETERS (ISOTROPIC PLATE IN IGA)
%% ============================================================
k  = 150;        % Thermal conductivity (W/m.K)

% BOUNDARY TEMPERATURES
T0 = 273;        % Cold (Left, Right, Bottom)
Tg = 773;        % Hot (Top)

% GEOMETRY: ASPECT RATIOS (a/b)
AR = [1 1.5 2 2.5 3];

% MESH PARAMETERS (NURBS)
nx = 24;         % Number of elements in x
ny = 24;         % Number of elements in y
p  = 3;          % Degree in u
q  = 3;          % Degree in v

% STORAGE FOR RESULTS
TA = zeros(1,length(AR)); 
TB = zeros(1,length(AR)); 
TC = zeros(1,length(AR)); 
TD = zeros(1,length(AR));

%% ============================================================
% LOOP OVER ASPECT RATIOS
%% ============================================================
fprintf('Starting IGA Analysis...\n');

for r = 1:length(AR)
    
    b = 10;          % Plate height
    a = AR(r)*b;     % Plate width
    
    %% 1. GEOMETRY GENERATION (Linear Basis)
    % ======================================
    CtrlPts = zeros(4,2,2);
    CtrlPts(1:3,1,1) = [0; 0; 0];
    CtrlPts(1:3,2,1) = [a; 0; 0];
    CtrlPts(1:3,1,2) = [0; b; 0];
    CtrlPts(1:3,2,2) = [a; b; 0];
    CtrlPts(4,:,:)   = 1; % Weights

    uKnoti = [0 0 a a];
    vKnoti = [0 0 b b];
    
    %% 2. K-REFINEMENT 
    % ================
    Surf = KRefine({uKnoti,vKnoti}, CtrlPts, [nx ny], [p q], [1 1]);
    
    % Evaluate Basis Functions
    [uKnot,vKnot,noCtrPts,noElems,~,~,~, ...
     elRangeU,~,elRangeV,~,element,index,...
     noGPs,controlPts,weights,p,q] = Evaluation_Nurbs_Par(Surf);

    Surf.ControlPoints = controlPts;
    Surf.Weights       = weights;

    %% 3. ASSEMBLY
    % ============
    K_global = sparse(noCtrPts, noCtrPts);
    
    % Gaussian Quadrature
    [W, Q] = quadrature(noGPs, 'GAUSS', 2); 

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
            
            [dRdxi, dRdeta] = NURBS2Dders([xi;eta], p, q, uKnot, vKnot, weights');
            
            J1    = [dRdxi; dRdeta] * pts;        
            detJ1 = det(J1);
            
            % Inverse Jacobian for physical derivatives
            invJ1 = inv(J1);
            dRdx  = invJ1 * [dRdxi; dRdeta];
            B     = dRdx;
            
            J2    = jacobianPaPaMapping(xiE,etaE);
            
            Ke = Ke + B' * k * B * detJ1 * J2 * W(gp);
        end
        K_global(sctr,sctr) = K_global(sctr,sctr) + Ke;
    end
    
    %% 4. BOUNDARY CONDITIONS (Dirichlet)
    % ===================================
    tol = 1e-8;
    x_coords = controlPts(:,1); 
    y_coords = controlPts(:,2);

    % Cold Edges (Left, Right, Bottom)
    cold = find(abs(x_coords) < tol | abs(x_coords - a) < tol | abs(y_coords) < tol);
    
    % Hot Edge (Top)
    hot  = find(abs(y_coords - b) < tol);

    fixed = unique([cold; hot]);
    free  = setdiff(1:noCtrPts, fixed);

    T_sol = zeros(noCtrPts, 1);
    T_sol(cold) = T0;
    T_sol(hot)  = Tg;
    
    % Solve for free DOFs
    T_sol(free) = K_global(free,free) \ (-K_global(free,fixed) * T_sol(fixed));
    
    %% 5. POINT EVALUATION
    % ====================
    Pt_A = [a/2  b/2];
    Pt_B = [a/2  5*b/8];
    Pt_C = [a/2  3*b/4];
    Pt_D = [a/2  7*b/8];
    
    TA(r) = eval_IGA_Temperature(Pt_A, T_sol, Surf);
    TB(r) = eval_IGA_Temperature(Pt_B, T_sol, Surf);
    TC(r) = eval_IGA_Temperature(Pt_C, T_sol, Surf);
    TD(r) = eval_IGA_Temperature(Pt_D, T_sol, Surf);

    fprintf('AR = %-4.1f Done.\n', AR(r));
    
    %% 6. PLOTTING 
    % ============
    plotTempField(Surf, T_sol, a, b, 0);
    title(sprintf('Orthotropic Temperature Field (AR=%.1f)', AR(r)));
end

%% ===========
% PRINT TABLE 
%% ===========
fprintf('\nTable 4. Nodal temperatures (K) Isotropic plate (IGA Results)\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Location ');
for r = 1:length(AR)
    fprintf('   AR=%-4.1f',AR(r));
end
fprintf('\n-------------------------------------------------------------\n');

fprintf('A       '); fprintf('%10.2f',TA); fprintf('\n');
fprintf('B       '); fprintf('%10.2f',TB); fprintf('\n');
fprintf('C       '); fprintf('%10.2f',TC); fprintf('\n');
fprintf('D       '); fprintf('%10.2f',TD); fprintf('\n');