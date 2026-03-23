clc, clear, close all
addpath(genpath('../Functions'));

% ====================================================================
% NURBS: ORTHOTROPIC SQUARE PLATE WITH CENTRAL HOLE
%        (Convergence Study)
%  Author : Dr. Rahmouni Faouzi  | Email : rahmounifaouzi01@gmail.com
%           Pr. Khennane Amar    | Email : a.khennane@adfa.edu.au
% ====================================================================

% ------------------------ 
% 1. SIMULATION PARAMETERS
% ------------------------ 
L   = 50;        % Physical dimension of the plate (Width & Height)
R   = 4;       % Radius of the central hole
kx  = 1003;       % Thermal conductivity in X-direction
ky  = 171;        % Thermal conductivity in Y-direction
D_mat = [kx 0;    % Material Matrix (Orthotropic properties)
         0 ky]; 

T_hole  = 773;    % Dirichlet BC: Fixed Temperature at the hole (Hot)
T_outer = 273;    % Dirichlet BC: Fixed Temperature at outer edge (Cold)

% Refine the mesh (h-refinement)
refinement_levels = [6];

% Plot at level 3
plotLevel = 1;

% Pre-allocate arrays for results
res_elems = zeros(length(refinement_levels), 1);
res_T     = zeros(length(refinement_levels), 4);
res_err   = zeros(length(refinement_levels), 1);

fprintf('Running NURBS Simulation...\n');

% -------------------------------------------------------------------------
% LOOP OVER MESH REFINEMENTS
% -------------------------------------------------------------------------
for i = 1:length(refinement_levels)
    refineCount = refinement_levels(i);
    fprintf('  Processing Refinement Level: %d\n', refineCount);

    % =====================================================================
    % 1. GEOMETRY GENERATION
    % =====================================================================
    fac = sind(45); % Factor to place control points exactly on the circle arc
    a   = fac * R;
    
    % Control Points
    controlPts = zeros(4,3,3); 
    
    % --- Row 1: The Inner Circular Arc ---
    controlPts(1:2,1,1) = [a; a];     % Start of arc
    controlPts(1:2,2,1) = [0; 2*a];   % Middle control point
    controlPts(1:2,3,1) = [-a; a];    % End of arc
    
    % --- Row 3: The Outer Square Boundary ---
    controlPts(1:2,1,3) = [L/2; L/2];
    controlPts(1:2,2,3) = [0; L/2];
    controlPts(1:2,3,3) = [-L/2; L/2];
    
    % --- Row 2: Mid-points ---
    controlPts(1:3,:,2) = (controlPts(1:3,:,1) + controlPts(1:3,:,3)) / 2;
    
    % Shift Geometry to Center (L/2, L/2)
    for r = 1:3
        for c = 1:3
            controlPts(1:2,r,c) = controlPts(1:2,r,c) + [L/2; L/2];
        end
    end
    
    % Apply Homogeneous Coordinates (Weights)
    controlPts(4,:,:)   = 1;          % Default weight
    controlPts(4,2,1)   = fac;        % Special weight for the circle arc
    controlPts(1:2,2,1) = controlPts(1:2,2,1) * fac; % Pre-multiply x,y by w
    
    % Knot Vectors: order p=2 (quadratic)
    uKnot = [0 0 0 1 1 1]; 
    vKnot = [0 0 0 1 1 1]; 
    
    % Create the first patch (North)
    solid1 = nrbmak(controlPts, {uKnot, vKnot});
    
    % Create the other 3 patches by rotating the first one 90 degrees 3 times.
    T_cen = [1 0 0 -L/2; 0 1 0 -L/2; 0 0 1 0; 0 0 0 1];
    R_90  = [0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
    T_back= [1 0 0 L/2; 0 1 0 L/2; 0 0 1 0; 0 0 0 1];
    M_trans = T_back * R_90 * T_cen;
    
    solid2 = nrbtform(solid1, M_trans);
    solid3 = nrbtform(solid2, M_trans);
    solid4 = nrbtform(solid3, M_trans);
    
    % Perform h-Refinement (Knot Insertion)
    if refineCount > 0
        solid1 = hRefineNURBS(solid1, refineCount);
        solid2 = hRefineNURBS(solid2, refineCount);
        solid3 = hRefineNURBS(solid3, refineCount);
        solid4 = hRefineNURBS(solid4, refineCount);
    end
    
    % Convert standard NURBS toolbox objects to a struct easier to loop over
    patches(1) = convert2DNurbsToPatch(solid1);
    patches(2) = convert2DNurbsToPatch(solid2);
    patches(3) = convert2DNurbsToPatch(solid3);
    patches(4) = convert2DNurbsToPatch(solid4);
    
    % =====================================================================
    % 2. CONNECTIVITY
    % =====================================================================
    all_CPs = [];
    for p=1:4
        cp = patches(p).controlPts;
        % Project to physical space (x = wx / w) to compare locations
        cp_geom = [cp(:,1)./cp(:,4), cp(:,2)./cp(:,4)]; 
        all_CPs = [all_CPs; cp_geom];
    end
    [unique_nodes, ~, ic] = unique(round(all_CPs, 5), 'rows');
    ndof = size(unique_nodes, 1); % Total number of unique unknowns
    
    % Store the Global ID pattern for each patch for easy assembly later
    global_node_patterns = cell(4,1);
    curr = 0;
    for p=1:4
        n_p = size(patches(p).controlPts,1);
        indices_local = (1:n_p)' + curr;
        global_indices = ic(indices_local);
        global_node_patterns{p} = reshape(global_indices, patches(p).noPtsY, patches(p).noPtsX);
        curr = curr + n_p;
    end
    res_elems(i) = sum([patches.noElemsU] .* [patches.noElemsV]);

    % =====================================================================
    % 3. ASSEMBLY
    % =====================================================================
    K = sparse(ndof, ndof); % Stiffness Matrix
    F = zeros(ndof, 1);     % Force Vector
    
    % Get Gauss Points
    [W, Q_gauss] = quadrature(3, 'GAUSS', 2);
    
    % Loop over every Patch -> Element -> Gauss Point
    for p = 1:4
        patch = patches(p);
        for e = 1:size(patch.element, 1)
            % Get global indices (DOFs) for this element
            local_sctr = patch.element(e, :);
            glob_sctr  = global_node_patterns{p}(local_sctr); 
            pts        = patch.controlPts(local_sctr, :); 
            
            % Get Parametric Range
            idx_u = patch.index(e,1); 
            idx_v = patch.index(e,2);
            xiE   = patch.elRangeU(idx_u, :); 
            etaE  = patch.elRangeV(idx_v, :);
            
            for gp = 1:length(W)
                wt = W(gp);
                
                % Map Parent -> Parametric
                Xi  = parent2ParametricSpace(xiE, Q_gauss(gp,1));
                Eta = parent2ParametricSpace(etaE, Q_gauss(gp,2));
                J2  = jacobianPaPaMapping(xiE,etaE); % Jacobian of this step

                % Calculate NURBS Basis functions (R) and derivatives (dR/dxi)
                [~, dRdxi, dRdeta] = NURBS2DBasisDers([Xi; Eta], ...
                    patch.p, patch.q, patch.uKnot, patch.vKnot, pts(:,4));
                
                % Calculate Jacobian: Parametric -> Physical
                CP_element_phys = [pts(:,1)./pts(:,4), pts(:,2)./pts(:,4)];
                Jac = CP_element_phys' * [dRdxi, dRdeta]; 
                
                detJ = det(Jac);
                % Transform derivatives to Physical Space
                dRdx = [dRdxi, dRdeta] * inv(Jac);
                B = dRdx';
                
                % Add to Stiffness Matrix (K = B' * D * B * detJ * weight)
                K(glob_sctr, glob_sctr) = K(glob_sctr, glob_sctr) + B' * D_mat * B * detJ * J2 * wt;
            end
        end
    end
    
    % =====================================================================
    % 4. BOUNDARY CONDITIONS (Dirichlet)
    % =====================================================================
    n_hole = []; n_outer = [];
    % Identify nodes on the hole and outer edge based on the patch topology
    for p = 1:4
        hole_nodes_p  = global_node_patterns{p}(1, :);   % Bottom edge of patch = Hole
        outer_nodes_p = global_node_patterns{p}(end, :); % Top edge of patches = Outer
        n_hole  = [n_hole; hole_nodes_p(:)];
        n_outer = [n_outer; outer_nodes_p(:)];
    end
    n_hole  = unique(n_hole); 
    n_outer = unique(n_outer);
    
    fixed    = unique([n_hole; n_outer]);
    free     = setdiff(1:ndof, fixed);
    
    % Apply Temperatures
    T = zeros(ndof, 1);
    T(n_hole)  = T_hole;
    T(n_outer) = T_outer;
    
    % Solve the linear system for free nodes: K_free * T_free = F - K_fixed * T_fixed
    T(free) = K(free, free) \ (F(free) - K(free, fixed)*T(fixed));
    
    % 5. EVALUATION AT EXACT POINTS (Inverse Mapping)
    % ==============================================
    patch_id = 1; % North Patch
    % Load target physical points generated by the FEM code
    pts      = load('Exact_Nodes.mat');
    pts_eval = pts.final_pts;
            
    u_fixed = 0.5; % The points lie on the symmetry line (middle of patch)
    
    for pt_idx = 1:4
        y_target = pts_eval(pt_idx, 2);
        
        % Use Newton-Raphson to solve the non-linear NURBS mapping equation
        v_guess = (y_target - (L/2+R)) / (L/2 - R); 
        v_guess = max(0, min(1, v_guess)); % Clamp to [0,1]
        
        v_sol = find_v_for_y(v_guess, u_fixed, y_target, patches(patch_id));
        
        % Once we have (u, v_sol), we can evaluate the Temperature field there
        [val, phys] = eval_nurbs_field(u_fixed, v_sol, patches(patch_id), global_node_patterns{patch_id}, T);
        res_T(i, pt_idx) = val;
    end
    
    % Calculate Flux Error (Balance of heat entering vs leaving)
    Q_full = K * T; 
    res_err(i) = abs(sum(Q_full(n_hole)) + sum(Q_full(n_outer)));
    
    % Visualize the result for the first refinement level
    if (i == plotLevel)
        plotHolesPlate(patches, global_node_patterns, T)
    end
end

% RESULTS TABLE ---
fprintf('Display IGA Results... \n');
fprintf('\n=============================================================\n');
fprintf(' RESULTS TABLE (NURBS)\n');
fprintf('=============================================================\n');
fprintf('%-10s %-10s %-10s %-10s %-10s %-10s\n', 'Elements', 'A', 'B', 'C', 'D', 'Error');
fprintf('-------------------------------------------------------------\n');
for i = 1:length(refinement_levels)
    fprintf('%-10d %-10.2f %-10.2f %-10.2f %-10.2f %-10.2e\n', ...
        res_elems(i), res_T(i,1), res_T(i,2), res_T(i,3), res_T(i,4), res_err(i));
end

