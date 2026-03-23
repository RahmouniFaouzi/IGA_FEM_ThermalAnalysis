function patches = holedPlateGeometry(L, R, p, q, refineCount)
    % Generates a multi-patch NURBS mesh for a square plate with a central hole.
    % Uses 4 patches rotated around the center to form the complete domain.
    % Note: The order is fixed to 2*2...
    
    % 1. DEFINE BASE PATCH (Quarter Section)
    % Factor for exact representation of a circular arc using quadratic NURBS
    fac = sind(45);  
    a   = fac * R;
    
    % Initialize 3x3 Control Points (9 points, 4 components: x, y, z, weight)
    controlPts = zeros(9, 4); 
    
    % Define physical coordinates for the base patch (before weighting)
    % This forms a transition from the circular hole (inner) to the square edge (outer)
    base_cp_coords = [ a, a; 
                       0, 2*a; 
                      -a, a; 
                     L/4, L/4; 
                       0, L/2; 
                    -L/4, L/4; 
                     L/2, L/2; 
                       0, L/2; 
                    -L/2, L/2 ];
                    
    % Interpolate middle row to ensure smooth mesh lines
    base_cp_coords(4:6, :) = (base_cp_coords(1:3,:) + base_cp_coords(7:9,:))/2;
    
    % Define weights (corner points of circle arc get weight 'fac')
    w = [1; fac; 1; 1; 1; 1; 1; 1; 1];
    
    % Apply homogeneous coordinates (x*w, y*w, z*w, w)
    controlPts(:, 1) = base_cp_coords(:,1) .* w; 
    controlPts(:, 2) = base_cp_coords(:,2) .* w;
    controlPts(:, 3) = 0; 
    controlPts(:, 4) = w;
    
    % Standard open knot vectors for quadratic patch (p=2)
    uKnot = [0 0 0 1 1 1]; 
    vKnot = [0 0 0 1 1 1];
    
    % Initialize structure to hold data for the 4 patches
    patches = struct('controlPts', [], ...
                     'uKnot', [], ...
                     'vKnot', [], ...
                     'n_u', 0, ...
                     'n_v', 0, ...
                     'element', [], ...
                     'index', [], ...
                     'spans_u', [], ...
                     'spans_v', []);
                     
    % Rotation angles to generate the full plate from the base patch
    angles = [0, 90, 180, 270];
    
    % 2. GENERATE AND REFINE 4 PATCHES
    for i = 1:4
        % Calculate rotation matrix
        theta = deg2rad(angles(i));
        Rot = [cos(theta) -sin(theta); 
               sin(theta) cos(theta)];
        
        controlPts_rot = controlPts;
        
        % Rotate physical coordinates
        coords = bsxfun(@rdivide, controlPts(:,1:2), controlPts(:,4)); % Project to physical
        coords_rot = (Rot * coords')';
        
        % Shift to center of the plate (L/2, L/2)
        coords_rot(:,1) = coords_rot(:,1) + L/2; 
        coords_rot(:,2) = coords_rot(:,2) + L/2;
        
        % Convert back to homogeneous coordinates
        controlPts_rot(:,1) = coords_rot(:,1) .* controlPts_rot(:,4); 
        controlPts_rot(:,2) = coords_rot(:,2) .* controlPts_rot(:,4);
        
        % Apply h-refinement (insert knots) to increase mesh density
        [controlPts_new, uK_new, vK_new] = refinePatchUniform(controlPts_rot, uKnot, vKnot, p, q, refineCount);
        
        % Store patch data
        patches(i).controlPts = controlPts_new; 
        patches(i).uKnot = uK_new; 
        patches(i).vKnot = vK_new;
        patches(i).n_u = length(uK_new) - p - 1; % Number of basis functions in U
        patches(i).n_v = length(vK_new) - q - 1; % Number of basis functions in V
        
        % Generate element connectivity for analysis
        [patches(i).element, patches(i).index, patches(i).spans_u, patches(i).spans_v] = buildConnectivity(patches(i).n_u, patches(i).n_v, p, q, uK_new, vK_new);
    end
end

function [controlPts, U, V] = refinePatchUniform(CP_in, U_in, V_in, p, q, level)
    % Recursively subdivides the mesh by inserting knots at midpoints.
    controlPts = CP_in; 
    U = U_in; 
    V = V_in; 
    nu = 3; % Initial number of control points in U (for 3x3 patch)
    nv = 3; % Initial number of control points in V
    
    for l = 1:level
        % Refine U-direction
        unique_u = unique(U); 
        new_knots = (unique_u(1:end-1) + unique_u(2:end)) / 2; % Find midpoints
        for k = 1:length(new_knots)
            [U, controlPts, nu] = insertKnot_u(U, controlPts, p, nu, nv, new_knots(k)); 
        end
        
        % Refine V-direction
        unique_v = unique(V); 
        new_knots_v = (unique_v(1:end-1) + unique_v(2:end)) / 2;
        for k = 1:length(new_knots_v)
            [V, controlPts, nv] = insertKnot_v(V, controlPts, q, nu, nv, new_knots_v(k)); 
        end
    end
end



function [U_new, CP_new, nu_new] = insertKnot_u(U, controlPts, p, nu, nv, u_val)
    % Inserts a single knot into the U-vector and updates Control Points.
    span  = FindSpan_(nu, p, u_val, U, 0); % Find where the knot fits
    U_new = [U(1:span), u_val, U(span+1:end)]; 
    nu_new = nu + 1;
    CP_new = zeros(nu_new * nv, 4);
    
    % Update control points row by row
    for j = 1:nv
        row_indices = (j-1)*nu + (1:nu); 
        P_row = controlPts(row_indices, :); 
        Q_row = zeros(nu_new, 4);
        
        % Boehm's algorithm for knot insertion
        for i = 1:nu_new
            if i <= span - p
                Q_row(i,:) = P_row(i,:);
            elseif i >= span + 1
                Q_row(i,:) = P_row(i-1,:);
            else
                % Linear interpolation of existing points
                alpha = (u_val - U(i)) / (U(i+p) - U(i)); 
                Q_row(i,:) = (1-alpha)*P_row(i-1,:) + alpha*P_row(i,:); 
            end
        end
        new_row_indices = (j-1)*nu_new + (1:nu_new); 
        CP_new(new_row_indices, :) = Q_row;
    end
end

function [V_new, CP_new, nv_new] = insertKnot_v(V, controlPts, q, nu, nv, v_val)
    % Inserts a single knot into the V-vector and updates Control Points.
    span   = FindSpan_(nv, q, v_val, V, 0);
    V_new  = [V(1:span), v_val, V(span+1:end)]; 
    nv_new = nv + 1;
    CP_new = zeros(nu * nv_new, 4);
    
    % Update control points column by column
    for i = 1:nu
        col_indices = i:nu:(nu*(nv-1)+i); 
        P_col = controlPts(col_indices, :); 
        Q_col = zeros(nv_new, 4);
        
        for j = 1:nv_new
            if j <= span - q
                Q_col(j,:) = P_col(j,:);
            elseif j >= span + 1
                Q_col(j,:) = P_col(j-1,:);
            else
                % Linear interpolation
                alpha      = (v_val - V(j)) / (V(j+q) - V(j)); 
                Q_col(j,:) = (1-alpha)*P_col(j-1,:) + alpha*P_col(j,:); 
            end
        end
        new_col_indices = i:nu:(nu*(nv_new-1)+i); 
        CP_new(new_col_indices, :) = Q_col;
    end
end

function [element, index, spans_u, spans_v] = buildConnectivity(nu, nv, p, q, U, V)
    % Generates the IGA element connectivity table (IEN).
    u_uniq = unique(U); 
    v_uniq = unique(V); 
    num_u = length(u_uniq) - 1; % Number of non-zero spans in U
    num_v = length(v_uniq) - 1; % Number of non-zero spans in V
    
    element = []; 
    index = []; 
    spans_u = zeros(num_u, 2); 
    spans_v = zeros(num_v, 2);
    
    % Identify valid knot spans
    for i=1:num_u, spans_u(i,:) = [u_uniq(i), u_uniq(i+1)]; end
    for i=1:num_v, spans_v(i,:) = [v_uniq(i), v_uniq(i+1)]; end
    
    count = 1;
    % Loop over 2D elements
    for j = 1:num_v
        for i = 1:num_u
            % Find knot span indices
            u_mid = sum(spans_u(i,:))/2; 
            v_mid = sum(spans_v(j,:))/2;
            spanU = FindSpan_(nu, p, u_mid, U, 0); 
            spanV = FindSpan_(nv, q, v_mid, V, 0);
            
            % Collect global indices of all basis functions supporting this element
            sctr = zeros(1, (p+1)*(q+1)); 
            cc = 1;
            for jj = (spanV-q):spanV 
                for ii = (spanU-p):spanU 
                    sctr(cc) = (jj-1)*nu + ii; 
                    cc = cc + 1; 
                end
            end
            element(count, :) = sctr; 
            index(count, :)   = [i, j]; 
            count = count + 1;
        end
    end
end