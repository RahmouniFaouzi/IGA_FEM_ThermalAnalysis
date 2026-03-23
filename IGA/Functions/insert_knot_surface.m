function [cp_new, U, V] = insert_knot_surface(cp, U, V, p, q, val, dir)
    % INSERT_KNOT_SURFACE Performs knot insertion on a NURBS surface patch.
    % This process increases the number of elements and control points (h-refinement)
    % without changing the physical shape or the continuity of the geometry.
    %
    % Inputs:
    %   cp   : [4 x nu x nv] Original control point grid (homogeneous coordinates).
    %   U, V : Current knot vectors in the u and v directions.
    %   p, q : Polynomial degrees in the u and v directions.
    %   val  : The parametric value [0, 1] where the new knot will be inserted.
    %   dir  : Direction of insertion (1 for U-direction, 2 for V-direction).
    %
    % Outputs:
    %   cp_new : Updated control point grid (size increases by 1 in the insertion direction).
    %   U, V   : Updated knot vectors containing the new value.

    % Get current dimensions of the control point grid
    [dim, nu, nv] = size(cp);

    if dir == 1  % --- INSERTION IN U-DIRECTION (Angular/Horizontal) ---
        
        % Pre-allocate new grid (adding one basis function in U)
        nu_new = nu + 1; 
        cp_new = zeros(dim, nu_new, nv);
        
        % Apply curve knot insertion to every "row" of control points in V
        for j = 1:nv
            % squeeze(cp(:,:,j)) extracts a 2D curve of CPs at constant V
            cp_new(:,:,j) = knot_insert_curve_4D(squeeze(cp(:,:,j)), p, U, val); 
        end
        
        % Update the U knot vector with the new value
        [U, ~] = insert_knot_vec(U, val);
        
    else         % --- INSERTION IN V-DIRECTION (Radial/Vertical) ---
        
        % Pre-allocate new grid (adding one basis function in V)
        nv_new = nv + 1; 
        cp_new = zeros(dim, nu, nv_new);
        
        % Apply curve knot insertion to every "column" of control points in U
        for i = 1:nu
            % squeeze(cp(:,i,:)) extracts a 2D curve of CPs at constant U
            cp_new(:,i,:) = knot_insert_curve_4D(squeeze(cp(:,i,:)), q, V, val); 
        end
        
        % Update the V knot vector with the new value
        [V, ~] = insert_knot_vec(V, val);
    end
end