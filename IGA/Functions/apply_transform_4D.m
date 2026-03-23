function cp_new = apply_transform_4D(cp, M)
    % APPLY_TRANSFORM_4D Applies a 3D transformation matrix to NURBS control points.
    %
    % Inputs:
    %   cp : 3D array [4 x nu x nv] containing control points in weighted 
    %        homogeneous coordinates [wx, wy, wz, w].
    %   M  : [4 x 4] Transformation matrix (Rotation, Translation, Scaling).
    %
    % Output:
    %   cp_new : Transformed control points in weighted homogeneous form.

    % Get dimensions of the control point grid
    [~, nu, nv] = size(cp); 
    cp_new = cp; % Initialize output array

    % Iterate through the grid of control points
    for j = 1:nv
        for i = 1:nu
            % Extract the current control point [wx; wy; wz; w]
            pt = cp(:,i,j); 
            
            % STEP 1: Perspective Projection (De-homogenization)
            % Divide by weight (pt(4)) to get physical 3D coordinates [x; y; z; 1]
            phys = [pt(1)/pt(4); pt(2)/pt(4); pt(3)/pt(4); 1];
            
            % STEP 2: Apply Transformation
            % Multiply by the 4x4 transformation matrix
            phys_new = M * phys; 
            
            % STEP 3: Re-homogenization
            % Multiply the new physical coordinates by the original weight 
            % to maintain the NURBS shape property.
            cp_new(1:3,i,j) = phys_new(1:3) * pt(4);
            
            % Note: pt(4) (the weight) remains unchanged in standard rigid transforms.
        end
    end
end