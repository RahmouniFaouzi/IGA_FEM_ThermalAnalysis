function Q = knot_insert_curve_4D(P, p, U, u_val)
    % KNOT_INSERT_CURVE_4D Implementation of Boehm's Algorithm for NURBS.
    % This function computes the new set of control points after inserting a 
    % knot into a NURBS curve. It operates in 4D space to handle weights (wx, wy, wz, w).
    %
    % Inputs:
    %   P     : [4 x num_cp] Matrix of control points in homogeneous space.
    %   p     : Polynomial degree of the curve.
    %   U     : The current knot vector.
    %   u_val : The parametric value where the knot is to be inserted.
    %
    % Output:
    %   Q     : [4 x (num_cp+1)] The new set of control points.

    % STEP 1: Find the affected knot span index 'k'
    k = find(U <= u_val, 1, 'last'); 

    % Get dimensions
    [dim, num_cp] = size(P); 
    Q = zeros(dim, num_cp + 1); % Knot insertion adds exactly one control point

    % STEP 2: Copy unaffected control points (Left side)
    for i = 1:(k-p)
        Q(:,i) = P(:,i); 
    end

    % STEP 3: Copy unaffected control points (Right side)
    for i = (k):num_cp
        Q(:,i+1) = P(:,i); 
    end

    % STEP 4: Calculate new control points (The blending zone)
    for i = (k-p+1):k
        % Calculate interpolation weight 'alpha' based on knot spans
        if (U(i+p) - U(i)) == 0
            alpha = 0; % Avoid division by zero at multiple knots
        else
            alpha = (u_val - U(i)) / (U(i+p) - U(i)); 
        end
        
        % Compute the new control point in 4D homogeneous space
        Q(:,i) = (1 - alpha) * P(:,i-1) + alpha * P(:,i);
    end
end