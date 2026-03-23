function [cp, U, V] = refine_uniform_dir(cp, U, V, p, q, dir)
    % REFINE_UNIFORM_DIR Performs uniform h-refinement in a specified direction.
    % This function identifies all existing knot spans (elements) and bisects 
    % them by inserting a new knot at each midpoint.
    %
    % Inputs:
    %   cp      : [4 x nu x nv] Control point grid (homogeneous).
    %   U, V    : Current knot vectors.
    %   p, q    : Polynomial degrees.
    %   dir     : Direction to refine (1 for U/Angular, 2 for V/Radial).
    %
    % Outputs:
    %   cp      : Updated control point grid with increased density.
    %   U, V    : Updated knot vectors with midpoints inserted.

    % STEP 1: Select the knot vector associated with the refinement direction
    if dir == 1
        vec = U; 
    else
        vec = V; 
    end

    % STEP 2: Identify unique knots to find element boundaries
    uniq = unique(vec); 

    % STEP 3: Calculate midpoints of all existing knot spans
    mid = (uniq(1:end-1) + uniq(2:end)) / 2;

    % STEP 4: Iteratively insert knots at calculated midpoints
    for k = 1:length(mid)
        % Knot insertion is performed sequentially
        [cp, U, V] = insert_knot_surface(cp, U, V, p, q, mid(k), dir); 
    end
end