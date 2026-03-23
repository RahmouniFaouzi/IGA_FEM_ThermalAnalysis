function [U_new, n_new] = insert_knot_vec(U, val)
    % INSERT_KNOT_VEC Updates a knot vector by inserting a new parametric value.
    % This is a fundamental step in h-refinement, providing the new parameter 
    % boundaries for the refined elements.
    %
    % Inputs:
    %   U   : Original knot vector (monotonically non-decreasing).
    %   val : The new parametric coordinate to be inserted (must be within U's range).
    %
    % Outputs:
    %   U_new : The expanded knot vector containing the new value.
    %   n_new : The new length of the knot vector.

    % STEP 1: Find the insertion index 'k'
    k = find(U <= val, 1, 'last'); 

    % STEP 2: Concatenate the vector
    U_new = [U(1:k), val, U(k+1:end)]; 

    % STEP 3: Return the new size
    n_new = length(U_new); 
end