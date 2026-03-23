function [node, global_node_patterns, ndof] = mergePatches(patches)
    % Merges separate NURBS patches into a single mesh by identifying shared control points.
    
    % 1. COLLECT ALL PHYSICAL CONTROL POINTS
    all_CPs = [];
    for k = 1:4 
        cp = patches(k).controlPts; 
        % Project from homogeneous (weighted) space to physical space (x = wx/w)
        phys = bsxfun(@rdivide, cp(:,1:2), cp(:,4));
        all_CPs = [all_CPs; phys]; 
    end
    
    % 2. IDENTIFY UNIQUE NODES (GEOMETRIC MERGE)
    tol = 1e-6; 
    % 'unique' finds distinct coordinates; 'ic' maps original points to the unique list
    [node, ~, ic] = unique(round(all_CPs/tol)*tol, 'rows');
    ndof = size(node, 1); % Total number of unique Degrees of Freedom
    
    % 3. BUILD LOCAL-TO-GLOBAL CONNECTIVITY MAP
    global_node_patterns = cell(4,1); 
    curr = 0;
    for k = 1:4
        n_pts = size(patches(k).controlPts, 1); 
        % Assign global IDs to the current patch based on the unique list
        global_node_patterns{k} = ic((1:n_pts)' + curr); 
        curr = curr + n_pts; 
    end
end