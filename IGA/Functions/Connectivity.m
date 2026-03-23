function [elRange, elConn] = Connectivity(p, knots)
% BUILDCONNECTIVITY Generates Element Ranges and Connectivity
%   Inputs:
%       p     : Polynomial degree
%       knots : Knot vector (e.g., [0 0 0 1 1 1])
%
%   Outputs:
%       elRange : [start, end] parametric bounds for each element
%       elConn  : Indices of the Control Points supporting each element

    % 1. IDENTIFY ELEMENTS
    unique_knots = unique(knots); 
    num_elems = length(unique_knots) - 1;
    
    % Initialize arrays
    elRange = zeros(num_elems, 2);    % Store [u_start, u_end]
    elConn  = zeros(num_elems, p+1);  % Store [idx1, idx2, ... idx_p+1]
    
    count = 1;
    
    % 2. LOOP THROUGH KNOT VECTOR
    for i = 1:length(knots)-1
        
        % Check for non-zero span (Element)
        if knots(i+1) > knots(i) + 1e-10
            
            % A. Store Parametric Range
            elRange(count, :) = [knots(i), knots(i+1)]; 
            
            % B. Determine Connectivity
            elConn(count, :) = (i-p):i; 
            
            count = count + 1;
        end
    end
end