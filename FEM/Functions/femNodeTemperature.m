function Tq = femNodeTemperature(node, XP, T, tol)
% INTERP_T  Interpolates nodal temperature at a given point
%
% INPUTS:
%   node - Nx2 array of nodal coordinates [x, y]
%   XP   - 1x2 vector, query point [x, y] where temperature is needed
%   T    - Nx1 vector of nodal temperatures
%   tol  - tolerance for matching coordinates (to account for numerical errors)
%
% OUTPUT:
%   Tq   - interpolated temperature at XP (NaN if point not found or ambiguous)

% Find node(s) whose coordinates match XP within tolerance
idXP = find( abs(node(:,1)-XP(1)) < tol & ...
             abs(node(:,2)-XP(2)) < tol );

% Handle cases based on number of matching nodes
if isempty(idXP)
    % No node found near the query point
    Tq = NaN;
elseif length(idXP) > 1
    % Multiple nodes found at the query point (ambiguous)
    warning('Multiple nodes matched point X.');
    Tq = NaN;
else
    % Single matching node found, return its temperature
    Tq = T(idXP(1));
end
end
