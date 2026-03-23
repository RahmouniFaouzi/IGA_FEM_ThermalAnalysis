function Me = T3_Mass(x,y,rho,c)
% -------------------------------------------------------------------------
% T3 CONSISTENT MASS MATRIX (2D HEAT CONDUCTION)
% -------------------------------------------------------------------------
% x, y : nodal coordinates of the triangular element (3 nodes)
% rho  : material density
% c    : specific heat capacity
% Me   : 3x3 element mass matrix
% -------------------------------------------------------------------------

A = polyarea(x,y);          
% Compute the area of the triangular element using its nodal coordinates

% Consistent mass matrix for a linear (T3) triangular element
% Derived from: ?(N^T * rho * c * N) dA
% The coefficient (A/12) comes from exact integration of linear shape
% functions over a triangle
Me = (rho*c*A/12) * ...
     [2 1 1;          % Diagonal terms represent self-capacitance
      1 2 1;          % Off-diagonal terms represent coupling between nodes
      1 1 2];

end