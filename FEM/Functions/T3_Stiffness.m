function Ke = T3_Stiffness(x,y,D)
% Computes the element conductivity (stiffness) matrix for a 2D
% linear triangular (T3) finite element with anisotropic conductivity

A = polyarea(x,y);          % Area of the triangular element

% Coefficients related to derivatives of shape functions w.r.t. y
b = [y(2)-y(3);             % b1
     y(3)-y(1);             % b2
     y(1)-y(2)];            % b3

% Coefficients related to derivatives of shape functions w.r.t. x
c = [x(3)-x(2);             % c1
     x(1)-x(3);             % c2
     x(2)-x(1)];            % c3

% Strain�displacement (temperature gradient) matrix
% Relates nodal temperatures to temperature gradients
B = 1/(2*A)*[b';            % dN/dx terms
             c'];           % dN/dy terms

% Element stiffness (conductivity) matrix
Ke = A*(B'*D*B);            % Integrated over the element area
end
