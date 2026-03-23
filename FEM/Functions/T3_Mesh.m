function [node,elem] = T3_Mesh(a,b,nx,ny)
%--------------------------------------------------------------------------
% Generates a structured 2D triangular (T3) mesh over a rectangular domain
%
% Inputs:
%   a  - length of the domain in x-direction
%   b  - length of the domain in y-direction
%   nx - number of divisions in x-direction
%   ny - number of divisions in y-direction
%
% Outputs:
%   node - nodal coordinates [x y] (Nnode x 2)
%   elem - element connectivity for T3 elements (Nelement x 3)
%--------------------------------------------------------------------------

% Create a structured grid of points over the rectangle [0,a] x [0,b]
[x,y] = meshgrid(linspace(0,a,nx+1), linspace(0,b,ny+1));

% Store nodal coordinates in a single array:
% each row corresponds to one node [xi, yi]
node = [x(:) y(:)];

% Initialize element connectivity array
elem = [];

% Loop over each rectangular cell in the structured grid
for j = 1:ny
    for i = 1:nx
        
        % Node numbering of the current rectangular cell
        % n1 ---- n2
        %  |       |
        % n3 ---- n4
        n1 = (j-1)*(nx+1) + i;     % bottom-left node
        n2 = n1 + 1;              % bottom-right node
        n3 = n1 + nx + 1;         % top-left node
        n4 = n3 + 1;              % top-right node
        
        % Split the rectangle into two triangular (T3) elements
        % Triangle 1: (n1, n2, n4)
        % Triangle 2: (n1, n4, n3)
        elem = [elem;
                n1 n2 n4;
                n1 n4 n3];
    end
end

end
