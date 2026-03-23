function [node,elem] = Q4_Mesh(a,b,nx,ny)
% Q4_Mesh
% Generates a structured quadrilateral (Q4) finite element mesh
% over a rectangular domain [0,a] x [0,b].
%
% INPUTS:
% a  : length of the domain in x-direction
% b  : length of the domain in y-direction
% nx : number of elements along x-direction
% ny : number of elements along y-direction
%
% OUTPUTS:
% node : (Nnode x 2) array of nodal coordinates [x y]
% elem : (Nelem x 4) array of element connectivity (Q4 elements)

% Generate a regular grid of nodal coordinates
% linspace creates equally spaced points in x and y directions
[x,y] = meshgrid(linspace(0,a,nx+1), linspace(0,b,ny+1));

% Store nodal coordinates in a single array
% Each row corresponds to one node: [x_i, y_i]
node = [x(:) y(:)];

% Initialize element connectivity matrix
elem = [];

% Loop over each rectangular cell in the structured grid
for j = 1:ny
    for i = 1:nx
        
        % Node numbering for the current quadrilateral element
        % n1 ---- n2
        %  |       |
        % n4 ---- n3
        %
        % Global node indices based on structured grid ordering
        n1 = (j-1)*(nx+1)+i;      % bottom-left node
        n2 = n1+1;                % bottom-right node
        n3 = n2+nx+1;             % top-right node
        n4 = n1+nx+1;             % top-left node
        
        % Append current element connectivity (Q4)
        elem = [elem; n1 n2 n3 n4];
    end
end
end
