function [TR, node, elem, ndof] = triMeshGen(n_circ, n_rad, bias, R, L)
% TRIMESHGEN Generates a structured triangular mesh for a square plate with a hole.
% Inputs: n_circ circumferential divisions
%         n_rad  radial divisions
%         bias   mesh packing
%         R      hole radius
%         L      plate width

% Parametric Space
theta = linspace(0, pi/4, n_circ/8 + 1);

% Radial geometric spacing (Biased towards 0)
xi = linspace(0, 1, n_rad)';
xi = (exp(bias*xi) - 1) / (exp(bias) - 1);

node = [];
elem = [];

% Generate 4 Quadrants
for quad = 1:4
    base_angle = (quad-1)*pi/2;
    for sect = 1:2
        if sect == 1
            th_loc = theta;
        else
            th_loc = fliplr(pi/2 - theta);
        end
        
        th_glob = base_angle + th_loc;
        
        for i = 1:length(th_glob)-1
            t1 = th_glob(i);
            t2 = th_glob(i+1);
            
            % Inner Points (Circle)
            P1_in = [R*cos(t1), R*sin(t1)];
            P2_in = [R*cos(t2), R*sin(t2)];
            
            % Outer Points (Square Boundary)
            P1_out = get_box_intersect(t1, L);
            P2_out = get_box_intersect(t2, L);
            
            % Generate Radial Lines
            for j = 1:n_rad-1
                r_a = xi(j);
                r_b = xi(j+1);
                
                % Four nodes of a logical quadrilateral
                N1 = (1-r_a)*P1_in + r_a*P1_out;
                N2 = (1-r_a)*P2_in + r_a*P2_out;
                N3 = (1-r_b)*P2_in + r_b*P2_out;
                N4 = (1-r_b)*P1_in + r_b*P1_out;
                
                % Add nodes to list (simple method, duplicates removed later)
                start_id = size(node,1);
                node = [node; N1; N2; N3; N4];
                
                % Split Quad into 2 Triangles
                elem = [elem; start_id+1, start_id+2, start_id+3;
                    start_id+1, start_id+3, start_id+4];
            end
        end
    end
end

% Shift Origin to Center of Plate (0..L -> L/2 is center)
node(:,1) = node(:,1) + L/2;
node(:,2) = node(:,2) + L/2;

% Cleanup: Remove Duplicate Nodes (Weld the patches)
[node, ~, idx] = unique(node, 'rows');
elem = idx(elem);
ndof = size(node,1);

% Create Triangulation Object
TR = triangulation(elem, node);
end

function P = get_box_intersect(theta, L)
% Ray casting to find square intersection
half_L = L/2;
t = mod(theta, 2*pi);

if (t >= 7*pi/4 || t < pi/4)       % Right Face
    x = half_L;
    y = half_L * tan(t);
elseif (t >= pi/4 && t < 3*pi/4)   % Top Face
    y = half_L;
    x = half_L / tan(t);
elseif (t >= 3*pi/4 && t < 5*pi/4) % Left Face
    x = -half_L;
    y = -half_L * tan(t);
else                               % Bottom Face
    y = -half_L;
    x = -half_L / tan(t);
end
P = [x, y];
end