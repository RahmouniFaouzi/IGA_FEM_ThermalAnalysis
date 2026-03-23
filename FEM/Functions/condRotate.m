function D = condRotate(kL,kT,theta)
% condRotate  Computes rotated thermal conductivity matrix
%
% INPUTS:
%   kL    - thermal conductivity along local (longitudinal) direction
%   kT    - thermal conductivity along transverse direction
%   theta - rotation angle of material axes w.r.t global axes (radians)
%
% OUTPUT:
%   D     - 2x2 thermal conductivity matrix in global coordinates

% Compute cosine and sine of rotation angle
c = cos(theta); 
s = sin(theta);

% Construct rotated conductivity matrix using standard transformation formula
D = [kL*c^2 + kT*s^2, (kL-kT)*s*c;
     (kL-kT)*s*c,     kL*s^2 + kT*c^2];
end