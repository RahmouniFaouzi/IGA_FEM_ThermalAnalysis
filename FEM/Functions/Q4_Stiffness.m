function Ke = Q4_Stiffness(x,y,kx,ky)
% Computes the element stiffness matrix for a 4-node quadrilateral (Q4)
% heat conduction element using 2x2 Gauss integration

gp = [-1 1]/sqrt(3);      % Gauss points for 2-point Gauss quadrature
Ke = zeros(4);            % Initialize 4x4 element stiffness matrix

for i = 1:2               % Loop over Gauss points in xi-direction
    for j = 1:2           % Loop over Gauss points in eta-direction
        
        xi  = gp(i);      % Natural coordinate xi
        eta = gp(j);      % Natural coordinate eta
        
        % Derivatives of shape functions w.r.t. xi and eta
        % Rows: shape functions N1..N4
        % Col 1: dN/dxi, Col 2: dN/deta
        dN = 1/4 * [
            -(1-eta)  -(1-xi)
             (1-eta)  -(1+xi)
             (1+eta)   (1+xi)
            -(1+eta)   (1-xi)];
        
        % Jacobian matrix J = [dx/dxi  dy/dxi; dx/deta  dy/deta]
        J = [dN(:,1)'*x   dN(:,1)'*y;
             dN(:,2)'*x   dN(:,2)'*y];
        
        % B matrix: temperature gradient matrix (dT/dx, dT/dy)
        B = inv(J) * dN';
        
        % Thermal conductivity matrix (anisotropic case)
        D = [kx  0;
             0  ky];
        
        % Add Gauss point contribution to element stiffness matrix
        Ke = Ke + (B' * D * B) * det(J);
    end
end
end
