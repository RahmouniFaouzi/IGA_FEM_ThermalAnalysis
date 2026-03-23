function [phi, dphi_mat, neighs] = MLSShapeFunc(x, y, nodes, dmax_all)
% MLS SHAPE FUNCTION (Full Derivative Implementation)
% -------------------------------------------------------------------------
% Inputs:
%   x, y      : Coordinates of the evaluation point (Gauss point or pixel)
%   nodes     : Global coordinate matrix of all nodes [N x 2]
%   dmax_all  : Vector of support radii (influence domains) for all nodes
%
% Outputs:
%   phi       : Shape function values for neighbor nodes [1 x n]
%   dphi_mat  : Derivatives [dphi/dx; dphi/dy] size [2 x n]
%   neighs    : Indices of the nodes in the support domain
% -------------------------------------------------------------------------

    %% 1. NEIGHBOR SEARCH
    % Calculate squared distance from evaluation point (x,y) to all nodes
    d2 = (nodes(:,1)-x).^2 + (nodes(:,2)-y).^2;
    
    % Optimization: First filter using the maximum possible support radius
    max_d = max(dmax_all);
    neighs = find(d2 < max_d^2);

    % Strict filter: Check if point is within the specific dmax of each neighbor
    if ~isempty(neighs)
        valid = d2(neighs) < dmax_all(neighs).^2;
        neighs = neighs(valid);
    end

    % Stability check: We need enough neighbors to invert the moment matrix A.
    % For a linear basis in 2D, min is 3, but 6+ is recommended for stability.
    n = length(neighs);
    if n < 6, phi=[]; dphi_mat=[]; return; end

    % Extract local data for neighbors
    xn = nodes(neighs, 1); 
    yn = nodes(neighs, 2);
    dm = dmax_all(neighs);

    %% 2. BASIS FUNCTIONS (Linear & Centered)
    % Linear Basis p = [1, x, y]. 
    % We shift the basis to the evaluation point (x,y) to improve conditioning.
    % p(x_I) = [1, x_I - x, y_I - y]
    p = [ones(n,1), xn-x, yn-y];
    
    % The basis evaluated at the point (x,y) itself becomes [1, 0, 0]
    p_ev = [1, 0, 0];

    % Spatial derivatives of the basis vector p w.r.t x and y
    % dp/dx = [0, -1, 0] and dp/dy = [0, 0, -1]
    dp_dx_vec = [0, -1, 0];
    dp_dy_vec = [0, 0, -1];

    %% 3. WEIGHT FUNCTIONS (Cubic Spline)
    w = zeros(n,1);
    dw_dx = zeros(n,1);
    dw_dy = zeros(n,1);

    for k = 1:n
        % Vector from evaluation point TO node k
        dx = xn(k)-x; dy = yn(k)-y;
        
        r2 = dx^2 + dy^2;
        r = sqrt(r2) / dm(k); % Normalized distance (0 to 1)
        
        % Cubic Spline Weight Function & its derivative w.r.t r (dv)
        if r <= 0.5
            val = 2/3 - 4*r^2 + 4*r^3;
            dv  = -8*r + 12*r^2;
        elseif r <= 1
            val = 4/3 - 4*r + 4*r^2 - (4/3)*r^3;
            dv  = -4 + 8*r - 4*r^2;
        else
            val=0; dv=0;
        end
        w(k) = val;
        
        % Chain Rule for spatial derivatives: dw/dx = (dw/dr) * (dr/dx)
        % dr/dx = -(xn - x) / (dm * dist)
        if r > 0
            term = dv / (dm(k) * sqrt(r2));
            dw_dx(k) = term * (-(xn(k)-x));
            dw_dy(k) = term * (-(yn(k)-y));
        end
    end

    %% 4. MOMENT MATRIX (A) AND SOLVE
    W = diag(w);      % Weight matrix
    A = p' * W * p;   % Moment matrix: Sum( w * p * p^T )

    % Tikhonov Regularization: Ensure A is invertible even if nodes are collinear
    if rcond(A) < 1e-12, A = A + eye(3)*1e-6; end

    B = p' * W;       % The 'B' matrix part of MLS
    iA = inv(A);      % Inverse of A

    % Standard MLS Shape Function: phi = p(0)^T * A^-1 * B
    phi = p_ev * iA * B;

    %% 5. FULL DERIVATIVES (dphi/dx, dphi/dy)
    % We compute the exact derivative including d(A^-1)/dx and dB/dx terms.
    % This is more accurate than the "Diffuse" approximation for heat transfer.
    
    A_x = zeros(3,3); A_y = zeros(3,3);
    B_x = zeros(3,n); B_y = zeros(3,n);

    for k=1:n
        pk = p(k,:)'; % Column vector of basis at node k
        
        % Basis derivatives at node k (constant for linear basis)
        dp_dx_k = [0; -1; 0];
        dp_dy_k = [0; 0; -1];
        
        % --- Derivative of Moment Matrix A ---
        % A = Sum(w * p * p^T)
        % dA/dx = Sum( dw/dx * p * p^T  +  w * dp/dx * p^T  +  w * p * dp/dx^T )
        
        term1_x = dw_dx(k) * (pk * pk');
        term1_y = dw_dy(k) * (pk * pk');
        
        term2_x = w(k) * (dp_dx_k * pk');
        term2_y = w(k) * (dp_dy_k * pk');
        
        term3_x = term2_x'; % Symmetry
        term3_y = term2_y';
        
        A_x = A_x + term1_x + term2_x + term3_x;
        A_y = A_y + term1_y + term2_y + term3_y;
        
        % --- Derivative of B Matrix ---
        % B = [w*p]
        % dB/dx = dw/dx * p + w * dp/dx
        B_x(:,k) = dp_dx_k * w(k) + pk * dw_dx(k);
        B_y(:,k) = dp_dy_k * w(k) + pk * dw_dy(k);
    end

    % Derivative of Inverse Matrix: d(A^-1) = -A^-1 * dA * A^-1
    iA_x = -iA * A_x * iA;
    iA_y = -iA * A_y * iA;

    % Final Chain Rule for Shape Function Derivatives
    % dphi = p_ev * ( d(A^-1)*B + A^-1*dB ) 
    % (Note: p_ev is constant [1,0,0], so its derivative is 0)
    dphi_dx = p_ev * (iA_x * B + iA * B_x);
    dphi_dy = p_ev * (iA_y * B + iA * B_y);

    dphi_mat = [dphi_dx; dphi_dy];
end