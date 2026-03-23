function [xi, eta, success] = invMapNewton(patch, pt_phys)
    % Solves for parametric (xi, eta) given a physical point (pt_phys)
    % using the Newton-Raphson iterative method.
    
    % Initialize guess at center of parametric domain
    xi  = 0.5; 
    eta = 0.5; 
    tol = 1e-4; 
    max_iter = 20; 
    success  = false;
    
    for iter = 1:max_iter
        % Locate current knot spans
        spanU = FindSpan_(patch.n_u, 2, xi, patch.uKnot, 0);
        spanV = FindSpan_(patch.n_v, 2, eta, patch.vKnot, 0);
        
        % Evaluate basis functions and derivatives
        [Nu, dNu] = HB_BasisFuns(spanU, xi, 2,  patch.uKnot);
        [Nv, dNv] = HB_BasisFuns(spanV, eta, 2,  patch.vKnot);
        
        % Collect local indices for the 3x3 control point support
        sctr = [];
        for jj = (spanV-2):spanV
            for ii = (spanU-2):spanU
                sctr = [sctr, (jj-1)*patch.n_u + ii];
            end
        end
        
        % Retrieve physical coordinates of local control points
        wts = patch.controlPts(sctr, 4);
        pts = bsxfun(@rdivide, patch.controlPts(sctr, 1:2), wts);
        
        % Initialize NURBS basis (R) and derivative arrays
        n_en  = 9; 
        R     = zeros(n_en, 1); 
        dRdxi = zeros(n_en, 1); 
        dRdeta = zeros(n_en, 1);
        sum_w  = 0; 
        sum_dw_xi  = 0; 
        sum_dw_eta = 0;
        
        % Compute NURBS basis functions and gradients
        k = 0;
        for j = 1:3
            for i = 1:3
                k = k+1;
                % Tensor product of B-splines
                val = Nu(i)*Nv(j); 
                dval_xi  = dNu(i)*Nv(j); 
                dval_eta = Nu(i)*dNv(j);
                
                % Weighted basis calculations
                w    = wts(k);
                R(k) = val*w; 
                dRdxi(k)  = dval_xi*w; 
                dRdeta(k) = dval_eta*w;
                
                % Accumulate sums for quotient rule
                sum_w = sum_w+R(k); 
                sum_dw_xi  = sum_dw_xi+dRdxi(k); 
                sum_dw_eta = sum_dw_eta+dRdeta(k);
            end
        end
        % Apply quotient rule for NURBS derivatives
        inv_w   = 1/sum_w;
        dR_dxi  = (dRdxi*sum_w - R*sum_dw_xi)*inv_w^2;
        dR_deta = (dRdeta*sum_w - R*sum_dw_eta)*inv_w^2;
        R = R*inv_w;
        
        % Calculate Residual: vector from current mapped point to target point
        cur_pt = (R' * pts)';
        res = cur_pt - pt_phys';
        
        % Convergence check
        if norm(res) < tol
            success = true; return;
        end
        
        % Form Jacobian matrix (dR/dxi * P, dR/deta * P)
        Jac = [dR_dxi'; dR_deta'] * pts;
        
        % Solve linear system for update step (J * delta = res)
        delta = Jac \ res; 
        xi = xi - delta(1);
        eta = eta - delta(2);
        
        % Keep solution within valid parametric domain [0,1]
        xi  = max(0, min(1, xi)); 
        eta = max(0, min(1, eta));
    end
end