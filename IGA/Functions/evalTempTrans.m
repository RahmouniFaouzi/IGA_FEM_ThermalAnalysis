function T_val = evalTempTrans(patch, glob_idx, T, xi, eta, p, q)
    % Evaluates the field variable T at parametric coordinates (xi, eta).
    
    % Find knot spans and evaluate basis functions
    spanU = FindSpan_(patch.n_u, p, xi, patch.uKnot, 0);
    spanV = FindSpan_(patch.n_v, q, eta, patch.vKnot, 0);
    Nu = BasisFuns(spanU, xi, p, patch.uKnot);
    Nv = BasisFuns(spanV, eta, q, patch.vKnot);
    
    num = 0; den = 0;
    % Loop over non-zero basis functions
    for jj = (spanV-q):spanV
        for ii = (spanU-p):spanU
            % Global index in the patch control point array
            idx = (jj-1)*patch.n_u + ii;
            
            % Local index for basis function array access
            loc_i = ii - (spanU-p) + 1; 
            loc_j = jj - (spanV-q) + 1;
            
            % Compute weighted basis function value
            val = Nu(loc_i) * Nv(loc_j);
            w = patch.controlPts(idx, 4);
            basis = val * w;
            
            % Accumulate numerator (basis * value) and denominator (basis weights)
            den = den + basis;
            gid = glob_idx(idx);
            num = num + basis * T(gid);
        end
    end
    % Normalize by the sum of weights (NURBS projection)
    T_val = num / den;
end