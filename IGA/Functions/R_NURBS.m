function R = R_NURBS(n_u, n_v, xi_edge, eta, p, q, uKnot,vKnot, weights, sctr)
% Calcuate the NURBS Rational shape Funcion
uspan = FSpan(n_u, p, xi_edge, uKnot);
vspan = FSpan(n_v, q, eta, vKnot);

Nu = BasisFuns(uspan, xi_edge, p, uKnot);
Nv = BasisFuns(vspan, eta, q, vKnot);

R = zeros(length(sctr), 1);
w_sum = 0;
c = 0;
% Tensor Product Loop
for j = 1:q+1
    for i = 1:p+1
        c = c + 1;
        val = Nu(i) * Nv(j) * weights(sctr(c));
        R(c) = val;
        w_sum = w_sum + val;
    end
end
R = R / w_sum;
end

% HELPER FUNCTION
function idx = FSpan(n, p, pts, KntVect)
% Based on Algorithm A2.1 [The NURBS BOOK, p.68]
idx = zeros(size(pts));
for i = 1 : numel(pts)
    if (pts(i) == KntVect(n + 1))
        idx(i) = n;
        return
    end
    low = p + 1;
    high = n + 1;
    mid = floor((low + high)/2);
    while(pts(i) < KntVect(mid) || pts(i) >= KntVect(mid + 1))
        if(pts(i) < KntVect(mid))
            high = mid;
        else
            low = mid;
        end
        mid = floor((low + high) / 2);
    end
    idx(i) = mid;
end
end