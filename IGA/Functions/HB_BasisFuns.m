function [N, dN] = HB_BasisFuns(i, u, p, U)
    N = zeros(1, p+1);
    left = zeros(1,p+1); right = zeros(1,p+1);
    ndu = zeros(p+1, p+1); ndu(1,1) = 1.0;
    for j = 1:p
        left(j+1) = u - U(i+1-j);
        right(j+1) = U(i+j) - u;
        saved = 0.0;
        for r = 0:(j-1)
            ndu(j+1, r+1) = right(r+2) + left(j-r+1);
            temp = ndu(r+1, j) / ndu(j+1, r+1);
            ndu(r+1, j+1) = saved + right(r+2)*temp;
            saved = left(j-r+1)*temp;
        end
        ndu(j+1, j+1) = saved;
    end
    for j = 0:p, N(j+1) = ndu(j+1, p+1); end
    
    eps = 1e-5;
    N_p = HB_BasisVal(i, u+eps, p, U); 
    N_m = HB_BasisVal(i, u-eps, p, U);
    dN  = (N_p-N_m)/(2*eps);
end