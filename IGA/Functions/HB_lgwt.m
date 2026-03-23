function [x,w] = HB_lgwt(N, a, b)
    N1 = N; 
    N2 = N+1; 
    y  = cos((2*(0:N-1)'+1)*pi/(2*N)); 
    y0 = 2;
    while max(abs(y-y0))>eps
        L = zeros(N1,N2); 
        L(:,1) = 1;
        L(:,2) = y;
        for k = 2:N1
            L(:,k+1) = ((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k; 
        end
        Lp = (N2)*(L(:,N1)-y.*L(:,N2))./(1-y.^2); 
        y0 = y; 
        y = y0-L(:,N2)./Lp;
    end
    x = (a*(1-y)+b*(1+y))/2; 
    w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end