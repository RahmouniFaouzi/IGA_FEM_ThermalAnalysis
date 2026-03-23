function S = Subdivision(U_coarse, p)
    % Oslo Algorithm / Knot Insertion
    nc_c = length(U_coarse) - p - 1;
    S = eye(nc_c);
    U_curr = U_coarse;
    
    % Identify knots to insert
    % Compare U_coarse and U_fine. U_fine has extra midpoints.
    to_insert = [];
    for k=1:(length(U_curr)-1)
        if U_curr(k) < U_curr(k+1)
            mid = (U_curr(k) + U_curr(k+1))/2;
            to_insert = [to_insert, mid];
        end
    end
    
    % Apply Insertion
    for u_ins = to_insert
        k = find(U_curr <= u_ins, 1, 'last');
        n = length(U_curr) - p - 1;
        T = zeros(n+1, n);
        
        for j = 1:(n+1)
            if j <= k-p 
                T(j,j)=1;
            elseif j >= k+1 
                T(j, j-1)=1;
            else
                denom = U_curr(j+p) - U_curr(j);
                if denom == 0 
                    a=0; 
                else
                    a = (u_ins - U_curr(j))/denom; 
                end
                T(j, j-1) = 1 - a;
                T(j, j)   = a;
            end
        end
        S = T * S;
        U_curr = sort([U_curr, u_ins]);
    end
end
