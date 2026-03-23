function span = FindSpan_(n, p, u, U, L)
% FINDSPAN_ Binary search to find the knot span index
%   Because knot vectors are sorted, Binary Search is O(log n), making it
%   much faster than a linear check for large meshes.

    if nargin == 4, L = 1;end
    % Special case: If u is exactly at the end of the domain
    if L
        if u >= U(n+2), span = n+1; return; end
        if u < U(p+2), span = p+1; return; end
        % Standard Binary Search
        low  = p+1; 
        high = n+2; 
        mid  = floor((low+high)/2);
    else
        if u >= U(n+1), span = n; return; end
        low  = p+1; 
        high = n+1; 
        mid = floor((low+high)/2);
    end
    while (u < U(mid) || u >= U(mid+1))
        if u < U(mid)
            high = mid; 
        else
            low = mid; 
        end
        mid = floor((low+high)/2); 
    end
    span = mid;
end
