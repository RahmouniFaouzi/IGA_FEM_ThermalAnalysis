function patch = convert2DNurbsToPatch(nurb)
% CONVERT2DNURBSTOPATCH Converts NURBS Toolbox object to Analysis Struct
%   Input:
%       nurb : Struct from nrbmak/nrbrefine (coefs is 4xNuxNv)
%   Output:
%       patch : Analysis-ready struct with .element, .controlPts, etc.

% 1. EXTRACT BASIC PROPERTIES
patch.p = nurb.order(1)-1;
patch.q = nurb.order(2)-1;

patch.uKnot = nurb.knots{1};
patch.vKnot = nurb.knots{2};

% 2. FLATTEN CONTROL POINTS
cp = nurb.coefs;
[~, nu, nv] = size(cp);
patch.noPtsX = nu;
patch.noPtsY = nv;

flat_cp = zeros(nu*nv, 4);
cnt = 1;
% We loop u (col) then v (row) to flatten
for i=1:nu
    for j=1:nv
        % Transpose to get row vector [x, y, z, w]
        flat_cp(cnt,:) = cp(:,i,j)';
        cnt=cnt+1;
    end
end
patch.controlPts = flat_cp;

% 3. GENERATE 1D CONNECTIVITY
[patch.elRangeU, elConnU] = Connectivity(patch.p, patch.uKnot);
[patch.elRangeV, elConnV] = Connectivity(patch.q, patch.vKnot);

patch.noElemsU = size(patch.elRangeU, 1);
patch.noElemsV = size(patch.elRangeV, 1);
noElems = patch.noElemsU * patch.noElemsV;

% 4. CREATE LOCAL GRID MAP
local_pattern = zeros(nv, nu);
for col = 1:nu
    for row = 1:nv
        % Column-major ordering formula (Standard MATLAB)
        local_pattern(row, col) = row + (col-1)*nv;
    end
end

% 5. BUILD 2D ELEMENT TABLE (Tensor Product)
patch.element = zeros(noElems, (patch.p+1)*(patch.q+1));
patch.index   = zeros(noElems, 2);

e = 1;
% Loop over every non-zero knot span grid (Elements)
for v = 1:patch.noElemsV
    vConn = elConnV(v, :); % Indices of basis functions in V direction
    
    for u = 1:patch.noElemsU
        uConn = elConnU(u, :); % Indices of basis functions in U direction
        c = 1;
        for i = 1:length(vConn)
            for j = 1:length(uConn)
                % Lookup the linear node ID for this combination
                patch.element(e, c) = local_pattern(vConn(i), uConn(j));
                c=c+1;
            end
        end
        
        % Store which knot span (index) this element belongs to
        patch.index(e, :) = [u, v];
        e = e + 1;
    end
end
end