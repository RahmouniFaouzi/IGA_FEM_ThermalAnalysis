function [R, dRdxi, dRdeta] = NURBS2DBasisDers(uv, p, q, U, V, weights)
% NURBS2DBASISDERS Calculates 2D NURBS Basis functions and Derivatives
%
%   Mathematically, a NURBS basis function R is defined as:
%       R(u,v) = ( N(u)*M(v)*w ) / sum( N(u)*M(v)*w )
%
%   Inputs:
%       uv      : [u, v] parametric coordinates
%       p, q    : Polynomial orders (u-dir, v-dir)
%       U, V    : Knot vectors
%       weights : Vector of weights for the active element
    
    u = uv(1); v = uv(2);
    
    % 1. FIND KNOT SPANS
    % NURBS are local. We first find which "element" (knot span) the point lies in.
    nU = length(U) - p - 2; 
    nV = length(V) - q - 2;
    spanU = FindSpan_(nU, p, u, U); 
    spanV = FindSpan_(nV, q, v, V);
    
    % 2. COMPUTE 1D B-SPLINE BASIS (Non-Rational)
    % We compute the univariate basis functions N(u) and M(v) separately.
    % 'dersU' contains [N values; dN/du values]
    dersU = BasisFunsDers(spanU, u, p, U, 1); 
    dersV = BasisFunsDers(spanV, v, q, V, 1);
    
    Nu = dersU(1,:); dNu = dersU(2,:); 
    Nv = dersV(1,:); dNv = dersV(2,:);
    
    % 3. TENSOR PRODUCT GENERATION
    % A 2D surface basis is the product of two 1D curve bases: N_2D = N(u) * M(v)
    num_basis = (p+1)*(q+1); 
    N = zeros(num_basis, 1); 
    dNdxi = zeros(num_basis, 1); 
    dNdeta = zeros(num_basis, 1);
    
    k = 1;
    % Loop order must match the control point ordering (usually v then u)
    for jj = 1:q+1 
        for ii = 1:p+1 
             % The tensor product rule:
             N(k)      = Nu(ii) * Nv(jj); 
             
             % Product rule for derivatives:
             % d/dxi depends only on u-terms, d/deta depends only on v-terms
             dNdxi(k)  = dNu(ii) * Nv(jj); 
             dNdeta(k) = Nu(ii) * dNv(jj); 
             k = k + 1;
        end 
    end
    
    % 4. APPLY WEIGHTS (The "Rational" part of NURBS)
    % Compute the Weight Function W(u,v) = sum( N_i * w_i )
    try
        W_sum = dot(N, weights);
    catch ME
        % code that runs if an error occurs
        disp('An error occurred:');
        disp(ME.message);
    end

    % Compute derivatives of the Weight Function dW/dxi and dW/deta
    dW_dxi = dot(dNdxi, weights); 
    dW_deta = dot(dNdeta, weights);
    
    % 5. COMPUTE RATIONAL BASIS AND DERIVATIVES
    % We must divide by W_sum. To differentiate (N*w)/W, we use the Quotient Rule:
    % (f/g)' = (f'g - fg') / g^2
    
    R = zeros(num_basis, 1); 
    dRdxi = zeros(num_basis, 1); 
    dRdeta = zeros(num_basis, 1);
    
    invW = 1/W_sum; % Pre-compute 1/W for speed
    
    for k = 1:num_basis
        % The Shape Function: R = (N*w) / W
        fac = weights(k) * invW; 
        R(k) = N(k) * fac;
        
        % The Derivative (Quotient Rule optimized):
        % dR/dxi = ( w*dN/dxi - R * dW/dxi ) / W
        dRdxi(k)  = (weights(k) * dNdxi(k)) * invW - R(k) * dW_dxi * invW;
        dRdeta(k) = (weights(k) * dNdeta(k)) * invW - R(k) * dW_deta * invW;
    end
end

function ders = BasisFunsDers(i, u, p, U, n_ders)
% BASISFUNSDERS Cox-de Boor Recursion Formula
%   Computes non-rational B-Spline basis functions and derivatives.
%   
%
%   It builds a triangle of values starting from p=0 (box functions)
%   up to order p, ensuring numerical stability.

    ders = zeros(n_ders+1, p+1); 
    ndu = zeros(p+1, p+1); 
    ndu(1,1)=1;
    
    left = zeros(p+1,1); right = zeros(p+1,1);
    
    % --- Step 1: Compute Basis Functions (Triangle Up) ---
    for j=1:p
        left(j+1) = u - U(i+1-j); 
        right(j+1) = U(i+j) - u; 
        saved = 0;
        
        % Recursive update based on previous order
        for r=0:j-1
            % Compute terms for the triangular recurrence
            ndu(j+1,r+1) = right(r+2) + left(j-r+1); 
            temp = ndu(r+1,j) / ndu(j+1,r+1);
            
            ndu(r+1,j+1) = saved + right(r+2)*temp; 
            saved = left(j-r+1)*temp;
        end
        ndu(j+1,j+1) = saved;
    end
    
    % Store function values (0-th derivative)
    ders(1,:) = ndu(:,p+1)'; 
    
    % --- Step 2: Compute Derivatives (Triangle Down) ---
    % Derivatives of B-splines are calculated using linear combinations
    % of B-splines of lower order (p-1).
    a = zeros(2,p+1);
    for r=0:p
        s1=0; s2=1; a(1,1)=1;
        for k=1:n_ders
            d=0; rk=r-k; pk=p-k;
            
            % Compute coefficients based on knot distances
            if(r>=k)
                a(s2+1,1) = a(s1+1,1) / ndu(pk+2,rk+1); 
                d = a(s2+1,1) * ndu(rk+1,pk+1); 
            end
            
            if(rk>=-1), j1=1; else, j1=-rk; end
            if(r-1<=pk), j2=k-1; else, j2=p-r; end
            
            % Loop to accumulate derivative terms
            for j=j1:j2
                a(s2+1,j+1) = (a(s1+1,j+1)-a(s1+1,j)) / ndu(pk+2,rk+j+1); 
                d = d + a(s2+1,j+1) * ndu(rk+j+1,pk+1); 
            end
            
            if(r<=pk)
                a(s2+1,k+1) = -a(s1+1,k) / ndu(pk+2,r+1); 
                d = d + a(s2+1,k+1) * ndu(r+1,pk+1); 
            end
            
            ders(k+1,r+1) = d; 
            j=s1; s1=s2; s2=j;
        end
    end
    
    % Multiply by p for the correct derivative scaling factor
    r=p; 
    for k=1:n_ders
        for j=0:p
            ders(k+1,j+1) = ders(k+1,j+1) * r; 
        end
        r = r * (p-k); 
    end
end