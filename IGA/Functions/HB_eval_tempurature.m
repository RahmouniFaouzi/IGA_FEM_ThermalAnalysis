function T_ev = HB_eval_tempurature(pts, uKnot, vKnot, p, q, a, b, ncp1_x, T)
% Evaluate tempurature in pts point
x0 = pts(1);
y0 = pts(2);

% Find Knot Spans
sp_u = find(uKnot <= x0/a, 1, 'last');
if (x0/a == 1), sp_u = length(uKnot)-p-1; end

sp_v = find(vKnot <= y0/b, 1, 'last');
if (y0/b == 1), sp_v = length(vKnot)-q-1; end

% Evaluate Basis Functions
[Nu,~] = HB_BasisFuns(sp_u, x0/a, p, uKnot);
[Nv,~] = HB_BasisFuns(sp_v, y0/b, q, vKnot);

% Compute Solution at Point
idx_u = (sp_u-p):sp_u;
idx_v = (sp_v-q):sp_v;
T_ev = 0;

for j = 1:length(idx_v)
    for i = 1:length(idx_u)
        dof = (idx_v(j)-1) * ncp1_x + idx_u(i);
        T_ev = T_ev + Nu(i)*Nv(j) * T(dof);
    end
end
end