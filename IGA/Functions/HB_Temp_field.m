function HB_Temp_field(pts, a, b, uKnot, vKnot, ncp1_x, ne_x, ne_y, dx, dy, p, q, T, RefinedFlags)
% -- 1. Temperature Contour Plot --
figure('Color','w', 'Name', 'Temperature Field');
hold on;
axis equal;
axis([0 a 0 b]);
xlabel('x'); ylabel('y');

res = 200; % Resolution for smooth plotting
Xp  = linspace(0, a, res);
Yp  = linspace(0, b, res);
[Xm,Ym] = meshgrid(Xp, Yp);
Tm = zeros(size(Xm));

fprintf('Generating Temperature Field (%dx%d)... \n', res, res);

% Interpolate solution onto plotting grid
for k = 1:numel(Xm)
    x0 = Xm(k); y0 = Ym(k);
    
    % Normalize coords for basis evaluation
    u_param = x0/a;
    v_param = y0/b;
    
    sp_u = find(uKnot <= u_param, 1, 'last');
    if (u_param >= 1), sp_u = length(uKnot) - p - 1; end
    
    sp_v = find(vKnot <= v_param, 1, 'last');
    if (v_param >= 1), sp_v = length(vKnot) - q - 1; end
    
    [Nu,~] = HB_BasisFuns(sp_u, u_param, p, uKnot);
    [Nv,~] = HB_BasisFuns(sp_v, v_param, q, vKnot);
    
    val = 0;
    idx_u = (sp_u-p):sp_u; idx_v = (sp_v-q):sp_v;
    
    for j = 1:length(idx_v)
        for i = 1:length(idx_u)
            dof = (idx_v(j) - 1) * ncp1_x + idx_u(i);
            val = val + Nu(i) * Nv(j) * T(dof);
        end
    end
    Tm(k) = val;
end

contourf(Xm,Ym,Tm, 40, 'LineStyle','none');
colormap jet;
c = colorbar; c.Label.String = 'Temperature [K]';
title(['HB-Spline Temperature Distribution (p=', num2str(p), ', q=', num2str(q), ')']);
hold off;

% -- 2. Mesh Topology Plot --
figure('Color','w', 'Name', 'Mesh Topology');
axis equal; 
axis([0 a 0 b]);
xlabel('x'); 
ylabel('y');
hold on;

for i = 1:ne_x
    for j = 1:ne_y
        x1 = (i-1) * dx;
        y1 = (j-1) * dy;
        
        if (RefinedFlags(i,j) == 0)
            % Draw Coarse Element
            rectangle('Position',[x1,y1,dx,dy], 'EdgeColor','k', 'LineWidth', 1.0);
        else
            % Draw 4 Fine Elements (Refined Block)
            rectangle('Position',[x1,      y1,      dx/2, dy/2], 'EdgeColor','b');
            rectangle('Position',[x1+dx/2, y1,      dx/2, dy/2], 'EdgeColor','b');
            rectangle('Position',[x1,      y1+dy/2, dx/2, dy/2], 'EdgeColor','b');
            rectangle('Position',[x1+dx/2, y1+dy/2, dx/2, dy/2], 'EdgeColor','b');
        end
    end
end

% -- 3. Plot Probe Points (ADDED THIS) --
% Plot red filled circles
plot(pts(:,1), pts(:,2), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

% -- Add text labels (A, B, C, D) --
lbls = {'A','B','C','D'};
for k = 1:min(4, size(pts,1))
    text(pts(k,1)+0.02*a, pts(k,2), lbls{k}, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
end

title('Hierarchical Mesh Topology (Black=Coarse, Blue=Refined)');
hold off;
fprintf('Done.\n');
end