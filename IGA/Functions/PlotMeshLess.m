function PlotMeshLess(nodes, hole_idx, L, R, num_nodes, dmax, T)

% --- FIG 1: NODES ---
figure('Name', 'Mesh', 'Color', 'w');
hold on;

plot(nodes(:,1), nodes(:,2), 'b.', 'MarkerSize', 5); 
plot(nodes(hole_idx,1), nodes(hole_idx,2), 'r.', 'MarkerSize', 8);

box on;
axis equal; 
axis([0 L 0 L]); 
title(['Node Distribution (N=' num2str(num_nodes) ')']);
xlabel('X'); ylabel('Y');

% --- FIG 2: SMOOTH FIELD ---
fprintf('Generating Temperature Field...\n');
figure('Name', 'Temperature Field', 'Color', 'w');

% Grid Resolution
N_grid = 120;
gx = linspace(0, L, N_grid);
gy = linspace(0, L, N_grid);
[Xq, Yq] = meshgrid(gx, gy);
Tq = nan(size(Xq));

% DIRECT LOOP: Compute T at every pixel
for row = 1:N_grid
    for col = 1:N_grid
        xp = Xq(row, col);
        yp = Yq(row, col);

        % Skip hole region
        if sqrt((xp-L/2)^2 + (yp-L/2)^2) < R
            continue;
        end

        [phi, ~, neighs] = MLSShapeFunc(xp, yp, nodes, dmax);
        if ~isempty(neighs)
            val = T(neighs);
            Tq(row, col) = sum(phi(:) .* val(:));
        end
    end
end

% Contour Plot
contourf(Xq, Yq, Tq, 80, 'LineColor', 'none');
hold on;

% Draw hole
theta_c = linspace(0, 2*pi, 200);
xc_c = L/2 + R * cos(theta_c);
yc_c = L/2 + R * sin(theta_c);
fill(xc_c, yc_c, 'w', 'EdgeColor', 'none');

colormap jet;
c = colorbar; 
c.Label.String = 'Temperature (K)';
caxis([273 773]);

axis equal; 
axis([0 L 0 L]); 
box on;
title('EFG Temperature Field');
xlabel('X'); 
ylabel('Y');
set(gca, 'FontSize', 12);

end
