function plot_results(pts_eval, elem, node, L, T)
%% ============================================================
%  3. PLOTTING (Mesh + Points + Field)
% ============================================================
figure('Color','w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.5]);

% Subplot 1: Mesh & Points
subplot(1,2,1); hold on;
triplot(elem, node(:,1), node(:,2), 'Color', [0.5 0.5 0.5]); % Grey Mesh
plot(pts_eval(:,1), pts_eval(:,2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8); % Red Points
labels = {'A', 'B', 'C', 'D'};
% Label Points
for i=1:4
    text(pts_eval(i,1)+1, pts_eval(i,2), labels{i}, 'FontWeight','bold', 'FontSize',12);
end
axis equal; 
axis([0 L 0 L]);
title('Mesh & Evaluation Points');
xlabel('X'); 
ylabel('Y');

% Subplot 2: Temperature Field
subplot(1,2,2);
trisurf(elem, node(:,1), node(:,2), T, 'EdgeColor', 'none');
view(2); 
shading interp; 
colormap jet; 
colorbar;
axis equal; 
axis([0 L 0 L]);
title('Temperature Field');
xlabel('X'); 
ylabel('Y');
end