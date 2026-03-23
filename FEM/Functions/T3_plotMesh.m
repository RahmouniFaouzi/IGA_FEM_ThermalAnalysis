function T3_plotMesh(node, elem, a, b, A, B, C, D, opts)
% ------------------------------------------------------------
% Plot T3 FEM mesh with full user control %
% opts.showMesh = true/false
% opts.showNodes = true/false
% opts.showABCD = true/false
% opts.showBoundary = true/false %
% Example:
% opts = struct('showMesh',1,'showNodes',1,'showABCD',1,'showBoundary',1);
% T3_plotMesh(node,elem,a,b,opts)
% ------------------------------------------------------------
figure; hold on; box on; axis equal;
%% ---- Plot mesh ----

if opts.showMesh
    for e = 1:size(elem,1)
        n = elem(e,:);
        fill(node(n,1), node(n,2), 'w', ...
            'EdgeColor','k','LineWidth',0.5);
    end
end

%% ---- Plot node numbers ----
if opts.showNodes
    for i = 1:size(node,1)
        text(node(i,1), node(i,2), num2str(i), ...
            'FontSize',8,'Color','r', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle');
    end
end

%% ---- Plot A, B, C, D ----
if opts.showABCD
    plot(A(1),A(2),'ko','MarkerFaceColor','k','MarkerSize',5)
    plot(B(1),B(2),'ko','MarkerFaceColor','k','MarkerSize',5)
    plot(C(1),C(2),'ko','MarkerFaceColor','k','MarkerSize',5)
    plot(D(1),D(2),'ko','MarkerFaceColor','k','MarkerSize',5)
    text(A(1),A(2),' A','FontSize',11,'FontWeight','bold')
    text(B(1),B(2),' B','FontSize',11,'FontWeight','bold')
    text(C(1),C(2),' C','FontSize',11,'FontWeight','bold')
    text(D(1),D(2),' D','FontSize',11,'FontWeight','bold')
end

%% ---- Plot boundaries ----
if opts.showBoundary
    tol = 1e-10;
    bottom = find(abs(node(:,2)) < tol);
    top = find(abs(node(:,2)-b) < tol);
    left = find(abs(node(:,1)) < tol);
    right = find(abs(node(:,1)-a) < tol);
    plot(node(bottom,1),node(bottom,2),'ro','MarkerFaceColor','r')
    plot(node(top,1),node(top,2),'bo','MarkerFaceColor','b')
    plot(node(left,1),node(left,2),'go','MarkerFaceColor','g')
    plot(node(right,1),node(right,2),'mo','MarkerFaceColor','m')
end

%% ---- Labels ----
xlabel('x');
ylabel('y');
title(sprintf('T3 FEM Mesh (a/b = %.2f)',a/b));
xlim([0 a]);
ylim([0 b]);
end
