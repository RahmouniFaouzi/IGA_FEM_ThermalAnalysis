function plotGRADEDIGA(patches, T, U, V, p, q)
% plotGRADEDIGA  Plot temperature field and computational mesh for IGA model
% ============================================================
% FIGURE 1: Temperature Field
% ============================================================
figure('Color','w', 'Name', 'Temperature Field');
hold on; axis equal; box on; axis off;

for k = 1:length(patches)
    patch = patches{k};
    
    uu = linspace(0,1,40);
    vv = linspace(0,1,40);
    [UU, VV] = meshgrid(uu, vv);
    
    XX = zeros(size(UU));
    YY = XX;
    TT = XX;
    
    for m = 1:numel(UU)
        u_pt = UU(m);
        v_pt = VV(m);
        
        spanU = FindSpan_(length(U)-p-2, p, u_pt, U);
        spanV = FindSpan_(length(V)-q-2, q, v_pt, V);
        
        idx_u = (spanU-p):spanU;
        idx_v = (spanV-q):spanV;
        
        glob_sctr = patch.g_map(idx_u, idx_v);
        glob_sctr = glob_sctr(:);
        
        [JJ, II] = meshgrid(idx_v, idx_u);
        lin_idx  = sub2ind([patch.nu, patch.nv], II(:), JJ(:));
        
        local_w  = patch.cp_list(lin_idx, 4);
        local_cp = patch.cp_list(lin_idx, 1:2) ./ repmat(local_w, 1, 2);
        
        [R_b, ~] = NURBS2DBasisDers([u_pt; v_pt], p, q, U, V, local_w);
        
        phys   = R_b' * local_cp;
        XX(m)  = phys(1);
        YY(m)  = phys(2);
        TT(m)  = dot(R_b, T(glob_sctr));
    end
    
    surf(XX, YY, zeros(size(XX)), TT, ...
        'EdgeColor', 'none', 'FaceColor', 'interp');
end

colormap jet;
colorbar;
title('Temperature Field');

% ============================================================
% FIGURE 2: Computational Mesh
% ============================================================
figure('Color','w', 'Name', 'Computational Mesh');
hold on; axis equal; box on; axis off;
title('Graded Computational Mesh');

for k = 1:length(patches)
    patch = patches{k};
    uniqueU = unique(U);
    uniqueV = unique(V);
    
    % ---- V lines ----
    for v_val = uniqueV
        u_tr = linspace(0,1,40);
        lx = zeros(size(u_tr));
        ly = lx;
        
        for i = 1:length(u_tr)
            u_pt = u_tr(i);
            v_pt = v_val;
            
            spanU = FindSpan_(length(U)-p-2, p, u_pt, U);
            spanV = FindSpan_(length(V)-q-2, q, v_pt, V);
            
            idx_u = (spanU-p):spanU;
            idx_v = (spanV-q):spanV;
            
            [JJ, II] = meshgrid(idx_v, idx_u);
            lin_idx  = sub2ind([patch.nu, patch.nv], II(:), JJ(:));
            
            local_w  = patch.cp_list(lin_idx, 4);
            local_cp = patch.cp_list(lin_idx, 1:2) ./ repmat(local_w, 1, 2);
            
            [R_b, ~] = NURBS2DBasisDers([u_pt; v_pt], p, q, U, V, local_w);
            
            phys  = R_b' * local_cp;
            lx(i) = phys(1);
            ly(i) = phys(2);
        end
        
        plot(lx, ly, 'k-', 'LineWidth', 0.5);
    end
    
    % ---- U lines ----
    for u_val = uniqueU
        v_tr = linspace(0,1,40);
        lx = zeros(size(v_tr));
        ly = lx;
        
        for i = 1:length(v_tr)
            u_pt = u_val;
            v_pt = v_tr(i);
            
            spanU = FindSpan_(length(U)-p-2, p, u_pt, U);
            spanV = FindSpan_(length(V)-q-2, q, v_pt, V);
            
            idx_u = (spanU-p):spanU;
            idx_v = (spanV-q):spanV;
            
            [JJ, II] = meshgrid(idx_v, idx_u);
            lin_idx  = sub2ind([patch.nu, patch.nv], II(:), JJ(:));
            
            local_w  = patch.cp_list(lin_idx, 4);
            local_cp = patch.cp_list(lin_idx, 1:2) ./ repmat(local_w, 1, 2);
            
            [R_b, ~] = NURBS2DBasisDers([u_pt; v_pt], p, q, U, V, local_w);
            
            phys  = R_b' * local_cp;
            lx(i) = phys(1);
            ly(i) = phys(2);
        end
        
        plot(lx, ly, 'k-', 'LineWidth', 1);
    end
end
end
