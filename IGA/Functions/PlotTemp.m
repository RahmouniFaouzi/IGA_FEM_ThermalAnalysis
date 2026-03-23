function PlotTemp(patches, global_indices, T, L, R)
    hold on;
    % Draw the central hole as a white filled circle
    th = linspace(0, 2*pi, 100); 
    fill(L/2+R*cos(th), L/2+R*sin(th), 'w');
    
    for k = 1:4
        P = patches(k); 
        g_map = global_indices{k};
        
        % Define parametric grid for plotting resolution
        nu_grid = 15; nv_grid = 15;
        us = linspace(0, 1, nu_grid); 
        vs = linspace(0, 1, nv_grid);
        [UU, VV] = meshgrid(us, vs);
        
        XX = zeros(size(UU)); 
        YY = zeros(size(UU)); 
        TT = zeros(size(UU));
        
        % Loop over grid points to compute physical coordinates and field values
        for r = 1:size(UU,1)
            for c = 1:size(UU,2)
                u = UU(r,c); 
                v = VV(r,c);
                
                % Evaluate Field (Temperature) at (u,v)
                val = evalTempTrans(P, g_map, T, u, v, 2, 2);
                
                % Compute Geometry Mapping (Parametric u,v -> Physical x,y)
                spanU = FindSpan_(P.n_u, 2, u, P.uKnot, 0); 
                spanV = FindSpan_(P.n_v, 2, v, P.vKnot, 0);
                Nu = BasisFuns(spanU, u, 2, P.uKnot); 
                Nv = BasisFuns(spanV, v, 2, P.vKnot);
                
                x_val=0; 
                y_val=0; 
                sum_w=0;
                
                % Sum over local basis support
                for jj=(spanV-2):spanV
                    for ii=(spanU-2):spanU
                        idx = (jj-1)*P.n_u + ii;
                        N_v = Nu(ii-(spanU-2)+1) * Nv(jj-(spanV-2)+1);
                        cp_w = P.controlPts(idx, 4);
                        
                        % Calculate weighted physical coordinates
                        x_val = x_val + N_v*(P.controlPts(idx,1)/cp_w); 
                        y_val = y_val + N_v*(P.controlPts(idx,2)/cp_w); 
                        sum_w = sum_w + N_v*cp_w;
                    end
                end
                % Project back to physical space
                XX(r,c) = (x_val*sum_w)/sum_w; 
                YY(r,c) = (y_val*sum_w)/sum_w; 
                TT(r,c) = val;
            end
        end
        % Plot the surface
        surf(XX, YY, TT, 'EdgeColor', 'none', 'FaceColor', 'interp');
    end
    % Manually expand the subplot
    pos = get(gca, 'Position');
    inset = 0.02; % Expansion factor
    set(gca, 'Position', [pos(1)-inset, pos(2)-inset, pos(3)+2*inset, pos(4)+2*inset]);
    colorbar;
end