%% ========================================================================
%  REGRESSION EQUATION OF TEMPERATURE FIELD EQUATION | IGA DATABASE
%  2D Steady-State Anisotropic Heat Conduction: Plate with Circular Hole
%
%  BCs:  T = Th  on inner circle r = R  (centred at L/2, L/2)
%        T = To  on outer square [0,L] x [0,L]
%
%  EQUATION FORM:
%    T(x,y) = To + (Th - To) * C(L, kx, ky, R, x, y)
%
%  where f is a polynomial in 6 or 7 physically-derived variables (depend on terms):
%    p                                     radial position
%    q                                     domain size
%    la                                    anisotropy ratio
%    RL                                    normalised source radius
%    c2                                    n = 2 angular harmonic
%    c4                                    n = 4 angular harmonic
%    c6                                    n = 6 angular harmonic
% =========================================================================

clear, clc, close all

%% =========================================================================
%  INPUT
% =========================================================================
L   = 50;
R   = 4;
kx  = 1003;
ky  = 171;
Th  = 773;
To  = 273;

x_pt  = 25;   % Test point x
y_pt  = 50;   % Test point y  
Ngrid = 0;    % grid for plots 

% ===================================
%  EVALUATE TEMPERATURE AT TEST POINT ------------------------------------
% ===================================
[T_pred, T_norm_pred, vars] = ML_temp_prediction(x_pt, y_pt, L, kx, ky, R, To, Th);

fprintf('=================================================================\n')
fprintf(' TEMPERATURE EQUATION RESULT\n')
fprintf('=================================================================\n')
fprintf('  Inputs:\n')
fprintf('    L = %g  kx = %g  ky = %g  R = %g\n',L,kx,ky,R)
fprintf('    x = %.4f  y = %.4f\n',x_pt,y_pt)
fprintf('    Th = %g K   To = %g K\n\n',Th,To)
fprintf('  Intermediate variables:\n')
fprintf('    xi  = (x-L/2)/sqrt(kx)    = %10.6f\n', vars.xi)
fprintf('    eta = (y-L/2)/sqrt(ky)    = %10.6f\n', vars.eta)
fprintf('    phi (angle in xi-eta)     = %10.6f rad (%.4f deg)\n', vars.phi, rad2deg(vars.phi))
fprintf('    rho_inner(phi)            = %10.6f\n', vars.rho_inner)
fprintf('    rho_outer(phi)            = %10.6f\n', vars.rho_outer)
fprintf('    p = ln(rho/rho_inner)     = %10.6f\n', vars.p)
fprintf('    q = ln(rho_outer/inner)   = %10.6f\n', vars.q)
fprintf('    lambda = ln(sqrt(kx/ky))  = %10.6f\n', vars.la)
fprintf('    RL = R/L                  = %10.6f\n', vars.RL)
fprintf('    cos(2*phi)                = %10.6f\n', vars.c2)
fprintf('    cos(4*phi)                = %10.6f\n', vars.c4)
fprintf('    cos(6*phi)                = %10.6f\n\n', vars.c6)
fprintf('  Result:\n')
fprintf('    T      = %.4f K  (%.4f C)\n\n', T_pred, T_pred-273.000)
fprintf('  Equation accuracy (from IGA training):\n')
fprintf('    R² = 0.99961   Mean |T| = 1.62 K\n')
fprintf('=================================================================\n')

if Ngrid ~= 0
    %% ===========================
    %  GENERATE 2D TEMPERATURE MAP -------------------------------------------
    % ============================
    fprintf('\nComputing 2D field...\n')
    xvec = linspace(0,L,Ngrid);
    yvec = linspace(0,L,Ngrid);
    [XX,YY] = meshgrid(xvec,yvec);
    
    for ii=1:Ngrid
        for jj=1:Ngrid
            xij = XX(ii,jj); yij = YY(ii,jj);
            rc = sqrt((xij-L/2)^2 + (yij-L/2)^2);
            if rc <= R
                TT(ii,jj) = Th;
            elseif xij<=0||xij>=L||yij<=0||yij>=L
                TT(ii,jj) = To;
            else
                TT(ii,jj) = ML_Temp(xij, yij, L, kx, ky, R, To, Th, sim_num);
            end
        end
    end
    
    % ======
    %  PLOT ----------------------------------------------------------------
    % ======
    
    figure('Color','w','Units','centimeters','Position',[5 5 11.5 9.6]);
    hold on;
    
    % ---- Mask hole first ----
    rc = sqrt((XX - L/2).^2 + (YY - L/2).^2);
    TT(rc <= R) = NaN;
    
    % ---- Field ----
    contourf(XX, YY, TT, Ngrid, 'LineColor','none');
    colormap(jet);
    shading interp;
    
    % ---- IMPORTANT: lock geometry BEFORE drawing circle ----
    axis equal tight;
    view(2);
    
    % ---- Clean circle (TRUE geometry) ----
    th_c = linspace(0,2*pi,Ngrid*2);   % more points = smoother circle
    xhole = L/2 + R*cos(th_c);
    yhole = L/2 + R*sin(th_c);
    
    plot(xhole, yhole, 'w', 'LineWidth', 1);
    
    % ---- Formatting ----
    ax = gca;
    set(ax, 'FontName','cambria', ...
        'FontSize',12, ...
        'TickDir','out', ...
        'Box','off');
    
    xlabel('x');
    ylabel('y');
    title('PIMLC 150 Temperature Field T(x,y)', 'FontWeight','bold');
    
    % ---- Colorbar ----
    Tmin = floor(min(TT(:),[],'omitnan'));
    Tmax = ceil(max(TT(:),[],'omitnan'));
    caxis([Tmin Tmax]);
    
    % ---- Colorbar ----
    cb = colorbar;
    set(cb, 'FontName','cambria', 'FontSize',12);
    
    % ---- EXACT 8 INTEGER TICKS ----
    ticks = linspace(Tmin, Tmax, 8);
    ticks = round(ticks);
    ticks = unique(ticks);
    
    % If duplicates reduce number, rebuild properly
    if numel(ticks) < 8
        ticks = Tmin : (Tmax - Tmin)/7 : Tmax;
        ticks = round(ticks);
        ticks = unique(ticks);
    end
    
    cb.Ticks = ticks;
    set(cb, 'FontName','cambria', 'FontSize',12);
    set(gcf,'Renderer','painters')
end