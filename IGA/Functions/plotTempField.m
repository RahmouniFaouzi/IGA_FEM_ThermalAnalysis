function plotTempField(Surf, T, a, b, DIS)
% plotTempField: Generates a high-res thermal contour plot (MATLAB 2015)
% Inputs:
%   Surf - NURBS Geometry Structure
%   T    - Solved temperature vector (control point DOFs)
%   a,b  - domain width and height
%   DIS  - optional verbosity flag (1 = print progress)

if nargin < 5, DIS = 1; end
if DIS, fprintf('Generating Temperature Field and Plot...\n'); end

% 1. SETUP GRID 
res = 200;                        % resolution of the plotted grid
xx  = linspace(0, a, res);
yy  = linspace(0, b, res);
[XX, YY] = meshgrid(xx, yy);

T_field = zeros(res, res);

% Evaluate temperature at each grid point using the provided evaluator
for i = 1:res
    for j = 1:res
        pt = [XX(i,j), YY(i,j)];
        T_field(i,j) = eval_IGA_Temperature(pt, T, Surf);
    end
end

if DIS, fprintf('Field evaluated (%d x %d). Preparing figure...\n', res, res); end

% ----------------------------
% 2. VISUALIZATION
% ----------------------------
figure('Name','NURBS-IGA Temperature Field','Color','w','Units','normalized');
% filled contour - hide contour lines for a smooth field
contourf(XX, YY, T_field, 100, 'LineStyle', 'none'); 
colormap(jet(256));

% Use actual field range for color scaling
min_T = min(T_field(:));
max_T = max(T_field(:));
if min_T == max_T
    caxis([min_T-1e-6, max_T+1e-6]);
else
    caxis([min_T, max_T]);
end

% Improve rendering
shading interp;                % smooth color transitions
set(gcf,'Renderer','opengl');  % better for smooth shading

% Colorbar
cb = colorbar('Location','eastoutside');
set(cb,'FontSize',16,'LineWidth',1);

% Select number of ticks
numTicks = 5;
ticks = linspace(min_T, max_T, max(2,round(numTicks)));
set(cb,'Ticks',ticks);

% Choose label format depending on the data range
rangeT = abs(max_T - min_T);
if rangeT >= 1e4
    fmt = '%.0e';        % use scientific if huge range
elseif rangeT >= 100
    fmt = '%.0f';
elseif rangeT >= 1
    fmt = '%.2f';
else
    fmt = '%.3f';
end

% Create tick labels (MATLAB 2015 safe)
tickLabels = cell(length(ticks),1);
for k = 1:length(ticks)
    tickLabels{k} = sprintf(fmt, ticks(k));
end
set(cb,'TickLabels',tickLabels);

% Axis labels and formatting
% --------------------------
% xlabel('x (mm)','FontSize',14);
% ylabel('y (mm)','FontSize',14);
% title('Temperature Field','FontSize',16,'Interpreter','none');

axis equal;
axis([0 a 0 b]);
set(gca,'FontSize',14,'LineWidth',1.2);
box on;

drawnow;

end
