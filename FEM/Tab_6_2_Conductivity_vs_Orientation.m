clc, clear, close all
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'Functions')))

%% Fiber and matrix properties
kf = 16;       % Fiber conductivity (W/m-K)
km = 171;      % Matrix conductivity (W/m-K)
vf = 0.5;      % Fiber volume fraction

%% On-axis ply conductivities
k1 = kf*vf + km*(1-vf);       % along fiber
k2 = km + ( (vf*(kf-km)*km)/ (0.5*(1-vf)*(kf-km)+km) );  % transverse

%% Ply orientations (degrees)
theta_deg = [30 45 60 90];

%% Storage
kx = zeros(size(theta_deg));
ky = zeros(size(theta_deg));

%% Compute [0/theta]s symmetric laminate
for i = 1:length(theta_deg)
    theta = theta_deg(i)*pi/180;  % convert to radians
    
    % Conductivity of theta ply in global x-y
    kx_theta = k1*cos(theta)^2 + k2*sin(theta)^2;
    ky_theta = k1*sin(theta)^2 + k2*cos(theta)^2;
    
    % Conductivity of 0-degree ply
    kx_0 = k1;
    ky_0 = k2;
    
    % Symmetric laminate [0/theta]s: average
    kx(i) = 0.5*(kx_0 + kx_theta);
    ky(i) = 0.5*(ky_0 + ky_theta);
end

%% Print Table 6
fprintf('\nTable 6. Thermal conductivity for symmetric laminates\n');
fprintf('-----------------------------------------------------\n');
fprintf('Ply orientation     kx (W/m-K)     ky (W/m-K)\n');
fprintf('-----------------------------------------------------\n');

labels = {'[0/30]s','[0/45]s','[0/60]s','[0/90]s'};
for i = 1:length(labels)
    fprintf('%8s            %7.2f          %7.2f\n', labels{i}, kx(i), ky(i));
end
