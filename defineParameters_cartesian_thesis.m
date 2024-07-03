% thesis version

global northLat Toffset Soffset
global kappax kappay kappaz kappaz_enh
global aH aV yHbar Htc deltax deltay deltaz yHvariation
global alph bet
global Tstar Sstar

% set diffusion parameters:
kappax = 5*10^-3; % diffusion in x direction
kappay = 5*10^-3; % diffusion in y direction
kappaz = 5*10^-4; % diffusion in z direction
kappaz_enh = 0.1; % enhanced diffusion in z direction @ density inversion events


northLat = (53-65)/(10-65);
kappaz_north = 1;

%% double-gyre surface flow parameters:
aH = 1; %strength scaling

%considering north atlantic from equator 0degN to iceland 70degN, the gulf
%stream border is at approx 40degN, i.e. 3/7 in our units:
yHbar = 3/7; % separatrix location
%the variation size does not correspond to the oscillation amplitude size,
%which is approx 2.5 degrees:
yHvariation = 0.02; % separatrix north-south oscillation extent

Htc = 0.1; %scales penetration depth of double gyre. 1 corresponds to 4km, 0.1 to 400 meters
deltax = 0.01; %scales westward intensification length


%% overturning vertical flow parameters: 
aV = 10^-2; %strength scaling

deltay = 0.1; %parameter for latitude at which overturning switches from upwelling to downwelling, where 0 corresponds to 0deg, 1 to 70degN
deltaz = 0.1; %parameter for depth at which overturning switches from northbound to southbound.


%% scale of variations in temperature, salinity around 15degC, 35psu:
%Tstar = 12;
Tstar = 20;
%Sstar = 2;
Sstar = 10;

% linear coefficients of T and S in setting the density, for density
% equation of state linearized around Toffset, Soffset
alph = 2.1e-4;
bet = 7.5e-4;
% Toffset, Soffset:
Toffset = 15.8;
Soffset = 35.3;

%% Grid parameters
alph = 2.2e-4;

%%%%%%%%%%%%%%%%%%%%%%
% The chosen resolution:
x1 = 0:0.005:0.1;
x2 = 0.1:0.02:1;
x = [x1,x2];
x_cell = unique(x);
y_cell = 0:0.02:1;
z1 = 0:0.05:0.7;
z2 = 0.7:0.02:0.98;
z3 = 0.98:0.01:1;
z = [z1,z2,z3];
z_cell = unique(z);
%%%%%%%%%%%%%%%%%%%%%

x_size = diff(x_cell);
y_size = diff(y_cell);
z_size = diff(z_cell);

%then the evaluation points are at the middle of each cell:
x_eval = (x_cell(1:end-1) + x_cell(2:end))/2;
y_eval = (y_cell(1:end-1) + y_cell(2:end))/2;
z_eval = (z_cell(1:end-1) + z_cell(2:end))/2;

%identify y-index at which to raise vertical diffusion kappaz for density inversions:
[~, northLat_ind] =  min(abs(y_eval - northLat));

Lx = length(x_eval);
Ly = length(y_eval);
Lz = length(z_eval);

lon_eval = (1-x_eval) * 60; %0degW to 60degW
lat_eval = (1-y_eval) * 70; %0degN(Eq) to 70degN

[Y_eval,Z_eval,X_eval] = meshgrid(y_eval,z_eval, x_eval);
[Y_size,Z_size,X_size] = meshgrid(y_size,z_size, x_size);
boxvolume = X_size .* Y_size .* Z_size;

X_size2 = diff(X_eval,1,3);
Y_size2 = diff(Y_eval,1,2);
Z_size2 = diff(Z_eval,1,1);

[Y_cell,Z_cell,X_cell] = meshgrid(y_cell, z_cell, x_cell);

dt = 0.0005;

lon_eval = (1-x_eval) * 60;
lat_eval = (1-y_eval) * (65-10) + 10;
depth_eval = (1-z_eval)*4000;

