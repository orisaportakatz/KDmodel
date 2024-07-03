function [aVsave, T, S] = modelRun_realistic_noMOW_SELECTdiff_MIXED_setGamma(Time,saveyesno,name, T, S, Gamma)
% run the KD model with selective vertical diffusion that is enhanced where
% there are density inversions, and original (no MOW) velocity advection
% field - an overturning component (uV=0, vV, wV) and a double-gyre surface
% component (uE, vE, wE=0)

% Input parameters:
% Time - length of simulation, in years
% saveyesno - set to 1 to save the final T and S fields under the name
% "name"
% T, S - initial condition temperature and salinity matrices, make sure
% their resolution fits simulation resolution.
% Gamma - scaling of overturning flow strength with density difference
% between north and south. Typical value - Gamma = 0.065

% General approach: After setting boundary conditions, velocity fields,
% parameters, run with a given Gamma until reaching steady state. Check
% overturning strength. Adjust Gamma and repeat to reach desired
% overturning strength.

%% prepare data
% load S and T from a previous run if you want to sgtart with a different initial condition 
%load('interim/modelRun_param19_eff_smooth_MIXED_IC19long_200yrs.mat', 'S', 'T')
%load('IC_cosine_tem_salt.mat', 'tem_cosine', 'salt_cosine', 'lataxis_model')

%% parameters:
defineParameters_cartesian_thesis

Gamma = 0.065; % default value for Gamma, comment this to set your own Gamma

%% grid:
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

%% Time parameters
dt = 0.0005;
saveevery = 100*dt; %can use this to save the T and S matrices every saveevery times

numberofsteps = length(dt:dt:Time); %calculate number of steps.
numberofsavedsteps = floor(Time/saveevery);
numberofstepsinayear = length(dt:dt:1);

%% boundary conditions - restoring:

%% Option 1: Use a simple cosine for temperature and salinity:
load('IC_cosine_tem_salt.mat', 'tem_cosine', 'salt_cosine', 'lataxis_model')

%evaluate salinity BC at y_eval:
S_atm_month = zeros(12, Ly);
T_atm_month = zeros(12, Ly);
for i = 1:12
    S_atm_month(i,:) = interp1(lataxis_model, fliplr(salt_cosine(i,:)), y_eval);
    T_atm_month(i,:) = interp1(lataxis_model, fliplr(tem_cosine(i,:)), y_eval);
end

% can plot this to see what is being used:

%figure; plot(T_atm_month, y_eval, 'x-')
%set(gca, 'Ydir', 'reverse')

%figure; plot(S_atm_month, y_eval, 'x-')
%set(gca, 'Ydir', 'reverse')

S_atm = zeros(numberofstepsinayear, Ly);
T_atm = zeros(numberofstepsinayear, Ly);

for i = 1:Ly
    temp1 = interp1(linspace(0,1,13), [S_atm_month(:,i); S_atm_month(1,i)], 0:dt:1);
    S_atm(1:numberofstepsinayear, i) = temp1(1:end-1)';
    temp2 = interp1(linspace(0,1,13), [T_atm_month(:,i); T_atm_month(1,i)], 0:dt:1);
    T_atm(1:numberofstepsinayear, i) = temp2(1:end-1)';
end

% move boundary conditions around zero, add the offset values at the end
% for visualization
T_atm = T_atm - Toffset;
S_atm = S_atm - Soffset;

%% Option 2: Use boundary conditions derived from restoring simulation with realistic SODA data:
load('SForcing_param19.mat', 'SForcing_IC19long_param19');
SForcing = SForcing_IC19long_param19;

load('SSTtypes.mat', 'SSTtypes');
T_atm_bare = SSTtypes.SSTzonal_5th_smooth;
T_atm_bare = T_atm_bare - 15.5;
X = SSTtypes.x;

%evaluate salinity BC at y_eval:
S_atm = zeros(2000,Ly);
T_atm = zeros(2000,Ly);
for i = 1:2000
    S_atm(i,:) = interp1(X, S_atm_bare(:,i), y_eval);
    S_atm(i,:) = flip(S_atm(i,:));
    T_atm(i,:) = interp1(X, T_atm_bare(:,i), y_eval);
    T_atm(i,:) = flip(T_atm(i,:));
end

%% Initial conditions:
S0 = S;
T0 = T;

% if you uncomment the saveevery inner loop in the main loop, the Smat, Tmat is
% the S, T saved matrix

%Smat = zeros([floor(size(S)/2) numberofsavedsteps]); %CAN LOOK AT PARTIAL DATA INSTEAD OF SIZE(S) TO SAVE SPACE
Smat = zeros([size(S) numberofsavedsteps]);

%Tmat = zeros([floor(size(T)/2) numberofsavedsteps]); %CAN LOOK AT PARTIAL DATA INSTEAD OF SIZE(S) TO SAVE SPACE
Smat = zeros([size(S) numberofsavedsteps]);

s = 1; %matrix saving index


%% overturning velocity strength calculations:

% The overturning strength is saved for each timestep in aVsave:
aVsave = zeros(numberofsteps+1,1); %  Multiply by 3.1e3 to obtain Sverdrup units

% set latitude of border between north and south boxes, at which density
% averages are calculated to set the overturning strength:
temp = abs(y_eval - ( yHbar - yHvariation)); 
%temp = abs(y_eval - 0.35);
%temp = abs(y_eval - 0.22);
y_switchy = find(temp == min(temp));
y_switchy = y_switchy(1);
%y_switchy = 17;

% matrix of volume of each infinetesimal box for the finite-volume
% approach:
boxvolume = X_size .* Y_size .* Z_size;

%calculate average density in north and south
density_north = sum(sum(sum((1026*(1-alph * T(:,1:y_switchy,:) + bet * S(:,1:y_switchy,:))) .* boxvolume(:,1:y_switchy,:))))/sum(sum(sum(boxvolume(:,1:y_switchy,:))));
density_south = sum(sum(sum((1026*(1-alph * T(:,y_switchy+1:end,:) + bet * S(:,y_switchy+1:end,:))).*boxvolume(:,y_switchy+1:end,:))))/sum(sum(sum(boxvolume(:,y_switchy+1:end,:))));


% set overturning circulation velocity strength as depends on Gamma and
% density difference between north and south:
s1 = 1;
%aV = Gamma * (density_south - density_north); <- tests 1,2
aV = -Gamma * (density_south - density_north);
aVsave(s1) = aV;
s1 = s1 + 1;

%% create velocity field:
x_eval = (x_cell(1:end-1) + x_cell(2:end))/2;
y_eval = (y_cell(1:end-1) + y_cell(2:end))/2;
z_eval = (z_cell(1:end-1) + z_cell(2:end))/2;

%u is defined on (X_cell, Y_eval, Z_eval), etc. for v and w.
[Yx,Zx,Xx] = meshgrid(y_eval,z_eval, x_cell(2:end-1));
[Yy,Zy,Xy] = meshgrid(y_cell(2:end-1), z_eval,x_eval);
[Yz,Zz,Xz] = meshgrid(y_eval, z_cell(2:end-1), x_eval);

%Period of oscillation:
%Th =  0.15;
Th = 1;

timevec = dt:dt:Th;

% clear uH vH wH

% uH = zeros([size(Xx), length(timevec)]);
% vH = zeros([size(Xy), length(timevec)]);
% wH = zeros([size(Xz), length(timevec)]);
%
% for i = 1:length(timevec)
%     tau = timevec(i);
% yHp = yHvariation * sin(2*pi*tau/Th);
% %yHp = 0;
%
% yH = yHbar + yHp;
%
% %These correspond to Mathematica file
% %"CartesianSimulationVelocityFieldCalculations"
%
uHtau = (-1).*deltax.^(-1).*exp(1).^((-1).*deltax.^(-1).*Xx+Htc.^(-1).*Zx).* ...
    ((-1)+exp(1).^(deltax.^(-1).*Xx)).*(exp(1).^(Htc.^(-1))+exp(1).^( ...
    Htc.^(-1).*Zx)).^(-1).*pi.*((-1)+Xx).*((-1)+lambertw(exp(1).^(1+ ...
    deltax.^(-1)))).^(-2).*lambertw(exp(1).^(1+deltax.^(-1)));

vHtau =deltax.^(-2).*exp(1).^((-1).*deltax.^(-1).*Xy+Htc.^(-1).*Zy).*(exp( ...
    1).^(Htc.^(-1))+exp(1).^(Htc.^(-1).*Zy)).^(-1).*((-1)+deltax.*((-1) ...
    +exp(1).^(deltax.^(-1).*Xy))+Xy).*((-1)+lambertw(exp(1).^(1+ ...
    deltax.^(-1)))).^(-2).*lambertw(exp(1).^(1+deltax.^(-1))).*sin( ...
    pi.*Yy);

wHtau = 0;

% uH(:,:,:,i) = uHtau;
% vH(:,:,:,i) = vHtau;
% wH(:,:,:,i) = wHtau;
% end

uV = 0;

vV = deltay.^(-1).*deltaz.^(-2).*exp(1).^((-1).*deltay.^(-1).*((-1)+( ...
    -2).*deltay+Yy+deltay.*(lambertw(exp(1).^(1+deltay.^(-1)))+ ...
    lambertw(exp(1).^(1+deltaz.^(-1)))))).*((-1)+exp(1).^(deltay.^( ...
    -1).*Yy)).*((-1)+Yy).*((-1).*deltaz.*exp(1).^(deltaz.^(-1))+exp(1) ...
    .^(deltaz.^(-1).*Zy).*(deltaz+Zy)).*((-1)+lambertw(exp(1).^(1+ ...
    deltay.^(-1)))).^(-2).*((-1)+lambertw(exp(1).^(1+deltaz.^(-1)))) ...
    .^(-2);

wV = deltay.^(-2).*deltaz.^(-1).*exp(1).^((-1).*deltay.^(-1).*((-1)+( ...
    -2).*deltay+Yz+deltay.*(lambertw(exp(1).^(1+deltay.^(-1)))+ ...
    lambertw(exp(1).^(1+deltaz.^(-1)))))).*(exp(1).^(deltaz.^(-1))+( ...
    -1).*exp(1).^(deltaz.^(-1).*Zz)).*((-1)+deltay.*((-1)+exp(1).^( ...
    deltay.^(-1).*Yz))+Yz).*Zz.*((-1)+lambertw(exp(1).^(1+deltay.^(-1)) ...
    )).^(-2).*((-1)+lambertw(exp(1).^(1+deltaz.^(-1)))).^(-2);

%% start main loop
tic;

timevec = dt:dt:Time;
Tave = 0;
Save = 0;
for t = dt:dt:Time
    
    %% set velocity field for this timestep:
    %taurel = floor(mod(t,Th)/dt)+1;
    
    yHp = yHvariation * sin(2*pi*t/Th);
    %yHp = 0;
    
    yH = yHbar + yHp;
    
    uH = uHtau .*sin(pi.*(2.*Yx+(-1).*yH));
    vH = vHtau .*sin(pi.*(Yy+(-1).*yH));
    
    u = aH * uH + aV * uV;
    v = aH * vH + aV * vV;
    w = aV * wV;
    
    
    %% calculate which month we're in for seasonal forcing:
    k = floor(mod(t * 2000, 2000));
    if k == 0
        k = 2000;
    end
    
    %% identify vertical density inversions
    rho_invert_indicator = sign(diff(1-alph * T + bet * S)) == 1;
    
    
    %% advance PDE
    
    RHS_S = diffusionEquation_oscBC_param19_mixed_SELECTDiff_mixed...
        (t, S, 1, Lx, Ly, Lz, y_eval, X_eval, Y_eval, Z_eval, X_cell, Y_cell, Z_cell,...
        S_atm(k,:),T_atm(k,:), northLat_ind, rho_invert_indicator);
    
    RHS_S = RHS_S + advectionEquation_multidim_FL77_SE_eff(S, u, v, w, dt, Lx, Ly, Lz, X_size2, Y_size2, Z_size2, X_cell, Y_cell, Z_cell);
    
    S = S + dt * RHS_S; %calculate the salinity field at the next time step.
    
    RHS_T = diffusionEquation_oscBC_param19_mixed_SELECTDiff_mixed...
        (t, T, 2, Lx, Ly, Lz, y_eval, X_eval, Y_eval, Z_eval, X_cell, Y_cell, Z_cell,...
        S_atm(k,:), T_atm(k,:), northLat_ind, rho_invert_indicator);
    
    RHS_T = RHS_T + advectionEquation_multidim_FL77_SE_eff(T, u, v, w, dt, Lx, Ly, Lz, X_size2, Y_size2, Z_size2, X_cell, Y_cell, Z_cell);
    
    T = T + dt * RHS_T; %calculate the salinity field at the next time step.
    
    
    %% new aV value:
    
    %param8:
    density_north = sum(sum(sum((1026*(1-alph * T(:,1:y_switchy,:) + bet * S(:,1:y_switchy,:))) .* boxvolume(:,1:y_switchy,:))))/sum(sum(sum(boxvolume(:,1:y_switchy,:))));
    density_south = sum(sum(sum((1026*(1-alph * T(:,y_switchy+1:end,:) + bet * S(:,y_switchy+1:end,:))).*boxvolume(:,y_switchy+1:end,:))))/sum(sum(sum(boxvolume(:,y_switchy+1:end,:))));
    
    
    %aV = Gamma * (density_south - density_north); <- tests 1,2
    
    aV = - Gamma * (density_south - density_north);
    aVsave(s1) = aV;
    
    s1 = s1 + 1;
    
    %% check we are still running:
    if sign(sum(sum(sum(isnan(T))))^2 + sum(sum(sum(isnan(S))))^2)
        temp = 1
    end
    
    %% can save every "saveevery" timestep
    if mod(t,saveevery) == 0 %save every "saveevery" time-steps
        t
    end
    
    Tave = Tave + T/length(timevec);
    Save = Save + S/length(timevec);
end

toc

if saveyesno
    save(name, 'aVsave', 'T', 'S');
end
