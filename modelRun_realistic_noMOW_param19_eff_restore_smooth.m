function modelRun_realistic_noMOW_param19_eff_restore_smooth(Time,saveyesno,name, T, S)

%IC19long:

%load('interim/modelRun_param19_eff_smooth_MIXED_IC19long_200yrs.mat', 'S', 'T')
load('SForcing_param19.mat', 'SForcing_IC19long_param19');
SForcing = SForcing_IC19long_param19;

load('SSTtypes.mat', 'SSTtypes');
T_atm_bare = SSTtypes.SSTzonal_5th_smooth;
T_atm_bare = T_atm_bare - 15.5;
X = SSTtypes.x;

load('SSStypes.mat', 'SSStypes');
S_atm_bare = SSStypes.SSSzonal_5th_smooth;
S_atm_bare = S_atm_bare - 15.5;
X = SSTtypes.x;

defineParameters_cartesian

%param14_eff:
%Gamma = 0.14;
%Gamma = 0.1;
Gamma = 0.065;
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

%evaluate salinity BC at y_eval:
S_atm = zeros(2000,Ly);
T_atm = zeros(2000,Ly);
for i = 1:2000
    S_atm(i,:) = interp1(X, S_atm_bare(:,i), y_eval);
    S_atm(i,:) = flip(S_atm(i,:));
    T_atm(i,:) = interp1(X, T_atm_bare(:,i), y_eval);
    T_atm(i,:) = flip(T_atm(i,:));
end

%all param except parma9:
dt = 0.0005;
%saveevery = 500*dt;
saveevery = 100*dt;

numberofsteps = length(dt:dt:Time); %calculate number of steps.
numberofsavedsteps = floor(Time/saveevery);

S0 = S;
%Smat = zeros([floor(size(S)/2) numberofsavedsteps]); %CAN LOOK AT PARTIAL DATA INSTEAD OF SIZE(S) TO SAVE SPACE
Smat = zeros([size(S) numberofsavedsteps]);

T0 = T;
%Tmat = zeros([floor(size(T)/2) numberofsavedsteps]); %CAN LOOK AT PARTIAL DATA INSTEAD OF SIZE(S) TO SAVE SPACE
Smat = zeros([size(S) numberofsavedsteps]);

s = 1; %matrix saving index

aVsave = zeros(numberofsteps+1,1);
%temp = abs(y_eval - ( yHbar - yHvariation));
%temp = abs(y_eval - 0.35);
temp = abs(y_eval - 0.22);
y_switchy = find(temp == min(temp));
y_switchy = y_switchy(1);
%y_switchy = 17;

boxvolume = X_size .* Y_size .* Z_size;

%param8:
density_north = sum(sum(sum((1026*(1-alph * T(:,1:y_switchy,:) + bet * S(:,1:y_switchy,:))) .* boxvolume(:,1:y_switchy,:))))/sum(sum(sum(boxvolume(:,1:y_switchy,:))));
density_south = sum(sum(sum((1026*(1-alph * T(:,y_switchy+1:end,:) + bet * S(:,y_switchy+1:end,:))).*boxvolume(:,y_switchy+1:end,:))))/sum(sum(sum(boxvolume(:,y_switchy+1:end,:))));


s1 = 1;
%aV = Gamma * (density_south - density_north); <- tests 1,2
aV = -Gamma * (density_south - density_north);
aVsave(s1) = aV;
s1 = s1 + 1;

%create velocity field:
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



% %MOW source:
% A = 2.5;
% sigmax = 0.2;
% sigmay = 0.1;
% sigmaz = 0.1;
% FS = abs(A*exp(...
%         -(X_eval-0.99).^2/(2*sigmax^2) ...
%         - (Y_eval-0.5).^2/(2*sigmay^2)...
%         - (Z_eval-0.75).^2/(2*sigmaz^2)));
% FS = FS - mean(mean(mean(FS)));
%
%
% A = 9;
% sigmax = 0.2;
% sigmay = 0.1;
% sigmaz = 0.1;
% FT = abs(A*exp(...
% -(X_eval-0.99).^2/(2*sigmax^2) ...
% - (Y_eval-0.5).^2/(2*sigmay^2)...
% - (Z_eval-0.75).^2/(2*sigmaz^2)));
% FT = FT - mean(mean(mean(FT))) - 12;
%

tic;
%
% f1 = figure;
% f2 = figure;
% f3 = figure;
% vSSS = VideoWriter('modelRun_param19_fullyear_SSS.avi');
% open(vSSS);
% vSST = VideoWriter('modelRun_param19_fullyear_SST.avi');
% open(vSST);
% vrho = VideoWriter('modelRun_param19_fullyear_rho.avi');
% open(vrho);
timevec = dt:dt:Time;
Tave = 0;
Save = 0;
for t = dt:dt:Time
    
    %taurel = floor(mod(t,Th)/dt)+1;
    
    yHp = yHvariation * sin(2*pi*t/Th);
    %yHp = 0;
    
    yH = yHbar + yHp;
    
    uH = uHtau .*sin(pi.*(2.*Yx+(-1).*yH));
    vH = vHtau .*sin(pi.*(Yy+(-1).*yH));
    
    u = aH * uH + aV * uV;
    v = aH * vH + aV * vV;
    w = aV * wV;
    
    RHS_S = 0;
    %calculate which timestep we're in:
    k = floor(mod(t * 2000, 2000));
    if k == 0
        k = 2000;
    end
    
    %RHS_S = diffusionEquation_oscBC_param19_mixed_northHighDiff(t, S, 1, Lx, Ly, Lz, X_eval, Y_eval, Z_eval, X_cell, Y_cell, Z_cell, SForcing(:,:,:,k),T_atm(k,:), northLat_ind);
    RHS_S = diffusionEquation_oscBC_param19_mixed_northHighDiff_restoring(t, S, 1, Lx, Ly, Lz, X_eval, Y_eval, Z_eval, X_cell, Y_cell, Z_cell, SForcing(:,:,:,k),T_atm(k,:), northLat_ind);
    %RHS_S = RHS_S + advectionEquation_multidim_FL77_SE(x_cell, y_cell, z_cell, t, S,u,v,w, dt);
    RHS_S = RHS_S + advectionEquation_multidim_FL77_SE_eff(S,u,v,w, dt, Lx, Ly, Lz, X_size2, Y_size2, Z_size2, X_cell, Y_cell, Z_cell);
    S = S + dt * RHS_S; %calculate the salinity field at the next time step.
    
    RHS_T = 0;
    RHS_T = diffusionEquation_oscBC_param19_mixed_northHighDiff_restoring(t, T, 2, Lx, Ly, Lz, X_eval, Y_eval, Z_eval, X_cell, Y_cell, Z_cell, SForcing(:,:,:,k),T_atm(k,:), northLat_ind);
    %RHS_T = RHS_T + advectionEquation_multidim_FL77_SE(x_cell, y_cell, z_cell, t, T,u,v,w, dt);
    RHS_T = RHS_T + advectionEquation_multidim_FL77_SE_eff(T,u,v,w, dt, Lx, Ly, Lz, X_size2, Y_size2, Z_size2, X_cell, Y_cell, Z_cell);
    T = T + dt * RHS_T; %calculate the salinity field at the next time step.
    
    %param8:
    density_north = sum(sum(sum((1026*(1-alph * T(:,1:y_switchy,:) + bet * S(:,1:y_switchy,:))) .* boxvolume(:,1:y_switchy,:))))/sum(sum(sum(boxvolume(:,1:y_switchy,:))));
    density_south = sum(sum(sum((1026*(1-alph * T(:,y_switchy+1:end,:) + bet * S(:,y_switchy+1:end,:))).*boxvolume(:,y_switchy+1:end,:))))/sum(sum(sum(boxvolume(:,y_switchy+1:end,:))));
    
    
    %aV = Gamma * (density_south - density_north); <- tests 1,2
    
    aV = - Gamma * (density_south - density_north);
    aVsave(s1) = aV;
    
    s1 = s1 + 1;
    %plot(aVsave(1:s1),'.')
    %     imagesc(lon_eval, lat_eval, 15.5+squeeze((T(end,:,:)))); colorbar
    %         set(gca, 'Ydir', 'normal')
    %         set(gca, 'Xdir', 'reverse')
    %         title(num2str(t))
    %         pause
    
    if sign(sum(sum(sum(isnan(T))))^2 + sum(sum(sum(isnan(S))))^2)
        temp = 1
    end
    
    
    if mod(t,saveevery) == 0 %save every "saveevery" time-steps
        t
        imagesc(lon_eval, lat_eval, 15.5+squeeze((T(end,:,:)))); colorbar
        set(gca, 'Ydir', 'normal')
        set(gca, 'Xdir', 'reverse')
        title(num2str(t))
        pause(0.01)
    end
    
    Tave = Tave + T/length(timevec);
    Save = Save + S/length(timevec);
end

toc

if saveyesno
    save(name, 'aVsave');
end
