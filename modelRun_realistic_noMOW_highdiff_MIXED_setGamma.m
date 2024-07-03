function [aVsave, T, S] = modelRun_realistic_noMOW_highdiff_MIXED_setGamma(Time,saveyesno,name, T, S, Gamma)
%param14_eff:
%Gamma = 0.065;
disp('updated 12:58')

%% prepare data
%load('interim/modelRun_param19_eff_smooth_MIXED_IC19long_200yrs.mat', 'S', 'T')
%load('IC_cosine_tem_salt.mat', 'tem_cosine', 'salt_cosine', 'lataxis_model')

%% parameters:
defineParameters_cartesian_thesis

%param14_eff:
%Gamma = 0.065;
alph = 2.2e-4;

Gamma

yHvariation

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


%all param except parma9:
dt = 0.0005;
%saveevery = 500*dt;
saveevery = 100*dt;

numberofsteps = length(dt:dt:Time); %calculate number of steps.
numberofsavedsteps = floor(Time/saveevery);
numberofstepsinayear = length(dt:dt:1);

%% boundary conditions - restoring:

load('IC_cosine_tem_salt.mat', 'tem_cosine', 'salt_cosine', 'lataxis_model')

%evaluate salinity BC at y_eval:
S_atm_month = zeros(12, Ly);
T_atm_month = zeros(12, Ly);
for i = 1:12
    S_atm_month(i,:) = interp1(lataxis_model, fliplr(salt_cosine(i,:)), y_eval);
    T_atm_month(i,:) = interp1(lataxis_model, fliplr(tem_cosine(i,:)), y_eval);
end

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

T_atm = T_atm - Toffset;
S_atm = S_atm - Soffset;

%% Initial conditions:

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
% bamba1 = figure;
% bamba2 = figure;
% bamba3 = figure;


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
    RHS_S = diffusionEquation_oscBC_param19_mixed_northHighDiff_mixed(t, S, 1, Lx, Ly, Lz, y_eval, X_eval, Y_eval, Z_eval, X_cell, Y_cell, Z_cell, S_atm(k,:),T_atm(k,:), northLat_ind);
    %RHS_S = RHS_S + advectionEquation_multidim_FL77_SE(x_cell, y_cell, z_cell, t, S,u,v,w, dt);
    RHS_S = RHS_S + advectionEquation_multidim_FL77_SE_eff(S,u,v,w, dt, Lx, Ly, Lz, X_size2, Y_size2, Z_size2, X_cell, Y_cell, Z_cell);
    S = S + dt * RHS_S; %calculate the salinity field at the next time step.
    
    RHS_T = 0;
    RHS_T = diffusionEquation_oscBC_param19_mixed_northHighDiff_mixed(t, T, 2, Lx, Ly, Lz, y_eval, X_eval, Y_eval, Z_eval, X_cell, Y_cell, Z_cell,S_atm(k,:),T_atm(k,:), northLat_ind);
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
        
%         figure(bamba1)
%         imagesc(lon_eval, lat_eval, 15.8+squeeze((T(end,:,:)))); colorbar
%         set(gca, 'Ydir', 'normal')
%         set(gca, 'Xdir', 'reverse')
%         title(num2str(t))
%         
%         figure(bamba2)
%         imagesc(lon_eval, lat_eval, 35.3+squeeze((S(end,:,:)))); colorbar
%         set(gca, 'Ydir', 'normal')
%         set(gca, 'Xdir', 'reverse')
%         title(num2str(t))
%         
%         figure(bamba3)
%         rho = 1026 * (1 - alph * T + bet * S + 1/(1026 * 1500^2) );
%         depth_eval = -(1-z_eval) * 4000;
%         contourf(lat_eval, depth_eval, squeeze(mean(rho(:,:,:),3))); colorbar
%         title(num2str(t))
%         
%         pause(0.01)
    end
    
    Tave = Tave + T/length(timevec);
    Save = Save + S/length(timevec);
end

toc

if saveyesno
    save(name, 'aVsave', 'T', 'S');
end
