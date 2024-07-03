function RHS = diffusionEquation_oscBC_param19_mixed_SELECTDiff_mixed...
    (t, S, BCtype, Lx, Ly, Lz, y_eval, X_eval, Y_eval, Z_eval, X_cell, Y_cell, Z_cell,...
    S_atm, T_atm, northLat_ind, rho_invert_indicator)

%BCtype = type of boundary conditions.
%0 for no flux, 1 for RESTORING for salinity at surface,
% 2 for restoring flux at surface.

%Restore to zonally averaged rough data

global kappax kappay kappaz kappaz_enh northLat
 
% x_eval = (x_cell(1:end-1) + x_cell(2:end))/2;
% y_eval = (y_cell(1:end-1) + y_cell(2:end))/2;
% z_eval = (z_cell(1:end-1) + z_cell(2:end))/2;
%  
% Lx = length(x_eval);
% Ly = length(y_eval);
% Lz = length(z_eval);
%  
% [Y_eval,Z_eval,X_eval] = meshgrid(y_eval, z_eval, x_eval);
% [Y_cell,Z_cell,X_cell] = meshgrid(y_cell, z_cell, x_cell);
 
%diffusion coefficient matrices.
%defined at internal half-points
 
 
%Fx/Fy/Fz is the flux matrix in the X/Y/Z direction, respectively, defined
%at half points, of size Lx+1-Ly-Lz/Lx-Ly+1-Lz/Lx-Ly-Lz+1, with its edges providing the
%relevant boundary conditions:
Fx_diffusion = zeros(Lz, Ly, Lx+1);
Fy_diffusion = zeros(Lz, Ly+1, Lx);
Fz_diffusion = zeros(Lz+1, Ly, Lx);
 
%assign flux values inside bulk:
 
Fx_diffusion(:,:,2:end-1) = kappax * diff(S,1,3)./(diff(X_eval,1,3));
 
Fy_diffusion(:,2:end-1,:) = kappay * diff(S,1,2)./(diff(Y_eval,1,2));

% High diffusion where there are density inversions:

diffZ_matrix = diff(S,1,1)./(diff(Z_eval,1,1));
kappaz_matrix = kappaz + zeros(size(diffZ_matrix));
kappaz_matrix(rho_invert_indicator>0) = kappaz_enh;
Fz_diffusion(2:end-1,:,:) = kappaz_matrix .* diffZ_matrix;

%Fz_diffusion(2:end-1,1:northLat_ind,:) = kappaz_enh * diffZ_matrix(:,1:northLat_ind,:);
%Fz_diffusion(2:end-1,(northLat_ind+1):end,:) = kappaz * diffZ_matrix(:,(northLat_ind+1):end, :);
 
%Only on surface the boundary flux can be non-zero. There are three
%options:
%1. Zero flux: leave it as is. When ICtype = 0 this happens.
 
if BCtype == 1 %RESTORING BC for SALINITY!!!
    
    S_atmosphere = repmat(S_atm, 1, 1, Lx);
    
    %S_atmosphere = repmat(S_atm, 1, 1, Lx);
    %S_atmosphere = (2 - 0.2 * cos(2*pi*t))* cos(pi*(1-Y_eval(end,:,:)));
    
    %param13:
    kappaSurface = 0.01; %(52 days restoring time)
    
    
    %Fz_diffusion(end,:,:) = kappaSurface * (S_atmosphere - S(end,:,:));
    Fz_diffusion(end,:,:) = kappaSurface * S_atmosphere;
 
elseif BCtype == 2 %3. Flux depends on surface value: dS/dz = kappa(S_atmosphere-S): Fits temperature, around 15.5degC.

    S_atmosphere = repmat(T_atm, 1, 1, Lx);
 
    %param13:
    kappaSurface = 0.07; %(52 days restoring time)

    Fz_diffusion(end,:,:) = kappaSurface * (S_atmosphere - S(end,:,:));


end
 
RHS = diff(Fx_diffusion,1,3)./diff(X_cell(2:end,2:end,:),1,3) + ...
    diff(Fy_diffusion,1,2)./diff(Y_cell(2:end,:,2:end),1,2) + ...
    diff(Fz_diffusion,1,1)./diff(Z_cell(:,2:end,2:end),1,1);
 
 
 
end
 

