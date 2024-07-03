function RHS = advectionEquation_multidim_FL77_SE_eff(S,u,v,w, dt, Lx, Ly, Lz, X_size, Y_size, Z_size, X_cell, Y_cell, Z_cell)
%BCtype = type of boundary conditions.
%0 for no flux, 1 for constant flux at surface,
% 2 for restoring flux at surface.


%x_eval = (x_cell(1:end-1) + x_cell(2:end))/2;
%y_eval = (y_cell(1:end-1) + y_cell(2:end))/2;
%z_eval = (z_cell(1:end-1) + z_cell(2:end))/2;



%Lx = length(x_eval);
%Ly = length(y_eval);
%Lz = length(z_eval);

%[Y_eval,Z_eval,X_eval] = meshgrid(y_eval, z_eval, x_eval);
%[Y_cell,Z_cell,X_cell] = meshgrid(y_cell, z_cell, x_cell);

%X_size = diff(X_eval,1,3);
%Y_size = diff(Y_eval,1,2);
%Z_size = diff(Z_eval,1,1);



%implement MITGCM's multidimensional scheme for RHS calculation:

%BCtype = type of boundary conditions.
%[~, Fx_advection_1, ~, ~] = advectionEquation_FL77_SE_eff(x_cell, y_cell, z_cell, t, S,u,v,w, dt);
[~, Fx_advection_1, ~, ~] = advectionEquation_FL77_SE_eff(S,u,v,w, dt, Lx, Ly, Lz, X_size, Y_size, Z_size, X_cell, Y_cell, Z_cell);

S_1ov3 = S - dt * (diff(Fx_advection_1,1,3)./diff(X_cell(2:end,2:end,:),1,3));

%[~, ~, Fy_advection_2, ~] = advectionEquation_FL77_SE_eff(x_cell, y_cell, z_cell, t, S_1ov3,u,v,w, dt);
[~, ~, Fy_advection_2, ~] = advectionEquation_FL77_SE_eff(S,u,v,w, dt, Lx, Ly, Lz, X_size, Y_size, Z_size, X_cell, Y_cell, Z_cell);

S_2ov3 = S_1ov3 - dt * diff(Fy_advection_2,1,2)./diff(Y_cell(2:end,:,2:end),1,2);

%[~, ~, ~, Fz_advection_3] = advectionEquation_FL77_SE_eff(x_cell, y_cell, z_cell, t, S_2ov3,u,v,w, dt);
[~, ~, ~, Fz_advection_3] = advectionEquation_FL77_SE_eff(S,u,v,w, dt, Lx, Ly, Lz, X_size, Y_size, Z_size, X_cell, Y_cell, Z_cell);

S_3ov3 = S_2ov3 - dt * diff(Fz_advection_3,1,1)./diff(Z_cell(:,2:end,2:end),1,1);

RHS = -1/dt * (S_3ov3 - S);