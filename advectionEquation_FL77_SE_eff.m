function [RHS, Fx_advection, Fy_advection, Fz_advection] = advectionEquation_FL77_SE_eff(S,u,v,w, dt, Lx, Ly, Lz, X_size, Y_size, Z_size, X_cell, Y_cell, Z_cell)
%BCtype = type of boundary conditions.
%0 for no flux, 1 for constant flux at surface,
% 2 for restoring flux at surface.

%low order scheme: upwind
%high order scheme: Lax-Wendorff
%Using Stack-exchange scheme
%"https://math.stackexchange.com/questions/3321872/how-do-i-properly-implement-a-hyperbolic-tvd-flux-limiter"

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


%advection coefficient matrices.
%defined at internal half-points


%Fx/Fy/Fz is the flux matrix in the X/Y/Z direction, respectively, defined
%at half points, of size Lx+1-Ly-Lz/Lx-Ly+1-Lz/Lx-Ly-Lz+1, with its edges providing the
%relevant boundary conditions:

Fx_advection = zeros(Lz, Ly, Lx+1);
Fy_advection = zeros(Lz, Ly+1, Lx);
Fz_advection = zeros(Lz+1, Ly, Lx);

% Fx_advection_U = zeros(Lz, Ly, Lx+1);
% Fy_advection_U = zeros(Lz, Ly+1, Lx);
% Fz_advection_U = zeros(Lz+1, Ly, Lx);

%find velocity fields: 
%MUCH MORE EFFICIENT TO INSERT AS INPUTS!
%[u, v, w] = velocityFields(x_cell, y_cell, z_cell, t);

%slope function:
rx = ones(size(Fx_advection));
ry = ones(size(Fy_advection));
rz = ones(size(Fz_advection));

%if u is a matrix, change to u(:,:,2:end-1):
% rx(:,:,3:end-2) =  1/2 * (sign(u(:,:,2:end-1)) + 1)... %1 if u>0, 0 if u<0
%     .* (S(:,:,2:end-2) - S(:,:,1:end-3))./(S(:,:,3:end-1) - S(:,:,2:end-2))...
%     + 1/2 * (1 - sign(u(:,:,2:end-1)))... %1 if u<0, 0 if u>0
%     .* (S(:,:,4:end) - S(:,:,3:end-1))./(S(:,:,3:end-1) - S(:,:,2:end-2));

rx(:,:,3:end-2) =  (u(:,:,2:end-1) > 0)... %1 if u>0, 0 if u<0
    .* (diff(S(:,:,1:end-2),1,3)./diff(S(:,:,2:end-1),1,3))...
    + (u(:,:,2:end-1)<0)... %1 if u<0, 0 if u>0
    .* (diff(S(:,:,3:end),1,3)./diff(S(:,:,2:end-1),1,3));

% ry(:,3:end-2,:) =  1/2 * (sign(v(:,2:end-1,:)) + 1)... %1 if u>0, 0 if u<0
%     .* (S(:,2:end-2,:) - S(:,1:end-3,:))./(S(:,3:end-1,:) - S(:,2:end-2,:))...
%     + 1/2 * (1 - sign(v(:,2:end-1,:)))... %1 if u<0, 0 if u>0
%     .* (S(:,4:end,:) - S(:,3:end-1,:))./(S(:,3:end-1,:) - S(:,2:end-2,:));

ry(:,3:end-2,:) =  (v(:,2:end-1,:) > 0)... %1 if u>0, 0 if u<0
     .* (diff(S(:,1:end-2,:),1,2)./diff(S(:,2:end-1,:),1,2))...
     + (v(:,2:end-1,:)<0)... %1 if u<0, 0 if u>0
     .* (diff(S(:,3:end,:),1,2)./diff(S(:,2:end-1,:),1,2));

% rz(3:end-2,:,:) =  1/2 * (sign(w(2:end-1,:,:)) + 1)... %1 if u>0, 0 if u<0
%     .* (S(2:end-2,:,:) - S(1:end-3,:,:))./(S(3:end-1,:,:) - S(2:end-2,:,:))...
%     + 1/2 * (1 - sign(w(2:end-1,:,:)))... %1 if u<0, 0 if u>0
%     .* (S(4:end,:,:) - S(3:end-1,:,:))./(S(3:end-1,:,:) - S(2:end-2,:,:));

rz(3:end-2,:,:) = (w(2:end-1,:,:)>0)... %1 if u>0, 0 if u<0
    .* (diff(S(1:end-2,:,:),1,1)./diff(S(2:end-1,:,:),1,1))...
    + (w(2:end-1,:,:)<0)... %1 if u<0, 0 if u>0
    .* (diff(S(3:end,:,:),1,1)./diff(S(2:end-1,:,:),1,1));


%slope limiter - minmod limiter:
% phix = min(ones(size(rx)), rx);
% phix = max(zeros(size(rx)), phix);
% 
% phiy = min(ones(size(ry)), ry);
% phiy = max(zeros(size(ry)), phiy);
% 
% phiz = min(ones(size(rz)), rz);
% phiz = max(zeros(size(rz)), phiz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%slope limiter - Sweby limiter: THE BEST SO FAR (with MD, beta = 1.5)
%betaSweby = 1.5; BEST

betaSweby = 1.5;

phix1 = min(ones(size(rx)), rx * betaSweby);
phix2 = min(rx, betaSweby * ones(size(rx)));
phix = max(phix1, phix2);
%phix = max(phix, zeros(size(rx)));
phix = (phix > 0) .* phix;

phiy1 = min(ones(size(ry)), ry * betaSweby);
phiy2 = min(ry, betaSweby * ones(size(ry)));
phiy = max(phiy1, phiy2);
%phiy = max(phiy, zeros(size(ry)));
phiy = (phiy > 0) .* phiy;

phiz1 = min(ones(size(rz)), rz * betaSweby);
phiz2 = min(rz, betaSweby * ones(size(rz)));
phiz = max(phiz1, phiz2);
%phiz = max(phiz, zeros(size(rz)));
phiz = (phiz > 0) .* phiz;

%For upwind:
% phix = 2*zeros(size(phix));
% phiy = 2*zeros(size(phiy));
% phiz = 2*zeros(size(phiz));

%assign flux values inside bulk using first order upwind scheme (disadvantage: anisotropic numerical diffusion depending on velocity and gradients):

% Fx_advection_U(:,:,2:end-1) = - (min(u,0) .* S(:,:,2:end) + max(u,0) .* S(:,:,1:end-1));
% 
% 
% testx = -(u .* (S(:,:,1:end-1) + S(:,:,2:end))/2 - abs(u) .* (S(:,:,2:end) - S(:,:,1:end-1))/2);
%     
% 
% 
% Fy_advection_U(:,2:end-1,:) = - (min(v,0) .* S(:,2:end,:) + max(v,0) .* S(:,1:end-1,:));
% 
% testy = -(v .* (S(:,1:end-1,:) + S(:,2:end,:))/2 - abs(v) .* (S(:,2:end,:) - S(:,1:end-1,:))/2);
% 
% Fz_advection_U(2:end-1,:,:) = - (min(w,0) .* S(2:end,:,:) + max(w,0) .* S(1:end-1,:,:));
% 
% testz = -(w .* (S(1:end-1,:,:) + S(2:end,:,:))/2 - abs(w) .* (S(2:end,:,:) - S(1:end-1,:,:))/2);
% 
% Fx_advection(:,:,2:end-1) = Fx_advection_U(:,:,2:end-1) +...
%     phix(:,:,2:end-1) .* (sign(u) - u .* dt./X_size) .* u .* (S(:,:,2:end) - S(:,:,1:end-1))/2;
% 
% Fy_advection(:,2:end-1,:) = Fy_advection_U(:,2:end-1,:) +...
%     phiy(:,2:end-1,:) .* (sign(v) - v .* dt./Y_size) .* v .* (S(:,2:end,:) - S(:,1:end-1,:))/2;
% 
% Fz_advection(2:end-1,:,:) = Fz_advection_U(2:end-1,:,:) +...
%     phiz(2:end-1,:,:) .* (sign(w) - w .* dt./Z_size) .* w .* (S(2:end,:,:) - S(1:end-1,:,:))/2;
% 


%  Fx_advection(:,:,2:end-1) = -(u .* (S(:,:,1:end-1) + S(:,:,2:end))/2 - abs(u) .* (S(:,:,2:end) - S(:,:,1:end-1))/2) ...
%      - phix(:,:,2:end-1) .* (sign(u) - u .* dt./X_size) .* u .* (S(:,:,2:end) - S(:,:,1:end-1))/2;
% 
% Fx_advection(:,:,2:end-1) = -(u .* (S(:,:,1:end-1) + S(:,:,2:end))/2 - abs(u) .* (S(:,:,2:end) - S(:,:,1:end-1))/2) ...
%     - phix(:,:,2:end-1) .* (sign(u) - u .* dt./X_size) .* u .* (S(:,:,2:end) - S(:,:,1:end-1))/2;

% 
% Fy_advection(:,2:end-1,:) = -(v .* (S(:,1:end-1,:) + S(:,2:end,:))/2 - abs(v) .* (S(:,2:end,:) - S(:,1:end-1,:))/2) ...
%     - phiy(:,2:end-1,:) .* (sign(v) - v .* dt./Y_size) .* v .* (S(:,2:end,:) - S(:,1:end-1,:))/2;
% 
% Fz_advection(2:end-1,:,:) = -(w .* (S(1:end-1,:,:) + S(2:end,:,:))/2 - abs(w) .* (S(2:end,:,:) - S(1:end-1,:,:))/2) ...
%     - phiz(2:end-1,:,:) .* (sign(w) - w .* dt./Z_size) .* w .* (S(2:end,:,:) - S(1:end-1,:,:))/2;

Fx_advection(:,:,2:end-1) = -(u .* (S(:,:,1:end-1) + S(:,:,2:end))/2 - abs(u) .* diff(S,1,3)/2) ...
- phix(:,:,2:end-1) .* (sign(u) - u .* dt./X_size) .* u .* (diff(S,1,3))/2;

Fy_advection(:,2:end-1,:) = -(v .* (S(:,1:end-1,:) + S(:,2:end,:))/2 - abs(v) .* (diff(S,1,2))/2) ...
    - phiy(:,2:end-1,:) .* (sign(v) - v .* dt./Y_size) .* v .* (diff(S,1,2))/2;

Fz_advection(2:end-1,:,:) = -(w .* (S(1:end-1,:,:) + S(2:end,:,:))/2 - abs(w) .* (diff(S,1,1))/2) ...
    - phiz(2:end-1,:,:) .* (sign(w) - w .* dt./Z_size) .* w .* (diff(S,1,1))/2;

%Calculate the RHS:
RHS = diff(Fx_advection,1,3)./diff(X_cell(2:end,2:end,:),1,3) + ...
    diff(Fy_advection,1,2)./diff(Y_cell(2:end,:,2:end),1,2) + ...
    diff(Fz_advection,1,1)./diff(Z_cell(:,2:end,2:end),1,1);
    
