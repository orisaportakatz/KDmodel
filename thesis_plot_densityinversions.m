% plot density figures for thesis density inversion section:

load('/Volumes/orika/data/MOWpathways/SODAdata2/monthlyData_MATfiles/density_climatological.mat')
load('/Volumes/orika/data/MOWpathways/SODAdata2/monthlyData_MATfiles/climatological_monthly_NorthAtlantic_salt.mat')
load('/Volumes/orika/data/MOWpathways/SODAdata2/monthlyData_MATfiles/climatological_monthly_NorthAtlantic_temperature.mat')

%% averaging over longitudes to obtain density vs. depth:

%% average over longitudes:
lat1_ind = find(latitude == 10);
lat2_ind = find(latitude == 65);

prho_climMon_ave = zeros(lat2_ind-lat1_ind+1, 50, 12);


% from 10°N to 65°N
% 10°N-20°N, 70°W-5°W
lat1_ind = find(latitude == 10);
lat2_ind = find(latitude == 20);
lon1_ind = find(longitude == -70);
lon2_ind = find(longitude == -5);
numind_p1 = lat2_ind-lat1_ind+1;

prho_climMon_ave(1:numind_p1,:, :) =...
    squeeze(mean(prho_climMon_nan(lon1_ind:lon2_ind, lat1_ind:lat2_ind, :, :), 1, 'omitnan'));

%
% 20°N-50°N, 80°W-5°W
lat1_ind = find(latitude == 20)+1;
lat2_ind = find(latitude == 50);
lon1_ind = find(longitude == -80);
lon2_ind = find(longitude == -5);
numind_p2 = lat2_ind-lat1_ind+1;

prho_climMon_ave(numind_p1+(1:numind_p2), :, :) = ...
    squeeze(mean(prho_climMon_nan(lon1_ind:lon2_ind, lat1_ind:lat2_ind, :, :), 1, 'omitnan'));

% 50°N-65°N, 65°W-5°W
lat1_ind = find(latitude == 50)+1;
lat2_ind = find(latitude == 65);
lon1_ind = find(longitude == -65);
lon2_ind = find(longitude == -5);
numind_p3 = lat2_ind-lat1_ind+1;

prho_climMon_ave(numind_p1+numind_p2+(1:numind_p3), :, :) =...
    squeeze(mean(prho_climMon_nan(lon1_ind:lon2_ind, lat1_ind:lat2_ind, :, :), 1, 'omitnan'));




%% creating the cosine-IC for thesis density inversion runs:

% prepare IC for kinematic-dynamic model:
addpath(genpath('/Users/admin/Documents/MATLAB/oceanography'))
mowstartup

load('/Volumes/orika/data/MOWpathways/SODAdata2/monthlyData_MATfiles/climatological_monthly_NorthAtlantic_temperature.mat', 'tem_all_surface')
load('/Volumes/orika/data/MOWpathways/SODAdata2/monthlyData_MATfiles/climatological_monthly_NorthAtlantic_salt.mat')
%%
salt_all_surface = squeeze(salt_climMon(:,:,1,:));

figure; contourf(longitude, latitude, tem_all_surface(:,:,1)', 200,'LineColor','none')
colorbar
grid on

figure; contourf(longitude, latitude, salt_climMon(:,:,1,1)', 200,'LineColor','none')
colorbar
caxis([32 37])

%% average over longitudes:
lat1_ind = find(latitude == 10);
lat2_ind = find(latitude == 65);

tem_lat_climMon = zeros(lat2_ind-lat1_ind+1, 12);
salt_lat_climMon = zeros(lat2_ind-lat1_ind+1, 12);

% from 10°N to 65°N
% 10°N-20°N, 70°W-5°W
lat1_ind = find(latitude == 10);
lat2_ind = find(latitude == 20);
lon1_ind = find(longitude == -70);
lon2_ind = find(longitude == -5);
numind_p1 = lat2_ind-lat1_ind+1;

tem_lat_climMon(1:numind_p1,:) =...
    squeeze(mean(tem_all_surface(lon1_ind:lon2_ind, lat1_ind:lat2_ind, :), 1, 'omitnan'));
salt_lat_climMon(1:numind_p1,:) =...
    squeeze(mean(salt_all_surface(lon1_ind:lon2_ind, lat1_ind:lat2_ind, :), 1, 'omitnan'));

%
% 20°N-50°N, 80°W-5°W
lat1_ind = find(latitude == 20)+1;
lat2_ind = find(latitude == 50);
lon1_ind = find(longitude == -80);
lon2_ind = find(longitude == -5);
numind_p2 = lat2_ind-lat1_ind+1;

tem_lat_climMon(numind_p1+(1:numind_p2),:) =...
    squeeze(mean(tem_all_surface(lon1_ind:lon2_ind, lat1_ind:lat2_ind, :), 1, 'omitnan'));
salt_lat_climMon(numind_p1+(1:numind_p2),:) =...
    squeeze(mean(salt_all_surface(lon1_ind:lon2_ind, lat1_ind:lat2_ind, :), 1, 'omitnan'));


% 50°N-65°N, 65°W-5°W
lat1_ind = find(latitude == 50)+1;
lat2_ind = find(latitude == 65);
lon1_ind = find(longitude == -65);
lon2_ind = find(longitude == -5);
numind_p3 = lat2_ind-lat1_ind+1;

tem_lat_climMon(numind_p1+numind_p2+(1:numind_p3),:) =...
    squeeze(mean(tem_all_surface(lon1_ind:lon2_ind, lat1_ind:lat2_ind, :), 1, 'omitnan'));
salt_lat_climMon(numind_p1+numind_p2+(1:numind_p3),:) =...
    squeeze(mean(salt_all_surface(lon1_ind:lon2_ind, lat1_ind:lat2_ind, :), 1, 'omitnan'));

tem_max = max(tem_lat_climMon);
tem_min = min(tem_lat_climMon);

salt_max = max(salt_lat_climMon);
salt_min = min(salt_lat_climMon);

%% temperature:
% create cosines for restoring boundary conditions, to be used in order to create the spin-up of the model:

lataxis_model = 0:0.01:1;
temp1 = repmat(1/2*(1+cos(lataxis_model*pi)),12,1);
temp2 = repmat((tem_max - tem_min),101,1)';
tem_cosine = temp1 .* temp2 + repmat(tem_min, 101,1)';
figure; plot(tem_cosine, lataxis_model, 'x-')

%% sanity check:
figure; plot(tem_lat_climMon, latitude(42:152))
hold on; plot(tem_cosine, linspace(latitude(42), latitude(152), 101), 'x-')


%% salinity:
% create cosines for restoring boundary conditions, to be used in order to create the spin-up of the model:

lataxis_model = 0:0.01:1;
temp1 = repmat(1/2*(1+cos(lataxis_model*pi)),12,1);
temp2 = repmat((salt_max - salt_min),101,1)';
salt_cosine = temp1 .* temp2 + repmat(salt_min, 101,1)';
figure; plot(salt_cosine, lataxis_model, 'x-')

%% sanity check:
figure; plot(salt_lat_climMon, latitude(42:152))
hold on; plot(salt_cosine, linspace(latitude(42), latitude(152), 101), 'x-')

%% save these restoring boundary conditions:
save('IC_cosine_tem_salt.mat', 'salt_cosine', 'tem_cosine', 'lataxis_model')
