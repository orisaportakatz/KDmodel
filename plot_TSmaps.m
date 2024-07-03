set(groot,'defaultAxesFontName','Cambria')
set(groot,'defaultAxesFontSize',12)
set(groot, 'defaultfigurecolor', [1,1,1])

defineParameters_cartesian_thesis
timevec = dt * (0:(len(aVsave)-1));
load('thesis_spinuptuning_clustertrial_03_mixed_SELECTdiff_t200_part2.mat')

figure; contourf(lon_eval, lat_eval, Toffset + squeeze(T(end, :,:)), 200, 'LineColor', 'none')
set(gca, 'xdir', 'reverse')
colormap(jet)
colorbar
caxis([0 30])

figure; contourf(lon_eval, lat_eval, Soffset + squeeze(S(end, :,:)), 200, 'LineColor', 'none')
set(gca, 'xdir', 'reverse')
colormap(jet)
colorbar
caxis([34.5 36.5])
%%
rho = 1026.3*(1-alph*T+bet*S);
figure; contourf(lat_eval, depth_eval, squeeze(mean(rho,3)), 30, 'LineColor', 'none')
set(gca, 'ydir', 'reverse')
colorbar
colormap(jet)
caxis([1025 1029])

%%
figure; plot(timevec, 3.1e3*aVsave, 'x-')