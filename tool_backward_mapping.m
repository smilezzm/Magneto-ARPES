if ~(exist('kx','var')&&exist('ky','var')&&exist('I_thetax_thetay','var')...
        &&exist('targetE','var'))
    load('fermi_surface.mat');
end
if ~(exist('Fkx','var')&&exist('Fky','var'))
    load('inverse_mapping.mat');
end

kx_initial = Fkx(kx,ky);
ky_initial = Fky(kx,ky);

figure('Color','w');
scatter(kx_initial(:), ky_initial(:), 10, I_thetax_thetay(:), 'filled');
axis equal tight;
xlabel('k_x [10^{10} m^{-1}]');
ylabel('k_y [10^{10} m^{-1}]');
title(sprintf('backward mapping @ E = %.3f eV', targetE));
colormap turbo; colorbar;