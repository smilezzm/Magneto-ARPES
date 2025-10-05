function [kx_corrected, ky_corrected, I_corrected]=backward_mapping(final_measured, Fkx, Fky)
    kx = final_measured.kx;
    ky = final_measured.ky;
    I_corrected = final_measured.I_thetax_thetay;

    kx_corrected = Fkx(kx,ky);
    ky_corrected = Fky(kx,ky);
    
    figure('Color','w');
    scatter(kx_corrected(:), ky_corrected(:), 10, I_corrected(:), 'filled');
    axis equal tight;
    xlabel('k_x [10^{10} m^{-1}]');
    ylabel('k_y [10^{10} m^{-1}]');
    title('backward mapping @ fermi surface');
    colormap turbo; colorbar;
