function corrected_band=backward_mapping(final_measured, Fkx, Fky, do_plot)
    % Don't forget to change the location of saving
    kx = final_measured.kx;
    ky = final_measured.ky;
    corrected_band.I = final_measured.I_thetax_thetay;

    corrected_band.kx = Fkx(kx,ky);  % there may be nan
    corrected_band.ky = Fky(kx,ky);
    valid_idx = ~isnan(corrected_band.kx) & ~isnan(corrected_band.ky);
    corrected_band.kx = corrected_band.kx(valid_idx);  % vector
    corrected_band.ky = corrected_band.ky(valid_idx);  % vector
    corrected_band.I = corrected_band.I(valid_idx);

    % save('./matdata/corrected_band_110mA.mat', '-struct', 'corrected_band');

    if do_plot
        figure('Color','w');
        scatter(corrected_band.kx, corrected_band.ky, 10, corrected_band.I, 'filled');
        axis equal tight;
        xlabel('k_x [10^{10} m^{-1}]');
        ylabel('k_y [10^{10} m^{-1}]');
        title('backward mapping @ fermi surface');
        colormap turbo; colorbar;
    end