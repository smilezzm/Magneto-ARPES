function read_igor_h5(filepath, outputpath, W, H, D, W_delta, H_delta, D_delta, ...
    W_offset, H_offset, D_offset, Ef)
    % igor generated hdf5 file has no headers.
    % Input
    % W: energy
    % H: -18~18 angle (thetax/rotation around x)
    % D: -14~14 angle (thetay/rotation around y)
    % Ef: E in this h5 file is already translated so that E≈0 corresponds to
    %   fermi level. Ef is the kinetic energy of E=0 in this dataset
    % Output
    % - an interactive band figure at each energy slice
    % - save the 3D matrix in ./matdata/

    info = h5info(filepath);
    dataset = info.Datasets.Name; % is a string
    I_thetax_thetay = h5read(filepath, ['/', dataset]);  % should be (W,H,D)
    if ~isequal(size(I_thetax_thetay), [W,H,D])
        error('Make sure the W(energy),H(thetax),D(thetay) inputs are correct');
    end

    mins = [min(I_thetax_thetay(:)), max(I_thetax_thetay(:))];
    fprintf('Data range: [%g, %g]\n', mins(1), mins(2));
    
    if ~all(isfinite(mins))
        warning('Data contains non-finite values. Endianness or dtype may be wrong.');
    end

    me   = 9.10938356e-31;     % electron mass [kg]
    hbar = 1.054571817e-34;    % reduced Planck constant [J*s]
    e   = 1.602176634e-19;    % 1 eV in joules

    E = W_offset + Ef + (0:double(W)-1) * W_delta;
    Tx = H_offset + (0:double(H)-1) * H_delta;
    Ty = D_offset + (0:double(D)-1) * D_delta;

    Tx_rad = deg2rad(Tx);
    Ty_rad = deg2rad(Ty);
    [ThetaX_grid, ThetaY_grid] = ndgrid(Tx_rad, Ty_rad);
    % Notice here thetay is the rotation around x axis!
    
    % --- Compute k magnitude ---
    Ekin_J = E * e; % kinetic energy in joules
    k_mag  = sqrt(2 * me * Ekin_J) / hbar; % [m^-1]  it's a row vector
    
    % --- Map to kx, ky ---
    % ky = k_mag .* cos(ThetaY_grid) .* sin(ThetaX_grid);
    % kx = k_mag .* sin(ThetaY_grid);
    ThetaX_grid = reshape(ThetaX_grid, [1, H, D]);
    ThetaY_grid = reshape(ThetaY_grid, [1, H, D]);
    kx = reshape(k_mag, [W,1,1]) .* tan(ThetaY_grid).*cos(ThetaX_grid)./sqrt(1+tan(ThetaY_grid).^2.*cos(ThetaX_grid).^2);
    ky = reshape(k_mag, [W,1,1]) .* tan(ThetaX_grid).*cos(ThetaY_grid)./sqrt(1+tan(ThetaX_grid).^2.*cos(ThetaY_grid).^2);
    % x,y correspond to tilt axis and polar axis
    
    % --- Scale to 1e10 m^-1 ---
    kx = kx / 1e10;
    ky = ky / 1e10;
    
    % interactively show each energy slice
    fig = figure('Name', 'Interactive Band Viewer', ...
                 'NumberTitle', 'off', ...
                 'Position', [100 100 1200 550]);
    
    % Create axes for the plot, leaving room at bottom for slider
    ax = axes('Parent', fig, 'Position', [0.1 0.15 0.75 0.8]);
    
    kx_init = kx(1,:,:);
    ky_init = ky(1,:,:);
    I_thetax_thetay_init = I_thetax_thetay(1,:,:);
    band = scatter(ax, kx_init(:), ky_init(:), 10, I_thetax_thetay_init(:), 'filled');
    axis(ax, 'equal', 'tight');
    xlabel(ax, 'k_x [10^{10} m^{-1}]');
    ylabel(ax, 'k_y [10^{10} m^{-1}]');
    title(ax, sprintf('kx-ky @ E = %.3f eV (slice 1/%d)', E(1), W));
    colormap(ax, turbo);
    cb = colorbar(ax);
    cb.Label.String = 'Intensity';
    clim(ax, [mins(1), mins(2)]);
    
    % Create slider control
    slider = uicontrol('Parent', fig, 'Style', 'slider', ...
                      'Units', 'normalized', ...
                      'Position', [0.1 0.05 0.65 0.03], ...
                      'Min', 1, 'Max', W, 'Value', 1, ...
                      'SliderStep', [1/(W-1), 10/(W-1)]);
    
    % Create text label for slider
    slider_label = uicontrol('Parent', fig, 'Style', 'text', ...
                            'Units', 'normalized', ...
                            'Position', [0.76 0.04 0.15 0.05], ...
                            'String', sprintf('Slice: 1/%d', W), ...
                            'HorizontalAlignment', 'left', ...
                            'FontSize', 10);
    
    % store needed data in figure appdata for the callback
    setappdata(fig, 'kx', kx);
    setappdata(fig, 'ky', ky);
    setappdata(fig, 'I_thetax_thetay', I_thetax_thetay);
    setappdata(fig, 'E', E);
    setappdata(fig, 'mins', mins);
    setappdata(fig, 'band', band);
    setappdata(fig, 'ax', ax);
    setappdata(fig, 'slider_label', slider_label);
    
    % assign callback for slider
    addlistener(slider, 'ContinuousValueChange', @(src, event) sliderCallback(fig, src));
    
    
    save(outputpath,'kx','ky','I_thetax_thetay','E');

end

function sliderCallback(fig, slider)
    % Callback for slider: change energy slice and update scatter colors
    try
        % Get slider value and round to nearest integer
        idx = round(slider.Value);
        
        % Get stored data
        E = getappdata(fig, 'E');
        I_all = getappdata(fig, 'I_thetax_thetay');
        band = getappdata(fig, 'band');
        mins = getappdata(fig, 'mins');
        ax = getappdata(fig, 'ax');
        slider_label = getappdata(fig, 'slider_label');
        W = numel(E);

        % Ensure index is within bounds
        idx = max(1, min(W, idx));

        % Extract slice and update scatter colors
        I_slice = squeeze(I_all(idx, :, :));
        band.CData = I_slice(:);

        % Update title with energy
        if isvalid(ax)
            title(ax, sprintf('kx-ky @ E = %.3f eV (slice %d/%d)', E(idx), idx, W));
            % keep color limits constant
            clim(ax, [mins(1), mins(2)]);
        end
        
        % Update slider label
        if isvalid(slider_label)
            slider_label.String = sprintf('Slice: %d/%d', idx, W);
        end

        drawnow limitrate;
    catch ME
        warning('Slider callback failed: %s', ME.message);
    end
end
