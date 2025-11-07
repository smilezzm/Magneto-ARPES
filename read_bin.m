% ARPES loader for SES 3D mapping (Ba122_fmap)
% Assumes .bin has no header and is float32 (4 bytes), little-endian.
% Dimensions and axes are taken from your .ini.
% 
% Don't forget to change the path of .bin before running the script and
% to change the .mat filename when saving the output data.
    function read_bin(bin_path, W, H, D, W_offset, H_offset, D_offset, W_delta, H_delta, D_delta)
    %% ---------- Parameters from your .ini ----------
    params = struct();
    
    % File
    params.bin_path   = bin_path; % 'D:/yao/Rice University/Academic/magneto ARPES/B0005/Spectrum_Ba122_fmap.bin';
    
    % Dimensions
    params.width      = W; % 618;     Energy points       
    params.height     = H; % 1000;    Thetax points       
    params.depth      = D; % 95;      Thetay points       
    
    % Data type
    params.bytes_per_point = 4; %  byteperpoint
    params.datatype  = 'single'; % float32
    params.byteorder = 'ieee-le'; % SES on Windows is little-endian
    
    % Axes calibration
    params.energy_offset = W_offset; % 15.825000;   % width_offset (Energy [eV])
    params.energy_delta  = W_delta; % 0.003000;    % width_delta
    params.thetax_offset = -24.942574;  % height_offset (Thetax [deg])
    params.thetax_delta  = 0.049935;    % height_delta
    params.thetay_offset = -14.000000;  % depth_offset (Thetay [deg])
    params.thetay_delta  = 0.300000;    % depth_delta
    params.Ef = 16.91;                   % Fermi Level ([eV])
    
    params.energy_label  = 'Energy [eV]';
    params.thetax_label  = 'Thetax [deg]';
    params.thetay_label  = 'Thetay [deg]';
    params.region_name   = 'map';
    
    
    %% ---------- Derived quantities ----------
    W = params.width;
    H = params.height;
    D = params.depth;
    N = W * H * D;
    Ef = params.Ef;
    
    expected_bytes = N * params.bytes_per_point;
    
    E  = params.energy_offset + (0:W-1) * params.energy_delta;
    Tx = params.thetax_offset + (0:H-1) * params.thetax_delta;
    Ty = params.thetay_offset + (0:D-1) * params.thetay_delta;
    
    %% ---------- File size sanity check ----------
    info = dir(params.bin_path);
    assert(~isempty(info), 'File not found: %s', params.bin_path);
    
    if info.bytes ~= expected_bytes
        warning('Byte count mismatch: file has %d bytes, expected %d. Check header/endianness/dimensions.', ...
            info.bytes, expected_bytes);
    end
    
    %% ---------- Read binary ----------
    fid = fopen(params.bin_path, 'r');
    assert(fid>0, 'Could not open file: %s', params.bin_path);
    
    cleaner = onCleanup(@() fclose(fid));
    
    raw = fread(fid, N, ['*' params.datatype], 0, params.byteorder);
    assert(numel(raw)==N, 'Read %d elements, expected %d.', numel(raw), N);
    
    %% ---------- Reshape ----------
    % Layout: [width (energy), height (Thetax), depth (Thetay)]
    data = reshape(raw, [W, H, D]);
    
    % Optional: move to [H x W x D] to make imagesc show Thetax vs Energy by default
    data_HWD = permute(data, [2 1 3]);  % [Thetax, Energy, Thetay]
    
    %% ---------- Quick sanity checks ----------
    mins = [min(data(:)), max(data(:))];
    fprintf('Data range: [%g, %g]\n', mins(1), mins(2));
    
    if ~all(isfinite(mins))
        warning('Data contains non-finite values. Endianness or dtype may be wrong.');
    end
    
    %% ---------- Save compact .mat ----------
    out = struct();
    out.data_HWD = data_HWD;   % preferred orientation for plotting
    out.data_raw = data;       % original orientation [W x H x D] [energy x thetax x thetay]
    out.energy   = E(:);
    out.thetax   = Tx(:);
    out.thetay   = Ty(:);
    out.labels   = struct('energy', params.energy_label, ...
                          'thetax', params.thetax_label, ...
                          'thetay', params.thetay_label);
    out.meta     = params;     % keep full metadata
    
    save_name = sprintf('%s_ARPES_%dx%dx%d.mat', params.region_name, W, H, D);
    save(save_name, '-struct', 'out', '-v7.3');
    fprintf('Saved: %s\n', save_name);
    
    %% ---------- Quick-look plots ----------
    figure('Name','Energy-Thetax slice (Thetay near 0)','Color','w');
    % Find Thetay index closest to 0 deg
    [~, iz0] = min(abs(Ty - 0));
    imagesc(E, Tx, data_HWD(:,:,iz0)); axis xy;
    xlabel(params.energy_label); ylabel(params.thetax_label);
    title(sprintf('%s | Thetay = %.3g deg (idx %d)', params.region_name, Ty(iz0), iz0));
    colorbar;
    
    figure('Name','Thetax-Thetay map at fermi surface','Color','w');
    % Pick energy near the Fermi level
    [~, iw_mid] = min(abs(E - Ef));
    % For Thetax-Thetay we want [Tx x Ty], so take data as [H x W x D] and squeeze along energy
    img_TxTy = squeeze(data_HWD(:, iw_mid, :));  % [H x D]
    imagesc(Ty, Tx, img_TxTy); axis xy;
    xlabel(params.thetay_label); ylabel(params.thetax_label);
    title(sprintf('%s | Energy = %.3f eV', params.region_name, E(iw_mid)));
    colorbar;
    
    %% Plot the fermi surface as (kx,ky,intensity) (it was (thetax,thetay,intensity) previously)
    % --- Constants ---
    me   = 9.10938356e-31;     % electron mass [kg]
    hbar = 1.054571817e-34;    % reduced Planck constant [J*s]
    e   = 1.602176634e-19;    % 1 eV in joules
    
    % --- Choose energy slice ---
    targetE = 16.91; % eV, example
    [~, iw] = min(abs(E - targetE)); % index in your energy axis
    
    % Extract Thetax-Thetay intensity slice
    I_thetax_thetay = squeeze(data_HWD(:, iw, :)); % [Thetax x Thetay]
    
    % --- Build angle grids ---
    [ThetaY_grid, ThetaX_grid] = ndgrid(Tx, Ty); % deg
    ThetaX_rad = deg2rad(ThetaX_grid);
    ThetaY_rad = deg2rad(ThetaY_grid);
    
    % --- Compute k magnitude ---
    Ekin_J = targetE * e; % kinetic energy in joules
    k_mag  = sqrt(2 * me * Ekin_J) / hbar; % [m^-1]
    
    % --- Map to kx, ky ---
    kx = k_mag .* cos(ThetaX_rad) .* sin(ThetaY_rad);
    ky = k_mag .* sin(ThetaX_rad);
    
    % --- Scale to 1e10 m^-1 ---
    kx = kx / 1e10;
    ky = ky / 1e10;
    
    % --- Plot ---
    figure('Color','w');
    scatter(kx(:), ky(:), 10, I_thetax_thetay(:), 'filled');
    axis equal tight;
    xlabel('k_x [10^{10} m^{-1}]');
    ylabel('k_y [10^{10} m^{-1}]');
    title(sprintf('kx-ky @ E = %.3f eV', Ef));
    colormap turbo; colorbar;
    
    save('fermi_surface_0mA.mat','kx','ky','I_thetax_thetay','Ef');