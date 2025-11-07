% Don't forget to change the path of .h5 before running the script and
% to change the .mat filename when saving the output data.
% This script reads all the slices of the .h5 file, while read_h5_Ef.m
% reads only the fermi surface slice.
function [E,kx,ky,I_thetax_thetay] = read_h5(filename, outputpath)
    % kx,ky being [W=751,H=876,D=113] matrices
    % E being a [W=751] vector
    % I_thetax_thetay being a [W, H, D] matrix
    % outputpath = './matdata/150K_bands_-110mA.mat'
    if nargin < 1 || isempty(filename)
        filename = 'D:\yao\Rice University\Academic\magneto ARPES\C_0111.h5';
    end
    %% ---------- Parameters from your .ini ----------
    params = struct();
    
    %% File Structure
    params.filename = filename;  
    info = h5info(params.filename);
    disp('=== File Structure ===');
    disp(info);

    print_h5_tree(filename);
    
    %% Dimensions
    if strcmp(info.Groups(2).Groups(1).Attributes(3).Name, 'Count') && ...
        strcmp(info.Groups(2).Groups(1).Attributes(5).Name, 'Label') && ...
        strcmp(info.Groups(2).Groups(1).Attributes(5).Value, 'Kinetic Energy')
        params.axes0      = info.Groups(2).Groups(1).Attributes(3).Value;    % Energy point number (751)
    else
        error('the structure of this .h5 is different from default')
    end
    if strcmp(info.Groups(2).Groups(2).Attributes(3).Name, 'Count') && ...
        strcmp(info.Groups(2).Groups(2).Attributes(5).Name, 'Label') && ...
        strcmp(info.Groups(2).Groups(2).Attributes(5).Value, 'ThetaX')
        params.axes1      = info.Groups(2).Groups(2).Attributes(3).Value;    % thetaX point number (876)
    else
        error('the structure of this .h5 is different from default')
    end
    if strcmp(info.Groups(2).Groups(3).Attributes(3).Name, 'Count') && ...
        strcmp(info.Groups(2).Groups(3).Attributes(5).Name, 'Label') && ...
        strcmp(info.Groups(2).Groups(3).Attributes(5).Value, 'ThetaY')
        params.axes2      = info.Groups(2).Groups(3).Attributes(3).Value;    % thetaY point number (113)
    else
        error('the structure of this .h5 is different from default')
    end
    
    %% Axes calibration
    params.energy_offset = info.Groups(2).Groups(1).Attributes(1).Value;   % axes0_offset (Energy [eV]) (109.2350)
    params.energy_delta  = info.Groups(2).Groups(1).Attributes(2).Value;    % axes0_delta (0.0020)
    params.thetax_offset = info.Groups(2).Groups(2).Attributes(1).Value;  % axes1_offset (ThetaX [deg]) (-17.9795)
    params.thetax_delta  = info.Groups(2).Groups(2).Attributes(2).Value;   % axes1_delta (0.0411)
    params.thetay_offset = info.Groups(2).Groups(3).Attributes(1).Value;  % axes2_offset (Thetay [deg]) (-14)
    params.thetay_delta  = info.Groups(2).Groups(3).Attributes(2).Value;    % axes2_delta (0.25)
    
    params.energy_label  = 'Energy [eV]';
    params.thetax_label  = 'Thetax [deg]';
    params.thetay_label  = 'Thetay [deg]';
    params.region_name   = 'Kagome';
    
    %% ---------- Derived quantities ----------
    W = params.axes0;   % Energy dimensions 
    H = params.axes1;   % ThetaX dimensions [-18,18]
    D = params.axes2;   % ThetaY dimensions [-14,14]
    N = W * H * D;
    
    E  = params.energy_offset + (0:double(W)-1) * params.energy_delta;
    Tx = params.thetax_offset + (0:double(H)-1) * params.thetax_delta;
    Ty = params.thetay_offset + (0:double(D)-1) * params.thetay_delta;
    
    %% ---------- Read .h5 ----------
    raw = h5read(params.filename, '/Data/Count');
    assert(numel(raw)==N, 'Read %d elements, expected %d.', numel(raw), N);
    disp('=== Dataset Info ===');
    disp(size(raw));
    disp(class(raw));
    
    %% ---------- Reshape ----------
    % raw Layout: [depth (Thetay), height (Thetax), width (energy)] !!!!
    
    %% ---------- Quick sanity checks ----------
    mins = [min(raw(:)), max(raw(:))];
    fprintf('Data range: [%g, %g]\n', mins(1), mins(2));
    
    if ~all(isfinite(mins))
        warning('Data contains non-finite values. Endianness or dtype may be wrong.');
    end
    
    
    %% Plot the fermi surface as (kx,ky,intensity) (it was (thetax,thetay,intensity) previously)
    % --- Constants ---
    me   = 9.10938356e-31;     % electron mass [kg]
    hbar = 1.054571817e-34;    % reduced Planck constant [J*s]
    e   = 1.602176634e-19;    % 1 eV in joules
    
    
    % Extract Thetax-Thetay intensity slice
    I_thetax_thetay = permute(raw, [3 2 1]);   % [Energy x ThetaX x ThetaY]
    
    % --- Build angle grids ---
    Tx_rad = deg2rad(Tx);
    Ty_rad = deg2rad(Ty);
    [ThetaX_grid, ThetaY_grid] = ndgrid(Tx_rad, Ty_rad); 
    % Notice here thetay is the rotation around polar axis y, which contributes to kx
    
    % --- Compute k magnitude ---
    Ekin_J = E * e; % kinetic energy in joules
    k_mag  = sqrt(2 * me * Ekin_J) / hbar; % [m^-1]  it's a row vector
    
    % --- Map to kx, ky ---
    % ky = k_mag .* cos(ThetaY_grid) .* sin(ThetaX_grid);
    % kx = k_mag .* sin(ThetaY_grid);
    ThetaX_grid = reshape(ThetaX_grid, [1, H, D]);
    ThetaY_grid = reshape(ThetaY_grid, [1, H, D]);
    ky = reshape(k_mag, [W,1,1]) .* tan(ThetaX_grid).*cos(ThetaY_grid)./sqrt(1+tan(ThetaX_grid).^2.*cos(ThetaY_grid).^2);
    kx = reshape(k_mag, [W,1,1]) .* tan(ThetaY_grid).*cos(ThetaX_grid)./sqrt(1+tan(ThetaY_grid).^2.*cos(ThetaX_grid).^2);
    % x,y correspond to tilt axis and polar axis
    
    % --- Scale to 1e10 m^-1 ---
    kx = kx / 1e10;
    ky = ky / 1e10;
    
    % interactively show each energy slice
    fig = figure('Name', 'Interactive Band Viewer', ...
                 'NumberTitle', 'off', ...
                 'Position', [100 100 1200 500]);
    kx_init = kx(1,:,:);
    ky_init = ky(1,:,:);
    I_thetax_thetay_init = I_thetax_thetay(1,:,:);
    band = scatter(kx_init(:), ky_init(:), 10, I_thetax_thetay_init(:), 'filled');
    axis(ax, 'equal', 'tight');
    xlabel(ax, 'k_x [10^{10} m^{-1}]');
    ylabel(ax, 'k_y [10^{10} m^{-1}]');
    title(ax, sprintf('kx-ky @ E = %.3f eV (slice 1/%d)', E(1), W));
    colormap(ax, turbo);
    cb = colorbar(ax);
    cb.Label.String = 'Intensity';
    clim(ax, [mins(1), mins(2)]);

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


function print_h5_tree(filename)
    info = h5info(filename);
    fprintf('📂 %s\n', filename);
    walk_h5(info, '', true);
end

function walk_h5(group, indent, isRoot)
    if ~isRoot
        fprintf('%s📁 %s\n', indent, group.Name);
    end

    % Print attributes
    for i = 1:length(group.Attributes)
        fprintf('%s  🏷️ Attribute: %s\n', indent, group.Attributes(i).Name);
    end

    % Print datasets
    for i = 1:length(group.Datasets)
        fprintf('%s  📄 Dataset: %s\n', indent, group.Datasets(i).Name);
    end

    % Recurse into subgroups
    for i = 1:length(group.Groups)
        walk_h5(group.Groups(i), [indent '  '], false);
    end
end
