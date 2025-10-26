% Don't forget to change the path of .h5 before running the script and
% to change the .mat filename when saving the output data.
function [kx,ky,I_thetax_thetay] = read_h5_Ef(filename, Ef)
    if nargin < 2 || isempty(Ef)
        Ef=110.56;
    end
    if nargin < 1 || isempty(filename)
        filename = 'D:\yao\Rice University\Academic\magneto ARPES\C_0004.h5';
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
    params.Ef = Ef;                   % Fermi Level ([eV]) (manually assigned)
    
    params.energy_label  = 'Energy [eV]';
    params.thetax_label  = 'Thetax [deg]';
    params.thetay_label  = 'Thetay [deg]';
    params.region_name   = 'Kagome';
    
    %% ---------- Derived quantities ----------
    W = params.axes0;
    H = params.axes1;
    D = params.axes2;
    N = W * H * D;
    Ef = params.Ef;
    
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
    
    %% ---------- Quick-look plots ----------
    figure('Color','w');
    % Pick energy near the Fermi level
    [~, iw_mid] = min(abs(E - Ef));
    img_TxTy = squeeze(raw(:, :, iw_mid))';  % transpose to [H × D]
    imagesc(Ty, Tx, img_TxTy);
    xlabel(params.thetay_label);
    ylabel(params.thetax_label);
    title(sprintf('%s | Energy = %.3f eV', params.region_name, E(iw_mid)));
    axis xy;
    colorbar; colormap turbo;
    
    
    %% Plot the fermi surface as (kx,ky,intensity) (it was (thetax,thetay,intensity) previously)
    % --- Constants ---
    me   = 9.10938356e-31;     % electron mass [kg]
    hbar = 1.054571817e-34;    % reduced Planck constant [J*s]
    e   = 1.602176634e-19;    % 1 eV in joules
    
    % --- Choose energy slice ---
    [~, iw] = min(abs(E - Ef)); % index in your energy axis
    
    % Extract Thetax-Thetay intensity slice
    I_thetax_thetay = squeeze(raw(:, :, iw))'; % [Thetax x Thetay]
    
    % --- Build angle grids ---
    Tx_rad = deg2rad(Tx);
    Ty_rad = deg2rad(Ty);
    [ThetaX_grid, ThetaY_grid] = ndgrid(Tx_rad, Ty_rad); % deg
    % Notice here thetay is the rotation around x axis!
    
    % --- Compute k magnitude ---
    Ekin_J = Ef * e; % kinetic energy in joules
    k_mag  = sqrt(2 * me * Ekin_J) / hbar; % [m^-1]
    
    % --- Map to kx, ky ---
    % ky = k_mag .* cos(ThetaY_grid) .* sin(ThetaX_grid);
    % kx = k_mag .* sin(ThetaY_grid);
    ky = k_mag .* tan(ThetaX_grid).*cos(ThetaY_grid)./sqrt(1+tan(ThetaX_grid).^2.*cos(ThetaY_grid).^2);
    kx = k_mag .* tan(ThetaY_grid).*cos(ThetaX_grid)./sqrt(1+tan(ThetaY_grid).^2.*cos(ThetaX_grid).^2);
    % x,y correspond to tilt axis and polar axis
    
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
    
    save('./matdata/fermi_surface_110mA.mat','kx','ky','I_thetax_thetay','Ef');

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
