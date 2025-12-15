% function build_h5(original_h5path, output_h5path, measured_bands, standard_field, parameters)
%     % INPUT
%     % Don't forget to modify the offsets and delta in the function
%     % original_h5path - path of .h5 from which measured_bands are obtained
%     % output_h5path - path of .h5 to store the corrected bands
%     % measured_bands - the struct produced by loading "3D_bands_110mA.mat"
%     %   .kx: [W, H, D]double
%     %   .ky: [W, H, D]double
%     %   .I_thetax_thetay: [W, H, D]single
%     %   .E: [1, W]
%     %      W = params.axes0;   % Energy dimensions 751
%     %      H = params.axes1;   % ThetaX dimensions [-18,18] 876
%     %      D = params.axes2;   % ThetaY dimensions [-14,14] 113
%     %      N = W * H * D;
%     %
%     % standard_field - the struct produced by loading "standard_field.mat"
%     %   describes the field without translation/rotation
%     %   .X .Y .Z - 3D matrix containing grid position
%     %   .Bx .By .Bz - 3D matrix containing field strength
%     %   .current - current used to simulate the field
%     %   .target - highest limit of the position grid
%     % 
%     % parameters - the struct containing
%     %   .transX .transY .transZ .thetaX .thetaY .thetaZ .current
%     %
%     % OUTPUT
%     % .h5 file
%     %   count, offset, and delta of axes0, axes1, and axes2 are all the same
%     %   only the dataset /Data/Count in the original .h5 is changed
%     %   size (D,H,W)
% 
%     %% Constant
%     me   = 9.10938356e-31;     % electron mass [kg]
%     hbar = 1.054571817e-34;    % reduced Planck constant [J*s]
%     e   = 1.602176634e-19;    % 1 eV in joules
%     energy = measured_bands.E;
%     parameters.thetaX = parameters.thetaX / 180 * pi;
%     parameters.thetaY = parameters.thetaY / 180 * pi;
%     parameters.thetaZ = parameters.thetaZ / 180 * pi;
% 
%     %% get the grid of the original .h5 file
%     info = h5info(original_h5path);
% 
%     params.thetax_offset = info.Groups(2).Groups(2).Attributes(1).Value;  % axes1_offset (ThetaX [deg]) (-17.9795)
%     params.thetax_delta  = info.Groups(2).Groups(2).Attributes(2).Value;   % axes1_delta (0.0411)
%     params.thetay_offset = info.Groups(2).Groups(3).Attributes(1).Value;  % axes2_offset (Thetay [deg]) (-14)
%     params.thetay_delta  = info.Groups(2).Groups(3).Attributes(2).Value;    % axes2_delta (0.25)
%     H = info.Groups(2).Groups(2).Attributes(3).Value;   % ThetaX dimensions [-18,18]  876
%     D = info.Groups(2).Groups(3).Attributes(3).Value;   % ThetaY dimensions [-14,14]  113
%     W = info.Groups(2).Groups(1).Attributes(3).Value;   % Energy dimensions [-109.235, 110.735]  751
%     Tx = params.thetax_offset + (0:double(H)-1) * params.thetax_delta;
%     Ty = params.thetay_offset + (0:double(D)-1) * params.thetay_delta;
%     [Tx, Ty] = ndgrid(Tx, Ty);
% 
%     %% write the corrected data
%     % Don't preallocate the entire array 
%     copyfile(original_h5path, output_h5path);
% 
%     % Process and write one energy slice at a time
%     for ii = 1:W    % energy is longer than W, we discard the outlying part
%         measured_band.kx = squeeze(measured_bands.kx(ii, :, :));
%         measured_band.ky = squeeze(measured_bands.ky(ii, :, :));
%         measured_band.I_thetax_thetay = squeeze(measured_bands.I_thetax_thetay(ii, :, :));
%         [Fkx, Fky] = calc_inverse_mapping(standard_field, energy(ii), parameters, false);
%         corrected_band = backward_mapping(measured_band, Fkx, Fky, false);
% 
%         k_mag  = sqrt(2 * me * e * energy(ii)) / hbar * 1e-10;
%         thetax = asind(corrected_band.ky ./ sqrt(k_mag^2 - corrected_band.kx.^2));
%         thetay = asind(corrected_band.kx ./ sqrt(k_mag^2 - corrected_band.ky.^2));
% 
%         Ftheta = scatteredInterpolant(thetax, thetay, double(corrected_band.I), 'linear', 'none');
%         new_data = single(Ftheta(Tx, Ty));
%         nan_idx = isnan(new_data);
%         new_data(nan_idx) = 0;
% 
%         % Write directly to HDF5 file (one slice at a time)
%         h5write(output_h5path, '/Data/Count', permute(new_data, [2,1]), double([1, 1, ii]), double([D, H, 1]));
% 
%         % Optional: display progress
%         if mod(ii, 50) == 0
%             fprintf('Processed %d/%d energy slices\n', ii, length(energy));
%         end
%     end
% 
%     %% Set the manipulator offset 0
%     h5writeatt(output_h5path, '/Manipulator', 'T', 0);
%     h5writeatt(output_h5path, '/Manipulator', 'F', 0);
%     h5writeatt(output_h5path, '/Manipulator', 'A', 0);
% end


function build_bin(output_binpath, measured_bands, standard_field, parameters, ...
    thetaxOffset, thetaxDelta, thetaxNum, thetayOffset, thetayDelta, thetayNum)
    % INPUT
    % output_binpath - path of .bin file to store the corrected bands
    % measured_bands - the struct produced by loading "3D_bands_110mA.mat"
    %   [W,H,D] 799,361,281
    % standard_field - the struct produced by loading "standard_field.mat"
    % parameters - the struct containing transformation parameters 
    %   angles in degree !!!
    %
    % OUTPUT
    % .bin file containing the 3D intensity data in single precision
    % Data layout: [D, H, W] = [ThetaY, ThetaX, Energy]
    % Separate .mat file with metadata (grid information, dimensions, etc.)

    %% Constant
    me   = 9.10938356e-31;     % electron mass [kg]
    hbar = 1.054571817e-34;    % reduced Planck constant [J*s]
    e   = 1.602176634e-19;    % 1 eV in joules
    energy = measured_bands.E;
    parameters.thetaX = parameters.thetaX / 180 * pi;
    parameters.thetaY = parameters.thetaY / 180 * pi;
    parameters.thetaZ = parameters.thetaZ / 180 * pi;

    %% Define output grid
    % params.thetax_offset = -17.9795;  % ThetaX [deg]
    % params.thetax_delta  = 0.1;
    % params.thetay_offset = -14;       % ThetaY [deg]
    % params.thetay_delta  = 0.1;
    H = thetaxNum;  % 361;   % ThetaX dimensions [-18,18]
    D = thetayNum;  % 281;   % ThetaY dimensions [-14,14]
    Tx = thetaxOffset + (0:double(H)-1) * thetaxDelta;
    Ty = thetayOffset + (0:double(D)-1) * thetayDelta;
    [Tx, Ty] = ndgrid(Tx, Ty);

    %% Open binary file for writing
    fid = fopen(output_binpath, 'wb');  % 'wb' = write binary
    if fid == -1
        error('Cannot open file for writing: %s', output_binpath);
    end

    % Process and write one energy slice at a time
    fprintf('Processing %d energy slices...\n', length(energy));
    for ii = 1:length(energy)
        measured_band.kx = squeeze(measured_bands.kx(ii, :, :));
        measured_band.ky = squeeze(measured_bands.ky(ii, :, :));
        measured_band.I_thetax_thetay = squeeze(measured_bands.I_thetax_thetay(ii, :, :));
        [Fkx, Fky] = calc_inverse_mapping(standard_field, energy(ii), parameters, false);
        corrected_band = backward_mapping(measured_band, Fkx, Fky, false);
        
        k_mag  = sqrt(2 * me * e * energy(ii)) / hbar * 1e-10;
        thetax = asind(corrected_band.ky ./ sqrt(k_mag^2 - corrected_band.kx.^2));
        thetay = asind(corrected_band.kx ./ sqrt(k_mag^2 - corrected_band.ky.^2));

        Ftheta = scatteredInterpolant(thetax(:), thetay(:), double(corrected_band.I(:)), 'linear', 'none');
        new_data = single(Ftheta(Tx, Ty));
        nan_idx = isnan(new_data);
        new_data(nan_idx) = 0;

        % Write slice to binary file (column-major order)
        % Transpose to [D, H] format
        count = fwrite(fid, permute(new_data, [2,1]), 'single');
        if count ~= D * H
            error('Failed to write complete data slice %d', ii);
        end

        % Clear variables to save memory
        clear measured_band corrected_band Ftheta new_data nan_idx thetax thetay

        % Display progress
        if mod(ii, 50) == 0
            fprintf('Processed %d/%d energy slices\n', ii, length(energy));
        end
    end

    % Close the file
    fclose(fid);
    fprintf('Binary file saved: %s\n', output_binpath);

    %% Save metadata to companion .mat file
    metadata.dimensions = [D, H, length(energy)];  % [ThetaY, ThetaX, Energy]
    metadata.thetax_offset = thetaxOffset;
    metadata.thetax_delta = thetaxDelta;
    metadata.thetax_count = H;
    metadata.thetay_offset = thetayOffset;
    metadata.thetay_delta = thetayDelta;
    metadata.thetay_count = D;
    metadata.energy_offset = energy(1);
    metadata.energy_delta = energy(2) - energy(1);
    metadata.energy_count = length(energy);
    metadata.datatype = 'single';
    metadata.description = 'Binary data layout: [ThetaY, ThetaX, Energy] in column-major order';

    [filepath, name, ~] = fileparts(output_binpath);
    metadata_path = fullfile(filepath, name+"_metadata.mat");
    save(metadata_path, 'metadata');
    fprintf('Metadata saved: %s\n', metadata_path);
end