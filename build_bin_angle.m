function build_bin_angle(output_binpath, measured_bands, standard_field, parameters, ...
    thetaxOffset, thetaxDelta, thetaxNum, thetayOffset, thetayDelta, thetayNum, Energy)
    % build the corrected data as .bin, with angle mapping instead of k
    % mapping
    % measured_bands contains .thetax, .thetay, which are of size (H,D)

    %% Prepare
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

    %% Map the angle to angle, this is approximately the same for all energy slices
    [Ftx, Fty] = calc_inverse_mapping_angle(standard_field, Energy, parameters, false);
    Tx_corrected = Ftx(Tx, Ty);   % (H,D) with nans
    Ty_corrected = Fty(Tx, Ty);
    valid_idx = ~isnan(Tx_corrected) & ~isnan(Ty_corrected);

    %% Open binary file for writing
    fid = fopen(output_binpath, 'wb');  % 'wb' = write binary
    if fid == -1
        error('Cannot open file for writing: %s', output_binpath);
    end


    % Process and write one energy slice at a time
    fprintf('Processing %d energy slices...\n', length(energy));
    for ii = 1:length(energy)
        % Directly map between angles instead of momentum
        I_thetax_thetay = squeeze(measured_bands.I_thetax_thetay(ii, :, :));
        Ftheta = scatteredInterpolant(Tx_corrected(valid_idx), Ty_corrected(valid_idx), double(I_thetax_thetay(valid_idx)), 'linear', 'none');
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
        clear Ftheta new_data nan_idx

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