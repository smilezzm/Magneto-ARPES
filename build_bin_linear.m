function build_bin_linear(movingPoints, fixedPoints, measured_angle, outPath, doPlot, E_fs)
    % INPUT
    % - movingPoints: the high-symmetry points in measured_angle 
    %       [x1 y1; x2 y2; x3 y3; x4 y4; x5 y5; x6 y6; x7 y7]
    % - fixedPoints: the corresponding high-symmetry points in 0-field band (angle coordinates)
    % - measured_angle: the measured band on a single energy-slice
    %       containing .thetax, .thetay, .I, .E
    %       .thetax, .thetay are (H,D), .I is (W, H, D), .E is (W,1)
    % - doPlot: plot the corrected pattern in angle space
    % - outPath: the file to store the corrected data as .bin
    % - E_fs: (optional, required if doPlot=true) energy value for Fermi surface plot
    % OUTPUT
    % - corrected: .thetax, .thetay with size of (M,N), .I with size of
    %       (W,H,D), .E with the size of (W,1)

    %% Validate inputs
    if doPlot && (nargin < 6 || isempty(E_fs))
        error('E_fs must be provided when doPlot is true');
    end

    %% linear transormation
    tform = fitgeotrans(movingPoints, fixedPoints, 'affine');
    [thetax_new, thetay_new] = transformPointsForward(tform, measured_angle.thetax, measured_angle.thetay);

    %% Open binary file for writing
    fid = fopen(outPath, 'wb');  % 'wb' = write binary
    if fid == -1
        error('Cannot open file for writing: %s', output_binpath);
    end
    [H,D] = size(measured_angle.thetax);
    energy = measured_angle.E;
    
    % If plotting, store the I_new for the slice closest to E_fs
    if doPlot
        [~, fs_idx] = min(abs(energy - E_fs));
        I_new_fs = [];
    end

    %% Writing slice by slice
    for ii = 1:length(energy)
        I_slice = double(squeeze(measured_angle.I_thetax_thetay(ii,:,:)));
        Ftheta = scatteredInterpolant(thetax_new(:), thetay_new(:), I_slice(:), 'linear', 'none');
        I_new = Ftheta(measured_angle.thetax, measured_angle.thetay);   % (H,D)
        idx = isnan(I_new);
        I_new(idx) = 0;
        I_new = single(I_new);
        
        % Store the slice for plotting if this is the Fermi surface energy
        if doPlot && ii == fs_idx
            I_new_fs = I_new;
        end
        
        count = fwrite(fid, permute(I_new, [2,1]), 'single');
        if count ~= D * H
            error('Failed to write complete data slice %d', ii);
        end
    end

    fclose(fid);
    fprintf('Binary file saved: %s\n', outPath);

    %% Save metadata to companion .mat file
    metadata.dimensions = [D, H, length(energy)];  % [ThetaY, ThetaX, Energy]
    metadata.thetax_offset = measured_angle.thetax(1,1);
    metadata.thetax_delta = measured_angle.thetax(2,1)-measured_angle.thetax(1,1);
    metadata.thetax_count = H;
    metadata.thetay_offset = measured_angle.thetay(1,1);
    metadata.thetay_delta = measured_angle.thetay(1,2)-measured_angle.thetay(1,1);
    metadata.thetay_count = D;
    metadata.energy_offset = energy(1);
    metadata.energy_delta = energy(2) - energy(1);
    metadata.energy_count = length(energy);
    metadata.datatype = 'single';
    metadata.description = 'Binary data layout: [ThetaY, ThetaX, Energy] in column-major order';

    [filepath, name, ~] = fileparts(outPath);
    metadata_path = fullfile(filepath, name+"_metadata.mat");
    save(metadata_path, 'metadata');
    fprintf('Metadata saved: %s\n', metadata_path);

    if doPlot
        figure('Color', 'w');
        imagesc(measured_angle.thetax(1,:), measured_angle.thetay(:,1), I_new_fs');
        axis xy equal tight;
        colorbar;
        xlabel('ThetaX (deg)');
        ylabel('ThetaY (deg)');
        title(sprintf('Corrected Intensity at E = %.3f eV ', ...
                      energy(fs_idx)));
        colormap(hot);
    end
end