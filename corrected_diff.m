function corrected_diff(k_r, corrected_plus, corrected_minus,further_transform)
    % k_r=1.4
    % further_transform: true or false
    intersection = 0.68; % cut a line to see the difference along it
    kx_plus = corrected_plus.kx;
    ky_plus = corrected_plus.ky;
    try
        I_plus = corrected_plus.I;
    catch
        I_plus = corrected_plus.I_thetax_thetay;
    end
    kx_minus = corrected_minus.kx;
    ky_minus = corrected_minus.ky;
    try
        I_minus = corrected_minus.I;
    catch 
        I_minus = corrected_minus.I_thetax_thetay;
    end
    if ~isa(I_plus, 'double') 
        I_plus = double(I_plus);
    end
    if ~isa(I_minus, 'double')
        I_minus = double(I_minus);
    end

    if ~further_transform
        center_mask_plus = (kx_plus.^2 + ky_plus.^2 < k_r^2);
        kx_plus_ROI = kx_plus(center_mask_plus);
        ky_plus_ROI = ky_plus(center_mask_plus);
        I_plus_ROI = I_plus(center_mask_plus);
        
        F_minus = scatteredInterpolant(kx_minus(:),ky_minus(:),I_minus(:),'linear','none');
        I_minus_ROI_intp = F_minus(kx_plus_ROI, ky_plus_ROI);
        clean_mask = ~isnan(I_minus_ROI_intp) & ~isnan(I_plus_ROI);
        figure
        scatter(kx_plus_ROI(clean_mask), ky_plus_ROI(clean_mask), 10, I_plus_ROI(clean_mask)-I_minus_ROI_intp(clean_mask),'filled');
        axis equal tight;
        xlabel('k_x [10^{10} m^{-1}]');
        ylabel('k_y [10^{10} m^{-1}]');
        title('Difference between corrected +field and -field fermi surface');
        colormap(turbo); colorbar;
    else
        % Interpolate scattered data onto a grid
        [xq,yq] = meshgrid(linspace(-k_r,k_r,100), linspace(-k_r,k_r,100));
        I_grid_plus = griddata(kx_plus, ky_plus, I_plus, xq, yq, 'linear');
        I_grid_minus = griddata(kx_minus, ky_minus, I_minus, xq, yq, 'linear');
        mask_plus = isnan(I_grid_plus);
        mask_minus = isnan(I_grid_minus);
        I_grid_plus(mask_plus) = 0.0;
        I_grid_minus(mask_minus) = 0.0;
        
        % Convert to images
        I_img_plus = mat2gray(I_grid_plus);
        I_img_minus = mat2gray(I_grid_minus);
        
        [optimizerConfig, metricConfig] = imregconfig('monomodal');
        % Use imregtform (requires Image Processing Toolbox)
        tform0 = affine2d([cosd(30) -sind(30) 0; sind(30) cosd(30) 0; 0.15 0 1]); % ~30° CCW rotation + shift
        tform = imregtform(I_img_minus, I_img_plus, 'rigid', optimizerConfig, metricConfig, ...
                   'InitialTransformation', tform0);
        
        % Apply transform
        I_img_minus_registered = imwarp(I_img_minus, tform, 'OutputView', imref2d(size(I_img_plus)));
        I_diff = I_img_plus - I_img_minus_registered;
        
        % plot
        figure;
        subplot(1,3,1);
        imagesc(xq(1,:), yq(:,1), I_img_plus); axis xy equal; title('Plus image'); colorbar;
        xline(intersection, 'w--', 'LineWidth', 1.5);   % white dashed line
        
        subplot(1,3,2);
        imagesc(xq(1,:), yq(:,1), I_img_minus_registered); axis xy equal; title('Minus registered'); colorbar;
        xline(intersection, 'w--', 'LineWidth', 1.5);   % white dashed line
        
        subplot(1,3,3);
        imagesc(xq(1,:), yq(:,1), I_diff); axis xy equal; title('Difference'); colorbar;
        xline(intersection, 'w--', 'LineWidth', 1.5);   % white dashed line
        
        % Find nearest column index to intersection
        [~, idx] = min(abs(xq(1,:) - intersection));
        
        % Define a window of columns around idx
        halfwin = 2;  % number of pixels on each side
        cols = max(1, idx-halfwin) : min(size(I_img_plus,2), idx+halfwin);
        
        % Average across those columns
        I_line_plus  = mean(I_img_plus(:, cols), 2, 'omitnan');
        I_line_minus = mean(I_img_minus(:, cols), 2, 'omitnan');
        
        % ky axis
        ky_line = yq(:,1);
        
        % Plot
        figure;
        plot(ky_line, I_line_plus, 'b-', 'LineWidth', 1.5);
        hold on;
        plot(ky_line, I_line_minus, 'r-', 'LineWidth', 1.5);
        xlabel('$k_y \, (\AA^{-1})$','Interpreter','latex');
        ylabel(sprintf('Intensity (avg over %d px) at k_x = %.2f', length(cols), intersection));
        title('Line cut of the corrected Fermi surface');
        legend('+ field', '- field');

    end
        
end