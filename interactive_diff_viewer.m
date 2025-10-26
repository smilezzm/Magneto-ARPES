function interactive_diff_viewer(corrected_plus, corrected_minus, k_r)
    % Interactive viewer for corrected fermi surface difference with vertical line cuts
    % 
    % Inputs:
    %   - corrected_plus: struct with .kx, .ky, .I (M×N matrices)
    %   - corrected_minus: struct with .kx, .ky, .I (M×N matrices)
    %   - k_r: radius of region of interest (optional, default 1.5)
    
    if nargin < 3
        k_r = 1.5;
    end
    
    % Extract data
    kx_plus = corrected_plus.kx;
    ky_plus = corrected_plus.ky;
    I_plus = corrected_plus.I;
    
    kx_minus = corrected_minus.kx;
    ky_minus = corrected_minus.ky;
    I_minus = corrected_minus.I;
    
    % Convert to double if needed
    if ~isa(I_plus, 'double')
        I_plus = double(I_plus);
    end
    if ~isa(I_minus, 'double')
        I_minus = double(I_minus);
    end
    
    % Apply circular mask to plus field data
    center_mask_plus = (kx_plus.^2 + ky_plus.^2 < k_r^2);
    kx_plus_ROI = kx_plus(center_mask_plus);
    ky_plus_ROI = ky_plus(center_mask_plus);
    I_plus_ROI = I_plus(center_mask_plus);
    
    % Interpolate minus field onto plus field points
    F_minus = scatteredInterpolant(kx_minus(:), ky_minus(:), I_minus(:), 'linear', 'none');
    I_minus_ROI_intp = F_minus(kx_plus_ROI, ky_plus_ROI);
    
    % Remove NaN values
    clean_mask = ~isnan(I_minus_ROI_intp) & ~isnan(I_plus_ROI);
    kx_clean = kx_plus_ROI(clean_mask);
    ky_clean = ky_plus_ROI(clean_mask);
    I_diff_clean = I_plus_ROI(clean_mask) - I_minus_ROI_intp(clean_mask);
    I_plus_clean = I_plus_ROI(clean_mask);
    I_minus_clean = I_minus_ROI_intp(clean_mask);
    
    % Create interpolants for line cuts
    F_diff = scatteredInterpolant(kx_clean, ky_clean, I_diff_clean, 'linear', 'none');
    F_plus_full = scatteredInterpolant(kx_clean, ky_clean, I_plus_clean, 'linear', 'none');
    F_minus_full = scatteredInterpolant(kx_clean, ky_clean, I_minus_clean, 'linear', 'none');
    
    % Create figure
    fig = figure('Name', 'Interactive Difference Viewer', ...
                 'NumberTitle', 'off', ...
                 'Position', [100 100 1200 500]);
    
    % Left subplot: difference map
    ax1 = subplot(1, 2, 1);
    scatter(kx_clean, ky_clean, 10, I_diff_clean, 'filled');
    axis equal tight;
    xlabel('k_x [10^{10} m^{-1}]');
    ylabel('k_y [10^{10} m^{-1}]');
    title('Difference (Plus - Minus)');
    colormap(ax1, turbo);
    cb1 = colorbar(ax1);
    ylabel(cb1, 'Intensity Difference');
    hold(ax1, 'on');
    
    % Initial vertical line
    kx_init = 0;
    ky_range = [min(ky_clean), max(ky_clean)];
    vline = plot(ax1, [kx_init kx_init], ky_range, 'w--', 'LineWidth', 2);
    
    % Right subplot: line profiles
    ax2 = subplot(1, 2, 2);
    ky_cuts = linspace(min(ky_clean), max(ky_clean), 200)';
    kx_cuts = kx_init * ones(size(ky_cuts));
    
    I_diff_cut = F_diff(kx_cuts, ky_cuts);
    I_plus_cut = F_plus_full(kx_cuts, ky_cuts);
    I_minus_cut = F_minus_full(kx_cuts, ky_cuts);
    
    % Plot initial profiles
    valid_cut = ~isnan(I_diff_cut);
    h_plus = plot(ax2, ky_cuts(valid_cut), I_plus_cut(valid_cut), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Plus field');
    hold(ax2, 'on');
    h_minus = plot(ax2, ky_cuts(valid_cut), I_minus_cut(valid_cut), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Minus field');
    xlabel(ax2, 'k_y [10^{10} m^{-1}]');
    ylabel(ax2, 'Intensity');
    title(ax2, sprintf('Vertical Cut at k_x = %.3f', kx_init));
    legend(ax2, 'Location', 'best');
    grid(ax2, 'on');
    
    % Store data in figure for callback access
    data = struct();
    data.kx_clean = kx_clean;
    data.ky_clean = ky_clean;
    data.F_diff = F_diff;
    data.F_plus = F_plus_full;
    data.F_minus = F_minus_full;
    data.vline = vline;
    data.ax2 = ax2;
    data.h_plus = h_plus;
    data.h_minus = h_minus;
    data.ky_range = ky_range;
    
    guidata(fig, data);
    
    % Set up mouse motion callback
    set(fig, 'WindowButtonMotionFcn', @(src, evt) updateLineCut(src, evt, ax1));
    
    fprintf('Move your mouse over the left plot to update the vertical cut line.\n');
end

function updateLineCut(fig, ~, ax1)
    % Get current point in axes coordinates
    cp = get(ax1, 'CurrentPoint');
    kx_current = cp(1, 1);
    ky_current = cp(1, 2);
    
    % Check if cursor is within axes
    xlims = xlim(ax1);
    ylims = ylim(ax1);
    if kx_current < xlims(1) || kx_current > xlims(2) || ...
       ky_current < ylims(1) || ky_current > ylims(2)
        return;
    end
    
    % Retrieve stored data
    data = guidata(fig);
    
    % Update vertical line position
    set(data.vline, 'XData', [kx_current kx_current]);
    
    % Generate new cut data
    ky_cuts = linspace(min(data.ky_clean), max(data.ky_clean), 200)';
    kx_cuts = kx_current * ones(size(ky_cuts));
    
    I_diff_cut = data.F_diff(kx_cuts, ky_cuts);
    I_plus_cut = data.F_plus(kx_cuts, ky_cuts);
    I_minus_cut = data.F_minus(kx_cuts, ky_cuts);
    
    % Update line plots
    valid_cut = ~isnan(I_diff_cut);
    if any(valid_cut)
        set(data.h_plus, 'XData', ky_cuts(valid_cut), 'YData', I_plus_cut(valid_cut));
        set(data.h_minus, 'XData', ky_cuts(valid_cut), 'YData', I_minus_cut(valid_cut));
        
        % Update title
        title(data.ax2, sprintf('Vertical Cut at k_x = %.3f', kx_current));
        
        % Auto-scale y-axis for better viewing
        if any(~isnan(I_plus_cut(valid_cut))) || any(~isnan(I_minus_cut(valid_cut)))
            all_vals = [I_plus_cut(valid_cut); I_minus_cut(valid_cut)];
            all_vals = all_vals(~isnan(all_vals));
            if ~isempty(all_vals)
                y_range = [min(all_vals), max(all_vals)];
                y_margin = 0.1 * (y_range(2) - y_range(1));
                if y_margin > 0
                    ylim(data.ax2, [y_range(1) - y_margin, y_range(2) + y_margin]);
                end
            end
        end
    end
    
    drawnow limitrate;
end
