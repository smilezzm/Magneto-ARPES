function best_p = exhaust_parameters_dual_correction(final_measured_plus, final_measured_minus,...
    k_r, standard_field, initial_guess, ax, energy)
    % angles are in unit of degrees
    
    % Extract and validate inputs
    [kx_plus, ky_plus, intensity_plus] = extract_data(final_measured_plus);
    [kx_minus, ky_minus, intensity_minus] = extract_data(final_measured_minus);
    validate_inputs_dual(kx_plus, ky_plus, intensity_plus, ...
                        kx_minus, ky_minus, intensity_minus);
    
    % Pre-compute field interpolants
    field_interpolants = create_field_interpolants(standard_field);

    % Find the most promising thetaz
    fprintf('With thetax, thetay fixed, find the most promising thetaz...')
    thetaz = linspace(initial_guess(6), initial_guess(6)+180, 100);
    score = zeros(100,1);
    for ii = 1:100
        parameters = [initial_guess(1:5); thetaz(ii); initial_guess(7)];
        score(ii) = min_delta(kx_plus, ky_plus, intensity_plus, ...
            kx_minus, ky_minus, intensity_minus, ...
            field_interpolants, k_r, energy, parameters);
    end

    % Run optimization with multiple starting points for robustness
    fprintf('Starting dual-field optimization...\n');
    

    % Return structured parameters
    parameters = struct(...
        'transX', best_p(1), 'transY', best_p(2), 'transZ', best_p(3), ...
        'thetaX', best_p(4), 'thetaY', best_p(5), 'thetaZ', best_p(6), ...
        'current', best_p(7));
    save('./matdata/optimized_p_dual_correction.mat','parameters');

    % Visualize results
    visualize_dual_results(kx_plus, ky_plus, intensity_plus, ...
                          kx_minus, ky_minus, intensity_minus, ...
                          field_interpolants, kx_ROI, ky_ROI, ...
                          best_p, history, ax1, ax2);

end

function validate_inputs_dual(kx_plus, ky_plus, intensity_plus, ...
                             kx_minus, ky_minus, intensity_minus)
    if ~(isequal(size(kx_plus), size(ky_plus), size(intensity_plus)) && ...
         isequal(size(kx_minus), size(ky_minus), size(intensity_minus)))
        error('Input arrays must have matching dimensions');
    end
end

function score = min_delta(kx_plus, ky_plus, intensity_plus, ...
            kx_minus, ky_minus, intensity_minus, ...
            field_interpolants, k_r, energy, parameters)
    parameters = [initial_guess(1:5); thetaz(ii); initial_guess(7)];
    BFcn = create_magnetic_field_function(field_interpolants, parameters);
    [Fkx_plus, Fky_plus] = inverse_mapping(BFcn, energy);
    parameters(7)=-parameters(7);
    BFcn = create_magnetic_field_function(field_interpolants, parameters);
    [Fkx_minus, Fky_minus] =inverse_mapping(BFcn, energy);
    kx_corrected_plus = Fkx_plus(kx_plus, ky_plus);
    ky_corrected_plus = Fky_plus(kx_plus, ky_plus);
    kx_corrected_minus = Fkx_minus(kx_minus, ky_minus);
    ky_corrected_minus = Fky_minus(kx_minus, ky_minus);   % (n,m) with nan

    valid_idx = ~isnan(kx_corrected_plus) & ~isnan(ky_corrected_plus);
    FI_plus = scatteredInterpolant(kx_corrected_plus(valid_idx), ky_corrected_plus(valid_idx), intensity_plus(valid_idx), 'linear', 'none');
    valid_idx = ~isnan(kx_corrected_minus) & ~isnan(ky_corrected_minus);
    FI_minus = scatteredInterpolant(kx_corrected_minus(valid_idx), ky_corrected_minus(valid_idx), intensity_minus(valid_idx), 'linear', 'none');
    
    [transX, transY] = determine_trans()
    FI_plus()
end

function visualize_dual_results(kx_plus, ky_plus, intensity_plus, ...
                          kx_minus, ky_minus, intensity_minus, ...
                          field_interpolants, kx_ROI, ky_ROI, ...
                          best_p, history, ax1, ax2)
    % Plot optimization history
    if ~isempty(history)
        plot(ax1, -history(1, :), 'LineWidth', 2);
        xlabel(ax1, 'Iteration');
        ylabel(ax1, 'Combined Similarity Score');
        title(ax1, 'Dual Optimization Progress');
        grid(ax1, 'on');
    end
    
    BFcn = create_magnetic_field_function(field_interpolants, best_p);
    [Fkx_plus, Fky_plus] = inverse_mapping(BFcn);
    best_p(7) = - best_p(7);
    BFcn = create_magnetic_field_function(field_interpolants, best_p);
    [Fkx_minus, Fky_minus] = inverse_mapping(BFcn);
    kx_corrected_plus = Fkx_plus(kx_plus, ky_plus);
    ky_corrected_plus = Fky_plus(kx_plus, ky_plus);
    kx_corrected_minus = Fkx_minus(kx_minus, ky_minus);
    ky_corrected_minus = Fky_minus(kx_minus, ky_minus);   % (n,m) with nan

    valid_idx = ~isnan(kx_corrected_plus);
    F_plus_ROI = scatteredInterpolant(kx_corrected_plus(valid_idx), ...
        ky_corrected_plus(valid_idx), ...
        intensity_plus(valid_idx), 'linear', 'none');
    intensity_plus_ROI = F_plus_ROI(kx_ROI, ky_ROI);
    valid_idx = ~isnan(kx_corrected_minus);
    F_minus_ROI = scatteredInterpolant(kx_corrected_minus(valid_idx), ...
        ky_corrected_minus(valid_idx), ...
        intensity_minus(valid_idx), 'linear', 'none');
    intensity_minus_ROI = F_minus_ROI(kx_ROI, ky_ROI);
    valid_idx = ~isnan(intensity_plus_ROI) & ~isnan(intensity_minus_ROI);
    scatter(ax2, kx_ROI(valid_idx), ky_ROI(valid_idx), 10, ...
        intensity_plus_ROI(valid_idx) - intensity_minus_ROI(valid_idx), 'filled');
    axis(ax2, 'equal');
    axis(ax2, 'tight');
    xlabel(ax2, 'k_x [10^{10} m^{-1}]');
    ylabel(ax2, 'k_y [10^{10} m^{-1}]');
    title(ax2, 'Difference between plus/minus-field fermi surface (plus-minus)');
    colormap(ax2, turbo);
    colorbar(ax2);

    % Print final parameters
    fprintf('\nOptimized Parameters:\n');
    fprintf('Translation: X=%.6f, Y=%.6f, Z=%.6f [m]\n', best_p(1), best_p(2), best_p(3));
    fprintf('Rotation: X=%.4f, Y=%.4f, Z=%.4f [rad]\n', best_p(4), best_p(5), best_p(6));
    fprintf('Current: %.6f [A]\n', -best_p(7));
end

% Include all the helper functions from the original script
function [kx, ky, intensity] = extract_data(data_struct)
    kx = data_struct.kx;
    ky = data_struct.ky;
    intensity = double(data_struct.I_thetax_thetay);
end

function k0_valid = get_valid_mask(kx_standard, ky_standard, k_mag)
    k0_valid = (kx_standard.^2 + ky_standard.^2) <= (0.4 * k_mag)^2;
end

function interpolants = create_field_interpolants(standard_field)
    interpolants = struct();
    interpolants.Bx = griddedInterpolant(standard_field.X, standard_field.Y, standard_field.Z, ...
                                        standard_field.Bx, 'linear', 'nearest');
    interpolants.By = griddedInterpolant(standard_field.X, standard_field.Y, standard_field.Z, ...
                                        standard_field.By, 'linear', 'nearest');
    interpolants.Bz = griddedInterpolant(standard_field.X, standard_field.Y, standard_field.Z, ...
                                        standard_field.Bz, 'linear', 'nearest');
end

function BFcn = create_magnetic_field_function(field_interpolants, parameters)
    % Extract parameters
    transX = parameters(1); transY = parameters(2); transZ = parameters(3);
    thetaX = parameters(4)*pi/180; thetaY = parameters(5)*pi/180; thetaZ = parameters(6)*pi/180;
    current = parameters(7);
    
    % Rotation matrices
    cx = cos(thetaX); sx = sin(thetaX);
    cy = cos(thetaY); sy = sin(thetaY);
    cz = cos(thetaZ); sz = sin(thetaZ);
    Rx = [1 0 0; 0 cx -sx; 0 sx cx];
    Ry = [cy 0 sy; 0 1 0; -sy 0 cy];
    Rz = [cz -sz 0; sz cz 0; 0 0 1];
    Rot = Rz * Rx * Ry;  % maybe Rz * Rx * Ry
    
    field_ratio = current / 0.2;
    
    % Create optimized field function
    BFcn = @(x, y, z) compute_field_at_point(x, y, z, field_interpolants, ...
                                           Rot, transX, transY, transZ, field_ratio);
end

function [Fkx, Fky] = inverse_mapping(BFcn, energy)
    % The function returns a mapping function that find the (kx,ky) map
    % when the electrons just came out at the surface of the sample, from
    % the (kx,ky) map received by ARPES. 
    % The parameters is an object that contains .transX, .transY, .transZ,
    % .thetaX, .thetaY, .thetaZ, current, and that is obtained from
    % optimized_parameters(...)
    % plot_grid = true/false to plot the grid of forward/inverse mapping
    % with provided parameters
    % Output:
    %   - Fkx: map (kx,ky) onto a kx
    %   - Fky: map (kx,ky) onto a ky
    %% Basic values
    hbar = 1.055e-34;
    me = 9.1093837015e-31; % Electron mass (kg)
    e = 1.602176634e-19;  % Elementary charge (C)
    z_target = 0.02;
    gridNum = 100;
    
    kx_i = linspace(-1.7,1.7,gridNum);   % in units 1e10 m^(-1)
    ky_i = linspace(-1.7,1.7,gridNum);   % Should be adjusted as needed
    [kx_i, ky_i] = ndgrid(kx_i, ky_i);
    k_mag = sqrt(2 * energy * e * me) / hbar / 1e10;
    k0_valid = get_valid_mask(kx_i, ky_i, k_mag);
    kx_i_valid = kx_i(k0_valid);
    ky_i_valid = ky_i(k0_valid);
    n_points = length(kx_i_valid);
    kx_f = zeros(n_points, 1);
    ky_f = zeros(n_points, 1);

    if n_points < 10
        error('few valid initial grids');
    end
    
    if exist('parfor', 'builtin') && n_points > 100
        parfor ii = 1:n_points
            vx0 = hbar * kx_i_valid(ii) * 1e10 / me;
            vy0 = hbar * ky_i_valid(ii) * 1e10 / me;
            v_final = getFinalVelocity(BFcn, energy * e, ...
                                               vx0, vy0, z_target);
            kx_f(ii) = v_final(1) * me / hbar * 1e-10;
            ky_f(ii) = v_final(2) * me / hbar * 1e-10;
        end
    else
        for ii = 1:n_points
            vx0 = hbar * kx_i_valid(ii) * 1e10 / me;
            vy0 = hbar * ky_i_valid(ii) * 1e10 / me;
            v_final = getFinalVelocity(BFcn, energy * e, ...
                                               vx0, vy0, z_target);
            kx_f(ii) = v_final(1) * me / hbar * 1e-10;
            ky_f(ii) = v_final(2) * me / hbar * 1e-10;
        end
    end

    try
        % NOTICE HERE I MAKE THE EXTRAPOLANT METHOD TO BE 'NONE'
        Fkx = scatteredInterpolant(kx_f, ky_f, kx_i_valid, 'linear', 'none');
        Fky = scatteredInterpolant(kx_f, ky_f, ky_i_valid, 'linear', 'none');
    catch ME
        error('Failed to interpolate for inverse mapping: %s', ME.message);
    end

end

function B = compute_field_at_point(x, y, z, interpolants, Rot, transX, transY, transZ, field_ratio)
    % Transform coordinates
    pos_global = [x - transX; y - transY; z - transZ];
    pos_local = Rot * pos_global;  % impose CW rotation thetaX, thetaY, thetaZ on the coil frame (during which the axes of thetaY, thetaZ are also rotating)
    % pos_global is the coordinates in lab frame, pos_local is the coordinates in coil frame
    
    % Interpolate field in local coordinates
    Bl = [interpolants.Bx(pos_local(1), pos_local(2), pos_local(3));
          interpolants.By(pos_local(1), pos_local(2), pos_local(3));
          interpolants.Bz(pos_local(1), pos_local(2), pos_local(3))];
    
    % Transform back to global coordinates
    B = Rot' * Bl * field_ratio;
end

function vxvy = getFinalVelocity(BFcn, Energy, vx0, vy0, z_target)
    me = 9.1093837015e-31; % Electron mass (kg)
    e = 1.602176634e-19;  % Elementary charge (C)
    % Adaptive time step Boris integration
    vz0 = sqrt(2 * Energy / me - vx0^2 - vy0^2);
    
    % Initial conditions
    x = [0; 0; 0];
    v_half = [vx0; vy0; vz0];
    
    % Adaptive time stepping
    B_initial = BFcn(0, 0, 0);
    omega_c = e * norm(B_initial) / me;
    dt = 2 * pi / omega_c / 200;
    
    maxSteps = 500;
    for n = 1:maxSteps
        v_prev = v_half;
        % Boris rotation
        B = BFcn(x(1), x(2), x(3));
        tvec = (- e * dt / (2 * me)) * B;
        t2 = dot(tvec, tvec);
        svec = 2 * tvec / (1 + t2);
        
        v_minus = v_half;
        v_prime = v_minus + cross(v_minus, tvec);
        v_plus = v_minus + cross(v_prime, svec);
        v_half = v_plus;
        
        % Position update
        x_new = x + dt * v_half; 
        
        % Check for target crossing
        if x_new(3) >= z_target
            % Linear interpolation for exact crossing point
            alpha = (z_target - x(3)) / (x_new(3) - x(3));
            % Interpolate velocity at exact crossing
            v_exact = v_prev + alpha * (v_half - v_prev);
            vxvy = v_exact(1:2);
            return;
        end
        
        x = x_new;
        
        % Adaptive time step based on field gradient
        if mod(n, 10) == 0
            B_norm = norm(B);
            if B_norm > 0
                omega_c = e * B_norm / me;
                dt = 2 * pi / omega_c / 200;
            end
        end
    end
    
    error('Target plane not reached within maxSteps');
end
