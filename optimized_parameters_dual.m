function parameters = optimized_parameters_dual(standard_measured, final_measured_plus, final_measured_minus,...
    k_r, standard_field, initial_guess, options)
    % This function returns the parameters that make the forward simulated
    % fermi surface figures match the best with both measured figures, for the
    % plus field and minus field
    %
    % Inputs:
    %   - standard_measured: the measured fermi surface of the standard sample with field off
    %       .kx, .ky, .I_thetax_thetay
    %   - final_measured_plus: the measured fermi surface with positive field
    %       .kx, .ky, .I_thetax_thetay
    %   - final_measured_minus: the measured fermi surface with negative field
    %       .kx, .ky, .I_thetax_thetay
    %   - k_r: the radius of k (in 1e10 m^(-1)) that we are interested at the center 
    %   - standard_field: an object containing the standard field values on grids
    %       .Bx, .By, .Bz, .X, .Y, .Z
    %   - initial_guess: a vector [transX; transY; transZ; thetaX; thetaY; thetaZ; current]
    %   - options: (optional) optimization options structure
    % Outputs:
    %   - parameters: optimized parameters structure
    
    if nargin < 7
        options = struct();
    end
    
    % Constants
    CONSTANTS = struct(...
        'hbar', 1.055e-34, ...
        'me', 9.1093837015e-31, ...
        'e', 1.602176634e-19, ...
        'target_E', 110.56, ...
        'z_target', 20e-3);
    
    % Extract and validate inputs
    [kx_standard, ky_standard, intensity_standard] = extract_data(standard_measured);
    [kx_plus, ky_plus, intensity_plus] = extract_data(final_measured_plus);
    [kx_minus, ky_minus, intensity_minus] = extract_data(final_measured_minus);
    validate_inputs_dual(kx_standard, ky_standard, intensity_standard, ...
                        kx_plus, ky_plus, intensity_plus, ...
                        kx_minus, ky_minus, intensity_minus);
    kx_ROI = linspace(-k_r,k_r,70);
    ky_ROI = linspace(-k_r,k_r,70);
    [kx_ROI,ky_ROI]=ndgrid(kx_ROI,ky_ROI);
    F_plus = scatteredInterpolant(kx_plus(:), ky_plus(:), intensity_plus(:), 'linear', 'none');
    F_minus = scatteredInterpolant(kx_minus(:), ky_minus(:), intensity_minus(:), 'linear', 'none');
    intensity_plus_ROI = F_plus(kx_ROI, ky_ROI);
    intensity_minus_ROI = F_minus(kx_ROI, ky_ROI);

    % Pre-compute k magnitude and valid mask
    k_mag = sqrt(2 * CONSTANTS.target_E * CONSTANTS.e * CONSTANTS.me) / CONSTANTS.hbar / 1e10;
    k0_valid = get_valid_mask(kx_standard, ky_standard, k_mag);
    
    % Pre-compute field interpolants
    field_interpolants = create_field_interpolants(standard_field);
    
    % Setup optimization
    history = [];
    objfun = @(p) -calc_dual_similarity_optimized(p, kx_standard, ky_standard, intensity_standard,...
        kx_ROI, ky_ROI, intensity_plus_ROI, intensity_minus_ROI, ...
        field_interpolants, CONSTANTS, k0_valid);
    
    % Set optimization options
    opt_options = get_optimization_options(options, @outfun);
    lb = [-1e-3, -1e-3, -1e-3, -0.1*pi, -0.1*pi, -pi, 0.068];
    ub = [1e-3, 1e-3, 1e-3, 0.1*pi, 0.1*pi, pi, 0.085];
    
    % Run optimization with multiple starting points for robustness
    fprintf('Starting dual-field optimization...\n');
    best_p = run_multi_start_optimization(objfun, initial_guess, lb, ub, opt_options);
    
    % Visualize results
    visualize_dual_results(kx_standard, ky_standard, intensity_standard, ...
                          kx_plus, ky_plus, intensity_plus, ...
                          kx_minus, ky_minus, intensity_minus, ...
                          field_interpolants, CONSTANTS, k0_valid, best_p, history);
    
    % Return structured parameters
    parameters = struct(...
        'transX', best_p(1), 'transY', best_p(2), 'transZ', best_p(3), ...
        'thetaX', best_p(4), 'thetaY', best_p(5), 'thetaZ', best_p(6), ...
        'current', best_p(7));
    save('optimized_p_dual.mat','parameters');
    
    function stop = outfun(p, optimValues, state)
        stop = false;
        if isequal(state, 'iter')
            history = [history, [optimValues.fval; p(:)]];
            if mod(optimValues.iteration, 10) == 0
                fprintf('Iteration %d: f = %.6f\n', optimValues.iteration, optimValues.fval);
            end
        end
    end
end

function validate_inputs_dual(kx_standard, ky_standard, intensity_standard, ...
                             kx_plus, ky_plus, intensity_plus, ...
                             kx_minus, ky_minus, intensity_minus)
    if ~(isequal(size(kx_standard), size(ky_standard), size(intensity_standard)) && ...
         isequal(size(kx_plus), size(ky_plus), size(intensity_plus)) && ...
         isequal(size(kx_minus), size(ky_minus), size(intensity_minus)))
        error('Input arrays must have matching dimensions');
    end
end

function score = calc_dual_similarity_optimized(parameters, kx_standard, ky_standard, intensity_standard,...
    kx_ROI, ky_ROI, intensity_plus_ROI, intensity_minus_ROI, ...
    field_interpolants, CONSTANTS, k0_valid)
    
    % Forward simulation with positive current
    [kx_sim_plus, ky_sim_plus, I_sim_plus] = forward_sim_vectorized(kx_standard, ky_standard, intensity_standard,...
        field_interpolants, CONSTANTS, k0_valid, parameters);
    
    % Forward simulation with negative current
    parameters_minus = parameters;
    parameters_minus(7) = -parameters(7);  % Flip current sign
    [kx_sim_minus, ky_sim_minus, I_sim_minus] = forward_sim_vectorized(kx_standard, ky_standard, intensity_standard,...
        field_interpolants, CONSTANTS, k0_valid, parameters_minus);
    
    % Calculate similarity for positive field
    score_plus = calc_single_similarity(kx_sim_plus, ky_sim_plus, I_sim_plus, ...
                                       kx_ROI, ky_ROI, intensity_plus_ROI);
    
    % Calculate similarity for negative field
    score_minus = calc_single_similarity(kx_sim_minus, ky_sim_minus, I_sim_minus, ...
                                        kx_ROI, ky_ROI, intensity_minus_ROI);
    
    % Combined score 
    score = score_plus * score_minus;
    
    % Debug output every 50 evaluations
    persistent eval_count;
    if isempty(eval_count)
        eval_count = 0;
    end
    eval_count = eval_count + 1;
    if mod(eval_count, 50) == 0
        fprintf('  Score breakdown: Plus=%.4f, Minus=%.4f, Combined=%.4f\n', ...
                score_plus, score_minus, score);
    end
end

function score = calc_single_similarity(kx_sim, ky_sim, I_sim, kx_ROI, ky_ROI, intensity_measured_ROI)
    % Handle case with insufficient simulation data
    if numel(kx_sim) < 3
        score = -Inf;
        return;
    end
    
    try
        % Interpolate simulation results to measured grid
        F_ROI = scatteredInterpolant(kx_sim, ky_sim, I_sim, 'linear', 'none');
        intensity_sim_ROI = F_ROI(kx_ROI, ky_ROI);

        % Remove NaN values
        valid_idx = ~isnan(intensity_sim_ROI) & ~isnan(intensity_measured_ROI);
        if sum(valid_idx(:)) < 0.8 * numel(intensity_sim_ROI)
            score = -Inf;
            return;
        end
        intensity_measured_ROI(~valid_idx)=0;
        intensity_sim_ROI(~valid_idx)=0;
        
        % Calculate normalized correlation
        if std(intensity_measured_ROI(:)) == 0 || std(intensity_sim_ROI(:)) == 0
            score = -Inf;
            return;
        end

        norm_measured = (intensity_measured_ROI - mean(intensity_measured_ROI(:))) / std(intensity_measured_ROI(:));
        norm_sim = (intensity_sim_ROI - mean(intensity_sim_ROI(:))) / std(intensity_sim_ROI(:));
        [~, ssim_map] = ssim(norm_sim, norm_measured);
        score = mean(ssim_map(valid_idx));
    catch
        score = -Inf;
    end
end

function visualize_dual_results(kx_standard, ky_standard, intensity_standard, ...
                               kx_plus, ky_plus, intensity_plus, ...
                               kx_minus, ky_minus, intensity_minus, ...
                               field_interpolants, CONSTANTS, k0_valid, best_p, history)
    % Plot optimization history
    if ~isempty(history)
        figure('Name', 'Dual Optimization Progress');
        plot(-history(1, :), 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Combined Similarity Score');
        title('Dual Optimization Progress');
        grid on;
    end
    
    % Generate final simulations for visualization
    [kx_sim_plus, ky_sim_plus, I_sim_plus] = forward_sim_vectorized(kx_standard, ky_standard, intensity_standard,...
        field_interpolants, CONSTANTS, k0_valid, best_p);
    
    best_p_minus = best_p;
    best_p_minus(7) = -best_p(7);
    [kx_sim_minus, ky_sim_minus, I_sim_minus] = forward_sim_vectorized(kx_standard, ky_standard, intensity_standard,...
        field_interpolants, CONSTANTS, k0_valid, best_p_minus);
    
    % Standard (no field)
    figure
    scatter(kx_standard(:), ky_standard(:), 10, intensity_standard(:), 'filled');
    axis equal tight;
    xlabel('k_x [10^{10} m^{-1}]');
    ylabel('k_y [10^{10} m^{-1}]');
    title('Standard (No Field)');
    colormap(turbo); colorbar;
    
    % % Measured plus field
    % subplot(2, 3, 2);
    % scatter(kx_plus(:), ky_plus(:), 10, intensity_plus(:), 'filled');
    % axis equal tight;
    % xlabel('k_x [10^{10} m^{-1}]');
    % ylabel('k_y [10^{10} m^{-1}]');
    % title(sprintf('Measured (+Field) @ E = %.3f eV', CONSTANTS.target_E));
    % colormap(turbo); colorbar;
    % 
    % % Simulated plus field
    % subplot(2, 3, 3);
    % if ~isempty(kx_sim_plus)
    %     scatter(kx_sim_plus, ky_sim_plus, 10, I_sim_plus, 'filled');
    % end
    % axis equal tight;
    % xlabel('k_x [10^{10} m^{-1}]');
    % ylabel('k_y [10^{10} m^{-1}]');
    % title(sprintf('Simulated (+Field) @ E = %.3f eV', CONSTANTS.target_E));
    % colormap(turbo); colorbar;
    
    % Create comparison plot
    figure
    % Here directly subtract them two
    k_r=1.4;
    center_mask_plus = (kx_plus.^2 + ky_plus.^2 < k_r^2);
    kx_plus_ROI = kx_plus(center_mask_plus);
    ky_plus_ROI = ky_plus(center_mask_plus);
    I_plus_ROI = intensity_plus(center_mask_plus);

    F_I = scatteredInterpolant(kx_sim_plus,ky_sim_plus,I_sim_plus,'linear','none');
    I_sim_plus_intp = F_I(kx_plus_ROI, ky_plus_ROI);
    clean_mask = ~isnan(I_sim_plus_intp) & ~isnan(I_plus_ROI);
    subplot(1,2,1)
    scatter(kx_plus_ROI(clean_mask), ky_plus_ROI(clean_mask), 10, I_plus_ROI(clean_mask)-I_sim_plus_intp(clean_mask),'filled');
    axis equal tight;
    xlabel('k_x [10^{10} m^{-1}]');
    ylabel('k_y [10^{10} m^{-1}]');
    title('Difference between measured and simulated +field fermi surface (measured-simulated)');
    colormap(turbo); colorbar;
    annotation('textbox', [0.45 0.3 0.1 0.4], 'String', ...
    sprintf(['The parameters are:\n' ...
             'transX = %.5f\ntransY = %.5f\ntransZ = %.5f\n' ...
             'thetaX = %.4f\nthetaY = %.4f\nthetaZ = %.4f\n' ...
             'current = %.3f'], ...
             best_p(1), best_p(2), best_p(3), ...
             best_p(4), best_p(5), best_p(6), ...
             best_p(7)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 12);

    % % Measured minus field
    % subplot(2, 3, 5);
    % scatter(kx_minus(:), ky_minus(:), 10, intensity_minus(:), 'filled');
    % axis equal tight;
    % xlabel('k_x [10^{10} m^{-1}]');
    % ylabel('k_y [10^{10} m^{-1}]');
    % title(sprintf('Measured (-Field) @ E = %.3f eV', CONSTANTS.target_E));
    % colormap(turbo); colorbar;
    % 
    % % Simulated minus field
    % subplot(2, 3, 6);
    % if ~isempty(kx_sim_minus)
    %     scatter(kx_sim_minus, ky_sim_minus, 10, I_sim_minus, 'filled');
    % end
    % axis equal tight;
    % xlabel('k_x [10^{10} m^{-1}]');
    % ylabel('k_y [10^{10} m^{-1}]');
    % title(sprintf('Simulated (-Field) @ E = %.3f eV', CONSTANTS.target_E));
    % colormap(turbo); colorbar;
    
    % Also, directly plot the subtraction
    center_mask_minus = (kx_minus.^2 + ky_minus.^2 < k_r^2);
    kx_minus_ROI = kx_minus(center_mask_minus);
    ky_minus_ROI = ky_minus(center_mask_minus);
    I_minus_ROI = intensity_minus(center_mask_minus);

    F_I = scatteredInterpolant(kx_sim_minus,ky_sim_minus,I_sim_minus,'linear','none');
    I_sim_minus_intp = F_I(kx_minus_ROI, ky_minus_ROI);
    clean_mask = ~isnan(I_sim_minus_intp) & ~isnan(I_minus_ROI);
    subplot(1,2,2)
    scatter(kx_minus_ROI(clean_mask), ky_minus_ROI(clean_mask), 10, I_minus_ROI(clean_mask)-I_sim_minus_intp(clean_mask),'filled');
    axis equal tight;
    xlabel('k_x [10^{10} m^{-1}]');
    ylabel('k_y [10^{10} m^{-1}]');
    title('Difference between measured and simulated -field fermi surface (measured-simulated)');
    colormap(turbo); colorbar;

    % Print final parameters
    fprintf('\nOptimized Parameters:\n');
    fprintf('Translation: X=%.6f, Y=%.6f, Z=%.6f [m]\n', best_p(1), best_p(2), best_p(3));
    fprintf('Rotation: X=%.4f, Y=%.4f, Z=%.4f [rad]\n', best_p(4), best_p(5), best_p(6));
    fprintf('Current: %.6f [A]\n', best_p(7));
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

function [kx_sim, ky_sim, I_sim] = forward_sim_vectorized(kx_standard, ky_standard, intensity_standard,...
    field_interpolants, CONSTANTS, k0_valid, parameters)
    
    % Extract valid points
    kx_valid = kx_standard(k0_valid);
    ky_valid = ky_standard(k0_valid);
    intensity_valid = intensity_standard(k0_valid);
    
    if isempty(kx_valid)
        kx_sim = []; ky_sim = []; I_sim = [];
        return;
    end
    
    % Create magnetic field function
    BFcn = create_magnetic_field_function(field_interpolants, parameters);
    
    % Vectorized trajectory calculation
    n_points = length(kx_valid);
    kx_final = zeros(n_points, 1);
    ky_final = zeros(n_points, 1);
    
    % Parallel processing if available
    if exist('parfor', 'builtin') && n_points > 100
        parfor i = 1:n_points
            try
                vx0 = CONSTANTS.hbar * kx_valid(i) * 1e10 / CONSTANTS.me;
                vy0 = CONSTANTS.hbar * ky_valid(i) * 1e10 / CONSTANTS.me;
                v_final = getFinalVelocity_optimized(BFcn, CONSTANTS.target_E * CONSTANTS.e, ...
                                                   vx0, vy0, CONSTANTS.z_target, CONSTANTS);
                kx_final(i) = v_final(1) * CONSTANTS.me / CONSTANTS.hbar * 1e-10;
                ky_final(i) = v_final(2) * CONSTANTS.me / CONSTANTS.hbar * 1e-10;
            catch
                kx_final(i) = NaN;
                ky_final(i) = NaN;
            end
        end
    else
        for i = 1:n_points
            try
                vx0 = CONSTANTS.hbar * kx_valid(i) * 1e10 / CONSTANTS.me;
                vy0 = CONSTANTS.hbar * ky_valid(i) * 1e10 / CONSTANTS.me;
                v_final = getFinalVelocity_optimized(BFcn, CONSTANTS.target_E * CONSTANTS.e, ...
                                                   vx0, vy0, CONSTANTS.z_target, CONSTANTS);
                kx_final(i) = v_final(1) * CONSTANTS.me / CONSTANTS.hbar * 1e-10;
                ky_final(i) = v_final(2) * CONSTANTS.me / CONSTANTS.hbar * 1e-10;
            catch
                kx_final(i) = NaN;
                ky_final(i) = NaN;
            end
        end
    end
    
    % Remove invalid results
    valid_results = ~isnan(kx_final) & ~isnan(ky_final);
    kx_sim = kx_final(valid_results);
    ky_sim = ky_final(valid_results);
    I_sim = intensity_valid(valid_results);
end

function BFcn = create_magnetic_field_function(field_interpolants, parameters)
    % Extract parameters
    transX = parameters(1); transY = parameters(2); transZ = parameters(3);
    thetaX = parameters(4); thetaY = parameters(5); thetaZ = parameters(6);
    current = parameters(7);
    
    % Rotation matrices
    cx = cos(thetaX); sx = sin(thetaX);
    cy = cos(thetaY); sy = sin(thetaY);
    cz = cos(thetaZ); sz = sin(thetaZ);
    Rx = [1 0 0; 0 cx -sx; 0 sx cx];
    Ry = [cy 0 sy; 0 1 0; -sy 0 cy];
    Rz = [cz -sz 0; sz cz 0; 0 0 1];
    Rot = Rz * Rx * Ry;
    
    field_ratio = current / 0.2;
    
    % Create optimized field function
    BFcn = @(x, y, z) compute_field_at_point(x, y, z, field_interpolants, ...
                                           Rot, transX, transY, transZ, field_ratio);
end

function B = compute_field_at_point(x, y, z, interpolants, Rot, transX, transY, transZ, field_ratio)
    % Transform coordinates
    pos_global = [x - transX; y - transY; z - transZ];
    pos_local = Rot * pos_global;
    
    % Interpolate field in local coordinates
    Bl = [interpolants.Bx(pos_local(1), pos_local(2), pos_local(3));
          interpolants.By(pos_local(1), pos_local(2), pos_local(3));
          interpolants.Bz(pos_local(1), pos_local(2), pos_local(3))];
    
    % Transform back to global coordinates
    B = Rot' * Bl * field_ratio;
end

function vxvy = getFinalVelocity_optimized(BFcn, Energy, vx0, vy0, z_target, CONSTANTS)
    % Adaptive time step Boris integration
    vz0 = sqrt(2 * Energy / CONSTANTS.me - vx0^2 - vy0^2);
    
    % Initial conditions
    x = [0; 0; 0];
    v_half = [vx0; vy0; vz0];
    
    % Adaptive time stepping
    B_initial = BFcn(0, 0, 0);
    omega_c = CONSTANTS.e * norm(B_initial) / CONSTANTS.me;
    dt = 0.02 * 2 * pi / omega_c;
    
    maxSteps = 500;
    for n = 1:maxSteps
        v_prev = v_half;
        % Boris rotation
        B = BFcn(x(1), x(2), x(3));
        tvec = (-CONSTANTS.e * dt / (2 * CONSTANTS.me)) * B;
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
                omega_c = CONSTANTS.e * B_norm / CONSTANTS.me;
                dt = 0.02 * 2 * pi / omega_c;
            end
        end
    end
    
    error('Target plane not reached within maxSteps');
end

function options = get_optimization_options(user_options, outfun)
    default_options = optimoptions('fmincon', ...
        'Display', 'iter-detailed', ...
        'OutputFcn', outfun, ...
        'MaxIterations', 500, ...
        'MaxFunctionEvaluations', 2000, ...
        'StepTolerance', 1e-8, ...
        'OptimalityTolerance', 1e-6, ...
        'UseParallel', true);
    
    if ~isempty(user_options)
        fn = fieldnames(user_options);
        for i = 1:length(fn)
            default_options.(fn{i}) = user_options.(fn{i});
        end
    end
    options = default_options;
end

function best_p = run_multi_start_optimization(objfun, initial_guess, lb, ub, options)
    % Single optimization run (can be extended to multi-start)
    try
        [best_p, fval, exitflag] = fmincon(objfun, initial_guess, [], [], [], [], lb, ub, [], options);
        fprintf('Dual optimization completed with exit flag: %d, final value: %.6f\n', exitflag, fval);
    catch ME
        warning('Optimization failed: %s', ME.message);
        best_p = initial_guess;
    end
end