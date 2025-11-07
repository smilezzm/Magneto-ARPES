function parameters = optimized_parameters(standard_measured, final_measured,...
    k_r, standard_field, initial_guess, similarity_function, options)
    % This function returns the parameters that make the forward simulated
    % fermi surface figure match the best with measured figure, for the
    % "calibrating" sample or "standard" sample. 
    %
    % Inputs:
    %   - standard_measured: the measured fermi surface of the standard sample with field off
    %       .kx, .ky, .I_thetax_thetay
    %   - final_measured: the measured fermi surface of the standard sample with field on
    %       .kx, .ky, .I_thetax_thetay
    %   - k_r: the radius of k (in 1e10 m^(-1)) that we are interested at the center 
    %   - standard_field: an object containing the standard field values on grids
    %       .Bx, .By, .Bz, .X, .Y, .Z
    %   - initial_guess: a vector [transX; transY; transZ; thetaX; thetaY; thetaZ; current]
    %   - options: (optional) optimization options structure
    %   - similarity_function: 'innerProduct' or 'ssim'， which specify the
    %       function that quantitize similarity
    % Outputs:
    %   - parameters: optimized parameters structure
    
    if nargin < 7
        options = struct();
    end
    if ~strcmp(similarity_function,'innerProduct') && ~strcmp(similarity_function,'ssim')
        error('Please assign a similarity_function in ''innerProduct'' and ''ssim'' ');
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
    [kx_final, ky_final, intensity_final] = extract_data(final_measured);
    validate_inputs(kx_standard, ky_standard, intensity_standard, ...
                   kx_final, ky_final, intensity_final);
    
    % Pre-compute k magnitude and valid mask
    k_mag = sqrt(2 * CONSTANTS.target_E * CONSTANTS.e * CONSTANTS.me) / CONSTANTS.hbar / 1e10;
    k0_valid = get_valid_mask(kx_standard, ky_standard, k_mag);
    if strcmp(similarity_function, 'innerProduct')
        % Pre-compute interpolation grid for measured data
        center_mask = (kx_final.^2 + ky_final.^2 < k_r^2);
        kx_final_ROI = kx_final(center_mask);
        ky_final_ROI = ky_final(center_mask);
        intensity_final_ROI = intensity_final(center_mask); % A vector
    else
        kx_ROI = linspace(-k_r,k_r,70);
        ky_ROI = linspace(-k_r,k_r,70);
        [kx_ROI,ky_ROI]=ndgrid(kx_ROI,ky_ROI);
        F_final = scatteredInterpolant(kx_final(:), ky_final(:), intensity_final(:), 'linear', 'none');
        intensity_final_ROI = F_final(kx_ROI, ky_ROI); % A 2D array
    end
    
    % Pre-compute field interpolants
    field_interpolants = create_field_interpolants(standard_field);
    
    % Setup optimization
    history = [];
    if strcmp(similarity_function, 'innerProduct')
        objfun = @(p) -calc_similarity_innerProduct(p, kx_standard, ky_standard, intensity_standard,...
            kx_final_ROI, ky_final_ROI, intensity_final_ROI, ...
            field_interpolants, CONSTANTS, k0_valid);
    else
        objfun = @(p) -calc_similarity_ssim(p, kx_standard, ky_standard, intensity_standard,...
            kx_ROI, ky_ROI, intensity_final_ROI, ...
            field_interpolants, CONSTANTS, k0_valid);
    end
    
    % Set optimization options
    opt_options = get_optimization_options(options, @outfun);
    lb = [-1e-3, -1e-3, -1e-3, -0.1*pi, -0.1*pi, -pi, 0.068];
    ub = [1e-3, 1e-3, 1e-3, 0.1*pi, 0.1*pi, pi, 0.08];
    
    % Run optimization with multiple starting points for robustness
    fprintf('Starting optimization...\n');
    best_p = run_multi_start_optimization(objfun, initial_guess, lb, ub, opt_options);
    
    % Visualize results
    visualize_results(kx_standard, ky_standard, intensity_standard, ...
                     kx_final, ky_final, intensity_final, ...
                     field_interpolants, CONSTANTS, k0_valid, best_p, history, similarity_function);
    
    % Return structured parameters
    parameters = struct(...
        'transX', best_p(1), 'transY', best_p(2), 'transZ', best_p(3), ...
        'thetaX', best_p(4), 'thetaY', best_p(5), 'thetaZ', best_p(6), ...
        'current', best_p(7));
    save('./matdata/optimized_p_plus_innerProduct.mat','parameters');
    
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

function [kx, ky, intensity] = extract_data(data_struct)
    kx = data_struct.kx;
    ky = data_struct.ky;
    intensity = double(data_struct.I_thetax_thetay);
end

function validate_inputs(kx_standard, ky_standard, intensity_standard, ...
                        kx_final, ky_final, intensity_final)
    if ~(isequal(size(kx_standard), size(ky_standard), size(intensity_standard)) && ...
         isequal(size(kx_final), size(ky_final), size(intensity_final)))
        error('Input arrays must have matching dimensions');
    end
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

function score = calc_similarity_innerProduct(parameters, kx_standard, ky_standard, intensity_standard,...
    kx_final_ROI, ky_final_ROI, intensity_final_ROI, ...
    field_interpolants, CONSTANTS, k0_valid)
    % Notice here kx_final_ROI, ky_final_ROI, intensity_final_ROI are
    % vectors. kx_sim, ky_sim, I_sim below are also vectors
    
    % Forward simulation
    [kx_sim, ky_sim, I_sim] = forward_sim_vectorized(kx_standard, ky_standard, intensity_standard,...
        field_interpolants, CONSTANTS, k0_valid, parameters);
    
    % Interpolate simulation results to measured grid
    if length(kx_sim) < 3
        score = -Inf; % Penalty for insufficient data
        return;
    end
    
    try
        F_ROI = scatteredInterpolant(kx_sim, ky_sim, I_sim, 'linear', 'none');
        intensity_sim_ROI = F_ROI(kx_final_ROI, ky_final_ROI);
        
        % Remove NaN values
        valid_idx = ~isnan(intensity_sim_ROI);
        if sum(valid_idx) < 0.8 * length(intensity_sim_ROI)
            score = -Inf;
            return;
        end
        
        intensity_final_clean = intensity_final_ROI(valid_idx);
        intensity_sim_clean = intensity_sim_ROI(valid_idx);
        
        % Calculate normalized correlation
        if std(intensity_final_clean) == 0 || std(intensity_sim_clean) == 0
            score = -Inf;
            return;
        end
        
        norm_final = (intensity_final_clean - mean(intensity_final_clean)) / std(intensity_final_clean);
        norm_sim = (intensity_sim_clean - mean(intensity_sim_clean)) / std(intensity_sim_clean);
        
        score = sum(norm_final .* norm_sim) / (length(norm_final) - 1);
        
    catch
        score = -Inf;
    end
end

function score = calc_similarity_ssim(parameters, kx_standard, ky_standard, intensity_standard,...
    kx_ROI, ky_ROI, intensity_final_ROI, ...
    field_interpolants, CONSTANTS, k0_valid)
    % Notice here kx_final_ROI, ky_final_ROI, intensity_final_ROI are
    % vectors. kx_sim, ky_sim, I_sim below are also vectors
    
    % Forward simulation
    [kx_sim, ky_sim, I_sim] = forward_sim_vectorized(kx_standard, ky_standard, intensity_standard,...
        field_interpolants, CONSTANTS, k0_valid, parameters);

    % Interpolate simulation results to measured grid
    if length(kx_sim) < 3
        score = -Inf; % Penalty for insufficient data
        return;
    end

    try
        F_ROI = scatteredInterpolant(kx_sim, ky_sim, I_sim, 'linear', 'none');
        intensity_sim_ROI = F_ROI(kx_ROI, ky_ROI);
        
        % Remove NaN values
        valid_idx = ~isnan(intensity_sim_ROI) & ~isnan(intensity_final_ROI);
        if sum(valid_idx(:)) < 0.7 * numel(intensity_sim_ROI)
            score = -Inf;
            return;
        end

        intensity_final_ROI(~valid_idx)=0;
        intensity_sim_ROI(~valid_idx)=0;

        % Calculate normalized correlation
        if std(intensity_final_ROI(:)) == 0 || std(intensity_sim_ROI(:)) == 0
            score = -Inf;
            return;
        end

        [~, ssim_map] = ssim(intensity_final_ROI, intensity_sim_ROI);
        score = mean(ssim_map(valid_idx));
        
    catch
        score = -Inf;
    end
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
        'UseParallel', false);
    
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
        fprintf('Optimization completed with exit flag: %d, final value: %.6f\n', exitflag, fval);
    catch ME
        warning('Optimization failed: %s', ME.message);
        best_p = initial_guess;
    end
end

function visualize_results(kx_standard, ky_standard, intensity_standard, ...
                          kx_final, ky_final, intensity_final, ...
                          field_interpolants, CONSTANTS, k0_valid, best_p, history, similarity_function)
    % Plot optimization history
    if ~isempty(history)
        figure('Name', 'Optimization Progress');
        plot(-history(1, :), 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Similarity Score');
        title('Optimization Progress');
        grid on;
    end
    
    % Generate final simulation for visualization
    [kx_sim, ky_sim, I_sim] = forward_sim_vectorized(kx_standard, ky_standard, intensity_standard,...
        field_interpolants, CONSTANTS, k0_valid, best_p);
    
    % Create comparison plot
    figure('Color', 'w');
    k_r=1.4;
    center_mask = (kx_final.^2 + ky_final.^2 < k_r^2);
    kx_final_ROI = kx_final(center_mask);
    ky_final_ROI = ky_final(center_mask);
    I_final_ROI = intensity_final(center_mask);

    F_I = scatteredInterpolant(kx_sim,ky_sim,I_sim,'linear','none');
    I_sim_intp = F_I(kx_final_ROI, ky_final_ROI);
    clean_mask = ~isnan(I_sim_intp) & ~isnan(I_final_ROI);
    subplot(1,2,1)
    scatter(kx_final_ROI(clean_mask), ky_final_ROI(clean_mask), 10, I_final_ROI(clean_mask)-I_sim_intp(clean_mask),'filled');
    axis equal tight;
    xlabel('k_x [10^{10} m^{-1}]');
    ylabel('k_y [10^{10} m^{-1}]');
    title(['Difference between measured and simulated +field fermi surface (measured-simulated)', ...
       newline, ...
       'Using similarity function: ', char(similarity_function)]);
    annotation('textbox', [0.6 0.3 0.1 0.4], 'String', ...
    sprintf(['The parameters are:\n' ...
             'transX = %.5f\ntransY = %.5f\ntransZ = %.5f\n' ...
             'thetaX = %.4f\nthetaY = %.4f\nthetaZ = %.4f\n' ...
             'current = %.3f'], ...
             best_p(1), best_p(2), best_p(3), ...
             best_p(4), best_p(5), best_p(6), ...
             best_p(7)), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 12);
    colormap(turbo); colorbar;
end
