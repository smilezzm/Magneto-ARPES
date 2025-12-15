function parameters = optimized_parameters_dual_correction(final_measured_plus, final_measured_minus,...
    k_r, standard_field, initial_guess, lb, ub, ax1, ax2, Ef, deltaE, options)
    % This function returns the parameters that make the corrected
    % fermi surface figures match the best from both measured figures, for the
    % plus field and minus field
    % Sometimes, we can fix the thetaX, thetaY, thetaZ and only tune the
    % four parameters left.
    % The final_measured_plus and final_measured_minus are .mat files that
    % contain 3D band information
    % 
    % Inputs:
    %   - final_measured_plus: the measured fermi surface with positive field
    %       .kx, .ky, .I_thetax_thetay, .E 
    %       [E,thetax,thetay]            [1, E]
    %   - final_measured_minus: the measured fermi surface with negative field
    %       .kx, .ky, .I_thetax_thetay, .E
    %   - k_r: the radius of k (in 1e10 m^(-1)) that we are interested at the center 
    %       In this case, from experience, k_r could be 1.1
    %   - standard_field: an object containing the standard field values on grids
    %       .Bx, .By, .Bz, .X, .Y, .Z
    %   - initial_guess: a vector [transX; transY; transZ; thetaX; thetaY; thetaZ; current]
    %   - lb, ub: lower bound and upper bound for the parameters
    %
    %   ANGLE ARE IN UNITS OF rad, LENGTH ARE IN UNITS OF mm, CURRENT IS IN
    %   UNITS OF A.     So that the parameters are of similar scales
    %
    %   - ax: the axes to which the figures are plotted
    %   - Ef: energy of the slice which you'd like to correct
    %   - deltaE: energy range around Ef (to take average over a few
    %   slices), in eV
    %   - options: Provide
    %       options.free_indices to choose which of the 7 parameters are
    %       optimized (indices follow the order of the initial_guess
    %       vector); the remaining entries stay fixed at their initial
    %       values.
    %       options.similarity_method: 'delta' or 'ssim'
    %           'delta' is recommended to finetune the already good
    %           parameters
    % Outputs:
    %   - parameters: optimized parameters structure

    if nargin < 10
        options = struct();
        options.similarity_method = 'delta';
    end

    free_idx = 1:7;
    if isstruct(options) && isfield(options, 'free_indices')
        free_idx = options.free_indices;
    end
    free_idx = free_idx(:)';
    if isempty(free_idx)
        error('At least one parameter must be selected for optimization.');
    end
    if any(free_idx < 1 | free_idx > 7)
        error('free_indices must reference positions within the 7-parameter vector.');
    end
    free_idx = unique(free_idx, 'stable');
    initial_guess = initial_guess(:);
    lb = lb(:);
    ub = ub(:);
    if numel(initial_guess) ~= 7 || numel(lb) ~= 7 || numel(ub) ~= 7
        error('initial_guess, lb, and ub must each have 7 elements.');
    end
    param_template = initial_guess;
    
    % Extract and validate inputs
    Ef_idx_plus = find(abs(final_measured_plus.E - Ef) < deltaE/2);
    Ef_idx_minus = find(abs(final_measured_minus.E - Ef) < deltaE/2);
    if isempty(Ef_idx_minus) | isempty(Ef_idx_plus)
        error('please make sure there are slices in range [Ef-deltaE/2, Ef+deltaE/2]');
    end

    [kx_plus, ky_plus, intensity_plus] = extract_data(final_measured_plus, Ef_idx_plus);
    [kx_minus, ky_minus, intensity_minus] = extract_data(final_measured_minus, Ef_idx_minus);
    validate_inputs_dual(kx_plus, ky_plus, intensity_plus, ...
                        kx_minus, ky_minus, intensity_minus);
    kx_ROI = linspace(-k_r,k_r,100);
    ky_ROI = linspace(-k_r,k_r,100);
    [kx_ROI,ky_ROI]=ndgrid(kx_ROI,ky_ROI);
    
    
    % Pre-compute field interpolants
    field_interpolants = create_field_interpolants(standard_field);
    
    % Setup optimization
    history = [];
    objfun_full = @(p) -calc_dual_similarity_optimized(p, ...
        kx_ROI, ky_ROI, kx_plus, ky_plus,  intensity_plus, ...
        kx_minus, ky_minus, intensity_minus, ...
        field_interpolants, k_r, options.similarity_method);
    objfun = @(p_free) objfun_full(expand_params(p_free));
    initial_guess_free = initial_guess(free_idx);
    lb_free = lb(free_idx);
    ub_free = ub(free_idx);

    % Set optimization options
    opt_options = get_optimization_options(@outfun);

    % Run optimization with multiple starting points for robustness
    fprintf('Starting dual-field optimization...\n');
    best_p_free = run_multi_start_optimization(objfun, initial_guess_free, lb_free, ub_free, opt_options);
    best_p = expand_params(best_p_free);

    % Return structured parameters
    parameters = struct(...
        'transX', best_p(1), 'transY', best_p(2), 'transZ', best_p(3), ...
        'thetaX', best_p(4), 'thetaY', best_p(5), 'thetaZ', best_p(6), ...
        'current', best_p(7));
    % save('./matdata/optimized_p_dual_correction.mat','parameters');

    % Visualize results
    visualize_dual_results(kx_plus, ky_plus, intensity_plus, ...
                          kx_minus, ky_minus, intensity_minus, ...
                          field_interpolants, kx_ROI, ky_ROI, ...
                          best_p, history, ax1, ax2);
    
    function p_full = expand_params(p_free)
        p_full = param_template;
        p_full(free_idx) = p_free;
        p_full = p_full(:);
    end

    function stop = outfun(p_free, optimValues, state)
        stop = false;
        if isequal(state, 'iter')
            p_full = expand_params(p_free);
            history = [history, [optimValues.fval; p_full(:)]];
            if mod(optimValues.iteration, 10) == 0
                fprintf('Iteration %d: f = %.6f\n', optimValues.iteration, optimValues.fval);
            end
        end
    end
end

function validate_inputs_dual(kx_plus, ky_plus, intensity_plus, ...
                             kx_minus, ky_minus, intensity_minus)
    if ~(isequal(size(kx_plus), size(ky_plus), size(intensity_plus)) && ...
         isequal(size(kx_minus), size(ky_minus), size(intensity_minus)))
        error('Input arrays must have matching dimensions');
    end
end

function score = calc_dual_similarity_optimized(parameters, ...
    kx_ROI, ky_ROI, kx_plus, ky_plus,  intensity_plus, ...
        kx_minus, ky_minus, intensity_minus, ...
    field_interpolants, k_r, method)
    
    parameters(1:3) = parameters(1:3)/1000; % convert to m
    BFcn = create_magnetic_field_function(field_interpolants, parameters);
    [Fkx_plus, Fky_plus] = inverse_mapping(BFcn);
    parameters(7) = - parameters(7);
    BFcn = create_magnetic_field_function(field_interpolants, parameters);
    [Fkx_minus, Fky_minus] = inverse_mapping(BFcn);
    
    kx_corrected_plus = Fkx_plus(kx_plus, ky_plus);
    ky_corrected_plus = Fky_plus(kx_plus, ky_plus);
    kx_corrected_minus = Fkx_minus(kx_minus, ky_minus);
    ky_corrected_minus = Fky_minus(kx_minus, ky_minus);   % (n,m) with nan
    
    % Calculate similarity for positive field
    score = calc_similarity(kx_corrected_plus, ky_corrected_plus, intensity_plus, ...
                               kx_corrected_minus, ky_corrected_minus, intensity_minus, ...
                               kx_ROI, ky_ROI, k_r, method);
    
    % Debug output every 50 evaluations
    persistent eval_count;
    if isempty(eval_count)
        eval_count = 0;
    end
    eval_count = eval_count + 1;
    if mod(eval_count, 50) == 0
        fprintf('  Score breakdown: %.4f\n', score);
    end
end

function score = calc_similarity(kx_corrected_plus, ky_corrected_plus, intensity_plus, ...
                                       kx_corrected_minus, ky_corrected_minus, intensity_minus, ...
                                       kx_ROI, ky_ROI, k_r, method)
    if strcmp(method, 'delta')
        ROI_mask = (kx_ROI.^2 + ky_ROI.^2 <= k_r^2);
        kx_ROI_vector = kx_ROI(ROI_mask);
        ky_ROI_vector = ky_ROI(ROI_mask);
        valid_idx = ~isnan(kx_corrected_plus);
        F_plus_ROI = scatteredInterpolant(kx_corrected_plus(valid_idx), ...
            ky_corrected_plus(valid_idx), ...
            intensity_plus(valid_idx), 'linear', 'none');
        valid_idx = ~isnan(kx_corrected_minus);
        F_minus_ROI = scatteredInterpolant(kx_corrected_minus(valid_idx), ...
            ky_corrected_minus(valid_idx), ...
            intensity_minus(valid_idx), 'linear', 'none');
        I_plus_ROI_vector = F_plus_ROI(kx_ROI_vector, ky_ROI_vector);
        I_minus_ROI_vector = F_minus_ROI(kx_ROI_vector, ky_ROI_vector);
        valid_idx = ~isnan(I_plus_ROI_vector) & ~isnan(I_minus_ROI_vector);
        if sum(valid_idx(:)) < 0.7 * numel(I_minus_ROI_vector)
            score = -Inf;
            return;
        end
        I_plus_clean = I_plus_ROI_vector(valid_idx);
        I_minus_clean = I_minus_ROI_vector(valid_idx);
        % norm_I_plus = (I_plus_clean - mean(I_plus_clean)) / std(I_plus_clean);
        % norm_I_minus = (I_minus_clean - mean(I_minus_clean)) / std(I_minus_clean);
        % score = - mean(abs(norm_I_plus-norm_I_minus));
        score = - mean(abs(I_plus_clean - I_minus_clean));

    elseif strcmp(method, 'ssim')
        try
            % Interpolate backmapping results to grid
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
    
            % Remove NaN values
            valid_idx = ~isnan(intensity_plus_ROI) & ~isnan(intensity_minus_ROI);
            if sum(valid_idx(:)) < 0.7 * numel(intensity_minus_ROI)
                score = -Inf;
                return;
            end
            intensity_plus_ROI(~valid_idx)=0;
            intensity_minus_ROI(~valid_idx)=0;
            
            % Calculate normalized correlation
            if std(intensity_plus_ROI(:)) == 0 || std(intensity_minus_ROI(:)) == 0
                score = -Inf;
                return;
            end
    
            [~, ssim_map] = ssim(intensity_plus_ROI, intensity_minus_ROI);
            score = mean(ssim_map(valid_idx));
        catch
            score = -Inf;
        end
    else
        error('input options.similarity_method must be either ''ssim'' or ''delta'' ');
    end
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

    best_p(1:3) = best_p(1:3)/1000;
    
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

    vals = intensity_plus_ROI(valid_idx) - intensity_minus_ROI(valid_idx);
    scatter(ax2, kx_ROI(valid_idx), ky_ROI(valid_idx), 10, vals, 'filled');
    axis(ax2, 'equal');
    axis(ax2, 'tight');
    xlabel(ax2, 'k_x [10^{10} m^{-1}]');
    ylabel(ax2, 'k_y [10^{10} m^{-1}]');
    title(ax2, 'Difference between plus/minus-field fermi surface (plus-minus)');
    % Build a smooth colormap
    n = 256; % number of steps in the gradient
    blue_to_white = [linspace(0,1,n/2)' linspace(0,1,n/2)' ones(n/2,1)];
    white_to_red   = [ones(n/2,1) linspace(1,0,n/2)' linspace(1,0,n/2)'];
    cmap = [blue_to_white; white_to_red];
    colormap(ax2, cmap);
    clim(ax2, [-max(abs(vals)), max(abs(vals))]);
    colorbar(ax2);

    % Print final parameters
    fprintf('\nOptimized Parameters:\n');
    fprintf('Translation: X=%.6f, Y=%.6f, Z=%.6f [mm]\n', best_p(1), best_p(2), best_p(3));
    fprintf('Rotation: X=%.4f, Y=%.4f, Z=%.4f [rad]\n', best_p(4), best_p(5), best_p(6));
    fprintf('Current: %.6f [A]\n', -best_p(7));
end

% Include all the helper functions from the original script
function [kx, ky, intensity] = extract_data(data_struct, Ef_idx)
    intensity = double(squeeze(mean(data_struct.I_thetax_thetay(Ef_idx,:,:),1)));
    Ef_idx = round(median(Ef_idx));
    kx = squeeze(data_struct.kx(Ef_idx,:,:));
    ky = squeeze(data_struct.ky(Ef_idx,:,:));
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

function [Fkx, Fky] = inverse_mapping(BFcn)
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
    Energy = 110.56;  % in eV
    z_target = 0.02;
    gridNum = 100;
    
    kx_i = linspace(-1.7,1.7,gridNum);   % in units 1e10 m^(-1)
    ky_i = linspace(-1.7,1.7,gridNum);   % Should be adjusted as needed
    [kx_i, ky_i] = ndgrid(kx_i, ky_i);
    k_mag = sqrt(2 * Energy * e * me) / hbar / 1e10;
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
            v_final = getFinalVelocity(BFcn, Energy * e, ...
                                               vx0, vy0, z_target);
            kx_f(ii) = v_final(1) * me / hbar * 1e-10;
            ky_f(ii) = v_final(2) * me / hbar * 1e-10;
        end
    else
        for ii = 1:n_points
            vx0 = hbar * kx_i_valid(ii) * 1e10 / me;
            vy0 = hbar * ky_i_valid(ii) * 1e10 / me;
            v_final = getFinalVelocity(BFcn, Energy * e, ...
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


function BFcn = create_magnetic_field_function(field_interpolants, parameters)
    % units of parameters: m and rad and A
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
    Rot = Rz * Rx * Ry;  % maybe Rz * Rx * Ry
    
    field_ratio = current / 0.2;
    
    % Create optimized field function
    BFcn = @(x, y, z) compute_field_at_point(x, y, z, field_interpolants, ...
                                           Rot, transX, transY, transZ, field_ratio);
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

function options = get_optimization_options(outfun)
    default_options = optimoptions('fmincon', ...
        'Display','iter-detailed', ...
        'OutputFcn',outfun, ...
        'MaxIterations',500, ...          % allow more iterations
        'MaxFunctionEvaluations',3000, ...  % allow more evaluations
        'StepTolerance',1e-3, ...         % smaller step stop criterion
        'FunctionTolerance',1e-4, ...     % keep iterating until objective change is tiny
        'OptimalityTolerance',1e-4, ...    % tighten first-order optimality tolerance
        'UseParallel',true);
    
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