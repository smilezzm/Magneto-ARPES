function [Fkx, Fky] = calc_inverse_mapping(standard_field, Energy, parameters, plot_grid)
    % The function returns a mapping function that find the (kx,ky) map
    % when the electrons just came out at the surface of the sample, from
    % the (kx,ky) map received by ARPES. 
    % The parameters is an object that contains .transX, .transY, .transZ,
    % .thetaX, .thetaY, .thetaZ, current. Degrees unit is rad!!
    % plot_grid = true/false to plot the grid of forward/inverse mapping
    % with provided parameters
    % Output:
    %   - Fkx: map (kx,ky) onto a kx
    %   - Fky: map (kx,ky) onto a ky
    %% Basic values
    hbar = 1.055e-34;
    me = 9.1093837015e-31; % Electron mass (kg)
    e = 1.602176634e-19;  % Elementary charge (C)
    % for fermi-surface, Energy = 110.56;  % in eV
    z_target = 0.02;
    gridNum = 100;
    
    kx_i = linspace(-1.72,1.72,gridNum);   % in units 1e10 m^(-1)
    ky_i = linspace(-1.72,1.72,gridNum);   % Should be adjusted as needed
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
    
    field_interpolants = create_field_interpolants(standard_field);
    BFcn = create_magnetic_field_function(field_interpolants, parameters);
    
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
        Fkx = scatteredInterpolant(kx_f, ky_f, kx_i_valid, 'linear', 'none');
        Fky = scatteredInterpolant(kx_f, ky_f, ky_i_valid, 'linear', 'none');
    catch ME
        error('Failed to interpolate for this energy slice: %s', ME.message);
    end
    fprintf('reverse grid mapping has been established\n');

    if plot_grid
        % Plot: initial kx-ky to final kx-ky
        figure('Color', 'w');
        subplot(1,2,1);
        h1 = scatter(kx_i_valid, ky_i_valid, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
        hold on
        h2 = scatter(kx_f, ky_f, 10, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
        xlabel('kx(10^{10} m^{-1})');
        ylabel('ky(10^{10} m^{-1})');
        legend([h1, h2], {'Intrinsic fermi surface', 'Predicted measured fermi surface'}, 'Location', 'best');
        title('Simulation of with-field fermi surface')
        grid on; axis equal;
        set(gcf, 'Color', 'w')  % gcf = get current figure
        
        % Plot: final kx-ky to initial kx-ky
        subplot(1,2,2);
        kx_test = linspace(-1, 1, 20);
        ky_test = kx_test;
        [kx_test,ky_test] = meshgrid(kx_test, ky_test);
        kx_eval = Fkx(kx_test, ky_test);
        ky_eval = Fky(kx_test, ky_test);
        valid_mask = ~isnan(kx_eval) & ~isnan(ky_eval);
        h1 = scatter(kx_test(:),ky_test(:),'r','filled','MarkerFaceAlpha',0.6);
        hold on
        h2 = scatter(kx_eval(valid_mask),ky_eval(valid_mask),'b','filled','MarkerFaceAlpha',0.6);
        xlabel('kx (10^{10} m^{-1})');
        ylabel('ky (10^{10} m^{-1})');
        legend([h1, h2], {'Measured fermi surface', 'Corrected fermi surface'}, 'Location', 'best');
        title('From measured with-field fermi surface to corrected fermi surface')
        set(gcf, 'Color', 'w')  % gcf = get current figure
        grid on; axis equal;
        
        % Add annotation to record the parameters
    annotation('textbox', [0.45 0.3 0.1 0.4], 'String', ...
        sprintf(['The parameters are:\n' ...
                 'transX = %.5f\ntransY = %.5f\ntransZ = %.5f\n' ...
                 'thetaX = %.4f\nthetaY = %.4f\nthetaZ = %.4f\n' ...
                 'current = %.3f'], ...
                 parameters.transX, parameters.transY, parameters.transZ, ...
                 parameters.thetaX, parameters.thetaY, parameters.thetaZ, ...
                 parameters.current), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize', 12);
    end
    
    % Store Fkx and Fky, which map (kx,ky) to kx0 and ky0 respectively
    % save('./matdata/inverse_mapping.mat','BFcn','Fkx','Fky','z_target');
end

function k0_valid = get_valid_mask(kx_standard, ky_standard, k_mag)
    k0_valid = (kx_standard.^2 + ky_standard.^2) <= (0.47 * k_mag)^2;
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
    transX = parameters.transX; transY = parameters.transY; transZ = parameters.transZ;
    thetaX = parameters.thetaX; thetaY = parameters.thetaY; thetaZ = parameters.thetaZ;
    current = parameters.current;
    
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
    pos_local = Rot * pos_global; % impose CW rotation thetaX, thetaY, thetaZ on the coil frame (during which the axes of thetaY, thetaZ are also rotating)
    % pos_global is the coordinates in lab frame, pos_local is the coordinates in coil frame
    
    % Interpolate field in local coordinates
    Bl = [interpolants.Bx(pos_local(1), pos_local(2), pos_local(3));
          interpolants.By(pos_local(1), pos_local(2), pos_local(3));
          interpolants.Bz(pos_local(1), pos_local(2), pos_local(3))];
    
    % Transform back to global coordinates
    B = Rot' * Bl * field_ratio; % Note here it's Rot not Rot' or Rz'Rx'Ry'
end

% Helper function to integrate the trajectory to get the final velocity
function vxvy = getFinalVelocity(BFcn, Energy, vx0, vy0, z_target)
    % Energy in (J), vx0 and vy0 in (m/s)
    e = 1.602176634e-19;  % Elementary charge (C)
    me = 9.1093837015e-31; % Electron mass (kg)
    B_mag_center = norm(BFcn(0,0,0));
    omega_c = e * B_mag_center / me;
    T_c = 2 * pi / omega_c;
    dt = T_c / 200;
    x = [0;0;0];
    q = -e;
    v_half = [vx0;vy0;sqrt(2*Energy/me-vx0^2-vy0^2)];
    maxSteps = 500;
    for n = 1:maxSteps
        v_prev = v_half;
        % Boris rotation (E=0)
        B = BFcn(x(1),x(2),x(3));              % [Bx; By; Bz] at position x
        tvec = (q*dt/(2*me)) * B;  % t vector
        t2 = dot(tvec, tvec);
        svec = 2*tvec/(1 + t2);
    
        v_minus = v_half;               % E-half kick skipped (E=0)
        v_prime = v_minus + cross(v_minus, tvec);
        v_plus  = v_minus + cross(v_prime, svec);
        v_half  = v_plus;               % velocity at (n+1/2)
    
        % Position update
        x_new = x + dt * v_half;
        % Check crossing of target plane
        if x_new(3) >= z_target
            % Linear interpolation for exact crossing point
            alpha = (z_target - x(3)) / (x_new(3) - x(3));
            % Interpolate velocity at exact crossing
            v_exact = v_prev + alpha * (v_half - v_prev);
            vxvy = v_exact(1:2);
            return
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
    error('Plane not reached within maxSteps.');
end
