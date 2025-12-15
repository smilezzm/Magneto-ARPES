function [Ftx, Fty] = calc_inverse_mapping_angle(standard_field, Energy, parameters, plot_grid)
    % The function returns a mapping function that find the (thetax,thetay) map
    % when the electrons just came out at the surface of the sample, from
    % the (thetax,thetay) map received by ARPES. 
    % The parameters is an object that contains .transX, .transY, .transZ,
    % .thetaX, .thetaY, .thetaZ, current. Degrees unit is rad!!
    % plot_grid = true/false to plot the grid of forward/inverse mapping
    % with provided parameters
    % Output:
    %   - Ftx: map (thetax,thetay) onto a thetax
    %   - Fty: map (thetax,thetay) onto a thetay

    %% Basic values
    hbar = 1.055e-34;
    me = 9.1093837015e-31; % Electron mass (kg)
    e = 1.602176634e-19;  % Elementary charge (C)
    % for fermi-surface, Energy = 110.56;  % in eV
    z_target = 0.02;
    gridNum = 200;
    
    %% Build angle-grid
    thetax_i = linspace(-19, 19, gridNum);   % in units deg
    thetay_i = linspace(-15,15,gridNum);   % Should be adjusted as needed
    [thetax_i, thetay_i] = ndgrid(thetax_i, thetay_i);
    thetax_i = thetax_i(:);
    thetay_i = thetay_i(:);

    %% Compute corresponding ks to simulate the trajectories
    kmag = sqrt(2 * Energy * e * me) / hbar / 1e10;
    kx_i = kmag * tand(thetay_i).*cosd(thetax_i)./sqrt(1+tand(thetay_i).^2.*cosd(thetax_i).^2);
    ky_i = kmag * tand(thetax_i).*cosd(thetay_i)./sqrt(1+tand(thetax_i).^2.*cosd(thetay_i).^2);

    %% Simulate the trajectories (in terms of kx,ky)
    kx_f = zeros(gridNum^2, 1);  % Column vector
    ky_f = zeros(gridNum^2, 1);  % Column vector
    
    field_interpolants = create_field_interpolants(standard_field);
    BFcn = create_magnetic_field_function(field_interpolants, parameters);
    
    if exist('parfor', 'builtin')
        parfor ii = 1:gridNum^2
            vx0 = hbar * kx_i(ii) * 1e10 / me;
            vy0 = hbar * ky_i(ii) * 1e10 / me;
            v_final = getFinalVelocity(BFcn, Energy * e, ...
                                               vx0, vy0, z_target);
            kx_f(ii) = v_final(1) * me / hbar * 1e-10;
            ky_f(ii) = v_final(2) * me / hbar * 1e-10;
        end
    else
        for ii = 1:gridNum^2
            vx0 = hbar * kx_i(ii) * 1e10 / me;
            vy0 = hbar * ky_i(ii) * 1e10 / me;
            v_final = getFinalVelocity(BFcn, Energy * e, ...
                                               vx0, vy0, z_target);
            kx_f(ii) = v_final(1) * me / hbar * 1e-10;
            ky_f(ii) = v_final(2) * me / hbar * 1e-10;
        end
    end

    %% Convert kx_f, ky_f into angles
    thetax_f = asind(ky_f ./ sqrt(kmag^2 - kx_f.^2));
    thetay_f = asind(kx_f ./ sqrt(kmag^2 - ky_f.^2));

    try
        Ftx = scatteredInterpolant(thetax_f, thetay_f, thetax_i, 'linear', 'none');
        Fty = scatteredInterpolant(thetax_f, thetay_f, thetay_i, 'linear', 'none');
    catch ME
        error('Failed to interpolate for this energy slice: %s', ME.message);
    end
    fprintf('reverse grid mapping has been established\n');

    if plot_grid
        % Plot: initial thetax-thetay to final thetax-thetay
        figure('Color', 'w');
        subplot(1,2,1);
        thetax_test = linspace(-18, 18, 37);
        thetay_test = linspace(-14, 14, 29);
        [thetax_test,thetay_test] = meshgrid(thetax_test, thetay_test);
        thetax_eval = Ftx(thetax_test, thetay_test);
        thetay_eval = Fty(thetax_test, thetay_test);
        h1 = scatter(thetax_eval(:), thetay_eval(:), 10, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
        hold on
        h2 = scatter(thetax_test(:), thetay_test(:), 10, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
        xlabel('thetax (deg)');
        ylabel('thetay (deg)');
        legend([h1, h2], {'Predicted intrinsic fermi surface', 'Measured fermi surface'}, 'Location', 'best');
        title('From measured with-field fermi surface to corrected fermi surface')
        grid on; axis equal;
        set(gcf, 'Color', 'w')  % gcf = get current figure
        
        % Plot: final kx-ky to initial kx-ky
        subplot(1,2,2);
        h1 = scatter(kx_i(:), ky_i(:),'b','filled','MarkerFaceAlpha',0.6);
        hold on
        h2 = scatter(kx_f(:), ky_f(:),'r','filled','MarkerFaceAlpha',0.6);
        xlabel('kx (10^{10} m^{-1})');
        ylabel('ky (10^{10} m^{-1})');
        legend([h1, h2], {'Intrinsic fermi surface', 'Predicted measured fermi surface'}, 'Location', 'best');
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
