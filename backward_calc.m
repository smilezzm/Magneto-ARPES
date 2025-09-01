function [vx0, vy0, vz0] = backward_calc(DA_pos, v_final, R, z0, options)
    % Honestly, this function is useless right now, due to difficulties
    % when finding the outcoming position of electrons.
    %
    % backward_calc finds the initial velocity of an electron comming out at (r,
    % theta, beta) with velocity (vx; vy; vz).
    %   (r, theta, beta) is the current position of electrons when entering
    %   electronic lens. However, we are not sure about the exact
    %   position of the electron coming in, because the entrance is a
    %   relatively large area. 
    % Inputs:
    %   DA_pos - position vector of somewhere at the DA entrance in arpes 
    %            coordinate (r, theta, beta)
    %   v_final - velocity vector of the electrons entering the lens (vx,
    %             vy, vz). And the +e 'electron' will be virtually emitted
    %             with -v_final
    %   R - an object containing information about the field
    %   z0 - z value of the sample, also the height of emitting point. This
    %        may be airHeight/2 + coilHeight/2
    %   options - sturct with additional parameters
    %             .BFcn - the pre-sample function that build B on a grid to
    %                     enhance speed
    %             .RelTol - relative tolerance of the ODE solver
    %             .AbsTol - absolute tolerance of the ODE solver
    
    if ~isfield(options, 'RelTol'), options.RelTol = 1e-9; end
    if ~isfield(options, 'AbsTol'), options.AbsTol = 1e-12; end
    
    %% Physical constants
    e = 1.602176634e-19;  % Elementary charge (C)
    me = 9.1093837015e-31; % Electron mass (kg)
    q_over_m = e/me;     % Charge-to-mass ratio for electron (positive)

    %% Define the ODE function
    function dydz = electronMotion(z, y)
        % y = [x, y, z, vx, vy, vz]
        pos = y(1:3);
        vel = y(4:6);
        % --- Use fast interpolant if provided ---
        if isfield(options,'BFcn')
            B = options.BFcn(pos(1), pos(2), pos(3));
        else
            intrpB = interpolateMagneticFlux(R, pos(1), pos(2), pos(3));
            B = [intrpB.Bx; intrpB.By; intrpB.Bz];
        end
        v_cross_B = cross(vel, B);
        accel = q_over_m * v_cross_B;
        dydz = [y(4:6) / y(6); accel / y(6)];
    end
    
    %% Define y0 (6,1), containing initial condition of ODE
    y0 = zeros(6,1);
    y0(1:3) = arpes2cart(DA_pos(1), DA_pos(2), DA_pos(3));
    y0(4:6) = - v_final;  
    z_max = z0 - y0(3); % negative
    
    %% Solve the differential equation
    ode_options = odeset('RelTol', options.RelTol, 'AbsTol', options.AbsTol);
    [~, y_solution] = ode45(@electronMotion, [0, z_max], y0, ode_options);

    vx0 = - y_solution(end, 4);
    vy0 = - y_solution(end, 5);
    vz0 = - y_solution(end, 6);

end

%% Helper function to convert spherical coordinates in ARPES to cartesian
function position = arpes2cart(r, theta, beta)
    x = r * sin(beta);
    y = r * cos(beta) * sin(theta);
    z = r * cos(beta) * cos(theta);
    position = [x;y;z];
end