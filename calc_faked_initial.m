function [faked_initial, trajectory, velocity, slit_entry] = calc_faked_initial(R, y0, z_max, options)
    % At first I thought it was important to consider whether electrons
    % will enter slit, but there is an electronic lens before the slit
    % which receive electrons with θ in [-15, 15].
    % So we can neglect anything about slit. Just set options.slit.enabled=false. 
    % Currently I only care about the starting and ending velocity, so this function is useless.
    %
    % Inputs:
    %   R - FEM result object containing magnetic field data
    %   y0 - Initial conditions [x0; y0; z0; vx0; vy0; vz0] (6x1 vector)
    %        Position in meters, velocity in m/s
    %   z_max - Maximum simulation height (from 0 to ..., not from the actual starting z value!)
    %   options - (optional) Structure with additional parameters:
    %             .plot - true/false to show trajectory plot (default: false)
    %             .reverseExtend - true/false to show the line that
    %                              extended from the final velocity reversely
    %             .coilGm - its type should be the discrete pde geometry.
    %                       If don't want it, then don't include it into
    %                       options
    %             .BFcn - the pre-sample function that build B on a grid to
    %                     enhance speed
    %             .slitGm - true/false to show slit geometry for visualization
    %             .RelTol - Relative tolerance for ODE solver (default: 1e-9)
    %             .AbsTol - Absolute tolerance for ODE solver (default: 1e-12)
    %             .verbose - true/false for console output (default: false)
    %             .slit - a struct containing the information of the position of the slit of ARPES
    %                     .enabled - true/false
    %                     .theta0 - the central position of the slit
    %                     .beta0 - the central position of the slit
    %                     .radius - the central position of the slit 
    %                     All of which are estimated in spherical coordinates for ARPES (slightly different from the normal spherical coordinates)
    %                     .deltaTheta - corresponding the 'width' of the slit
    %                     .deltaBeta - corresponding the 'length' of the slit
    %                     .tolerance - the range of radius so that particles entering a volume are all considered getting into the slit 
    % Outputs:
    %   faked_initial - Initial position and direction backwards asymptote,
    %                   if the particle doesn't enter the slit, then it's calculated from
    %                   the final position (at z_max)
    %                   if the particle enter the slit, then it's calculated from the
    %                   position&velocity when entering slit
    %   trajectory - Trajectory matrix [x; y; z] (3xN)
    %   velocity - Velocity matrix [vx; vy; vz] (3xN)
    %   slit_entry - a struct containing 
    %                .reached - ture or false
    %                .position - position when entering the slit, or []
    %                .velocity - velocity when entering the slit, or []
    %
    % Example:
    %   r0 = [0; 0; 0.01];  % Initial position (m)
    %   v0 = [1e5; 0; 0];   % Initial velocity (m/s)
    %   y0 = [r0; v0];
    %   [init, traj, vel, entry] = calc_faked_initial(R, y0, 2e-3);
    
    %% Input validation and default options
    if nargin < 4
        options = struct();
    end
    
    % Set default options
    if ~isfield(options, 'plot'), options.plot = false; end
    if ~isfield(options, 'reverseExtend'), options.reverseExtend = false; end
    if ~isfield(options, 'RelTol'), options.RelTol = 1e-9; end
    if ~isfield(options, 'AbsTol'), options.AbsTol = 1e-12; end
    if ~isfield(options, 'verbose'), options.verbose = false; end
    if ~isfield(options, 'slitGm'), options.slitGm = false; end
    % Add slit parameters to options if not present
    if ~isfield(options, 'slit'), options.slit = struct(); end
    if ~isfield(options.slit, 'enabled'), options.slit.enabled = false; end
    if ~isfield(options.slit, 'theta0'), options.slit.theta0 = 0.2; end
    if ~isfield(options.slit, 'beta0'), options.slit.beta0 = 0.2; end
    if ~isfield(options.slit, 'radius'), options.slit.radius = 0.007; end % Distance from origin
    if ~isfield(options.slit, 'deltaTheta'), options.slit.deltaTheta = 0.04; end
    if ~isfield(options.slit, 'deltaBeta'), options.slit.deltaBeta = 0.2; end
    if ~isfield(options.slit, 'tolerance'), options.slit.tolerance = 6e-4; end 
    % the last one is "thickness" in radius, which means electrons entering a
    % volume of space are all considered entering into slit
    
    % Validate input dimensions
    if length(y0) ~= 6
        error('Initial conditions y0 must be a 6-element vector [x0; y0; z0; vx0; vy0; vz0]');
    end
    
    %% Physical constants
    e = 1.602176634e-19;  % Elementary charge (C)
    me = 9.1093837015e-31; % Electron mass (kg)
    q_over_m = -e/me;     % Charge-to-mass ratio for electron (negative)
    
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

    %% Determine whether the electron has entered the slit.
    function [value, isterminal, direction] = slitEvent(z, y)
        if ~options.slit.enabled
            value = 1;        % No event
            isterminal = 0;   % Don't stop
            direction = 0;    % Any direction
            return;
        end
        
        % Convert to spherical coordinates for ARPES (different from the normal spherical coordinates)
        [r, theta, beta] = cart2arpes(y(1), y(2), z); 
        % do not use y(3) because y(3) stands for actual z component, here
        % z is calculated from the starting plane
        
        % Check if electron is within slit angular range
        v_theta = abs(wrapToPi(theta - options.slit.theta0)) - options.slit.deltaTheta/2;
        v_beta = abs(wrapToPi(beta - options.slit.beta0)) - options.slit.deltaBeta/2;
        v_radius = abs(r - options.slit.radius) - options.slit.tolerance/2;
        % value crosses 0 when electron enters slit
        value = max([v_theta, v_beta, v_radius]); % so that value is continuous, preventing potential error when solving ODE
        isterminal = 1;   % stop integrating when entering slit
        direction = -1;   % only entering
    end
    
    %% Solve the ODE
    if options.verbose
        fprintf('Starting electron trajectory calculation...\n');
        fprintf('Initial position: [%.3e, %.3e, %.3e] m\n', y0(1), y0(2), y0(3));
        fprintf('Initial velocity: [%.3e, %.3e, %.3e] m/s\n', y0(4), y0(5), y0(6));
        fprintf('Simulation z range: 0 - %.3e m\n', z_max);
        if options.slit.enabled
            fprintf('Slit detection enabled: θ = %.3f, Δθ = %.3f, β = %.3f, Δβ = %.3f, r = %.3e, tolerance = %.3e\n', ...
                options.slit.theta0, options.slit.deltaTheta, ...
                options.slit.beta0, options.slit.deltaBeta, options.slit.radius, options.slit.tolerance)
        end
    end
    
    % Set ODE solver options and solve ODEs
    if options.slit.enabled
        ode_options = odeset('RelTol', options.RelTol, 'AbsTol', options.AbsTol, ...
                        'Events', @slitEvent);
        [z_solution, y_solution, z_event, y_event, ~] = ode45(@electronMotion, [0, z_max], y0, ode_options);
    else
        ode_options = odeset('RelTol', options.RelTol, 'AbsTol', options.AbsTol);
        [z_solution, y_solution] = ode45(@electronMotion, [0, z_max], y0, ode_options);
    end
    

    % Check if slit was reached
    slit_reached = ~isempty(z_event) && options.slit.enabled;
    if slit_reached && options.verbose
        fprintf('Electron reached slit at z = %.6e m\n', z_event(end));
    end
    
    % Record slit entry information
    if slit_reached
        slit_entry.reached = true;
        % I hope there is at most one element in y_event and z_event
        % check it
        if length(z_event)>1
            error("length of z_event larger than 1")
        end
        % check over
            
        slit_entry.position = y_event(1, 1:3).';
        slit_entry.velocity = y_event(1, 4:6).';
        
        if options.verbose
            fprintf('Slit entry - Absolute Position: [%.3e, %.3e, %.3e] m\n', slit_entry.position);
            fprintf('Slit entry - Velocity: [%.3e, %.3e, %.3e] m/s\n', slit_entry.velocity);
        end
    else
        slit_entry.reached = false;
        slit_entry.position = [];
        slit_entry.velocity = [];
    end

    % Extract trajectory and velocity
    trajectory = y_solution(:, 1:3)';  % [x; y; z] vs z
    velocity = y_solution(:, 4:6)';    % [vx; vy; vz] vs z
    faked_initial = calculate_backwards_asymptote(y_solution(end, 1:3), y_solution(end, 4:6), y0(3));


    if options.verbose
        fprintf('Trajectory calculation completed.\n');
        fprintf('Total z steps: %d\n', length(z_solution));
        fprintf('Final position: [%.3e, %.3e, %.3e] m\n', trajectory(1,end), trajectory(2,end), trajectory(3,end));
        fprintf('Final speed: %.3e m/s\n', norm(velocity(:,end)));
    end
    
    %% Optional plotting
    if options.plot
        figure;
        plot3(trajectory(1,:), trajectory(2,:), trajectory(3,:), 'r-', 'LineWidth', 2);
        hold on;
        
        % Mark start and end points
        plot3(trajectory(1,1), trajectory(2,1), trajectory(3,1), 'go', ...
            'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
        hold on;
        plot3(trajectory(1,end), trajectory(2,end), trajectory(3,end), 'ro', ...
            'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'End');
        hold on;
        % Plot coil geometry if provided
        if isfield(options, 'coilGm')
            pdegplot(options.coilGm, 'FaceAlpha', 0.2);
        end
        if options.slitGm
            P1 = arpes2cart(options.slit.radius, options.slit.theta0-options.slit.deltaTheta/2, options.slit.beta0-options.slit.deltaBeta/2);
            P2 = arpes2cart(options.slit.radius, options.slit.theta0+options.slit.deltaTheta/2, options.slit.beta0-options.slit.deltaBeta/2);
            P3 = arpes2cart(options.slit.radius, options.slit.theta0+options.slit.deltaTheta/2, options.slit.beta0+options.slit.deltaBeta/2);
            P4 = arpes2cart(options.slit.radius, options.slit.theta0-options.slit.deltaTheta/2, options.slit.beta0+options.slit.deltaBeta/2);
            X = [P1(1), P2(1), P3(1), P4(1)];
            Y = [P1(2), P2(2), P3(2), P4(2)];
            Z = [P1(3), P2(3), P3(3), P4(3)] + y0(3);
            hold on
            fill3(X, Y, Z, 'c');
            grid on;
        end
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        title('Electron Trajectory in Magnetic Field');
        legend('Trajectory', 'Start', 'End');
        grid on;
        axis equal;
    end
    if options.reverseExtend
        figure;
        plot3([faked_initial(1), trajectory(1,end)], [faked_initial(2), trajectory(2, end)], [faked_initial(3), trajectory(3,end)], '--', 'Color', 'r', 'Linewidth', 2);
        if isfield(options, 'coilGm')
            hold on
            pdegplot(options.coilGm, 'FaceAlpha', 0.2);
        end
        if options.slitGm
            P1 = arpes2cart(options.slit.radius, options.slit.theta0-options.slit.deltaTheta/2, options.slit.beta0-options.slit.deltaBeta/2);
            P2 = arpes2cart(options.slit.radius, options.slit.theta0+options.slit.deltaTheta/2, options.slit.beta0-options.slit.deltaBeta/2);
            P3 = arpes2cart(options.slit.radius, options.slit.theta0+options.slit.deltaTheta/2, options.slit.beta0+options.slit.deltaBeta/2);
            P4 = arpes2cart(options.slit.radius, options.slit.theta0-options.slit.deltaTheta/2, options.slit.beta0+options.slit.deltaBeta/2);
            X = [P1(1), P2(1), P3(1), P4(1)];
            Y = [P1(2), P2(2), P3(2), P4(2)];
            Z = [P1(3), P2(3), P3(3), P4(3)] + y0(3); % Notice to translate the coordinate relative to the starting point to absolute coordinate
            hold on
            fill3(X, Y, Z, 'c');
            grid on;
        end
        grid on;
        title('Extended reversely from the final velocity (faked initial position)')
        xlabel('X'); ylabel('Y'); zlabel('Z');
    end
end


%% Helper function to find the faked initial position and velocity, 
% obtained from extending the final velocity vector in the reversed direction to intersect with the xy plane at the initial position.
function faked_initial = calculate_backwards_asymptote(final_pos, final_vel, z_initial)
    % CALCULATE_BACKWARDS_ASYMPTOTE Find intersection of backwards velocity vector with xy plane
    %
    % Inputs:
    %   final_pos - Final position [x_f; y_f; z_f] (3x1 vector)
    %   final_vel - Final velocity [vx_f; vy_f; vz_f] (3x1 vector)
    %   z_initial - Z coordinate of the xy plane to intersect with (scalar)
    %
    % Output:
    %   faked_initial - [x_init; y_init; z_init; vx_f; vy_f; vz_f] (6x1 vector)
    %                   Position where backwards line intersects xy plane + final velocity
    
    % Extract components
    x_f = final_pos(1);
    y_f = final_pos(2);
    z_f = final_pos(3);
    
    vx_f = final_vel(1);
    vy_f = final_vel(2);
    vz_f = final_vel(3);
    
    % Check if velocity has z-component (avoid division by zero)
    if abs(vz_f) < 1e-15
        warning('Final velocity has no z-component, cannot find intersection with xy plane');
        faked_initial = [final_pos; final_vel];
        return;
    end
    
    % Parametric line equation: r(t) = r_final - t * v_final
    % We want to find t such that z(t) = z_initial
    % z(t) = z_f - t * vz_f = z_initial
    % Solving for t: t = (z_f - z_initial) / vz_f
    
    t_intersect = (z_f - z_initial) / vz_f;
    
    % Calculate intersection point
    x_init = x_f - t_intersect * vx_f;
    y_init = y_f - t_intersect * vy_f;
    z_init = z_initial;  % By definition
    
    % Return faked initial conditions
    faked_initial = [x_init; y_init; z_init; vx_f; vy_f; vz_f];
end

%% Helper function to convert cartesian to spherical coordinates 
% (notice the spherical coordinates in ARPES is different from the normal one)
function [r, theta, beta] = cart2arpes(x,y,z)
    r = sqrt(x^2 + y^2 + z^2);
    theta = atan(y/z);
    beta = asin(x/r);
end

%% Helper function to convert spherical coordinates in ARPES to cartesian
function position = arpes2cart(r, theta, beta)
    x = r * sin(beta);
    y = r * cos(beta) * sin(theta);
    z = r * cos(beta) * cos(theta);
    position = [x,y,z];
end