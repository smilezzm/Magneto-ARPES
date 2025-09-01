%% Electron Trajectory Calculation in Magnetic Field
% This code calculates the trajectory of electrons emitted from the coil surface

%% Constants
e = 1.602176634e-19;  % Elementary charge (C)
me = 9.1093837015e-31; % Electron mass (kg)
q_over_m = -e/me;     % Charge-to-mass ratio for electron (negative)

%% Initial conditions
% Starting position: central point of top surface of coil
r0 = [0; 0; airHeight/2 + coilHeight/2];  % Initial position [x; y; z]

% Initial velocity (example: electron emitted with some initial energy)
v0_magnitude = 1e5;  % Initial speed in m/s (adjust as needed)
% Direction: can be vertical, angled, or any direction you want
theta = 0.4;  % Angle from vertical (radians)
phi = 0;    % Azimuthal angle (radians)
v0 = v0_magnitude * [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];

%% Integration parameters (here use z instead of t)
dz = 2e-5;  % Time step (seconds) - very small for accurate integration
z_max = 8e-3; % Maximum simulation time (seconds)
z = 0:dz:z_max;
n_steps = length(z);

%% Initialize trajectory arrays
trajectory = zeros(3, n_steps);
velocity = zeros(3, n_steps);
trajectory(:,1) = r0;
velocity(:,1) = v0;

% use ODE solver method
trajectory_ode_solver(R, r0, v0, z_max, q_over_m, coilGm);

function trajectory_ode_solver(R, r0, v0, z_max, q_over_m, coilGm)
    electronMotion = @(z, y) localelectronMotion(z, y, q_over_m, R);
    % Define ODE function
    function dydz = localelectronMotion(z, y, q_over_m, R)
        % y = [x, y, z, vx, vy, vz]
        pos = y(1:3);
        vel = y(4:6);
        intrpB = interpolateMagneticFlux(R, pos(1), pos(2), pos(3));
        B = [intrpB.Bx; intrpB.By; intrpB.Bz];
        v_cross_B = cross(vel, B) / y(6);
        accel = q_over_m * v_cross_B;
        dydz = [y(4:6)/y(6); accel];
    end
    
    % Initial conditions for ODE solver
    y0 = [r0; v0];
    
    % Solve ODE
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
    [t_ode, y_ode] = ode45(electronMotion, [0, z_max], y0, options);
    
    % Extract trajectory
    trajectory_ode = y_ode(:, 1:3)';
    
    % Plot ODE solution
    figure;
    plot3(trajectory_ode(1,:), trajectory_ode(2,:), trajectory_ode(3,:), 'r-', 'LineWidth', 2);
    hold on;
    plot3(r0(1), r0(2), r0(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    hold on;
    pdegplot(coilGm, FaceAlpha=0.2);
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('Electron Trajectory (ODE Solver)');
    grid on;
    axis equal;
end