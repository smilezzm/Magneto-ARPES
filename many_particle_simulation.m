%% Set the velocity distribution
deltaVx = 2e6;
deltaVy = 2e6;
vxElementNum = 30;
vyElementNum = 30;
v0Range = [1e6, 6e6];
vzElementNum = 8; % also the number of energies for each (kx, ky)
v0_distr = createv0distr(deltaVx, deltaVy, vxElementNum, vyElementNum, v0Range, vzElementNum);
v0_valid = true(vxElementNum, vyElementNum, vzElementNum); % v0 is valid if theta is not too big
for ii = 1:vxElementNum
    for jj = 1:vyElementNum
        for kk = 1:vzElementNum
            if sqrt(v0_distr(ii, jj, kk, 1)^2+v0_distr(ii, jj, kk, 2)^2) > v0_distr(ii, jj, kk, 3)*0.8
                v0_valid(ii, jj, kk) = false;
            end
        end
    end
end
final_velocity = zeros(vxElementNum, vyElementNum, vzElementNum, 3); % the same size as v0_distr
z_max = 1e-2; % upper limit during trajectory integration 

%% Set the displacement of the coil (translation and rotation)
coilRin = 0.004;
coilRout = 0.007;
coilHeight = 0.004;
airR = 0.05;
airHeight = 0.05;
currentDensity = 2E5;
geometry.coilRin = coilRin;
geometry.coilRout = coilRout;
geometry.coilHeight = coilHeight;
geometry.airR = airR;
geometry.airHeight = airHeight;
fieldOptions.plot = true;
fieldOptions.rotation.theta = 0 * 180 / pi ;
fieldOptions.rotation.beta = 0 * 180 / pi;
fieldOptions.rotation.phi = 0 * 180 / pi;
[R, coilGm, airGm] = calc_field(currentDensity, geometry, fieldOptions);

%% --- Build magnetic field interpolants (choose region of interest) ---
% When rerunning this program, there is no need to build interpolants
% When running for the first time, please uncomment this part
xg = linspace(-3e-3, 3e-3, 20);
yg = linspace(-3e-3, 3e-3, 20);
zg = linspace( (airHeight/2+coilHeight/2), (airHeight/2+coilHeight/2)+1.1*z_max, 20);
[Xg,Yg,Zg] = ndgrid(xg,yg,zg);
Bsample = R.interpolateMagneticFlux(Xg(:), Yg(:), Zg(:));
BxI = griddedInterpolant(Xg, Yg, Zg, reshape(Bsample.Bx, size(Xg)), 'linear', 'nearest');
ByI = griddedInterpolant(Xg, Yg, Zg, reshape(Bsample.By, size(Xg)), 'linear', 'nearest');
BzI = griddedInterpolant(Xg, Yg, Zg, reshape(Bsample.Bz, size(Xg)), 'linear', 'nearest');
BFcn = @(x,y,z) [BxI(x,y,z); ByI(x,y,z); BzI(x,y,z)];
fprintf('BFcn constructed successfully!\n');

%% prepare for trajectory integration
options = struct();
options.RelTol = 5e-7;
options.AbsTol = 5e-10; % slightly relaxed for speed

%% --- Parallel loop over all (ii,jj,kk) ---
% strangely, there are numerical errors when using parallel loops ???
% totalN = vxElementNum * vyElementNum * vzElementNum;
% entry_initial_lin = false(totalN,1);   % linear workspace variable
% parfor linearIdx = 1:totalN
%     [ii,jj,kk] = ind2sub([vxElementNum, vyElementNum, vzElementNum], linearIdx);
%     y0 = [0; 0; airHeight/2+coilHeight/2; ...
%           v0_distr(ii,jj,kk,1); v0_distr(ii,jj,kk,2); v0_distr(ii,jj,kk,3)];
%     [~, ~, ~, slit_entry] = calc_faked_initial(R, y0, z_max, options);
%     entry_initial_lin(linearIdx) = slit_entry.reached;
%     fprintf('vx0=%.3e, vy0=%.3e, vz0=%.3e', y0(4), y0(5), y0(6));
% end
% entry_initial = reshape(entry_initial_lin, vxElementNum, vyElementNum, vzElementNum);


%% Integrate the path of each electron
for ii = 1:vxElementNum
    for jj = 1:vyElementNum
        for kk = 1:vzElementNum
            if ~v0_valid(ii, jj, kk)
                continue
            end
            y0 = [0;0;airHeight/2+coilHeight/2;v0_distr(ii,jj,kk,1);v0_distr(ii,jj,kk,2);v0_distr(ii,jj,kk,3)];
            final_velocity(ii,jj,kk,:) = calc_velocity(BFcn, y0, z_max, options);
        end
    end
end

%% draw the comparison between initial (vx,vy) and final (vx, vy)
% for z_num = 1:vzElementNum
%     initial_vx = v0_distr(:,:,z_num,1);
%     initial_vy = v0_distr(:,:,z_num,2);
%     initial_vx_flat = initial_vx(:);
%     initial_vy_flat = initial_vy(:);
%     final_vx = final_velocity(:,:,z_num,1);
%     final_vy = final_velocity(:,:,z_num,2);
%     final_vx_flat = final_vx(:);
%     final_vy_flat = final_vy(:);
%     figure
%     scatter(initial_vx_flat, initial_vy_flat, 10, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
%     hold on
%     scatter(final_vx_flat, final_vy_flat, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
%     xlabel('vx (m/s)');
%     ylabel('vy (m/s)');
%     title(sprintf('Velocity distribution (vx, vy) for slice z_n = %d', z_num));
%     grid on;
% end


%% Build mapping with interpolation (First fix the energy, then interpolate in 2D)
% If directly interpolating in 3D, the relation may not be one on one
% Create separate interpolants for each vz slice
invMapCell = cell(vzElementNum, 1);

for z_slice = 1:vzElementNum
    % Extract data for this energy slice
    mask_2d = v0_valid(:, :, z_slice);
    
    if sum(mask_2d(:)) < 4  % Need at least 4 points for interpolation
        continue;
    end
    
    initial_vx_2d = v0_distr(:, :, z_slice, 1);
    initial_vy_2d = v0_distr(:, :, z_slice, 2);
    final_vx_2d = final_velocity(:, :, z_slice, 1);
    final_vy_2d = final_velocity(:, :, z_slice, 2);
    
    % Flatten and filter valid points
    initial_vx_flat = initial_vx_2d(mask_2d);
    initial_vy_flat = initial_vy_2d(mask_2d);
    final_vx_flat = final_vx_2d(mask_2d);
    final_vy_flat = final_vy_2d(mask_2d);
    
    % Create 2D interpolants (final -> initial)
    try
        Gx_2d = scatteredInterpolant(final_vx_flat, final_vy_flat, initial_vx_flat, 'linear', 'none');
        Gy_2d = scatteredInterpolant(final_vx_flat, final_vy_flat, initial_vy_flat, 'linear', 'none');
        
        invMapCell{z_slice} = @(vx_query, vy_query) [Gx_2d(vx_query, vy_query), Gy_2d(vx_query, vy_query)];
    catch ME
        warning('Failed to create interpolant for z_slice %d: %s', z_slice, ME.message);
        invMapCell{z_slice} = [];
    end
end










%% Function that integrate the trajectories and find the final velocity
function final_velocity = calc_velocity(BFcn, y0, z_max, options)
    e = 1.602176634e-19;  % Elementary charge (C)
    me = 9.1093837015e-31; % Electron mass (kg)
    q_over_m = -e/me;     % Charge-to-mass ratio for electron (negative)
    function dydz = electronMotion(z, y)
        % y = [x, y, z, vx, vy, vz]
        pos = y(1:3);
        vel = y(4:6);
        B = BFcn(pos(1), pos(2), pos(3));
        v_cross_B = cross(vel, B);
        accel = q_over_m * v_cross_B;
        dydz = [y(4:6) / y(6); accel / y(6)];
    end
    ode_options = odeset('RelTol', options.RelTol, 'AbsTol', options.AbsTol);
    [~, y_solution] = ode45(@electronMotion, [0, z_max], y0, ode_options);
    final_velocity = y_solution(end, 4:6)';
end

%% Helper function used when plotting the slit
function position = arpes2cart(r, theta, beta)
    x = r * sin(beta);
    y = r * cos(beta) * sin(theta);
    z = r * cos(beta) * cos(theta);
    position = [x,y,z];
end