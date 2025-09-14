if exist('coilDout','var')&&exist('coilDin','var')&&...
        exist('coilHeight','var')&&exist('current','var')&&...
        exist('X','var')&&exist('Y','var')&&exist('Z','var')&&...
        exist('Bx','var')&&exist('By','var')&&exist('Bz','var')
else
    try
        load('standard_field.mat');
    catch
        error('You should first calculate a "standard" field with calc_standard_field.m');
    end
end
%% Basic values
hbar = 1.055e-34;
me = 9.1093837015e-31; % Electron mass (kg)
e = 1.602176634e-19;  % Elementary charge (C)
photon_energy = 21.2;
work_function = 4.29;
binding_energy = 0.0;
Energy = photon_energy - work_function - binding_energy;  % in eV
z_target = 0.01;
real_current = 0.025;
field_ratio = real_current / current;  % The basic field is calculated with no translation, no rotation and current of 0.2A.
gridNum = 100;
thetax = 0/180*pi;  % radiant
thetay = 5/180*pi;
thetaz = 0;
translationX = 0;
translationY = 0;
translationZ = 0;

%% Construct interpolant from the grid value of standard field
Bxl = griddedInterpolant(X,Y,Z,Bx,'linear','nearest');
Byl = griddedInterpolant(X,Y,Z,By,'linear','nearest');
Bzl = griddedInterpolant(X,Y,Z,Bz,'linear','nearest');

%% Build the magnetic field with rotated&translated coil from the field without rotation&translation (stored in R)
% Also, a field_ratio term is multiplied for the current applied here. 
% The standard field is calculated based on I=0.2A
if evalin('base', 'exist(''BFcn'', ''var'')')
    fprintf('interpolant of real field has been established\n');
else
    cx = cos(thetax); sx = sin(thetax);
    cy = cos(thetay); sy = sin(thetay);
    cz = cos(thetaz); sz = sin(thetaz);
    Rx = [1 0 0; 0 cx -sx; 0 sx cx];
    Ry = [cy 0 sy; 0 1 0; -sy 0 cy];
    Rz = [cz -sz 0; sz cz 0; 0 0 1];
    Rot = Rz*Ry*Rx;

    xg = linspace(-coilDout, coilDout, 100);      % global(lab) coordinate range
    yg = xg;
    zg = linspace(-coilHeight, z_target, 100);    % mainly care the part that electrons may pass
    [Xg,Yg,Zg] = ndgrid(xg,yg,zg);  % lab coordinates
    
    % Stack into a 3 x (N*M*L) matrix of vectors
    posStack = [(Xg(:)-translationX)'; (Yg(:)-translationY)'; (Zg(:)-translationZ)'];   % size 3 x (N*M*L)
    % Apply the 3x3 rotation/transform matrix R
    positionStackl = Rot' * posStack;                % still 3 x (N*M*L)
    % Reshape back to N x M x L
    Xl = reshape(positionStackl(1,:), size(Xg));  % local coordinates (N,M,L)
    Yl = reshape(positionStackl(2,:), size(Xg));
    Zl = reshape(positionStackl(3,:), size(Xg));
    Blx_sample = Bxl(Xl(:),Yl(:),Zl(:));  % Blx_sample size: (N*M*L,1)
    Bly_sample = Byl(Xl(:),Yl(:),Zl(:));
    Blz_sample = Bzl(Xl(:),Yl(:),Zl(:));
    BlStack = [Blx_sample'; Bly_sample'; Blz_sample'];   
    BgStack = Rot * BlStack * field_ratio;
    Bxg = griddedInterpolant(Xg, Yg, Zg, reshape(BgStack(1,:), size(Xg)), 'linear', 'nearest');
    Byg = griddedInterpolant(Xg, Yg, Zg, reshape(BgStack(2,:), size(Xg)), 'linear', 'nearest');
    Bzg = griddedInterpolant(Xg, Yg, Zg, reshape(BgStack(3,:), size(Xg)), 'linear', 'nearest');
    BFcn = @(x,y,z) [Bxg(x,y,z); Byg(x,y,z); Bzg(x,y,z)];
    fprintf('interpolant of real field has been established\n')
end
   
%% build the grid mapping
fprintf('start building the reverse grid mapping\n')
if evalin('base', 'exist(''Fkx'', ''var'') && exist(''Fky'', ''var'')')
    fprintf('reverse grid mapping has been established\n');
else
    kx_i = linspace(-0.9,0.9,gridNum);   % in units 1e10 m^(-1)
    ky_i = linspace(-0.9,0.9,gridNum);
    k_mag = sqrt(2*me*Energy*e)/hbar/1e10;
    [kx_i, ky_i] = ndgrid(kx_i, ky_i);
    k0_valid = true(size(kx_i));
    for ii = 1:gridNum
        for jj = 1:gridNum
            if sqrt(kx_i(ii,jj)^2+ky_i(ii,jj)^2) > 0.4 * k_mag
                k0_valid(ii,jj) = false;
            end
        end
    end
    kx_f = zeros(size(kx_i));
    ky_f = zeros(size(kx_i));
    
    for ii = 1:gridNum
        for jj = 1:gridNum
            if ~k0_valid(ii,jj)
                kx_f(ii,jj) = NaN;
                ky_f(ii,jj) = NaN;
            else
                vx0 = hbar * kx_i(ii,jj) * 1e10 / me;
                vy0 = hbar * ky_i(ii,jj) * 1e10 / me;
                v_final = getFinalVelocity(BFcn, Energy * e, vx0, vy0, z_target);
                kx_f(ii,jj) = v_final(1) * me / hbar * 1e-10;
                ky_f(ii,jj) = v_final(2) * me / hbar * 1e-10;
            end
        end
    end
    try
        Fkx = scatteredInterpolant(kx_f(k0_valid), ky_f(k0_valid), kx_i(k0_valid), 'linear', 'nearest');
        Fky = scatteredInterpolant(kx_f(k0_valid), ky_f(k0_valid), ky_i(k0_valid), 'linear', 'nearest');
    catch ME
        error('Failed to interpolate for this energy slice: %s', ME.message);
    end
    fprintf('reverse grid mapping has been established\n');
end

%% Plot: initial kx-ky to final kx-ky
figure
h1 = scatter(kx_i(k0_valid), ky_i(k0_valid), 10, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
hold on
h2 = scatter(kx_f(k0_valid), ky_f(k0_valid), 10, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('kx(10^{10} m^{-1})');
ylabel('ky(10^{10} m^{-1})');
legend([h1, h2], {'Initial kx-ky', 'Predicted final kx-ky'}, 'Location', 'best');
title('from initial kx-ky to final kx-ky')
grid on
set(gcf, 'Color', 'w')  % gcf = get current figure

%% Plot: final kx-ky to initial kx-ky
figure
kx_test = linspace(-0.6, 0.6, 15);
ky_test = kx_test;
[kx_test,ky_test] = meshgrid(kx_test, ky_test);
kx_eval = Fkx(kx_test, ky_test);
ky_eval = Fky(kx_test, ky_test);
h1 = scatter(kx_test(:),ky_test(:),'r','filled','MarkerFaceAlpha',0.6);
hold on
h2 = scatter(kx_eval(:),ky_eval(:),'b','filled','MarkerFaceAlpha',0.6);
xlabel('kx (10^{10} m^{-1})');
ylabel('ky (10^{10} m^{-1})');
legend([h1, h2], {'Final kx-ky', 'Predicted initial kx-ky'}, 'Location', 'best');
title('from final kx-ky to initial kx-ky')
set(gcf, 'Color', 'w')  % gcf = get current figure
grid on

%% Store Fkx and Fky, which map (kx,ky) to kx0 and ky0 respectively
save('inverse_mapping.mat','BFcn','Fkx','Fky','z_target','real_current');

%% Helper function to integrate the trajectory to get the final velocity
function vxvy = getFinalVelocity(BFcn, Energy, vx0, vy0, z_target)
    % Energy in (J), vx0 and vy0 in (m/s)
    e = 1.602176634e-19;  % Elementary charge (C)
    me = 9.1093837015e-31; % Electron mass (kg)
    B_mag_center = norm(BFcn(0,0,0));
    omega_c = e * B_mag_center / me;
    T_c = 2 * pi / omega_c;
    dt = T_c / 100;
    x = [0;0;0];
    q = -e;
    v_half = [vx0;vy0;sqrt(2*Energy/me-vx0^2-vy0^2)];
    t = 0;
    maxSteps = 5000;
    for n = 1:maxSteps
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
        x = x + dt * v_half;
        t = t + dt;
        % Check crossing of target plane
        if x(3) >= z_target
            vxvy = v_half(1:2);
            return
        end
    end
    error('Plane not reached within maxSteps.');
end
