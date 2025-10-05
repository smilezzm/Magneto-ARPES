function [kx_f,ky_f,I_f] = forward_mapping(standard_measured, inverse_mapping)
% - standard_measured: the fermi surface from ARPES data 0mA(without field)
%   .kx, .ky, .I_thetax_thetay, .Ef
%   kx,ky are (N,M) obtained from ARPES data
%   I_thetax_thetay is the intensity of the same size
% - inverse_mapping: 
%   Fkx, Fky : inverse mapping
%   BFcn : interpolant of the field for the chosen parameters
%   target_z: the upper limit of z. It is directly assigned in
%   calc_inverse_mapping.m
    
    hbar = 1.055e-34;
    me = 9.1093837015e-31; % Electron mass (kg)
    e = 1.602176634e-19;  % Elementary charge (C)
    BFcn = inverse_mapping.BFcn;
    z_target = inverse_mapping.z_target;
    kx = standard_measured.kx;
    ky = standard_measured.ky;
    I_f = standard_measured.I_thetax_thetay;

    k_mag = sqrt(2 * Ef * e * me) / hbar / 1e10;
    k0_valid = true(size(kx));
    [N1,N2] = size(kx);
    for ii = 1:N1
        for jj = 1:N2
            if sqrt(kx(ii,jj)^2+ky(ii,jj)^2) > 0.5 * k_mag
                k0_valid(ii,jj) = false;
            end
        end
    end
    kx_f = zeros(size(kx));
    ky_f = zeros(size(kx));
    
    for ii = 1:N1
        for jj = 1:N2
            if ~k0_valid(ii,jj)
                kx_f(ii,jj) = NaN;
                ky_f(ii,jj) = NaN;
            else
                vx0 = hbar * kx(ii,jj) * 1e10 / me;
                vy0 = hbar * ky(ii,jj) * 1e10 / me;
                v_final = getFinalVelocity(BFcn, Ef * e, vx0, vy0, z_target);
                kx_f(ii,jj) = v_final(1) * me / hbar * 1e-10;
                ky_f(ii,jj) = v_final(2) * me / hbar * 1e-10;
            end
        end
    end
    
    figure('Color','w');
    scatter(kx_f(:), ky_f(:), 10, I_f(:), 'filled');
    axis equal tight;
    xlabel('k_x [10^{10} m^{-1}]');
    ylabel('k_y [10^{10} m^{-1}]');
    title(sprintf('forward mapping @ E = %.3f eV', Ef));
    colormap turbo; colorbar;
end

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
