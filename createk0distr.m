%% Helper function to build a series of initial v0
function v0_distr = createv0distr(deltaVx, deltaVy, vxElementNum, vyElementNum, v0Range, vzElementNum)
    % v0_distr is a 4-D matrix 
    % deltaVx = 1e4;
    % deltaVy = 1e4;
    % vxElementNum = 100;
    % vyElementNum = 100;
    % v0Range = [4e4, 3e5];
    % vzElementNum = 50; % also the number of energies for each (kx, ky)
    vx = linspace(-deltaVx/2, +deltaVx/2, vxElementNum);
    vy = linspace(-deltaVy/2, +deltaVy/2, vyElementNum);
    v0_distr = zeros(vxElementNum, vyElementNum, vzElementNum, 3);
    
    for ii = 1:vxElementNum
        for jj = 1:vyElementNum
            v0_distr(ii,jj,:,3) = createvzFromXY(vx(ii), vy(jj), v0Range, vzElementNum);
            v0_distr(ii,jj,:,1) = vx(ii);
            v0_distr(ii,jj,:,2) = vy(jj);
        end
    end
end


%% Helper function to create vz given a certain (vx, vy)
function vz_series = createvzFromXY(vx, vy, vrange, N)
% createvzFromXY  Generate N velocity vectors sharing px, py (vx, vy)
% The total speed squared is linearly sampled between v_min^2 and v_max^2.
%
% Inputs:
%   vx, vy   - Fixed x and y velocity components (scalars)
%   vrange   - 1x2 vector [v_min, v_max]; total speed limits (m/s)
%   N        - Number of samples (positive integer)
%
% Output:
%   vz_series - 1*N vector; [vz1 vz2 ... vzN]
%
% Example:
%   vz_series = createInitialVelocity2(1e5, 2e5, [3e5, 6e5], 5);

    % Basic input checks
    if numel(vrange) ~= 2
        error('vrange must be a 2-element vector [v_min, v_max].');
    end
    if ~isscalar(N) || N < 1 || N ~= floor(N)
        error('N must be a positive integer.');
    end

    v_min = vrange(1);
    v_max = vrange(2);
    if v_min > v_max
        error('v_min must be <= v_max.');
    end

    horiz2 = vx^2 + vy^2;
    if v_min^2 < horiz2
        error('v_min is too small: v_min^2 must be >= vx^2 + vy^2 to keep vz real.');
    end

    % Total speed squared samples
    v2_series = linspace(v_min^2, v_max^2, N);

    % Compute vz^2 ensuring numerical non-negativity
    vz2_series = v2_series - horiz2;
    if any(vz2_series < -1e-12)
        error('Computed vz^2 became negative (check inputs).');
    end
    vz2_series(vz2_series < 0) = 0; % clamp tiny negative due to rounding

    % Fill result
    vz_series = sqrt(vz2_series);
end