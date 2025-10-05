function transformed_points = rigid_transform(original_points, tx, ty, theta)
    % theta is the CCW angle in degree 
    kx = original_points.kx;
    ky = original_points.ky;
    I_f = original_points.I_thetax_thetay;
    transformed_points = struct();

    kx_f = kx + tx;
    ky_f = ky + ty;
    Rot = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
    k_stacked = Rot * [kx_f(:)'; ky_f(:)'];
    kx_f = reshape(k_stacked(1,:), size(kx));
    ky_f = reshape(k_stacked(2,:), size(ky));
    transformed_points.kx = kx_f;
    transformed_points.ky = ky_f;
    transformed_points.I = I_f;
end