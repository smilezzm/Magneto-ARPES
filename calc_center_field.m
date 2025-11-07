function field = calc_center_field(standard_field, parameters, position)
    field_interpolants = create_field_interpolants(standard_field);
    BFcn = create_magnetic_field_function(field_interpolants, parameters);
    field = BFcn(position(1), position(2), position(3));
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