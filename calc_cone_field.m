function field = calc_cone_field(z_target, doPlot)
% Roughly calculate the field for cone-like permalloy geometry.
% There is no current here, instead, the Ni cylinder is magnetized at . 
% Notice z_target should be smaller than airHeight/2

    %% Geometry parameters
    airD = 50e-3;
    airHeight = 50e-3;
    coneHeight = 6e-3;
    coneDtop = 1e-3;
    coneDbottom = 6e-3;
    coreHeight = 6e-3;
    permeability = 83000;  % of permalloy

    %% Create cone geometry using STL
    % Generate truncated cone as STL file
    stl_path = './matdata/cone_core.stl';
    create_cone_stl(stl_path, coneDtop, coneDbottom, coneHeight);
    
    % Import the STL geometry
    conegm = importGeometry(stl_path);

    %% Build air geometry and core geometry
    airgm = multicylinder(airD/2, airHeight);
    airgm = translate(airgm, [0, 0, -airHeight/2]);

    coregm = multicylinder(coneDbottom/2, coreHeight);
    coregm = translate(coregm, [0, 0, - coneHeight - coreHeight - 1e-4]);

    % Combine geometries
    gm = addCell(airgm, conegm);
    gm = addCell(gm, coregm);
    gm = fegeometry(gm);

    %% Set the model
    model = femodel(AnalysisType='magnetostatic', Geometry=gm);
    model.VacuumPermeability = 1.2566370614E-6;
    model.MaterialProperties = materialProperties(RelativePermeability=1);

    % Find permalloy material cells (cone region)
    cone_ids = findCell(gm, [0,0,-coneHeight/2]);
    model.MaterialProperties(cone_ids) = materialProperties(RelativePermeability=permeability);

    % Assign the Ni with magnetization of 0.25T/mu_0
    magnetization = 0.2 / model.VacuumPermeability;
    core_ids = findCell(gm, [0,0,-coneHeight-coreHeight/2]);
    model.CellLoad(core_ids) = cellLoad("Magnetization",[0, 0, magnetization]);

    % boundary condition
    model.FaceBC(1:3) = faceBC(MagneticPotential=[0;0;0]);  % 1:3 are airgm face indices

    % Generate mesh with reasonable parameters
    model = generateMesh(model, 'Hmax', 4e-3, 'Hmin', 0.1e-3);

    if doPlot
        figure;
        pdegplot(model, 'CellLabels', 'on', 'FaceAlpha', 0.5);
        title('Geometry: Cone + Core + Coil');
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        axis equal; view(3);

        figure;
        pdemesh(model, FaceAlpha=0.3);  % it's hard to show the inner mesh, so you'll see only the mesh at airgm surface
        title('Finite Element Mesh');   
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        axis equal; view(3);
    end

    %% Solve
    result = solve(model);

    %% Extract magnetic field at target plane
    x = linspace(-1.5 * coneDbottom, 1.5 * coneDbottom, 120);
    y = x;
    z = linspace(- 1.5 * coreHeight, z_target, 140);
    [X, Y, Z] = ndgrid(x, y, z);
    Bdata = interpolateMagneticFlux(result, X, Y, Z); % struct with fields Bx, By, Bz

    % Shape to grid and package
    Bx = reshape(Bdata.Bx, size(X));
    By = reshape(Bdata.By, size(Y));
    Bz = reshape(Bdata.Bz, size(Z));

    field = struct('X', X, 'Y', Y, 'Z', Z, 'Bx', Bx, 'By', By, 'Bz', Bz);
    save('./matdata/cone_field.mat', 'X','Y','Z','Bx','By','Bz','z_target');

    if doPlot
        %% Plot field magnitude at z=1e-4 and 2e-3
        % z=1e-4
        figure;
        subplot(1,2,1);
        [X, Y] = ndgrid(x,y);
        Z_low = 1e-4 * ones(size(X));
        B_data = interpolateMagneticFlux(result, X, Y, Z_low);
        Br = sqrt(B_data.Bx.^2+B_data.By.^2);
        Bz_max = max(B_data.Bz);
        Br_max = max(Br);
        Br = reshape(Br, size(X));
        Bz = reshape(B_data.Bz, size(X));

        surf(X*1e3, Y*1e3, -ones(size(X)), Bz*1e3, 'FaceAlpha', 0.9);
        shading interp; colorbar;
        clim([0, Bz_max*1e3]);
        title(sprintf('Bz (mT) at z = %.2f mm', Z_low(1,1)*1000));
        xlabel('X (mm)'); ylabel('Y (mm)');
        hold on;
        % Add circle at center with diameter coneDtop
        theta_circle = linspace(0, 2*pi, 100);
        x_circle = (coneDtop/2) * cos(theta_circle) * 1e3;
        y_circle = (coneDtop/2) * sin(theta_circle) * 1e3;
        plot3(x_circle, y_circle, ones(size(x_circle))*10, 'r-', 'LineWidth', 2);
        hold off;
        view(2);
        xlim([-coneDbottom/2*1e3, coneDbottom/2*1e3]);
        ylim([-coneDbottom/2*1e3, coneDbottom/2*1e3]);


        subplot(1,2,2);
        surf(X*1e3, Y*1e3, -ones(size(X)), Br*1e3, 'FaceAlpha', 0.9);
        shading interp; colorbar;
        clim([0, Br_max*1e3]);
        title(sprintf('Br (mT) at z = %.2f mm', Z_low(1,1)*1000));
        xlabel('X (mm)'); ylabel('Y (mm)');
        hold on;
        % Add circle at center with diameter coneDtop
        plot3(x_circle, y_circle, ones(size(x_circle))*10, 'r-', 'LineWidth', 2);
        hold off;
        view(2);
        xlim([-coneDbottom/2*1e3, coneDbottom/2*1e3]);
        ylim([-coneDbottom/2*1e3, coneDbottom/2*1e3]);

        % z=2e-3
        figure;
        subplot(1,2,1);
        Z_high = 2e-3 * ones(size(X));
        B_data = interpolateMagneticFlux(result, X, Y, Z_high);
        Bz = reshape(B_data.Bz, size(X));
        Br = reshape(sqrt(B_data.Bx.^2+B_data.By.^2), size(X));
        surf(X*1e3, Y*1e3, -ones(size(X)), Bz*1e3, 'FaceAlpha', 0.9);
        shading interp; colorbar;
        clim([0, Bz_max*1e3]);
        title(sprintf('Bz (mT) at z = %.2f mm', Z_high(1,1)*1000));
        xlabel('X (mm)'); ylabel('Y (mm)');
        hold on;
        % Add circle at center with diameter coneDtop
        plot3(x_circle, y_circle, ones(size(x_circle))*10, 'r-', 'LineWidth', 2);
        hold off;
        view(2);
        xlim([-coneDbottom/2*1e3, coneDbottom/2*1e3]);
        ylim([-coneDbottom/2*1e3, coneDbottom/2*1e3]);

        subplot(1,2,2);
        surf(X*1e3, Y*1e3, -ones(size(X)), Br*1e3, 'FaceAlpha', 0.9);
        shading interp; colorbar;
        clim([0, Br_max*1e3]);
        title(sprintf('Br (mT) at z = %.2f mm', Z_high(1,1)*1000));
        xlabel('X (mm)'); ylabel('Y (mm)');
        hold on;
        % Add circle at center with diameter coneDtop
        plot3(x_circle, y_circle, ones(size(x_circle))*10, 'r-', 'LineWidth', 2);
        hold off;
        view(2);
        xlim([-coneDbottom/2*1e3, coneDbottom/2*1e3]);
        ylim([-coneDbottom/2*1e3, coneDbottom/2*1e3]);

        %% x=0 intersection plot
        figure;
        y_slice = linspace(-coneDbottom, coneDbottom, 80);
        z_slice = linspace(-(coneHeight+coreHeight), 10e-3, 100);
        [Y_slice, Z_slice] = meshgrid(y_slice, z_slice);
        X_slice = zeros(size(Y_slice));

        % Interpolate field at x=0 plane
        B_slice = interpolateMagneticFlux(result, X_slice, Y_slice, Z_slice);
        B_mag_slice = sqrt(B_slice.Bx.^2 + B_slice.By.^2 + B_slice.Bz.^2);
        B_mag_slice = reshape(B_mag_slice, size(Y_slice));

        % Plot field magnitude as colormap
        contourf(Y_slice*1e3, Z_slice*1e3, B_mag_slice*1e3, 20, 'LineStyle', 'none');
        hold on;
        colorbar;
        xlabel('Y (mm)'); ylabel('Z (mm)');
        title('Field at X=0 Plane');
        axis equal tight;

        % Add sparsely sampled field direction arrows (in-plane: By, Bz)
        skip = 4;  % Sample every 5th point
        By = reshape(B_slice.By, size(Y_slice));
        Bz = reshape(B_slice.Bz, size(Y_slice));
        quiver(Y_slice(1:skip:end, 1:skip:end)*1e3, Z_slice(1:skip:end, 1:skip:end)*1e3, ...
               By(1:skip:end, 1:skip:end), Bz(1:skip:end, 1:skip:end), ...
               'k', 'LineWidth', 1.2, 'AutoScale', 'on', 'AutoScaleFactor', 1.5);

        % Draw cone/core contour
        % Cone profile
        z_cone = [0, -coneHeight];
        r_cone_profile = [coneDtop/2, coneDbottom/2];
        plot(r_cone_profile*1e3, z_cone*1e3, 'r-', 'LineWidth', 2);
        plot(-r_cone_profile*1e3, z_cone*1e3, 'r-', 'LineWidth', 2); 

        plot([-coneDbottom/2, coneDbottom/2], [-coneHeight, -coneHeight], 'r-', 'LineWidth', 2);

        % Core profile (cylinder)
        z_core = [-coneHeight-1e-4, -(coneHeight+coreHeight+1e-4)];
        r_core = [coneDbottom/2, coneDbottom/2];
        plot(r_core*1e3, z_core*1e3, 'r-', 'LineWidth', 2); 
        plot(-r_core*1e3, z_core*1e3, 'r-', 'LineWidth', 2); 

        % Core top
        plot([-(coneDbottom/2), (coneDbottom/2)]*1e3, ...
             [-(coneHeight+1e-4), -(coneHeight+1e-4)]*1e3, 'r-', 'LineWidth', 2);

        % Core bottom
        plot([-(coneDbottom/2), (coneDbottom/2)]*1e3, ...
             [-(coneHeight+coreHeight+1e-4), -(coneHeight+coreHeight+1e-4)]*1e3, 'r-', 'LineWidth', 2);

        hold off;
        legend('Field Magnitude', 'Field Direction', 'Cone/Core/Coil', 'Location', 'best');
    end
end

function create_cone_stl(filename, top_diameter, bottom_diameter, cone_height)
    % Create a truncated cone + cylindrical core as STL
    
    % Parameters
    n_sides = 32;  % Number of sides for approximation
    theta = linspace(0, 2*pi, n_sides+1);
    
    % Cone vertices
    r_top = top_diameter / 2;
    r_bottom = bottom_diameter / 2;
    
    % Top circle (z = 0)
    x_top = r_top * cos(theta);
    y_top = r_top * sin(theta);
    z_top = zeros(size(theta));
    
    % Bottom circle of cone (z = -cone_height)
    x_bottom = r_bottom * cos(theta);
    y_bottom = r_bottom * sin(theta);
    z_bottom = -cone_height * ones(size(theta));
    
    % Open file
    fid = fopen(filename, 'w');
    fprintf(fid, 'solid cone\n');
    
    % Write cone side faces
    for i = 1:n_sides
        % each stripe is divided into two triangles
        % Triangle 1
        fprintf(fid, '  facet normal 0 0 0\n');
        fprintf(fid, '    outer loop\n');
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', x_top(i), y_top(i), z_top(i));
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', x_bottom(i), y_bottom(i), z_bottom(i));
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', x_bottom(i+1), y_bottom(i+1), z_bottom(i+1));
        fprintf(fid, '    endloop\n');
        fprintf(fid, '  endfacet\n');
        
        % Triangle 2
        fprintf(fid, '  facet normal 0 0 0\n');
        fprintf(fid, '    outer loop\n');
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', x_top(i), y_top(i), z_top(i));
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', x_bottom(i+1), y_bottom(i+1), z_bottom(i+1));
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', x_top(i+1), y_top(i+1), z_top(i+1));
        fprintf(fid, '    endloop\n');
        fprintf(fid, '  endfacet\n');
    end
    
    % Top cap (triangular fan)
    for i = 1:n_sides
        fprintf(fid, '  facet normal 0 0 1\n');
        fprintf(fid, '    outer loop\n');
        fprintf(fid, '      vertex 0 0 0\n');
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', x_top(i), y_top(i), z_top(i));
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', x_top(i+1), y_top(i+1), z_top(i+1));
        fprintf(fid, '    endloop\n');
        fprintf(fid, '  endfacet\n');
    end
    
    % Bottom cap
    for i = 1:n_sides
        fprintf(fid, '  facet normal 0 0 -1\n');
        fprintf(fid, '    outer loop\n');
        fprintf(fid, '      vertex 0 0 %.6e\n', z_bottom(1));
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', x_bottom(i+1), y_bottom(i+1), z_bottom(i+1));
        fprintf(fid, '      vertex %.6e %.6e %.6e\n', x_bottom(i), y_bottom(i), z_bottom(i));
        fprintf(fid, '    endloop\n');
        fprintf(fid, '  endfacet\n');
    end
    
    fprintf(fid, 'endsolid cone\n');
    fclose(fid);
    
    fprintf('STL file created: %s\n', filename);
end