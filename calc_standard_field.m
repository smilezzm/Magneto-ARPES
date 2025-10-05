function standard_field = calc_standard_field(z_target, doPlot)
%CALC_STANDARD_FIELD Build PDE model, solve magnetostatics at 0.2 A, and sample B on a grid.
%   [geometry, standard_field] = calc_standard_field(z_target, doPlot)
%   Outputs:
%     - geometry: PDE geometry used in the model (model.Geometry)
%     - standard_field: struct with fields X, Y, Z, Bx, By, Bz
%       The range of the grid X,Y,Z is:
%           x = linspace(-coilDout, coilDout, 100); 
%           y = x;
%           z = linspace(-0.5 * coilHeight, z_target, 100)
%   Args:
%     - z_target: target plane z (m), upper bound of sampling grid. It'd
%       better be larger than the expected upper bound of real field, since
%       the standard field is to be translated/rotated/scaled to get the real
%       field.
%     - doPlot (optional): true to plot geometry and field; default false

    if nargin < 2
        doPlot = false;
    end

    % Geometry parameters
    airHeight     = 0.05;
    airD          = 0.05;
    coilDin       = 6.3e-3;
    coilDout      = 12e-3;
    coilHeight    = 5.6e-3;
    shieldDout    = 13e-3;
    shieldDin     = 4e-3;
    shieldHeight  = 0.3e-3;
    coreDout1     = 5.334e-3;
    coreDin       = 2.286e-3;
    coreDout2     = 3.5e-3;
    coreHeight1   = 2.794e-3;
    coreHeight2   = 3.048e-3;
    gapCoilShield = 0.3e-3;
    gapCoilCore   = 1.28e-3;
    grooveLength  = 1.016e-3;
    wireD         = 0.165e-3;
    current = 0.2;               % A

    % No rotation/translation (standard configuration)
    thetax = 0; thetay = 0; thetaz = 0;
    translationX = 0; translationY = 0; translationZ = 0;

    % Build geometry
    shieldGm = multicylinder([shieldDin/2, shieldDout/2], shieldHeight, Void=[1 0]);
    shieldGm = translate(shieldGm, [0,0,-shieldHeight+gapCoilShield]);

    coilGm = multicylinder([coilDin/2, coilDout/2], coilHeight, Void=[1 0]);
    coilGm = translate(coilGm, [0,0,-shieldHeight-coilHeight]);

    airGm = multicylinder(airD/2, airHeight);
    airGm = translate(airGm, [0,0,-airHeight/2-coilHeight/2]);

    coreGm1 = multicylinder([coreDin/2, coreDout1/2], coreHeight1, Void=[1 0]);
    coreGm2 = multicylinder([coreDin/2, coreDout2/2], coreHeight2, Void=[1 0]);
    coreGm1 = translate(coreGm1, [0,0,-shieldHeight-coreHeight1-gapCoilCore]);
    coreGm2 = translate(coreGm2, [0,0,-shieldHeight-coreHeight1-coreHeight2-gapCoilCore+2e-5]);
    cuttingBox = multicuboid(grooveLength, coreDout1, grooveLength);
    cuttingBox = translate(cuttingBox, [0, 0, -shieldHeight-gapCoilCore-grooveLength/2]);

    coreGm1 = fegeometry(coreGm1);
    coreGm2 = fegeometry(coreGm2);
    shieldGm = fegeometry(shieldGm);
    coilGm   = fegeometry(coilGm);
    cuttingBox = fegeometry(cuttingBox);
    airGm    = fegeometry(airGm);
    groovedCoreGm = subtract(union(coreGm1, coreGm2), cuttingBox);

    objs = {shieldGm, groovedCoreGm, coilGm};
    for k = 1:numel(objs)
        objs{k} = rotate(objs{k}, thetax*180/pi, [0, 0, 0], [1e-3, 0, 0]);
        objs{k} = rotate(objs{k}, thetay*180/pi, [0, 0, 0], [0, 1e-3, 0]);
        objs{k} = rotate(objs{k}, thetaz*180/pi, [0, 0, 0], [0, 0, 1e-3]);
        objs{k} = translate(objs{k}, [translationX, translationY, translationZ]);
    end
    gm = addCell(airGm, objs{1});
    gm = addCell(gm, objs{2});
    gm = addCell(gm, objs{3});
    fprintf('geometry built\n');

    if doPlot
        figure; pdegplot(gm, FaceAlpha=0.2, CellLabels="on");
        title('Geometry'); axis equal
    end

    % Magnetostatic model
    model = femodel(AnalysisType="magnetostatic", Geometry=gm);
    model.VacuumPermeability = 1.2566370614E-6;
    currentDensity = current / wireD^2; % A/m^2

    % Materials
    model.MaterialProperties = materialProperties(RelativePermeability=1);
    model.MaterialProperties(2) = materialProperties(RelativePermeability=83000);   % shield
    model.MaterialProperties(3) = materialProperties(RelativePermeability=100000);  % core

    % Coil current on Cell 4 (air=1, shield=2, core=3, coil=4 in this build)
    model.CellLoad(4) = cellLoad(CurrentDensity=@(region, state) ...
        windingCurrent3D(region, state, thetax, thetay, thetaz, ...
                         [translationX, translationY, translationZ], currentDensity));

    % Mesh (balanced for speed/accuracy)
    model = applyOptimizedMesh(model, shieldHeight, gapCoilShield, gapCoilCore, coilHeight, airD);
    fprintf('mesh built\n');

    % Solver options and BC
    model.SolverOptions.RelativeTolerance = 5e-5;
    model.SolverOptions.AbsoluteTolerance = 1e-7;
    model.FaceBC(1:3) = faceBC(MagneticPotential=[0;0;0]); % air cylinder outer faces

    % Solve
    fprintf('start solving the model\n');
    R = solve(model);
    fprintf('model solved\n');

    % Sample B on a 3D grid
    x = linspace(-coilDout, coilDout, 100);
    y = x;
    z = linspace(-0.5 * coilHeight, z_target, 100);
    [X, Y, Z] = ndgrid(x, y, z);
    Bdata = interpolateMagneticFlux(R, X, Y, Z); % struct with fields Bx, By, Bz

    % Shape to grid and package
    Bx = reshape(Bdata.Bx, size(X));
    By = reshape(Bdata.By, size(Y));
    Bz = reshape(Bdata.Bz, size(Z));

    standard_field = struct('X', X, 'Y', Y, 'Z', Z, 'Bx', Bx, 'By', By, 'Bz', Bz);
    save('standard_field.mat', 'X','Y','Z','Bx','By','Bz','current','z_target');

    if doPlot
        % Plot B field in 3D together with the geometry
        figure
        stepX = 4;
        X_sparse = X(1:stepX:end, 1:stepX:end, 1:stepX:end);
        Y_sparse = Y(1:stepX:end, 1:stepX:end, 1:stepX:end);
        Z_sparse = Z(1:stepX:end, 1:stepX:end, 1:stepX:end);
        Bx_sparse = Bx(1:stepX:end, 1:stepX:end, 1:stepX:end);
        By_sparse = By(1:stepX:end, 1:stepX:end, 1:stepX:end);
        Bz_sparse = Bz(1:stepX:end, 1:stepX:end, 1:stepX:end);
        quiver3(X_sparse,Y_sparse,Z_sparse,Bx_sparse,By_sparse,Bz_sparse,Color="r")
        hold on
        pdegplot(gm,FaceAlpha=0.2);
        xlim([-coilDout, coilDout]);
        ylim([-coilDout, coilDout]);
        zlim([-coilHeight, 3e-3]);
        
        
        % Plot B slice at y=0
        % --- Step 1: Field slice ---
        Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2);
        [~, iy] = min(abs(y));  % index of y closest to 0
        Xslice = squeeze(X(:, iy, :));  
        Zslice = squeeze(Z(:, iy, :));
        Bslice = squeeze(Bmag(:, iy, :));
        figure
        surf(Xslice, Zslice, -ones(size(Bslice)), Bslice, 'EdgeColor', 'none', 'FaceAlpha', 0.9)
        view(2)
        axis equal tight
        colormap('turbo')
        colorbar
        xlabel('X [m]')
        ylabel('Z [m]')
        title('|B| at y ≈ 0')
        hold on
        % --- Step 2: Add B-field vectors (projections onto slice plane) ---
        BxSlice = squeeze(Bx(:, iy, :));
        BzSlice = squeeze(Bz(:, iy, :));
        stepX = 2;
        stepZ = 2;
        Xq = Xslice(1:stepX:end, 1:stepZ:end);
        Zq = Zslice(1:stepX:end, 1:stepZ:end);
        Bxq = BxSlice(1:stepX:end, 1:stepZ:end);
        Bzq = BzSlice(1:stepX:end, 1:stepZ:end);
        quiver(Xq, Zq, Bxq, Bzq, 1.2, 'r', 'AutoScale', 'on', 'LineWidth', 1);
        hold on
        % --- Step 3: Geometry slice ---
        rectangle('Position', [-shieldDout/2, -shieldHeight+gapCoilShield, (shieldDout-shieldDin)/2, shieldHeight], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [shieldDin/2, -shieldHeight+gapCoilShield, (shieldDout-shieldDin)/2, shieldHeight], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [-coilDout/2, -shieldHeight-coilHeight, (coilDout-coilDin)/2, coilHeight], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [coilDin/2, -shieldHeight-coilHeight, (coilDout-coilDin)/2, coilHeight], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [-coreDout1/2, -shieldHeight-gapCoilCore-coreHeight1, (coreDout1-coreDin)/2, coreHeight1], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [coreDin/2, -shieldHeight-gapCoilCore-coreHeight1, (coreDout1-coreDin)/2, coreHeight1], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [-coreDout2/2, -shieldHeight-gapCoilCore-coreHeight1-coreHeight2, (coreDout2-coreDin)/2, coreHeight2], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [coreDin/2, -shieldHeight-gapCoilCore-coreHeight1-coreHeight2, (coreDout2-coreDin)/2, coreHeight2], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
    
        
        % Plot B slice at x=0
        % --- Step 1: Field slice ---
        [~, ix] = min(abs(x));  % index of x closest to 0
        Yslice = squeeze(Y(ix, :, :));  % Notice X (length(y),length(x), length(z)) !!!!
        Zslice = squeeze(Z(ix, :, :));
        Bslice = squeeze(Bmag(ix, :, :));
        figure
        surf(Yslice, Zslice, -ones(size(Bslice)), Bslice, 'EdgeColor', 'none', 'FaceAlpha', 0.9)
        view(2)
        axis equal tight
        colormap('turbo')
        colorbar
        xlabel('y [m]')
        ylabel('Z [m]')
        title('|B| at x ≈ 0')
        hold on
        % --- Step 2: Add B-field vectors (projections onto slice plane) ---
        BySlice = squeeze(By(ix, :, :));
        BzSlice = squeeze(Bz(ix, :, :));
        stepY = 2;
        stepZ = 2;
        Yq = Yslice(1:stepY:end, 1:stepZ:end);
        Zq = Zslice(1:stepY:end, 1:stepZ:end);
        Byq = BySlice(1:stepY:end, 1:stepZ:end);
        Bzq = BzSlice(1:stepY:end, 1:stepZ:end);
        quiver(Yq, Zq, Byq, Bzq, 1.2, 'r', 'AutoScale', 'on', 'LineWidth', 1);
        hold on
        % --- Step 3: Geometry slice ---
        rectangle('Position', [-shieldDout/2, -shieldHeight+gapCoilShield, (shieldDout-shieldDin)/2, shieldHeight], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [shieldDin/2, -shieldHeight+gapCoilShield, (shieldDout-shieldDin)/2, shieldHeight], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [-coilDout/2, -shieldHeight-coilHeight, (coilDout-coilDin)/2, coilHeight], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [coilDin/2, -shieldHeight-coilHeight, (coilDout-coilDin)/2, coilHeight], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [-coreDout1/2, -shieldHeight-gapCoilCore-coreHeight1, (coreDout1-coreDin)/2, coreHeight1], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [coreDin/2, -shieldHeight-gapCoilCore-coreHeight1, (coreDout1-coreDin)/2, coreHeight1], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [-coreDout2/2, -shieldHeight-gapCoilCore-coreHeight1-coreHeight2, (coreDout2-coreDin)/2, coreHeight2], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
        hold on
        rectangle('Position', [coreDin/2, -shieldHeight-gapCoilCore-coreHeight1-coreHeight2, (coreDout2-coreDin)/2, coreHeight2], 'EdgeColor', 'r', 'LineWidth', 1, 'Clipping', 'off');
    end
end

function f3D = windingCurrent3D(region, ~, thetax, thetay, thetaz, Tvec, currentDensity)
    % Don't forget to consider rotation and translation
    cx = cos(thetax); sx = sin(thetax);
    cy = cos(thetay); sy = sin(thetay);
    cz = cos(thetaz); sz = sin(thetaz);
    Rx = [1 0 0; 0 cx -sx; 0 sx cx];
    Ry = [cy 0 sy; 0 1 0; -sy 0 cy];
    Rz = [cz -sz 0; sz cz 0; 0 0 1];
    Rot = Rz*Ry*Rx;

    xg = region.x - Tvec(1);
    yg = region.y - Tvec(2);
    zg = region.z - Tvec(3);
    Pglob = [xg; yg; zg];
    Ploc = Rot.' * Pglob;

    xl = Ploc(1,:); yl = Ploc(2,:);
    phi = atan2(yl, xl);
    I_local = [-sin(phi); cos(phi); zeros(size(phi))];

    f3D = Rot * I_local * currentDensity;
end

function model = applyOptimizedMesh(model, shieldHeight, gapCoilShield, ...
                                        gapCoilCore, coilHeight, airD)
    %APPLYOPTIMIZEDMESH  Apply tuned mesh parameters to a PDE model
    %
    %   model = applyOptimizedMesh(model, shieldHeight, gapCoilShield, ...
    %                              gapCoilCore, coilHeight, airD)
    %
    %   This function encapsulates the mesh sizing logic so that intermediate
    %   variables do not clutter the base workspace.
    
    % --- Balance between accuracy and speed ---
    smallFeature   = min(shieldHeight, gapCoilShield);     
    hSmallTarget   = smallFeature / 3;       % Reduced from 5 (coarser mesh)
    shieldTarget   = hSmallTarget * 1.1;     % Slightly increase for faster calc
    HminVal        = 1.1 * hSmallTarget;     % Increased minimum element size
    coilGapTarget  = gapCoilCore / 3;        % Reduced from 5
    coilBulkTarget = coilHeight / 5;         % Reduced from 7
    airCoarse      = airD / 6;               % Coarser air mesh

    % Increased growth rate for faster transition between fine and coarse
    HmaxVal  = min(airCoarse, 22 * HminVal); % Increased from 16
    HgradVal = 1.8;                          % Increased from 1.65

    % Safety factor function
    sf = @(h) max(h, 1.05 * HminVal);

    % Face groups
    shieldFaces = cellFaces(model.Geometry, 2);
    coreFaces   = cellFaces(model.Geometry, 3);
    coilFaces   = cellFaces(model.Geometry, 4);

    % --- Generate mesh ---
    model = generateMesh(model, ...
        GeometricOrder = "linear", ...
        Hmin  = HminVal, ...
        Hmax  = HmaxVal, ...
        Hgrad = HgradVal, ...
        Hface = { ...
            shieldFaces, sf(shieldTarget), ...
            coilFaces,   sf(min(coilGapTarget, coilBulkTarget)), ...
            coreFaces,   sf(coilGapTarget)});

    fprintf('mesh created\n');
end