function [vx0, vy0] = getInitialVelocity(vx_final, vy_final, v_total, invMapCell, v0Range)
    % v0Range is a row vector containing all values of |v| in the simulation
    % Find closest energy slice (|v| slice)
    vzElementNum = size(v0Range, 2);
    best_slice = 1;
    min_diff = inf;
    for z_slice = 1:vzElementNum
        diff = abs(v_total - v0Range(z_slice));
        if diff < min_diff
            min_diff = diff;
            best_slice = z_slice;
        end
    end
    
    if isempty(invMapCell{best_slice})
        vx0 = NaN;
        vy0 = NaN;
        return;
    end
    
    result = invMapCell{best_slice}(vx_final, vy_final);
    vx0 = result(1);
    vy0 = result(2);
end