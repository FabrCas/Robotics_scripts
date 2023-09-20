function j_analytical = analytical_jacobian(f_r,q_in)
% f_r -> function for the direct kinematics, equations should be in columns
% rx = N*cos(q2)-q3*sin(q2);
% ry = N*sin(q2)+q3*cos(q2);
% rz = q1;
% f_r = [rx;ry;rz]

% q_in -> joint variables of the manipulator, should be in columns
    [temp ,numCols_f] = size(f_r)
    [temp, numCols_q] = size(q_in)
    if numCols_f > 1
        f_r = f_r.'
    end

    if numCols_q > 1
        q_in = q_in.'
    end

    j_analytical = jacobian(f_r,q_in)
end 