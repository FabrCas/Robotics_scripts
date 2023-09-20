% phi dot is not equal omega but we have the following relation omega = T(phi) phi_dot
% this function is used to find T(phi)
% the parameters are 2 vectors: omega 1x3 and angle dot 1x3 usually
% containing symbolic expression of angle dot.

function res = decomposition_omega(omega, angles_dot)
    res = sym(zeros(length(angles_dot)));
    for i=1:length(angles_dot)
        for j=1:length(angles_dot)
            res(i,j) = diff(omega(i), angles_dot(j));
        end
    end
end
