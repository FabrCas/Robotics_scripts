
% 
% torque by the dynamic model like (M*qdd + c
% a -> vector of dynamic coefficients

function Y = u_to_Y(u,a)
    [Y, b] = equationsToMatrix(u,a);
end