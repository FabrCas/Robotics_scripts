%% different from christoffel_symbol.m, this function returns The lower c multiplied by qdot^2 
% in fact we have the following relation: c(q, qdot) = qdot^T * C(q) * qdot 
% M -> the inertia matrix
% q -> the list of joint variable
% q_dot = vector (column format) of time derivative for joint variables. 
function C = M_to_C(M, q, q_dot)
    C = sym([]);
    n = length(q);
    
    for i=1:n
        M_k = M(:, i);
        c_k = 0.5*(jacobian(M_k, q)+(jacobian(M_k, q)') - diff(M, q(i)));
        %c_k = simplify(c_k);
        c_k = simplify(expand(c_k),5);
        C = [C; simplify((q_dot')*c_k*q_dot)];
    end
end
