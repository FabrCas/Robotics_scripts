% U -> same of the potential energy for all the links
% q_in -> joint variables used for the differentiation
function g_q = U_to_g(U, q_in)
    g_q = simplify(transpose(jacobian(U,q_in)));
end