%{
this function uses the jacobian pseudoinverse + the projection in the null
space of q0_dot (correction) to compute q_dot.

q0_dot is the preferred joint velocity.

correction direction is orthognoal to the one with minimum norm (jacobian
based), so obviously the final solution is no more the one with minimum
norm.
freedom to choose q0_dot, to have a solution with other nice properties.
%}

function q_dot = null_space_method(J, r_dot, q0_dot)
    I = eye(length(J));
    projector = (I-pinv(J)*J);
    q_dot = pinv(J)*r_dot+ projector*q0_dot;
    q_dot = simplify(q_dot);
end


% an alternative version 
%{
function q_dot = null_space_method2(J, r_dot, q0_dot)
    q_dot = q0_dot + pinv(J)*(r_dot - J*q0_dot);
    q_dot = simplify(q_dot);
end
%}



