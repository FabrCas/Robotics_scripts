%{
H -> optimization function H(q)
q_in -> joint variables in H(q)

example of H(q) used -> sin(q2)^2 + sin(q3)^3

%}

% gradient is for a single function, while Jacobian for multiple functions
% like in the case of direct kinematics.

function q0_dot = projected_gradient(H,q_in)
    q0_dot = gradient(H,q_in);
end

% combine this function with null space method for the complete results of
% qdot