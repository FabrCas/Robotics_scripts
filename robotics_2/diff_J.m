% simple function for the differentiation of J

function J_dot = diff_J(J, q, qdot)
    % compute n
    n = length(q);

    % initialize J_dot
    J_dot = zeros(size(J,1),size(J,2));

    for i=1:n
        J_dot = J_dot + simplify(diff(J,q(i))*qdot(i));
    end 

    J_dot = simplify(J_dot,15);
end