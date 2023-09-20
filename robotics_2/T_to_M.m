function M = T_to_M(T, q_dot)
    % could be useful to colled q dot? 
    r = length(q_dot);
    M = sym(zeros(r));
    c = r;
    for i=1:r
        for j=1:c
            % differentiate to extract the elements, on the diagonal we
            % have a second order derivative
            tmp = simplify(diff(T, q_dot(i)));
            M(i, j) = simplify(diff(tmp, q_dot(j)), 10);
        end
    end
end