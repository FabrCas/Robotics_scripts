% if the matrix is not full rank i cannot do the pseudoinverse
% this is used to minimize the weighted norm

% i.e. of weight matrix is the inertia one M(q)

function inverse = wpinv(matrix, weight_matrix)
    n = size(matrix,2);
    m = size(matrix, 1);

    if n>m
        inverse = simplify(inv(weight_matrix)*transpose(matrix)*inv(matrix*inv(weight_matrix)*transpose(matrix)));
    elseif m>n
        disp("not found expression for the case m<n")
    else 
        disp("matrix is squared used normal inverse function")
    end 
end