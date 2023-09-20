% function used to check the property M dot - 2S is skew


function isSkewSymmetric = check_factorization(M, q, qdot, S, print_info)

% compute M dot
n = size(M,1);

% initialize dM
dM = zeros(n,n);

for i=1:n
    dM = dM + simplify(diff(M,q(i))*qdot(i));
end 

if print_info
    disp("derivate of inertia matrix")
    disp(dM)
end

test_matrix = simplify(expand(dM - 2*S),15);

if print_info
    disp("test_matrix => dM - 2S")
    disp(test_matrix)
end


% Check if the matrix is skew-symmetric
isSkewSymmetric = isequal(test_matrix, -test_matrix');

if print_info
    disp("check result:")
    disp(isSkewSymmetric)
end

end