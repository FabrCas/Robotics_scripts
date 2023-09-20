%script to find a suitable D(q) matrix that makes the augmented one
%[A(q)D(q)]' non-singular


%function [outputArg1,outputArg2] = red_dym_D(inputArg1,inputArg2)
%outputArg1 = inputArg1;
%outputArg2 = inputArg2;
%end

syms l1 q1 l2 q2 real
% Example 2x3 matrix
A = [1 1 1; -(l1*sin(q1) + l2* sin(q1+q2)), -l2*sin(q1+q2), 0];

disp(A)


% Check if the matrix is non-square
if size(A, 1) ~= size(A, 2)
    % Compute the null space of the matrix
    null_space = null(A);
    
    % Check if the null space is non-empty
    if ~isempty(null_space)
        % Extract a non-zero vector from the null space
        additional_row = simplify(null_space(:, 1)');
        disp(additional_row);
    else
        disp("No additional row found.");
    end
else
    disp("The matrix is already square.");
end

