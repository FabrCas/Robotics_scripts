%todo

function Y = Y_independent(Y, pi_vector)
    % Compute the reduced row echelon form
    rref_Y = rref(Y);
    disp(rref_Y)

    % Find the indices of the linearly independent columns
    dependent_columns = find(all(Y == 0, 1));
    disp(dependent_columns)
    
    independent_matrix = [];
    depedent_matrix = [];
    pi_independent = [];
    pi_dependent = [];

    for i=1:size(Y,2)
        if not any(dependent_columns(:) == i)
                

        end
    
    end


end