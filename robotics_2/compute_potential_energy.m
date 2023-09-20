% cell_r0_ci -> cell array containing all the r vector from RF0 origin to
%               each  CoM
% g_vector -> gravity vector expressed in referene frame 0, this should be
%               in column format

function [U, array_Ui] = compute_potential_energy(array_m_i, cell_r0_ci, g_vector, print_info)
    n_links = length(array_m_i); %excluding the base
    disp(n_links)
    % initialization
    U = 0;
    array_Ui = sym([]);
    
    if size(g_vector,1) == 1 % row vector
        disp("gravity vector should be transposed for a column format, transposing...")
        g_vector = g_vector';
    end


    disp(cell_r0_ci)

    for i=1:n_links
        Ui = - array_m_i(i) * transpose(g_vector) * cell_r0_ci{i};
        array_Ui = [array_Ui, Ui];
        U = simplify(U + Ui);

        if print_info
            fprintf("Link %d------------------------------------\n", i);
            display(Ui);
            fprintf("\n---------------------------------------------\n\n");
        end
    end

end