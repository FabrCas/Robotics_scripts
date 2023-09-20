% alpha = [pi/2, 0,0];
% a=[0,a2,a3];
% d=[d1 0 0];
% theta=[q1,q2,q3];
% 
% joints = 'RRR'
% table=[alpha',a',d',theta']
% 
% [T, A] = DHMatrix(table);
% A0_1=A{1};
% A1_2=A{2};
% A2_3=A{3};
% 
% f_r = get_f_r(T);
% f_r_ = f_r(1:3) %hide last row
% 
% % manual computation
% fprintf('manual computation\n')
% 
% Jl = jacobian(f_r(1:3), [q1,q2,q3])
% 
% z_i = [0 0 1]';
% R0_1 = get_rotation_mat(A0_1);
% R1_2 = get_rotation_mat(A1_2);
% Ja = [z_i, R0_1*z_i, R0_1*R1_2*z_i]
% 
% % Using function directly
% fprintf('Using function directly\n')
% 
% [jl, ja] = geometric_jacobian(f_r, joints, [q1 q2 q3], table)

function [Jl, Ja] = geometric_jacobian(f_r, sequence, q_in, params)
% geoJ = geometric_jacobian() takes as inputs:
%   -f_r: mapping from joint space to cartesian space
%   -sequence: A string containing the sequence of 'r's and 'p's
%   -q_in: The list of the variables we want to derive wrt
%   -params: List of params needed to get the DHMatrix
% and outputs:
%   -geoJ: The resulting geometric jacobian

    [~, A] = DHMatrix(params);
    sequence = char(lower(sequence));
    n_dof = strlength(sequence);
    
    f_r = f_r(1:3, :); %extract p, to derivate 
    
    %to calculate Jl, we use analytic jacobian
    Jl = jacobian(f_r, q_in);
    
    Ja = sym(zeros(3, n_dof));
    
    cells = cell(1, n_dof);
    cells{1} = get_rotation_mat(A{1});
    for i = 2:n_dof
        cells{i} = cells{i-1}*get_rotation_mat(A{i});
    end
    
    for i = 1:n_dof
        c = sequence(i);
        
        if c ~= 'r' && c ~= 'p'
            disp("Sequence must be composed of only 'r's and 'p's")
            geoJ = -1;
            return
        end
        
        %if the joint is prismatic then the element in Ja is equals to zero
        %otherwise:
        %   Ja(:, i) = z_(i-1) = R0_1 * ... * R(i-2)_(i-1)*[0; 0; 1]
        if c == 'r'
            if i == 1
                Ja(:, i) = [0; 0; 1];
            else
                Ja(:, i) = cells{i-1}*[0; 0; 1]; %[0;0;1] to take only the z component
            end
        end
    end
    Jl = simplify(Jl);
    Ja = simplify(Ja);
    %geoJ = simplify([Jl; Ja]);

    end