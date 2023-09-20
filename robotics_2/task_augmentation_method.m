%to test 

%{
JA_k -> cell array with the jacobians sapareted by priority 
        i.e. JK  ={[J1,J2] , J3}
q_in -> array of symbolic joint variables
rd_k -> r dot, cell array with the task velocities sapareted by priority
%}

function qd_k = task_augmentation_method(JA_k, q_in, rd_k)

n_k = size(JA_k,2);
disp("number of priorities for the tasks -> " + num2str(n_k))

% initialization
qd_0 = zeros(1, length(q_in))';
PA_0 =  eye(length(q_in));

qd_k = qd_0;
PA_k = PA_0;

for k = 1:n_k
    % each JAk{k} is m_ki x N, with m_ki the task associted to the priority
    % level k
    disp("J" + num2str(k) + ":")
    disp(JA_k{k})
    
    % update step for q dot
    qd_k = qd_k + pinv(JA_k{k}*PA_k) * (rd_k{k} - JA_k{k}* qd_k);

    % update step for projection matrix
    PA_k = PA_k - pinv(JA_k{k}*PA_k) * (JA_k{k}*PA_k);
end


    % I = eye(length(Jak{k}));
    % Pak = (I-pinv(Jak{k})*Jak{k});
