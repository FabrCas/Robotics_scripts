function J_dls = dls_method(J,lambda)
    m = size(J,1);
    I = eye(m,m);
    J_dls = transpose(J)*inv( (lambda*I)+(J*transpose(J)) );
end 