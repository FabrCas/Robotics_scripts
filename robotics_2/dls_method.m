%to test 

%{
instead of using jabobian pseudoinverse for q dot, which both reduce the
velocity norm and the error (if exists, so not totally feasible), we can
use DLS method which minimize the error with a tolerance and minimize the
joint velocity norm.
this method is more robust closer the singularities of J (rank leak), price
to pay is that we have always an error.

lambda -> integer >= 0, if == 0 (deactivated) you just minimze the norm of the error.
            so this value is a weight to give more importance to the
            velocity norm.
            suitable choice of lambda, using SVD, is choose according to
            minimum singular value sigma(q).

we can think of deactive DLS when we are "enough" far from singularities.

lambda or mu as symbolic definition of the formula
%}

function J_dls = dls_method(J,lambda)
    m = size(J,1);
    I = eye(m,m);
    J_dls = transpose(J)*inv( (lambda^2*I)+(J*transpose(J)) );
end 