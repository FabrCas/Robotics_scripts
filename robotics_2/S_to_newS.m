% function used to find a different valid factorization that satisfy:
% M dot - 2S is skew matrix.
% to Works the starting S matrix should already satisfy this condition.
% qdot -> vector (column format) of joint velocities

function newS = S_to_newS(S,qdot)
    S_qdot = vector_to_skew(qdot);
    newS = simplify(S + S_qdot);
end