% v can be either a vector of R3 or a scalar

function S = vector_to_skew(v)

% old version, no result for v scalar
%    x = vector(1);
%    y = vector(2);
%    z = vector(3);
%    skew_symmetric_matrix = [0 -z y; z 0 -x; -y x 0];
%end

% If V (1x1) then S =
%
%           | 0  -v |
%           | v   0 |
%
if length(v) == 3
    % SO(3) case
    S = [  0   -v(3)  v(2)
          v(3)  0    -v(1)
         -v(2) v(1)   0];
else
    % SO(2) case
    S = [0 -v; v 0];

end