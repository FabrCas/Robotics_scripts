function B = licols2(A)
    A = A.';
    N = size(A,1) ;                  % number of rows
    IncludeTF = false(N,1) ;         % by default, exclude all rows, except ...
    IncludeTF(1) = true ;            % first row which can always be included
    R0 = rank(A) ;                   % the original rank
    for k = 2:N,                     % loop over all rows
       B = A(IncludeTF,:) ;          % select the currently included rows of A
       IncludeTF(k) = rank(B) < R0 ; % include in B when the rank is less
    end
    isequal(rank(B), R0)             % check!
    B = B.';
end 