function det = det3(matrix, useRow, n)
    if not(size(matrix,1) == size(matrix,2) )
        error("mica è quadrata eh")
    end 
    det = 0;
    if useRow
        for i=1:3
            det = det + matrix(n,i)*(-1)^(n+i)* det2(reduceMatrix(matrix,n,i));
        end 
    else
        for i=1:3
            det = det + matrix(i,n)*(-1)^(n+i)* det2(reduceMatrix(matrix,i,n));
        end 
    end 
end 


function det2 = det2(matrix)
    det2 = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1);
end 


function mat = reduceMatrix(matrix,m,n)
    mat = matrix;
    mat(:,n) = [];
    mat(m,:) = [];
end 