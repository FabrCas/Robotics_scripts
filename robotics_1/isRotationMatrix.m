function res = isRotationMatrix(R)
    try 
        RR = eval(simplify(R'*R)); 
        det_R = eval(simplify(det(R))); 
    catch
        RR = R'*R; 
        det_R = det(R); 
    end 
    % orthonormal check
%     disp(det_R);
%     disp(RR);

    if isequal(round(RR,5), eye(3)) && isequal(round(det_R,5), 1)
        res = true;
    else
        res = false;
    end
end

