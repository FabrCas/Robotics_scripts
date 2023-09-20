function R = R2homo(R_e)
    R = [R_e(1,:) 0; R_e(2,:) 0; R_e(3,:) 0; 0 0 0 1];
end