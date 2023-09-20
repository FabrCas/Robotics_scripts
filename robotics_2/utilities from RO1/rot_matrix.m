% agle in radians !
% choose axis between: 'x','y','z'

function [matrix] = rot_matrix(angle ,axis)
%ROT_MATRIX Summary of this function goes here
%   Detailed explanation goes here
if axis == 'x'
    matrix = [1 0 0; 0 cos(angle) sin(angle); 0 -sin(angle) cos(angle)];
elseif axis == 'y'
    matrix = [cos(angle) 0 -sin(angle); 0 1 0; sin(angle) 0 cos(angle)];
elseif axis == 'z'
    matrix = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];
end

