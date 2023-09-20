function [extreme_points,extreme_values,flection_points] = maxminVal_f(f,unknown)
% f-> function for the calculus of the max
% unknown -> symbolic unknown used to derivate respect it

%example:
% f = (pi*((6*t)/T^2 - (6*t^2)/T^3))/4
% unknown = t

% how to use:
% given a function f(x) the results describe 
% [1] the extreme values of x for which are obtained the max/min values
% [2] the extreme values of f(x) using as x the extreme values of [1]
% [3] inflection points: use this to determine whether a point is a max or
% a min.

% if the unkwnown found after the first derivative are zero, this means
% haveing a linear fuction, so 

%-------------------- first derivative -------------------------

fprintf("function:")
disp(f)

f2 = diff(f) == 0;
disp(f2)
extreme_points = solve(f2,unknown);

n_extreme_points = length(extreme_points);

extreme_values = subs(f, unknown, extreme_points);

%-------------------- second derivative -------------------------
f3 = diff(f,2);
flection_points = subs(f3, unknown, extreme_points);

% ------------------- return values ---------------------------
try
    extreme_points = eval(extreme_points);
    extreme_values = eval(extreme_values);
    flection_points = eval(flection_points);
catch
    try
        extreme_points = simplify(extreme_points);
    catch
        extreme_points = double(extreme_points);
    end
    try    
        extreme_values = simplify(extreme_values);
    catch
        extreme_values = double(extreme_values);
    end
    try
        flection_points = simplify(flection_points);
    catch
        flection_points = double(flection_points);
    end
end 
if n_extreme_points >= 1
    fprintf("extreme points: f'(x) == 0\n")
    disp(extreme_points)
    fprintf("extreme values in f(x) given the extreme points:\n")
    disp(extreme_values)
    fprintf("inflection points sign: (- abs max) (+ abs min) :\n")
    disp(flection_points)
else
    fprintf("linear function, max and min at the interval of the domain use subs with initial time and final time\n")
end
end
