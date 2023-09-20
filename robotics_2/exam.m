%% ------------------------------ ex 1 ------------------------------------
clear all
close all
clc

% kinematic variables
syms q1 q2 q3 q4 real                         % q
syms q1d q2d q3d q4d real                       % q dot 
syms q1dd q2dd q3dd q4dd real                     % q dot dot

% dynamic variables
syms I1 I2 I3 I4 real                     % baricentral inertia elements  (also I1 I2 ... generic is fine)
syms l real                             % link length
syms m1 m2 real                         % link masses

q = [q1 q2 q3 q4]';                           % vector of joint variables
q_dot = [q1d q2d q3d q4d]';                     % vector of joint velocities
q_dotdot = [q1dd, q2dd q3dd q4dd]';               % vector of joint accelerations

%%

% i define the DH parameters
alpha = [pi/2, -pi/2, pi/2, 0]';
a = [0, 0, 0, 0]';
d = [q1, q2, q3,q4]';
theta = [0,0,0,0]';

dhtable = [alpha, a, d, theta];
%% 

syms dc2 dc4 g0 real 
syms m1 m2 m3 m4 real

rc1y = 0;
rc2y = q2 + dc2 ;
rc3y = q2;
rc4y = q2+ q4 +dc4;

U1 = 0;
U2 = -m2*g0*rc2y;
U3 = -m3*g0*rc3y;
U4 = -m4*g0*rc4y;


U = U1 +U2 + U3 +U4 ;
disp(U)
q_in = [q1 q2 q3 q4]';


g_q = simplify(transpose(jacobian(U,q_in)));
disp(g_q)

%% 
M = [m1 + m2 + m3 + m4, 0, m3 + m4, 0;0, m2+m3+m4, 0, m4;m3+m4, 0, m3+m4, 0;0, m4,0,m4];
disp(M)

q = [q1 q2 q3 q4]';                
q_dot = [q1d q2d q3d q4d]';  

C = M_to_C(M, q, q_dot)

dyn_model = M*q_dotdot+ g_q;

disp("Dynamic model lsh");
disp(dyn_model)
%% 

J = [1 0 1 0; 0 1 0 1];
J_T = J';
disp(J)
disp(J_T)
disp(M)

M_inv = inv(M);
disp(M_inv)

% cartesian inertia martix
M_c = inv(J*M_inv*J_T);
disp(M_c)
%% 

M_c_val = subs(M_c,[m3 m4],[1 1]);
disp(M_c_val)

%% 
r_dd = [0 0]';

JM = J * M_inv;
disp(JM)

JM_inv = pinv(JM);
disp(JM_inv)

tau0 = JM_inv * r_dd
disp(tau0)

%% 

r_dd = [2 -3]';
JM = J * M_inv;
disp(JM)

JM_inv = pinv(JM);
disp(JM_inv)

JM_inv_val = subs(JM_inv,[m1 m2 m3 m4],[1 1 1 1]);
disp(JM_inv_val)

%% 
g_val = [0 -29.4 0 -9.8]';
JM_val = subs(JM, [m1 m2 m3 m4],[1 1 1 1])
second_term = r_dd- JM_val*g_val
%% 
tau = JM_inv_val * second_term


%% ------------------------------ ex 2 ------------------------------------
clear all
close all
clc



%% 
clear all
close all
clc



%% ------------------------------ ex 3 ------------------------------------
clear all
close all
clc




%% ------------------------------ ex 4 ------------------------------------
clear all
close all
clc




%% ------------------------------ ex 5 ------------------------------------
clear all
close all
clc





%% ------------------------------ ex 6 ------------------------------------
clear all
close all
clc



