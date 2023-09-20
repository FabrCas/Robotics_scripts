%% exercise first part

% April 14, 2021 exercise 1

syms q1 q1d l1 rc1x rc1y real

pc1_1 = [rc1x, rc1y, 0, 1]';

% transformation matrix

DH_params =[[0],[l1],[0],[q1]]; %alpha,a,d,theta
[T, A] = DHMatrix(DH_params);
disp(T)

pc1_0 = T * pc1_1;
pc1_0 = pc1_0(1:3);
disp(pc1_0)

%% 

pc1_0d = diff(pc1_0,q1)*q1d;
disp(pc1_0d)

vc1_squared_norm = simplify(pc1_0d'*pc1_0d, 15);
disp(vc1_squared_norm)

%% reset everything
clear all
close all
clc
%% defintion of symbolic variables

% kinematic variables
syms q1 q2 real                         % q
syms q1d q2d real                       % q dot 
syms q1dd q2dd real                     % q dot dot

% dynamic variables
syms I1zz I2zz real                     % baricentral inertia elements  (also I1 I2 ... generic is fine)
syms l real                             % link length
syms m1 m2 real                         % link masses

q = [q1 q2]';                           % vector of joint variables
q_dot = [q1d q2d]';                     % vector of joint velocities
q_dotdot = [q1dd, q2dd]';               % vector of joint accelerations
%% DH parameters 

% i define the DH parameters
alpha = [0, 0]';
a = [l, l]';
d = [0, 0]';
theta = [q1, q2]';
dhtable = [alpha, a, d, theta];   %  [[alpha]' [a]' [d]' [theta]']

% actually is not required the computation of the transformation matrix
if false
    [T, A] = DHMatrix(dhtable);
    disp("transformation matrix")
    disp(T)
end
% 

%% define the CoM position in the same RF (i.e. CoM 1 in RF1)
% even if planar is convenient to define all in 3D
pc1_1 = [-l/2, 0, 0]';
pc2_2 = [-l/2, 0, 0]';


%% % computation of the kinetic energy
% we define the linear and angular velocity for the base (fixed) 
v0_0 = [0 0 0]';
omega0_0 = [0 0 0]';
string_joints_types = 'rr';
array_q = q;
array_q_dot = q_dot;
array_m_i = [m1 m2];
cell_ri_ci = {pc1_1, pc2_2};
array_I_ci = [I1zz, I2zz];


[T, array_Ti] = compute_kinetic_energy_T_Ti(string_joints_types, dhtable, array_q, array_q_dot, cell_ri_ci, array_m_i, array_I_ci, v0_0, omega0_0, true);
disp("Full Kinetic energy")
disp(T)
%% M(q)

M = T_to_M(T, q_dot);
disp("Inertia matrix")
disp(M)
%% dynamic coefficients for M    (not mandatory)
syms a1 a2 a3 real 

a1_exp = (m2*l^2)/4 + I2zz;
a2_exp = l^2*m2;
a3_exp = (l^2*m1)/4 + I1zz + I2zz;


M = subs(M,[a1_exp,a2_exp,a3_exp],[a1,a2,a3]);
disp("M(q) after parametrization via dynamic coefficients");
disp(M)

%% Coriolis and Centrifugal

c = M_to_C(M, q, q_dot);
disp("Coriolis and Centrifugal term");
disp(c)

%% Potential energy computation
syms g0 real

%define the gravity vector
g_vector = [g0 0 0];

%position vector in RF0 of CoM (since gravity term is just in x coordinate
%here y and z position can be whatever, so i left zero (but actually is
%not)
pc1_0 = [l*cos(q1), 0, 0]';   
pc2_0 = [l*cos(q1) + l*cos(q1+q2), 0, 0]';

array_m_i = [m1 m2];
cell_r0_ci ={pc1_0, pc2_0};

[U, array_Ui] = compute_potential_energy(array_m_i, cell_r0_ci, g_vector, true);
disp("Potential energy U(q)");
disp(U)
%% gravity term

g_q = U_to_g(U, q);
disp("Gravity term g(q)");
disp(g_q);
%% dynamic coefficients for g(q)    (not mandatory)
% can be useful to expand the gravity term for substitution
%g_q = expand(g_q);

%first substitution of coefficient from M
g_q = subs(g_q,[a1_exp,a2_exp,a3_exp],[a1,a2,a3]);
disp("")
disp(g_q)
disp("")

% define new dynamic coefficients
syms a4 a5 real
a4_exp = g0*l*m1;
a5_exp = g0*l*m2;

g_q = subs(g_q,[a4_exp,a5_exp],[a4,a5]);
disp("g(q) after parametrization via dynamic coefficients");
disp(g_q)
%% define the whole dynamic model

dyn_model = M*q_dotdot + c + g_q;

disp("Dynamic model lsh");
disp(dyn_model)
%% 
syms q1 q2 l real
syms dq1 dq2 real

f_q = [l*cos(q1) + l*cos(q1 +q2); l*sin(q1)+ l*sin(q1+q2)];  %direct kinematics as a column vector
disp(f_q)

q = [q1 q2]';
qdot = [dq1 dq2]';

J_q = jacobian(f_q, q);
disp(J_q)

J_dot = diff_J(J_q, q, qdot);
disp(J_dot)








