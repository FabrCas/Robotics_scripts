clc
clear all 

%% ------------------- same exercise of dyn model using moving frame algorithm
syms q1 q2 q3 dc1 dc2 dc3 I1 I2 I3 l1 l2 l3 real
syms m1 m2 m3 real
syms q1d q2d q3d real

%% define dh parameters
alpha = [pi/2, -pi/2, 0];
a=[0, 0, l3];
d=[l1, q2, 0];
theta=[q1, pi/2, q3];

dhtable=[alpha',a',d',theta'];

%% define rci

% define this vector locally
r1_c1 = [0, -dc1, 0]';
r2_c2 = [0, dc2, 0]';
r3_c3 = [-l3+dc3, 0, 0]';


%% 

string_joints_types = "RPR";
array_q = [q1 q2 q3]';
array_q_dot = [q1d q2d q3d]';
array_m_i = [m1 m2 m3]';
array_I_ci = [I1 I2 I3]';
cell_ri_ci = {r1_c1, r2_c2, r3_c3};
v0_0     = [0,0,0]';
omega0_0 = [0,0,0]';
print_info = true;


[T, array_Ti] = compute_kinetic_energy_T_Ti(string_joints_types, dhtable, array_q, array_q_dot, cell_ri_ci, array_m_i, array_I_ci, v0_0, omega0_0, print_info)
