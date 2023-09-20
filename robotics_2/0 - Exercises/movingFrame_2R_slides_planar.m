% define the symbolic variables for the 2R robot under gravity (vertical
% plane)

syms q1 q2 real
syms d1 d2 real 
syms l1 l2 real
syms m1 m2 real 
syms Ic1_xx Ic1_yy Ic1_zz real  %symbolic variables explicitly defined as real
syms Ic2_xx Ic2_yy Ic2_zz real  %symbolic variables explicitly defined as real
% derivative of joint variables
syms q1d q2d real
syms q1dd q2dd real
% other terms of the dynamic model 
syms u1 u2 g0 real


% let's derive the DH parameters and matrices.

alpha = transpose([0, 0]);
a = transpose([l1, l2]);
d = transpose([0, 0]);
theta = transpose([q1, q2]);
DH = [alpha, a, d , theta];

% to use the function for the rotatiton matrices the parameters should be
% in this order [alpha a d theta]

[T,A] = DHMatrix(DH);   % T-> matrix, A -> cell array

disp("full transformation matrix by DH parameters")
disp(T)

% extraction of rotation matrices
R1 = A{1}(1:3,1:3);
R2 = A{2}(1:3,1:3);


% i get the direct kinematics from the full t
f_r = T(1:3,4);
disp(R1)
disp(R2)
disp(f_r)
%% kinetic energy
% i deifine the vectors going from link rf to com
rc1 = [-l1 + d1, 0 ,0]';
rc2 = [-l2 + d2, 0 ,0]';

% initialization of the moving frame algorithm

omega0_0 = [0, 0, 0]';
v0_0 = [0,0,0]';
print_info= true;

% i compute the kinetic energy via script
[T, array_Ti] = compute_kinetic_energy_T_Ti("RR", DH, [q1,q2], [q1d, q2d], {rc1, rc2}, [m1,m2], {Ic1_zz Ic2_zz}, v0_0, omega0_0, 1); % last variable is to print info

disp("Kinetic energy:")
disp(T)
%% inertia matrix

% next step is the derivation of the inertia matrix from the total kinetic
% energy from konig theorem

qdot = [q1d, q2d]';
M = T_to_M(T,qdot);
disp("Inertia matrix")
disp(M)
pretty(M(1,1))
pretty(M(1,2))
pretty(M(2,1))
pretty(M(2,2))


%% here is possible to do regrouping extracting the dynamic coefficients ai
a1_v = m1*d1^2 + m2*d2^2 + m2*l1^2 + Ic1_zz + Ic2_zz;
a2_v = 2*m2*d2*l1;
a3_v = m2*d2^2 + Ic2_zz;

syms a1 a2 a3 real
% substitute to get the new matrix

M_a = subs(M, a1_v, a1);
M_a = subs(M_a, a2_v, a2);
M_a = subs(M_a, a3_v, a3);
disp(M_a)
%% coriolis and centrifugal term exploting christoffel symbols
% Follows the computation of the Christoffel symbols for the coriolis and
% centrifugal term

q = [q1, q2]';
qdot = [q1d, q2d]';

c = M_to_C(M, q, qdot);
disp("coriolis and centrifugal term c(q, qdot)")
disp(c)

% even there i can substitute a2
%% now the computation for the gravity term (-mgh)

% first is needed the computation of vector going from base frame the each
% CoM by ispection
r0_c1 = [d1 * cos(q1), d1 * sin(q1), 0]';
r0_c2 = [l1 * cos(q1) + d2*cos(q1+q2), l1*sin(q1) + d2*sin(q1+q2), 0]';
%btw it's just really important the y component of these vectors since the
%gravity acceleration 
g = [0, -g0, 0]';
% then U1 and U2 are strightforward

U1 = - m1 * transpose(g) * r0_c1;   %alternatively U1 = m1 * -g0 * r0_c1(1)
U2 = - m2 * transpose(g) * r0_c2;
U = U1 + U2;
disp(U1)
disp(U2)
disp("potential energy:")
disp(U);

% gravity term
g_q  = U_to_g(U,[q1 q2]');
disp("gravity term:")
disp(g_q)

%% (alternatively)

% or using directly the function for the potential energy
[U, Ui] = compute_potential_energy([m1,m2], {r0_c1, r0_c2}, g, true);

disp("potential energy:")
disp(U)

g_q = U_to_g(U, [q1 q2]);
disp("gravity term:")
disp(g_q)
%% 

% i can also substitute in the g term the following dynamic coefficients
a5 = m2 * d2 * g0;
a4 = d1*g0*m1*cos(q1) + l1*cos(q1)*g0*m2;  
%% 

tmp = g_q; %M*[q1dd q2dd]';
disp(tmp)
disp(size(tmp))
%% the dynamic model
dyn_model = M*[q1dd q2dd]'+ c + g_q == [u1, u2]';

disp("the dynamic model of the 2R robot vertical plane")
disp(dyn_model);









