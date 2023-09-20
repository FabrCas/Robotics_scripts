%% exercise 2 
% i define the DH parameters
syms a1 a2 a3 real % actually are the same of l1 l2 l3
syms theta1 theta2 theta3 real % in this case are the joint variables


% kinematic parameters
syms l1 l2 l3 real
syms q1 q2 q3 real


% dynamic paramters
syms m1 m2 m3 real
syms rcx_1 rcy_1 rcx_2 rcy_2 rcx_3 rcy_3 real % planar value for the vector that locates the CoM
syms g0 real

% given the reference frame 0 with x axis pointing downward
g = [g0,0,0]'; 


% i define the full vector for the CoM location

r1_c1 = [rcx_1, rcy_1, 0, 1]';
r2_c2 = [rcx_2, rcy_2, 0, 1]';
r3_c3 = [rcx_3, rcy_3, 0, 1]';
%% 

% i compute the trasformation matrix

alpha = [0,0,0]';
a = [l1 l2 l3]';
d = [0,0,0]';
theta = [q1,q2,q3]';

DH_params = [alpha,a,d,theta];

[T, A] = DHMatrix(DH_params);

disp(A{1})
disp(A{2})
disp(A{3})

% the CoM vector in RF0 (actually we are interested only on the first
% element)
r0_c1 = simplify(A{1}* r1_c1);
disp(r0_c1)
r0_c2 = simplify(A{1}*A{2}* r2_c2);
disp(r0_c2)
r0_c3 = simplify(A{1}*A{2}*A{3}*r3_c3);
disp(r0_c3)


% the x-component
disp(r0_c1(1,:))
disp(r0_c2(1,:))
disp(r0_c3(1,:))
%% 
% we compute the potential energy

[U, array_Ui] = compute_potential_energy([m1 m2 m3], {r0_c1(1:3), r0_c2(1:3), r0_c3(1:3)}, g, true);

disp(U);

%% 
% we compute the gravity term
g_q = U_to_g(U, [q1, q2, q3]);
g_q = collect(g_q, g0);

%%
% now we want to have g_q = [0,0,0], let's start from the third equation
%eqn3 = rcy_3*cos(q1 + q2 + q3) + l3*sin(q1 + q2 + q3) + rcx_3*sin(q1 + q2 + q3) == 0;
eqn3 = collect(g_q(3), sin(q1 + q2 + q3)) == 0;
pretty(eqn3)

%sol = solve(eqn3, rcy_3,'Real',true);
%disp(sol)

% better to do this inspection by hand and we find
% rcy_3 = 0;
% rcx_3 = -l3

%% 
% second equation is 

%lseq2 = collect(g_q(2), g0);
%lseq2 = collect(lseq2, sin(q1 + q2));
%lseq2 = collect(lseq2, sin(q1 + q2 + q3));

eqn2 = g_q(2) == 0;
% i substitute the results from the first equation
eqn2 = subs(eqn2,rcy_3,0);
eqn2 = subs(eqn2,rcx_3,-l3); 
eqn2 = simplify(eqn2);

pretty(eqn2)

% finding:
% rcy_2 = 0
% rcx_2 = -((m2+m3)/m2) * l2
%% 
% third equation is 
eqn1 = g_q(1) ==  m1*g0*rcy_1*cos(q1);


disp(eqn1)

% i substitute the results from the first equation and second equation
eqn1 = subs(eqn1,rcy_3,0);
eqn1 = subs(eqn1,rcx_3,-l3); 
eqn1 = subs(eqn1,rcy_2,0);
eqn1 = subs(eqn1,rcx_2,-((m2+m3)/m2) * l2); 
eqn1 = simplify(eqn1);


pretty(eqn1)
 % finding final condition = m1*rcx_1 = -l1*(m1+m2+m3)

%% exercise 3

% DH params 
syms q1 q2 q3 q4 real
syms dq1 dq2 dq3 dq4 real

alpha = [pi/2, -pi/2, pi/2, 0]';
a = [0,0,0,0]';
d = [q1, q2, q3,q4]';
theta = [0,0,0,0]';

DH_params = [alpha,a,d,theta];
[T, A] = DHMatrix(DH_params);
disp(T)
disp(A{1})
disp(A{2})
disp(A{3})
disp(A{4})

%% 
% we compute v,vc and omega
v0_0 = [0 0 0]';
omega0_0 = [0 0 0]';

%rotation matrix for rf0 to world frame
matrix = rot_matrix(pi/2 ,'y');
disp(matrix)


syms r1_c1z r2_c2z r3_c3z r4_c4z real;
r1_c1 = [0 0 r1_c1z]';

R0_1 = affine_get_R(A{1});
r0_1 = affine_get_translation(A{1});


R0_1 = round(matrix*R0_1);
disp(R0_1)


[omega1_1, v1_c1, v1_1] = compute_omega_vc_v(r1_c1, r0_1, R0_1, dq1, omega0_0, v0_0, true);

disp(omega1_1)
disp(v1_c1)
disp(v1_1)


% i can continue in this way but since by inspection i can obtain vc_i i
% define
%% 

syms m1 m2 m3 m4 real
syms Ici_zz real

vc1 = [dq1 0 0]';
vc2 = [dq1 dq2 0]';
vc3 = [dq1 + dq3  dq2 0]';
vc4 = [dq1 + dq3  dq2 + dq4 0]';

omega_i = [0,0,0]';

% i compute all the kinetic energies

T1 = konig_theorem(m1, vc1, omega_i, Ici_zz);
disp(T1)
T2 = konig_theorem(m2, vc2, omega_i, Ici_zz);
disp(T2)
T3 = konig_theorem(m3, vc3, omega_i, Ici_zz);
disp(T3)
T4 = konig_theorem(m4, vc4, omega_i, Ici_zz);
disp(T4)

T = simplify(T1 + T2 + T3 +T4);
disp(T)
%% 

q_dot = [dq1, dq2, dq3, dq4];

M = T_to_M(T, q_dot);
disp(M)
%% 
% the direct kineamatics and the jacobian
px = q1 + q3;
py = q2+ q4;

f_r = [px, py];

J = jacobian(f_r, [q1 q2 q3 q4]);
disp(J)

%% 


syms vdx vdy real

vd = [vdx, vdy]';
% i compute q dot first by simpling using the pseudoinverse and then using
% the weighted pseudoinverse.

q_dot = pinv(J) * vd;
disp(q_dot)

%% 

% then the solution that tries to minimize also the kinetic energy involved
% using the inertia matrix as weight matrix

q_dot_w = wpinv(J,M) * vd;
disp(q_dot_w)

%% exercise 4

syms m1 m2 m3 real
syms l1 l2 l3 real
syms dc1 dc2 dc3 real
syms q1 q2 q3 real 
syms Ic1_yy Ic2_yy Ic3_zz real

% we define the position vector for the CoM

r1_c1 = [0 -dc1 0]';
r2_c2 = [0 dc2 0]';
r3_c3 = [-l3 + dc3 0 0]';

% get the DH table
alpha = [pi/2 -pi/2 0]';
a = [0 0 l3]';
d = [l1 q2 0]';
theta = [q1 pi/2 q3]';
DH_params = [alpha,a,d,theta];

[T, A] = DHMatrix(DH_params);
disp(T)

%% 
% computation of pc2_0 

pc2_2 = [0 dc2-q2 0]';

% rotation matrix from RF2 to RF0
T0_2 = A{1}*A{2};
R0_2 = T0_2(1:3,1:3);

disp(R0_2)
pc2_0 = R0_2* pc2_2;
disp(pc2_0)

% we differntiate respect to time, for doing this we introduce the t inside
% q
%% 

% code neeeded to perform derivation over time
syms t real;
q1_t = str2sym('q1(t)');
q2_t = str2sym('q2(t)');


R0_2 = subs(R0_2,[q1 q2], [q1_t q2_t]);
pc2_2 =  subs(pc2_2,[q1 q2], [q1_t q2_t]);
disp(R0_2);
disp(pc2_2);
%% 

pc2_0 = R0_2* pc2_2;
disp(pc2_0)
vc2_0 = simplify(diff(pc2_0,t));
disp(vc2_0)

norm_vc2 = simplify(expand(norm(vc2_0)^2))

%% 
% computation of pc3_0 

pc3_3 = [-l3 + dc3, 0 ,0, 1]';

pc3_0 = simplify(T * pc3_3)

%% 

syms q1d q2d q3d real

% try to solve this exercise using the moving frame algorithm and scripts
% how define array I_ci, i look to the reference frame
% for the first (RF1) the direction interested is the y1, so Ic1_yy, since
% link rotates around the y1 axis.
% for the second link the rotation is still the one of the previous one,
% now in RF2 the direction is x, so Ic2_xx.
% for the third omega3 rotates around z axis, so Ic3_zz
array_I_ci = [Ic1_yy Ic2_yy Ic3_zz];

[T, array_Ti] = compute_kinetic_energy_T_Ti("rpr", DH_params, [q1 q2 q3], [q1d q2d q3d], {r1_c1, r2_c2, r3_c3}, [m1, m2, m3], array_I_ci, [0,0,0]',  [0,0,0]', true);

%% 

pretty(T)

M = T_to_M(T, [q1d q2d q3d]);
M = simplify(expand(M),15);
disp(M)
%%  exercise 5
clear all
clc

syms a1 a2 a3 a4 a5 a6 real
syms q1 q2 q3 real
syms q1d q2d q3d real

M = sym(zeros(3,3));
M(1,1) = a1+ 2*a2*q2 + a3*q2^2 + 2*a4*q2*sin(q3) + a5*sin(q3)^2;
M(2,2) = a3;
M(2,3) = a4*cos(q3);
M(3,2) = a4*cos(q3);
M(3,3) = a6;

disp(M)

%% 
% first request is to compute the coriolis and centrifugal terms 

c= M_to_C(M,[q1,q2,q3],[q1d,q1d,q3d]');
disp(c)

%%
% second request
% we can compute easily the first factorization or using the coded function
% to get the c(q, q dot) or using another function for just the
% Coriolis-Centrifugal terms. let's do with both.

% first factorization S
q = [q1,q2,q3]';
qd = [q1d,q2d,q3d]';

 S_s = M_to_S_skew_sym(M, q, qd, true);

 %%
 % % first factorization S by hand

 C = christoffel_symbols(M, q);
 S_s_1= qd' * C{1};
 S_s_2= qd' * C{2};
 S_s_3= qd' * C{3};

 S_s2 = [S_s_1;S_s_2;S_s_3];

 disp(S_s2)
%% 
% test the first factorization S 

isSkewSymmetric = check_factorization(M, q, qd, S_s, true);


%% 
% let's find the second factorization 

 s_prime = S_to_newS(S_s,qd);
 disp(s_prime)


%% 

isSkewSymmetric = check_factorization(M, q, qd, s_prime, true);

%% 

% third factorization matrix
factor_3 = sym(zeros(3,3));
factor_3(3,2) = q3d;
factor_3(3,3) = -q2d;
disp(factor_3)
S_3 = S_s +  factor_3;
disp(S_3)

% the be a factorization we have to satisfy this:

condition = simplify(expand(S_3*qd),15);
c = simplify(expand(c),15);

%disp(isequal(expand(condition), expand(c)))
disp(condition)
disp(c)

isSkewSymmetric = check_factorization(M, q, qd, S_3, true);
%% 
% find the regressor matrix
disp(M) 
disp(c)

pi_a = [a1 a2 a3 a4 a5 a6]';
disp(pi_a);

% find the regressor matrix
syms q1dd q2dd q3dd real

qdd = [q1dd q2dd q3dd]';
u = simplify(M*qdd + c);
disp(u)


Y = u_to_Y(u,pi_a);

disp(Y)
%% exercise 6









