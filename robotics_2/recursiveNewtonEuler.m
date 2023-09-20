function RNE()
clc
clear all
close all

VISUALIZE = false;

% -------------------------------------------------------------------------
% Define the DH parameters
N_DOFS = 3;
dh.theta = [0 0 0];
dh.alpha = [0 0 0];
dh.offset = [0 0 0];        % always zero vector for us
dh.d = [0 0 0];
dh.a = [0.5 0.5 0.5];
dh.type = ['r' 'r' 'r'];

% -------------------------------------------------------------------------
% Rigid body paramaters: inertia, mass, and center of mass for each link
% (in this case 3 links)
rb.I =  repmat(eye(3), 1, 1, N_DOFS);
rb.m = [1 1 1];
% In standard DH, COM is mesured respect to its end-effector (using local
% frame). When COM is [0 0 0]', it means the COM is located exactly at the
% end-effector. Therefore, COM usually has negative values, which means the
% COM is behind the end-effector
rb.r = [-0.25 0 0; -0.25 0 0; -0.25 0 0]'; % this is the ri_ci in frame i vector

% -------------------------------------------------------------------------
% Arbitrary trajectory as the inputs: joint position, velocity, and 
% acceleration
ts = 0.001;
time_span = 0:ts:1;
qc = [pi/3*sin(2*pi*1*time_span)' pi/3*sin(2*pi*1*time_span)' pi/3*sin(2*pi*1*time_span)'];
qcdot = gradient(qc', ts)';
qcddot = gradient(qcdot', ts)';

% -------------------------------------------------------------------------
% Kinematic simulation, optional, for visualization purpose!   
% End-effector position, form the base which is located at [0 0 0]'

if VISUALIZE
    figure;
    hold on;
    h1 = plot3(0,0,0, 'b');   % The link
    h2 = plot3(0,0,0, 'or');  % The joint/end-effector
    h3 = plot3(0,0,0, 'xk');  % The center of mass
    xlabel('x')
    ylabel('y')
    zlabel('z')
    view([0 0 1]);
    axis equal;
    xlim([-2 2]);
    ylim([-2 2]);
    zlim([-0.5 0.5]);
    legend('Link', 'Joint', 'Center of Mass')
    
    EE = zeros(3, N_DOFS+1);
    COM = zeros(3, N_DOFS);    
    for k = 1 : length(time_span)  
        for i = 1 : 1 : N_DOFS
            T = calc_transformation(0, i, dh, qc(k,:));
            EE(:,i+1) = T(1:3,4);
            COM(:,i) = EE(:,i+1) + T(1:3,1:3) * rb.r(:,i);
        end
        
        % Draw the robot
        set(h1, 'XData', EE(1, :), 'YData', EE(2, :),'ZData', EE(3, :));
        set(h2, 'XData', EE(1, :), 'YData', EE(2, :),'ZData', EE(3, :));
        set(h3, 'XData', COM(1, :), 'YData', COM(2, :),'ZData', COM(3, :));
        drawnow;
    end
end 

% -------------------------------------------------------------------------
% get the torque for the trajectory
u_n = invdyn(dh, rb, qc, qcdot, qcddot, [0; -9.8; 0]);

disp(u_n)

% ------------------------------------------------------------------------
% Plotting
figure;
hold on;
plot(time_span, u_n(1,:), 'b');
plot(time_span, u_n(2,:), 'r');
plot(time_span, u_n(3,:), 'g');
end


% ---------------------------aux function for NER -------------------------
function Q = invdyn(dh, rb, qc, qcdot, qcddot, grav)    % we use the gravity trick
% Inverse dynamic with recursive Newton-Euler

%if nargin < 6
%    grav = [0;0;0];
%end

z0 = [0; 0; 1];

for k = 1 : length(qc)  
    q = qc(k, :);
    qdot = qcdot(k, :);
    qddot = qcddot(k, :);

    N_DOFS = length(q);
    
    % ---------------------------------------------------------------------
    % 1) Forward recursion (acceleration recursion)
    for i = 1 : N_DOFS
        T = calc_transformation(i-1, i, dh, q);
        R = T(1:3, 1:3);
        p = [dh.a(i); dh.d(i)*sin(dh.alpha(i)); dh.d(i)*cos(dh.alpha(i))]; % represents ri-1_i in reference frame i
        
        if i > 1
            % compute omega_i in frame i
            w(:, i) =  R'*(w(:, i-1) + z0.*qdot(i));
            % compute omegad_i in frame i
            wdot(:, i) = R'*(wdot(:, i-1) +  z0.*qddot(i) + ...
                cross(w(:, i-1), z0.*qdot(i)));
            % compute a_i in frame i
            vdot(:,i) = R'*vdot(:,i-1) + cross(wdot(:,i), p) + ...
                cross(w(:,i), cross(w(:,i),p));

            % a_ci in frame i is computed directly on the backward
            % recursion with the name vcdot
        
        else % if first joint we have omega_0 = 0, omegad_0 = 0 and a_0 = gravity vector
            w(:, i) =  R'*(z0.*qdot(i));
            wdot(:, i) = R'*(z0.*qddot(i));
            vdot(:,i) = R'*grav + cross(wdot(:,i), p) + ...
                cross(w(:,i), cross(w(:,i),p));
        end
    end
    
    % Dynamic simulation
    % 2) Backward recursion (force/torque computation)
    for i = N_DOFS:-1:1
        p = [dh.a(i); dh.d(i)*sin(dh.alpha(i)); dh.d(i)*cos(dh.alpha(i))];    % represents ri-1_i in reference frame i
        
        % this is the computation of a_ci in frame i, in the slides is
        % modelled in the forward recursion but it's the same
        vcdot = vdot(:,i) + cross(wdot(:,i),rb.r(:,i)) + ...
            cross(w(:,i),cross(w(:,i),rb.r(:,i)));
        
        % the update term 
        F = rb.m(i)*vcdot;                                              % force
        N = rb.I(:,:,i)*wdot(:,i) + cross(w(:,i),rb.I(:,:,i)*w(:,i));   % torque

        % to be added external forces initialization here (so f_n+1 and
        % n_n+2
        
        % summing contribution of i+1 generalized force
        % p -> ri-1_i in frame i
        % rb.r ->  ri_ci in frame i
        if i < N_DOFS
            T = calc_transformation(i, i+1, dh, q);
            R = T(1:3, 1:3);

            f(:,i) = R*f(:,i+1) + F;
            %n(:,i) = R*(n(:,i+1) + cross(R'*p, f(:,i+1))) + ...
            %    cross(rb.r(:,i)+p,F) + N;
            n(:,i) = R*(n(:,i+1)) + cross(R*f(:,i+1), rb.r(:,i)) - ...
                cross(f(:,i), p + rb.r(:,i)) + N;
            
        else
            n(:,i) = cross(rb.r(:,i)+p,F) + N;
            f(:,i) = F;
        end
        
        T = calc_transformation(i-1, i, dh, q);
        R = T(1:3, 1:3);
        
        
        % 3) Force projection
        % insert here dissipative term if needed like viscous friction
        % adding to Q(i,k) -> Fv(i) * qcdot
        
        if dh.type(i) == 't'
            Q(i,k) = f(:,i)'*R'*z0;
        elseif dh.type(i) == 'r'
            Q(i,k) = n(:,i)'*R'*z0;
        end
    end
end
end

% -------------------------------------------------------------------------
% get DH Transformation matrix
function  T = calc_transformation(from, to, dh, qc)

% qc -> position from the trajectory 
% Transformation from one joint to another joint
% 0<=from<N_DOFS
% 0<to<=N_DOFS

T = eye(4);
N_DOFS = length(qc);

% Sanity check
if (from >= N_DOFS) || (from < 0) || (to <= 0) || (to >  N_DOFS)
    return;
end

for i = from+1 : to
    % the sampled joint variables have a different interpretation in DH
    % based on the joint type: 
    if dh.type(i) == 'r'
        dh.theta(i) = qc(i);
    elseif dh.type(i) == 'p'
        dh.d(i) = qc(i);
    end
    
    ct = cos(dh.theta(i) + dh.offset(i));
    st = sin(dh.theta(i) + dh.offset(i));
    ca = cos(dh.alpha(i));
    sa = sin(dh.alpha(i));
    
    T = T * [ ct    -st*ca   st*sa     dh.a(i)*ct ; ...
              st    ct*ca    -ct*sa    dh.a(i)*st ; ...
              0     sa       ca        dh.d(i)    ; ...
              0     0        0         1          ];
end

end