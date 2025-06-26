clc
clear variables
close all
yalmip('clear');

clear defaultNode
clear newDomainNode
%% General settings

nx = 4;
nu = 2;

% Radius
r_ego = 0.5;

Ts = 0.05;

% Initial condition
x0 = [0, 0, 0, 0]';

% Reference
x_r = [0, 0, 0, 0]';

n_ref = size(x_r,2);

ref_tol = [0.1, 0.1, inf, inf]';

u_o = [0;0];

%% Koopman

% CT system equations
x = sym('x',[nx,1],'real');
u = sym('u',[nu,1],'real');

f_sym = [x(3)*cos(x(4));
	x(3)*sin(x(4));
	u(1);
	u(2)];

% Initial set of observables
phi_init = [x; x(1)^2; x(2)^2];

N_obs = 16;

% Compute reduced Koopman lifted system
[Ap,B0p,Bp,As,B0s,Bs,~,phi_sym,Jz_sym,Ju_sym,b_sym] = koopman_calc_red_fun(x,u,f_sym,phi_init,N_obs);

fprintf('-------------------------\n')

nz = length(phi_sym);

z = sym('z',[nz,1],'real');

phi = matlabFunction(phi_sym,'vars',{[x]});

Jz_fun = matlabFunction(Jz_sym,'vars',{[z],[u]});
Ju_fun = matlabFunction(Ju_sym,'vars',{[z],[u]});
b_fun = matlabFunction(b_sym,'vars',{[z],[u]});

A = [Ap; As];
B0 = [B0p; B0s];
B{1} = [Bp{1}; Bs{1}];
B{2} = [Bp{2}; Bs{2}];


%% MPC data

Np = 15;
Nc = 5;

Tp = Ts;


Q = diag([30, 30, 1, 1e-03]);
R = diag([1e-03, 1e-03]);
Qd = diag([1e-01, 1e-01, 10, 10]);
Rd = diag([1, 1]);
P = 1*Q;
S = diag(1e06*ones(nx,1));
S_nl = 1e06;

u_lb = [-10, -pi]';
u_ub = [10, pi]';

x_lb = [-0.5, -0.5, -10, -inf]';
x_ub = [3.5, 3.5, 10, inf]';

n_obst = 1;
r_obst = 0.5;

%% MPC controller

% Yalmip variables

% Obstacles position and velocity
c_obst = sdpvar(2,n_obst,'full');
v_obst = sdpvar(2,n_obst,'full');

% Lifted states
z = sdpvar(nz*ones(1,Np+1),1*ones(1,Np+1));
z_in = sdpvar(nz,1);

% Inputs
u = sdpvar(nu*ones(1,Nc),1*ones(1,Nc));
if Nc == 1
	u = {u};
end

% Reference
x_ref = sdpvar(nx,1);

% Additional terms for reduced lifting
Ad = sdpvar(nz,nz,'full');
Bd = sdpvar(nz,nu,'full');
bd = sdpvar(nz,1);

% Slack variables
e = sdpvar(nx,1);
e_nl = sdpvar(1,1);

% Input parametrization

ind_u = zeros(Np,1);


j = 0;
for i = 1:1:Np
	if mod(i-1, floor(Np/Nc)) == 0
		if j < Nc
			j = j+1;
		end
	end
	ind_u(i) = j;
end


% Cost function

cost = 0;

% Stage cost
for i = 1:1:Np
	cost = cost + quad_form(z{i}(1:nx)-x_ref, Q);
	cost = cost + quad_form(u{ind_u(i)}, R);
end

for i = 2:1:Np
	cost = cost + quad_form(z{i}(1:nx)-z{i-1}(1:nx), Qd);
	cost = cost + quad_form(u{ind_u(i)}-u{ind_u(i-1)}, Rd);
end

cost = cost + quad_form(z{i}(1:nx)-x_ref, P);

% Slack variables penalty
cost = cost + quad_form(e, S);
cost = cost + quad_form(e_nl, S_nl);

% Constraints

constr = [];

% Initial condition
constr = [constr; z{1} == z_in];

% Prediction model
for i = 1:1:Np
	constr = [constr; z{i+1} == Ad*z{i} + Bd*u{ind_u(i)} + bd];
end

% State constraints
for i = 1:1:Np
	constr = [constr; x_lb + [r_ego; r_ego; 0; 0] - e <= z{i}(1:nx) <= x_ub - [r_ego; r_ego; 0; 0] + e];
	for j=1:1:n_obst
		constr = [constr; z{i}(5) - 2*(c_obst(1,j) + v_obst(1,j)*(i-1)*Ts)*z{i}(1) + (c_obst(1,j) + v_obst(1,j)*(i-1)*Ts)^2 + ...
			z{i}(6) - 2*(c_obst(2,j) + v_obst(2,j)*(i-1)*Ts)*z{i}(2) + (c_obst(2,j) + v_obst(2,j)*(i-1)*Ts)^2 >= ...
			(r_obst+r_ego)^2 - e_nl];
	end
end

% Input bounds
for i = 1:1:Nc
	constr = [constr; u_lb <= u{i} <= u_ub];
end

% Slack variables constr.
constr = [constr; e >= 0; e_nl >= 0];

% Optimizer object (MPC controller)

params_in = {z_in, x_ref, Ad, Bd, bd, c_obst, v_obst};

sol_out = {[u{:}], [z{:}]};

options = sdpsettings('verbose',0,'solver','mosek');

% MOSEK settings
options.mosek.MSK_DPAR_INTPNT_QO_TOL_DFEAS = 1e-06;
options.mosek.MSK_DPAR_INTPNT_QO_TOL_PFEAS = 1e-06;

knmpc = optimizer(constr,cost,options,params_in,sol_out);

%%
% Reference robot1
% ref1 = [3;3;0;0];
ref1 = x_r;
ref1_reached = false;
ref1_tol = 3e-1;

% Reference robot2
ref2 = [0;0];
ref2_reached = false;
ref2_tol = 3e-1;


%% ROS2 Network
setenv('ROS_DOMAIN_ID','30');
ros2 topic list

% Creating a ROS2 node
node_robot1=ros2node('/node_robot1', 30);
node_robot2=ros2node("/node_robot2", 30);


% Definition of:
%   - subscribers: for odometry (both robots);
%   - publishers: for velocity commands (both robots);
odom_sub_1 = ros2subscriber(node_robot1, "/odom","nav_msgs/Odometry"); 
odom_sub_2 = ros2subscriber(node_robot2, "/turtlebot2/odom","nav_msgs/Odometry"); 

vel_pub_1 = ros2publisher(node_robot1, "/cmd_vel", "geometry_msgs/TwistStamped"); 
vel_pub_2 = ros2publisher(node_robot1, "/turtlebot2/cmd_vel", "geometry_msgs/TwistStamped"); 


%% Control loop
% f1 = figure(1);
% 
% axis equal, grid on, box on
% set(gca,'TickLabelInterpreter','latex','fontsize',12)
% xlim([-5,7])
% xlabel('$x$ [m]','interpreter','latex')
% ylabel('$y$ [m]','interpreter','latex')
% title('Planar motion','interpreter','latex')
counter = 0;
while ~ref1_reached

    % hold on
    % Receiving odometry (both robots)
    odom_msg_rob1 = receive(odom_sub_1);
    % odom_msg_rob2 = receive(odom_sub_2);

    % Extracting information for robot1 (K-NMPC)
    pose_rob1 = odom_msg_rob1.pose.pose;
    twist_rob1 = odom_msg_rob1.twist.twist;

    position_rob1 = [pose_rob1.position.x; pose_rob1.position.y];
    eul_angles_rob1 = quat2eul([pose_rob1.orientation.w, pose_rob1.orientation.x, pose_rob1.orientation.y, pose_rob1.orientation.z]); %from quaternion to Euler angles (ZXY)
    orientation_rob1 = eul_angles_rob1(1);
    lin_vel_rob1 = twist_rob1.linear.x;
    ang_vel_rob1 = twist_rob1.angular.z;


    
    % Extracting information for robot2
    % pose_rob2 = odom_msg_rob2.pose.pose;
    % twist_rob2 = odom_msg_rob2.twist.twist;
    % 
    % position_rob2 = [pose_rob2.position.x; pose_rob2.position.y];
    % eul_angles_rob2 = quat2eul([pose_rob2.orientation.w, pose_rob2.orientation.x, pose_rob2.orientation.y, pose_rob2.orientation.z]); %from quaternion to Euler angles (ZXY)
    % orientation_rob2 = eul_angles_rob2(1);
    % lin_vel_rob2 = twist_rob2.linear.x;
    % ang_vel_rob2 = twist_rob2.angular.z;

    position_rob2 = [6; 6];
    orientation_rob2 = 0;
    lin_vel_rob2 = 0;
    ang_vel_rob2 = 0;
    % fprintf("Ciao")
    % Computing the current robot distance from corresponding reference
    ref1_dist = sqrt(sum((ref1(1:2)-position_rob1).^2));
    % ref2_dist = sqrt(sum((ref2(1:2)-position_rob2).^2));
    
    % NMPC control for robot1
    if ref1_dist<=ref1_tol && ~ref1_reached
        % Reference of robot1 reached
        velMsg_rob1 = ros2message(vel_pub_1);
        velMsg_rob1.linear.x = 0;
        velMsg_rob1.angular.z = 0;
        send(vel_pub_1, velMsg_rob1)
        fprintf("Robot 1 has reached its reference\n")
        ref1_reached = true;
    else
        % KNMPC
        x_rob1 = position_rob1(1);
        y_rob1 = position_rob1(2);
        v_lin = lin_vel_rob1;
        theta = orientation_rob1;

        x_mpc = [x_rob1; y_rob1; v_lin; theta]
        z_mpc = phi(x_mpc);

	    % Pred. model lin.
	    A_lin = Jz_fun(z_mpc,u_o) + B{1}*u_o(1) + B{2}*u_o(2);
	    B_lin = Ju_fun(z_mpc,u_o) + [B{1}*z_mpc, B{2}*z_mpc];
	    b_lin = b_fun(z_mpc,u_o) - [B{1}*z_mpc, B{2}*z_mpc]*u_o;
    
	    Ad = (eye(nz) + Tp*(A + A_lin) + Tp^2/2*(A + A_lin)^2 + Tp^3/6*(A + A_lin)^3 + Tp^4/24*(A + A_lin)^4);
	    Bd = (Tp*eye(nz) + Tp^2/2*(A + A_lin) + Tp^3/6*(A + A_lin)^2 + Tp^4/24*(A + A_lin)^3)*(B0 + B_lin);
	    bd = (Tp*eye(nz) + Tp^2/2*(A + A_lin) + Tp^3/6*(A + A_lin)^2 + Tp^4/24*(A + A_lin)^3)*b_lin;
        
        % Absolute linear velocity of robot2 [Vx,Vy]
        absolute_lin_vel_rob2 = [lin_vel_rob2*cos(orientation_rob2); lin_vel_rob2*sin(orientation_rob2)];
	   
        % Compute optimal input    
	    sol = knmpc({z_mpc, ref1, Ad, Bd, bd, position_rob2, absolute_lin_vel_rob2});
	    u_mpc = sol{1}(:,1);
        
        % Extraction of angular velocity to apply
        ctrl_inp_ang_vel = u_mpc(2)

	    % Kinematic model
	    ctrl_inp_lin_vel_next = x_mpc(3) + Ts*u_mpc(1)

        if abs(ctrl_inp_lin_vel_next)>0.5
            ctrl_inp_lin_vel_next = sign(ctrl_inp_lin_vel_next)*0.5; % Limit linear velocity to max 0.5
        end

        if abs(ctrl_inp_ang_vel)>1.5
            ctrl_inp_ang_vel = sign(ctrl_inp_ang_vel)*1.5;
        end
    
	    % Velocity message
        velMsg_rob1 = ros2message(vel_pub_1);
        velMsg_rob1.twist.linear.x = ctrl_inp_lin_vel_next;
        velMsg_rob1.twist.angular.z = ctrl_inp_ang_vel;
        velMsg_rob1.header.stamp = ros2time(node_robot1,"now");
        velMsg_rob1.header.frame_id = 'odom';
        
        send(vel_pub_1, velMsg_rob1);
        pause(Ts)
    end
    
    
    % Movement of robot2
    % if ref2_dist<=ref2_tol && ~ref2_reached
    %     % Reference of robot2 reached
    %     velMsg_rob2.linear.x = 0;
    %     velMsg_rob2.angular.z = 0;
    %     send(vel_pub_2, velMsg_rob2)
    %     fprintf("Robot 2 has reached its reference\n")
    %     ref2_reached = true;
    % else
    %     % linear movement
    %     velMsg_rob2.linear.x = 0;
    %     velMsg_rob2.angular.z = 0;
    % end

    % plot(position_rob1(1), position_rob1(2), '.', 'LineWidth', 2,'color', 'b')
    % plot(position_rob2(1), position_rob2(2), '.', 'LineWidth', 2, 'color', 'r')
    % 
    % 
    % hold off
	% drawnow
	% pause(0.1)
	% if ~ref1_reached 
	% 	cla
    % end
    

end