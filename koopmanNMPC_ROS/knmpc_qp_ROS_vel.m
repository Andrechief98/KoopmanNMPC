clc
clear variables
close all
yalmip('clear');

clear defaultNode
clear newDomainNode
%% General settings

nx = 3;
nu = 2;

% Radius
r_ego = 0.25;

Ts = 0.2;

% % Initial condition
% x0 = [0, 0, 0]';

% Reference
x_r = [2, 2, 0]';

n_ref = size(x_r,2);

ref_tol = [0.1, 0.1, inf]';


%% Koopman

% CT system equations
x = sym('x',[nx,1],'real');
u = sym('u',[nu,1],'real');

f_sym = [u(1)*cos(x(3));
	u(1)*sin(x(3));
	u(2)];

% Initial set of observables
phi_init = [x; x(1)^2; x(2)^2];

N_obs = 20;

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

Np = 10;
Nc = 5;

Tp = Ts;


Q = diag([50, 50, 1e-03]);
R = diag([1e-03, 1]);
Qd = diag([1e-01, 1e-01, 1e-01]);
Rd = diag([0.1, 10]);
P = 1*Q;
S = diag(1e06*ones(nx,1));
S_nl = 1e06;

max_lin_vel = 0.2;
max_ang_vel = 1.5;

u_lb = [-max_lin_vel, -max_ang_vel]';
u_ub = [max_lin_vel, max_ang_vel]';

% Bound su state variables
x_lb = [-0.25, -0.25, -inf]';
x_ub = [2.25, 2.25, inf]';

% Obstacles [box, robot_2, bin] 
obst.x_c = [0.75, 1.75, 1.5];
obst.y_c = [0.4, 0.75, 1.75];

obst.lx = [0.4/sqrt(2), 0.3, 0.15];
obst.ly = [0.6/sqrt(2), 0.3, 0.15];

% obst.v_x = [0, -0.1, 0];
% obst.v_y = [0, 0.1, 0];

obst.alpha = [1.1, 1.1, 1.1];

obst.type = {"r", "c", "c"};

n_obst = length(obst.x_c);
r_obst = 0.3;

%% MPC controller

% Yalmip variables

% Obstacles position and velocity[
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

u_p = sdpvar(nu,1);

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
cost = cost + quad_form(u{ind_u(1)}-u_p, Rd);

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
	constr = [constr; x_lb + [r_ego; r_ego; 0] - e <= z{i}(1:nx) <= x_ub - [r_ego; r_ego; 0] + e];
	for j=1:1:n_obst
		% constr = [constr; z{i}(4) - 2*(c_obst(1,j) + v_obst(1,j)*(i-1)*Ts)*z{i}(1) + (c_obst(1,j) + v_obst(1,j)*(i-1)*Ts)^2 + ...
		% 	z{i}(5) - 2*(c_obst(2,j) + v_obst(2,j)*(i-1)*Ts)*z{i}(2) + (c_obst(2,j) + v_obst(2,j)*(i-1)*Ts)^2 >= ...
		% 	(obst.r(j)+r_ego)^2 - e_nl];

		a_o = (obst.lx(j) + obst.alpha(j)*r_ego)^2;
		b_o = (obst.ly(j) + obst.alpha(j)*r_ego)^2;

		constr = [constr; b_o*( z{i}(4) - 2*(c_obst(1,j) + v_obst(1,j)*(i-1)*Ts)*z{i}(1) + (c_obst(1,j) + v_obst(1,j)*(i-1)*Ts)^2 ) + ...
			a_o*( z{i}(5) - 2*(c_obst(2,j) + v_obst(2,j)*(i-1)*Ts)*z{i}(2) + (c_obst(2,j) + v_obst(2,j)*(i-1)*Ts)^2 ) >= ...
			a_o*b_o - e_nl];
	end
end

% Input bounds
for i = 1:1:Nc
	constr = [constr; u_lb <= u{i} <= u_ub];
end

% Slack variables constr.
constr = [constr; e >= 0; e_nl >= 0];

% Optimizer object (MPC controller)

params_in = {z_in, x_ref, Ad, Bd, bd, c_obst, v_obst, u_p};

sol_out = {[u{:}], [z{:}]};

options = sdpsettings('verbose',0,'solver','mosek');

% MOSEK settings
options.mosek.MSK_DPAR_INTPNT_QO_TOL_DFEAS = 1e-06;
options.mosek.MSK_DPAR_INTPNT_QO_TOL_PFEAS = 1e-06;

knmpc = optimizer(constr,cost,options,params_in,sol_out);

%%
% Reference robot1
ref1 = x_r;
ref1_reached = false;
ref1_tol = 2e-2;

% Reference robot2
ref2 = [0.25; 2.25];
ref2_reached = false;
ref2_tol = 2e-1;
ref2_angle_tol = 1e-1;
ref2_angle = deg2rad(-135); 
max_lin_vel_rob2 = 0.13;


%% ROS2 Network
setenv('ROS_DOMAIN_ID','30');
ros2 topic list

% Creating a ROS2 node
node_robot1=ros2node('/node_robot1', 30);
node_robot2=ros2node("/node_robot2", 30);


% Definition of:
%   - subscribers: for odometry (both robots);
%   - publishers: for velocity commands (both robots);
odom_sub_1 = ros2subscriber(node_robot1, "/tb3_3/odom","nav_msgs/Odometry"); 
odom_sub_2 = ros2subscriber(node_robot2, "/tb3_4/odom","nav_msgs/Odometry"); 

vel_pub_1 = ros2publisher(node_robot1, "/tb3_3/cmd_vel", "geometry_msgs/TwistStamped", "Reliability", "besteffort"); 
vel_pub_2 = ros2publisher(node_robot1, "/tb3_4/cmd_vel", "geometry_msgs/TwistStamped", "Reliability", "besteffort"); 


%% Control loop

counter = 0;
u_o = zeros(nu,1);
u_p = zeros(nu,1);
velMsg_rob1 = ros2message(vel_pub_1);
velMsg_rob1.header.frame_id = 'tb3_3/odom';

velMsg_rob2 = ros2message(vel_pub_2);
velMsg_rob2.header.frame_id = 'tb3_4/odom';

while ~ref1_reached

    % Receiving odometry (both robots)
    odom_msg_rob1 = receive(odom_sub_1);
    odom_msg_rob2 = receive(odom_sub_2);

    % Extracting information for robot1 (K-NMPC)
    pose_rob1 = odom_msg_rob1.pose.pose;
    twist_rob1 = odom_msg_rob1.twist.twist;

    position_rob1 = [pose_rob1.position.x; pose_rob1.position.y];
    eul_angles_rob1 = quat2eul([pose_rob1.orientation.w, pose_rob1.orientation.x, pose_rob1.orientation.y, pose_rob1.orientation.z]); %from quaternion to Euler angles (ZXY)
    orientation_rob1 = eul_angles_rob1(1);
    lin_vel_rob1 = twist_rob1.linear.x;
    ang_vel_rob1 = twist_rob1.angular.z;
    
    % Extracting information for robot2
    pose_rob2 = odom_msg_rob2.pose.pose;
    twist_rob2 = odom_msg_rob2.twist.twist;

    position_rob2 = [pose_rob2.position.x; pose_rob2.position.y];
    eul_angles_rob2 = quat2eul([pose_rob2.orientation.w, pose_rob2.orientation.x, pose_rob2.orientation.y, pose_rob2.orientation.z]); %from quaternion to Euler angles (ZXY)
    orientation_rob2 = eul_angles_rob2(1);
    lin_vel_rob2 = twist_rob2.linear.x;
    ang_vel_rob2 = twist_rob2.angular.z;

    
    % Computing the current robot distance from corresponding reference
    ref1_dist = sqrt(sum((ref1(1:2)-position_rob1).^2));
    ref2_dist = sqrt(sum((ref2(1:2)-position_rob2).^2));


    % NMPC control for robot1
    if ref1_dist<=ref1_tol && ~ref1_reached
        % Reference of robot1 reached
        velMsg_rob1.linear.x = 0;
        velMsg_rob1.angular.z = 0;
        velMsg_rob1.header.stamp = ros2time(node_robot1,"now");
        send(vel_pub_1, velMsg_rob1);
        send(vel_pub_1, velMsg_rob1)
        fprintf("Robot 1 has reached its reference\n")
        ref1_reached = true;
    else
        % KNMPC
        x_rob1 = position_rob1(1);
        y_rob1 = position_rob1(2);
        v_lin = lin_vel_rob1;
        theta = orientation_rob1;

        x_mpc = [x_rob1; y_rob1; theta];
        z_mpc = phi(x_mpc);

	    % Pred. model lin.
	    A_lin = Jz_fun(z_mpc,u_o) + B{1}*u_o(1) + B{2}*u_o(2);
	    B_lin = Ju_fun(z_mpc,u_o) + [B{1}*z_mpc, B{2}*z_mpc];
	    b_lin = b_fun(z_mpc,u_o) - [B{1}*z_mpc, B{2}*z_mpc]*u_o;
    
	    Ad = (eye(nz) + Tp*(A + A_lin) + Tp^2/2*(A + A_lin)^2 + Tp^3/6*(A + A_lin)^3 + Tp^4/24*(A + A_lin)^4);
	    Bd = (Tp*eye(nz) + Tp^2/2*(A + A_lin) + Tp^3/6*(A + A_lin)^2 + Tp^4/24*(A + A_lin)^3)*(B0 + B_lin);
	    bd = (Tp*eye(nz) + Tp^2/2*(A + A_lin) + Tp^3/6*(A + A_lin)^2 + Tp^4/24*(A + A_lin)^3)*b_lin;
        
        % Absolute linear velocity of robot2 [Vx,Vy]
        absolute_lin_vel_obstacles = [
            0, lin_vel_rob2*cos(orientation_rob2), 0; 
            0, lin_vel_rob2*sin(orientation_rob2), 0
            ];
        position_obstacles = [
            [obst.x_c(1); obst.y_c(1)], position_rob2, [obst.x_c(3); obst.y_c(3)]
            ];
	   
        % Compute optimal input    
	    sol = knmpc({z_mpc, ref1, Ad, Bd, bd, position_obstacles, absolute_lin_vel_obstacles, u_p});
	    u_mpc = sol{1}(:,1);
        u_p = u_mpc;
        u_o = 0.5*u_mpc;
        
        % Extraction of linear velocity to apply
	    ctrl_inp_lin_vel = u_mpc(1)

        % Extraction of angular velocity to apply
        ctrl_inp_ang_vel = u_mpc(2)

	    % Velocity message
        velMsg_rob1.twist.linear.x = ctrl_inp_lin_vel;
        velMsg_rob1.twist.angular.z = ctrl_inp_ang_vel;
        velMsg_rob1.header.stamp = ros2time(node_robot1,"now");        
        send(vel_pub_1, velMsg_rob1);
        % pause(0.001)
    end 
    
    % Movement of robot2
    if ref2_dist<=ref2_tol && ~ref2_reached
        % Reference of robot2 reached
        fprintf("Robot 2: reference position reached\n")
        
        velMsg_rob2.twist.linear.x = 0;
        
        ref2_angle_dist = ref2_angle - orientation_rob2;

        if abs(ref2_angle_dist) < ref2_angle_tol
            velMsg_rob2.twist.angular.z = 0;
            ref2_reached = true;
            fprintf("Robot 2: orientation reached\n")
        else
            velMsg_rob2.twist.angular.z = ref2_angle_dist;
        end
        
        velMsg_rob2.header.stamp = ros2time(node_robot2, "now");
        send(vel_pub_2, velMsg_rob2)
        % pause(0.001)

    else
        % proportional linear movement
        angle_diff_2 = wrapToPi(atan2(ref2(2) - position_rob2(2), ref2(1) - position_rob2(1)) - orientation_rob2);

        if abs(angle_diff_2) > ref2_angle_tol
            if abs(angle_diff_2)>max_ang_vel
                velMsg_rob2.twist.angular.z = sign(angle_diff_2)*max_ang_vel;
            else
                velMsg_rob2.twist.angular.z = angle_diff_2;
            end
        else
            if abs(ref2_dist)> max_lin_vel_rob2
                velMsg_rob2.twist.linear.x = max_lin_vel_rob2;
            else 
                velMsg_rob2.twist.linear.x = ref2_dist;
            end
        end

        velMsg_rob2.twist.linear.y = 0.0; 
        velMsg_rob2.twist.linear.z = 0.0; 
        
        velMsg_rob2.twist.angular.x = 0.0; 
        velMsg_rob2.twist.angular.y = 0.0; 
        velMsg_rob2.header.stamp = ros2time(node_robot2, "now");
        send(vel_pub_2, velMsg_rob2);
        % pause(0.001)
        
    end

       
end

%%
velMsg_rob1.twist.linear.x = 0.0; 
velMsg_rob1.twist.linear.y = 0.0; 
velMsg_rob1.twist.linear.z = 0.0; 

velMsg_rob1.twist.angular.x = 0.0; 
velMsg_rob1.twist.angular.y = 0.0; 
velMsg_rob1.twist.angular.z = 0.0; 
velMsg_rob1.header.stamp = ros2time(node_robot1,"now");
velMsg_rob1.header.frame_id = 'tb3_3/odom';
send(vel_pub_1, velMsg_rob1);

%%
velMsg_rob2.twist.linear.x = 0.0; 
velMsg_rob2.twist.linear.y = 0.0; 
velMsg_rob2.twist.linear.z = 0.0; 

velMsg_rob2.twist.angular.x = 0.0; 
velMsg_rob2.twist.angular.y = 0.0; 
velMsg_rob2.twist.angular.z = 0.0; 
velMsg_rob2.header.stamp = ros2time(node_robot2,"now");
velMsg_rob2.header.frame_id = 'tb3_4/odom';
send(vel_pub_2, velMsg_rob2);