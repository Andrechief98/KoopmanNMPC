clc
clear

%% DATA EXTRACTION
bagPath = "./tests/test_4"; 
bag = ros2bag(bagPath);

% We extract the relevant topics
topics = ["/tb3_3/odom", "/tb3_4/odom", "/tb3_3/cmd_vel", "/tb3_4/cmd_vel"];
bagSel = select(bag, "Topic", topics);

% We extract the messages from each topic
odom_tb3_3 = readMessages(select(bagSel, "Topic", "/tb3_3/odom"));
odom_tb3_4 = readMessages(select(bagSel, "Topic", "/tb3_4/odom"));
cmd_tb3_3 = readMessages(select(bagSel, "Topic", "/tb3_3/cmd_vel"));
cmd_tb3_4 = readMessages(select(bagSel, "Topic", "/tb3_4/cmd_vel"));

% Function to extract robot's position and velocity from odometry
extractOdom = @(msgs) struct( ...
    'x', cellfun(@(m) m.pose.pose.position.x, msgs), ...
    'y', cellfun(@(m) m.pose.pose.position.y, msgs), ...
    'vx', cellfun(@(m) m.twist.twist.linear.x, msgs), ...
    'omega_z', cellfun(@(m) m.twist.twist.angular.z, msgs), ...
    'sec', cellfun(@(m) m.header.stamp.sec, msgs), ...
    'nanosec', cellfun(@(m) m.header.stamp.nanosec, msgs), ...
    'theta', cellfun(@(m) ...
        quat2yaw([ ...
            m.pose.pose.orientation.w, ...
            m.pose.pose.orientation.x, ...
            m.pose.pose.orientation.y, ...
            m.pose.pose.orientation.z]), msgs) ...
);
% Function to extract robot's velocity from velocity topic
extractVel = @(msgs) struct( ...
    'vx', cellfun(@(m) m.twist.linear.x, msgs), ...
    'omega_z', cellfun(@(m) m.twist.angular.z, msgs), ...
    'sec', cellfun(@(m) m.header.stamp.sec, msgs), ...
    'nanosec', cellfun(@(m) m.header.stamp.nanosec, msgs) ...
);

% Velocity data extraction
cmdData.tb3_3 = extractVel(cmd_tb3_3);
cmdData.tb3_4 = extractVel(cmd_tb3_4);

cmdData.tb3_3.time = double(cmdData.tb3_3.sec) + double(cmdData.tb3_3.nanosec)*10^-9;
cmdData.tb3_4.time = double(cmdData.tb3_4.sec) + double(cmdData.tb3_4.nanosec)*10^-9;

cmdData.tb3_3.time = cmdData.tb3_3.time(:) - cmdData.tb3_3.time(1);
cmdData.tb3_4.time = cmdData.tb3_4.time(:) - cmdData.tb3_4.time(1);

% Odometry data extraction
odomData.tb3_3 = extractOdom(odom_tb3_3);
odomData.tb3_4 = extractOdom(odom_tb3_4);

odomData.tb3_3.time = double(odomData.tb3_3.sec) + double(odomData.tb3_3.nanosec)*10^-9;
odomData.tb3_4.time = double(odomData.tb3_4.sec) + double(odomData.tb3_4.nanosec)*10^-9;


% Rearrange the x-axis of the odometry data from robot_1 (K-NMPC
% controlled)
odomData.tb3_3.time = odomData.tb3_3.time(:) - odomData.tb3_3.time(1);

start_idx_tb3_3_time = find(odomData.tb3_3.x < 0.01);
odomData.tb3_3.time(start_idx_tb3_3_time) = [];
odomData.tb3_3.x(start_idx_tb3_3_time) = [];
odomData.tb3_3.y(start_idx_tb3_3_time) = [];
odomData.tb3_3.theta(start_idx_tb3_3_time) = [];
odomData.tb3_3.time = odomData.tb3_3.time - odomData.tb3_3.time(1);

% end_idx_tb3_3_time = find(odomData.tb3_3.x >=2);
% odomData.tb3_3.time(end_idx_tb3_3_time) = [];
% odomData.tb3_3.x(end_idx_tb3_3_time) = [];
% odomData.tb3_3.y(end_idx_tb3_3_time) = [];
% odomData.tb3_3.theta(end_idx_tb3_3_time) = [];
% 
% % Refining the final data:
% odomData.tb3_3.time = [odomData.tb3_3.time; 20];
% odomData.tb3_3.x = [odomData.tb3_3.x; odomData.tb3_3.x(end)];
% odomData.tb3_3.y = [odomData.tb3_3.y; odomData.tb3_3.y(end)];
% odomData.tb3_3.theta = [odomData.tb3_3.theta ; odomData.tb3_3.theta(end)];
% 
% cmdData.tb3_3.time = [cmdData.tb3_3.time; 20];
% cmdData.tb3_3.vx = [cmdData.tb3_3.vx ; cmdData.tb3_3.vx(end)];
% cmdData.tb3_3.omega_z = [cmdData.tb3_3.omega_z ; cmdData.tb3_3.omega_z(end)];

save("extracted_rosbags_data.mat")
%% PLOTS
% TRAJECTORIES PLOTS
clc
clear
close all
load("extracted_rosbags_data.mat")

% Radius of the two robots
r_rob1 = 0.25;
r_rob2 = 0.3;

% Obstacles [box, bin] 
obst.x_c = [0.75, 1.5];
obst.y_c = [0.4, 1.75];
obst.lx = [0.4/sqrt(2), 0.15];
obst.ly = [0.6/sqrt(2), 0.15];


% Bound su state variables
x_lb = [-0.25, -0.25, -inf]';
x_ub = [2.25, 2.25, inf]';

f1 = figure(1); set(f1,'WindowStyle','normal','color','w'); f1.Position = [200   200   400   400];

axis equal, grid on, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([x_lb(1), x_ub(1)]), ylim([x_lb(2), x_ub(2)])
xlabel('$x$ [m]','interpreter','latex')
ylabel('$y$ [m]','interpreter','latex')
title('Planar motion','interpreter','latex')
hold on

% Bounds
plot([x_lb(1), x_ub(1), x_ub(1), x_lb(1), x_lb(1)], [x_lb(2), x_lb(2), x_ub(2), x_ub(2), x_lb(2)], '-', 'linewidth', 1.25, 'color', [0.35, 0.35, 0.35])

% Box
plot_rect([obst.x_c(1); obst.y_c(1)], [obst.lx(1)*sqrt(2); obst.ly(1)*sqrt(2)], '-', 1.25, [0.35, 0.35, 0.35])
plot_obst([obst.x_c(1); obst.y_c(1)], [obst.lx(1); obst.ly(1)], '--', 1.25, [0.35 0.35 0.35])

% Bin
plot_obst([obst.x_c(2); obst.y_c(2)], [obst.lx(2); obst.ly(2)], '--', 1.25, [0.35 0.35 0.35])

step = 5;
% Robot 2 (tb3_4)
plot(odomData.tb3_4.x(1:step:end), odomData.tb3_4.y(1:step:end),  '-', 'linewidth', 1.25, 'color', [0.35 0.35 0.35])
plot(odomData.tb3_4.x(1:step:end), odomData.tb3_4.y(1:step:end), '.', 'markersize', 15, 'color', [0.35 0.35 0.35])

for i=1:step:length(odomData.tb3_4.x)
    plot_obst([odomData.tb3_4.x(i), odomData.tb3_4.y(i)], [r_rob2; r_rob2], '-', 1.5, [0.75 0.75 0.75])
end

% Robot 1 (tb3_3)
plot(odomData.tb3_3.x(1:step:end), odomData.tb3_3.y(1:step:end), '-', 'linewidth', 1.5, 'color', [0.5 0.5 1])
plot(odomData.tb3_3.x(1:step:end), odomData.tb3_3.y(1:step:end), 'b.', 'markersize', 15)

for i=1:step:length(odomData.tb3_3.x)
    plot_obst([odomData.tb3_3.x(i), odomData.tb3_3.y(i)], [r_rob1; r_rob1], '-', 1.5, [0.5 0.5 1])
end



% Reference
ref = [2,2];
plot(ref(1), ref(2) , '.', 'markersize', 20, 'color', [0 0.75 0])

hold off



% ODOMETRY AND CONTROL INPUT PLOTS
f2 = figure(2); set(f2,'WindowStyle','normal','color','w'); f2.Position = [600   200   700   400];

tiledlayout(2,3,'tilespacing','tight','padding','tight')

nexttile(1), hold on


plot(odomData.tb3_3.time, odomData.tb3_3.x, 'b-')

hold off
grid on, grid minor, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([0, odomData.tb3_3.time(end)])
xlabel('$t$ [s]','interpreter','latex')
title('$x$ [m]','interpreter','latex')

nexttile(2), hold on

plot(odomData.tb3_3.time, odomData.tb3_3.y, 'b-')

hold off
grid on, grid minor, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([0, odomData.tb3_3.time(end)])
xlabel('$t$ [s]','interpreter','latex')
title('$y$ [m]','interpreter','latex')

nexttile(3), hold on

plot(odomData.tb3_3.time, rad2deg(odomData.tb3_3.theta), 'b-')

hold off
grid on, grid minor, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([0, odomData.tb3_3.time(end)])
xlabel('$t$ [s]','interpreter','latex')
title('Heading $\theta$ [deg]','interpreter','latex')

nexttile(4), hold on

plot(cmdData.tb3_3.time, cmdData.tb3_3.vx, 'r.-')

hold off
grid on, grid minor, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([0, cmdData.tb3_3.time(end)])
xlabel('$t$ [s]','interpreter','latex')
title('Velocity $v$ [m/s]','interpreter','latex')

nexttile(5), hold on

plot(cmdData.tb3_3.time, rad2deg(cmdData.tb3_3.omega_z), 'r.-')

hold off
grid on, grid minor, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([0, cmdData.tb3_3.time(end)])
xlabel('$t$ [s]','interpreter','latex')
title('Angular velocity $\omega$ [deg/s]','interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_obst(c, size, line_style, line_width, color)
fplot(@(t) c(1) + size(1)*sin(t), @(t) c(2) + size(2)*cos(t), line_style, 'linewidth', line_width, 'color', color);
end

function yaw = quat2yaw(q)
% q = [w, x, y, z] (Matlab quaterion format)
w = q(1); x = q(2); y = q(3); z = q(4);

% Conversion from quaternion to euler angles (ZXY)
eul_angles = quat2eul([w, x, y, z]);
yaw = eul_angles(1);
end

function plot_rect(c, size, line_style, line_width, color)
plot([c(1)-size(1)/2; c(1)+size(1)/2; c(1)+size(1)/2; c(1)-size(1)/2; c(1)-size(1)/2], ...
	[c(2)-size(2)/2; c(2)-size(2)/2; c(2)+size(2)/2; c(2)+size(2)/2; c(2)-size(2)/2], line_style, 'color', color, 'linewidth', line_width)
end






