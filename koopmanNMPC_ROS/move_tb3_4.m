clear
clc
clear defaultNode
clear newDomainNode

%%
setenv('ROS_DOMAIN_ID','30');

% === CONFIGURAZIONE INIZIALE ===
node = ros2node('/matlab_controller_node', 30);

odomSub = ros2subscriber(node, "/tb3_4/odom", "nav_msgs/Odometry");
cmdVelPub = ros2publisher(node, "/tb3_4/cmd_vel", "geometry_msgs/TwistStamped", "Reliability", "besteffort"); 


ref_pos = [1, 1];        
ref_angle = deg2rad(-135);     

ref_reached = false;

ref_tol = 1e-1;
ref_angle_tol = 1e-1;
max_lin_vel = 0.15;
max_ang_vel = 1;

%%

velMsg = ros2message(cmdVelPub);

while ~ref_reached
    
    % Receive the odom msg
    odomMsg = receive(odomSub);

    % Extract the position
    position_rob = [odomMsg.pose.pose.position.x, odomMsg.pose.pose.position.y];
    
    % Extract the orientation and convert from quaternion to euler
    quat = [odomMsg.pose.pose.orientation.w, odomMsg.pose.pose.orientation.x, odomMsg.pose.pose.orientation.y, odomMsg.pose.pose.orientation.z];
    eul = quat2eul(quat);
    orientation_rob = eul(1);
    
    % Compute the reference distance
    ref_dist = sqrt(sum((ref_pos(1:2)-position_rob).^2));
    
    % Creation of the velocity message
    
    

    if ref_dist<=ref_tol && ~ref_reached
        fprintf("Robot 2: position reached\n")

        velMsg.twist.linear.x = 0;

        ref_angle_dist = ref_angle - orientation_rob;

        if abs(ref_angle_dist) < ref_angle_tol
            velMsg.twist.angular.z = 0;
            ref_reached = true;
            fprintf("Robot 2: orientation reached\n")
        else
            velMsg.twist.angular.z = ref_angle_dist;
        end
        send(cmdVelPub, velMsg)
        pause(0.001)
        
        
    else
        
        angle_diff = wrapToPi(atan2(ref_pos(2) - position_rob(2), ref_pos(1) - position_rob(1)) - orientation_rob)
        
        if abs(angle_diff) > ref_angle_tol
            if abs(angle_diff)>max_ang_vel
                velMsg.twist.angular.z = sign(angle_diff)*max_ang_vel;
            else
                velMsg.twist.angular.z = angle_diff;
            end
        else
            if abs(ref_dist)> max_lin_vel
                velMsg.twist.linear.x = max_lin_vel;
            else 
                velMsg.twist.linear.x = ref_dist;
            end
        end


        velMsg.header.stamp = ros2time(node, "now");
        send(cmdVelPub, velMsg);
        pause(0.001)
        
    end

end

%% Stop the robot
% 

velMsg.twist.linear.x = 0.0; 
velMsg.twist.linear.y = 0.0; 
velMsg.twist.linear.z = 0.0; 

velMsg.twist.angular.x = 0.0; 
velMsg.twist.angular.y = 0.0; 
velMsg.twist.angular.z = 0.0; 
velMsg.header.stamp = ros2time(node, "now");
velMsg.header.frame_id = 'tb3_4/odom';
send(cmdVelPub, velMsg);
