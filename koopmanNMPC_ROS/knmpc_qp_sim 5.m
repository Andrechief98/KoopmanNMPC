clc
clear variables
close all
yalmip('clear');

%% General settings

% MPC
u_param_type = "piecewise";
pred_model_discr = "rk4";

% Simulation
plant_integr = "rk4";
disp_curr_exec_time = true;

% Plots
plot_pred = true;
N_pred_plot = 1;

% Data
save_data = false;
save_path = 'sim_data';
save_name = 'knmpc_qp_sim.mat';

%% Plant data

nx = 3;
nu = 2;

f = @(x,u) [u(1)*cos(x(3));
	u(1)*sin(x(3));
	u(2)];

% Radius
r_ego = 0.25;

%% Simulation data

T_sim = 20;
Ts = 0.2;

N = floor(T_sim/Ts);
t = 0:Ts:(N-1)*Ts;

% Initial condition
x0 = [0, 0, 0]';

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

u_lb = [-0.2, -1.5]';
u_ub = [0.2, 1.5]';

x_lb = [-0.5, -0.5, -inf]';
x_ub = [2.5, 2.5, inf]';

% Obstacles
obst.x_c = [0.75, 1.75, 1.5];
obst.y_c = [0.4, 0.75, 1.75];

obst.lx = [0.4/sqrt(2), 0.3, 0.15];
obst.ly = [0.6/sqrt(2), 0.3, 0.15];

obst.v_x = [0, -0.1, 0];
obst.v_y = [0, 0.1, 0];

obst.alpha = [1.1, 1, 1];

obst.type = {"r", "c", "c"};

n_obst = length(obst.x_c);

%% Plot scenario

f1 = figure(1); set(f1,'WindowStyle','normal','color','w'); f1.Position = [200   200   400   400];

hold on

% Obstacles
for i=1:1:n_obst
	plot_obst([obst.x_c(i); obst.y_c(i)], [obst.lx(i); obst.ly(i)], '--', 1.25, [0.5, 0.5, 0.5])

	if obst.type{i} == "r"
		plot_rect([obst.x_c(i); obst.y_c(i)], [obst.lx(i)*sqrt(2); obst.ly(i)*sqrt(2)], '-', 1.25, [0.35, 0.35, 0.35])
	end
end

% Bounds
plot([x_lb(1), x_ub(1), x_ub(1), x_lb(1), x_lb(1)], [x_lb(2), x_lb(2), x_ub(2), x_ub(2), x_lb(2)], ...
	'-', 'linewidth', 1.25, 'color', [0.35, 0.35, 0.35])

% Reference
plot(x_r(1,:), x_r(2,:), '.', 'markersize', 20, 'color', [0 0.75 0])

% Ego robot
plot_obst([x0(1); x0(2)], [r_ego; r_ego], '-', 1.25, [0 0 1])

hold off

axis equal, grid on, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([x_lb(1), x_ub(1)]), ylim([x_lb(2), x_ub(2)])
xlabel('$x$ [m]','interpreter','latex')
ylabel('$y$ [m]','interpreter','latex')
title('Planar motion','interpreter','latex')

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

if u_param_type == "piecewise"
	% Piecewise constant over Np
	j = 0;
	for i = 1:1:Np
		if mod(i-1, floor(Np/Nc)) == 0
			if j < Nc
				j = j+1;
			end
		end
		ind_u(i) = j;
	end
end

if u_param_type == "first"
	% First Nc steps, then constant
	for i = 1:1:Np
		if i <= Nc
			ind_u(i) = i;
		else
			ind_u(i) = Nc;
		end
	end
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

%% Simulation

x_mpc = zeros(nx,N+1);
u_mpc = zeros(nu,N);
z_mpc = zeros(nz,N+1);

x_pred_mpc = cell(N,1);

z0 = phi(x0);

x_mpc(:,1) = x0;
z_mpc(:,1) = z0;

z_o = z0;
u_o = zeros(nu,1);

exec_time = [];

Ad = zeros(nz,nz);
Bd = zeros(nz,nu);
bd = zeros(nz,1);

i_ref = 1;

c_obst_curr = [obst.x_c; obst.y_c];
v_obst_curr = [obst.v_x; obst.v_y];

u_p = zeros(nu,1);

fprintf('K-NMPC (QP) simulation...\n');
for k=1:1:N

	z_o = z_mpc(:,k);

	% Pred. model lin.
	A_lin = Jz_fun(z_o,u_o) + B{1}*u_o(1) + B{2}*u_o(2);
	B_lin = Ju_fun(z_o,u_o) + [B{1}*z_o, B{2}*z_o];
	b_lin = b_fun(z_o,u_o) - [B{1}*z_o, B{2}*z_o]*u_o;

	% Pred. model discr.
	if pred_model_discr == "euler"
		% Euler
		Ad = eye(nz) + Tp*(A + A_lin);
		Bd = Tp*(B0 + B_lin);
		bd = Tp*b_lin;
	elseif pred_model_discr == "rk4"
		% RK4
		Ad = (eye(nz) + Tp*(A + A_lin) + Tp^2/2*(A + A_lin)^2 + Tp^3/6*(A + A_lin)^3 + Tp^4/24*(A + A_lin)^4);
		Bd = (Tp*eye(nz) + Tp^2/2*(A + A_lin) + Tp^3/6*(A + A_lin)^2 + Tp^4/24*(A + A_lin)^3)*(B0 + B_lin);
		bd = (Tp*eye(nz) + Tp^2/2*(A + A_lin) + Tp^3/6*(A + A_lin)^2 + Tp^4/24*(A + A_lin)^3)*b_lin;
	end

	% Compute optimal input
	tim1 = tic;

	sol = knmpc({z_mpc(:,k), x_r(:,i_ref), Ad, Bd, bd, c_obst_curr, v_obst_curr, u_p});
	
	curr_exec_time = toc(tim1);
	exec_time = [exec_time, curr_exec_time];

	if disp_curr_exec_time == true
		fprintf('K-NMPC (QP) step solved (Sim. time: t = %.3f s / %.3f s) (Exec. time: %.4f ms)\n', t(k), t(end), curr_exec_time*1e03)
	end

	u_mpc(:,k) = sol{1}(:,1);

	u_p = u_mpc(:,k);
	u_o = 0.5*u_mpc(:,k);

	x_pred_mpc{k} = sol{2}(1:nx,:);

	% Plant model integration
	if plant_integr == "euler"
		x_mpc(:,k+1) = x_mpc(:,k) + Ts*f(x_mpc(:,k), u_mpc(:,k));
	end
	if plant_integr == "rk4"
		k1 = f(x_mpc(:,k), u_mpc(:,k));
		k2 = f(x_mpc(:,k) + k1*Ts/2, u_mpc(:,k));
		k3 = f(x_mpc(:,k) + k2*Ts/2, u_mpc(:,k));
		k4 = f(x_mpc(:,k) + k3*Ts, u_mpc(:,k));
		x_mpc(:,k+1) = x_mpc(:,k) + Ts/6 * ( k1 + 2*k2 + 2*k3 + k4 );
	end

	z_mpc(:,k+1) = phi(x_mpc(:,k+1));

	% Move obstacles
	c_obst_curr = c_obst_curr + Ts*v_obst_curr;

	if all(abs(x_mpc(:,k)-x_r(:,i_ref)) < ref_tol(:,i_ref)) & i_ref < n_ref
		i_ref = i_ref + 1;
	elseif all(abs(x_mpc(:,k)-x_r(:,i_ref)) < ref_tol(:,i_ref)) & i_ref == n_ref
		N = k; T_sim = Ts*N;
		t(k+1:end) = [];
		x_mpc(:,k+2:end) = [];
		u_mpc(:,k+1:end) = [];
		break
	end
	
end

fprintf('-------------------------\n')
fprintf('Exec. time: [%.4f, %.4f] ms; mean = %.4f ms\n\n', min(exec_time)*1e03, max(exec_time)*1e03, mean(exec_time)*1e03)

%% Save data

if save_data == true
	save([save_path,'/',save_name], ...
		't', 'x_mpc', 'u_mpc', 'x_r', 'x_lb', 'x_ub', 'u_lb', 'u_ub', 'obst', 'exec_time')
end

%% Animated plot

f1 = figure(1); set(f1,'WindowStyle','normal','color','w'); f1.Position = [200   200   400   400];

axis equal, grid on, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([x_lb(1), x_ub(1)]), ylim([x_lb(2), x_ub(2)])
xlabel('$x$ [m]','interpreter','latex')
ylabel('$y$ [m]','interpreter','latex')
title('Planar motion','interpreter','latex')

for k=1:1:N
	hold on

	% Obstacles
	for i=1:1:n_obst
		plot_obst([obst.x_c(i) + obst.v_x(i)*(k-1)*Ts; obst.y_c(i) + obst.v_y(i)*(k-1)*Ts], ...
			[obst.lx(i); obst.ly(i)], '--', 1.25, [0.35 0.35 0.35])
	end

	% Bounds
	plot([x_lb(1), x_ub(1), x_ub(1), x_lb(1), x_lb(1)], [x_lb(2), x_lb(2), x_ub(2), x_ub(2), x_lb(2)], ...
		'-', 'linewidth', 1.25, 'color', [0.35, 0.35, 0.35])

	% Reference
	plot(x_r(1,:), x_r(2,:), '.', 'markersize', 20, 'color', [0 0.75 0])

	% Ego robot
	plot_obst([x_mpc(1,k); x_mpc(2,k)], [r_ego; r_ego], '-', 1.25, [0 0 1])

	% Trajectory
	plot(x_mpc(1,1:1:k), x_mpc(2,1:1:k), '-', 'linewidth', 1.5, 'color', [0.5 0.5 1])
	plot(x_mpc(1,k), x_mpc(2,k), 'b.', 'markersize', 15)

	hold off
	drawnow
	pause(0.1)
	if k < N
		cla
	end
end

%% Plots

f2 = figure(2); set(f2,'WindowStyle','normal','color','w'); f2.Position = [400   300   400   400];

hold on

% Obstacles
for i=1:1:n_obst
	for k=1:1:N-1
		plot_obst([obst.x_c(i) + obst.v_x(i)*(k-1)*Ts; obst.y_c(i) + obst.v_y(i)*(k-1)*Ts], ...
			[obst.lx(i); obst.ly(i)], '--', 1, [0.75 0.75 0.75])
	end
	plot_obst([obst.x_c(i) + obst.v_x(i)*(N-1)*Ts; obst.y_c(i) + obst.v_y(i)*(N-1)*Ts], ...
		[obst.lx(i); obst.ly(i)], '--', 1.25, [0.35 0.35 0.35])
end

% Bounds
plot([x_lb(1), x_ub(1), x_ub(1), x_lb(1), x_lb(1)], [x_lb(2), x_lb(2), x_ub(2), x_ub(2), x_lb(2)], ...
	'-', 'linewidth', 1.25, 'color', [0.35, 0.35, 0.35])

% Reference
plot(x_r(1,:), x_r(2,:), '.', 'markersize', 20, 'color', [0 0.75 0])

% MPC predictions
if plot_pred == true
	for k=1:N_pred_plot:N
		plot(x_pred_mpc{k}(1,:), x_pred_mpc{k}(2,:), '.-', 'linewidth', 0.5, 'color', [0.35 0.35 1])
	end
end

% Ego robot
for k=1:1:N-1
	plot_obst([x_mpc(1,k); x_mpc(2,k)], [r_ego; r_ego], '-', 1, [0.5 0.5 1])
end
plot_obst([x_mpc(1,N); x_mpc(2,N)], [r_ego; r_ego], '-', 1.25, [0 0 1])

% Trajectory
plot(x_mpc(1,:), x_mpc(2,:), 'b.-', 'linewidth', 1.25, 'markersize', 10)

hold off

axis equal, grid on, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([x_lb(1), x_ub(1)]), ylim([x_lb(2), x_ub(2)])
xlabel('$x$ [m]','interpreter','latex')
ylabel('$y$ [m]','interpreter','latex')
title('Planar motion','interpreter','latex')

f3 = figure(3); set(f3,'WindowStyle','normal','color','w'); f3.Position = [600   200   700   400];

tiledlayout(2,3,'tilespacing','tight','padding','tight')

nexttile(1), hold on

plot(t, x_mpc(1,1:end-1), 'b-')

hold off
grid on, grid minor, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([0, t(end)])
xlabel('$t$ [s]','interpreter','latex')
title('$x$ [m]','interpreter','latex')

nexttile(2), hold on

plot(t, x_mpc(2,1:end-1), 'b-')

hold off
grid on, grid minor, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([0, t(end)])
xlabel('$t$ [s]','interpreter','latex')
title('$y$ [m]','interpreter','latex')

nexttile(3), hold on

plot(t, rad2deg(x_mpc(3,1:end-1)), 'b-')

hold off
grid on, grid minor, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([0, t(end)])
xlabel('$t$ [s]','interpreter','latex')
title('Heading $\theta$ [deg]','interpreter','latex')

nexttile(4), hold on

plot(t, u_mpc(1,:), 'r.-')

hold off
grid on, grid minor, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([0, t(end)])
xlabel('$t$ [s]','interpreter','latex')
title('Velocity $v$ [$\mathrm{m \, s^{-1}}$]','interpreter','latex')

nexttile(5), hold on

plot(t, rad2deg(u_mpc(2,:)), 'r.-')

hold off
grid on, grid minor, box on
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlim([0, t(end)])
xlabel('$t$ [s]','interpreter','latex')
title('Angular velocity $\omega$ [deg/s]','interpreter','latex')

% ====================

function y = quad_form(x,M)
y = x'*M*x;
end

function plot_obst(c, size, line_style, line_width, color)
fplot(@(t) c(1) + size(1)*sin(t), @(t) c(2) + size(2)*cos(t), line_style, 'linewidth', line_width, 'color', color);
end

function plot_rect(c, size, line_style, line_width, color)
plot([c(1)-size(1)/2; c(1)+size(1)/2; c(1)+size(1)/2; c(1)-size(1)/2; c(1)-size(1)/2], ...
	[c(2)-size(2)/2; c(2)-size(2)/2; c(2)+size(2)/2; c(2)+size(2)/2; c(2)-size(2)/2], line_style, 'color', color, 'linewidth', line_width)
end









