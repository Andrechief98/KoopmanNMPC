clc
clear variables
close all

% Nonlinear system
nx = 4;
nu = 2;

x = sym('x',[nx,1],'real');
u = sym('u',[nu,1],'real');

f = [x(3)*cos(x(4));
	x(3)*sin(x(4));
	u(1);
	u(2)];

% Initial set of observables

% phi = x;
phi = [x; x(1)^2; x(2)^2];

%% Compute observables

N_obs = 18; % Max n. of observ. to be generated

i = 1;
while length(phi) < N_obs && i <= length(phi)

	z_dot_curr = gradient(phi(i),x)'*f; % Koopman infinitesimal generator

	z_dot_curr = expand(subs(z_dot_curr, sin(x(4))^2, 1-cos(x(4))^2)); % Replace sin(x)^2 = 1-cos(x)^2

	[~, phi_curr] = coeffs(z_dot_curr);
	phi_curr = subs(phi_curr,u,ones(nu,1));
	
	for j=1:1:length(phi_curr)
		new_obs = logical(phi == phi_curr(j));
	
		if max(new_obs) == 0 && ~isSymType(phi_curr(j),'number')
			phi = [phi; phi_curr(j)];
		end
	end

	i=i+1;

end

if length(phi) > N_obs
	phi = phi(1:N_obs);
	i = i-1;
end

z_dot_nl = []; % Nonlinear lifted equations

while i <= N_obs && i <= length(phi)

	z_dot_nl_curr = gradient(phi(i),x)'*f; % Koopman infinitesimal generator

	z_dot_nl_curr = expand(subs(z_dot_nl_curr, sin(x(4))^2, 1-cos(x(4))^2)); % Replace sin(x)^2 = 1-cos(x)^2

	z_dot_nl = [z_dot_nl; z_dot_nl_curr];

	i=i+1;

end

phi

z_dot_nl

if isempty(z_dot_nl)
	fprintf('Finite set of observables found (%d observables). No nonlinear lifted equations\n', length(phi))
else
	fprintf('%d observables found. %d nonlinear lifted equations\n', length(phi), length(z_dot_nl))
end

%% Compute lifted system matrices

nz = length(phi);
nz_nl = length(z_dot_nl);

nz_lin = nz - nz_nl;

Ap = zeros(nz_lin,nz);

B0p = zeros(nz_lin,nu);
Bp = cell(nu,1);
for i=1:1:nu
	Bp{i} = zeros(nz_lin,nz);
end

for i=1:1:nz_lin

	z_dot_curr = gradient(phi(i),x)'*f; % Koopman infinitesimal generator

	z_dot_curr = expand(subs(z_dot_curr, sin(x(4))^2, 1-cos(x(4))^2)); % Replace sin(x)^2 = 1-cos(x)^2

	[coeffs_curr, phi_curr] = coeffs(z_dot_curr);

	% State matrix A
	for j=1:1:nz
		ind_phi = logical(phi_curr == phi(j));

		if max(ind_phi) == 1
			Ap(i,j) = double(coeffs_curr(ind_phi));
		end
	end

	% Input matrix B0
	for j=1:1:nu
		ind_u = logical(phi_curr == u(j));

		if max(ind_u) == 1
			B0p(i,j) = double(coeffs_curr(ind_u));
		end
	end

	% Bilinear matrices Bk, k=1,...,nu
	for k=1:1:nu
		phi_u = phi*u(k);

		for j=1:1:nz
			ind_phi_u = logical(phi_curr == phi_u(j));

			if max(ind_phi_u) == 1
				Bp{k}(i,j) = double(coeffs_curr(ind_phi_u));
			end
		end
	end

end

%% Compute lifted system matrices (nonlinear part)

As = zeros(nz_nl,nz);

B0s = zeros(nz_nl,nu);
Bs = cell(nu,1);
for i=1:1:nu
	Bs{i} = zeros(nz_nl,nz);
end

z = sym('z',[nx,1],'real');

z_dot_nl_res = sym(zeros(nz_nl,1));

for i=1:1:nz_nl

	z_dot_curr = z_dot_nl(i);

	[coeffs_curr, phi_curr] = coeffs(z_dot_curr);

	% State matrix Ap
	for j=1:1:nz
		ind_phi = logical(phi_curr == phi(j));

		if max(ind_phi) == 1
			As(i,j) = double(coeffs_curr(ind_phi));
			z_dot_curr = z_dot_curr - As(i,j)*phi(j);
		end
	end

	% Input matrix B0
	for j=1:1:nu
		ind_u = logical(phi_curr == u(j));

		if max(ind_u) == 1
			B0s(i,j) = double(coeffs_curr(ind_u));
			z_dot_curr = z_dot_curr - B0s(i,j)*u(j);
		end
	end

	% Bilinear matrices Bk, k=1,...,nu
	for k=1:1:nu
		phi_u = phi*u(k);

		for j=1:1:nz
			ind_phi_u = logical(phi_curr == phi_u(j));

			if max(ind_phi_u) == 1
				Bs{k}(i,j) = double(coeffs_curr(ind_phi_u));
				z_dot_curr = z_dot_curr - Bs{k}(i,j)*phi_u(j);
			end
		end
	end

	z_dot_nl_res(i) = z_dot_curr;

end

z_dot_nl_res = subs(z_dot_nl_res,x,z)

%% Linearize nonlinear residual

Jz_nl = jacobian(z_dot_nl_res,z);
Jz(z) = [sym(zeros(nz_lin,nx)), sym(zeros(nz_lin,nz-nx));
	subs(Jz_nl,u,zeros(nu,1)), sym(zeros(nz_nl,nz-nx))]

Ju_nl = jacobian(z_dot_nl_res,u);
Ju(z) = [sym(zeros(nz_lin,nu));
	subs(Ju_nl,u,zeros(nu,1))]

bs_nl = z_dot_nl_res - Jz_nl*z;
bs(z) = [sym(zeros(nz_lin,1));
	subs(bs_nl,u,zeros(nu,1))]















