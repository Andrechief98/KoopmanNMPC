function [A,B0,B,phi] = koopman_calc_fun(x,u,f,phi,N_obs)

% Nonlinear system
nu = length(u);

% Compute observables

i = 1;
while length(phi) < N_obs && i <= length(phi)

	z_dot_curr = gradient(phi(i),x)'*f; % Koopman infinitesimal generator

	z_dot_curr = expand(subs(z_dot_curr, sin(x(4))^2, 1-cos(x(4))^2)); % Replace sin(x)^2 = 1-cos(x)^2

	[~,phi_curr] = coeffs(z_dot_curr);
	phi_curr = subs(phi_curr,u,ones(nu,1));
	
	phi_new = [];
	for j=1:1:length(phi_curr)
		new_obs = logical(phi == phi_curr(j));
	
		if max(new_obs) == 0 && ~isSymType(phi_curr(j),'number')
			phi_new = [phi_new; phi_curr(j)];
		end
	end

	phi = [phi; phi_new];

	i=i+1;

end

if i == length(phi)+1
	phi
	fprintf('Finite set of observables found (%d observables)\n', length(phi))
elseif length(phi) >= N_obs
	phi
	error('%d observables found (finite set not found)\n', length(phi))
end

% Compute lifted system matrices

nz = length(phi);

A = zeros(nz,nz);

B0 = zeros(nz,nu);
B = cell(nu,1);
for i=1:1:nu
	B{i} = zeros(nz,nz);
end

for i=1:1:length(phi)

	z_dot_curr = gradient(phi(i),x)'*f; % Koopman infinitesimal generator

	z_dot_curr = expand(subs(z_dot_curr, sin(x(4))^2, 1-cos(x(4))^2)); % Replace sin(x)^2 = 1-cos(x)^2

	[phi_temp,input_terms] = coeffs(z_dot_curr,u);

	% State matrix
	ind_z = logical(input_terms == 1);
	if sum(ind_z) >= 1
		[coeffs_z,phi_curr_z] = coeffs(phi_temp(ind_z));

		ind_obs = logical(phi == phi_curr_z);

		A(i,:) = double(coeffs_z)*ind_obs'; % State terms
	end

	% Input matrices
	for j=1:1:nu
		ind_u = logical(input_terms == u(j));
		if sum(ind_u) >= 1
			[coeffs_u,phi_curr_u] = coeffs(phi_temp(ind_u));

			ind_obs = logical(phi == phi_curr_u);

			B{j}(i,:) = double(coeffs_u)*ind_obs'; % Input bilinear terms

			if any(phi_curr_u == 1)
				B0(i,j) = 1; % Input linear terms
			end

		end
	end

end

end
















