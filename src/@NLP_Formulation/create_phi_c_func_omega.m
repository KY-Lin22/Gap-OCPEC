function phi_c_func = create_phi_c_func_omega(self, OCPEC, param_c)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

import casadi.*
% strongly convex function d and its derivative
[d_func, d_grad, ~] = self.create_strongly_convex_func(OCPEC);
% variable
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);
eta = SX.sym('eta', OCPEC.Dim.lambda, 1); 
omega_c = SX.sym('omega_c', OCPEC.Dim.lambda, 1);
% formulate phi_c
p_c = d_func(lambda) - d_func(omega_c) + d_grad(lambda) * (omega_c - lambda);
phi_c = eta' * (lambda - omega_c) + param_c * p_c;
phi_c_func = Function('phi_c_func', {lambda, eta, omega_c}, {phi_c}, {'lambda', 'eta', 'omega_c'}, {'phi_c'});

end