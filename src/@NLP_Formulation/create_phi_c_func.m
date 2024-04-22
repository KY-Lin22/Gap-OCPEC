function phi_c_func = create_phi_c_func(self, OCPEC, param_c)
% substitute omega expression into phi_c expression, current only for box-cstr and nonnega-orth

import casadi.*
% strongly convex function d and its derivative
[d_func, d_grad, ~] = self.create_strongly_convex_func(OCPEC);
% variable
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);
eta = SX.sym('eta', OCPEC.Dim.lambda, 1); 
% omega_c expression
if strcmp(OCPEC.VISetType, 'finitely_representable') || strcmp(OCPEC.VISetType, 'polyhedral')
    % TODO omega_solver is a function object about a convex or QP optimization solver 
elseif strcmp(OCPEC.VISetType, 'box_constraint') || strcmp(OCPEC.VISetType, 'nonnegative_orthant')
    omega_solver = self.create_omega_solver(OCPEC, param_c);
    omega_c = omega_solver(lambda, eta);
end
% formulate phi_c
p_c = d_func(lambda) - d_func(omega_c) + d_grad(lambda) * (omega_c - lambda);
phi_c = eta' * (lambda - omega_c) + param_c * p_c;
phi_c_func = Function('phi_c_func', {lambda, eta}, {phi_c}, {'lambda', 'eta'}, {'phi_c'});

end