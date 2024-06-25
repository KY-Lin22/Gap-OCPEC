function FuncObj = create_FuncObj(self)
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
FuncObj = struct();

%% initialize other MX-type symbolic variable and parameter
% dual variable
gamma_h = MX.sym('gamma_h', self.NLP.Dim.h, 1);
gamma_c = MX.sym('gamma_c', self.NLP.Dim.c, 1);
% smooth FB parameter
sigma = MX.sym('sigma', 1, 1);
% vector that collects variable and parameter
Y = [self.NLP.z; gamma_h; gamma_c];
p = [self.NLP.s, sigma];
% dz, dgamma_h, dgamma_c
dz = MX.sym('dz', self.NLP.Dim.z, 1);
dgamma_h = MX.sym('dgamma_h', self.NLP.Dim.h, 1);
dgamma_c = MX.sym('dgamma_c', self.NLP.Dim.c, 1);
dY = [dz; dgamma_h; dgamma_c];

%% IPOPT solver for solving the first parameterized NLP
NLP_Prob = struct('x', self.NLP.z, 'f', self.NLP.J, 'g', [self.NLP.h; self.NLP.c], 'p', self.NLP.s); 
IPOPT_Option = self.Option.IPOPT_Solver;
FuncObj.IPOPT_Solver = nlpsol('Solver', 'ipopt', NLP_Prob, IPOPT_Option);

%% NLP function, Jacobian, and Hessian
% function
FuncObj.J = Function('J', {self.NLP.z}, {self.NLP.J}, {'z'}, {'J'});
FuncObj.h = Function('h', {self.NLP.z}, {self.NLP.h}, {'z'}, {'h'});
FuncObj.c = Function('c', {self.NLP.z, self.NLP.s}, {self.NLP.c}, {'z', 's'}, {'c'});
% Jacobian
J_grad = jacobian(self.NLP.J, self.NLP.z); 
h_grad = jacobian(self.NLP.h, self.NLP.z);
c_grad = jacobian(self.NLP.c, self.NLP.z);
% Hessian
[J_hessian, ~] = hessian(self.NLP.J, self.NLP.z);

%% NLP Lagrangian
% function
LAG = self.NLP.J + gamma_h' * self.NLP.h - gamma_c' * self.NLP.c;
% jacobian
LAG_grad = jacobian(LAG, self.NLP.z); % or LAG_grad = J_grad + gamma_h' * h_grad - gamma_c' * c_grad;
% Hessian
switch self.Option.KKT.HessianApproximation
    case 'Exact'
        % exact Hessian        
        [LAG_hessian, ~] = hessian(LAG, self.NLP.z);
    case 'Gauss_Newton'
        % Gauss-Newton Hessian approximation
        LAG_hessian = J_hessian;
    case 'Quasi_Newton'
        % quasi-Newton Hessian approximation (TODO)
        error('quasi-Newton Hessian approximation method currently is not supported')
    otherwise
        error('specified Hessian approximation method is not supported')
end

%% perturbed system of equation for complementarity between inequality constraint and its dual variable 
% smooth FB function
a = SX.sym('a', self.NLP.Dim.c, 1);
b = SX.sym('b', self.NLP.Dim.c, 1);
sigma_SX = SX.sym('sigma', 1, 1);
psi = sqrt(a.^2 + b.^2 + sigma_SX.^2) - a - b;
psi_FuncObj = Function('psi', {a, b, sigma_SX}, {psi}, {'a', 'b', 'sigma'}, {'psi'});
% PSI function
PSI = psi_FuncObj(self.NLP.c, gamma_c, sigma);
% PSI jacobian
PSI_grad_z = jacobian(PSI, self.NLP.z);
PSI_grad_gamma_c = jacobian(PSI, gamma_c);

%% line search in non-interior-point method
% constraint violation M
M = [self.NLP.h; PSI];
FuncObj.M = Function('M', {Y, p}, {M}, {'Y', 'p'}, {'M'});
% directional derivative of cost 
J_grad_times_dz = J_grad * dz;
FuncObj.J_grad_times_dz = Function('J_grad_times_dz', {Y, dY}, {J_grad_times_dz}, {'Y', 'dY'}, {'J_grad_times_dz'});

%% KKT residual and matrix
% KKT residual
KKT_residual = [LAG_grad'; self.NLP.h; PSI];
FuncObj.KKT_residual = Function('KKT_residual', {Y, p}, {KKT_residual}, {'Y', 'p'}, {'KKT_residual'});
% KKT matrix
nu_h = self.Option.KKT.RegParam.nu_h;
nu_c = self.Option.KKT.RegParam.nu_c;
nu_H = self.Option.KKT.RegParam.nu_H;
KKT_matrix = ...
    [LAG_hessian + nu_H * MX.eye(self.NLP.Dim.z), h_grad',                            -c_grad';...
    h_grad,                                       -nu_h * MX.eye(self.NLP.Dim.h),     MX(self.NLP.Dim.h, self.NLP.Dim.c);...
    PSI_grad_z,                                   MX(self.NLP.Dim.c, self.NLP.Dim.h), PSI_grad_gamma_c - nu_c * MX.eye(self.NLP.Dim.c)];
FuncObj.KKT_matrix = Function('KKT_matrix', {Y, p}, {KKT_matrix}, {'Y', 'p'}, {'KKT_matrix'});

%% sensitivity matrix (w.r.t. parameter)
% lagrangian gradient (w.r.t. s)
LAG_grad_sensitivity = jacobian(LAG_grad, self.NLP.s);
% PSI (w.r.t. s and sigma)
PSI_sensitivity_s = jacobian(PSI, self.NLP.s);
PSI_sensitivity_sigma = jacobian(PSI, sigma);
% formulate sensitivity matrix
sensitivity_matrix = ...
    [LAG_grad_sensitivity,  MX(self.NLP.Dim.z, 1);...
    MX(self.NLP.Dim.h, 1),  MX(self.NLP.Dim.h, 1);...
    PSI_sensitivity_s,      PSI_sensitivity_sigma];
FuncObj.sensitivity_matrix = Function('sensitivity_matrix', {Y, p}, {sensitivity_matrix}, {'Y', 'p'}, {'sensitivity_matrix'});

end