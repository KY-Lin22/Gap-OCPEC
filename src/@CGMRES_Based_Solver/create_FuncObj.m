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

%% NLP function, Jacobian, and Hessian
% function
FuncObj.J = Function('J', {self.NLP.z}, {self.NLP.J}, {'z'}, {'J'});
FuncObj.h = Function('h', {self.NLP.z}, {self.NLP.h}, {'z'}, {'h'});
FuncObj.c = Function('c', {self.NLP.z, self.NLP.s}, {self.NLP.c}, {'z', 's'}, {'c'});
% Jacobian
J_grad = jacobian(self.NLP.J, self.NLP.z); 
h_grad = jacobian(self.NLP.h, self.NLP.z);
c_grad = jacobian(self.NLP.c, self.NLP.z);
FuncObj.J_grad = Function('J_grad', {self.NLP.z}, {J_grad}, {'z'}, {'J_grad'});
FuncObj.h_grad = Function('h_grad', {self.NLP.z}, {h_grad}, {'z'}, {'h_grad'});
FuncObj.c_grad = Function('c_grad', {self.NLP.z, self.NLP.s}, {c_grad}, {'z', 's'}, {'c_grad'});
% Hessian
[J_hessian, ~] = hessian(self.NLP.J, self.NLP.z);
FuncObj.J_hessian = Function('J_hessian', {self.NLP.z}, {J_hessian}, {'z'}, {'J_hessian'});

%% NLP Lagrangian
% function
LAG = self.NLP.J + gamma_h' * self.NLP.h - gamma_c' * self.NLP.c;
FuncObj.LAG = Function('LAG', {self.NLP.z, gamma_h, gamma_c, self.NLP.s}, {LAG},...
    {'z', 'gamma_h', 'gamma_c', 's'}, {'LAG'});
% jacobian
LAG_grad = jacobian(LAG, self.NLP.z); % or LAG_grad = J_grad + gamma_h' * h_grad - gamma_c' * c_grad;
FuncObj.LAG_grad = Function('LAG_grad', {self.NLP.z, gamma_h, gamma_c, self.NLP.s}, {LAG_grad},...
    {'z', 'gamma_h', 'gamma_c', 's'}, {'LAG_grad'});
% Hessian
switch self.Option.NIP.HessianApproximation
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
FuncObj.LAG_hessian = Function('LAG_hessian', {self.NLP.z, gamma_h, gamma_c, self.NLP.s}, {LAG_hessian},...
    {'z', 'gamma_h', 'gamma_c', 's'}, {'LAG_hessian'});

%% perturbed system of equation for complementarity between inequality constraint and its dual variable 
% smooth FB function
a = SX.sym('a', self.NLP.Dim.c, 1);
b = SX.sym('b', self.NLP.Dim.c, 1);
sigma_SX = SX.sym('sigma', 1, 1);
psi = sqrt(a.^2 + b.^2 + sigma_SX.^2) - a - b;
psi_FuncObj = Function('psi', {a, b, sigma_SX}, {psi}, {'a', 'b', 'sigma'}, {'psi'});
% PSI function
PSI = psi_FuncObj(self.NLP.c, gamma_c, sigma);
FuncObj.PSI = Function('PSI', {self.NLP.z, gamma_c, self.NLP.s, sigma}, {PSI},...
    {'z', 'gamma_c', 's', 'sigma'}, {'PSI'});
% PSI jacobian
PSI_grad_z = jacobian(PSI, self.NLP.z);
PSI_grad_gamma_c = jacobian(PSI, gamma_c);

%% line search
% constraint violation M
M = [self.NLP.h; PSI];
FuncObj.M = Function('M', {Y, p}, {M}, {'Y', 'p'}, {'M'});

%% KKT residual, matrix
% KKT residual
KKT_residual = [LAG_grad'; self.NLP.h; PSI];
FuncObj.KKT_residual = Function('KKT_residual', {Y, p}, {KKT_residual}, {'Y', 'p'}, {'KKT_residual'});
% KKT matrix
nu_h = self.Option.NIP.RegParam.nu_h;
nu_c = self.Option.NIP.RegParam.nu_c;
nu_H = self.Option.NIP.RegParam.nu_H;
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