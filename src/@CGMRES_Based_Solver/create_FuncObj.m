function FuncObj = create_FuncObj(self, NLP)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
FuncObj = struct();

%% function
FuncObj.J = Function('J', {NLP.z}, {NLP.J}, {'z'}, {'J'});
FuncObj.h = Function('h', {NLP.z}, {NLP.h}, {'z'}, {'h'});
FuncObj.c = Function('c', {NLP.z, NLP.s}, {NLP.c}, {'z', 's'}, {'c'});

%% Jacobian
J_grad = jacobian(NLP.J, NLP.z); 
h_grad = jacobian(NLP.h, NLP.z);
c_grad = jacobian(NLP.c, NLP.z);
FuncObj.J_grad = Function('J_grad', {NLP.z}, {J_grad}, {'z'}, {'J_grad'});
FuncObj.h_grad = Function('h_grad', {NLP.z}, {h_grad}, {'z'}, {'h_grad'});
FuncObj.c_grad = Function('c_grad', {NLP.z, NLP.s}, {c_grad}, {'z', 's'}, {'c_grad'});

%% Hessian
[J_hessian, ~] = hessian(NLP.J, NLP.z);
FuncObj.J_hessian = Function('J_hessian', {NLP.z}, {J_hessian}, {'z'}, {'J_hessian'});

%% quasi-Newton Hessian approximation (TODO)


%% smooth FB function (element-wise)
% formulate smooth FB function (scalar)
a = SX.sym('a', 1, 1);
b = SX.sym('b', 1, 1);
sigma = SX.sym('sigma', 1, 1);
psi = sqrt(a^2 + b^2 + sigma^2) - a - b;
psi_FuncObj = Function('psi', {a, b, sigma}, {psi}, {'a', 'b', 'sigma'}, {'psi'});

% formulate smooth FB function (element-wise) for inequality constraint c
gamma_c = SX.sym('gamma_c', NLP.Dim.c, 1);
xi_c = SX.sym('xi_c', NLP.Dim.c, 1);
psi_FuncObj_map_c = psi_FuncObj.map(NLP.Dim.c);
PSI = (psi_FuncObj_map_c(gamma_c', xi_c', sigma))';

% formulate smooth FB jacobian (element-wise) for inequality constraint c
PSI_grad_gamma_c = jacobian(PSI, gamma_c);
PSI_grad_xi_c = jacobian(PSI, xi_c);

% create function object
FuncObj.PSI = Function('PSI', {gamma_c, xi_c, sigma}, {PSI}, {'gamma_c', 'xi_c', 'sigma'}, {'PSI'});
FuncObj.PSI_grad_gamma_c = Function('PSI_grad_gamma_c', {gamma_c, xi_c, sigma}, {PSI_grad_gamma_c},...
    {'gamma_c', 'xi_c', 'sigma'}, {'PSI_grad_gamma_c'});
FuncObj.PSI_grad_xi_c = Function('PSI_grad_xi_c', {gamma_c, xi_c, sigma}, {PSI_grad_xi_c},...
    {'gamma_c', 'xi_c', 'sigma'}, {'PSI_grad_xi_c'});

end