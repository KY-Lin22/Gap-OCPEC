function FuncObj = create_FuncObj(self)
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
FuncObj = struct();

%% initialize other MX-type symbolic variable and parameter
% auxiliary variable
v_c = MX.sym('v_c', self.NLP.Dim.c, 1);
% dual variable
gamma_h = MX.sym('gamma_h', self.NLP.Dim.h, 1);
gamma_c = MX.sym('gamma_c', self.NLP.Dim.c, 1);
% vector that collects variable
Y = [self.NLP.z; v_c; gamma_h; gamma_c];

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
%J_grad = jacobian(self.NLP.J, self.NLP.z); 
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
switch self.Option.KKT.Hessian_approximation
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

%% system of equation for complementarity between inequality constraint auxiliary variable and dual variable 
% FB function
a_SX = SX.sym('a', self.NLP.Dim.c, 1);
b_SX = SX.sym('b', self.NLP.Dim.c, 1);
psi_SX = sqrt(a_SX.^2 + b_SX.^2) - a_SX - b_SX;
psi_FuncObj = Function('psi', {a_SX, b_SX}, {psi_SX}, {'a', 'b'}, {'psi'});
% PSI function
PSI = psi_FuncObj(v_c, gamma_c);
% PSI jacobian
PSI_grad_v_c = jacobian(PSI, v_c);
PSI_grad_gamma_c = jacobian(PSI, gamma_c);

%% KKT residual and matrix
% KKT residual
KKT_residual = [LAG_grad'; self.NLP.h; self.NLP.c - v_c; PSI];
FuncObj.KKT_residual = Function('KKT_residual', {Y, self.NLP.s}, {KKT_residual}, {'Y', 's'}, {'KKT_residual'});
% KKT matrix
nu_h = self.Option.KKT.RegParam.nu_h;
nu_c = self.Option.KKT.RegParam.nu_c;
nu_H = self.Option.KKT.RegParam.nu_H;
KKT_matrix = ...
    [LAG_hessian + nu_H * MX.eye(self.NLP.Dim.z), MX(self.NLP.Dim.z, self.NLP.Dim.c),           h_grad',                                -c_grad';...
    h_grad,                                       MX(self.NLP.Dim.h, self.NLP.Dim.c),           -nu_h * MX.eye(self.NLP.Dim.h),     MX(self.NLP.Dim.h, self.NLP.Dim.c);...
    c_grad,                                       -MX.eye(self.NLP.Dim.c),                      MX(self.NLP.Dim.c, self.NLP.Dim.h),     MX(self.NLP.Dim.c, self.NLP.Dim.c);...
    MX(self.NLP.Dim.c, self.NLP.Dim.z),           PSI_grad_v_c - nu_c * MX.eye(self.NLP.Dim.c), MX(self.NLP.Dim.c, self.NLP.Dim.h),     PSI_grad_gamma_c - nu_c * MX.eye(self.NLP.Dim.c)];

%% sensitivity matrix (w.r.t. parameter)
% lagrangian gradient (w.r.t. s)
LAG_grad_sensitivity = jacobian(LAG_grad, self.NLP.s);
% inequality constraint (w.r.t. s)
c_sensitivity_s = jacobian(self.NLP.c, self.NLP.s);
% formulate sensitivity matrix
sensitivity_matrix = ...
    [LAG_grad_sensitivity;...
    MX(self.NLP.Dim.h, 1);...
    c_sensitivity_s;...
    MX(self.NLP.Dim.c, 1)];

%% differential equation for s
epsilon_s = self.Option.Continuation.epsilon_s;
s_End = self.Option.Continuation.s_End;

s_dot = -epsilon_s*(self.NLP.s - s_End);
FuncObj.s_dot = Function('s_dot', {self.NLP.s}, {s_dot}, {'s'}, {'s_dot'});

%% differential equation for Y
epsilon_T = self.Option.Continuation.epsilon_T;
Y_dot = KKT_matrix\(-epsilon_T * KKT_residual - sensitivity_matrix * s_dot);
FuncObj.Y_dot = Function('Y_dot', {Y, self.NLP.s}, {Y_dot}, {'Y', 's'}, {'Y_dot'});

end