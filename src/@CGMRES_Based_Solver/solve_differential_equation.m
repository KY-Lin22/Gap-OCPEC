function [Y_dot, Info] = solve_differential_equation(self, Y, p, p_dot, epsilon)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
timeStart = tic;
% Y node (z, gamma_h, gamma_c)
Y_Node = cumsum([self.NLP.Dim.z, self.NLP.Dim.h, self.NLP.Dim.c]);
% extract variable and parameter
z       = Y(            1 : Y_Node(1), 1);
gamma_c = Y(Y_Node(2) + 1 : Y_Node(3), 1);
s = p(1);
sigma = p(2);

% KKT residual, 
KKT_Residual = self.evaluate_KKT_Residual(Y, p);

% sensitivity matrix
sensitivity_Matrix = self.evaluate_sensitivity_Matrix(Y, p);

% KKT matrix 
c = full(self.FuncObj.c(z, s));
h_grad = sparse(self.FuncObj.h_grad(z));
c_grad = sparse(self.FuncObj.c_grad(z, s));
% Lagrangian Hessian
switch self.Option.HessianApproximation
    case 'Gauss_Newton'
        LAG_hessian = sparse(self.FuncObj.J_hessian(z));
    otherwise
        error('specified Hessian approximation method is not supported')
end
PSI_grad_c = sparse(self.FuncObj.PSI_grad_c(c, gamma_c, sigma));
PSI_grad_gamma_c = sparse(self.FuncObj.PSI_grad_gamma_c(c, gamma_c, sigma));
KKT_Matrix = self.evaluate_KKT_Matrix(h_grad, c_grad, LAG_hessian, PSI_grad_c, PSI_grad_gamma_c);

% Y_dot
Y_dot = KKT_Matrix\(-epsilon * KKT_Residual - sensitivity_Matrix * p_dot);

% Info
Info.time = toc(timeStart);
Info.rho = 0;
end