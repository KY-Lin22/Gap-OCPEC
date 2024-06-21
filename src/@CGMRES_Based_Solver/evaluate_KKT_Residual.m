function KKT_Residual = evaluate_KKT_Residual(self, Y, p)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

% Y node (z, gamma_h, gamma_c)
Y_Node = cumsum([self.NLP.Dim.z, self.NLP.Dim.h, self.NLP.Dim.c]);
% extract variable and parameter
z       = Y(            1 : Y_Node(1), 1);
gamma_h = Y(Y_Node(1) + 1 : Y_Node(2), 1);
gamma_c = Y(Y_Node(2) + 1 : Y_Node(3), 1);
s = p(1);
sigma = p(2);

% evaluate function and Jacobian
h = full(self.FuncObj.h(z));
c = full(self.FuncObj.c(z, s));
J_grad = full(self.FuncObj.J_grad(z));
h_grad = sparse(self.FuncObj.h_grad(z));
c_grad = sparse(self.FuncObj.c_grad(z, s));
PSI = full(self.FuncObj.PSI(c, gamma_c, sigma));
LAG_grad_z = J_grad + gamma_h' * h_grad - gamma_c' * c_grad;

% KKT residual
KKT_Residual = [LAG_grad_z'; h; PSI];
end