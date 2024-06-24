function sensitivity_Matrix = evaluate_sensitivity_Matrix(self, Y, p)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here
% Y node (z, gamma_h, gamma_c)
Y_Node = cumsum([self.NLP.Dim.z, self.NLP.Dim.h, self.NLP.Dim.c]);
p_Node = cumsum([1, 1]);
% extract variable and parameter
z       = Y(            1 : Y_Node(1), 1);
gamma_c = Y(Y_Node(2) + 1 : Y_Node(3), 1);
s = p(1);
sigma = p(2);

% evaluate function, jacobian and sensitivity
c = full(self.FuncObj.c(z, s));

PSI_grad_c = sparse(self.FuncObj.PSI_grad_c(c, gamma_c, sigma));

LAG_grad_z_T_sensitivity = sparse(self.FuncObj.LAG_grad_z_T_sensitivity(z, gamma_c, s));
c_sensitivity = sparse(self.FuncObj.c_sensitivity(z, s));
PSI_sensitivity = sparse(self.FuncObj.PSI_sensitivity(c, gamma_c, sigma));

%% extract nonzero index and value
% LAG_grad_z w.r.t. s
[i_LAG_grad_z_T_sensitivity, j_LAG_grad_z_T_sensitivity, s_LAG_grad_z_T_sensitivity] = find(LAG_grad_z_T_sensitivity);

% PSI w.r.t. s
[i_PSI_wrt_s, j_PSI_wrt_s, s_PSI_wrt_s] = find(PSI_grad_c * c_sensitivity);
i_PSI_wrt_s = i_PSI_wrt_s + Y_Node(2);

% PSI w.r.t. sigma
[i_PSI_sensitivity, j_PSI_sensitivity, s_PSI_sensitivity] = find(PSI_sensitivity);
i_PSI_sensitivity = i_PSI_sensitivity + Y_Node(2);
j_PSI_sensitivity = j_PSI_sensitivity + p_Node(1);

%% assemble sensitivity matrix
i_sensitivity = [i_LAG_grad_z_T_sensitivity; i_PSI_wrt_s; i_PSI_sensitivity];
j_sensitivity = [j_LAG_grad_z_T_sensitivity; j_PSI_wrt_s; j_PSI_sensitivity];
s_sensitivity = [s_LAG_grad_z_T_sensitivity; s_PSI_wrt_s; s_PSI_sensitivity];

sensitivity_Matrix = sparse(i_sensitivity, j_sensitivity, s_sensitivity, Y_Node(3), p_Node(2), length(s_sensitivity));
end