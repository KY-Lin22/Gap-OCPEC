function KKT_Matrix = evaluate_KKT_Matrix(self, h_grad, c_grad, LAG_hessian, PSI_grad_c, PSI_grad_gamma_c)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
%%
% load regularization parameter
nu_h = self.Option.RegParam.nu_h;
nu_c = self.Option.RegParam.nu_c;
nu_H = self.Option.RegParam.nu_H;
% Y node (z, gamma_h, gamma_c)
Y_Node = cumsum([self.NLP.Dim.z, self.NLP.Dim.h, self.NLP.Dim.c]);

%% extract nonzero index and value
% Hessian + nu_H * eye(Dim.z)
[i_hessian, j_hessian, s_hessian] = find(LAG_hessian + nu_H * speye(self.NLP.Dim.z));
% h_grad'
[i_h_grad_T, j_h_grad_T, s_h_grad_T] = find(h_grad');
j_h_grad_T = j_h_grad_T + Y_Node(1);
% -c_grad'
[i_nega_c_grad_T, j_nega_c_grad_T, s_nega_c_grad_T] = find(-c_grad');
j_nega_c_grad_T = j_nega_c_grad_T + Y_Node(2);
% h_grad
[i_h_grad, j_h_grad, s_h_grad] = find(h_grad);
i_h_grad = i_h_grad + Y_Node(1);
% -nu_h * eye(Dim.h)
i_nu_h = Y_Node(1) + 1 : Y_Node(2);
j_nu_h = Y_Node(1) + 1 : Y_Node(2);
s_nu_h = - nu_h * ones(1, self.NLP.Dim.h);
% PSI_grad_c * c_grad
[i_PSI_grad_c_c_grad, j_PSI_grad_c_c_grad, s_PSI_grad_c_c_grad] = find(PSI_grad_c * c_grad);
i_PSI_grad_c_c_grad = i_PSI_grad_c_c_grad + Y_Node(2);
% PSI_grad_gamma_c - nu_c * eye(Dim.c)
[i_PSI_grad_gamma_c, j_PSI_grad_gamma_c, s_PSI_grad_gamma_c] = find(PSI_grad_gamma_c - nu_c * speye(self.NLP.Dim.c));
i_PSI_grad_gamma_c = i_PSI_grad_gamma_c + Y_Node(2);
j_PSI_grad_gamma_c = j_PSI_grad_gamma_c + Y_Node(2);

%% assemble KKT matrix
i_KKT = [i_hessian; i_h_grad_T; i_nega_c_grad_T; i_h_grad; i_nu_h'; i_PSI_grad_c_c_grad; i_PSI_grad_gamma_c];
j_KKT = [j_hessian; j_h_grad_T; j_nega_c_grad_T; j_h_grad; j_nu_h'; j_PSI_grad_c_c_grad; j_PSI_grad_gamma_c];
s_KKT = [s_hessian; s_h_grad_T; s_nega_c_grad_T; s_h_grad; s_nu_h'; s_PSI_grad_c_c_grad; s_PSI_grad_gamma_c];
KKT_Matrix = sparse(i_KKT, j_KKT, s_KKT, Y_Node(3), Y_Node(3), length(s_KKT));

end