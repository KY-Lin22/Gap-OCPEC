function [KKT_error_primal, KKT_error_dual, KKT_error_complementary] = evaluate_KKT_error(self, Y, LAG_grad_T, h, c)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here

% extract dual variable
Y_Node = cumsum([self.NLP.Dim.z, self.NLP.Dim.h, self.NLP.Dim.c]);
gamma_h = Y(Y_Node(1) + 1 : Y_Node(2), 1);
gamma_c = Y(Y_Node(2) + 1 : Y_Node(3), 1);

% scaling parameter
KKT_scaling_max = self.Option.NIP.KKT_scaling_max;
scaling_dual = max([KKT_scaling_max, norm([gamma_h; gamma_c], 1)/(self.NLP.Dim.h + self.NLP.Dim.c)])/KKT_scaling_max;
scaling_complementary = max([KKT_scaling_max, norm(gamma_c, 1)/self.NLP.Dim.c])/KKT_scaling_max;

% KKT error
KKT_error_primal = norm([h; min([zeros(self.NLP.Dim.c, 1), c], [], 2)], inf);
KKT_error_dual = norm([LAG_grad_T; min([zeros(self.NLP.Dim.c, 1), gamma_c], [], 2)], inf)/scaling_dual;
KKT_error_complementary = norm(c .* gamma_c, inf)/scaling_complementary;

end