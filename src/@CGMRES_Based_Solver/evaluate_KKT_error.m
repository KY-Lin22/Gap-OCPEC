function [KKT_error_primal, KKT_error_dual, KKT_error_complementary] = ...
    evaluate_KKT_error(self, gamma_h, gamma_c, h, c, LAG_grad_z)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
% scaling parameter
scaling_dual = max([self.Option.KKT_scaling_max,...
    norm([gamma_h; gamma_c], 1)/(self.NLP.Dim.h + self.NLP.Dim.c)])/self.Option.KKT_scaling_max;
scaling_complementary = max([self.Option.KKT_scaling_max,...
    norm(gamma_c, 1)/self.NLP.Dim.c])/self.Option.KKT_scaling_max;

% KKT error
KKT_error_primal = norm([h; min([zeros(self.NLP.Dim.c, 1), c], [], 2)], inf);
KKT_error_dual = norm([LAG_grad_z'; min([zeros(self.NLP.Dim.c, 1), gamma_c], [], 2)], inf)/scaling_dual;
KKT_error_complementary = norm(c .* gamma_c, inf)/scaling_complementary;

end