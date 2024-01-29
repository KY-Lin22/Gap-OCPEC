function [d_func, d_grad, d_hessian] = create_strongly_convex_func(self, OCPEC)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
% vector variable of this convex function
varDim = OCPEC.Dim.lambda;
x = SX.sym('x', varDim, 1); 
% define strongly convex function
switch self.strongly_convex_func
    case 'quadratic'
        d = 0.5 * x' * eye(varDim) * x; 
    case 'general'
        % TODO
end
% gradient and hessian
dx = jacobian(d, x);
[dxx, ~] = hessian(d, x);
% function object
d_func = Function('d_func', {x}, {d}, {'x'}, {'d'});
d_grad = Function('d_grad', {x}, {dx}, {'x'}, {'dx'});
d_hessian = Function('d_hessian', {x}, {dxx}, {'x'}, {'dxx'});
end