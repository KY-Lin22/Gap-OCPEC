function [d_func, d_grad, d_hessian] = createStronglyConvexFunction(self, OCPEC)
% create function object about a scalar strongly convex function d(x) and its derivative
% where x is a vector variable having the VI variable dimension
import casadi.*
% vector variable of this convex function
varDim = OCPEC.Dim.lambda;
x = SX.sym('x', varDim, 1); 
% define strongly convex function
switch self.stronglyConvexFuncType
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

