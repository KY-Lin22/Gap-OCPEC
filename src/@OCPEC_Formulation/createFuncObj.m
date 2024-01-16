function FuncObj = createFuncObj(self)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
%% Jacobian
glambda = jacobian(self.g, self.lambda);

%% function object
% cost function
FuncObj.L_S = Function('L_S', {self.x, self.u, self.lambda}, {self.L_S}, {'x', 'u', 'lambda'}, {'L_S'});
FuncObj.L_T = Function('L_T', {self.x}, {self.L_T}, {'x'}, {'L_T'});

% DVI
FuncObj.f = Function('f', {self.x, self.u, self.lambda}, {self.f}, {'x', 'u', 'lambda'}, {'f'});
FuncObj.g = Function('g', {self.lambda}, {self.g}, {'lambda'}, {'g'});
FuncObj.F = Function('F', {self.x, self.u, self.lambda}, {self.F}, {'x', 'u', 'lambda'}, {'F'});

FuncObj.glambda = Function('glambda', {self.lambda}, {glambda}, {'lambda'}, {'glambda'});

% inequality and equality constraint 
FuncObj.G = Function('G', {self.x, self.u}, {self.G}, {'x', 'u'}, {'G'});
FuncObj.C = Function('C', {self.x, self.u}, {self.C}, {'x', 'u'}, {'C'});

end

