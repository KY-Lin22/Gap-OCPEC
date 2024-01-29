function discre_state_equ_func = create_discre_state_equ_func(self, OCPEC)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

import casadi.*

% symbolic variable
xPrev = SX.sym('xPrev', OCPEC.Dim.x, 1);
x = SX.sym('x', OCPEC.Dim.x, 1);
u = SX.sym('u', OCPEC.Dim.u, 1);
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);

% discretized state equation: discre_state_equ = 0
switch self.state_equation_discretization
    case 'implicit_Euler'
        discre_state_equ = xPrev - x + OCPEC.timeStep * OCPEC.FuncObj.f(x, u, lambda);
end

% function object
discre_state_equ_func = Function('discre_state_equ_func', {xPrev, x, u, lambda}, {discre_state_equ},...
    {'xPrev', 'x', 'u', 'lambda'}, {'discre_state_equ'});

end