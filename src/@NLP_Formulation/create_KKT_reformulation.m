function [KKT_stationarity_func, KKT_complementarity_func] = create_KKT_reformulation(self, OCPEC)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

import casadi.*
x = SX.sym('x', OCPEC.Dim.x, 1);
u = SX.sym('u', OCPEC.Dim.u, 1);
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);
zeta = SX.sym('zeta', OCPEC.Dim.g, 1); % dual variable for VI set K function g >= 0

s = SX.sym('s', 1, 1); % relaxation parameter for complementarity condition

%% KKT stationarity (equality)
F = OCPEC.FuncObj.F(x, u, lambda);
g_grad = OCPEC.FuncObj.g_grad(lambda);

KKT_stationarity = F - g_grad' * zeta;

KKT_stationarity_func = Function('KKT_stationarity_func',...
    {x, u, lambda, zeta}, {KKT_stationarity},...
    {'x', 'u', 'lambda', 'zeta'}, {'KKT_stationarity'});

%% KKT complementarity (inequality)
g = OCPEC.FuncObj.g(lambda);

switch self.KKT_complementarity_relaxation_strategy
    case 'Scholtes'
        KKT_complementarity = [zeta; g; s - zeta .* g];
    case 'Lin_Fukushima'
        KKT_complementarity = [s^2 - zeta .* g; (zeta + s) .* (g + s) - s^2];
    case 'Kadrani'
        KKT_complementarity = [zeta + s; g + s; - (zeta - s) .* (g - s)];
end
KKT_complementarity_func = Function('KKT_complementarity_func',...
    {lambda, zeta, s}, {KKT_complementarity},  {'lambda', 'zeta', 's'}, {'KKT_complementarity'});

end