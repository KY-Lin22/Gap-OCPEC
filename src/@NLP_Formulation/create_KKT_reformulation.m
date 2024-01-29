function [KKT_stationarity_func, KKT_complementarity_func] = create_KKT_reformulation(self, OCPEC)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

import casadi.*
x = SX.sym('x', OCPEC.Dim.x, 1);
u = SX.sym('u', OCPEC.Dim.u, 1);
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);
zeta = SX.sym('zeta', OCPEC.Dim.g, 1); % dual variable for VI set K function g >= 0
w = SX.sym('w', OCPEC.Dim.g, 1); % auxiliary variable for VI set K function g >= 0
s = SX.sym('s', 1, 1); % relaxation parameter for complementarity condition

%% KKT stationarity (equality)
F = OCPEC.FuncObj.F(x, u, lambda);
glambda = OCPEC.FuncObj.glambda(lambda);

KKT_stationarity = F - glambda' * zeta;

KKT_stationarity_func = Function('KKT_stationarity_func',...
    {x, u, lambda, zeta}, {KKT_stationarity},...
    {'x', 'u', 'lambda', 'zeta'}, {'KKT_stationarity'});

%% KKT complementarity (inequality)
switch self.KKT_complementarity_relaxation_strategy
    case 'Scholtes'
        KKT_complementarity = [zeta; w; s - zeta .* w];
    case 'Lin_Fukushima'
        KKT_complementarity = [s^2 - zeta .* w; (zeta + s) .* (w + s) - s^2];
    case 'Kadrani'
        KKT_complementarity = [zeta + s; w + s; - (zeta - s) .* (w - s)];
end
KKT_complementarity_func = Function('KKT_complementarity_func',...
    {zeta, w, s}, {KKT_complementarity},  {'zeta', 'w', 's'}, {'KKT_complementarity'});


end