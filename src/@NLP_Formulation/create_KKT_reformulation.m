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
g_grad = OCPEC.FuncObj.g_grad(lambda);

KKT_stationarity = F - g_grad' * zeta;

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
    case 'Steffensen_Ulbrich'
        % ref: A new relaxation scheme for mathematical programs with equilibrium constraints, Steffensen and Ulbrich
        % theta function
        z = SX.sym('z', OCPEC.Dim.g, 1);
        theta = 2/pi * sin(z*pi/2 + 3*pi/2) + 1;
        theta_func = Function('theta_func', {z}, {theta});
        % phi_SU function
        phi_SU = if_else(abs(z) - s >= 0, abs(z), s * theta_func(z/s));
        phi_SU_func = Function('phi_SU_func', {z, s}, {phi_SU});
        % KKT complementarity
        KKT_complementarity = [zeta; w; phi_SU_func(zeta - w, s) - zeta - w];
    case 'Kanzow_Schwartz'
        % ref: A new regularization method for mathematical programs with complementarity constraints with strong convergence properties, Kanzow, Christian and Schwartz, Alexandra
        % phi_KS
        a = SX.sym('a', 1, 1);
        b = SX.sym('b', 1, 1);
        phi_KS = if_else(a + b >= 0, a*b, -1/2*(a^2 + b^2));
        phi_KS_func = Function('phi_KS_func', {a, b}, {phi_KS});
        % KKT complementarity
        KKT_complementarity = [zeta; w; - phi_KS_func(zeta - s, w - s)];
end
KKT_complementarity_func = Function('KKT_complementarity_func',...
    {zeta, w, s}, {KKT_complementarity},  {'zeta', 'w', 's'}, {'KKT_complementarity'});

end