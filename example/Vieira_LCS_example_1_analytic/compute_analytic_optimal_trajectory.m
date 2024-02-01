function [X_analytic_Opt, U_analytic_Opt, LAMBDA_analytic_Opt] = compute_analytic_optimal_trajectory(OCPEC)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
param.a = 3;
param.b = -0.5;
param.d = 1;
param.e = -2;
param.f = 3;
timeAxis = 0 : OCPEC.timeStep : OCPEC.nStages * OCPEC.timeStep;

% compute gamma based on x_0
if param.e * param.f * OCPEC.x0 <= 0
    gamma = param.a^2 + param.f^2;
else
    gamma = param.a^2 + (param.f - param.b * param.e / param.d)^2;
end

% compute p_0 based on gamma and x_0
p_0 = - (OCPEC.x0 * (exp(2*sqrt(gamma)*OCPEC.TimeHorizon) - 1)) / ...
    ( (sqrt(gamma) - param.a) * exp(2*sqrt(gamma) * OCPEC.TimeHorizon) + sqrt(gamma) + param.a);

% symbolic time variable t
t = SX.sym('t', 1, 1);

p_first_term = (sqrt(gamma) - param.a) * exp(sqrt(gamma) * t) ...
             + (sqrt(gamma) + param.a) * exp(- sqrt(gamma) * t);
p_second_term = exp(sqrt(gamma)*t) - exp(-sqrt(gamma)*t);

p = 1/(2*sqrt(gamma)) * (p_first_term * p_0 + p_second_term * OCPEC.x0);

dpdt = jacobian(p, t);

p_func = Function('p_func', {t}, {p}, {'t'}, {'p'});
dpdt_func = Function('dpdt_func', {t}, {dpdt}, {'t'}, {'dpdt'});

p_func_map = p_func.map(OCPEC.nStages);
dpdt_func_map = dpdt_func.map(OCPEC.nStages);

dpdt_value = full(dpdt_func_map(timeAxis(2:end)));
p_value = full(p_func_map(timeAxis(2:end)));

% anlytical optimal state trajectory
X_analytic_Opt = dpdt_value + param.a * p_value;

% anlytical optimal control trajectory
if param.e * param.f * OCPEC.x0 <= 0
    U_analytic_Opt = param.f * p_value;
else
    U_analytic_Opt = (param.f - param.e * param.b/ param.d) * p_value;
end

% anlytical optimal algebraic trajectory
LAMBDA_analytic_Opt = 1/param.d*max(zeros(1, length(U_analytic_Opt)), -param.e * U_analytic_Opt);
end