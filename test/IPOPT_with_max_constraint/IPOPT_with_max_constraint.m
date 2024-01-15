clear all
clc
addpath('E:\GitHub\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

% variable
lambda = SX.sym('lambda', 1, 1);
eta = SX.sym('eta', 1, 1);
s = SX.sym('s', 1, 1);

% problem cost and constraint
J = 0.5 * ((lambda - 1)^2 + (eta - 1)^2);

% primal gap
% phi = 0.5 * (eta^2 - (max(0, eta - lambda))^2);
% g = [lambda;...
%     s - phi];

% D gap
b = 2;
a = 0.3;
phi = (b-a)/(2*a*b)*eta^2 - 1/(2*a)*(max(0, eta-a*lambda))^2 + 1/(2*b)*(max(0, eta-b*lambda))^2;
g = s - phi;

% create solver
Prob = struct('x', [lambda; eta], 'f', J, 'g', g, 'p', s);

Option = struct;
Option.ipopt.tol = 1e-6;

x_Init = [0; 1];
s_Init = 1e-6;

solver = nlpsol('solver', 'ipopt', Prob, Option);

% 
solution = solver('x0', x_Init, 'p', s_Init,...
    'lbg', zeros(size(g, 1), 1), 'ubg', inf*ones(size(g, 1), 1));

%
x_Opt = full(solution.x);
J_opt = full(solution.f);