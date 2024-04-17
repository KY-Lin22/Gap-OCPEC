function OCPEC = Stewart_FDI_Init_Val()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%ref: equ (7), section 2 in ''Optimal control of systems with discontinuous 
%   differential equations'', 2010, Stewart, D.E. and Anitescu, M.

import casadi.*
%%
% time parameter
timeHorizon = 2; % time horizon T
nStages = 100; % number of discretized stages
timeStep = timeHorizon ./ nStages; % discretization time step

% initial state
x0 = -1; 

% variable 
xDim = 1;
uDim = 0;
lambdaDim = 1;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

% cost function
L_S = x^2;
L_T = (x - 5/3);

% DVI
f_1 = 1; % switch function > 0 
f_2 = 3; % switch function < 0
f = f_1*(1 - lambda) + f_2*lambda; % state equation f

lambdaMax = 1;
lambdaMin = 0;
g = [lambda - lambdaMin;...
    lambdaMax - lambda]; % VI set g >= 0

F = x; % VI function F

VISetType = 'box_constraint'; 
bl = lambdaMin;
bu = lambdaMax;

% inequality constraint G >= 0
G = SX(0,1);
% equality constraint C = 0
C = SX(0,1);

%% create OCPEC instant
OCPEC = OCPEC_Formulation(...
    timeHorizon, nStages, timeStep,...
    x0, ...
    x, u, lambda,...
    L_T, L_S,...
    f, g, F, VISetType, bl, bu,...
    G, C);
end