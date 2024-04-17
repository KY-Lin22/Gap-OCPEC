function OCPEC = Vieira_LCS_State_Jump_1()
% ref: example 7 in ''Quadratic Optimal Control of Linear Complementarity 
%      Systems : First order necessary conditions and numerical analysis''
%      2018, A. Vieira, et. al,
import casadi.*
%%
% time parameter
timeHorizon = 10; % time horizon T
nStages = 1000; % number of discretized stages
timeStep = timeHorizon ./ nStages; % discretization time step

% initial state
x0 = [-2; 1; -1];

% variable 
xDim = 3;
uDim = 1;
lambdaDim = 1;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

% cost function
alpha = 10;
L_S = x' * x + u^2 + alpha * lambda^2;
L_T = 0;

% DVI
f = [0, 1, 0; 0, 0, 1; 0, 0, 0] * x + [0; 0; 1] * u + [0; 0; 1] * lambda; % state equation f
g = lambda;
F = [1, 0, 0] * x + u; % VI function F
VISetType = 'nonnegative_orthant'; 
bl = 0;
bu = inf;

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
    f, g, F, VISetType, bl, bu, ...
    G, C);
end