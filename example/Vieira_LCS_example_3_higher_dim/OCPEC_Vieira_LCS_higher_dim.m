function OCPEC = OCPEC_Vieira_LCS_higher_dim()
% ref: example 3 in ''Quadratic Optimal Control of Linear Complementarity 
%      Systems : First order necessary conditions and numerical analysis''
%      2018, A. Vieira, et. al,
import casadi.*
%%
% time parameter
TimeHorizon = 1; % time horizon T
nStages = 100; % number of discretized stages
timeStep = TimeHorizon ./ nStages; % discretization time step

% initial and reference state
x0 = [-0.5; 1]; % initial state
xRef = [0; 0]; % ref state

% variable and their bounds
xDim = 2;
uDim = 2;
lambdaDim = 2;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

xMax = [inf; inf];
xMin = [-inf; -inf];
uMax = [inf; inf];
uMin = [-inf; -inf];
lambdaMax = [inf; inf];
lambdaMin = [0; 0];

% cost function
L_S = x' * x + 25* (u' * u);
L_T = 0;

% DVI
f = [1, 2; 2, 1] * x + [1, 3; 2, 1] * u + [-1, 1; -1, 1] * lambda; % state equation f
g = lambda;
F = [3, -1; -2, 0] * x + [1, -1; -1, 2] * u + lambda; % VI function F
VISetType = 'nonnegative_orthant'; 
% inequality constraint G >= 0
G = SX(0,1);
% equality constraint C = 0
C = SX(0,1);
%% create OCPEC instant
OCPEC = OCPEC_Formulation(...
    TimeHorizon, nStages, timeStep,...
    x0, xRef,...
    x, u, lambda,...
    xMax, xMin, uMax, uMin, lambdaMax, lambdaMin,...
    L_T, L_S,...
    f, g, F, VISetType,...
    G, C);
end