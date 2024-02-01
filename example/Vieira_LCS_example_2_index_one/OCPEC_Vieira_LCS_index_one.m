function OCPEC = OCPEC_Vieira_LCS_index_one()
% ref: example 2 in ''Quadratic Optimal Control of Linear Complementarity 
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
uDim = 1;
lambdaDim = 1;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

xMax = [inf; inf];
xMin = [-inf; -inf];
uMax = inf;
uMin = -inf;
lambdaMax = inf;
lambdaMin = 0;

% cost function
L_S = x'*x + u^2;
L_T = 0;

% DVI
f = [0, 1; 0, 0]*x + [0; 1] * u + [-1; 1] * lambda; % state equation f
g = lambda;
F = [-1, 1] * x + u; % VI function F
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