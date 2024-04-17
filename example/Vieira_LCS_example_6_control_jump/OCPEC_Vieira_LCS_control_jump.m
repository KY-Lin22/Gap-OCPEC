function OCPEC = OCPEC_Vieira_LCS_control_jump()
% ref: example 6 in ''Quadratic Optimal Control of Linear Complementarity 
%      Systems : First order necessary conditions and numerical analysis''
%      2018, A. Vieira, et. al,
import casadi.*
%%
% time parameter
timeHorizon = 1; % time horizon T
nStages = 100; % number of discretized stages
timeStep = timeHorizon ./ nStages; % discretization time step

% initial state
x0 = [-0.5; -1]; 

% variable 
xDim = 2;
uDim = 1;
lambdaDim = 1;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

% cost function
L_S = x' * x + u^2;
L_T = 0;

% DVI
f = [1, -3; -8, 10] * x + [4; 8] * u + [-3; -1] * lambda; % state equation f
g = lambda;
F = [1, -3] * x + 3 * u + 5 * lambda; % VI function F
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
    x0,...
    x, u, lambda,...
    L_T, L_S,...
    f, g, F, VISetType, bl, bu,...
    G, C);
end

