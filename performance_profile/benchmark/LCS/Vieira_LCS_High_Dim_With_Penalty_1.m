function OCPEC = Vieira_LCS_High_Dim_With_Penalty_1()
% ref: example 3 in ''Quadratic Optimal Control of Linear Complementarity 
%      Systems : First order necessary conditions and numerical analysis''
%      2018, A. Vieira, et. al,
import casadi.*
%%
% time parameter
timeHorizon = 1; % time horizon T
nStages = 100; % number of discretized stages
timeStep = timeHorizon ./ nStages; % discretization time step

% initial state
x0 = [-0.5; 1]; 

% variable 
xDim = 2;
uDim = 2;
lambdaDim = 2;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

% cost function
alpha = 10;
L_S = x'*x + 25 * (u' * u) + alpha*  (lambda' * lambda);
L_T = 0;

% DVI
f = [1, 2; 2, 1] * x + [1, 3; 2, 1] * u + [-1, 1; -1, 1] * lambda; % state equation f
g = lambda;
F = [3, -1; -2, 0] * x + [1, -1; -1, 2] * u + lambda; % VI function F
VISetType = 'nonnegative_orthant'; 
bl = [0; 0];
bu = [inf; inf];

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