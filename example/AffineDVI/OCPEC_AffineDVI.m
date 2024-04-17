function OCPEC = OCPEC_AffineDVI()
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
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
xWeight = [20; 20];
uWeight = 1;
L_S = 0.5 * x'*diag(xWeight)*x ...
    + 0.5 * u'*diag(uWeight)*u;
L_T = 0;

% DVI
lambdaMax = 1;
lambdaMin = -1;

f = [1, -3; -8, 10] * x + [4; 8] * u + [-3; -1] * lambda; % state equation f
g = [lambda - lambdaMin;...
    lambdaMax - lambda]; % VI set g >= 0
F = [1, -3] * x + 3 * u + 5 * lambda; % VI function F
VISetType = 'box_constraint'; 
bl = lambdaMin;
bu = lambdaMax;

% inequality constraint G >= 0
xMax = [2; 2];
xMin = [-2; -2];
uMax = 2;
uMin = -2;

G = [xMax - x;...
    x - xMin;...
    uMax - u;...
    u - uMin];
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

