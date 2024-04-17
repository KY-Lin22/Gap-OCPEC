function OCPEC = OCPEC_CartPoleWithFriction()
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
%%
% time parameter
timeHorizon = 3; % time horizon T
nStages = 300; % number of discretized stages
timeStep = timeHorizon ./ nStages; % discretization time step

% initial state
x0 = [1; 0/180*pi; 0; 0]; 

% variable
xDim = 4;
uDim = 1;
lambdaDim = 1;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

% cost function
xRef = [1; 180/180*pi; 0; 0]; % ref state
xWeight_S = [1; 100; 1; 1];
uWeight_S = 1;
lambdaWeight_S = 1;
L_S = 0.5 * (x - xRef)'*diag(xWeight_S)*(x - xRef)...
    + 0.5 * u'*diag(uWeight_S)*u ...
    + 0.5 * lambda'*diag(lambdaWeight_S)*lambda;

xWeight_T = [1; 100; 10; 20];
L_T = 0.5 * (x - xRef)'*diag(xWeight_T)*(x - xRef);

% DVI
mass = [1; 0.1];
linkLength = 1;
gravityParam = 9.8;
M = [mass(1) + mass(2),                 mass(2) * linkLength * cos(x(2));...
     mass(2) * linkLength * cos(x(2)),  mass(2) * linkLength^2];
C_Mat = [0,   -mass(2) * linkLength * x(4) * sin(x(2));...
     0,   0]; 
G_Mat = [0;...
     -mass(2) * gravityParam * linkLength * sin(x(2))];
Bu = [u(1);...
      0]; 
LAMBDA = [lambda(1);...
     0]; % friction bewteen cart and ground  
H = G_Mat + Bu + LAMBDA - C_Mat * [x(3); x(4)];

lambdaMax = 2;
lambdaMin = -2;

f = [x(3:4);...
    inv(M)*H];% xDot = f(x, u, lambda)
g = [lambda - lambdaMin;...
    lambdaMax - lambda]; % VI set g >= 0
F = x(3);
VISetType = 'box_constraint'; 
bl = lambdaMin;
bu = lambdaMax;

% inequality constraint G >= 0
xMax = [5; 240/180*pi; 20; 20];
xMin = [0; -240/180*pi; -20; -20];
uMax = 30;
uMin = -30;
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

