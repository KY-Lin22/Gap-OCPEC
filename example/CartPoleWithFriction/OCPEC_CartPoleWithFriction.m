function OCPEC = OCPEC_CartPoleWithFriction()
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
%%
% time parameter
TimeHorizon = 3; % time horizon T
nStages = 300; % number of discretized stages
timeStep = TimeHorizon ./ nStages; % discretization time step

% initial and reference state
x0 = [0; 0/180*pi; 0; 0]; % initial state
xRef = [0; 180/180*pi; 0; 0]; % ref state

% variable and their bounds
xDim = 4; % cart position, pole angle, car vel, pole vel
uDim = 1;
lambdaDim = 1;
x = SX.sym('x', xDim, 1);
u = SX.sym('u', uDim, 1);
lambda = SX.sym('lambda', lambdaDim, 1);

xMax = [5; 240/180*pi; 20; 20];
xMin = [0; -240/180*pi; -20; -20];
uMax = 30;
uMin = -30;

% cost function
xWeight_S = [1; 1; 1; 1];
uWeight_S = 1;
lambdaWeight = 0.001;
L_S = 0.5 * xWeight_S(1) * (x(1) - xRef(1))^2 ...
    + 0.5 * xWeight_S(2) * (1 + cos(x(2)))^2 ...
    + 0.5 * xWeight_S(3) * (x(3) - xRef(3))^2 ...
    + 0.5 * xWeight_S(4) * (x(4) - xRef(4))^2 ...
    + 0.5 * uWeight_S * u^2 ...
    + 0.5 * lambdaWeight * lambda^2;

xWeight_T = [1; 1; 1; 1];
L_T = 0.5 * xWeight_T(1) * (x(1) - xRef(1))^2 ...
    + 0.5 * xWeight_T(2) * (1 + cos(x(2)))^2 ...
    + 0.5 * xWeight_T(3) * (x(3) - xRef(3))^2 ...
    + 0.5 * xWeight_T(4) * (x(4) - xRef(4))^2;

% DVI
mass = [1; 0.1];
linkLength = 1;
g = 9.8;
M = [mass(1) + mass(2),                 mass(2) * linkLength * cos(x(2));...
     mass(2) * linkLength * cos(x(2)),  mass(2) * linkLength^2];
C_Mat = [0,   -mass(2) * linkLength * x(4) * sin(x(2));...
     0,   0]; 
G_Mat = [0;...
     -mass(2) * g * linkLength * sin(x(2))];
Bu = [u(1);...
      0]; 
LAMBDA = [lambda(1);...
     0]; % friction bewteen cart and ground  
H = G_Mat + Bu + LAMBDA - C_Mat * [x(3); x(4)];

f = [x(3:4);...
    inv(M)*H];% xDot = f(x, u, lambda)
VISetType = 'box_constraint'; 
bl = -2;
bu = 2;
F = x(3);

% inequality constraint G >= 0
G = [xMax - x;...
    x - xMin;...
    uMax - u;...
    u - uMin];
% equality constraint C = 0
C = SX(0,1);

%% create OCPEC instant
OCPEC = OCPEC_Formulation(...
    TimeHorizon, nStages, timeStep,...
    x0, xRef,...
    x, u, lambda,...
    L_T, L_S,...
    f, g, F, VISetType, bl, bu,...
    G, C);
end

