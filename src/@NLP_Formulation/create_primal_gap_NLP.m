function nlp = create_primal_gap_NLP(self, OCPEC)
% formulate a relax generalized primal gap constraint based NLP (omega_c is a implicit variable)
%
% OCPEC has the form:
%  min  L_T(x) + int_0^T L_S(x, u, lambda) dt,
%  s.t. Dot{x} = f(x, u, lambda)
%       lambda \in SOL(K, F(x, u, lambda))
%       K = {lambda | g(lambda) >= 0}
%       G(x, u) >= 0,
%       C(x, u) = 0,
%
% NLP has the form:
%  min  J(z, p),
%  s.t. h(z, p) = 0,
%       c(z, p) >= 0
% where: (1) z: collects all the variables to be optimized and arranged in a stagewise manner
%            z = [z_1;...z_n;...z_N] and z_n = [x_n; u_n; lambda_n; eta_n] 
%            with x_n:      system state
%                 u_n:      control                 
%                 lambda_n: algebraic variable   
%                 eta_n:    auxiliary variable for VI function F in gap function
%        (2) p: collects all the NLP problem parameters p = s
%            with s:  nonnegative relax parameter for gap constraints
%        (3) J: cost function J = sum(J_stage_n)*dt + J_terminal
%            with J_stage_n:   stage cost defined in OCPEC
%                 J_terminal:  terminal cost defined in OCPEC
%        (4) h: equality constraint arranged in a stagewise manner
%            h = [h_1;...h_n;...h_N] and h_n = [f_n; F_n - eta_n; C_n]     
%            with f_n:   discretized state equation f defined in OCPEC
%                 F_n:   VI function
%                 C_n:   path equality constraints C defined in OCPEC 
%        (5) c: inequality constraint arranged in a stagewise manner
%            c = [c_1;...c_n;...c_N] and c_n = [g_n; s - phi_n; G_n] 
%            with g_n: inequality constraints g formulating VI set K  
%                 phi_n: generalized primal gap function (function of lambda, eta)
%                 G_n: path inequality constraints G defined in OCPEC
% output: nlp is a structure with fields:
%         z: variable
%         p: parameter
%         J, h, c: cost and constraint function,  
%         Dim: problem size

import casadi.*

%% initialize NLP variable (stagewise, capital)
% initialize problem variable 
X = SX.sym('X', OCPEC.Dim.x, OCPEC.nStages); 
XPrev = [OCPEC.x0, X(:, 1 : end - 1)];
U = SX.sym('U', OCPEC.Dim.u, OCPEC.nStages);
LAMBDA = SX.sym('LAMBDA', OCPEC.Dim.lambda, OCPEC.nStages);
ETA = SX.sym('ETA', OCPEC.Dim.lambda, OCPEC.nStages); 
% initialize problem parameter
s = SX.sym('s', 1, 1);

%% mapping function object
% stage cost
L_S_map = OCPEC.FuncObj.L_S.map(OCPEC.nStages);
% DVI
f_map = OCPEC.FuncObj.f.map(OCPEC.nStages);
F_map = OCPEC.FuncObj.F.map(OCPEC.nStages);
g_map = OCPEC.FuncObj.g.map(OCPEC.nStages);
% path inequality constraints G and equality constraint C defined in OCPEC
G_map = OCPEC.FuncObj.G.map(OCPEC.nStages);
C_map = OCPEC.FuncObj.C.map(OCPEC.nStages);
% generalized primal gap function phi_c (function of lambda and eta)
param_c = 1;
phi_c_func = self.create_phi_c_func(OCPEC, param_c);
phi_c_func_map = phi_c_func.map(OCPEC.nStages);

%% formulate NLP function (stagewise)
% cost
L_S_stage = L_S_map(X, U, LAMBDA);
L_T = OCPEC.FuncObj.L_T(X(:, end));
% DVI
f_stage = f_map(X, U, LAMBDA);
F_stage = F_map(X, U, LAMBDA);
g_stage = g_map(LAMBDA);
% path inequality constraint G and equality constraint C
G_stage = G_map(X, U);
C_stage = C_map(X, U);
% gap function
phi_stage = phi_c_func_map(LAMBDA, ETA);

%% reshape NLP variable and function (column, lowercase)
% variable
Z = [X; U; LAMBDA; ETA];
z = reshape(Z, (OCPEC.Dim.x + OCPEC.Dim.u + OCPEC.Dim.lambda + OCPEC.Dim.lambda) * OCPEC.nStages, 1);
Dim.z_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.u, OCPEC.Dim.lambda, OCPEC.Dim.lambda]); 
Dim.z = size(z, 1);

% problem parameter
p = s;
Dim.p = size(p, 1);

% cost function
J = sum(L_S_stage) * OCPEC.timeStep + L_T;

% equality constraint h = 0
h_stage = [...
    XPrev - X + f_stage * OCPEC.timeStep;...
    F_stage - ETA;...
    C_stage];
h = reshape(h_stage, (OCPEC.Dim.x + OCPEC.Dim.lambda + OCPEC.Dim.C) * OCPEC.nStages, 1);
Dim.h_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.lambda, OCPEC.Dim.C]);
Dim.h = size(h, 1);

% inequality constraint c >= 0
c_stage = [...
    g_stage;...
    s - phi_stage;...
    G_stage];
c = reshape(c_stage, (OCPEC.Dim.g + 1 + OCPEC.Dim.G) * OCPEC.nStages, 1);
Dim.c_Node = cumsum([OCPEC.Dim.g, 1, OCPEC.Dim.G]);
Dim.c = size(c, 1);

%% create output struct
nlp = struct('z', z, 'p', p,...
    'J', J, 'h', h, 'c', c,...
    'Dim', Dim);

end