function nlp = create_KKT_based_NLP(self, OCPEC)
% formulate a relax KKT based NLP
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
%            z = [z_1;...z_n;...z_N] and z_n = [x_n; u_n; lambda_n; zeta_n, w_n] 
%            with x_n:      system state
%                 u_n:      control                 
%                 lambda_n: algebraic variable   
%                 zeta_n:   dual variable for inequality g formulating VI set K   
%                 w_n:      auxiliary variable to transfer inequality g into equality
%        (2) p: collects all the NLP problem parameters p = [s, mu]
%            with s:  nonnegative relax parameter for gap constraints
%                 mu: nonnegative penalty parameter for penalty function
%        (3) J: cost function J = J_ocp + J_penalty
%            J_ocp = sum(J_stage_n) + J_terminal 
%            J_penalty = sum(J_penalty_n)
%            with J_stage_n:   stage cost defined in OCPEC
%                 J_terminal:  terminal cost defined in OCPEC
%                 J_penalty_n: penalty cost (KKT based NLP sets penalty term as zero)
%        (4) h: equality constraint arranged in a stagewise manner
%            h = [h_1;...h_n;...h_N] and h_n = [f_n; KKT_stationarity_n; g_n - w_n; C_n]     
%            with f_n:                discretized state equation f defined in OCPEC
%                 KKT_stationarity_n: KKT stationarity condition for VI reformulation
%                 g_n:                inequality constraints g formulating VI set K
%                 C_n:                path equality constraints C defined in OCPEC 
%        (5) c: inequality constraint arranged in a stagewise manner
%            c = [c_1;...c_n;...c_N] and c_n = [KKT_complementarity_n; G_n] 
%            with KKT_complementarity_n: KKT complementarity condition for VI reformulation
%                 G_n:                   path inequality constraints G defined in OCPEC
% output: nlp is a structure with fields:
%         z: variable
%         p: parameter
%         J, h, c: cost and constraint function,  
%         J_ocp, J_penalty: cost term
%         Dim: problem size

import casadi.*


%% initialize NLP variable (stagewise, capital)
% define dual and auxiliary variable dimension
zeta_Dim = OCPEC.Dim.g;
w_Dim = OCPEC.Dim.g;
% initialize problem variable 
X = SX.sym('X', OCPEC.Dim.x, OCPEC.nStages); 
XPrev = [OCPEC.x0, X(:, 1 : end - 1)];
U = SX.sym('U', OCPEC.Dim.u, OCPEC.nStages);
LAMBDA = SX.sym('LAMBDA', OCPEC.Dim.lambda, OCPEC.nStages);
ZETA = SX.sym('ETA', zeta_Dim, OCPEC.nStages); % dual variable for g
W = SX.sym('W', w_Dim, OCPEC.nStages); % auxiliary variable for g
% initialize problem parameter
s = SX.sym('s', 1, 1);
mu = SX.sym('mu', 1, 1);

%% mapping function object
% stage cost
L_S_map = OCPEC.FuncObj.L_S.map(OCPEC.nStages);
% discretization state equation
f_map = self.discre_state_equ_func.map(OCPEC.nStages);
% VI set K inequality constraint g>=0
g_map = OCPEC.FuncObj.g.map(OCPEC.nStages);
% path inequality constraints G and equality constraint C defined in OCPEC
G_map = OCPEC.FuncObj.G.map(OCPEC.nStages);
C_map = OCPEC.FuncObj.C.map(OCPEC.nStages);
% KKT reformulation for VI
KKT_stationarity_map = self.KKT_stationarity_func.map(OCPEC.nStages);
KKT_complementarity_map = self.KKT_complementarity_func.map(OCPEC.nStages);

%% formulate NLP function (stagewise, capital)
% cost
L_S_stage = L_S_map(X, U, LAMBDA);
L_T = OCPEC.FuncObj.L_T(X(:, end));
penalty_stage = mu * 0;
% discretized state equation
f_stage = f_map(XPrev, X, U, LAMBDA);
% VI function and set
g_stage = g_map(LAMBDA);
% path inequality constraint G and equality constraint C
G_stage = G_map(X, U);
C_stage = C_map(X, U);
% KKT reformulation for VI
KKT_stationarity_stage = KKT_stationarity_map(X, U, LAMBDA,ZETA);
KKT_complementarity_stage = KKT_complementarity_map(ZETA, W, s);

%% reshape NLP variable and function (column, lowercase)
% variable
Z = [X; U; LAMBDA; ZETA; W];
z = reshape(Z, (OCPEC.Dim.x + OCPEC.Dim.u + OCPEC.Dim.lambda + zeta_Dim + w_Dim) * OCPEC.nStages, 1);
Dim.z_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.u, OCPEC.Dim.lambda, zeta_Dim, w_Dim]); 
Dim.z = size(z, 1);

% problem parameter
p = [s; mu];

% cost function
J_ocp = OCPEC.timeStep*sum(L_S_stage) + L_T;
J_penalty = OCPEC.timeStep*sum(penalty_stage);
J = J_ocp + J_penalty;

% equality constraint h = 0
h_stage = [f_stage; KKT_stationarity_stage; g_stage - W; C_stage];
h = reshape(h_stage, (OCPEC.Dim.x + OCPEC.Dim.lambda + OCPEC.Dim.g + OCPEC.Dim.C) * OCPEC.nStages, 1);
Dim.h_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.lambda, OCPEC.Dim.g, OCPEC.Dim.C]);
Dim.h = size(h, 1);

% inequality constraint c >= 0
c_stage = [KKT_complementarity_stage; G_stage];
c = reshape(c_stage, (size(KKT_complementarity_stage, 1) + OCPEC.Dim.G) * OCPEC.nStages, 1);
Dim.c_Node = cumsum([size(KKT_complementarity_stage, 1), OCPEC.Dim.G]);
Dim.c = size(c, 1);

%% create output struct
nlp = struct('z', z, 'p', p,...
    'J', J, 'h', h, 'c', c,...
    'J_ocp', J_ocp, 'J_penalty', J_penalty,...
    'Dim', Dim);

end