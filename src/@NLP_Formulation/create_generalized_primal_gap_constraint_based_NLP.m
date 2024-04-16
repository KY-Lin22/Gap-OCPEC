function nlp = create_generalized_primal_gap_constraint_based_NLP(self, OCPEC)
% formulate a relax generalized primal gap constraint based NLP
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
%            z = [z_1;...z_n;...z_N] and z_n = [x_n; u_n; lambda_n; eta_n; w_n] 
%            with x_n:      system state
%                 u_n:      control                 
%                 lambda_n: algebraic variable   
%                 eta_n:    auxiliary variable for VI function F in gap function
%                 w_n:      auxiliary variable for generalized primal gap function 
%        (2) p: collects all the NLP problem parameters p = [s, mu]
%            with s:  nonnegative relax parameter for gap constraints
%                 mu: nonnegative penalty parameter for penalty function
%        (3) J: cost function J = J_ocp + J_penalty
%            J_ocp = sum(J_stage_n)*dt + J_terminal 
%            J_penalty = mu*sum(J_penalty_n)*dt
%            with J_stage_n:   stage cost defined in OCPEC
%                 J_terminal:  terminal cost defined in OCPEC
%                 J_penalty_n: penalty cost for auxiliary variable w_n
%        (4) h: equality constraint arranged in a stagewise manner
%            h = [h_1;...h_n;...h_N] and h_n = [f_n; F_n - eta_n; phi_n - w_n; C_n]     
%            with f_n:   discretized state equation f defined in OCPEC
%                 F_n:   VI function
%                 phi_n: generalized primal gap function
%                 C_n:   path equality constraints C defined in OCPEC 
%        (5) c: inequality constraint arranged in a stagewise manner
%            c = [c_1;...c_n;...c_N] and c_n = [g_n; s - w_n; G_n] 
%            with g_n: inequality constraints g formulating VI set K  
%                 G_n: path inequality constraints G defined in OCPEC
% output: nlp is a structure with fields:
%         z: variable
%         p: parameter
%         J, h, c: cost and constraint function,  
%         omega: variable used when it can not be explicitly represented by lambda and eta (To Be Done)
%         J_ocp, J_penalty: cost term
%         Dim: problem size

import casadi.*

%% initialize NLP variable (stagewise, capital)
% define auxiliary variable dimension
eta_Dim = OCPEC.Dim.lambda;
w_Dim = 1;
% initialize problem variable 
X = SX.sym('X', OCPEC.Dim.x, OCPEC.nStages); 
XPrev = [OCPEC.x0, X(:, 1 : end - 1)];
U = SX.sym('U', OCPEC.Dim.u, OCPEC.nStages);
LAMBDA = SX.sym('LAMBDA', OCPEC.Dim.lambda, OCPEC.nStages);
ETA = SX.sym('ETA', eta_Dim, OCPEC.nStages); 
W = SX.sym('W', w_Dim, OCPEC.nStages); 
% initialize problem parameter
s = SX.sym('s', 1, 1);
mu = SX.sym('mu', 1, 1);

%% mapping function object
% stage cost
L_S_map = OCPEC.FuncObj.L_S.map(OCPEC.nStages);
% penalty cost
penalty_func = self.create_penalty_func();
penalty_func_map = penalty_func.map(OCPEC.nStages);
% DVI
f_map = OCPEC.FuncObj.f.map(OCPEC.nStages);
F_map = OCPEC.FuncObj.F.map(OCPEC.nStages);
g_map = OCPEC.FuncObj.g.map(OCPEC.nStages);
% path inequality constraints G and equality constraint C defined in OCPEC
G_map = OCPEC.FuncObj.G.map(OCPEC.nStages);
C_map = OCPEC.FuncObj.C.map(OCPEC.nStages);
% phi_c to construct generalized primal gap function phi
phi_c_func = self.create_weighting_generalized_primal_gap_function(OCPEC);
phi_c_func_map = phi_c_func.map(OCPEC.nStages);

%% formulate NLP function (stagewise)
% cost
L_S_stage = L_S_map(X, U, LAMBDA);
L_T = OCPEC.FuncObj.L_T(X(:, end));
penalty_func_stage = penalty_func_map(W);
% DVI
f_stage = f_map(X, U, LAMBDA);
F_stage = F_map(X, U, LAMBDA);
g_stage = g_map(LAMBDA);
% path inequality constraint G and equality constraint C
G_stage = G_map(X, U);
C_stage = C_map(X, U);
% gap function
primal_gap_param_c = 1;
if strcmp(OCPEC.VISetType, 'box_constraint') || strcmp(OCPEC.VISetType, 'nonnegative_orthant')
    phi_stage = phi_c_func_map(LAMBDA, ETA, primal_gap_param_c);
elseif strcmp(OCPEC.VISetType, 'finitely_representable') || strcmp(OCPEC.VISetType, 'polyhedral')
    % TODO, phi is a function of lambda, eta, and omega_c
end

%% reshape NLP variable and function (column, lowercase)
% variable
Z = [X; U; LAMBDA; ETA; W];
z = reshape(Z, (OCPEC.Dim.x + OCPEC.Dim.u + OCPEC.Dim.lambda + eta_Dim + w_Dim) * OCPEC.nStages, 1);
Dim.z_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.u, OCPEC.Dim.lambda, eta_Dim, w_Dim]); 
Dim.z = size(z, 1);

% problem parameter
p = [s; mu];

% cost function
J_ocp = sum(L_S_stage) * OCPEC.timeStep + L_T;
J_penalty = mu * sum(penalty_func_stage) * OCPEC.timeStep;
J = J_ocp + J_penalty;

% equality constraint h = 0
h_stage = [...
    XPrev - X + f_stage * OCPEC.timeStep;...
    F_stage - ETA;...
    phi_stage - W;...
    C_stage];
h = reshape(h_stage, (OCPEC.Dim.x + eta_Dim + w_Dim + OCPEC.Dim.C) * OCPEC.nStages, 1);
Dim.h_Node = cumsum([OCPEC.Dim.x, eta_Dim, w_Dim, OCPEC.Dim.C]);
Dim.h = size(h, 1);

% inequality constraint c >= 0
c_stage = [...
    g_stage;...
    s - W;...
    G_stage];
c = reshape(c_stage, (OCPEC.Dim.g + w_Dim + OCPEC.Dim.G) * OCPEC.nStages, 1);
Dim.c_Node = cumsum([OCPEC.Dim.g, w_Dim, OCPEC.Dim.G]);
Dim.c = size(c, 1);

%% create output struct
nlp = struct('z', z, 'p', p,...
    'J', J, 'h', h, 'c', c,...
    'J_ocp', J_ocp, 'J_penalty', J_penalty,...
    'Dim', Dim);

end