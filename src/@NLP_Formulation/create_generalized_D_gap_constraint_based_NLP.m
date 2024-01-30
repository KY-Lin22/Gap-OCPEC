function nlp = create_generalized_D_gap_constraint_based_NLP(self, OCPEC)
% formulate a relax generalized D gap constraint based NLP
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
%            z = [z_1;...z_n;...z_N] and z_n = [x_n; u_n; lambda_n; eta_n; w_n; v_n] 
%            with x_n:      system state
%                 u_n:      control                 
%                 lambda_n: algebraic variable   
%                 eta_n:    auxiliary variable for VI function F in gap function
%                 w_n:      auxiliary variable for phi_a_n in generalized D gap function phi_ab = phi_a - phi_b 
%                 v_n:      auxiliary variable for phi_b_n in generalized D gap function phi_ab = phi_a - phi_b 
%        (2) p: collects all the NLP problem parameters p = [s, mu]
%            with s:  nonnegative relax parameter for gap constraints
%                 mu: nonnegative penalty parameter for penalty function
%        (3) J: cost function J = J_ocp + J_penalty
%            J_ocp = sum(J_stage_n) + J_terminal 
%            J_penalty = sum(J_penalty_n)
%            with J_stage_n:   stage cost defined in OCPEC
%                 J_terminal:  terminal cost defined in OCPEC
%                 J_penalty_n: penalty cost for auxiliary variable w_n and v_n
%        (4) h: equality constraint arranged in a stagewise manner
%            h = [h_1;...h_n;...h_N] and h_n = [f_n; F_n - eta_n; phi_a_n - w_n; phi_b_n - v_n; C_n]     
%            with f_n:     discretized state equation f defined in OCPEC
%                 F_n:     VI function
%                 phi_a_n: first part of generalized D gap function phi_ab = phi_a - phi_b
%                 phi_b_n: second part of generalized D gap function phi_ab = phi_a - phi_b
%                 C_n:     path equality constraints C defined in OCPEC 
%        (5) c: inequality constraint arranged in a stagewise manner
%            c = [c_1;...c_n;...c_N] and c_n = [s - w_n + v_n; G_n] 
%            with G_n: path inequality constraints G defined in OCPEC
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
v_Dim = w_Dim;
% initialize problem variable 
X = SX.sym('X', OCPEC.Dim.x, OCPEC.nStages); 
XPrev = [OCPEC.x0, X(:, 1 : end - 1)];
U = SX.sym('U', OCPEC.Dim.u, OCPEC.nStages);
LAMBDA = SX.sym('LAMBDA', OCPEC.Dim.lambda, OCPEC.nStages);
ETA = SX.sym('ETA', eta_Dim, OCPEC.nStages); % auxiliary variable for VI function F
W = SX.sym('W', w_Dim, OCPEC.nStages); % auxiliary variable for phi_a in generalized D gap function
V = SX.sym('V', v_Dim, OCPEC.nStages); % auxiliary variable for phi_b in generalized D gap function
% initialize problem parameter
s = SX.sym('s', 1, 1);
mu = SX.sym('mu', 1, 1);

%% mapping function object
% stage cost
L_S_map = OCPEC.FuncObj.L_S.map(OCPEC.nStages);
% penalty cost
penalty_map = self.penalty_func.map(OCPEC.nStages);
% discretization state equation
f_map = self.discre_state_equ_func.map(OCPEC.nStages);
% VI function F and set K inequality constraint g
F_map = OCPEC.FuncObj.F.map(OCPEC.nStages);
% path inequality constraints G and equality constraint C defined in OCPEC
G_map = OCPEC.FuncObj.G.map(OCPEC.nStages);
C_map = OCPEC.FuncObj.C.map(OCPEC.nStages);
% phi_c to construct generalized D gap function phi_ab = phi_a - phi_b
phi_c_map = self.phi_c_func.map(OCPEC.nStages);

%% formulate NLP function (stagewise)
% cost
L_S_stage = L_S_map(X, U, LAMBDA);
L_T = OCPEC.FuncObj.L_T(X(:, end));
penalty_stage = penalty_map(W, mu);
% discretized state equation
f_stage = f_map(XPrev, X, U, LAMBDA);
% VI function and set
F_stage = F_map(X, U, LAMBDA);
% path inequality constraint G and equality constraint C
G_stage = G_map(X, U);
C_stage = C_map(X, U);
% gap function
D_gap_param_a = self.D_gap_param_a;
D_gap_param_b = self.D_gap_param_b;
if strcmp(OCPEC.VISetType, 'box_constraint') || strcmp(OCPEC.VISetType, 'nonnegative_orthant')
    phi_a_stage = phi_c_map(LAMBDA, ETA, D_gap_param_a);
    phi_b_stage = phi_c_map(LAMBDA, ETA, D_gap_param_b);
elseif strcmp(OCPEC.VISetType, 'finitely_representable') || strcmp(OCPEC.VISetType, 'polyhedral')
    % TODO, phi_a is a function of lambda, eta, and omega_a
    % TODO, phi_b is a function of lambda, eta, and omega_b
end

%% reshape NLP variable and function (column, lowercase)
% variable
Z = [X; U; LAMBDA; ETA; W; V];
z = reshape(Z, (OCPEC.Dim.x + OCPEC.Dim.u + OCPEC.Dim.lambda + eta_Dim + w_Dim + v_Dim) * OCPEC.nStages, 1);
Dim.z_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.u, OCPEC.Dim.lambda, eta_Dim, w_Dim, v_Dim]); 
Dim.z = size(z, 1);

% problem parameter
p = [s; mu];

% cost function
J_ocp = OCPEC.timeStep*sum(L_S_stage) + L_T;
J_penalty = OCPEC.timeStep*sum(penalty_stage);
J = J_ocp + J_penalty;

% equality constraint h = 0
h_stage = [f_stage; F_stage - ETA; phi_a_stage - W; phi_b_stage - V; C_stage];
h = reshape(h_stage, (OCPEC.Dim.x + eta_Dim + w_Dim + v_Dim + OCPEC.Dim.C) * OCPEC.nStages, 1);
Dim.h_Node = cumsum([OCPEC.Dim.x, eta_Dim, w_Dim, v_Dim, OCPEC.Dim.C]);
Dim.h = size(h, 1);

% inequality constraint c >= 0
c_stage = [s - W + V; G_stage];
c = reshape(c_stage, (w_Dim + OCPEC.Dim.G) * OCPEC.nStages, 1);
Dim.c_Node = cumsum([w_Dim, OCPEC.Dim.G]);
Dim.c = size(c, 1);

%% create output struct
nlp = struct('z', z, 'p', p,...
    'J', J, 'h', h, 'c', c,...
    'J_ocp', J_ocp, 'J_penalty', J_penalty,...
    'Dim', Dim);
end