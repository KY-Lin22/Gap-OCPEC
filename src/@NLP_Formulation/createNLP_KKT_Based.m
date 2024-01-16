function nlp = createNLP_KKT_Based(self, OCPEC)
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
%            z = [z_1;...z_n;...z_N] and z_n = [x_n; u_n; lambda_n; zeta_n] 
%            with x_n:      system state
%                 u_n:      control                 
%                 lambda_n: algebraic variable   
%                 zeta_n:     dual variable for inequality g formulating VI set K   
%        (2) p: collects all the NLP problem parameters
%            p = s
%        (3) J: cost function J = sum(J_n) and J_n = J_stage_n
%            with J_stage_n:   stage cost defined in OCPEC
%        (4) h: equality constraint arranged in a stagewise manner
%            h = [h_1;...h_n;...h_N] and h_n = [f_n; F_n - zeta_n' * glambda; C_n]     
%            with f_n:     discretized state equation f defined in OCPEC
%                 F_n:     VI function
%                 glambda: jacobian of inequality constraint g formulating VI set K
%                 C_n:     path equality constraints C defined in OCPEC  
%        (5) c: inequality constraint arranged in a stagewise manner
%            c = [c_1;...c_n;...c_N] and c_n = [zeta_n; g_n; s - zeta*g_n; G_n] 
%            with g_n: inequality constraints g formulating VI set K  
%                 G_n: path inequality constraints G defined in OCPEC
% output: nlp is a structure with fields:
%         z, omega: variable
%         p: parameter
%         J, h, c: cost and constraint function,                                 
%         Dim: problem size
import casadi.*

%% initialize NLP variable and function (stagewise, capital)
% define auxiliary variable dimension
zeta_Dim = OCPEC.Dim.g;

% initialize problem variable 
X = SX.sym('X', OCPEC.Dim.x, OCPEC.nStages); 
U = SX.sym('U', OCPEC.Dim.u, OCPEC.nStages);
LAMBDA = SX.sym('LAMBDA', OCPEC.Dim.lambda, OCPEC.nStages);
ZETA = SX.sym('ETA', zeta_Dim, OCPEC.nStages);

Z = [X; U; LAMBDA; ZETA]; % NLP variable in stagewise manner

% initialize problem parameter
s = SX.sym('s', 1, 1);

% initialize problem cost and constraint function
J = 0;
H = SX(OCPEC.Dim.x + OCPEC.Dim.lambda + OCPEC.Dim.C, OCPEC.nStages); 
C = SX(3 * OCPEC.Dim.g + OCPEC.Dim.G, OCPEC.nStages);

%% formulate NLP function
for n = 1 : OCPEC.nStages 
    % load variable
    if n == 1
        x_nPrev = OCPEC.x0;
    else
        x_nPrev = X(:, n - 1);
    end
    x_n = X(:, n);
    u_n = U(:, n);
    lambda_n = LAMBDA(:, n);
    zeta_n = ZETA(:, n);    

    % NLP stage cost function
    J_stage_n = OCPEC.timeStep * OCPEC.FuncObj.L_S(x_n, u_n, lambda_n);
    if (n == OCPEC.nStages) && (size(OCPEC.L_T, 1) == 1)
        J_stage_n = J_stage_n + OCPEC.FuncObj.L_T(x_n);
    end    
    % summarize NLP cost function
    J = J + J_stage_n; 
    
    % discretized state equation f defined in OCPEC (implicit Euler method)
    f_n = x_nPrev - x_n + OCPEC.timeStep * OCPEC.FuncObj.f(x_n, u_n, lambda_n); 
    
    % inequality constraint g formulating VI set K
    g_n = OCPEC.FuncObj.g(lambda_n);    
    glambda_n = OCPEC.FuncObj.glambda(lambda_n);
    
    % VI function
    F_n = OCPEC.FuncObj.F(x_n, u_n, lambda_n);
    
    % path inequality constraints G defined in OCPEC   
    G_n = OCPEC.FuncObj.G(x_n, u_n);  
    
    % path equality constraints C defined in OCPEC
    C_n = OCPEC.FuncObj.C(x_n, u_n);
    
    % summarize NLP constraint
    h_n = [f_n; F_n - glambda_n' * zeta_n; C_n];
    H(:, n) = h_n;    
    
    c_n = [zeta_n; g_n; s - zeta_n .* g_n; G_n];
    C(:, n) = c_n;     
end

%% formulate NLP variable and function (column, lowercase)
% problem variable
z = reshape(Z, (OCPEC.Dim.x + OCPEC.Dim.u + OCPEC.Dim.lambda + zeta_Dim) * OCPEC.nStages, 1);
% problem parameter
p = s;
% constraint
h = reshape(H, (OCPEC.Dim.x + OCPEC.Dim.lambda + OCPEC.Dim.C) * OCPEC.nStages, 1);
c = reshape(C, (3 * OCPEC.Dim.g + OCPEC.Dim.G) * OCPEC.nStages, 1);

% size and node point of variable and function (NLP)
Dim.z = size(z, 1);
Dim.h = size(h, 1);
Dim.c = size(c, 1);
Dim.z_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.u, OCPEC.Dim.lambda, zeta_Dim]); % node point after reshaping z into stagewise Z
Dim.h_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.lambda, OCPEC.Dim.C]); % node point after reshaping h into stagewise H
Dim.c_Node = cumsum([OCPEC.Dim.g, OCPEC.Dim.g, OCPEC.Dim.g, OCPEC.Dim.G]); % node point after reshaping c into stagewise C

%% create output struct
nlp = struct('z', z, 'p', p,...
    'J', J, 'h', h, 'c', c,...
    'Dim', Dim);

end
