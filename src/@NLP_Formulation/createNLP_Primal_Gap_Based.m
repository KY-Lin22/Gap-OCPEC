function nlp = createNLP_Primal_Gap_Based(self, OCPEC)
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
%                 w_n:      scalar auxiliary variable to transfer gap function inequality into equality
%        (2) p: collects all the NLP problem parameters
%            p = s
%            with s:  nonnegative relax parameter for gap constraints
%        (3) J: cost function J = sum(J_n) and J_n = J_stage_n
%            with J_stage_n:   stage cost defined in OCPEC
%        (4) h: equality constraint arranged in a stagewise manner
%            h = [h_1;...h_n;...h_N] and h_n = [f_n; F_n - eta_n; phi_n - w_n; C_n]     
%            with f_n:   discretized state equation f defined in OCPEC
%                 F_n:   VI function
%                 phi_n: primal gap constraint (variable: lambda_n, eta_n)
%                 C_n:   path equality constraints C defined in OCPEC   
%        (5) c: inequality constraint arranged in a stagewise manner
%            c = [c_1;...c_n;...c_N] and c_n = [g_n; s - w_n; G_n] 
%            with g_n: inequality constraints g formulating VI set K  
%                 G_n: path inequality constraints G defined in OCPEC
%
% output: nlp is a structure with fields:
%         z, omega: variable
%         p: parameter
%         J, h, c: cost and constraint function,                                 
%         Dim: problem size
import casadi.*
%
%% initialize NLP variable and function (stagewise, capital)
% define auxiliary variable dimension
eta_Dim = OCPEC.Dim.lambda;
w_Dim = 1;

% initialize problem variable 
X = SX.sym('X', OCPEC.Dim.x, OCPEC.nStages); 
U = SX.sym('U', OCPEC.Dim.u, OCPEC.nStages);
LAMBDA = SX.sym('LAMBDA', OCPEC.Dim.lambda, OCPEC.nStages);
ETA = SX.sym('ETA', eta_Dim, OCPEC.nStages);
W = SX.sym('W', w_Dim, OCPEC.nStages); % auxiliary variable for gap function

Z = [X; U; LAMBDA; ETA; W]; % NLP variable in stagewise manner

% initialize problem parameter
s = SX.sym('s', 1, 1);

% initialize problem cost and constraint function
J = 0;
H = SX(OCPEC.Dim.x + eta_Dim + w_Dim + OCPEC.Dim.C, OCPEC.nStages); 
C = SX(OCPEC.Dim.g + w_Dim + OCPEC.Dim.G, OCPEC.nStages);

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
    eta_n = ETA(:, n);   
    w_n = W(:, n);
    
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
    
    % VI function
    F_n = OCPEC.FuncObj.F(x_n, u_n, lambda_n);
    
    % path inequality constraints G defined in OCPEC   
    G_n = OCPEC.FuncObj.G(x_n, u_n);  
    
    % path equality constraints C defined in OCPEC
    C_n = OCPEC.FuncObj.C(x_n, u_n);
    
    % primal gap function phi (function of lambda, eta, but it needs an intermedia variable omega which is a function of lambda, eta)   
    switch OCPEC.VISetType
        case 'box_constraint'
            %
            switch self.stronglyConvexFuncType
                case 'quadratic'
                    % omega has an explicit expression: projection of the stationary point into [bl, bu]
                    omega_n = min(max(OCPEC.bl, lambda_n - eta_n), OCPEC.bu);  
                    phi_n = self.d_func(lambda_n) - self.d_func(omega_n) ...
                        + (eta_n' - self.d_grad(lambda_n)) * (lambda_n - omega_n);
                case 'general'
                    % To Be Done 
                    % omega also has an explicit expression: projection of the stationary point into [bl, bu]
                    omega_n = []; 
                    phi_n = self.d_func(lambda_n) - self.d_func(omega_n) ...
                        + (eta_n' - self.d_grad(lambda_n)) * (lambda_n - omega_n);
            end
        case 'nonnegative_orthant'
            %
            switch self.stronglyConvexFuncType
                case 'quadratic'
                    % explicit expression of phi (ref: P940, Finite dim VI and CP, 2003, snet c = 1 in the third equation)
                    max_operator = max(zeros(OCPEC.Dim.lambda, 1), eta_n - lambda_n);
                    phi_n = 1/2 * ((eta_n' * eta_n) - (max_operator' * max_operator));                    
                case 'general'
                    % To Be Done 
                    % omega has a explicit expression: projection of the stationary point into [0, inf]
                    omega_n = []; 
                    phi_n = self.d_func(lambda_n) - self.d_func(omega_n) ...
                        + (eta_n' - self.d_grad(lambda_n)) * (lambda_n - omega_n);
            end

        case 'finitely_representable'
            % To Be Done (use symbolic variable omega like the SGCL solver)
            switch self.stronglyConvexFuncType
                case 'quadratic'
                case 'general'
            end            
        case 'polyhedral'    
            % To Be Done (use symbolic variable omega like the SGCL solver)
            switch self.stronglyConvexFuncType
                case 'quadratic'
                case 'general'
            end            
    end
    
    % summarize NLP constraint
    h_n = [f_n; F_n - eta_n; phi_n - w_n; C_n];
    H(:, n) = h_n;    
    
    c_n = [g_n; s - w_n; G_n];
    C(:, n) = c_n;        
end

%% formulate NLP variable and function (column, lowercase)
% problem variable
z = reshape(Z, (OCPEC.Dim.x + OCPEC.Dim.u + OCPEC.Dim.lambda + eta_Dim + w_Dim) * OCPEC.nStages, 1);
% problem parameter
p = s;
% constraint
h = reshape(H, (OCPEC.Dim.x + eta_Dim + w_Dim + OCPEC.Dim.C) * OCPEC.nStages, 1);
c = reshape(C, (OCPEC.Dim.g + w_Dim + OCPEC.Dim.G) * OCPEC.nStages, 1);
% size and node point of variable and function (NLP)
Dim.z = size(z, 1);
Dim.h = size(h, 1);
Dim.c = size(c, 1);
Dim.z_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.u, OCPEC.Dim.lambda, eta_Dim, w_Dim]); % node point after reshaping z into stagewise Z
Dim.h_Node = cumsum([OCPEC.Dim.x, eta_Dim, w_Dim, OCPEC.Dim.C]); % node point after reshaping h into stagewise H
Dim.c_Node = cumsum([OCPEC.Dim.g, w_Dim, OCPEC.Dim.G]); % node point after reshaping c into stagewise C

%% create output struct
nlp = struct('z', z, 'p', p,...
    'J', J, 'h', h, 'c', c,...
    'Dim', Dim);

end

