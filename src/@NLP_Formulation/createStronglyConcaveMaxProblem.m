function OmegaEvalProb = createStronglyConcaveMaxProblem(self, OCPEC)
% create a strongly concave maximization problem to evaluate omega in gap function
% (1) generalized primal gap constraint based: one problem with the form of:
%  max  f(omega),
%  s.t. lbg <= g(omega) <= ubg,
% (2) generalized D gap constraint based: two problems respectively with the form of:
%  max  f_a(omega),
%  s.t. lbg <= g(omega) <= ubg,
%  max  f_b(omega),
%  s.t. lbg <= g(omega) <= ubg,

import casadi.*
%%
% initialize problem variable (omega)
omega = SX.sym('omega', OCPEC.Dim.lambda, 1);

% initialize problem parameter (lambda, eta)
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);
eta = SX.sym('eta', OCPEC.Dim.lambda, 1); % auxiliary variable for VI function F(x, u, lambda)

% constraint (VI set K): omega \in K := {omega | g(omega) >= 0}
g = OCPEC.FuncObj.g(omega); 
lbg = zeros(OCPEC.Dim.g, 1);
ubg = inf*ones(OCPEC.Dim.g, 1);

% cost function (max) 
switch self.relaxProbType
    case 'generalized_primal_gap_constraint_based'       
        f = - self.d_func(omega) - (eta' - self.d_grad(lambda)) * omega;
    case 'generalized_D_gap_constraint_based' 
        a = self.D_gap_param.a;
        b = self.D_gap_param.b;
        f.a = - a * self.d_func(omega) - (eta' - a * self.d_grad(lambda)) * omega;
        f.b = - b * self.d_func(omega) - (eta' - b * self.d_grad(lambda)) * omega;
end

%% create output struct
% problem variable, parameter, constraint function, and constraint bounds
OmegaEvalProb = struct(...
    'x', omega,...
    'p', [lambda; eta],...
    'g', g,...
    'lbg', lbg,...
    'ubg', ubg);
% problem cost function
switch self.relaxProbType
    case 'generalized_primal_gap_constraint_based'        
        OmegaEvalProb.f = f;
    case 'generalized_D_gap_constraint_based'
        OmegaEvalProb.f.a = f.a;
        OmegaEvalProb.f.b = f.b;
end

end

