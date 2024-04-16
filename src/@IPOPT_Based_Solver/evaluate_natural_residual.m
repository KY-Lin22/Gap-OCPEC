function natRes = evaluate_natural_residual(self, z_Opt)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*

%% extract solution
Z_Opt = reshape(z_Opt, self.NLP.Dim.z_Node(end), self.OCPEC.nStages);

X_Opt = Z_Opt(1 : self.NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(self.NLP.Dim.z_Node(1) + 1 : self.NLP.Dim.z_Node(2), :);
LAMBDA_Opt = Z_Opt(self.NLP.Dim.z_Node(2) + 1 : self.NLP.Dim.z_Node(3), :);

%% problem data for Euclidean projector
F_FuncObj_map = self.OCPEC.FuncObj.F.map(self.OCPEC.nStages);
F_Opt = full(F_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));
w_Opt = LAMBDA_Opt - F_Opt; 

%% construct qp solver to evaluate Euclidean projector
% variable (projector output) and parameter (projector input)
y = SX.sym('y', self.OCPEC.Dim.lambda, 1);
w = SX.sym('w', self.OCPEC.Dim.lambda, 1);
% cost function and constraint
J = 0.5 * y'*diag(ones(self.OCPEC.Dim.lambda, 1))*y - w'*y;
g = self.OCPEC.FuncObj.g(y); % y \in K := {y: g(y) >= 0}
lbg = zeros(self.OCPEC.Dim.g, 1);
ubg = inf*ones(self.OCPEC.Dim.g, 1);
% problem struct
Prob = struct('x', y, 'f', J, 'g', g, 'p', w);
% option
Option = struct(...
    'printLevel', 'none',...% 'none', 'low', 'medium', 'high' (see qpoases manual sec 5.2)
    'hessian_type', 'posdef',...% 'unknown', 'posdef', 'semidef', 'indef', 'zero', 'identity' (see qpoases manual sec 4.4, 4.5)
    'error_on_fail', false);
solver_singleStage = qpsol('EuclideanProjector', 'qpoases', Prob, Option);
% solver
EuclideanProjector = solver_singleStage.map(self.OCPEC.nStages);

%% evaluate Euclidean projector
solution = EuclideanProjector(...
    'lbg', repmat(lbg, 1, self.OCPEC.nStages),...
    'ubg', repmat(ubg, 1, self.OCPEC.nStages),...
    'p', w_Opt);
proj_w_Opt = full(solution.x);

%% evaluate natural residual
natRes = max(max(abs(LAMBDA_Opt - proj_w_Opt)));

end

