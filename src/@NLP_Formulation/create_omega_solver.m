function omega_solver = create_omega_solver(self, OCPEC, param_c)
% create omega_solver to evaluate omega_c in the weighting generalized primal gap function
% omega_c is the solution to the strongly concave maximization problem with the form of:
%  max  f_c(omega_c),
%  s.t. lbg <= g(omega_c) <= ubg,
% where f_c is parameterized by a given positive constant c
%
% omega_solver is a function object either being a optimization solver parameterized by lambda and omega 
%                                             or a symbolic function with input lambda and omega
import casadi.*

% initialize problem parameter (lambda, eta)
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);
eta = SX.sym('eta', OCPEC.Dim.lambda, 1);

%% create omega solver based on the VI set type
switch OCPEC.VISetType
    case 'finitely_representable'
        %% the strictly concave maximization problem is a general convex problem
        % utilize a convex optimization solver cvx or yalmz(TODO)

    case 'polyhedral'
        %% the strictly concave maximization problem is a QP
        % utilize a convex active-set QP solver qpoases (density pattern)
        
        % problem variable (omega_c)
        omega_c = SX.sym('omega_c', OCPEC.Dim.lambda, 1);
        % problem parameter vector (lambda, omega)
        p = [lambda; eta];
        % strongly convex function d and its derivative
        [d_func, d_grad, ~] = self.create_strongly_convex_func(OCPEC);
        % problem cost function (max)        
        f_c = - param_c * d_func(omega_c) - (eta' - param_c * d_grad(lambda)) * omega_c;
        % problem constraint (VI set K): omega \in K := {omega | g(omega) >= 0}
        g = OCPEC.FuncObj.g(omega_c);
        % create problem: transfer strictly concave maximization problem 
        % into strictgly convex minimization problem to utilize QP solver
        StrictlyConvexMinProb = struct(...
            'x', omega_c,...
            'f', -f_c,...% note: minimize cost!
            'g', g,...
            'p', p);
        % create option
        qpoases_Option = struct(...
            'printLevel', 'none',...% 'none', 'low', 'medium', 'high' (see qpoases manual sec 5.2)
            'hessian_type', 'posdef',...% 'unknown', 'posdef', 'semidef', 'indef', 'zero', 'identity' (see qpoases manual sec 4.4, 4.5)
            'error_on_fail', false); 
        % create function object: QP solver 
        omega_solver = qpsol('omega_solver', 'qpoases', StrictlyConvexMinProb, qpoases_Option);        

    case 'box_constraint'
        %% strictly concave maximization problem has explicit solution 
        % omega_c has the explicit expression: projection of the stationary point into [bl, bu]      
        switch self.strongly_convex_func
            case 'quadratic'
                stationary_point_c = lambda - (1/param_c) * eta;           
        end
        % smoothing omega_c' expression: omega_c = mid(bl, bu, stationary_point_c)
        epsilon = self.gap_func_smooth_param;
        if epsilon == 0
            % do not smoothing 
            omega_c = max(OCPEC.bl, min(stationary_point_c, OCPEC.bu));
        else
            % smoothing by CHKS function: refer to example 2 in ''A new look at smoothing Newton methods for NCPs and BVI'', Mathematical Programming, 2000
            omega_c = 0.5*(OCPEC.bl + sqrt((OCPEC.bl - stationary_point_c).^2 + 4*epsilon^2)) ...
                + 0.5*(OCPEC.bu - sqrt((OCPEC.bu - stationary_point_c).^2 + 4*epsilon^2));
        end
        % create a single stage solver 
        omega_solver = Function('omega_solver', {lambda, eta}, {omega_c}, {'lambda', 'eta'}, {'omega_c'});

    case 'nonnegative_orthant'
        %% strictly concave maximization problem has explicit solution
        % omega has the explicit expression: projection of the stationary point into [0, inf]
        switch self.strongly_convex_func
            case 'quadratic'
                stationary_point_c = lambda - (1/param_c) * eta;
        end
        % smoothing omega_c's expression: omega_c = max(0, stationary_point_c)
        epsilon = self.gap_func_smooth_param;
        if epsilon == 0
            % do not smoothing
            omega_c = max(zeros(OCPEC.Dim.lambda, 1), stationary_point_c);
        else
            % smoothing by CHKS function
            omega_c = 0.5*(sqrt(stationary_point_c.^2 + 4 * epsilon^2) + stationary_point_c);
        end
        % create a single stage solver 
        omega_solver = Function('omega_solver', {lambda, eta}, {omega_c}, {'lambda', 'eta'}, {'omega_c'});

end

end