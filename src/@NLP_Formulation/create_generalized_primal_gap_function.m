function phi_func = create_generalized_primal_gap_function(self, OCPEC)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);
eta = SX.sym('eta', OCPEC.Dim.lambda, 1); % auxiliary variable for VI function F in gap function
eps = self.gap_func_smooth_param;
%%
switch OCPEC.VISetType
    case 'box_constraint'
        % omega has an explicit expression: projection of the stationary point into [bl, bu]
        switch self.strongly_convex_func
            case 'quadratic'
                stationary_point = lambda - eta; 
            case 'general'
                stationary_point = []; % To Be Done                
        end
        % To be done: smooth omega = mid(bl, bu, stationary_point) by CHKS function in elementent-wise
        omega = min(max(OCPEC.bl, stationary_point), OCPEC.bu); 
        phi = self.d_func(lambda) - self.d_func(omega) + (eta' - self.d_grad(lambda)) * (lambda - omega);
        phi_func = Function('phi_func', {lambda, eta}, {phi}, {'lambda', 'eta'}, {'phi'});

    case 'nonnegative_orthant'
        % omega has an explicit expression: projection of the stationary point into [0, inf]
        switch self.strongly_convex_func
            case 'quadratic'
                stationary_point = lambda - eta;
            case 'general'
                stationary_point = []; % To Be Done
        end
        % smoothing omega = max(0, stationary_point) by CHKS function
        omega = 0.5*(sqrt(stationary_point.^2 + 4 * eps^2) + stationary_point); 
        phi = self.d_func(lambda) - self.d_func(omega) + (eta' - self.d_grad(lambda)) * (lambda - omega);
        phi_func = Function('phi_func', {lambda, eta}, {phi}, {'lambda', 'eta'}, {'phi'});

    case 'finitely_representable'
        % omega does not have an explicit expression
        omega = SX.sym('omega', OCPEC.Dim.lambda, 1);
        phi = self.d_func(lambda) - self.d_func(omega) + (eta' - self.d_grad(lambda)) * (lambda - omega);
        phi_func = Function('phi_func', {lambda, eta, omega}, {phi}, {'lambda', 'eta', 'omega'}, {'phi'});

    case 'polyhedral'
        % omega does not have an explicit expression
        omega = SX.sym('omega', OCPEC.Dim.lambda, 1);
        phi = self.d_func(lambda) - self.d_func(omega) + (eta' - self.d_grad(lambda)) * (lambda - omega);
        phi_func = Function('phi_func', {lambda, eta, omega}, {phi}, {'lambda', 'eta', 'omega'}, {'phi'});
        
end

end