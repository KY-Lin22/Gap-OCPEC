function phi_func  = create_generalized_D_gap_function(self, OCPEC)
import casadi.*
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);
eta = SX.sym('eta', OCPEC.Dim.lambda, 1); % auxiliary variable for VI function F in gap function
epsilon = self.gap_func_smooth_param;
a = self.D_gap_param_a;
b = self.D_gap_param_b;

%%
switch OCPEC.VISetType
    case 'box_constraint'
        % omega has an explicit expression: projection of the stationary point into [bl, bu]
        switch self.strongly_convex_func
            case 'quadratic'
                stationary_point_a = lambda - (1/a) * eta;
                stationary_point_b = lambda - (1/b) * eta;
            case 'general'
                stationary_point_a = [];% To Be Done
                stationary_point_b = [];
        end
        % smoothing omega = mid(bl, bu, stationary_point) by CHKS function 
        % ref: example 2 in ''A new look at smoothing Newton methods for NCPs and BVI'', Mathematical Programming, 2000
        omega_a = 0.5*(OCPEC.bl + sqrt((OCPEC.bl - stationary_point_a).^2 + 4*epsilon^2)) ...
            + 0.5*(OCPEC.bu - sqrt((OCPEC.bu - stationary_point_a).^2 + 4*epsilon^2));
        omega_b = 0.5*(OCPEC.bl + sqrt((OCPEC.bl - stationary_point_b).^2 + 4*epsilon^2)) ...
            + 0.5*(OCPEC.bu - sqrt((OCPEC.bu - stationary_point_b).^2 + 4*epsilon^2));
        q_a = self.d_func(lambda) - self.d_func(omega_a) + self.d_grad(lambda) * (omega_a - lambda);
        q_b = self.d_func(lambda) - self.d_func(omega_b) + self.d_grad(lambda) * (omega_b - lambda);
        phi_ab = (omega_b - omega_a)' * eta + a * q_a - b * q_b;
        phi_func = Function('phi_func', {lambda, eta}, {phi_ab}, {'lambda', 'eta'}, {'phi_ab'});

    case 'nonnegative_orthant'
        % omega has an explicit expression: projection of the stationary point into [0, inf]
        switch self.strongly_convex_func
            case 'quadratic'
                stationary_point_a = lambda - (1/a) * eta;
                stationary_point_b = lambda - (1/b) * eta;                
            case 'general'
                stationary_point_a = []; % To Be Done
                stationary_point_b = []; 
        end
        % smoothing omega = max(0, stationary_point) by CHKS function
        omega_a = 0.5*(sqrt(stationary_point_a.^2 + 4 * epsilon^2) + stationary_point_a);
        omega_b = 0.5*(sqrt(stationary_point_b.^2 + 4 * epsilon^2) + stationary_point_b);
        q_a = self.d_func(lambda) - self.d_func(omega_a) + self.d_grad(lambda) * (omega_a - lambda);
        q_b = self.d_func(lambda) - self.d_func(omega_b) + self.d_grad(lambda) * (omega_b - lambda);
        phi_ab = (omega_b - omega_a)' * eta + a * q_a - b * q_b;
        phi_func = Function('phi_func', {lambda, eta}, {phi_ab}, {'lambda', 'eta'}, {'phi_ab'});

    case 'finitely_representable'
        % omega does not have an explicit expression
        omega_a = SX.sym('omega_a', OCPEC.Dim.lambda, 1); 
        omega_b = SX.sym('omega_b', OCPEC.Dim.lambda, 1); 
        q_a = self.d_func(lambda) - self.d_func(omega_a) + self.d_grad(lambda) * (omega_a - lambda);
        q_b = self.d_func(lambda) - self.d_func(omega_b) + self.d_grad(lambda) * (omega_b - lambda);
        phi_ab = (omega_b - omega_a)' * eta + a * q_a - b * q_b;
        phi_func = Function('phi_func', {lambda, eta, omega_a, omega_b}, {phi_ab},...
            {'lambda', 'eta', 'omega_a', 'omega_b'}, {'phi_ab'});

    case 'polyhedral'
        % omega does not have an explicit expression
        omega_a = SX.sym('omega_a', OCPEC.Dim.lambda, 1); 
        omega_b = SX.sym('omega_b', OCPEC.Dim.lambda, 1); 
        q_a = self.d_func(lambda) - self.d_func(omega_a) + self.d_grad(lambda) * (omega_a - lambda);
        q_b = self.d_func(lambda) - self.d_func(omega_b) + self.d_grad(lambda) * (omega_b - lambda);
        phi_ab = (omega_b - omega_a)' * eta + a * q_a - b * q_b;
        phi_func = Function('phi_func', {lambda, eta, omega_a, omega_b}, {phi_ab},...
            {'lambda', 'eta', 'omega_a', 'omega_b'}, {'phi_ab'});
        
end

end