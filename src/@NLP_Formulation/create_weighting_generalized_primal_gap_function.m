function phi_c_func = create_weighting_generalized_primal_gap_function(self, OCPEC)
%create weighting generalized primal gap function phi_c
%   phi_c is a function of lambda, eta, and an intermedia variable omega, 
%   where omega also is a function of lambda and eta.
%   Therefore, according to the type of VI set in OCPEC, phi_c may be:
%   (1) a function of lambda and eta,
%   where omega CAN be explicitly represented by lambda and eta
%   (2) a function of lambda, eta, and omega
%   where omega CAN NOT be explicitly represented by lambda and eta

import casadi.*
% variable
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);
eta = SX.sym('eta', OCPEC.Dim.lambda, 1); 
% parameter
c = SX.sym('c', 1, 1); % weighting nonnegative scalar parameter 
epsilon = self.gap_func_smooth_param;
% strongly convex function d and its derivative
[d_func, d_grad, ~] = self.create_strongly_convex_func(OCPEC);

%% formulate phi_c
switch OCPEC.VISetType
    case 'box_constraint'
        % omega has an explicit expression: projection of the stationary point into [bl, bu]
        switch self.strongly_convex_func
            case 'quadratic'
                stationary_point_c = lambda - (1/c) * eta;           
        end
        % smoothing omega = mid(bl, bu, stationary_point) by CHKS function 
        % ref: example 2 in ''A new look at smoothing Newton methods for NCPs and BVI'', Mathematical Programming, 2000
        omega_c = 0.5*(OCPEC.bl + sqrt((OCPEC.bl - stationary_point_c).^2 + 4*epsilon^2)) ...
            + 0.5*(OCPEC.bu - sqrt((OCPEC.bu - stationary_point_c).^2 + 4*epsilon^2));
        p_c = d_func(lambda) - d_func(omega_c) + d_grad(lambda) * (omega_c - lambda);
        phi_c = eta' * (lambda - omega_c) + c * p_c;
        phi_c_func = Function('phi_c_func', {lambda, eta, c}, {phi_c}, {'lambda', 'eta', 'c'}, {'phi_c'});

    case 'nonnegative_orthant'
        % omega has an explicit expression: projection of the stationary point into [0, inf]
        switch self.strongly_convex_func
            case 'quadratic'
                stationary_point_c = lambda - (1/c) * eta;
        end
        % smoothing omega = max(0, stationary_point) by CHKS function
        omega_c = 0.5*(sqrt(stationary_point_c.^2 + 4 * epsilon^2) + stationary_point_c);
        p_c = d_func(lambda) - d_func(omega_c) + d_grad(lambda) * (omega_c - lambda);
        phi_c = eta' * (lambda - omega_c) + c * p_c;
        phi_c_func = Function('phi_c_func', {lambda, eta, c}, {phi_c}, {'lambda', 'eta', 'c'}, {'phi_c'});

    case 'finitely_representable'
        % omega does not have an explicit expression
        omega_c = SX.sym('omega_c', OCPEC.Dim.lambda, 1);
        p_c = d_func(lambda) - d_func(omega_c) + d_grad(lambda) * (omega_c - lambda);
        phi_c = eta' * (lambda - omega_c) + c * p_c;
        phi_c_func = Function('phi_c_func', {lambda, eta, omega_c, c}, {phi_c}, {'lambda', 'eta', 'omega_c', 'c'}, {'phi_c'});

    case 'polyhedral'
        % omega does not have an explicit expression
        omega_c = SX.sym('omega_c', OCPEC.Dim.lambda, 1);
        p_c = d_func(lambda) - d_func(omega_c) + d_grad(lambda) * (omega_c - lambda);
        phi_c = eta' * (lambda - omega_c) + c * p_c;
        phi_c_func = Function('phi_c_func', {lambda, eta, omega_c, c}, {phi_c}, {'lambda', 'eta', 'omega_c', 'c'}, {'phi_c'});
        
end

end