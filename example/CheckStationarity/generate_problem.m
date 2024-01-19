function Prob = generate_problem(relax_prob_type)
%UNTITLED3 Summary of this function goes here
%   generate a two dimension toy problem in the form of:
%   min J = (lambda - 1)^2 + (eta - 1)^2
%   s.t. 0 <= lambda \perp eta >= 0
% ref to C.Kanzow et,al:  A new regularization method for mathematical programs with complementarity 
%                         constraints with strong convergence properties (2010, report version)
import casadi.*

%% problem variable, parameter, and cost
% variables
lambda = SX.sym('lambda', 1, 1);
eta = SX.sym('eta', 1, 1);
% parameter
s = SX.sym('s', 1, 1);
% cost
J = (lambda - 1)^2 + (eta - 1)^2;

%% problem constraint: g >= 0
switch relax_prob_type
    case 'primal_gap'
        phi = 0.5 * (eta^2 - (max(0, eta - lambda))^2);
        g = [lambda;...
            s - phi];       
    case 'D_gap'    
        b = 1.1;
        a = 0.9;
        phi_ab = (b-a)/(2*a*b)*eta^2 - 1/(2*a)*(max(0, eta-a*lambda))^2 + 1/(2*b)*(max(0, eta-b*lambda))^2;       
        g = s - phi_ab;
    case 'Scholtes'
        g = [lambda;...
            eta;...
            s - lambda * eta];
    case 'Lin_Fukushima'
        g = [s^2 - lambda * eta;...
            (lambda + s) * (eta + s) - s^2];
    case 'Kadrani'
        g = [lambda + s;...
            eta + s;...
            -(lambda - s)*(eta - s)];
    case 'Kanzow_Schwartz'
        
    case 'Steffensen_Ulbrich'
        
end

%% output
Prob = struct('x', [lambda; eta], 'f', J, 'g', g, 'p', s);

end

