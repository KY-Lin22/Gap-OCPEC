function penalty_func = create_penalty_func(self)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here

import casadi.*

% symbolic variable
w = SX.sym('w', 1, 1); % scalar auxiliabry variable for gap function
% penalty function formulation
switch self.gap_func_auxiliary_variable_penalty
    case 'none'
        penalty = 0;
    case 'L2'
        penalty = w^2;
end
% output
penalty_func = Function('penalty_func', {w}, {penalty}, {'w'}, {'penalty'});

end