function penalty_func = create_penalty_func(self)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here

import casadi.*

% symbolic variable
w = SX.sym('w', 1, 1); % auxiliabry variable for gap function
mu = SX.sym('mu', 1, 1); % penalty parameter

% L1 penalty function (To Do: L2, barrier penalty function)
penalty = mu * w;

% output
penalty_func = Function('penalty_func', {w, mu}, {penalty}, {'w', 'mu'}, {'penalty'});

end