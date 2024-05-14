
clear all
clc
% create function object using callback
X = casadi.MX.sym('X', 2, 1);
Y = casadi.MX.sym('Y', 2, 1);
SQ = (X - 5)' * (X - 5) + (Y - 5)' * (Y - 5);
SQ_FuncObj = casadi.Function('Square_FuncObj', {X, Y}, {SQ}, {'X', 'Y'}, {'SQ'});
opts = struct('enable_fd',true);
f = NLP_callback('f', SQ_FuncObj, opts);

J = casadi.Function('J', {X, Y}, {jacobian(f(X, Y), X)});
disp(J)

%% nlp
x = casadi.MX.sym('x', 2, 1);
y = casadi.MX.sym('y', 2, 1);
cost = f(x, y) + 10;
disp(cost)
nlp = struct('x', [x; y], 'f', cost);
% jac_f = jacobian(y, x);

Option = struct();
% Option.ipopt.max_iter = 1000;
% Option.ipopt.hessian_approximation = 'limited-memory';
solver = casadi.nlpsol('solver', 'ipopt', nlp, Option);

%% solve
solution = solver('x0', [0; 0; 0; 0], 'lbx', [-100; -100; -100; -100], 'ubx', [100; 100; 100; 100]);
disp('optimal solution: ')
disp(full(solution.x))
disp('optimal cost: ')
disp(full(solution.f))

