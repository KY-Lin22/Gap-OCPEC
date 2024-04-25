function Rec = run_solver_test(solver_set, param_set)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% check input
if numel(param_set) ~= size(solver_set, 2)
    error('number of element in param_set should equal to the number of method (column) in solver_set ')
end

% init record
Rec.cost = zeros(size(solver_set));
Rec.time = zeros(size(solver_set));
Rec.dual_var_inf_norm = zeros(size(solver_set));
% solve
for i = 1 : size(solver_set, 1)
    for j = 1 : size(solver_set, 2)
        % show info
        disp(['-----------------------------',...
            'prob: ', num2str(i), ' / ', num2str(size(solver_set, 1)), '; ', ...
            'solver: ',  num2str(j), ' / ', num2str(size(solver_set, 2)),...
            '-----------------------------'])
        % load solver and parameter
        solver_i_j = solver_set{i, j};
        z_Init_i_j = ones(solver_i_j.NLP.Dim.z, 1);      
        param_j = param_set{j};
        p_Init_j = param_j.p_Init;
        p_End_j = param_j.p_End;
        % solve and record
        [~, Info_i_j] = solver_i_j.solve_NLP(z_Init_i_j, p_Init_j, p_End_j);
        if Info_i_j.terminalStatus == 1
            Rec.cost(i, j) = Info_i_j.cost;
            Rec.time(i, j) = Info_i_j.time;
            Rec.dual_var_inf_norm(i, j) = norm(Info_i_j.dual_var, inf);
        else
            % set fail case as inf
            Rec.cost(i, j) = inf;
            Rec.time(i, j) = inf;
            Rec.dual_var_inf_norm(i, j) = inf;
        end
    end
end

end