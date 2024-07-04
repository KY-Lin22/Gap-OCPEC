function Rec = run_solver_test(solver_set, param_set)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% check input
if numel(param_set) ~= size(solver_set, 2)
    error('number of element in param_set should equal to the number of method (column) in solver_set ')
end

% init record
Rec.cost = zeros(size(solver_set));
Rec.KKT_error = zeros(size(solver_set));
Rec.VI_nat_res = zeros(size(solver_set));
Rec.time = zeros(size(solver_set));
Rec.time_1 = zeros(size(solver_set));
Rec.continuationStepNum = zeros(size(solver_set));

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
        s_Init_j = param_j.s_Init;
        s_End_j = param_j.s_End;
        % solve and record
        [~, Info_i_j] = solver_i_j.solve_NLP(z_Init_i_j, s_Init_j, s_End_j);
        if Info_i_j.terminal_status == 1
            Rec.cost(i, j) = Info_i_j.cost;
            Rec.KKT_error(i, j) = Info_i_j.KKT_error;
            Rec.VI_nat_res(i, j) = Info_i_j.VI_natural_residual;

            time_i_j = Info_i_j.time;
            timeLog_i_j = Info_i_j.Log.time;
            continuationStepNum_i_j = Info_i_j.continuationStepNum;

            Rec.time(i, j) = time_i_j;         
            Rec.time_1(i, j) = timeLog_i_j(1);
            Rec.continuationStepNum(i, j) = continuationStepNum_i_j;                     
        else
            % set fail case as inf
            Rec.cost(i, j) = inf;
            Rec.KKT_error(i, j) = inf;
            Rec.VI_nat_res(i, j) = inf;

            Rec.time(i, j) = inf;
            Rec.time_1(i, j) = inf;
            Rec.continuationStepNum(i, j) = inf;
        end
    end
end

end