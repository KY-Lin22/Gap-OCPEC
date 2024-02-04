function Rec = run_test(solver_set)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% parameter
s_Init = 1e0;
s_End = 1e-8;
mu_Init = 1e0;
mu_End = 1e4;
p_Init = [s_Init; mu_Init];
p_End = [s_End; mu_End];
% record cost and time
Rec.cost = zeros(size(solver_set));
Rec.time = zeros(size(solver_set));
% solve
for i = 1 : size(solver_set, 1)
    for j = 1 : size(solver_set, 2)
        solver_i_j = solver_set{i, j};
        z_Init_i_j = zeros(solver_i_j.NLP.Dim.z, 1);
        [~, Info_i_j] = solver_i_j.solve_NLP(z_Init_i_j, p_Init, p_End);
        if Info_i_j.terminalStatus == 1
            Rec.cost(i, j) = Info_i_j.cost.ocp;
            Rec.time(i, j) = Info_i_j.time;
        else
            % set fail case as inf
            Rec.cost(i, j) = inf;
            Rec.time(i, j) = inf;
        end
    end
end

end