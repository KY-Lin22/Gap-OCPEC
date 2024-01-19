clear all
clc

%% generate problem and solver
relax_prob_type = 'Scholtes'; % 'primal_gap', 'D_gap'
                                % 'Scholtes', 'Lin_Fukushima', 'Kadrani'
                                % 'Kadrani', 'Kanzow_Schwartz'
                                % 'Steffensen_Ulbrich'
% two-dimension toy problem
Prob = generate_problem(relax_prob_type);

% option
Option = struct;
Option.print_time = false;
Option.record_time = true;
Option.ipopt.print_level = 0;
% Option.ipopt.tol = 1e-8; % default 1e-8
% Option.ipopt.max_iter = 3000; % default 3000
solver = casadi.nlpsol('solver', 'ipopt', Prob, Option);

%% single solve
x_0 = [-1; -1];
s_0 = 1e-8;
output = solve_problem_homotopy(solver, Prob, x_0, s_0);

%% solve problem and check stationarity
% generate a set of initial guess
lambda_sequ = -1 : 0.5 : 2;  
eta_sequ = -1 : 0.5 : 2;  
Rec.x_Init = cell(length(lambda_sequ), length(eta_sequ));
for i = 1 : length(lambda_sequ)
    for j = 1 : length(eta_sequ)
        Rec.x_Init{i, j} = [lambda_sequ(i); eta_sequ(j)];
    end   
end
Rec.x_Opt = cell(length(lambda_sequ), length(eta_sequ));
Rec.x_Opt_type = cell(length(lambda_sequ), length(eta_sequ)); 
Rec.x_Opt_color = cell(length(lambda_sequ), length(eta_sequ)); 

%% solve problem in a homotopy manner with a given initial guess
s_0 = 1e-1;
for i = 1 : length(lambda_sequ)
    for j = 1 : length(eta_sequ)
        x_0 = Rec.x_Init{i, j};        
        output = solve_problem_homotopy(solver, Prob, x_0, s_0);
        Rec.x_Opt{i, j} = output.x_Opt;
        Rec.x_Opt_type{i, j} = output.x_Opt_type;
        Rec.x_Opt_color{i, j} = output.x_Opt_color;
    end   
end

%% save data
% save('Data_Scholtes.mat', 'Rec')
save('Data_primal_gap.mat', 'Rec')
% save('Data_D_gap.mat', 'Rec')
% save('Data_Lin_Fukushima', 'Rec')
% save('Data_Kadrani', 'Rec')