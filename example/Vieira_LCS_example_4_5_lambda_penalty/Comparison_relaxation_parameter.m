clear all
clc

%% create OCPEC, NLP, and solver
% construct OCPEC problem
OCPEC = OCPEC_Vieira_LCS_lambda_penalty();

% discretize OCPEC into a NLP problem
Option.relaxation_problem = 'gap_constraint_based'; % 'gap_constraint_based', 'KKT_based'
Option.gap_constraint_relaxation_strategy = 'generalized_D_gap'; % 'generalized_primal_gap', 'generalized_D_gap'
Option.KKT_complementarity_relaxation_strategy = 'Scholtes'; % 'Scholtes', 'Lin_Fukushima', 'Kadrani'
Option.strongly_convex_func = 'quadratic';
Option.gap_func_smooth_param = 0.001;
Option.D_gap_param_a = 0.8;
Option.D_gap_param_b = 1.2;
Option.penalty_gap_func_auxiliary_variable = 'L2'; % 'none', 'L1', 'L2'
NLP = NLP_Formulation(OCPEC, Option);

% create solver
solver = IPOPT_Based_Solver(OCPEC, NLP);

% create initial guess
% z_Init = solver.create_initial_guess();

%% parameter and problem solve
% parameter
solver.Option.Homotopy.kappa_s_times = 0.8;
solver.Option.Homotopy.kappa_s_exp = 1.2;

s_Init = 1e-1;
s_End_sequ = [...
    9e-2, 7e-2, 5e-2, 3e-2, 1e-2,...
    9e-3, 7e-3, 5e-3, 3e-3, 1e-3,... 
    9e-4, 7e-4, 5e-4, 3e-4, 1e-4,...
    9e-5, 7e-5, 5e-5, 3e-5, 1e-5,...
    9e-6, 7e-6, 5e-6, 3e-6, 1e-6];
mu_Init = 1e0;
mu_End = 1e3;

Rec_size = numel(s_End_sequ);
Rec.s_End = s_End_sequ;
Rec.cost = zeros(1, Rec_size);
Rec.nat_res = zeros(1, Rec_size);
Rec.time = zeros(1, Rec_size);

z_Init = zeros(NLP.Dim.z, 1);
p_Init = [s_Init; mu_Init];

num_test = 5;
for j = 1 : Rec_size
    cost_test = 0;
    nat_res_test = 0;
    time_test = 0;
    s_End_j = s_End_sequ(j);
    p_End_j = [s_End_j; mu_End];
    for k = 1 : num_test
        [z_Opt, Info] = solver.solve_NLP(z_Init, p_Init, p_End_j);
        if ~strcmp(Info.terminalMsg, 'Solve_Succeeded')
            error('solver fail')
        end
        cost_test = cost_test + Info.cost.ocp;
        nat_res_test = nat_res_test + Info.VI_natural_residual;
        time_test = time_test + Info.time;
    end
    Rec.cost(j) = cost_test/num_test;
    Rec.nat_res(j) = nat_res_test/num_test;
    Rec.time(j) = time_test/num_test;
    disp(['s_End: ', num2str(s_End_j, '%10.2e'), '; ',...
        'cost: ', num2str(Rec.cost(j), '%10.2e'), '; ',...
        'natRes: ', num2str(Rec.nat_res(j), '%10.2e'), '; ',...
        'time: ', num2str(Rec.time(j), '%10.2f'), ' s'])    
end
%%
save('Data_primal_gap.mat', 'Rec')
