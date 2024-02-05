function solver_set = create_solver_set(OCPEC_func_handle, nStages_sequ, NLP_option_set, solver_Option)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%% generate OCPEC problem set
OCPEC_problem_set = cell(numel(OCPEC_func_handle), numel(nStages_sequ));
for i = 1 : numel(OCPEC_func_handle)
    for j = 1 : numel(nStages_sequ)
        OCPEC_i_j = OCPEC_func_handle{i}();
        OCPEC_i_j.nStages = nStages_sequ{j};
        OCPEC_problem_set{i, j} = OCPEC_i_j;
    end
end
OCPEC_problem_set = reshape(OCPEC_problem_set, [], 1);

%% generate relaxed NLP problem set
NLP_problem_set = cell(numel(OCPEC_problem_set), numel(NLP_option_set));
for i = 1 : numel(OCPEC_problem_set)
    for j = 1 : numel(NLP_option_set)
        OCPEC_i = OCPEC_problem_set{i};
        NLP_option_j = NLP_option_set{j};
        NLP_i_j = NLP_Formulation(OCPEC_i, NLP_option_j);
        NLP_problem_set{i, j} = NLP_i_j;
    end
end

%% generate IPOPT based solver set
solver_set = cell(size(NLP_problem_set));
for i = 1 : size(NLP_problem_set, 1)
    for j = 1 : size(NLP_problem_set, 2)
        OCPEC_i = OCPEC_problem_set{i};
        NLP_i_j = NLP_problem_set{i, j};
        solver_i_j = IPOPT_Based_Solver(OCPEC_i, NLP_i_j);
        % option
        solver_i_j.Option.Homotopy.kappa_s_times = solver_Option.kappa_s_times;
        solver_i_j.Option.Homotopy.kappa_s_exp = solver_Option.kappa_s_exp;
        solver_i_j.Option.Homotopy.VI_nat_res_tol = solver_Option.VI_nat_res_tol;
        solver_set{i, j} = solver_i_j;
    end
end

end