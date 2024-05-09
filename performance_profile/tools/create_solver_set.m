function solver_set = create_solver_set(OCPEC_func_handle, nStages_sequ, NLP_option_set, solver_Option_set)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% check input
if numel(NLP_option_set) ~= numel(solver_Option_set)
    error('NLP_option_set and solver_Option_set should have the same number of element')
end

% generate OCPEC problem set
OCPEC_problem_set = cell(numel(OCPEC_func_handle), numel(nStages_sequ));
for i = 1 : numel(OCPEC_func_handle)
    for j = 1 : numel(nStages_sequ)
        OCPEC_i_j = OCPEC_func_handle{i}();
        OCPEC_i_j.nStages = nStages_sequ{j};
        OCPEC_i_j.timeStep = OCPEC_i_j.timeHorizon ./ OCPEC_i_j.nStages;
        OCPEC_problem_set{i, j} = OCPEC_i_j;
    end
end
OCPEC_problem_set = reshape(OCPEC_problem_set, [], 1);

% generate relaxed NLP problem and IPOPT based solver set
solver_set = cell(numel(OCPEC_problem_set), numel(NLP_option_set));
for i = 1 : numel(OCPEC_problem_set)
    for j = 1 : numel(NLP_option_set)
        % load OCPEC problem, NLP and solver option
        OCPEC_i = OCPEC_problem_set{i};
        NLP_option_j = NLP_option_set{j};
        solver_option_j = solver_Option_set{j};
        % create NLP
        NLP_i_j = NLP_Formulation(OCPEC_i, NLP_option_j);
        % create solver
        solver_i_j = IPOPT_Based_Solver(OCPEC_i, NLP_i_j);      
        solver_i_j.Option.Homotopy.kappa_s_times = solver_option_j.kappa_s_times;
        solver_i_j.Option.Homotopy.VI_nat_res_tol = solver_option_j.VI_nat_res_tol;
        solver_set{i, j} = solver_i_j;
    end
end


end