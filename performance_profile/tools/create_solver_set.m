function solver_set = create_solver_set(OCPEC_func_handle, nStages_sequ, NLP_option_set, solver_option_set, solver_type)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% check input
if numel(NLP_option_set) ~= numel(solver_option_set)
    error('NLP_option_set and Solver_option_set should have the same number of element')
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
    % problem
    OCPEC_i = OCPEC_problem_set{i};
    for j = 1 : numel(NLP_option_set)
        % method       
        NLP_option_j = NLP_option_set{j};
        Solver_option_j = solver_option_set{j};
        % create NLP
        NLP_i_j = NLP_Formulation(OCPEC_i, NLP_option_j);
        % create solver
        Solver_type_j = solver_type{j};
        switch Solver_type_j
            case 'IPOPT-Based'
                solver_i_j = IPOPT_Based_Solver(OCPEC_i, NLP_i_j, Solver_option_j);
            case 'DynSys-Based'
                solver_i_j = DynSys_Based_Solver(OCPEC_i, NLP_i_j, Solver_option_j);
            otherwise
                error('specified solver type is not supported')
        end
              
        solver_set{i, j} = solver_i_j;
    end
end


end