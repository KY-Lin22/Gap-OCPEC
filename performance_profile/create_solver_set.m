function solver_set = create_solver_set(OCPEC_problem_set, NLP_problem_set)
%UNTITLED28 Summary of this function goes here
%   Detailed explanation goes here

solver_set = cell(size(NLP_problem_set));

for i = 1 : size(NLP_problem_set, 1)
    for j = 1 : size(NLP_problem_set, 2)
        OCPEC_i = OCPEC_problem_set{i};
        NLP_i_j = NLP_problem_set{i, j};
        solver_i_j = IPOPT_Based_Solver(OCPEC_i, NLP_i_j);
        % option
        solver_i_j.Option.Homotopy.kappa_s_times = 0.8;
        solver_i_j.Option.Homotopy.kappa_s_exp = 1.2;
        solver_i_j.Option.Homotopy.VI_nat_res_tol = 1e-2;
        solver_set{i, j} = solver_i_j;
    end
end

end