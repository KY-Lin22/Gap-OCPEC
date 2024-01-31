function z_Init = create_initial_guess(self)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% z_Init = randn(NLP.Dim.z, 1);

OCPEC = self.OCPEC;
NLP = self.NLP;
% x, u, lambda
X_Init = randn(OCPEC.Dim.x, OCPEC.nStages);
U_Init = randn(OCPEC.Dim.u, OCPEC.nStages);
LAMBDA_Init = randn(OCPEC.Dim.lambda, OCPEC.nStages);
% auxiliary variable
switch NLP.relaxation_problem
    case 'gap_constraint_based'
        F_FuncObj_map = OCPEC.FuncObj.F.map(OCPEC.nStages);
        ETA_Init = full(F_FuncObj_map(X_Init, U_Init, LAMBDA_Init));
        switch NLP.gap_constraint_relaxation_strategy
            case 'generalized_primal_gap'
                W_Init = zeros(1, OCPEC.nStages);
                z_Init = reshape([X_Init; U_Init; LAMBDA_Init; ETA_Init; W_Init], [], 1);
            case 'generalized_D_gap'
                W_Init = zeros(1, OCPEC.nStages);
                V_Init = zeros(1, OCPEC.nStages);
                z_Init = reshape([X_Init; U_Init; LAMBDA_Init; ETA_Init; W_Init; V_Init], [], 1);
        end
    case 'KKT_based'
        ZETA_Init = zeros(OCPEC.Dim.g, OCPEC.nStages);
        g_FuncObj_map = OCPEC.FuncObj.g.map(OCPEC.nStages);
        W_Init = full(g_FuncObj_map(LAMBDA_Init));
        z_Init = reshape([X_Init; U_Init; LAMBDA_Init; ZETA_Init; W_Init], [], 1);
end

end