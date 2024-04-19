function nlp = create_gap_constraint_based_NLP(self, OCPEC)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch self.gap_constraint_relaxation_strategy
    case 'generalized_primal_gap'

        switch self.auxiliary_variable_strategy
            case 'none'
                nlp = self.create_primal_gap_NLP(OCPEC);
            case 'primal_gap_func'
                nlp = self.create_primal_gap_NLP_w(OCPEC);
            case 'omega'
                % TO DO
        end

    case 'generalized_D_gap'

        switch self.auxiliary_variable_strategy
            case 'none'
                nlp = self.create_D_gap_NLP(OCPEC);
            case 'primal_gap_func'
                nlp = self.create_D_gap_NLP_w_v(OCPEC);
            case 'omega'
                % TO DO
        end
end


end