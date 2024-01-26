classdef NLP_Formulation < handle
    % formulate a NLP based on the given OCPEC and formulation option
    %
    % OCPEC has the form:
    %  min  L_T(x) + int_0^T L_S(x, u, lambda) dt,
    %  s.t. Dot{x} = f(x, u, lambda)
    %       lambda \in SOL(K, F(x, u, lambda))
    %       K = {lambda | g(lambda) >= 0}
    %       G(x, u) >= 0,
    %       C(x, u) = 0,    
    %
    % NLP has the form:
    %  min  J(z, p),
    %  s.t. h(z, p) = 0,
    %       c(z, p) >= 0
    %
    properties
        relaxProbType char {mustBeMember(relaxProbType, {...
            'generalized_primal_gap_constraint_based',...
            'generalized_D_gap_constraint_based'...
            'KKT_based'...
            })} = 'generalized_primal_gap_constraint_based' 
        D_gap_param = struct('a', 0.9, 'b', 1.1) % D gap function parameters: b > a > 0 (a ref value: a = 0.9. b = 1.1)
                                                 % Ref: Theoretical and numerical investigation of the D-gap function for BVI, 
                                                 % 1998, Mathematical Programming, C.Kanzow & M. Fukushima        
        stronglyConvexFuncType char {mustBeMember(stronglyConvexFuncType, {...
            'quadratic',...
            'general'...
            })} = 'quadratic' 
        KKT_relaxation_strategy char{mustBeMember(KKT_relaxation_strategy, {...
            'Scholtes',...
            'Lin_Fukushima',...
            'Kadrani'})} = 'Scholtes'
        d_func % function object, strongly convex function d  
        d_grad % function object, gradient of the strongly convex function d
        d_hessian % function object, hessian of the strongly convex function d  
        
        OmegaEvalProb % struct, the strongly concave maximization problem to evaluate omega
        
        z % symbolic variable, includes all the variable to be optimized,
        p % symbolic variable, including all the problem parameter
        J % symbolic function, cost function 
        h % symbolic function, equality constraint    
        c % symbolic function, inequality constraint   
        Dim % struct, problem dimension record
        
        FuncObj % structure, function object  
    end
    
    %% Constructor method        
    methods
        function self = NLP_Formulation(OCPEC, Option)
            %NLP_Formulation: Construct an instance of this class
            %   Detailed explanation goes here
            %% specify properties
            if ~isempty(Option.relaxProbType)
                self.relaxProbType = Option.relaxProbType;
            end
            if ~isempty(Option.stronglyConvexFuncType)
                self.stronglyConvexFuncType = Option.stronglyConvexFuncType;
            end
            if ~isempty(Option.D_gap_param.a)
                self.D_gap_param.a = Option.D_gap_param.a;
            end
            if ~isempty(Option.D_gap_param.b)
                self.D_gap_param.b = Option.D_gap_param.b;
            end            
            if ~isempty(Option.KKT_relaxation_strategy)
                self.KKT_relaxation_strategy = Option.KKT_relaxation_strategy;
            end

            %% create a Strongly Convex Function d and its derivative used in the gap functions PHI (or PHIab)
            [d_func, d_grad, d_hessian] = self.createStronglyConvexFunction(OCPEC);
            self.d_func = d_func;
            self.d_grad = d_grad;
            self.d_hessian = d_hessian;
            
            %% create a Strongly Concave Maximization Problem to evaluate variable omega (or omega_a, omega_b) used in the gap function PHI (or PHIab)                     
            if strcmp(self.relaxProbType, 'KKT_based')
                self.OmegaEvalProb = [];
            else
                self.OmegaEvalProb = self.createStronglyConcaveMaxProblem(OCPEC);
            end
            
            %% discretize OCPEC into NLP 
            switch self.relaxProbType
                case 'generalized_primal_gap_constraint_based'
                    nlp = self.createNLP_Primal_Gap_Based(OCPEC);
                case 'generalized_D_gap_constraint_based'
                    nlp = self.createNLP_D_Gap_Based(OCPEC);
                case 'KKT_based'
                    nlp = self.createNLP_KKT_Based(OCPEC);
            end
            self.z = nlp.z;   
            self.p = nlp.p;
            self.J = nlp.J;    
            self.h = nlp.h;  
            self.c = nlp.c;        
            self.Dim = nlp.Dim;   
            
            %% create function object
            self.FuncObj = self.createFuncObj(nlp);
            
        end
    end
    
    %% Other method
    methods
        [d_func, d_grad, d_hessian] = createStronglyConvexFunction(self, OCPEC)
        
        OmegaEvalProb = createStronglyConcaveMaxProblem(self, OCPEC)
        
        nlp = createNLP_Primal_Gap_Based(self, OCPEC);
        
        nlp = createNLP_D_Gap_Based(self, OCPEC);
        
        nlp = createNLP_KKT_Based(self, OCPEC);
        
        FuncObj = createFuncObj(self, nlp) 
    end
    
end

