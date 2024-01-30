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
        relaxation_problem char {mustBeMember(relaxation_problem, {...
            'gap_constraint_based',...
            'KKT_based'...
            })} = 'gap_constraint_based' 
        gap_constraint_relaxation_strategy char {mustBeMember(gap_constraint_relaxation_strategy, {...
            'generalized_primal_gap',...
            'generalized_D_gap'})} = 'generalized_primal_gap'
        KKT_complementarity_relaxation_strategy char {mustBeMember(KKT_complementarity_relaxation_strategy, {...
            'Scholtes',...
            'Lin_Fukushima',...
            'Kadrani'})} = 'Scholtes'

        strongly_convex_func char {mustBeMember(strongly_convex_func, {...
            'quadratic',...
            'general'...
            })} = 'quadratic'  % strongly convex function in Auchmuty's saddle function       
        gap_func_smooth_param double {mustBeNonnegative} = 0.001 % used in CHKS smoothing function for max(0, x)
        D_gap_param_a double {mustBeNonnegative} = 0.9; % D gap function parameters: b > a > 0 (a ref value: a = 0.9. b = 1.1)
        D_gap_param_b double {mustBeNonnegative} = 1.1; % Ref: Theoretical and numerical investigation of the D-gap function   
                                                        % for BVI, 1998, Mathematical Programming, C.Kanzow & M. Fukushima                  
        penalty_gap_func_auxiliary_variable char {mustBeMember(penalty_gap_func_auxiliary_variable, {...
            'none',...
            'L1'})} = 'none' % To Do: L2, barrier function

        state_equation_discretization char {mustBeMember(state_equation_discretization, {...
            'implicit_Euler'})} = 'implicit_Euler'
    end

    properties
        discre_state_equ_func % function object, discretization state equation

        d_func % function object, strongly convex function d  
        d_grad % function object, gradient of the strongly convex function d
        d_hessian % function object, hessian of the strongly convex function d   
        OmegaEvalProb % struct, strongly concave maximization problem to evaluate omega
        penalty_func % function object, penalty function for certain variables
        gap_func % function object, gap function for VI

        KKT_stationarity_func % function object, KKT stationarity condition for VI
        KKT_complementarity_func % function object, KKT complementarity condition for VI 
    end

    properties
        z % symbolic variable, includes all the variable to be optimized,
        p % symbolic variable, including all the problem parameter
        J % symbolic function, cost function 
        h % symbolic function, equality constraint    
        c % symbolic function, inequality constraint   
        Dim % struct, problem dimension record
        
        FuncObj % structure, NLP function object  
    end
    
    %% Constructor method        
    methods
        function self = NLP_Formulation(OCPEC, Option)
            %NLP_Formulation: Construct an instance of this class
            %   Detailed explanation goes here
            %% specify properties based on Option
            if isfield(Option, 'relaxation_problem')
                self.relaxation_problem = Option.relaxation_problem;
            end
            if isfield(Option, 'gap_constraint_relaxation_strategy')
                self.gap_constraint_relaxation_strategy = Option.gap_constraint_relaxation_strategy;
            end
            if isfield(Option, 'KKT_complementarity_relaxation_strategy')
                self.KKT_complementarity_relaxation_strategy = Option.KKT_complementarity_relaxation_strategy;
            end
            if isfield(Option, 'strongly_convex_func')
                self.strongly_convex_func = Option.strongly_convex_func;
            end
            if isfield(Option, 'gap_func_smooth_param')
                self.gap_func_smooth_param = Option.gap_func_smooth_param;
            end
            if isfield(Option, 'D_gap_param_a')
                self.D_gap_param_a = Option.D_gap_param_a;
            end
            if isfield(Option, 'D_gap_param_b')
                self.D_gap_param_b = Option.D_gap_param_b;
            end
            if self.D_gap_param_b <= self.D_gap_param_a
                error('D gap function parameter should satisfy: b > a > 0')
            end
            if isfield(Option, 'penalty_gap_func_auxiliary_variable')
                self.penalty_gap_func_auxiliary_variable = Option.penalty_gap_func_auxiliary_variable;
            end
            if isfield(Option, 'state_equation_discretization')
                self.state_equation_discretization = Option.state_equation_discretization;
            end

            %% specify properties about function object used in NLP reformulation
            % create discretization state equation
            self.discre_state_equ_func = self.create_discre_state_equ_func(OCPEC);
                    
            switch self.relaxation_problem
                case 'gap_constraint_based'
                    % create a strongly convex function d and its derivative used in gap function
                    [d_func, d_grad, d_hessian] = self.create_strongly_convex_func(OCPEC);
                    self.d_func = d_func;
                    self.d_grad = d_grad;
                    self.d_hessian = d_hessian;
                    % create a strongly concave maximization problem to evaluate variable omega
                    % used in gap function
                    self.OmegaEvalProb = self.create_strongly_concave_max_prob(OCPEC);
                    % create penalty function for auxiliary variable
                    self.penalty_func = self.create_penalty_func();
                    % create function object to evaluate gap function phi, 
                    % phi is a function of lambda, eta, and an intermedia variable omega
                    % omega also is a function of lambda and eta.
                    % Therefore, according to the type of the OCPEC VI set, the return 
                    % function object about phi may be:
                    % (1) a function of lambda and eta,
                    % where omega CAN be explicitly represented by lambda and eta
                    % (2) a function of lambda, eta, and omega
                    % where omega CAN NOT be explicitly represented by lambda and eta
                    switch self.gap_constraint_relaxation_strategy
                        case 'generalized_primal_gap'
                            self.gap_func = self.create_generalized_primal_gap_function(OCPEC);
                        case 'generalized_D_gap'
                            self.gap_func = self.create_generalized_D_gap_function(OCPEC);
                    end

                case 'KKT_based'
                    % create function object for KKT based VI reformulation
                    [KKT_stationarity_func, KKT_complementarity_func] = self.create_KKT_reformulation(OCPEC);
                    self.KKT_stationarity_func = KKT_stationarity_func;
                    self.KKT_complementarity_func = KKT_complementarity_func;
            end

            %% discretize OCPEC into NLP
            switch self.relaxation_problem
                case 'gap_constraint_based'
                    nlp = self.create_gap_constraint_based_NLP(OCPEC);
                case 'KKT_based'
                    nlp = self.create_KKT_based_NLP(OCPEC);
            end
            self.z = nlp.z;   
            self.p = nlp.p;
            self.J = nlp.J;    
            self.h = nlp.h;  
            self.c = nlp.c;        
            self.Dim = nlp.Dim;   
            
            %% create function object
            self.FuncObj = self.create_FuncObj(nlp);
            
        end
    end
    
    %% Other method
    methods
        discre_state_equ_func = create_discre_state_equ_func(self, OCPEC)

        [d_func, d_grad, d_hessian] = create_strongly_convex_func(self, OCPEC)
        
        OmegaEvalProb = create_strongly_concave_max_prob(self, OCPEC)

        penalty_func = create_penalty_func(self)

        phi_func = create_generalized_primal_gap_function(self, OCPEC)

        phi_func = create_generalized_D_gap_function(self, OCPEC)

        [KKT_stationarity_func, KKT_complementarity_func] = create_KKT_reformulation(self, OCPEC)

        nlp = create_gap_constraint_based_NLP(self, OCPEC)

        nlp = create_KKT_based_NLP(self, OCPEC)
        
        FuncObj = create_FuncObj(self, nlp) 
    end
    
end

