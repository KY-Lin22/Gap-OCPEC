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
            'generalized_D_gap'...
            })} = 'generalized_primal_gap'
        KKT_complementarity_relaxation_strategy char {mustBeMember(KKT_complementarity_relaxation_strategy, {...
            'Scholtes',...
            'Lin_Fukushima',...
            'Kadrani'...
            })} = 'Scholtes'

        strongly_convex_func char {mustBeMember(strongly_convex_func, {...
            'quadratic',...
            })} = 'quadratic'  % strongly convex function in Auchmuty's saddle function       
        gap_func_smooth_param double {mustBeNonnegative} = 1e-4 % used in CHKS smoothing function for max(0, x) and mid(bl, bu, x)
        D_gap_param_a double {mustBeNonnegative} = 0.9; % D gap function parameters: b > a > 0 (a ref value: a = 0.9, b = 1.1)
        D_gap_param_b double {mustBeNonnegative} = 1.1; % Ref: Theoretical and numerical investigation of the D-gap function   
                                                        % for BVI, 1998, Mathematical Programming, C.Kanzow & M. Fukushima  
    end

    properties
        z % symbolic variable, includes all the variable to be optimized,
        p % symbolic variable, including all the problem parameter
        J % symbolic function, cost function 
        h % symbolic function, equality constraint    
        c % symbolic function, inequality constraint   
        Dim % struct, problem dimension record
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

            %% discretize OCPEC into NLP
            switch self.relaxation_problem
                case 'gap_constraint_based'
                    switch self.gap_constraint_relaxation_strategy
                        case 'generalized_primal_gap'
                            nlp = self.create_primal_gap_NLP(OCPEC);
                        case 'generalized_D_gap'
                            nlp = self.create_D_gap_NLP(OCPEC);
                    end
                case 'KKT_based'
                    nlp = self.create_KKT_based_NLP(OCPEC);
            end
            self.z = nlp.z;   
            self.p = nlp.p;
            self.J = nlp.J;    
            self.h = nlp.h;  
            self.c = nlp.c;        
            self.Dim = nlp.Dim;               

            %% display NLP information
            disp('*----------------------------------- NLP Information ------------------------------------*')
            disp('1. equilibrium constraint reformulation')
            disp(['relaxation problem: .................................. ', self.relaxation_problem])
            switch self.relaxation_problem
                case 'gap_constraint_based'
                    disp(['relaxation strategy: ................................. ', self.gap_constraint_relaxation_strategy])
                    disp(['strongly convex function: ............................ ', self.strongly_convex_func])
                    disp(['CHKS smoothing parameter: ............................ ', num2str(self.gap_func_smooth_param, '%10.3e')])
                    switch self.gap_constraint_relaxation_strategy
                        case 'generalized_primal_gap'

                        case 'generalized_D_gap'
                            disp(['D gap function parameter (a / b): .................... ', ...
                                num2str(self.D_gap_param_a), ' / ', num2str(self.D_gap_param_b)])
                    end
                case 'KKT_based'
                    disp(['relaxation strategy: ................................. ', self.KKT_complementarity_relaxation_strategy])
            end
            disp('2. Problem Size')
            disp(['number of decision variable (z): ..................... ', num2str(self.Dim.z)])
            disp(['number of parameter (p): ............................. ', num2str(self.Dim.p)])
            disp(['number of equality constraint (h): ................... ', num2str(self.Dim.h)])
            disp(['number of inequality constraint (c): ................. ', num2str(self.Dim.c)])
            
        end
    end
    
    %% Other method
    methods
        % create gap constraint based NLP
        nlp = create_primal_gap_NLP(self, OCPEC)      

        nlp = create_D_gap_NLP(self, OCPEC) 

        [d_func, d_grad, d_hessian] = create_strongly_convex_func(self, OCPEC) 
        
        omega_solver = create_omega_solver(self, OCPEC, param_c)

        phi_c_func = create_phi_c_func(self, OCPEC, param_c)        

        % create KKT based NLP
        nlp = create_KKT_based_NLP(self, OCPEC) 

        [KKT_stationarity_func, KKT_complementarity_func] = create_KKT_reformulation(self, OCPEC)
        
    end
    
end

