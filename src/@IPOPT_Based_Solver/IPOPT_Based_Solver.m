classdef IPOPT_Based_Solver < handle
    %Implementation of IPOPT with continuation method which solves the following NLP problem:
    %  min  J(z, p),
    %  s.t. h(z, p) = 0,
    %       c(z, p) >= 0
    % where: z is the variable, 
    %        p is the parameter, 
    %        J is the cost, and h, c are the constraint arranged in a stagewise manner.
    %
    properties
        OCPEC % struct, optimal control problem with equalibrium constraints
        NLP % struct, nonlinear programming problem (discretized OCPEC)
        Option % struct, IPOPT solver option
        Solver % function object, IPOPT solver
    end
    
    %% Constructor method       
    methods
        function self = IPOPT_Based_Solver(OCPEC, NLP)
            %IPOPT_Based_Solver Construct an instance of this class
            %   Detailed explanation goes here
            import casadi.*
            disp('creating solver...')
            % properties: OCPEC and NLP
            self.OCPEC = OCPEC;
            self.NLP = NLP;  
            
            % properties: solver option
            self.Option = self.create_Option();
            
            % properties: solver
            Prob = struct('x', NLP.z, 'f', NLP.J, 'g', [NLP.h; NLP.c], 'p', NLP.p);
            
            NLP_Solver_Option = self.Option.NLP_Solver;
            self.Solver = nlpsol('Solver', 'ipopt', Prob, NLP_Solver_Option);

            disp('Done!')
        end
        
    end
    
    %% Other method
    methods
        % create solver option
        Option = create_Option(self)   

        % create initial guess
        z_Init = create_initial_guess(self)
        
        % solve a sequence of NLP in a homotopy manner from p_Init to p_End 
        [z_Opt, Info] = solve_NLP(self, z_Init, p_Init, p_End)

        % evaluate natural residual
        natRes = evaluate_natural_residual(self, z_Opt)

        % show result (TO DO)
        show_result(self, Info)        
    end
    
end

