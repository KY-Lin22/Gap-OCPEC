classdef IPOPT_Based_Solver < handle
    %Implementation of IPOPT with continuation method which solves the following NLP problem:
    %  min  J(z),
    %  s.t. h(z) = 0,
    %       c(z, s) >= 0
    % where: z is the variable, 
    %        s is the parameter, 
    %        J is the cost, and h, c are the constraint arranged in a stagewise manner.
    %
    properties
        OCPEC % struct, optimal control problem with equalibrium constraints
        NLP % struct, nonlinear programming problem (discretized OCPEC)
        Option % struct, solver option
        FuncObj % struct, function object
    end
    
    %% Constructor method       
    methods
        function self = IPOPT_Based_Solver(OCPEC, NLP, Option)
            %IPOPT_Based_Solver Construct an instance of this class
            %   Detailed explanation goes here
            import casadi.*
            disp('creating solver...')
            % initialize properties: OCPEC, NLP, Option, FuncObj
            self.OCPEC = OCPEC;
            self.NLP = NLP;           
            self.Option = Option; 
            self.FuncObj = self.create_FuncObj(); 
            disp('Done!')
        end
        
    end
    
    %% Other method
    methods
        % create function Object
        FuncObj = create_FuncObj(self) 

        % solve a sequence of parameterized NLP in a homotopy manner from s_Init to s_End 
        [z_Opt, Info] = solve_NLP(self, z_Init, s_Init, s_End)

        % create parameter sequence
        [S, l_Max] = create_parameter_sequence(self, s_Init, s_End)

        % evaluate natural residual
        natRes = evaluate_natural_residual(self, z_Opt)  
    end
    
    methods(Static)
        % create solver option
        Option = create_Option()  
    end
end

