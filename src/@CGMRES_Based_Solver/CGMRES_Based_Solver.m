classdef CGMRES_Based_Solver
    %Implementation of solver that combines non-interior-point and C/GMRES method which solves the following NLP problem:
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
        Option % struct, C/GMRES solver option
        FuncObj % struct, function object
    end

    %% Constructor method
    methods
        function self = CGMRES_Based_Solver(OCPEC, NLP, Option)
            %CGMRES_Based_Solver Construct an instance of this class
            %   Detailed explanation goes here
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
        %% basic method
        % create function Object
        FuncObj = create_FuncObj(self) 

        % main function of solver that combines non-interior-point method and CGMRES method

        % evaluate natural residual
        natRes = evaluate_natural_residual(self, z_Opt)

        %% stage 1: Non-Interior-Point Method
        % main function of non-interior-point method
        [z_Opt, Info] = non_interior_point_method(self, z_Init, p) 

        % evaluate KKT error
        [KKT_error_primal, KKT_error_dual, KKT_error_complementary] = evaluate_KKT_error(self, Y, LAG_grad_T, h, c)

        % merit line search
        [Y_k, Info] = LineSearch_Merit(self, beta, Y, p, dY);

        %% stage 2: C/GMRES Method 
        

        %% backup method for C/GMRES Method (through the lens of dynamical system)


    end

    methods(Static)
        % create solver option
        Option = create_Option()
    end

end