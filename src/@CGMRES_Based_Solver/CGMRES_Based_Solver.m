classdef CGMRES_Based_Solver < handle
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
        [z_Opt, Info] = solve_NLP(self, z_Init, s_Init, s_End)

        % create parameter and its time derivative sequence
        [P, P_dot, l_Max] = create_parameter_sequence(self, s_Init, s_End);

        % evaluate natural residual
        natRes = evaluate_natural_residual(self, z_Opt)

        %% stage 1: evaluate first iterate by solving first parameterized NLP
        % main function of solving first parameterized NLP
        [Y, Info] = solve_first_NLP(self, z_Init, p)

        % non-interior-point method
        [z_Opt, Info] = non_interior_point_method(self, z_Init, p) 

        % merit line search in non-interior-point method
        [Y_k, Info] = LineSearch_Merit(self, beta, Y, p, dY);

        %% stage 2: evaluate new iterate by integrating a differential equation
        % evaluate new iterate by integration
        [Y_l, Info] = integrate_differential_equation(self, Y, Y_dot)

        % solving differential equation
        [Y_dot, Info] = solve_differential_equation(self, Y, p, p_dot, Y_dot_Init, epsilon)

        % FDGMRES method



    end

    methods(Static)
        % create solver option
        Option = create_Option()
    end

end