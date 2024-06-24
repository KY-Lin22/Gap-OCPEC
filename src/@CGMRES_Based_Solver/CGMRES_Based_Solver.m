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
            self.FuncObj = self.create_FuncObj(NLP); 
            disp('Done!')
        end
    end

    %% Other method
    methods
        % create function Object
        FuncObj = create_FuncObj(self, NLP) 

        % main function of solver combining non-interior-point and CGMRES method
        [z_Opt, Info] = solve_NLP(self, z_Init, s_Init, s_End)

        % evaluate natural residual
        natRes = evaluate_natural_residual(self, z_Opt)

        %% stage 1: Non-Interior-Point Method
        % main function of non-interior-point method
        [z_Opt, Info] = non_interior_point_method(self, z_Init, p)        

        % evaluate KKT matrix 
        KKT_Matrix = evaluate_KKT_Matrix(self, h_grad, c_grad, LAG_hessian, PSI_grad_c, PSI_grad_gamma_c)

        % evaluate KKT error
        [KKT_error_primal, KKT_error_dual, KKT_error_complementary] = ...
            evaluate_KKT_error(self, gamma_h, gamma_c, h, c, LAG_grad_z)

        % merit line search
        [z_k, gamma_h_k, gamma_c_k, Info] = LineSearch_Merit(self,...
            beta, s, sigma, ...
            z, gamma_h, gamma_c, dz, dgamma_h, dgamma_c, ...
            J, h, PSI, J_grad)

        %% stage 2: CGMRES Method
        % evaluate KKT residual
        KKT_Residual = evaluate_KKT_Residual(self, Y, p)

        % evaluate sensitivity matrix
        sensitivity_Matrix = evaluate_sensitivity_Matrix(self, Y, p);

        % CGMRES method
        [Y_dot, Info] = CGMRES_method(self, Y, p, p_dot, Y_dot_Init, h_FD, k_max, epsilon)

        % directly solve differential equation
        [Y_dot, Info] = solve_differential_equation(self, Y, p, p_dot, epsilon)

    end

    methods(Static)
        % create solver option
        Option = create_Option()
    end

end