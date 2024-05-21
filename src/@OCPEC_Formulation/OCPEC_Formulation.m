classdef OCPEC_Formulation < handle
    %create an OCPEC with the form of:
    %  min  L_T(x) + int_0^T L_S(x, u, lambda) dt,
    %  s.t. Dot{x} = f(x, u, lambda)
    %       lambda \in SOL(K, F(x, u, lambda))
    %       K = {lambda | g(lambda) >= 0}
    %       G(x, u) >= 0,
    %       C(x, u) = 0,
    
    properties
        timeHorizon % time horizon
        nStages % number of discretized stage
        timeStep % discretization time step
        
        x0 % initial state

        x % differentiable state
        u % control input
        lambda % algebraic variable  
        
        L_T % terminal cost
        L_S % stage cost   
        
        f % ODE r.h.s function
        g % convex inequality formulating the VI set K
        F % VI function        
        VISetType char {mustBeMember(VISetType,{...
            'finitely_representable',...
            'polyhedral',...
            'box_constraint',...
            'nonnegative_orthant'...
            })} = 'box_constraint' % type of VI set K    
        bl % lower bound of box-constraint VI set K (used to construct the explicit expression of gap functions)
        bu % upper bound of box-constraint VI set K (used to construct the explicit expression of gap functions)
        
        G % path inequality constraint
        C % path equality constraint 
        
        Dim % variable dimemsion record
        FuncObj % CasADi function object           
    end
    
    %% Constructor method     
    methods
        function self = OCPEC_Formulation(...
                timeHorizon, nStages, timeStep,...
                x0, ...
                x, u, lambda,...
                L_T, L_S,...
                f, g, F, VISetType, bl, bu,...
                G, C)
            %OCPEC_Formulation: Construct an instance of this class
            %   Detailed explanation goes here
            disp('creating OCPEC...')
            %% formulate OCPEC
            % time parameter
            self.timeHorizon = timeHorizon;
            self.nStages = nStages;
            self.timeStep = timeStep;
            % initial state
            self.x0 = x0;            
            % variable
            self.x = x;
            self.u = u;
            self.lambda = lambda;             
            % cost function
            self.L_T = L_T;
            self.L_S = L_S;
            % DVI
            self.f = f;
            self.g = g;
            self.F = F;             
            self.VISetType = VISetType;
            self.bl = bl;
            self.bu = bu;
            % inequality and equality path constraint
            self.G = G;
            self.C = C;
            % dim record
            self.Dim = struct(...
                'x', size(x, 1), 'u', size(u, 1), 'lambda', size(lambda, 1),...
                'g', size(g, 1), 'G', size(G, 1), 'C', size(C, 1));             
            % function object
            self.FuncObj = self.create_FuncObj(); 

            %% display OCPEC informulation
            disp('*---------------------------------- OCPEC Information -----------------------------------*')
            disp('1. time parameter')
            disp(['time horizon: ............................... ', num2str(self.timeHorizon)])
            disp(['discretization stage: ....................... ', num2str(self.nStages)])
            disp(['time step: .................................. ', num2str(self.timeStep)])            
            disp('2. problem size')
            disp(['VI set type: ................................ ', self.VISetType])
            disp(['number of state variable (x): ............... ', num2str(self.Dim.x)])
            disp(['number of control variable (u): ............. ', num2str(self.Dim.u)])
            disp(['number of algebraic variable (lambda): ...... ', num2str(self.Dim.lambda)])
            disp(['number of VI set inequality (g): ............ ', num2str(self.Dim.g)])
            disp(['number of path inequality constraint (G): ... ', num2str(self.Dim.G)])
            disp(['number of path equality constraint (C): ..... ', num2str(self.Dim.C)])

            disp('Done!')

        end

    end
    %% Other method
    methods
        FuncObj = create_FuncObj(self)
    end
        
end

