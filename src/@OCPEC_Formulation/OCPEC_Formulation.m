classdef OCPEC_Formulation < handle
    %create an OCPEC with the form of:
    %  min  L_T(x) + int_0^T L_S(x, u, lambda) dt,
    %  s.t. Dot{x} = f(x, u, lambda)
    %       lambda \in SOL(K, F(x, u, lambda))
    %       K = {lambda | g(lambda) >= 0}
    %       G(x, u) >= 0,
    %       C(x, u) = 0,
    
    properties
        TimeHorizon % time horizon
        nStages % number of discretized stage
        timeStep % discretization time step
        
        x0 % initial state
        xRef % reference state     
        
        x % differentiable state
        u % control input
        lambda % algebraic variable  

        xMax % x upper bound
        xMin % x lower bound
        uMax % u upper bound
        uMin % u lower bound
        lambdaMax % lambda upper bound
        lambdaMin % lambda lower bound
        
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
        
        G % path inequality constraint (only including bound constraint for x and u, other types of constraints should be transferred into equality constraint)
        C % path equality constraint 
        
        Dim % variable dimemsion record
        FuncObj % CasADi function object           
    end
    
    %% Constructor method     
    methods
        function self = OCPEC_Formulation(...
                TimeHorizon, nStages, timeStep,...
                x0, xRef,...
                x, u, lambda,...
                xMax, xMin, uMax, uMin, lambdaMax, lambdaMin,...
                L_T, L_S,...
                f, g, F, VISetType,...
                G, C)
            %OCPEC_Formulation: Construct an instance of this class
            %   Detailed explanation goes here
            % time parameter
            self.TimeHorizon = TimeHorizon;
            self.nStages = nStages;
            self.timeStep = timeStep;
            % initial and reference state
            self.x0 = x0;
            self.xRef = xRef;             
            % variable and their bounds
            self.x = x;
            self.u = u;
            self.lambda = lambda;   
            self.xMax = xMax;
            self.xMin = xMin;
            self.uMax = uMax;
            self.uMin = uMin;
            self.lambdaMax = lambdaMax;
            self.lambdaMin = lambdaMin;            
            % cost function
            self.L_T = L_T;
            self.L_S = L_S;
            % DVI
            self.f = f;
            self.g = g;
            self.F = F;             
            self.VISetType = VISetType;
            if strcmp(self.VISetType, 'box_constraint')
                self.bl = lambdaMin;
                self.bu = lambdaMax;
            end
            % inequality and equality path constraint
            self.G = G;
            self.C = C;
            % dim record
            self.Dim = struct(...
                'x', size(x, 1), 'u', size(u, 1), 'lambda', size(lambda, 1),...
                'g', size(g, 1), 'G', size(G, 1), 'C', size(C, 1));             
            % function object
            self.FuncObj = self.createFuncObj();                        
        end

    end
    %% Other method
    methods
        FuncObj = createFuncObj(self)
    end
        
end
