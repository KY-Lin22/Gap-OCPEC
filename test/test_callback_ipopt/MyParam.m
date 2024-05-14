classdef MyParam < casadi.Callback
    properties
        f
    end
    %% constructor method
    methods       
        function self = MyParam(name, opts)
            self@casadi.Callback();
            self.construct(name, opts);
            self.f = 0;
        end
    end
    %% overrides other method
    methods
        % Number of inputs and outputs
        function output = get_n_in(self)
            output = 1;
        end
        function output = get_n_out(self)
            output = 1;
        end
        % Initialize the object
        function init(self)
            disp('initializing object')
        end
        % Evaluate numerically
        function output = eval(self)
            output = casadi.DM(self.f);
        end
        function set(self, value)
            self.f = Call_Toolbox(value);
        end
        function output = has_jacobian(self)
            output = true;
        end
        function J = get_jacobian(self, name, inames, onames, opts)
            x = casadi.SX.sym('x', 1, 1);
            dummy = casadi.SX.sym('dummy', casadi.Sparsity(1,1));
            J = casadi.Function(name, [x, dummy], [0], inames, onames, opts);
        end
        function output = nlp_callback(param, iter_x)
            param.set(iter_x);
            output = 0;
        end


    end
    


end