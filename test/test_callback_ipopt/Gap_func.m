classdef Gap_func < casadi.Callback
    %% constructor method
    properties
        phi_c_func
    end
    methods
        function self = Gap_func(name, phi_c, opts)
            self@casadi.Callback();
            self.phi_c_func = phi_c;
            self.construct(name, opts);
        end
    end

    %% other methods
    methods
        % Number of inputs and outputs
        function v = get_n_in(self)
            v = 2;
        end
        function v = get_n_out(self)
            v = 1;
        end
        % Set sparsity of input
        function out = get_sparsity_in(self, i)
            if i==0 % Note: index-0 based
                out = casadi.Sparsity.dense(1,1);
            elseif i == 1
                out = casadi.Sparsity.dense(1,1);
            end
        end
         % Set sparsity of output
        function out = get_sparsity_out(self, i)
            if i==0 % Note: index-0 based
                out = casadi.Sparsity.dense(1,1);
            end
        end       
        % Initialize the object
        function init(self)
            disp('initializing object')
        end
        % finalize the object
        function finalize(self)
            disp('Finalizing object construction')
        end        
        % Evaluate numerically
        function output = eval(self, arg)
            lambda = full(arg{1});
            eta = full(arg{2});
            phi_c = self.phi_c_func(lambda, eta);
            output = {phi_c};
        end
    end


end