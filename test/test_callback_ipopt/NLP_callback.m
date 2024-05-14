classdef NLP_callback < casadi.Callback
    % implement the python code of the developer answer in https://groups.google.com/g/casadi-users/c/nboKTSx2RRw
    properties        
        SQ_FuncObj
    end
    %% constructor method
    methods
        function self = NLP_callback(name, SQ_FuncObj, opts)
            self@casadi.Callback();      
            self.SQ_FuncObj = SQ_FuncObj;
            self.construct(name, opts);
        end
    end

    %% other method
    methods
        % Number of inputs and outputs
        function output = get_n_in(self)
            output = 2;
        end
        function output = get_n_out(self)
            output = 1;
        end
        % Initialize the object
        function init(self)
            disp('initializing object')
        end
        % finalize the object
        function finalize(self)
            disp('Finalizing object construction')
        end
        % Set sparsity of input
        function out = get_sparsity_in(self, i)
            if i==0 % Note: index-0 based
                out = casadi.Sparsity.dense(2,1);
            elseif i == 1
                out = casadi.Sparsity.dense(2,1);
            end
        end    
        % Set sparsity of output
        function out = get_sparsity_out(self, i)
            if i==0 % Note: index-0 based
                out = casadi.Sparsity.dense(1,1);
            end
        end         

        % Evaluate numerically
        function [output] = eval(self, arg)
            x = full(arg{1});
            y = full(arg{2});
            f = self.SQ_FuncObj(x, y);
            % disp(x)
            output = {f};
        end

        % jacobian
        % function output = has_jacobian(self)
        %     output = true;
        % end
        % function output = get_jacobian(self, name, options)
        %     output = casadi.Function(name, {self.X}, {jacobian(self.SQ, self.X)}, options);
        % end

    end

end