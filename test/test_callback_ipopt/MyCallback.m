classdef MyCallback < casadi.Callback
  properties
    d
  end
  methods
    function self = MyCallback(name, d)
      self@casadi.Callback();
      self.d = d;
      construct(self, name);
    end

    % Number of inputs and outputs
    function output = get_n_in(self)
      output = 1;
    end
    function output =get_n_out(self)
      output = 1;
    end

    % Initialize the object
    function init(self)
      disp('initializing object')
    end

    % Evaluate numerically
    function arg = eval(self, arg)
      x = arg{1};
      f = sin(self.d * x);
      arg = {f};
    end
  end
end