classdef MyGamma < casadi.Callback
  methods
    function self = MyGamma(name, opts)
      self@casadi.Callback();
      self.construct(name, opts);
    end
    function [results] = eval(self, arg)
      % A cell array of DMs comes in
      arg = full(arg{1});
      % Some code that is not supported by CasADi
      res = gamma(arg);     
      % A cell array of DMs should go out
      results = {res};
    end
  end
end