function [Y_l, Info] = integrate_differential_equation(self, Y, Y_dot)
%UNTITLED29 Summary of this function goes here
%   Detailed explanation goes here
timeStart = tic;

dtau = self.Option.Continuation.dtau;
switch self.Option.Continuation.integration_method
    case 'explitic_Euler'
        Y_l = Y + dtau * Y_dot;
    case 'RK4'

    otherwise
        error('specified integration method is not supported')
end

timeElapsed = toc(timeStart);
Info.time = timeElapsed;
Info.terminal_status = 1;
Info.terminal_msg = ['- Solver succeeds: ', 'because a new iterate found by integrating a differential equation'];
end