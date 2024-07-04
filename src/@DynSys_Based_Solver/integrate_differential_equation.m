function [Y_l, Info] = integrate_differential_equation(self, Y, Y_dot, p, p_dot, p_l, p_dot_l)
%UNTITLED29 Summary of this function goes here
%   approximate Y_l at time tau_l based on known quantities at time tau and tau_l
%  tau: previous time
%  tau_l: current time
%  tau_l - tau = dtau
timeStart = tic;

dtau = self.Option.Continuation.dtau;
if norm(Y_dot, inf) > 1e10
    % set as fail case: extremely large Y_dot may lead to divergence
    Y_l = Y;
    terminal_status = 0;
    terminal_msg = ['- Solver fails: ', 'because Y_dot is extremely large, which may lead to divergence']; 
else
    % set as success case: integrate to obtain a new iterate
    switch self.Option.Continuation.integration_method
        case 'explitic_Euler'
            % update iterate
            Y_l = Y + dtau * Y_dot;
        case 'RK4'
            % approximate the value of p and p_dot at the middle of interval [tau, tau_l]
            p_m = 0.5 * (p + p_l);
            p_dot_m = 0.5 * (p_dot + p_dot_l);
            % compute k_1, ..., k_4
            k_1 = Y_dot;
            k_2 = full(self.FuncObj.Y_dot(Y + (dtau/2)*k_1, p_m, p_dot_m));
            k_3 = full(self.FuncObj.Y_dot(Y + (dtau/2)*k_2, p_m, p_dot_m));
            k_4 = full(self.FuncObj.Y_dot(Y + dtau*k_3,     p_l, p_dot_l));
            % update iterate
            Y_l = Y + (dtau/6) * (k_1 + 2*k_2 + 2*k_3 + k_4);
        otherwise
            error('specified integration method is not supported')
    end
    terminal_status = 1;
    terminal_msg = ['- Solver succeeds: ', 'because a new iterate found by integrating a differential equation'];
end

timeElapsed = toc(timeStart);
Info.time = timeElapsed;
Info.terminal_status = terminal_status;
Info.terminal_msg = terminal_msg;
end