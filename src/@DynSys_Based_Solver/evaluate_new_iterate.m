function [Y_l, Info] = evaluate_new_iterate(self, Y, p, dtau)
%UNTITLED29 Summary of this function goes here
%   evaluate Y_l at time tau_l based on quantities Y, p at time tau
%  tau: previous time
%  tau_l: current time
%  tau_l - tau = dtau
timeStart = tic;

% evaluate time derivative Y_dot
Y_dot = full(self.FuncObj.Y_dot(Y, p));

%% evaluate new iterate Y_l
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
            % evaluate the value of p at the middle and end of interval [tau, tau_l]
            p_m = self.evaluate_new_parameter(p, dtau/2);
            p_l = self.evaluate_new_parameter(p, dtau);
            % compute k_1, ..., k_4
            k_1 = Y_dot;
            k_2 = full(self.FuncObj.Y_dot(Y + (dtau/2)*k_1, p_m));
            k_3 = full(self.FuncObj.Y_dot(Y + (dtau/2)*k_2, p_m));
            k_4 = full(self.FuncObj.Y_dot(Y + dtau*k_3,     p_l));
            % update iterate
            Y_l = Y + (dtau/6) * (k_1 + 2*k_2 + 2*k_3 + k_4);
        otherwise
            error('specified integration method is not supported')
    end
    terminal_status = 1;
    terminal_msg = ['- Solver succeeds: ', 'because a new iterate found by integrating a differential equation'];
end

%% Info
timeElapsed = toc(timeStart);

Info.time = timeElapsed;
Info.terminal_status = terminal_status;
Info.terminal_msg = terminal_msg;

end