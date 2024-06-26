function [P, P_dot, l_Max] = create_parameter_sequence(self, s_Init, s_End)
%UNTITLED27 Summary of this function goes here
%   Detailed explanation goes here

% load option parameter 
kappa_s_times = self.Option.Continuation.kappa_s_times;
sigma_Init = self.Option.Continuation.sigma_Init;
sigma_End = self.Option.Continuation.sigma_End;
kappa_sigma_times = self.Option.Continuation.kappa_sigma_times;
kappa_sigma_exp = self.Option.Continuation.kappa_sigma_exp;
dtau = self.Option.Continuation.dtau;
% formulate init and end parameter vector
p_Init = [s_Init; sigma_Init];
p_End = [s_End; sigma_End];
% evaluate parameter sequence and the max number of continuation step
p_l = p_Init;
P = p_l;
l_Max = 0;
while true
    if all(p_l == p_End)
        break
    else
        p_trial = [kappa_s_times*p_l(1);...
            min([kappa_sigma_times*p_l(2), p_l(2)^kappa_sigma_exp])];
        p_l = max([p_trial, p_End], [], 2);
        P = [P, p_l];
        l_Max = l_Max + 1;
    end
end
% evaluate parameter time derivative by finite difference
P_dot = zeros(size(P));
if l_Max ~= 0
    P_dot(:, 1 : end - 1) = (P(:, 2 : end) - P(:, 1 : end - 1))/dtau;
    P_dot(:, end) = P_dot(:, end - 1);
end

end