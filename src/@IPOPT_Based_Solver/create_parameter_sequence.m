function [S, l_Max] = create_parameter_sequence(self, s_Init, s_End)
%UNTITLED30 Summary of this function goes here
%   Detailed explanation goes here

% load parameter
kappa_s_times = self.Option.Continuation.kappa_s_times;
kappa_s_exp = self.Option.Continuation.kappa_s_exp;
% evaluate parameter sequence and the max number of continuation step
s_l = s_Init;
S = s_l;
l_Max = 0;
while true
    if s_l == s_End
        break
    else        
        s_trial = min([kappa_s_times*s_l, s_l^kappa_s_exp]);
        s_l = max([s_trial, s_End], [], 2);
        S = [S, s_l];
        l_Max = l_Max + 1;
    end
end

end