function [S, l_Max] = create_parameter_sequence(self)
%UNTITLED30 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
% load parameter
s_Init = self.Option.Continuation.s_Init;
s_End = self.Option.Continuation.s_End;
% evaluate parameter sequence and the max number of continuation step
s_l = s_Init;
switch self.Option.Continuation.update_rule
    case 'times_exp'
        S = s_l;
        l_Max = 0;
        kappa_s_times = self.Option.Continuation.kappa_s_times;
        kappa_s_exp = self.Option.Continuation.kappa_s_exp;
        while true
            if s_l == s_End
                break
            else
                s_trial = min([kappa_s_times*s_l, s_l^kappa_s_exp]);
                s_l = max([s_trial, s_End]);
                S = [S, s_l];
                l_Max = l_Max + 1;
            end
        end
    case 'dynamics'
        epsilon_s = self.Option.Continuation.epsilon_s;
        dtau = self.Option.Continuation.dtau;
        l_Max = self.Option.Continuation.l_Max;
        S = zeros(1, l_Max + 1);
        S(1) = s_l;
        % create dynamics
        s_SX = SX.sym('s', 1, 1);
        s_dot_SX = -epsilon_s*(s_SX - s_End);
        s_dot_FuncObj = Function('s_dot', {s_SX}, {s_dot_SX}, {'s'}, {'s_dot'});
        for l = 1 : l_Max
            % compute next parameter by integrating dynamics using RK4 method
            k_1 = full(s_dot_FuncObj(s_l));
            k_2 = full(s_dot_FuncObj(s_l + (dtau/2)*k_1));
            k_3 = full(s_dot_FuncObj(s_l + (dtau/2)*k_2));
            k_4 = full(s_dot_FuncObj(s_l + dtau*k_3));
            s_l = s_l + (dtau/6) * (k_1 + 2*k_2 + 2*k_3 + k_4);
            S(l + 1) = s_l;
        end
end

end