function p_l = evaluate_new_parameter(self, p, dtau)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% compute next parameter by integrating ODE using RK4 method
k_1 = full(self.FuncObj.p_dot(p));
k_2 = full(self.FuncObj.p_dot(p + (dtau/2)*k_1));
k_3 = full(self.FuncObj.p_dot(p + (dtau/2)*k_2));
k_4 = full(self.FuncObj.p_dot(p + dtau*k_3));
p_l = p + (dtau/6) * (k_1 + 2*k_2 + 2*k_3 + k_4);

end