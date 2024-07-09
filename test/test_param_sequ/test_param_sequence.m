clear all
clc
import casadi.*

% init and end parameter
s_Init = 1e0;
s_End = 1e-12;
sigma_Init = 1e-2;
sigma_End = 1e-6;

p_Init = [s_Init; sigma_Init];
p_End = [s_End; sigma_End];

% time parameter
dtau = 0.01;
l_Max = 500;
% parameter dynamical system
eplison_p = 5;
p_sym = SX.sym('p_sym', 2, 1);
p_dot_sym = -eplison_p*(p_sym - p_End);
p_dot_FuncObj = Function('p_dot', {p_sym}, {p_dot_sym}, {'p'}, {'p_dot'});

% evaluate parameter sequence 
P = zeros(2, l_Max);
P_dot = zeros(2, l_Max);

p = p_Init;
for l = 1 : l_Max
    p_dot = full(p_dot_FuncObj(p));
    P(:, l) = p;
    P_dot(:, l) = p_dot;
    p_l = p + dtau * p_dot;    
    p = p_l;
end

%%
figure(1)
semilogy(P(1, :), 'LineWidth', 2)
hold on
semilogy(P(2, :), 'LineWidth', 2)
legend('$ s $', '$ \sigma $', 'Interpreter','latex')

figure(2)
semilogy(P_dot(1, :), 'LineWidth', 2)
hold on
semilogy(P_dot(2, :), 'LineWidth', 2)
legend('$ \dot{s} $', '$ \dot{\sigma} $', 'Interpreter','latex')