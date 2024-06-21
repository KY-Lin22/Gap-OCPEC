function [Y_dot, Info] = CGMRES_method(self, Y, p, p_dot, Y_dot_Init, h_FD, k_max, epsilon)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
timeStart_CGMRES = tic;
%% step 1: initialization
% KKT residual
T = self.evaluate_KKT_Residual(Y, p);
% forward difference approximation for T_grad_p * p_dot
DhT_Y_p_0_pDot = (self.evaluate_KKT_Residual(Y, p + h_FD * p_dot) - T)/h_FD;
% compute b term
b = -epsilon * T - DhT_Y_p_0_pDot;
% forward difference approximation for T_grad_Y * Y_dot 
DhT_Y_pNext_YDotInit_0 = (self.evaluate_KKT_Residual(Y + h_FD * Y_dot_Init, p + h_FD * p_dot) ...
    - self.evaluate_KKT_Residual(Y, p + h_FD * p_dot))/h_FD;
% r_hat, v_1, rho, beta, k
r_hat = b - DhT_Y_pNext_YDotInit_0;
v_1 = r_hat/(norm(r_hat, 2));
rho = norm(r_hat, 2);
beta = rho;
k = 0;
disp(['k = ', num2str(k, '%10.1d'), '; ', 'rho = ', num2str(rho, '%10.2e')])
%% step 2: compute V_k and y_k
Y_Dim = self.NLP.Dim.z + self.NLP.Dim.h + self.NLP.Dim.c;
V_k = zeros(Y_Dim, k_max);
V_k(:, 1) = v_1;
while k < k_max
    % (a) update counter k
    k = k + 1;
    % (b) 
    % initialize v_kNext
    v_k = V_k(:, k);
    DhT_Y_pNext_vk_0 = (self.evaluate_KKT_Residual(Y + h_FD * v_k, p + h_FD * p_dot) ...
        - self.evaluate_KKT_Residual(Y, p + h_FD * p_dot))/h_FD;
    v_kNext = DhT_Y_pNext_vk_0;
    % compute H_k
    H_k = zeros(k + 1, k);
    for j = 1 : k
        v_j = V_k(:, j);
        h_jk = v_kNext' * v_j;
        H_k(j, k) = h_jk;
        v_kNext = v_kNext - h_jk * v_j;
    end
    % (c)
    h_kNext_k = norm(v_kNext, 2);
    H_k(k + 1, k) = h_kNext_k;
    % (d)
    v_kNext = v_kNext/(norm(v_kNext, 2));
    V_k(:, k + 1) = v_kNext;
    % (e)
    e_1 = [1; zeros(k, 1)];
    y_k = lsqlin(-H_k, beta * e_1);
    rho = norm(beta * e_1 - H_k * y_k, 2);
    disp(['k = ', num2str(k, '%10.1d'), '; ', 'rho = ', num2str(rho, '%10.2e')])
end

%% step 3: compute Y_dot
Y_dot = Y_dot_Init + V_k(:, 1 : k_max) * y_k;

% output
timeElasped_CGMRES = toc(timeStart_CGMRES);
Info.rho = rho;
Info.time = timeElasped_CGMRES;
end