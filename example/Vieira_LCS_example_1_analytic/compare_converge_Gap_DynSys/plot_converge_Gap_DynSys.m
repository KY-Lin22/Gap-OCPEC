clear all
clc

%%
timeHorizon = 1;
nStages = 100;
OCPEC = OCPEC_Vieira_LCS_analytic();
OCPEC.timeHorizon = timeHorizon;
OCPEC.nStages = nStages;
OCPEC.timeStep = OCPEC.timeHorizon ./ OCPEC.nStages;

Data_converge_Gap_DynSys = load('Data_test_converge_Gap_DynSys');
method_name = Data_converge_Gap_DynSys.rec.name;
lineStyles = {'-','--',':','-.'};

%% VI nature residual, KKT error, time/step v.s. continuation step
l_Max = 500;
% VI nature residual
figure(1)
for i = 1 : numel(method_name)
    Info_i = Data_converge_Gap_DynSys.rec.Info{i};
    semilogy(Info_i.Log.VI_nat_res(1 : l_Max + 1), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})      
    hold on
end
grid on
xlabel('Continuation step $l$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$ E_{VI} $', 'Interpreter','latex', 'FontSize', 11)
legend(method_name, 'Location','northeast', 'FontSize', 8)
xlim([0, l_Max])

% KKT error
figure(2)
for i = 1 : numel(method_name) 
    Info_i = Data_converge_Gap_DynSys.rec.Info{i};
    semilogy(Info_i.Log.KKT_error(1 : l_Max + 1), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})
    hold on
end
grid on
xlabel('Continuation step $l$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$ E_{KKT} $', 'Interpreter','latex', 'FontSize', 11)
legend(method_name, 'Location','northeast', 'FontSize', 8)
xlim([0, l_Max])

% KKT residual
figure(3)
for i = 1 : numel(method_name) 
    Info_i = Data_converge_Gap_DynSys.rec.Info{i};
    semilogy(Info_i.Log.KKT_res(1 : l_Max + 1), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})
    hold on
end
grid on
xlabel('Continuation step $l$', 'Interpreter','latex', 'FontSize', 11)
ylabel('KKT residual $ \| T \|_2 / N $', 'Interpreter','latex', 'FontSize', 11)
legend(method_name, 'Location','northeast', 'FontSize', 8)
xlim([0, l_Max])

% time
figure(4)
for i = 1 : numel(method_name) 
    Info_i = Data_converge_Gap_DynSys.rec.Info{i};
    stairs(Info_i.Log.time(2 : l_Max + 1), 'LineWidth', 1,  'LineStyle', lineStyles{mod(i,4) + 1})
    hold on
end
grid on
xlabel('Continuation step $l$', 'Interpreter','latex', 'FontSize', 11)
ylabel('Time [s]', 'Interpreter','latex', 'FontSize', 11)
legend(method_name, 'Location','northeast', 'FontSize', 8)

%% compare with KKT
% figure(5)
% time_mat = zeros(l_Max, numel(method_name));
% for i = 1 : numel(method_name)
%     Info_i = Data_converge_Gap_DynSys.rec.Info{i};
%     time_mat(:, i) = Info_i.Log.time(2 : l_Max + 1);
%     disp(['name: ' method_name{i}, '; ', 'avetime: ', num2str(sum(Info_i.Log.time(2 : l_Max + 1))/l_Max, '%10.4f'), 's' ])
% end
% boxplot(time_mat, 'Labels', method_name)

%% analytic optimal trajectory
timeAxis = 0 : OCPEC.timeStep : OCPEC.nStages * OCPEC.timeStep;
F_FuncObj_map = OCPEC.FuncObj.F.map(OCPEC.nStages);
[X_analytic_Opt, U_analytic_Opt, LAMBDA_analytic_Opt] = compute_analytic_optimal_trajectory(OCPEC);
F_analytic_Opt = full(F_FuncObj_map(X_analytic_Opt, U_analytic_Opt, LAMBDA_analytic_Opt));

figure(6)
subplot(4,1,1)
plot(timeAxis, [OCPEC.x0, X_analytic_Opt], 'r', 'LineWidth',1.5)
legend('$ x^{*}_{ana}(t)$', 'Interpreter','latex', 'FontSize', 11)
subplot(4,1,2)
stairs(timeAxis(2:end), U_analytic_Opt, 'b',  'LineWidth', 1.5)
legend('$ u^{*}_{ana}(t)$', 'Interpreter','latex', 'FontSize', 11)
subplot(4,1,3)
plot(timeAxis(2:end), LAMBDA_analytic_Opt, 'k',  'LineWidth', 1.5)
legend('$ \lambda^{*}_{ana}(t)$', 'Interpreter','latex', 'FontSize', 11)
subplot(4,1,4)
plot(timeAxis(2:end), F_analytic_Opt, 'g',  'LineWidth', 1.5)
legend('$ F^{*}_{ana}(t)$', 'Interpreter','latex', 'FontSize', 11)
xlabel('Time [s]', 'Interpreter','latex', 'FontSize', 11)

%% error between obtained numerical optimal trajectory and analytic optimal trajectory 
z_Node = cumsum([OCPEC.Dim.x, OCPEC.Dim.u, OCPEC.Dim.lambda, OCPEC.Dim.lambda]);
figure(7)
for i = 1 : numel(method_name)
    z_Opt_i = Data_converge_Gap_DynSys.rec.z_Opt{i};
    Z_Opt_i = reshape(z_Opt_i, [], OCPEC.nStages);
    
    X_Opt_i      = Z_Opt_i(            1 : z_Node(1), :);
    U_Opt_i      = Z_Opt_i(z_Node(1) + 1 : z_Node(2), :);    
    LAMBDA_Opt_i = Z_Opt_i(z_Node(2) + 1 : z_Node(3), :);
    F_Opt_i = full(F_FuncObj_map(X_Opt_i, U_Opt_i, LAMBDA_Opt_i));

    subplot(3,1,1)
    plot(timeAxis(2:end), abs(X_Opt_i - X_analytic_Opt), 'LineWidth',1.5, 'LineStyle', lineStyles{mod(i,4) + 1})
    title('$ |x^{*}_{ana}(t) - x^{*}_{num}(t) |$', 'Interpreter','latex', 'FontSize', 11)
    hold on

    subplot(3,1,2)
    stairs(timeAxis(2:end), abs(U_Opt_i - U_analytic_Opt), 'LineWidth', 1.5, 'LineStyle', lineStyles{mod(i,4) + 1})
    title('$ |u^{*}_{ana}(t) - u^{*}_{num}(t) |$', 'Interpreter','latex', 'FontSize', 11)
    hold on

    subplot(3,1,3)
    plot(timeAxis(2:end), abs(LAMBDA_Opt_i - LAMBDA_analytic_Opt), 'LineWidth', 1.5, 'LineStyle', lineStyles{mod(i,4) + 1})
    title('$ |\lambda^{*}_{ana}(t) - \lambda^{*}_{num}(t) |$', 'Interpreter','latex', 'FontSize', 11)
    hold on

    % subplot(4,1,4)
    % plot(timeAxis(2:end), abs(F_Opt_i - F_analytic_Opt), 'LineWidth', 1.5, 'LineStyle', lineStyles{mod(i,4) + 1})   
    % title('$ |F^{*}_{ana}(t) - F^{*}_{num}(t) |$', 'Interpreter','latex', 'FontSize', 11)
    % xlabel('time [s]')
    % hold on
end
xlabel('Time [s]', 'Interpreter','latex', 'FontSize', 11)
legend(method_name, 'Location','northeast', 'FontSize', 8)
