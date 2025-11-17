clear all
clc

%%
Data_converge_Gap_DynSys = load('Data_test_converge_Gap_DynSys');
method_name = Data_converge_Gap_DynSys.rec.name;
OCPEC = OCPEC_Vieira_LCS_analytic();
OCPEC.timeHorizon = Data_converge_Gap_DynSys.rec.OCPEC.timeHorizon;
OCPEC.nStages = Data_converge_Gap_DynSys.rec.OCPEC.nStages;
OCPEC.timeStep = Data_converge_Gap_DynSys.rec.OCPEC.timeStep;
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
