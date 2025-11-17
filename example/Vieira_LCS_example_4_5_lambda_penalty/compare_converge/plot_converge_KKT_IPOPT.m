clear all
clc

%%
Data_converge_KKT_IPOPT = load('Data_test_converge_KKT_IPOPT');
rec = Data_converge_KKT_IPOPT.rec;
method_name = rec.name;
lineStyles = {'-','--',':','-.'};

%% VI nature residual, KKT error, time/step v.s. continuation step
l_Max = 25;
N = 2000;
% VI nature residual
figure(1)
for i = 1 : numel(method_name)
    Info_i = rec.Info{i};
    semilogy(Info_i.Log.VI_nat_res(1 : l_Max + 1), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})      
    hold on
end
grid on
xlabel('Continuation step $l$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$ E_{VI} $', 'Interpreter','latex', 'FontSize', 11)
legend(method_name, 'Location','northeast', 'FontSize', 11)
xlim([0, l_Max])

% KKT error
figure(2)
subplot(2, 1, 1)
for i = 1 : numel(method_name) 
    Info_i = rec.Info{i};
    semilogy(Info_i.Log.KKT_error(1 : l_Max + 1)/N, 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})
    hold on
end
grid on
ylabel('KKT error', 'Interpreter','latex', 'FontSize', 11)
legend(method_name, 'Location','northeast', 'FontSize', 9)
xlim([0, l_Max])
ylim([1e-15, 1e5])

subplot(2, 1, 2)
for i = 1 : numel(method_name) 
    Info_i = rec.Info{i};
    % semilogy(Info_i.Log.time(2 : l_Max + 1), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})
    stairs(Info_i.Log.time(2 : l_Max + 1), 'LineWidth', 1,  'LineStyle', lineStyles{mod(i,4) + 1})
    hold on
end
grid on
xlabel('Continuation step $l$', 'Interpreter','latex', 'FontSize', 11)
ylabel('Computation time [s]', 'Interpreter','latex', 'FontSize', 11)
xlim([0, l_Max])

