clear all
clc

%%
Data_D_gap_param = load('Data_test_D_gap_param');

%%
lineStyles = {'-','--',':','-.'};
figure(1)
for i = 1 : numel(Data_D_gap_param.rec.param_a)
    continuationStepNum_i = Data_D_gap_param.rec.Info{i}.continuationStepNum;
    semilogy(Data_D_gap_param.rec.Info{i}.Log.VI_nat_res(1 : continuationStepNum_i), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})
    hold on
end
grid on
legend(Data_D_gap_param.rec.name, 'Location','northeast', 'FontSize', 8)
xlabel('Continuation step', 'FontSize', 11)
ylabel('VI natural residual (max)', 'FontSize', 11)

figure(2)
for i = 1 : numel(Data_D_gap_param.rec.param_a)
    continuationStepNum_i = Data_D_gap_param.rec.Info{i}.continuationStepNum;
    semilogy(Data_D_gap_param.rec.Info{i}.Log.KKT_error(1 : continuationStepNum_i), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})
    hold on
end
grid on
legend(Data_D_gap_param.rec.name, 'Location','northeast', 'FontSize', 8)
xlabel('Continuation step', 'FontSize', 11)
ylabel('KKT error', 'FontSize', 11)

figure(3)
for i = 1 : numel(Data_D_gap_param.rec.param_a)
    continuationStepNum_i = Data_D_gap_param.rec.Info{i}.continuationStepNum;
    % stairs(Data_D_gap_param.rec.Info{i}.Log.time(2 : continuationStepNum_i), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})
    semilogy(Data_D_gap_param.rec.Info{i}.Log.time(2 : continuationStepNum_i), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})
    hold on
end
grid on
legend(Data_D_gap_param.rec.name, 'Location','northeast', 'FontSize', 8)
xlabel('Continuation step', 'FontSize', 11)
ylabel('Time Elapsed [s]', 'FontSize', 11)
