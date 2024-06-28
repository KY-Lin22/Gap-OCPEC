clear all
clc

%%
Data_IPOPT_DynSys = load('Data_test_IPOPT_DynSys');

%%
lineStyles = {'-','--',':','-.'};
fast_index = 15;
figure(1)
for i = 1 : numel(Data_IPOPT_DynSys.rec.name) - fast_index
    continuationStepNum_i = Data_IPOPT_DynSys.rec.Info{fast_index + i}.continuationStepNum;
    semilogy(Data_IPOPT_DynSys.rec.Info{fast_index + i}.Log.VI_nat_res(1 : continuationStepNum_i + 1), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})
    hold on
end
grid on
legend(Data_IPOPT_DynSys.rec.name(fast_index + 1 : end), 'Location','northeast', 'FontSize', 8)
xlabel('Continuation step', 'FontSize', 11)
ylabel('VI natural residual (max)', 'FontSize', 11)

%%
figure(2)
for i = 1 : numel(Data_IPOPT_DynSys.rec.name) - fast_index
    continuationStepNum_i = Data_IPOPT_DynSys.rec.Info{fast_index + i}.continuationStepNum;
    semilogy(Data_IPOPT_DynSys.rec.Info{fast_index + i}.Log.KKT_error(1 : continuationStepNum_i + 1), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})
    hold on
end
grid on
legend(Data_IPOPT_DynSys.rec.name(fast_index + 1 : end), 'Location','northeast', 'FontSize', 8)
xlabel('Continuation step', 'FontSize', 11)
ylabel('KKT error', 'FontSize', 11)

%%
figure(3)
for i = 1 : numel(Data_IPOPT_DynSys.rec.name) - fast_index
    continuationStepNum_i = Data_IPOPT_DynSys.rec.Info{fast_index + i}.continuationStepNum;
    semilogy(Data_IPOPT_DynSys.rec.Info{fast_index + i}.Log.time(2 : continuationStepNum_i + 1), 'LineWidth', 2,  'LineStyle', lineStyles{mod(i,4) + 1})
    hold on
end
grid on
legend(Data_IPOPT_DynSys.rec.name(fast_index + 1 : end), 'Location','northeast', 'FontSize', 8)
xlabel('Continuation step', 'FontSize', 11)
ylabel('Time Elapsed [s]', 'FontSize', 11)

%% display time
disp('-----------------------------------------------------------------------------------------------------------------------------------------')
disp(' num | step | time(ave, 2:end)[s] | time(T)[s] | time(1)[s] | time(2:end)[s] |  method name ')
for i = 1 : numel(Data_IPOPT_DynSys.rec.name)
    time_i    = Data_IPOPT_DynSys.rec.Info{i}.time;
    timeLog_i = Data_IPOPT_DynSys.rec.Info{i}.Log.time;
    step_i    = Data_IPOPT_DynSys.rec.Info{i}.continuationStepNum;
    disp([num2str(i, '%10.4d'),  ' | ', ...        
         num2str(step_i, '%10.4d'), ' |       ', ...
         num2str((time_i - timeLog_i(1))/step_i, '%10.4f'), '        |   ', ...
         num2str(time_i, '%10.4f'), '   |   ',...
         num2str(timeLog_i(1), '%10.4f'), '   |     ',...
         num2str(time_i - timeLog_i(1), '%10.4f'), '     |   ',...
         Data_IPOPT_DynSys.rec.name{i}])
end
