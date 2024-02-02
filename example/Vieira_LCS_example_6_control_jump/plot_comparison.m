clear all
clc

%%
Data_primal_gap = load('Data_primal_gap.mat');
Data_primal_gap_Rec = Data_primal_gap.Rec;

Data_D_gap = load('Data_D_gap.mat');
Data_D_gap_Rec = Data_D_gap.Rec;

Data_KKT = load('Data_KKT.mat');
Data_KKT_Rec = Data_KKT.Rec;

figure(1)
semilogx(Data_primal_gap_Rec.s_End(1 : end), Data_primal_gap_Rec.cost(1 : end), 'b+',......
    Data_D_gap_Rec.s_End(1 : end), Data_D_gap_Rec.cost(1 : end), 'r*',...
    Data_KKT_Rec.s_End(1 : end), Data_KKT_Rec.cost(1 : end), 'go')
grid on
legend('primal gap', 'D gap', 'KKT')
xlabel('relaxation parameter')
ylabel('cost')

figure(2)
semilogx(Data_primal_gap_Rec.s_End(1 : end), Data_primal_gap_Rec.nat_res(1 : end), 'b+',......
    Data_D_gap_Rec.s_End(1 : end), Data_D_gap_Rec.nat_res(1 : end), 'r*',...
    Data_KKT_Rec.s_End(1 : end), Data_KKT_Rec.nat_res(1 : end), 'go')
grid on
legend('primal gap', 'D gap', 'KKT')
xlabel('relaxation parameter')
ylabel('equilibrium constraint violation')


figure(3)
semilogx(Data_primal_gap_Rec.s_End(1 : end), Data_primal_gap_Rec.time(1 : end), 'b+',......
    Data_D_gap_Rec.s_End(1 : end), Data_D_gap_Rec.time(1 : end), 'r*',...
    Data_KKT_Rec.s_End(1 : end), Data_KKT_Rec.time(1 : end), 'go')
grid on
legend('primal gap', 'D gap', 'KKT')
xlabel('relaxation parameter')
ylabel('time [s]')