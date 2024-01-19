clear all
clc

%%
Data_Scholtes = load('Data_Scholtes');
Data_primal_gap = load('Data_primal_gap');
Data_D_gap = load('Data_D_gap');
% Data_Lin_Fukushima = load('Data_Lin_Fukushima');
% Data_Kadrani = load('Data_Kadrani');

x_Init = Data_Scholtes.Rec.x_Init;

%% Scholtes
stationarity_type_Scholtes = Data_Scholtes.Rec.x_Opt_type;
color_Scholtes = Data_Scholtes.Rec.x_Opt_color;
figure(1)
for i = 1 : size(x_Init, 1)
    for j = 1 : size(x_Init, 2)
        x_Init_ij = x_Init{i, j};
        plot(x_Init_ij(1), x_Init_ij(2), stationarity_type_Scholtes{i, j}, 'MarkerSize', 10, 'MarkerEdgeColor', color_Scholtes{i, j})
        hold on
    end    
end
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
axis([-1.5, 2.5, -1.5, 2.5])

%% primal gap
stationarity_type_primal_gap = Data_primal_gap.Rec.x_Opt_type;
color_primal_gap = Data_primal_gap.Rec.x_Opt_color;
figure(2)
for i = 1 : size(x_Init, 1)
    for j = 1 : size(x_Init, 2)
        x_Init_ij = x_Init{i, j};
        plot(x_Init_ij(1), x_Init_ij(2), stationarity_type_primal_gap{i, j}, 'MarkerSize', 10, 'MarkerEdgeColor', color_primal_gap{i, j})
        hold on
    end    
end
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
axis([-1.5, 2.5, -1.5, 2.5])

%% D gap
stationarity_type_D_gap = Data_D_gap.Rec.x_Opt_type;
color_D_gap = Data_D_gap.Rec.x_Opt_color;
figure(3)
for i = 1 : size(x_Init, 1)
    for j = 1 : size(x_Init, 2)
        x_Init_ij = x_Init{i, j};
        plot(x_Init_ij(1), x_Init_ij(2), stationarity_type_D_gap{i, j}, 'MarkerSize', 10, 'MarkerEdgeColor', color_D_gap{i, j})
        hold on
    end    
end
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
axis([-1.5, 2.5, -1.5, 2.5])

%% Lin_Fukushima
stationarity_type_Lin_Fukushima = Data_Lin_Fukushima.Rec.x_Opt_type;
color_Lin_Fukushima = Data_Lin_Fukushima.Rec.x_Opt_color;
figure(4)
for i = 1 : size(x_Init, 1)
    for j = 1 : size(x_Init, 2)
        x_Init_ij = x_Init{i, j};
        plot(x_Init_ij(1), x_Init_ij(2), stationarity_type_Lin_Fukushima{i, j}, 'MarkerSize', 10, 'MarkerEdgeColor', color_Lin_Fukushima{i, j})
        hold on
    end    
end
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
axis([-1.5, 2.5, -1.5, 2.5])

%% Kadrani
stationarity_type_Kadrani = Data_Kadrani.Rec.x_Opt_type;
color_Kadrani = Data_Kadrani.Rec.x_Opt_color;

figure(5)
for i = 1 : size(x_Init, 1)
    for j = 1 : size(x_Init, 2)
        x_Init_ij = x_Init{i, j};
        plot(x_Init_ij(1), x_Init_ij(2), stationarity_type_Kadrani{i, j}, 'MarkerSize', 10, 'MarkerEdgeColor', color_Kadrani{i, j})
        hold on
    end    
end
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
axis([-1.5, 2.5, -1.5, 2.5])

