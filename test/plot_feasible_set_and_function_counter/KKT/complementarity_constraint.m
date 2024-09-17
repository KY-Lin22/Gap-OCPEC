%% feasible set formed by Scholtes relaxation stragegy
clear all
clc
s = 0.1;
plot_Scholtes_feasible_set(s)

%% feasible set formed by Lin-Fukushima relaxation stragegy
clear all
clc
s = 0.2;
plot_Lin_Fukushima_feasible_set(s)

%% feasible set formed by Kadrani relaxation stragegy
clear all
clc
s = 0.1;
plot_Kadrani_feasible_set(s)

%% feasible set formed by Steffensen-Ulbrich relaxation stragegy
clear all
clc
s = 0.4;
plot_Steffensen_Ulbrich_feasible_set(s)

%% feasible set formed by Kanzow-Schwartz relaxation stragegy
clear all
clc
s = 0.1;
plot_Kanzow_Schwartz_feasible_set(s)

%% sub function
function plot_Scholtes_feasible_set(s)
% parameter
stepsize = 0.01;
% node
origin_x = 0;
origin_y = 0;
% figure limit
x_lb = origin_x - 0.2;
x_ub = origin_x + 1;
y_lb = origin_y - 0.2;
y_ub = origin_y + 1;
% curve
curve_x = origin_x + stepsize : stepsize : x_ub;
curve_y = s ./ curve_x;
% region
reg_x = [origin_x, origin_x,   curve_x, curve_x(end)];
reg_y = [origin_y, curve_y(1), curve_y, origin_y];
% boundary
bound_xAxis_x = origin_x : stepsize : x_ub;
bound_xAxis_y = zeros(1, length(bound_xAxis_x));
bound_yAxis_y = origin_y : stepsize : y_ub;
bound_yAxis_x = zeros(1, length(bound_yAxis_y));
% plot
figure(1)
% region
patch(reg_x, reg_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% boundary
plot(bound_xAxis_x, bound_xAxis_y, 'k', 'LineWidth', 3)
hold on
plot(bound_yAxis_x, bound_yAxis_y, 'k', 'LineWidth', 3)
hold on
plot(curve_x, curve_y, 'k', 'LineWidth', 3)
axis([x_lb, x_ub, y_lb, y_ub])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
xticks([])
yticks([])
xlabel('$\zeta_i$', 'Interpreter','latex', 'FontSize', 30)
ylabel('$g_i$', 'Interpreter','latex', 'FontSize', 30)
end

function plot_Lin_Fukushima_feasible_set(s)
% parameter
stepsize = 0.01;
% node
origin_x = 0;
origin_y = 0;
% figure limit
x_lb = origin_x - 0.2;
x_ub = origin_x + 1;
y_lb = origin_y - 0.2;
y_ub = origin_y + 1;
% curve_1
curve_1_x = origin_x + stepsize : stepsize : x_ub;
curve_1_y = s^2 ./ curve_1_x;
% curve 2
curve_2_x = -s + stepsize: stepsize : x_ub;
curve_2_y = s^2 ./ (curve_2_x + s) - s;
% region
reg_x = [curve_1_x, curve_1_x(end), flip(curve_2_x), curve_2_x(1)];
reg_y = [curve_1_y, curve_2_y(end), flip(curve_2_y), curve_1_y(1)];
% plot
figure(1)
% region
patch(reg_x, reg_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% boundary
hold on
plot(curve_1_x, curve_1_y, 'k', 'LineWidth', 3)
hold on
plot(curve_2_x, curve_2_y, 'k', 'LineWidth', 3)
axis([x_lb, x_ub, y_lb, y_ub])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
xticks([])
yticks([])
xlabel('$\zeta_i$', 'Interpreter','latex', 'FontSize', 30)
ylabel('$g_i$', 'Interpreter','latex', 'FontSize', 30)
end

function plot_Kadrani_feasible_set(s)
% node
origin_x = 0;
origin_y = 0;
% figure limit
x_lb = origin_x - 0.2;
x_ub = origin_x + 1;
y_lb = origin_y - 0.2;
y_ub = origin_y + 1;
% region
reg_1_x = [-s,   -s, s, s];
reg_1_y = [y_ub, s,  s, y_ub];
reg_2_x = [x_ub, s, s,  x_ub];
reg_2_y = [s,    s, -s, -s];
% plot
figure(1)
% region
patch(reg_1_x, reg_1_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_2_x, reg_2_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% boundary
plot(reg_1_x, reg_1_y, 'k', 'LineWidth', 3)
hold on
plot(reg_2_x, reg_2_y, 'k', 'LineWidth', 3)
axis([x_lb, x_ub, y_lb, y_ub])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
xticks([])
yticks([])
xlabel('$\zeta_i$', 'Interpreter','latex', 'FontSize', 30)
ylabel('$g_i$', 'Interpreter','latex', 'FontSize', 30)
end

function plot_Steffensen_Ulbrich_feasible_set(s)
% parameter
stepsize = 0.01;
% node
origin_x = 0;
origin_y = 0;
% figure limit
x_lb = origin_x - 0.2;
x_ub = origin_x + 1;
y_lb = origin_y - 0.2;
y_ub = origin_y + 1;
% curve
curve_x = origin_x: stepsize : s;
curve_y = zeros(1, length(curve_x));
for i = 1 : length(curve_x)
    fun = @(y) (curve_x(i) + y - s*(1/8)*(-((curve_x(i) - y)/s)^4 + 6*((curve_x(i) - y)/s)^2 + 3));
    y0 = 1;
    [y_i,~] = fsolve(fun, y0);
    curve_y(1, i) = y_i;
end
% region
reg_x = [origin_x, origin_x,   curve_x, curve_x(end)];
reg_y = [origin_y, curve_y(1), curve_y, origin_y];
% boundary
bound_xAxis_x = origin_x : stepsize : x_ub;
bound_xAxis_y = zeros(1, length(bound_xAxis_x));
bound_yAxis_y = origin_y : stepsize : y_ub;
bound_yAxis_x = zeros(1, length(bound_yAxis_y));
% plot
figure(1)
% region
patch(reg_x, reg_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% boundary
plot(bound_xAxis_x, bound_xAxis_y, 'k', 'LineWidth', 3)
hold on
plot(bound_yAxis_x, bound_yAxis_y, 'k', 'LineWidth', 3)
hold on
plot(curve_x, curve_y, 'k', 'LineWidth', 3)
axis([x_lb, x_ub, y_lb, y_ub])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
xticks([])
yticks([])
xlabel('$\zeta_i$', 'Interpreter','latex', 'FontSize', 30)
ylabel('$g_i$', 'Interpreter','latex', 'FontSize', 30)
end

function plot_Kanzow_Schwartz_feasible_set(s)
% node
origin_x = 0;
origin_y = 0;
% figure limit
x_lb = origin_x - 0.2;
x_ub = origin_x + 1;
y_lb = origin_y - 0.2;
y_ub = origin_y + 1;
% region
reg_x = [origin_x, origin_x, x_ub,     x_ub, s, s];
reg_y = [y_ub,     origin_y, origin_y, s,    s  y_ub];
% boundary
bound_1_x = [origin_x, origin_x, x_ub];
bound_1_y = [y_ub,     origin_y, origin_y];
bound_2_x = [s,    s, x_ub];
bound_2_y = [y_ub, s, s];
% plot
figure(1)
% region
patch(reg_x, reg_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% boundary
plot(bound_1_x, bound_1_y, 'k', 'LineWidth', 3)
hold on
plot(bound_2_x, bound_2_y, 'k', 'LineWidth', 3)
axis([x_lb, x_ub, y_lb, y_ub])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
xticks([])
yticks([])
xlabel('$\zeta_i$', 'Interpreter','latex', 'FontSize', 30)
ylabel('$g_i$', 'Interpreter','latex', 'FontSize', 30)
end
