%% feasible set formed by relaxed generalized primal gap constraint for complementarity constraint
clear all
clc

c = 1;
s = 0.05;
plot_primal_gap_constraint_feasible_set(c, s)

%% feasible set formed by relaxed generalized D gap constraint for complementarity constraint
clear all
clc

a = 0.5;
b = 2; % b > a > 0
s = 0.05;
plot_D_gap_constraint_feasible_set(a, b, s)

%% sub function
function plot_primal_gap_constraint_feasible_set(c, s)
% parameter
stepsize = 0.0001;
% node
origin_x = 0;
origin_y = 0;
nodePoint1_x = sqrt(2*s/c);
nodePoint1_y = sqrt(2*c*s);
nodePoint2_x = 0;
nodePoint2_y = -sqrt(2*c*s);
ycxPoint1_x = 0.001;
ycxPoint1_y = c*ycxPoint1_x;
% figure limit
x_lb = origin_x - 0.2;
x_ub = nodePoint1_x + 1;
y_lb = -nodePoint1_y - 0.2;
y_ub = nodePoint1_y + 1;
% region
reg_curve_x = ycxPoint1_x : stepsize : nodePoint1_x;
reg_curve_y = zeros(1, length(reg_curve_x));
for i = 1 : length(reg_curve_x)
    reg_curve_y(i) = s/reg_curve_x(i) + c/2*reg_curve_x(i);
end
reg_x = [reg_curve_x, nodePoint1_x, x_ub,         x_ub,          nodePoint2_x];
reg_y = [reg_curve_y, nodePoint1_y, nodePoint1_y, -nodePoint1_y, nodePoint2_y];
% boundary
bound_x_geq_0_y = -nodePoint1_y : stepsize : reg_curve_y(1);
bound_x_geq_0_x = zeros(1, length(bound_x_geq_0_y));
% axis
bound_y_cx_x = x_lb : stepsize : x_ub;
bound_y_cx_y = c.* bound_y_cx_x;
% plot
figure(1)
% region
patch(reg_x, reg_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% node
plot(nodePoint1_x, nodePoint1_y, 'ko', 'LineWidth', 3)
hold on
plot(nodePoint2_x, nodePoint2_y, 'ko', 'LineWidth', 3)
hold on
% boundary
plot(bound_x_geq_0_x, bound_x_geq_0_y, 'k', 'LineWidth', 3)
hold on
plot(reg_curve_x, reg_curve_y, 'k', 'LineWidth', 3)
hold on
plot([nodePoint1_x, x_ub], [nodePoint1_y, nodePoint1_y], 'k', 'LineWidth', 3)
hold on
plot([origin_x, x_ub], [-nodePoint1_y, -nodePoint1_y], 'k', 'LineWidth', 3)
hold on
plot(bound_y_cx_x, bound_y_cx_y, 'k', 'LineStyle', '--', 'LineWidth', 3 )
hold on
text(0.8, 0.8, '$\eta = \lambda$', 'Interpreter','latex', 'FontSize', 20)
%text(0.5, 0.5, '$(\sqrt{2s}, \sqrt{2s})$', 'Interpreter','latex', 'FontSize', 20)
%text(0, -0.5, '$(0, -\sqrt{2s})$', 'Interpreter','latex', 'FontSize', 20)
axis([x_lb, x_lb + 1.5, y_lb, y_lb + 1.5])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
xticks([])
yticks([])
xlabel('$\lambda$', 'Interpreter','latex', 'FontSize', 20)
ylabel('$\eta$', 'Interpreter','latex', 'FontSize', 20)
end

function plot_D_gap_constraint_feasible_set(a, b, s)
% parameter
stepsize = 0.0001;
% node
nodePoint1_x = sqrt(2*s/(b-a));
nodePoint1_y = b*nodePoint1_x;
nodePoint2_x = sqrt(2*b*s/(a*(b-a)));
nodePoint2_y = sqrt(2*a*b*s/(b-a));
nodePoint3_x = -sqrt(2*s/(b-a));
nodePoint3_y = a*nodePoint3_x;
nodePoint4_x = -sqrt(2*a*s/(b*(b-a)));
nodePoint4_y = -sqrt(2*a*b*s/(b-a));
% figure limit
x_lb = nodePoint3_x - 1;
x_ub = nodePoint2_x + 1;
y_lb = nodePoint4_y - 1;
y_ub = nodePoint1_y + 1;
% region 
reg_2_curve_x = nodePoint2_x: -stepsize: nodePoint1_x;
reg_2_curve_y = zeros(1, length(reg_2_curve_x));
for i = 1 : length(reg_2_curve_x)
    reg_2_curve_y(i) = b*reg_2_curve_x(i) - b*sqrt((1-a/b)*(reg_2_curve_x(i))^2 - 2*s/b);
end
reg_4_curve_x = nodePoint3_x : stepsize : nodePoint4_x;
reg_4_curve_y = zeros(1, length(reg_4_curve_x));
for i = 1 : length(reg_4_curve_x)
    reg_4_curve_y(i) = a * reg_4_curve_x(i) - a * sqrt((1 - b/a)*(reg_4_curve_x(i))^2 + 2*s/a);
end
reg_x = [nodePoint3_x, nodePoint3_x, reg_4_curve_x, x_ub,...
    x_ub,         nodePoint2_x, reg_2_curve_x, nodePoint1_x];
reg_y = [y_ub,         nodePoint3_y, reg_4_curve_y, nodePoint4_y,...
    nodePoint2_y, nodePoint2_y, reg_2_curve_y, y_ub];
% axis
bound_y_ax_x = x_lb: stepsize : x_ub;
bound_y_ax_y = a.* bound_y_ax_x;
bound_y_bx_x = x_lb: stepsize : x_ub;
bound_y_bx_y = b.* bound_y_bx_x;
% plot
figure(2)
% region
patch(reg_x, reg_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% node
plot(nodePoint1_x, nodePoint1_y, 'ko', 'LineWidth', 3)
hold on
plot(nodePoint2_x, nodePoint2_y, 'ko', 'LineWidth', 3)
hold on
plot(nodePoint3_x, nodePoint3_y, 'ko', 'LineWidth', 3)
hold on
plot(nodePoint4_x, nodePoint4_y, 'ko', 'LineWidth', 3)
hold on
% boundary
plot([nodePoint3_x, nodePoint3_x], [y_ub, nodePoint3_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint1_x, nodePoint1_x], [y_ub, nodePoint1_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint2_x, x_ub], [nodePoint2_y, nodePoint2_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint4_x, x_ub], [nodePoint4_y, nodePoint4_y], 'k', 'LineWidth', 3)
hold on
plot(reg_2_curve_x, reg_2_curve_y, 'k', 'LineWidth', 3)% region 2 curve
hold on
plot(reg_4_curve_x, reg_4_curve_y, 'k', 'LineWidth', 3)% region 4 curve
hold on
% axis
plot(bound_y_ax_x, bound_y_ax_y, 'k', 'LineStyle', '--', 'LineWidth', 3 )
hold on
plot(bound_y_bx_x, bound_y_bx_y, 'k', 'LineStyle', '--', 'LineWidth', 3 )
hold on
text(0.6, 1.0, '$\eta = b\lambda$', 'Interpreter','latex', 'FontSize', 20)
text(1.0, 0.6, '$\eta = a\lambda$', 'Interpreter','latex', 'FontSize', 20)
%text(-0.4, 0.8, '$(\sqrt{\frac{2s}{b-a}}, \sqrt{\frac{2b^2s}{b-a}})$', 'Interpreter','latex', 'FontSize', 20)
%text(0.8, -0.4, '$(\sqrt{\frac{2bs}{a(b-a)}}, \sqrt{\frac{2abs}{b-a}})$', 'Interpreter','latex', 'FontSize', 20)
%text(-0.4, -0.4, '$(-\sqrt{\frac{2s}{b-a}}, \sqrt{\frac{2a^2s}{b-a}})$', 'Interpreter','latex', 'FontSize', 20)
%text(0.8, -0.8, '$(-\sqrt{\frac{2as}{b(b-a)}}, -\sqrt{\frac{2abs}{b-a}})$', 'Interpreter','latex', 'FontSize', 20)
axis([nodePoint3_x-0.2, nodePoint2_x+1.0, nodePoint4_y-0.2, nodePoint1_y+1.0])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
xticks([])
yticks([])
xlabel('$\lambda$', 'Interpreter','latex', 'FontSize', 20)
ylabel('$\eta$', 'Interpreter','latex', 'FontSize', 20)
end