%% feasible set formed by relaxed generalized primal gap constraint for BVI
clear all
clc

c = 1;
s = 0.1;
bl = -1;
bu = 1;
plot_primal_gap_constraint_feasible_set(c, s, bl, bu)

%% feasible set formed by relaxed generalized D gap constraint for BVI
clear all
clc

a = 0.5;
b = 2;
s = 0.1;
bl = -1;
bu = 1;
plot_D_gap_constraint_feasible_set(a, b, s, bl, bu)

%% sub function
function plot_primal_gap_constraint_feasible_set(c, s, bl, bu)
% parameter
stepsize = 0.0001;
% node
nodePoint1_x = bl;
nodePoint1_y = -sqrt(2*c*s);
nodePoint2_x = bu - sqrt(2*s/c);
nodePoint2_y = -sqrt(2*c*s);
nodePoint3_x = bu;
nodePoint3_y = sqrt(2*c*s);
nodePoint4_x = sqrt(2*s/c) + bl;
nodePoint4_y = sqrt(2*c*s);
% figure limit
x_lb = bl - 0.3;
x_ub = bu + 0.3;
y_lb = nodePoint1_y - 1;
y_ub = nodePoint3_y + 1;
% region
reg_1_curve_x = nodePoint1_x + stepsize : stepsize : nodePoint4_x;
reg_1_curve_y = zeros(1, length(reg_1_curve_x));
for i = 1 : length(reg_1_curve_x)
     reg_1_curve_y(i) = s/(reg_1_curve_x(i) - bl) + c/2*(reg_1_curve_x(i) - bl);
end
reg_3_curve_x = nodePoint3_x - stepsize : -stepsize: nodePoint2_x;
reg_3_curve_y = zeros(1, length(reg_3_curve_x));
for i = 1 : length(reg_3_curve_x)
    reg_3_curve_y(i) = s/(reg_3_curve_x(i) - bu) + c/2*(reg_3_curve_x(i) - bu);
end
reg_x = [reg_1_curve_x, nodePoint3_x, reg_3_curve_x, nodePoint1_x];
reg_y = [reg_1_curve_y, nodePoint3_y, reg_3_curve_y, nodePoint1_y];
% boundary
boundary_bl_y =  -sqrt(2*c*s): stepsize : reg_1_curve_y(1);
boundary_bl_x = bl * ones(1, length(boundary_bl_y));
boundary_bu_y = reg_3_curve_y(1): stepsize :  sqrt(2*c*s);
boundary_bu_x = bu * ones(1, length(boundary_bu_y));
% axis
boundary_c_bl_x = bl - 1 : stepsize : bu + 1;
boundary_c_bl_y = c*(boundary_c_bl_x - bl);
boundary_c_bu_x = bl - 1 : stepsize : bu + 1;
boundary_c_bu_y = c*(boundary_c_bu_x - bu);
% plot
figure(1)
% region
patch(reg_x, reg_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% node
plot(nodePoint4_x, nodePoint4_y, 'ko', 'LineWidth', 3)
hold on
plot(nodePoint2_x, nodePoint2_y, 'ko', 'LineWidth', 3)
hold on
% boundary
plot(boundary_bl_x, boundary_bl_y, 'k', 'LineWidth', 3)
hold on
plot(boundary_bu_x, boundary_bu_y, 'k', 'LineWidth', 3)
hold on
plot(reg_1_curve_x, reg_1_curve_y, 'k', 'LineWidth', 3)
hold on
plot(reg_3_curve_x, reg_3_curve_y, 'k', 'LineWidth', 3)
hold on
plot([nodePoint4_x, nodePoint3_x], [nodePoint4_y, nodePoint3_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint1_x, nodePoint2_x], [nodePoint1_y, nodePoint2_y], 'k', 'LineWidth', 3)
hold on
% axes
plot(boundary_c_bl_x, boundary_c_bl_y, 'k', 'LineStyle', '--', 'LineWidth', 3 )
hold on
plot(boundary_c_bu_x, boundary_c_bu_y, 'k', 'LineStyle', '--', 'LineWidth', 3 )
text(-1, 1, '$\eta = \lambda - b_l$', 'Interpreter','latex', 'FontSize', 20)
text(1, -1, '$\eta = \lambda - b_u$', 'Interpreter','latex', 'FontSize', 20)
text(-1, -1, '$\lambda = b_l$', 'Interpreter','latex', 'FontSize', 20)
text(1, 1, '$\lambda = b_u$', 'Interpreter','latex', 'FontSize', 20)
axis([x_lb, x_ub, y_lb, y_ub])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
xticks([])
yticks([])
xlabel('$\lambda$', 'Interpreter','latex', 'FontSize', 20)
ylabel('$\eta$', 'Interpreter','latex', 'FontSize', 20)

end

function plot_D_gap_constraint_feasible_set(a, b, s, bl, bu)
% parameter
stepsize = 0.0001;
% nodes
nodePoint1_x = bl + sqrt(2*s/(b-a));
nodePoint1_y = b * (nodePoint1_x - bl);
nodePoint2_y = sqrt(2*a*b*s/(b-a));
nodePoint2_x = 1/a * nodePoint2_y + bl;
nodePoint3_x = bl - sqrt(2*s/(b-a));
nodePoint3_y = a * (nodePoint3_x - bl);
nodePoint4_y = - sqrt(2*a*b*s/(b-a));
nodePoint4_x = 1/b * nodePoint4_y + bl;
nodePoint5_y = sqrt(2*a*b*s/(b-a));
nodePoint5_x = 1/b * nodePoint5_y + bu;
nodePoint6_x = bu + sqrt(2*s/(b-a));
nodePoint6_y = a * (nodePoint6_x - bu);
nodePoint7_y = -sqrt(2*a*b*s/(b-a));
nodePoint7_x = 1/a * nodePoint7_y + bu;
nodePoint8_x = bu - sqrt(2*s/(b-a));
nodePoint8_y = b * (nodePoint8_x - bu);
% figure limit
x_lb = bl - 0.6;
x_ub = bu + 0.6;
y_lb = nodePoint8_y - 1.2;
y_ub = nodePoint1_y + 1.2;
% region
reg_2_curve_x = nodePoint1_x: stepsize: nodePoint2_x;
reg_2_curve_y = zeros(1, length(reg_2_curve_x));
for i = 1 : length(reg_2_curve_x)
    reg_2_curve_y(i) = b * (reg_2_curve_x(i) - bl) - b * sqrt((1 - a/b)*(reg_2_curve_x(i) - bl)^2 - 2*s/b);
end
reg_3_curve_x = nodePoint4_x : -stepsize : nodePoint3_x;
reg_3_curve_y = zeros(1, length(reg_3_curve_x));
for i = 1 : length(reg_3_curve_x)
    reg_3_curve_y(i) = a * (reg_3_curve_x(i) - bl) - a * sqrt((1 - b/a)*(reg_3_curve_x(i) - bl)^2 + 2*s/a);
end
reg_6_curve_x = nodePoint5_x : stepsize : nodePoint6_x;
reg_6_curve_y = zeros(1, length(reg_6_curve_x));
for i = 1 : length(reg_6_curve_x)
    reg_6_curve_y(i) = a * (reg_6_curve_x(i) - bu) + a * sqrt((1 - b/a)*(reg_6_curve_x(i) - bu)^2 + 2*s/a);
end
reg_5_curve_x = nodePoint8_x : -stepsize : nodePoint7_x;
reg_5_curve_y = zeros(1, length(reg_5_curve_x));
for i = 1 : length(reg_5_curve_x)
    reg_5_curve_y(i) = b * (reg_5_curve_x(i) - bu) + b * sqrt((1 - a/b)*(reg_5_curve_x(i) - bu)^2 - 2*s/b);
end

reg_x = [nodePoint1_x, reg_2_curve_x, reg_6_curve_x, nodePoint6_x,...
    nodePoint8_x, reg_5_curve_x, reg_3_curve_x, nodePoint3_x];
reg_y = [y_ub,         reg_2_curve_y, reg_6_curve_y, y_lb,...
    y_lb,         reg_5_curve_y, reg_3_curve_y, y_ub];

% axis
axis_y_bx_bl_x = x_lb: stepsize : (x_lb + x_ub)/2;
axis_y_bx_bl_y = b.* (axis_y_bx_bl_x - bl);
axis_y_ax_bl_x = x_lb: stepsize : (x_lb + x_ub)/2;
axis_y_ax_bl_y = a.* (axis_y_ax_bl_x - bl);

axis_y_bx_bu_x = (x_lb + x_ub)/2: stepsize : x_ub;
axis_y_bx_bu_y = b.* (axis_y_bx_bu_x - bu);
axis_y_ax_bu_x = (x_lb + x_ub)/2: stepsize : x_ub;
axis_y_ax_bu_y = a.* (axis_y_ax_bu_x - bu);

axis_bl_y = 0 : stepsize : y_ub;
axis_bl_x = bl * ones(1, length(axis_bl_y));

axis_bu_y = y_lb : stepsize : 0;
axis_bu_x = bu * ones(1, length(axis_bu_y));

% plot
figure(4)
% region
patch(reg_x, reg_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% node
plot([nodePoint1_x, nodePoint2_x, nodePoint3_x, nodePoint4_x, nodePoint5_x, nodePoint6_x, nodePoint7_x, nodePoint8_x],...
    [nodePoint1_y, nodePoint2_y, nodePoint3_y, nodePoint4_y, nodePoint5_y, nodePoint6_y, nodePoint7_y, nodePoint8_y],...
    'ko', 'LineWidth', 3, 'LineStyle', 'none')
hold on
% boundary
plot([nodePoint3_x, nodePoint3_x], [y_ub, nodePoint3_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint1_x, nodePoint1_x], [y_ub, nodePoint1_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint2_x, nodePoint5_x], [nodePoint2_y, nodePoint5_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint4_x, nodePoint7_x], [nodePoint4_y, nodePoint7_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint8_x, nodePoint8_x], [nodePoint8_y, y_lb], 'k', 'LineWidth', 3)
hold on
plot([nodePoint6_x, nodePoint6_x], [nodePoint6_y, y_lb], 'k', 'LineWidth', 3)
hold on
plot(reg_2_curve_x, reg_2_curve_y, 'k', 'LineWidth', 3)
hold on
plot(reg_3_curve_x, reg_3_curve_y, 'k', 'LineWidth', 3)
hold on
plot(reg_5_curve_x, reg_5_curve_y, 'k', 'LineWidth', 3)
hold on
plot(reg_6_curve_x, reg_6_curve_y, 'k', 'LineWidth', 3)
hold on
% axis
plot(axis_y_bx_bl_x, axis_y_bx_bl_y, 'k', 'LineStyle', '--', 'LineWidth', 3 )
hold on
plot(axis_y_ax_bl_x, axis_y_ax_bl_y, 'k', 'LineStyle', '--', 'LineWidth', 3 )
hold on
plot(axis_y_bx_bu_x, axis_y_bx_bu_y, 'k', 'LineStyle', '--', 'LineWidth', 3 )
hold on
plot(axis_y_ax_bu_x, axis_y_ax_bu_y, 'k', 'LineStyle', '--', 'LineWidth', 3 )
hold on
plot(axis_bl_x, axis_bl_y, 'k', 'LineStyle', '--', 'LineWidth', 1.5)
hold on
plot(axis_bu_x, axis_bu_y, 'k', 'LineStyle', '--', 'LineWidth', 1.5)
text(0, 1.7, '$\eta = b(\lambda - b_l)$', 'Interpreter','latex', 'FontSize', 20)
text(0, 0.7, '$\eta = a(\lambda - b_l)$', 'Interpreter','latex', 'FontSize', 20)
text(-1.35, 1.2, '$\lambda = b_l$', 'Interpreter','latex', 'FontSize', 20)
text(-1.2, -1.7, '$\eta = b(\lambda - b_u)$', 'Interpreter','latex', 'FontSize', 20)
text(-1.2, -0.7, '$\eta = a(\lambda - b_u)$', 'Interpreter','latex', 'FontSize', 20)
text(0.65, -1.2, '$\lambda = b_u$', 'Interpreter','latex', 'FontSize', 20)
axis([x_lb, x_ub, y_lb, y_ub])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
xticks([])
yticks([])
xlabel('$\lambda$', 'Interpreter','latex', 'FontSize', 20)
ylabel('$\eta$', 'Interpreter','latex', 'FontSize', 20)
end