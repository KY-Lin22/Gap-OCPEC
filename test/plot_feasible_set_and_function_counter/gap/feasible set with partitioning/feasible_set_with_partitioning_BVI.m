
%% feasible set formed by relaxed generalized primal gap constraint 
clear all
clc

c = 1;
s = 0.1;
bl = -1;
bu = 1;
plotPrimalGapConstraintSet(c, s, bl, bu)

%% feasible set formed by relaxed generalized D gap constraint 
clear all
clc

a = 0.3;
b = 2;
s = 0.1;
bl = -1;
bu = 1;
plotDGapConstraintSet(a, b, s, bl, bu)



%% sub function
function plotPrimalGapConstraintSet(c, s, bl, bu)
% parameter
stepsize = 0.01;
% node
bl_Point_x = bl + 0.001;
bl_Point_y = c*(bl_Point_x - bl);
bu_Point_x = bu - 0.001;
bu_Point_y = c*(bu_Point_x - bu);
nodePoint1_x = bl;
nodePoint1_y = 0;
nodePoint2_x = sqrt(2*s/c) + bl;
nodePoint2_y = sqrt(2*c*s);
nodePoint3_x = bu;
nodePoint3_y = sqrt(2*c*s);
nodePoint4_x = bu;
nodePoint4_y = 0;
nodePoint5_x = bu - sqrt(2*s/c);
nodePoint5_y = -sqrt(2*c*s);
nodePoint6_x = bl;
nodePoint6_y = -sqrt(2*c*s);
% region 1
reg_1_curve_x = bl_Point_x : stepsize : nodePoint2_x;
reg_1_curve_y = zeros(1, length(reg_1_curve_x));
for i = 1 : length(reg_1_curve_x)
     reg_1_curve_y(i) = s/(reg_1_curve_x(i) - bl) + c/2*(reg_1_curve_x(i) - bl);
end
reg_1_x = [reg_1_curve_x, bl_Point_x];
reg_1_y = [reg_1_curve_y, bl_Point_y];
% region 2
reg_2_x = [nodePoint1_x, nodePoint2_x, nodePoint3_x, nodePoint4_x, nodePoint5_x, nodePoint6_x];
reg_2_y = [nodePoint1_y, nodePoint2_y, nodePoint3_y, nodePoint4_y, nodePoint5_y, nodePoint6_y];
% region 3
reg_3_curve_x = nodePoint5_x : stepsize: bu_Point_x;
reg_3_curve_y = zeros(1, length(reg_3_curve_x));
for i = 1 : length(reg_3_curve_x)
    reg_3_curve_y(i) = s/(reg_3_curve_x(i) - bu) + c/2*(reg_3_curve_x(i) - bu);
end
reg_3_x = [bu_Point_x, reg_3_curve_x];
reg_3_y = [bu_Point_y, reg_3_curve_y];
% boundary
boundary_bl_y =  -sqrt(2*c*s): stepsize : reg_1_curve_y(1);
boundary_bl_x = bl * ones(1, length(boundary_bl_y));
boundary_bu_y = reg_3_curve_y(end): stepsize :  sqrt(2*c*s);
boundary_bu_x = bu * ones(1, length(boundary_bu_y));
boundary_c_bl_x = bl - 1 : stepsize : bu + 1;
boundary_c_bl_y = c*(boundary_c_bl_x - bl);
boundary_c_bu_x = bl - 1 : stepsize : bu + 1;
boundary_c_bu_y = c*(boundary_c_bu_x - bu);
% plot
figure(1)
% original region

% relaxed region
patch(reg_1_x, reg_1_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_2_x, reg_2_y, 'red', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_3_x, reg_3_y, 'green', 'FaceAlpha', 0.1, 'LineStyle', 'none')
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
plot([nodePoint2_x, nodePoint3_x], [nodePoint2_y, nodePoint3_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint6_x, nodePoint5_x], [nodePoint6_y, nodePoint5_y], 'k', 'LineWidth', 3)
hold on
% axes
plot(boundary_c_bl_x, boundary_c_bl_y, 'k', 'LineStyle', '--', 'LineWidth', 1 )
hold on
plot(boundary_c_bu_x, boundary_c_bu_y, 'k', 'LineStyle', '--', 'LineWidth', 1 )
% node
% plot([nodePoint1_x, nodePoint2_x, nodePoint3_x, nodePoint4_x, nodePoint5_x, nodePoint6_x],...
%     [nodePoint1_y, nodePoint2_y, nodePoint3_y, nodePoint4_y, nodePoint5_y, nodePoint6_y],  'r*', 'LineWidth', 5)
axis([bl - 0.3, bu + 0.3, nodePoint5_y - 0.8, nodePoint2_y + 0.8])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
end



function plotDGapConstraintSet(a, b, s, bl, bu)
% parameter
stepsize = 0.01;

% node
nodePoint1_x = bl + sqrt(2*s/(b-a));
nodePoint1_y = b * (nodePoint1_x - bl);
nodePoint2_y = sqrt(2*a*b*s/(b-a));
nodePoint2_x = 1/a * nodePoint2_y + bl; 
nodePoint3_x = bl;
nodePoint3_y = 0;
nodePoint4_x = bl - sqrt(2*s/(b-a));
nodePoint4_y = a * (nodePoint4_x - bl);
nodePoint5_y = - sqrt(2*a*b*s/(b-a));
nodePoint5_x = 1/b * nodePoint5_y + bl;
nodePoint6_y = sqrt(2*a*b*s/(b-a));
nodePoint6_x = 1/b * nodePoint6_y + bu;
nodePoint7_x = bu + sqrt(2*s/(b-a));
nodePoint7_y = a * (nodePoint7_x - bu);
nodePoint8_x = bu;
nodePoint8_y = 0;
nodePoint9_y = -sqrt(2*a*b*s/(b-a));
nodePoint9_x = 1/a * nodePoint9_y + bu;
nodePoint10_x = bu - sqrt(2*s/(b-a));
nodePoint10_y = b * (nodePoint10_x - bu);

% figure boundary
x_lb = nodePoint4_x - 1;
x_ub = nodePoint7_x + 1;
y_lb = nodePoint10_y - 1;
y_ub = nodePoint1_y + 1;

% region 1
reg_1_x = [nodePoint4_x, nodePoint4_x, nodePoint3_x, nodePoint1_x, nodePoint1_x];
reg_1_y = [y_ub, nodePoint4_y, nodePoint3_y, nodePoint1_y, y_ub];

% region 2
reg_2_curve_x = nodePoint1_x: stepsize: nodePoint2_x;
reg_2_curve_y = zeros(1, length(reg_2_curve_x));
for i = 1 : length(reg_2_curve_x)
    reg_2_curve_y(i) = b * (reg_2_curve_x(i) - bl) - b * sqrt((1 - a/b)*(reg_2_curve_x(i) - bl)^2 - 2*s/b);
end
reg_2_x = [reg_2_curve_x, nodePoint3_x];
reg_2_y = [reg_2_curve_y, nodePoint3_y];

% region 3
reg_3_curve_x = nodePoint4_x : stepsize : nodePoint5_x;
reg_3_curve_y = zeros(1, length(reg_3_curve_x));
for i = 1 : length(reg_3_curve_x)
    reg_3_curve_y(i) = a * (reg_3_curve_x(i) - bl) - a * sqrt((1 - b/a)*(reg_3_curve_x(i) - bl)^2 + 2*s/a);
end
reg_3_x = [reg_3_curve_x, nodePoint3_x];
reg_3_y = [reg_3_curve_y, nodePoint3_y];

% region 4
reg_4_x = [nodePoint3_x, nodePoint2_x, nodePoint6_x, nodePoint8_x, nodePoint9_x, nodePoint5_x];
reg_4_y = [nodePoint3_y, nodePoint2_y, nodePoint6_y, nodePoint8_y, nodePoint9_y, nodePoint5_y];

% region 5
reg_5_curve_x = nodePoint10_x : -stepsize : nodePoint9_x;
reg_5_curve_y = zeros(1, length(reg_5_curve_x));
for i = 1 : length(reg_5_curve_x)
    reg_5_curve_y(i) = b * (reg_5_curve_x(i) - bu) + b * sqrt((1 - a/b)*(reg_5_curve_x(i) - bu)^2 - 2*s/b);
end
reg_5_x = [reg_5_curve_x, nodePoint8_x];
reg_5_y = [reg_5_curve_y, nodePoint8_y];

% region 6
reg_6_curve_x = nodePoint6_x : stepsize : nodePoint7_x;
reg_6_curve_y = zeros(1, length(reg_6_curve_x));
for i = 1 : length(reg_6_curve_x)
    reg_6_curve_y(i) = a * (reg_6_curve_x(i) - bu) + a * sqrt((1 - b/a)*(reg_6_curve_x(i) - bu)^2 + 2*s/a);
end
reg_6_x = [reg_6_curve_x, nodePoint8_x];
reg_6_y = [reg_6_curve_y, nodePoint8_y];

% region 7
reg_7_x = [nodePoint7_x, nodePoint7_x, nodePoint8_x, nodePoint10_x, nodePoint10_x];
reg_7_y = [y_lb, nodePoint7_y, nodePoint8_y, nodePoint10_y, y_lb];

% axis
axis_y_bx_bl_x = nodePoint4_x: stepsize : nodePoint2_x;
axis_y_bx_bl_y = b.* (axis_y_bx_bl_x - bl);
axis_y_ax_bl_x = x_lb: stepsize : (x_lb + x_ub)/2;
axis_y_ax_bl_y = a.* (axis_y_ax_bl_x - bl);

axis_y_bx_bu_x = nodePoint9_x: stepsize : nodePoint7_x;
axis_y_bx_bu_y = b.* (axis_y_bx_bu_x - bu);
axis_y_ax_bu_x = (x_lb + x_ub)/2: stepsize : x_ub;
axis_y_ax_bu_y = a.* (axis_y_ax_bu_x - bu);


axis_bl_y = y_lb : stepsize : y_ub;
axis_bl_x = bl * ones(1, length(axis_bl_y));

axis_bu_y = y_lb : stepsize : y_ub;
axis_bu_x = bu * ones(1, length(axis_bu_y));

% plot
figure(4)
% region
patch(reg_1_x, reg_1_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_2_x, reg_2_y, 'red', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_3_x, reg_3_y, 'yellow', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_4_x, reg_4_y, 'green', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_5_x, reg_5_y, 'red', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_6_x, reg_6_y, 'yellow', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_7_x, reg_7_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% boundary
plot([nodePoint4_x, nodePoint4_x], [y_ub, nodePoint4_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint1_x, nodePoint1_x], [y_ub, nodePoint1_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint2_x, nodePoint6_x], [nodePoint2_y, nodePoint6_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint5_x, nodePoint9_x], [nodePoint5_y, nodePoint9_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint10_x, nodePoint10_x], [nodePoint10_y, y_lb], 'k', 'LineWidth', 3)
hold on
plot([nodePoint7_x, nodePoint7_x], [nodePoint7_y, y_lb], 'k', 'LineWidth', 3)
hold on
f2 = @(x,y) -a./2 .* (x - bl).^2 + (x - bl).* y - 1./(2.* b) .* y.^2 - s;
fimplicit(f2, [nodePoint1_x, nodePoint2_x, nodePoint2_y, nodePoint1_y], 'k', 'LineWidth', 3) % region 2 curve
hold on
f3 = @(x,y) b./2 .* (x - bl).^2 - (x - bl).* y + 1./(2.* a) .* y.^2 - s;
fimplicit(f3, [nodePoint4_x, nodePoint5_x, nodePoint5_y, nodePoint4_y], 'k', 'LineWidth', 3) % region 3 curve
hold on
f5 = @(x,y) -a./2 .* (x - bu).^2 + (x - bu).* y - 1./(2.* b) .* y.^2 - s;
fimplicit(f5, [nodePoint9_x, nodePoint10_x, nodePoint10_y, nodePoint9_y], 'k', 'LineWidth', 3) % region 5 curve
hold on
f6 = @(x,y) b./2 .* (x - bu).^2 - (x - bu).* y + 1./(2.* a) .* y.^2 - s;
fimplicit(f6, [nodePoint6_x, nodePoint7_x, nodePoint7_y, nodePoint6_y], 'k', 'LineWidth', 3) % region 6 curve
hold on
% axis
plot(axis_y_bx_bl_x, axis_y_bx_bl_y, 'k', 'LineStyle', '--', 'LineWidth', 1 )
hold on
plot(axis_y_ax_bl_x, axis_y_ax_bl_y, 'k', 'LineStyle', '--', 'LineWidth', 1 )
hold on
plot(axis_y_bx_bu_x, axis_y_bx_bu_y, 'k', 'LineStyle', '--', 'LineWidth', 1 )
hold on
plot(axis_y_ax_bu_x, axis_y_ax_bu_y, 'k', 'LineStyle', '--', 'LineWidth', 1 )
hold on
plot(axis_bl_x, axis_bl_y, 'k', 'LineStyle', '--', 'LineWidth', 1)
hold on
plot(axis_bu_x, axis_bu_y, 'k', 'LineStyle', '--', 'LineWidth', 1)
axis([nodePoint4_x-0.4, nodePoint7_x+0.4, nodePoint10_y-0.4, nodePoint1_y+0.4])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
end



