
%% feasible set formed by Scholtes relaxation strategy


%% feasible set formed by relaxed generalized primal gap constraint 
clear all
clc

c = 1;
s = 0.1;
plotPrimalGapConstraintSet(c, s)

%% counter of generalized primal gap function
clear all
clc

c = 1;
plotPrimalGapFunContour(c)

%% feasible set formed by relaxed generalized D gap constraint 
clear all
clc

a = 0.3;
b = 2; % b > a > 0
s = 0.1;
plotDGapConstraintSet(a, b, s)

%% counter of generalized D gap function
clear all
clc

a = 0.3;
b = 2;
plotDGapFunContour(a, b)

%% sub function
function plotPrimalGapConstraintSet(c, s)
% parameter
stepsize = 0.01;

% node
origin_x = 0;
origin_y = 0;
nodePoint1_x = sqrt(2*s/c);
nodePoint1_y = sqrt(2*c*s);
ycxPoint1_x = 0.001;
ycxPoint1_y = c*ycxPoint1_x;
% figure boundary
x_lb = origin_x - 0.2;
x_ub = nodePoint1_x + 1;
y_lb = -nodePoint1_y - 0.2;
y_ub = nodePoint1_y + 1;

% region 1
reg_1_x = [origin_x, nodePoint1_x, x_ub, x_ub, origin_x];
reg_1_y = [origin_y, nodePoint1_y, nodePoint1_y, -nodePoint1_y, -nodePoint1_y];

% region 2
reg_2_curve_x = ycxPoint1_x : stepsize : nodePoint1_x;
reg_2_curve_y = zeros(1, length(reg_2_curve_x));
for i = 1 : length(reg_2_curve_x)
    reg_2_curve_y(i) = s/reg_2_curve_x(i) + c/2*reg_2_curve_x(i);
end
reg_2_x = [reg_2_curve_x, ycxPoint1_x];
reg_2_y = [reg_2_curve_y, ycxPoint1_y];

% boundary
bound_x_geq_0_y = -nodePoint1_y : stepsize : reg_2_curve_y(1);
bound_x_geq_0_x = zeros(1, length(bound_x_geq_0_y));
bound_y_cx_x = x_lb : stepsize : x_ub;
bound_y_cx_y = c.* bound_y_cx_x;

% plot
figure(1)
% region
patch(reg_1_x, reg_1_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_2_x, reg_2_y, 'red', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
% boundary
plot(bound_x_geq_0_x, bound_x_geq_0_y, 'k', 'LineWidth', 3)
hold on
plot(reg_2_curve_x, reg_2_curve_y, 'k', 'LineWidth', 3)
hold on
plot([nodePoint1_x, x_ub], [nodePoint1_y, nodePoint1_y], 'k', 'LineWidth', 3)
hold on
plot([origin_x, x_ub], [-nodePoint1_y, -nodePoint1_y], 'k', 'LineWidth', 3)
hold on
plot(bound_y_cx_x, bound_y_cx_y, 'k', 'LineStyle', '--', 'LineWidth', 1 )
hold on
axis([x_lb, x_ub, y_lb, y_ub])
xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
end

function plotPrimalGapFunContour(c)
% parameter
stepsize = 0.01;
% data
x = -1 : stepsize : 9;
y = -5 : stepsize : 5;
Z_surf = zeros(length(x), length(y));
for i = 1 : length(x)
    for j = 1 : length(y)
        Z_surf(j, i) = 1/(2*c)*(y(j)^2 - (max([0, y(j) - c*x(i)]))^2);
    end
end
% mesh and colour
[X, Y] = meshgrid(x, y);
C = 1.*Z_surf;
figure(2)
sfc_handle = surfc(X, Y, Z_surf, C, 'FaceAlpha',0.5, 'EdgeColor', 'none');
contourProperty = sfc_handle(2);
contourProperty.ContourZLevel = -max(max(Z_surf));
contourProperty.LineWidth = 1;
contourProperty.LevelList = [1, 5, 10];
contourProperty.ShowText = 'on';
hold on
colorbar
% contour
% surf
% fsurf
lighting gouraud;% best lighting algoithm for curved surfaces
material shiny
% sets the axis limits equal to the range of the data
axis([x(1)-0.5, x(end)+0.5, y(1)-0.5, y(end)+0.5, contourProperty.ContourZLevel, max(max(Z_surf)) + 1] )
% title('generalized primal Gap Function')
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
zlabel('$\varphi_{Au}$', 'Interpreter','latex')
end

function plotDGapConstraintSet(a, b, s)
% parameter
stepsize = 0.01;

% node
origin_x = 0;
origin_y = 0;
nodePoint1_x = sqrt(2*s/(b-a));
nodePoint1_y = b*nodePoint1_x;
nodePoint2_x = sqrt(2*b*s/(a*(b-a)));
nodePoint2_y = sqrt(2*a*b*s/(b-a));
nodePoint3_x = -sqrt(2*s/(b-a));
nodePoint3_y = a*nodePoint3_x;
nodePoint4_x = -sqrt(2*a*s/(b*(b-a)));
nodePoint4_y = -sqrt(2*a*b*s/(b-a));

% figure boundary
x_lb = nodePoint3_x - 1;
x_ub = nodePoint2_x + 1;
y_lb = nodePoint4_y - 1;
y_ub = nodePoint1_y + 1;

% region 1
reg_1_x = [nodePoint3_x, nodePoint3_x, origin_x, nodePoint1_x, nodePoint1_x];
reg_1_y = [y_ub, nodePoint3_y, origin_y, nodePoint1_y, y_ub];

% region 2
reg_2_curve_x = nodePoint1_x: stepsize: nodePoint2_x;
reg_2_curve_y = zeros(1, length(reg_2_curve_x));
for i = 1 : length(reg_2_curve_x)
    reg_2_curve_y(i) = b*reg_2_curve_x(i) - b*sqrt((1-a/b)*(reg_2_curve_x(i))^2 - 2*s/b);
end
reg_2_x = [reg_2_curve_x, origin_x];
reg_2_y = [reg_2_curve_y, origin_y];
% region 3
reg_3_x = [x_ub, nodePoint2_x, origin_x, nodePoint4_x, x_ub];
reg_3_y = [nodePoint2_y, nodePoint2_y, origin_y, nodePoint4_y, nodePoint4_y];
% region 4
reg_4_curve_x = nodePoint3_x : stepsize : nodePoint4_x;
reg_4_curve_y = zeros(1, length(reg_4_curve_x));
for i = 1 : length(reg_4_curve_x)
    reg_4_curve_y(i) = a * reg_4_curve_x(i) - a * sqrt((1 - b/a)*(reg_4_curve_x(i))^2 + 2*s/a);
end
reg_4_x = [reg_4_curve_x, origin_x];
reg_4_y = [reg_4_curve_y, origin_y];

% axis
bound_y_ax_x = x_lb: stepsize : x_ub;
bound_y_ax_y = a.* bound_y_ax_x;
bound_y_bx_x = x_lb: stepsize : x_ub;
bound_y_bx_y = b.* bound_y_bx_x;

% plot
figure(2)
patch(reg_1_x, reg_1_y, 'blue', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_2_x, reg_2_y, 'red', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_3_x, reg_3_y, 'green', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
patch(reg_4_x, reg_4_y, 'yellow', 'FaceAlpha', 0.1, 'LineStyle', 'none')
hold on
plot([nodePoint3_x, nodePoint3_x], [y_ub, nodePoint3_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint1_x, nodePoint1_x], [y_ub, nodePoint1_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint2_x, x_ub], [nodePoint2_y, nodePoint2_y], 'k', 'LineWidth', 3)
hold on
plot([nodePoint4_x, x_ub], [nodePoint4_y, nodePoint4_y], 'k', 'LineWidth', 3)
hold on
f2 = @(x,y) -a./2 .* x.^2 + x .* y - 1./(2.* b) .* y.^2 - s;
fimplicit(f2, [nodePoint1_x, nodePoint2_x, nodePoint2_y, nodePoint1_y], 'k', 'LineWidth', 3) % region 2 curve
hold on
f4 = @(x,y) b./2 .* x.^2 - x .* y + 1./(2.* a) .* y.^2 - s;
fimplicit(f4, [nodePoint3_x, nodePoint4_x, nodePoint4_y, nodePoint3_y], 'k', 'LineWidth', 3) % region 4 curve
hold on
plot(bound_y_ax_x, bound_y_ax_y, 'k', 'LineStyle', '--', 'LineWidth', 1 )
hold on
plot(bound_y_bx_x, bound_y_bx_y, 'k', 'LineStyle', '--', 'LineWidth', 1 )
axis([nodePoint3_x-0.6, nodePoint2_x+0.5, nodePoint4_y-0.4, nodePoint1_y+0.5])

xline(0, 'LineWidth', 1); 
yline(0, 'LineWidth', 1);
end

function plotDGapFunContour(a, b)
% parameter
stepsize = 0.01;
% data
x = -4 : stepsize : 6;
y = -4 : stepsize : 6;
Z_surf = zeros(length(x), length(y));
for i = 1 : length(x)
    for j = 1 : length(y)
        Z_surf(j, i) = (b-a)/(2*a*b)*y(j)^2 - 1/(2*a)*(max([0, y(j)-a*x(i)]))^2 + 1/(2*b)*(max([0,y(j)-b*x(i)]))^2;
    end
end
% mesh and colour
[X, Y] = meshgrid(x, y);
C = 1.*Z_surf;
figure(3)
sfc_handle = surfc(X, Y, Z_surf, C, 'FaceAlpha',0.5, 'EdgeColor', 'none');
contourProperty = sfc_handle(2);
contourProperty.ContourZLevel = -max(max(Z_surf));
contourProperty.LineWidth = 1;
contourProperty.LevelList = [1, 5, 10];
contourProperty.ShowText = 'on';
hold on

colorbar
% contour
% surf
% fsurf
lighting gouraud;% best lighting algoithm for curved surfaces
material shiny
box on
% sets the axis limits equal to the range of the data
axis([x(1)-0.5, x(end)+0.5, y(1)-0.5, y(end)+0.5, contourProperty.ContourZLevel, max(max(Z_surf)) + 5] )
% title('generalized D Gap Function')
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
zlabel('$\varphi^{ab}_{Au}$', 'Interpreter','latex')
end

