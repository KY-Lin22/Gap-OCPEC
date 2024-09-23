%% counter of generalized primal gap function for BVI
clear all
clc

c = 1;
bl = -1;
bu = 1;
counter_level = [1, 5, 10];
plot_primal_gap_fun_contour(c, bl, bu, counter_level)

%% counter of generalized D gap function for BVI
clear all
clc

a = 0.5;
b = 2;
bl = -1;
bu = 1;
counter_level = [1, 5, 10];
plot_D_gap_fun_contour(a, b, bl, bu, counter_level)


%% sub function
function plot_primal_gap_fun_contour(c, bl, bu, counter_level)
% parameter
stepsize = 0.01;
% data
x = -6 : stepsize : 6;
y = -6 : stepsize : 6;
Z_surf = zeros(length(x), length(y));
for i = 1 : length(x)
    for j = 1 : length(y)
        omega_i_j = min([max([bl, x(i)-1/c*y(j)]), bu]);
        Z_surf(j, i) = y(j)*(x(i) - omega_i_j) - c/2*(x(i) - omega_i_j)^2;
    end
end
% mesh and colour
[X, Y] = meshgrid(x, y);
C = 1.*Z_surf;
figure(1)
sfc_handle = surfc(X, Y, Z_surf, C, 'FaceAlpha',0.5, 'EdgeColor', 'none');

contourProperty = sfc_handle(2);
contourProperty.ContourZLevel = -15;
contourProperty.LineWidth = 1.5;
contourProperty.LevelList = counter_level;
contourProperty.ShowText = 'on';
hold on

colorbar

lighting gouraud;% best lighting algoithm for curved surfaces
material shiny
box on
% sets the axis limits equal to the range of the data
axis([x(1)-0.5, x(end)+0.5, y(1)-0.5, y(end)+0.5, contourProperty.ContourZLevel, 15] )
set(gca, 'FontSize', 20)
xlabel('$\lambda$', 'Interpreter','latex', 'FontSize', 20)
ylabel('$\eta$', 'Interpreter','latex', 'FontSize', 20)
zlabel('$\varphi_{Au}$', 'Interpreter','latex', 'FontSize', 20)
end

function plot_D_gap_fun_contour(a, b, bl, bu, counter_level)
% parameter
stepsize = 0.01;
% data
x = -6 : stepsize : 6;
y = -6 : stepsize : 6;
Z_surf = zeros(length(x), length(y));
for i = 1 : length(x)
    for j = 1 : length(y)
        omega_a_i_j = min([max([bl, x(i)-1/a*y(j)]), bu]);
        omega_b_i_j = min([max([bl, x(i)-1/b*y(j)]), bu]);
        Z_surf(j, i) = y(j)*(omega_b_i_j - omega_a_i_j) ...
            - a/2*(x(i) - omega_a_i_j)^2 ...
            +b/2*(x(i) - omega_b_i_j)^2;
    end
end
% mesh and colour
[X, Y] = meshgrid(x, y);
C = 1.*Z_surf;
figure(2)
sfc_handle = surfc(X, Y, Z_surf, C, 'FaceAlpha',0.5, 'EdgeColor', 'none');

contourProperty = sfc_handle(2);
contourProperty.ContourZLevel = -15;
contourProperty.LineWidth = 1.5;
contourProperty.LevelList = counter_level;
contourProperty.ShowText = 'on';
hold on

colorbar
lighting gouraud;% best lighting algoithm for curved surfaces
material shiny
box on
% sets the axis limits equal to the range of the data
axis([x(1)-0.5, x(end)+0.5, y(1)-0.5, y(end)+0.5, contourProperty.ContourZLevel, 15] )
set(gca, 'FontSize', 20)
xlabel('$\lambda$', 'Interpreter','latex', 'FontSize', 20)
ylabel('$\eta$', 'Interpreter','latex', 'FontSize', 20)
zlabel('$\varphi^{ab}_{Au}$', 'Interpreter','latex', 'FontSize', 20)
end