%% counter of smooth generalized primal gap function
clear all
clc
stepSize = 0.01;
lambda = -1 : stepSize : 9;
eta = -5 : stepSize : 5;
c = 1;
eps = 0.1; % smooth parameter
plotSmoothPrimalGapFunContour(c, eps, lambda, eta)

plotSmoothPrimalGapFunGrad(c, eps, lambda, eta)


%% counter of smooth generalized D gap function
clear all
clc
stepSize = 0.01;
lambda = -4 : stepSize : 6;
eta = -4 : stepSize : 6;
a = 0.3;
b = 2;
eps = 0.1;
plotSmoothD_GapFunContour(a, b, eps, lambda, eta)

plotSmoothD_GapFunGrad(a, b, eps, lambda, eta)


%% sub function
function plotSmoothPrimalGapFunContour(c, eps, lambda, eta)
Fun_surf = zeros(length(lambda), length(eta));
for i = 1 : length(lambda)
    for j = 1 : length(eta)
        smooth_max = 0.5*(sqrt((eta(j) - c*lambda(i))^2 + 4 * eps^2) + (eta(j) - c*lambda(i)));
        Fun_surf(i, j) = 1/(2*c)*(eta(j)^2 - (smooth_max)^2);
    end
end
Fun_surf = Fun_surf';
% mesh and colour
[X, Y] = meshgrid(lambda, eta);
C = 1.*Fun_surf;
figure(1)
sfc_handle = surfc(X, Y, Fun_surf, C, 'FaceAlpha',0.5, 'EdgeColor', 'none');
contourProperty = sfc_handle(2);
contourProperty.ContourZLevel = -max(max(Fun_surf));
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
axis([lambda(1)-0.5, lambda(end)+0.5, eta(1)-0.5, eta(end)+0.5, contourProperty.ContourZLevel, max(max(Fun_surf)) + 1] )
% title('generalized primal Gap Function')
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
zlabel('$\varphi_{Au}$', 'Interpreter','latex')
end

function plotSmoothPrimalGapFunGrad(c, eps, lambda, eta)
Grad_lambda = zeros(length(lambda), length(eta));
Grad_eta = zeros(length(lambda), length(eta));

for i = 1 : length(lambda)
    for j = 1 : length(eta)
        smooth_max = 0.5*(sqrt((eta(j) - c*lambda(i))^2 + 4 * eps^2) + (eta(j) - c*lambda(i)));
        Grad_lambda(i, j) = smooth_max;
        Grad_eta (i, j) = 1/c*eta(j) - 1/c * smooth_max;
    end
end
Grad_lambda = Grad_lambda';
Grad_eta = Grad_eta';

% mesh and colour
[X, Y] = meshgrid(lambda, eta);
color_lambda = 1.*Grad_lambda;
color_eta = 1 .* Grad_eta;
% grad_lambda
figure(2)
surf(X, Y, Grad_lambda, color_lambda, 'FaceAlpha',0.5, 'EdgeColor', 'none');
hold on
colorbar
lighting gouraud;% best lighting algoithm for curved surfaces
material shiny
axis([lambda(1)-0.5, lambda(end)+0.5, eta(1)-0.5, eta(end)+0.5, -max(max(Grad_lambda)), max(max(Grad_lambda))+1])
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
zlabel('$\nabla_{\lambda} \varphi_{Au}$', 'Interpreter','latex')
% grad_eta
figure(3)
surf(X, Y, Grad_eta, color_eta, 'FaceAlpha',0.5, 'EdgeColor', 'none');
hold on
colorbar
lighting gouraud;% best lighting algoithm for curved surfaces
material shiny
axis([lambda(1)-0.5, lambda(end)+0.5, eta(1)-0.5, eta(end)+0.5, -max(max(Grad_eta)), max(max(Grad_eta))+1])
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
zlabel('$\nabla_{\eta} \varphi_{Au}$', 'Interpreter','latex')
end

function plotSmoothD_GapFunContour(a, b, eps, lambda, eta)
Fun_surf = zeros(length(lambda), length(eta));
for i = 1 : length(lambda)
    for j = 1 : length(eta)
        smooth_max_a = 0.5*(sqrt((eta(j) - a*lambda(i))^2 + 4 * eps^2) + (eta(j) - a*lambda(i)));
        smooth_max_b = 0.5*(sqrt((eta(j) - b*lambda(i))^2 + 4 * eps^2) + (eta(j) - b*lambda(i)));
        Fun_surf(i, j) = 1/(2*a)*(eta(j)^2 - (smooth_max_a)^2) - 1/(2*b)*(eta(j)^2 - (smooth_max_b)^2);
    end
end
Fun_surf = Fun_surf';
% mesh and colour
[X, Y] = meshgrid(lambda, eta);
colour_ab = 1.*Fun_surf;
figure(1)
sfc_handle = surfc(X, Y, Fun_surf, colour_ab, 'FaceAlpha',0.5, 'EdgeColor', 'none');
contourProperty = sfc_handle(2);
contourProperty.ContourZLevel = -max(max(Fun_surf));
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
axis([lambda(1)-0.5, lambda(end)+0.5, eta(1)-0.5, eta(end)+0.5, contourProperty.ContourZLevel, max(max(Fun_surf)) + 1] )
% title('generalized primal Gap Function')
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
zlabel('$\varphi^{ab}_{Au}$', 'Interpreter','latex')

end

function plotSmoothD_GapFunGrad(a, b, eps, lambda, eta)
Grad_lambda = zeros(length(lambda), length(eta));
Grad_eta = zeros(length(lambda), length(eta));


for i = 1 : length(lambda)
    for j = 1 : length(eta)
        smooth_max_a = 0.5*(sqrt((eta(j) - a*lambda(i))^2 + 4 * eps^2) + (eta(j) - a*lambda(i)));
        smooth_max_b = 0.5*(sqrt((eta(j) - b*lambda(i))^2 + 4 * eps^2) + (eta(j) - b*lambda(i)));
        Grad_lambda(i, j) = smooth_max_a - smooth_max_b;
        Grad_eta (i, j) = (1/a*eta(j) - 1/a * smooth_max_a) - (1/b*eta(j) - 1/b * smooth_max_b);
    end
end

Grad_lambda = Grad_lambda';
Grad_eta = Grad_eta';

% mesh and colour
[X, Y] = meshgrid(lambda, eta);
color_lambda = 1.*Grad_lambda;
color_eta = 1 .* Grad_eta;
% grad_lambda
figure(2)
surf(X, Y, Grad_lambda, color_lambda, 'FaceAlpha',0.5, 'EdgeColor', 'none');
hold on
colorbar
lighting gouraud;% best lighting algoithm for curved surfaces
material shiny
axis([lambda(1)-0.5, lambda(end)+0.5, eta(1)-0.5, eta(end)+0.5, -max(max(Grad_lambda)), max(max(Grad_lambda))+1])
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
zlabel('$\nabla^{ab}_{\lambda} \varphi_{Au}$', 'Interpreter','latex')
% grad_eta
figure(3)
surf(X, Y, Grad_eta, color_eta, 'FaceAlpha',0.5, 'EdgeColor', 'none');
hold on
colorbar
lighting gouraud;% best lighting algoithm for curved surfaces
material shiny
axis([lambda(1)-0.5, lambda(end)+0.5, eta(1)-0.5, eta(end)+0.5, -max(max(Grad_eta)), max(max(Grad_eta))+1])
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('$\eta$', 'Interpreter','latex')
zlabel('$\nabla^{ab}_{\eta} \varphi_{Au}$', 'Interpreter','latex')
end