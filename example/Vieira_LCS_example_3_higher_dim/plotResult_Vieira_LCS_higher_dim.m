function plotResult_Vieira_LCS_higher_dim(OCPEC, NLP, z_Opt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCPEC.nStages);

X_Opt = Z_Opt(1 : NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :);
LAMBDA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);

timeAxis = 0 : OCPEC.timeStep : OCPEC.nStages * OCPEC.timeStep;

F_FuncObj_map = OCPEC.FuncObj.F.map(OCPEC.nStages);
F_Opt = full(F_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));

figure(111)
subplot(2,2,1)
plot(timeAxis, [OCPEC.x0(1), X_Opt(1, :)], 'r',...
     timeAxis, [OCPEC.x0(2), X_Opt(2, :)], 'g', 'LineWidth',1.2)
legend('x_1', 'x_2')
xlabel('time [s]')
title('differential state')

subplot(2,2,2)
plot(timeAxis(2:end), LAMBDA_Opt(1, :), 'r',...
    timeAxis(2:end), LAMBDA_Opt(2, :), 'g','LineWidth', 1.2, 'LineStyle', 'none', 'Marker', '.')
legend('\lambda_1', '\lambda_2')
xlabel('time [s]')
title('algebraic variable')

subplot(2,2,3)
plot(timeAxis(2:end), U_Opt(1, :), 'r',...
    timeAxis(2:end), U_Opt(2, :), 'g','LineWidth', 1.2, 'LineStyle', 'none', 'Marker', '.')
legend('u_1', 'u_2')
xlabel('time [s]')
title('control input')

subplot(2,2,4)
plot(timeAxis(2:end), F_Opt(1, :), 'r',...
    timeAxis(2:end), F_Opt(2, :), 'g','LineWidth', 1.2, 'LineStyle', 'none', 'Marker', '.')
legend('F_1', 'F_2')
xlabel('time [s]')
title('VI function')

end