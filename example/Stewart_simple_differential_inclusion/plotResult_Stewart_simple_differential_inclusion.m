function plotResult_Stewart_simple_differential_inclusion(OCPEC, NLP, z_Opt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Z_Opt = reshape(z_Opt, NLP.Dim.z_Node(end), OCPEC.nStages);

X_Opt = Z_Opt(1 : NLP.Dim.z_Node(1), :);
U_Opt = Z_Opt(NLP.Dim.z_Node(1) + 1 : NLP.Dim.z_Node(2), :);
LAMBDA_Opt = Z_Opt(NLP.Dim.z_Node(2) + 1 : NLP.Dim.z_Node(3), :);

timeAxis = 0 : OCPEC.timeStep : OCPEC.nStages * OCPEC.timeStep;

F_FuncObj_map = OCPEC.FuncObj.F.map(OCPEC.nStages);
F_Opt = full(F_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));
f_FuncObj_map = OCPEC.FuncObj.f.map(OCPEC.nStages);
f_Opt = full(f_FuncObj_map(X_Opt, U_Opt, LAMBDA_Opt));

figure(112)
subplot(3,1,1)
plot(timeAxis, [OCPEC.x0(1), X_Opt(1, :)], 'r', 'LineWidth',1.2)
legend('x')
xlabel('time [s]')
title('system state')

subplot(3,1,2)
plot(timeAxis(2:end), f_Opt(1, :), 'b', 'LineWidth', 1.2)
legend('f') 
xlabel('time(s)')
title('state equation')

subplot(3,1,3)
plot(timeAxis(2:end), F_Opt(1, :), 'k',...
    timeAxis(2:end), LAMBDA_Opt(1, :), 'b', 'LineWidth', 1.2)
legend('F', '\lambda') 
xlabel('time [s]')
title('check VI')
end