function [z_Opt, Info] = solve_NLP(self, z_Init, s_Init, s_End)
% solve parameterized NLP by solver using dynamical system method
% NLP has the form:
%  min  J(z),
%  s.t. h(z) = 0,
%       c(z, s) >= 0,
% where: z is the variable,
%        s is the parameter,
%        J is the cost, and h, c are the constraints
% Syntax:
%          [z_Opt, Info] = solve_NLP(self, z_Init, s_Init, s_End)
%          [z_Opt, Info] = self.solve_NLP(z_Init, s_Init, s_End)
% Argument:
%          z_Init: double, NLP.Dim.z X 1, initial guess
%          s_Init: double, relaxation parameter (initial)
%          s_End: double, relaxation parameter (end)
% Output:
%          z_Opt: double, NLP.Dim.z X 1, optimal solution found by solver
%          Info: struct, record the iteration information
import casadi.*

%% check input
% check input z_Init
if ~all(size(z_Init) == [self.NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end
% check relaxation parameter 
if (~isscalar(s_Init)) || (s_Init < 0)
    error('s_Init should be a nonnegative scalar')
end
if (~isscalar(s_End)) || (s_End < 0)
    error('s_End should be a nonnegative scalar')
end
if s_Init < s_End
    error('s_Init should not smaller than s_End')
end

%% Initialization
% Y node (z, gamma_h, gamma_c)
Y_Node = cumsum([self.NLP.Dim.z, self.NLP.Dim.h, self.NLP.Dim.c]);
% create parameter sequence
[P, P_dot, l_Max] = self.create_parameter_sequence(s_Init, s_End);
% create record
Log.param      = zeros(l_Max + 1, 2); % [s, sigma]
Log.p_dot      = zeros(l_Max + 1, 1);
Log.Y_dot      = zeros(l_Max + 1, 1);
Log.GMRES_res  = zeros(l_Max + 1, 1);
Log.cost       = zeros(l_Max + 1, 1);
Log.KKT_res    = zeros(l_Max + 1, 1);
Log.KKT_error  = zeros(l_Max + 1, 1); 
Log.VI_nat_res = zeros(l_Max + 1, 1);
Log.time       = zeros(l_Max + 1, 1); 

%% continuation loop (l: continuation step counter, Y_l: current iterate, Y: previous iterate)
Y = zeros(Y_Node(3), 1);
Y_dot = zeros(Y_Node(3), 1);
p = zeros(2, 1);
p_dot = zeros(2, 1);
l = 0;
while true
    %% step 1: evaluate iterate at current continuation step
    % specify parameter p_l and its time derivative p_dot_l
    p_l = P(:, l + 1);    
    p_dot_l = P_dot(:, l + 1);
    % evaluate iterate Y_l
    if l == 0
        % evaluate first iterate by solving first parameterized NLP
        [Y_l, Info_firstNLP] = self.solve_first_NLP(z_Init, p_l);
        terminal_status_l = Info_firstNLP.terminal_status;
        terminal_msg_l = Info_firstNLP.terminal_msg;
        timeElasped_Y = Info_firstNLP.time;
    else
        % evaluate new iterate by integrating a differential equation
        [Y_l, Info_integrator] = self.integrate_differential_equation(Y, Y_dot, p, p_dot, p_l, p_dot_l);        
        terminal_status_l = Info_integrator.terminal_status;
        terminal_msg_l = Info_integrator.terminal_msg; 
        timeElasped_Y = Info_integrator.time;
    end
    % evaluate time derivative of iterate Y_dot_l
    Y_dot_l_Init = Y_dot; % initial guess for iterate solver (GMRES)
    [Y_dot_l, Info_solve_Y_dot] = self.solve_differential_equation(Y_l, p_l, p_dot_l, Y_dot_l_Init);
    GMRES_res_l = Info_solve_Y_dot.GMRES_res;
    timeElasped_Y_dot = Info_solve_Y_dot.time;

    %% step 2: record and print information of the current continuation step
    % extract variable and parameter
    z_l       = Y_l(            1 : Y_Node(1), 1);
    gamma_c_l = Y_l(Y_Node(2) + 1 : Y_Node(3), 1);
    s_l       = p_l(1);
    % some quantities of current iterate
    J_l = full(self.FuncObj.J(z_l));
    c_l = full(self.FuncObj.c(z_l, s_l));
    KKT_residual_l = full(self.FuncObj.KKT_residual(Y_l, p_l));    
    LAG_grad_l = KKT_residual_l(            1 : Y_Node(1), 1);
    h_l        = KKT_residual_l(Y_Node(1) + 1 : Y_Node(2), 1); 
    KKT_res_l = norm(KKT_residual_l, inf);
    KKT_error_l = norm([...
        LAG_grad_l;...
        h_l;...
        min([zeros(self.NLP.Dim.c, 1), c_l], [], 2);...
        min([zeros(self.NLP.Dim.c, 1), gamma_c_l], [], 2);...
        c_l .* gamma_c_l], inf);
    VI_nat_res_l = self.evaluate_natural_residual(Y_l(1 : Y_Node(1), 1));
    % record
    Log.param(l + 1, :) = p_l';  
    Log.p_dot(l + 1, :) = norm(p_dot_l, inf);
    Log.Y_dot(l + 1, :) = norm(Y_dot_l, inf);
    Log.GMRES_res(l + 1, :) = GMRES_res_l;
    Log.cost(l + 1, :) = J_l;
    Log.KKT_res(l + 1, :) = KKT_res_l;
    Log.KKT_error(l + 1, :) = KKT_error_l;
    Log.VI_nat_res(l + 1, :) = VI_nat_res_l;
    Log.time(l + 1, :) = timeElasped_Y + timeElasped_Y_dot;
    % print
    if mod(l, 10) ==  0
        disp('--------------------------------------------------------------------------------------------------------------------------------')
        headMsg = '  StepNum |    s     |   sigma  |   p_dot  |   Y_dot  | GMRES_res |   cost   | KKT_res  | KKT_error | VI_nat_res | time[s] ';
        disp(headMsg)
    end
    continuation_Step_Msg = ['  ',...
        num2str(l,'%10.3d'), '/', num2str(l_Max,'%10.3d'),' | ',...
        num2str(Log.param(l + 1, 1), '%10.2e'),' | ', num2str(Log.param(l + 1, 2), '%10.2e'),' | ',...
        num2str(Log.p_dot(l + 1),'%10.2e'), ' | ',...
        num2str(Log.Y_dot(l + 1),'%10.2e'), ' | ',...
        num2str(Log.GMRES_res(l + 1),'%10.2e'), '  | ',...
        num2str(Log.cost(l + 1),'%10.2e'), ' | ',...
        num2str(Log.KKT_res(l + 1), '%10.2e'), ' | ',...
        num2str(Log.KKT_error(l + 1), '%10.2e'), '  |  ',...
        num2str(Log.VI_nat_res(l + 1), '%10.2e'),'  | ',...
        num2str(Log.time(l + 1), '%10.4f')];
    disp(continuation_Step_Msg)

    %% step 3: check ternimation based on the current continuation step
    if terminal_status_l && (VI_nat_res_l <= self.Option.Continuation.tol.VI_nat_res) && (KKT_error_l <= self.Option.Continuation.tol.KKT_error)
        % this continuation step finds the desired optimal solution
        terminal_status = 1;
        terminal_msg = terminal_msg_l;
        break
    elseif ~terminal_status_l
        % this continuation step fails to find the optimal solution
        terminal_status = 0;
        terminal_msg = terminal_msg_l;
        break
    elseif l == l_Max
        % final continuation step still can not find the desired optimal solution
        terminal_status = -1;
        terminal_msg = 'solver can not find the optimal solution satisfying the desired VI natural residual';
        break
    else
        % this continuation step finds the optimal solution, prepare for next step
        Y = Y_l;
        Y_dot = Y_dot_l;
        p = p_l;
        p_dot = p_dot_l;
        l = l + 1;
    end

end

%% return optimal solution and create information
% extract primal and dual variable
z_l       = Y_l(            1 : Y_Node(1), 1);
gamma_h_l = Y_l(Y_Node(1) + 1 : Y_Node(2), 1);
gamma_c_l = Y_l(Y_Node(2) + 1 : Y_Node(3), 1);
% return the current homotopy iterate as the optimal solution
z_Opt = z_l;
% create Info
Info.continuationStepNum = l;
Info.terminal_status = terminal_status;
Info.terminal_msg = terminal_msg;
Info.gamma_h = gamma_h_l;
Info.gamma_c = gamma_c_l;
Info.cost = J_l;
Info.KKT_error = KKT_error_l;
Info.VI_natural_residual = VI_nat_res_l;
Info.Log = Log;
Info.time = sum(Log.time);
% display result
disp('*------------------------------------------------- Solution Information --------------------------------------------------*')
disp('1. Terminal Message')
disp(Info.terminal_msg)
disp('2. Continuation Step Message')
disp(['- TimeElapsed: ................................................... ', num2str(Info.time,'%10.4f'), ' s'])
disp(['- Total time for solving first parameterized NLP: ................ ', num2str(Info.Log.time(1),'%10.4f'), ' s'])
disp(['- Total time for solving subsequent parameterized NLP: ........... ', num2str((Info.time - Info.Log.time(1)),'%10.4f'), ' s'])
disp(['- Continuation Step: ............................................. ', num2str(Info.continuationStepNum)])
if Info.continuationStepNum ~= 0
    disp(['- Average time for each step using dynamical system method: ...... ', num2str((Info.time - Info.Log.time(1))/Info.continuationStepNum,'%10.4f'), ' s'])
end
disp('3. Solution Message')
disp(['- Cost: .......................................................... ', num2str(Info.cost,'%10.3e'), '; '])
disp(['- KKT error: ..................................................... ', num2str(Info.KKT_error,'%10.3e'), '; '])
disp(['- VI natural residual: ........................................... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])

end