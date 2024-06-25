function [z_Opt, Info] = solve_NLP(self, z_Init, s_Init, s_End)
% solve parameterized NLP by solver combines non-interior-point method and CGMRES method
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
% Y node (z, gamma_h, gamma_c)
Y_Node = cumsum([self.NLP.Dim.z, self.NLP.Dim.h, self.NLP.Dim.c]);

% load option parameter (for parameter trajectory)
kappa_s_times = self.Option.Continuation.kappa_s_times;
sigma_Init = self.Option.Continuation.sigma_Init;
sigma_End = self.Option.Continuation.sigma_End;
kappa_sigma_times = self.Option.Continuation.kappa_sigma_times;
kappa_sigma_exp = self.Option.Continuation.kappa_sigma_exp;

% fictitious time step and stabilization parameter
dtau = self.Option.Continuation.dtau;
epsilon = 1/dtau;

% formulate parameter vector
p_Init = [s_Init; sigma_Init];
p_End = [s_End; sigma_End];

%% create record
% evaluate the max number of continuation step
p_test = p_Init;
j_Max = 0;
while true
    if all(p_test == p_End)
        break
    else
        p_trial = [kappa_s_times*p_test(1);...
            min([kappa_sigma_times*p_test(2), p_test(2)^kappa_sigma_exp])];
        p_test = max([p_trial, p_End], [], 2);
        j_Max = j_Max + 1;
    end
end
% log for each step's information
Log.param      = zeros(j_Max + 1, 2); % [s, sigma]
Log.Y_dot      = zeros(j_Max + 1, 1);
Log.GMRES_res  = zeros(j_Max + 1, 1);
Log.cost       = zeros(j_Max + 1, 1);
Log.KKT_error  = zeros(j_Max + 1, 1); 
Log.VI_nat_res = zeros(j_Max + 1, 1);
Log.time       = zeros(j_Max + 1, 1); 

%% continuation loop (j: continuation step counter, Y_j: current iterate, Y: previous iterate)
j = 0;
while true
    %% step 1: evaluate iterate
    if j == 0
        % specify parameter
        p_j = p_Init;
        % solve the first parameterized NLP
        [Y_j, Info_firstNLP] = self.solve_first_NLP(z_Init, p_j);
        terminal_status_j = Info_firstNLP.terminal_status;
        terminal_msg_j = Info_firstNLP.terminal_msg;
        time_j = Info_firstNLP.time;
        % initialize quantities
        Y_dot = zeros(Y_Node(3), 1);
        GMRES_res_j = 0;
    else
        % specify new parameter
        p_j_trial = [kappa_s_times*p(1); min([kappa_sigma_times*p(2), p(2)^kappa_sigma_exp])];
        p_j = max([p_j_trial, p_End], [], 2);
        % compute time derivative p_dot in previous fictitious time tau by finite difference 
        p_dot = (p_j - p)/dtau; 
        % compute time derivative Y_dot in previous fictitious time tau
        [Y_dot, Info_differential_equation] = self.solve_differential_equation(Y, p, p_dot, Y_dot_Init, epsilon);
        % evaluate new iterate by integrating a differential equation with explicit Euler method
        Y_j = Y + dtau * Y_dot;
        % terminal info and time
        terminal_status_j = 1;
        terminal_msg_j = ['- Solver succeeds: ', 'because a new iterate found by solving a differential equation']; 
        time_j = Info_differential_equation.time;
        GMRES_res_j = Info_differential_equation.GMRES_res;
    end

    %% step 2: record and print information of the current continuation step
    % some quantities of current iterate
    J_j = full(self.FuncObj.J(Y_j(1 : Y_Node(1), 1)));
    KKT_error_j = norm(full(self.FuncObj.KKT_residual(Y_j, p_j)), inf);
    VI_nat_res_j = self.evaluate_natural_residual(Y_j(1 : Y_Node(1), 1));
    % record
    Log.param(j + 1, :) = p_j';
    Log.Y_dot(j + 1, :) = norm(Y_dot, inf);
    Log.GMRES_res(j + 1, :) = GMRES_res_j;
    Log.cost(j + 1, :) = J_j;
    Log.KKT_error(j + 1, :) = KKT_error_j;
    Log.VI_nat_res(j + 1, :) = VI_nat_res_j;
    Log.time(j + 1, :) = time_j;
    % print
    if mod(j, 10) ==  0
        disp('---------------------------------------------------------------------------------------------------')
        headMsg = ' StepNum |    s     |   sigma  |   Y_dot  | GMRES_res |   cost   | KKT_error | VI_nat_res | time(s) ';
        disp(headMsg)
    end
    continuation_Step_Msg = ['  ',...
        num2str(j,'%10.2d'), '/', num2str(j_Max,'%10.2d'),'  | ',...
        num2str(Log.param(j + 1, 1), '%10.2e'),' | ', num2str(Log.param(j + 1, 2), '%10.2e'),' | ',...
        num2str(Log.Y_dot(j + 1),'%10.2e'), ' | ',...
        num2str(Log.GMRES_res(j + 1),'%10.2e'), '  | ',...
        num2str(Log.cost(j + 1),'%10.2e'), ' | ',...
        num2str(Log.KKT_error(j + 1), '%10.2e'), '  |  ',...
        num2str(Log.VI_nat_res(j + 1), '%10.2e'),'  | ',...
        num2str(Log.time(j + 1), '%10.4f')];
    disp(continuation_Step_Msg)

    %% step 3: check ternimation based on the current continuation step
    if terminal_status_j && (VI_nat_res_j <= self.Option.Continuation.tol.VI_nat_res) && (KKT_error_j <= self.Option.Continuation.tol.KKT_error)
        % IPOPT at this continuation step finds the optimal solution
        % satisfying the desired VI natural residual and KKT residual
        terminal_status = 1;
        terminal_msg = terminal_msg_j;
        break
    elseif ~terminal_status_j
        % IPOPT at this continuation step fails to find the optimal solution
        terminal_status = 0;
        terminal_msg = terminal_msg_j;
        break
    elseif j == j_Max
        % IPOPT still can not find the optimal solution in the final continuation step
        terminal_status = -1;
        terminal_msg = 'solver can not find the optimal solution satisfying the desired VI natural residual';
        break
    else
        % IPOPT at this homotopy iteration (not the final) finds the optimal solution, prepare for next homotopy iteration
        Y = Y_j;
        p = p_j;
        Y_dot_Init = Y_dot;
        j = j + 1;
    end

end

%% return optimal solution and create information
% extract primal and dual variable
z_j       = Y_j(            1 : Y_Node(1), 1);
gamma_h_j = Y_j(Y_Node(1) + 1 : Y_Node(2), 1);
gamma_c_j = Y_j(Y_Node(2) + 1 : Y_Node(3), 1);
% return the current homotopy iterate as the optimal solution
z_Opt = z_j;
% create Info
Info.continuationStepNum = j;
Info.terminal_status = terminal_status;
Info.terminal_msg = terminal_msg;
Info.gamma_h = gamma_h_j;
Info.gamma_c = gamma_c_j;
Info.cost = J_j;
Info.KKT_error = KKT_error_j;
Info.VI_natural_residual = VI_nat_res_j;
Info.Log = Log;
Info.time = sum(Log.time);
% display result
disp('*------------------------------------------------- Solution Information --------------------------------------------------*')
disp(['1. Terminal Message: ', Info.terminal_msg])
disp('2. Continuation Step Message')
disp(['- TimeElapsed: .............................................. ', num2str(Info.time,'%10.4f'), ' s'])
disp(['- Continuation Step: ........................................ ', num2str(Info.continuationStepNum)])
disp(['- Total time for solving first parameterized NLP: ........... ', num2str(Log.time(1),'%10.4f'), ' s'])
disp(['- Total time for solving differential equation: ............. ', num2str((Info.time - Log.time(1)),'%10.4f'), ' s'])
if Info.continuationStepNum ~= 0
    disp(['- Average time for solving each differential equation: ...... ', num2str((Info.time - Log.time(1))/Info.continuationStepNum,'%10.4f'), ' s'])
end
disp('3. Solution Message')
disp(['- Cost: ..................................................... ', num2str(Info.cost,'%10.3e'), '; '])
disp(['- KKT error: ................................................ ', num2str(Info.KKT_error,'%10.3e'), '; '])
disp(['- VI natural residual: ...................................... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])

end