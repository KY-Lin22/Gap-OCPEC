function [z_Opt, Info] = solve_NLP(self, z_Init, s_Init, s_End)
% solve parameterized NLP by solver combining non-interior-point and CGMRES 
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
kappa_s_times = self.Option.Homotopy.kappa_s_times;
sigma_Init = self.Option.Homotopy.sigma_Init;
sigma_End = self.Option.Homotopy.sigma_End;
kappa_sigma_times = self.Option.Homotopy.kappa_sigma_times;
kappa_sigma_exp = self.Option.Homotopy.kappa_sigma_exp;

% load option parameter (for CGMRES method)
dtau = self.Option.CGMRES.dtau;
h_FD = self.Option.CGMRES.h_FD;
k_max = self.Option.CGMRES.k_max;
% compute stabilization parameter
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
Log.param        = zeros(j_Max + 1, 2); % [s, sigma]
Log.GMRES_res    = zeros(j_Max + 1, 1);
Log.cost         = zeros(j_Max + 1, 1);
Log.KKT_res_norm = zeros(j_Max + 1, 1); 
Log.VI_nat_res   = zeros(j_Max + 1, 1);
Log.timeElapsed  = zeros(j_Max + 1, 1); 

%% continuation loop (j: continuation step counter, Y_j: current iterate, Y: previous iterate)
j = 0;
while true
    %% step 1: evaluate iterate
    if j == 0
        % specify parameter
        p_j = p_Init; 
        % provide initial guess of Y_dot for next iteration
        Y_dot = zeros(Y_Node(3), 1);
        % solve the first parameterized NLP by non-interior-point method        
        [z_j, Info_NIP] = self.non_interior_point_method(z_Init, p_j);
        gamma_h_j = Info_NIP.gamma_h;
        gamma_c_j = Info_NIP.gamma_c;
        Y_j = [z_j; gamma_h_j; gamma_c_j];
        % terminal info and time
        GMRES_res_j = 0;
        terminalStatus_j = Info_NIP.terminalStatus;
        terminalMsg_j = Info_NIP.terminalMsg;
        time_j = Info_NIP.Time.total;
    else
        % specify parameter
        p_j_trial = [kappa_s_times*p(1);...
            min([kappa_sigma_times*p(2), p(2)^kappa_sigma_exp])];
        p_j = max([p_j_trial, p_End], [], 2);        
        % compute time derivative p_dot in previous fictitious time tau by finite difference
        p_dot = (p_j - p)/dtau;
        % compute time derivative Y_dot in previous fictitious time tau by CGMRES method 
        [Y_dot, Info_CGMRES] = self.CGMRES_method(Y, p, p_dot, Y_dot_Init, h_FD, k_max, epsilon);
        % evaluate new iterate by integrating a differential equation with explicit Euler method
        Y_j = Y + dtau * Y_dot;
        % extract primal and dual variable
        z_j       = Y_j(            1 : Y_Node(1), 1);
        gamma_h_j = Y_j(Y_Node(1) + 1 : Y_Node(2), 1);
        gamma_c_j = Y_j(Y_Node(2) + 1 : Y_Node(3), 1);
        % terminal info and time
        GMRES_res_j = Info_CGMRES.rho;
        terminalStatus_j = 1;
        terminalMsg_j = ['- Solver succeeds: ', 'because a new iterate found by CGMRES method'];
        time_j = Info_CGMRES.time;
    end

    %% step 2: evaluate iterate information: cost, KKT residual, VI natural residual
    J_j = full(self.FuncObj.J(z_j));
    KKT_Residual_j = self.evaluate_KKT_Residual(Y_j, p_j);
    KKT_res_norm_j = norm(KKT_Residual_j, inf); % option: L2 norm scaled by dim_Y
    VI_nat_res_j = self.evaluate_natural_residual(z_j);

    %% step 3: record and print information of the current continuation iterate
    % record
    Log.param(j + 1, :) = p_j';
    Log.GMRES_res(j + 1, :) = GMRES_res_j;
    Log.cost(j + 1, :) = J_j;
    Log.KKT_res_norm(j + 1, :) = KKT_res_norm_j;
    Log.VI_nat_res(j + 1, :) = VI_nat_res_j;
    Log.timeElapsed(j + 1, :) = time_j;

    % print
    if mod(j, 10) ==  0
        disp('---------------------------------------------------------------------------------------------------')
        headMsg = ' StepNum |    s     |   sigma  | GMRES_res |   cost   |  KKT_res | VI_nat_res | time(s) ';
        disp(headMsg)
    end
    continuation_Step_Msg = ['  ',...
        num2str(j,'%10.2d'), '/', num2str(j_Max,'%10.2d'),'  | ',...
        num2str(Log.param(j + 1, 1), '%10.2e'),' | ', num2str(Log.param(j + 1, 2), '%10.2e'),' | ',...
        num2str(Log.GMRES_res(j + 1),'%10.2e'), '  | ',...
        num2str(Log.cost(j + 1),'%10.2e'), ' | ',...
        num2str(Log.KKT_res_norm(j + 1), '%10.2e'), ' |  ',...
        num2str(Log.VI_nat_res(j + 1), '%10.2e'),'  | ',...
        num2str(Log.timeElapsed(j + 1), '%10.4f')];
    disp(continuation_Step_Msg)

    %% step 4: check ternimation based on the current homotopy iterate
    if terminalStatus_j && (VI_nat_res_j <= self.Option.Homotopy.VI_nat_res_tol)
        % IPOPT at this homotopy iteration finds the optimal solution satisfying the desired VI natural residual
        exitFlag = true;
        terminalStatus = 1;
        terminalMsg = terminalMsg_j;
    elseif ~terminalStatus_j
        % IPOPT at this homotopy iteration fails to find the optimal solution
        exitFlag = true;
        terminalStatus = 0;
        terminalMsg = terminalMsg_j;
    elseif j == j_Max
        % IPOPT still can not find the optimal solution in the final homotopy iteration
        exitFlag = true;
        terminalStatus = -1;
        terminalMsg = 'solver can not find the optimal solution satisfying the desired VI natural residual';
    else
        % IPOPT at this homotopy iteration (not the final) finds the optimal solution, prepare for next homotopy iteration
        exitFlag = false;
        % update Y, p, Y_dot_Init and continuation step counter
        Y = Y_j;
        p = p_j;
        Y_dot_Init = Y_dot;
        j = j + 1;
    end

    %% step 5: check exitFlag and return optimal solution
    if exitFlag
        % return the current homotopy iterate as the optimal solution
        z_Opt = z_j;
        % create Info
        Info.continuationStepNum = j;
        Info.terminalStatus = terminalStatus; 
        Info.terminalMsg = terminalMsg; 
        Info.dual_var = [gamma_h_j; gamma_c_j];
        Info.cost = J_j;
        Info.KKT_residual_norm = KKT_res_norm_j;
        Info.VI_natural_residual = VI_nat_res_j;
        Info.Log = Log;
        Info.time = sum(Log.timeElapsed);
        % display homotopy terminal result and then break        
        disp('*--------------------------------------------- Solution Information ----------------------------------------------*')
        disp(['1. Terminal Message: ', Info.terminalMsg]) 
        disp('2. Continuation Step Message')
        disp(['- TimeElapsed: ................................. ', num2str(Info.time,'%10.4f'), ' s'])
        disp(['- Continuation Step: ........................... ', num2str(Info.continuationStepNum)])        
        disp(['- Time for non interior point method: .......... ', num2str(Log.timeElapsed(1),'%10.4f'), ' s'])
        disp(['- Time for CGMRES method: ...................... ', num2str((Info.time - Log.timeElapsed(1)),'%10.4f'), ' s'])
        disp('3. Solution Message')
        disp(['- Cost: ........................................ ', num2str(Info.cost,'%10.3e'), '; '])
        disp(['- KKT residual (Norm): ......................... ', num2str(Info.KKT_residual_norm,'%10.3e'), '; '])
        disp(['- equilibrium constraint(natural residual): .... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])

        break
    end
end

end