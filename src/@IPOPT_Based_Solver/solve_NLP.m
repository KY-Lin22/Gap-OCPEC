function [z_Opt, Info] = solve_NLP(self, z_Init, s_Init, s_End)
% solve parameterized NLP by IPOPT using continuation method
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
% create parameter sequence
[S, l_Max] = self.create_parameter_sequence(s_Init, s_End);
% create record 
Log.param      = zeros(l_Max + 1, 1);
Log.cost       = zeros(l_Max + 1, 1);
Log.KKT_error  = zeros(l_Max + 1, 1); % max([primal, dual])
Log.VI_nat_res = zeros(l_Max + 1, 1);
Log.iterNum    = zeros(l_Max + 1, 1);
Log.time       = zeros(l_Max + 1, 1); 

%% continuation loop (l: continuation step counter, z_l: current iterate, z: previous iterate)
z = z_Init;
l = 0;
while true
    %% step 1: evaluate iterate at current continuation step
    % specify parameter s_l
    s_l = S(:, l + 1);
    % evaluate iterate z_l by solving a parameterized NLP using IPOPT
    solution_l = self.Solver('x0', z, 'p', s_l,...
        'lbg', [zeros(self.NLP.Dim.h, 1); zeros(self.NLP.Dim.c, 1)],...
        'ubg', [zeros(self.NLP.Dim.h, 1); inf*ones(self.NLP.Dim.c, 1)]);
    % extract solution and information
    z_l = full(solution_l.x);
    gamma_h_l = full(solution_l.lam_g(1 : self.NLP.Dim.h, 1));
    gamma_c_l = full(solution_l.lam_g(self.NLP.Dim.h + 1 : end, 1));

    J_l = full(solution_l.f);
    KKT_error_primal_l = self.Solver.stats.iterations.inf_pr(end);
    KKT_error_dual_l = self.Solver.stats.iterations.inf_du(end); 
    VI_nat_res_l = self.evaluate_natural_residual(z_l);

    stepSize_primal_l = self.Solver.stats.iterations.alpha_pr(2:end);
    stepSize_dual_l = self.Solver.stats.iterations.alpha_du(2:end); 

    iterNum_l = self.Solver.stats.iter_count;
    solver_status_l = (strcmp(self.Solver.stats.return_status, 'Solve_Succeeded')) ...
        || (strcmp(self.Solver.stats.return_status, 'Solved_To_Acceptable_Level'))...
        || (strcmp(self.Solver.stats.return_status, 'Feasible_Point_Found'));
    terminal_status_l = solver_status_l;
    terminal_msg_l = self.Solver.stats.return_status;
    time_l = self.Solver.stats.t_wall_total; % self.Solver.t_proc_total;

    %% step 2: record and print information of the current continuation step
    % record
    Log.param(l + 1, :) = s_l;
    Log.cost(l + 1, :) = J_l;
    Log.KKT_error(l + 1, :) = max([KKT_error_primal_l, KKT_error_dual_l]);
    Log.VI_nat_res(l + 1, :) = VI_nat_res_l;
    Log.iterNum(l + 1, :) = iterNum_l;
    Log.time(l + 1, :) = time_l;
    if mod(l, 10) ==  0
        disp('----------------------------------------------------------------------------------------------------------------')
        headMsg = '  stepNum |    s    |   cost   | KKT(primal/dual)| alpha_p(min/avg)| alpha_d(min/avg)| VI_nat_res | iterNum | time[s] ';
        disp(headMsg)
    end
    continuation_Step_Msg = ['  ',...
        num2str(l,'%10.3d'), '/', num2str(l_Max,'%10.3d'),' | ',...
        num2str(Log.param(l + 1, 1), '%10.1e'), ' | ',...
        num2str(Log.cost(l + 1), '%10.2e'), ' | ',...
        num2str(KKT_error_primal_l, '%10.1e'), ' ', num2str(KKT_error_dual_l, '%10.1e'),' | ',...
        num2str(min(stepSize_primal_l), '%10.1e'), ' ', num2str(sum(stepSize_primal_l)/iterNum_l, '%10.1e'), ' | ',...
        num2str(min(stepSize_dual_l), '%10.1e'), ' ', num2str(sum(stepSize_dual_l)/iterNum_l, '%10.1e'),' |  ',...
        num2str(Log.VI_nat_res(l + 1), '%10.1e'),'   |   ',...
        num2str(Log.iterNum(l + 1, :), '%10.4d'),'  | ',...
        num2str(Log.time(l + 1, :), '%10.4f')];
    disp(continuation_Step_Msg)

    %% step 3: check ternimation based on the current continuation step
    if terminal_status_l && (VI_nat_res_l <= self.Option.Continuation.tol.VI_nat_res)
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
        z = z_l;
        l = l + 1;
    end

end

%% return optimal solution and create information
% return the current homotopy iterate as the optimal solution
z_Opt = z_l;
% create Info
Info.continuationStepNum = l;
Info.terminal_status = terminal_status;
Info.terminal_msg = terminal_msg;
Info.gamma_h = gamma_h_l;
Info.gamma_c = gamma_c_l;
Info.cost = J_l;
Info.KKT_error = max([KKT_error_primal_l, KKT_error_dual_l]);
Info.VI_natural_residual = VI_nat_res_l;
Info.Log = Log;
Info.time = sum(Log.time);
Info.iterNum = sum(Log.iterNum);
% display result
disp('*------------------------------------------------- Solution Information --------------------------------------------------*')
disp('1. Terminal Message')
disp(['- ', Info.terminal_msg])
disp('2. Continuation Step Message')
disp(['- TimeElapsed: ................................................... ', num2str(Info.time,'%10.4f'), ' s'])
disp(['- Iterations: .................................................... ', num2str(Info.iterNum)])
disp(['- TimeElapsed Per Iteration: ..................................... ', num2str(1000 * Info.time / Info.iterNum,'%10.2f'), ' ms/Iter'])
disp(['- Total time for solving first parameterized NLP: ................ ', num2str(Info.Log.time(1),'%10.4f'), ' s'])
disp(['- Total time for solving subsequent parameterized NLP: ........... ', num2str((Info.time - Info.Log.time(1)),'%10.4f'), ' s'])
disp(['- Continuation Step: ............................................. ', num2str(Info.continuationStepNum)])
if Info.continuationStepNum ~= 0
    disp(['- Average time for solving each subsequent parameterized NLP: .... ', num2str((Info.time - Info.Log.time(1))/Info.continuationStepNum,'%10.4f'), ' s'])
end
disp('3. Solution Message')
disp(['- Cost: .......................................................... ', num2str(Info.cost,'%10.3e'), '; '])
disp(['- KKT error: ..................................................... ', num2str(Info.KKT_error,'%10.3e'), '; '])
disp(['- VI natural residual: ........................................... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])

end

