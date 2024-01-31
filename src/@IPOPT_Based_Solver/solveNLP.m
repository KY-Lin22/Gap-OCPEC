function [z_Opt, Info] = solveNLP(self, z_Init, p_Init, p_End)
% solve NLP with given z_Init, p_Init, and p_End by IPOPT using continuation method
% NLP has the form:
%  min  J(z, p),
%  s.t. h(z, p) = 0,
%       c(z, p) >= 0,
% where: z is the variable,
%        p is the parameter,
%        J is the cost, and h, c are the constraints
% Syntax:
%          [z_Opt, Info] = solveNLP(self, z_Init, p_Init, p_End)
%          [z_Opt, Info] = self.solveNLP(z_Init, p_Init, p_End)
% Argument:
%          z_Init: double, NLP.Dim.z X 1, initial guess
%          p_Init: double, problem parameter (initial) p = [s; mu]
%          p_End: double, problem parameter (end) p = [s; mu]
% Output:
%          z_Opt: double, NLP.Dim.z X 1, optimal solution found by solver
%          Info: struct, record the iteration information

import casadi.*

%% check option and input
% check option (TODO)

NLP = self.NLP;
Option = self.Option;

% check input z_Init
if ~all(size(z_Init) == [NLP.Dim.z, 1])
    error('z_Init has wrong dimension')
end
% check relaxation parameter
if (p_Init(1) < 0) || (p_End(1) < 0)
    error('relax parameter s (i.e., p_1) should be nonnegative')
end
if p_Init(1) < p_End(1)
    error('s_Init should not smaller than s_End')
end
% check penalty parameter
if (p_Init(2) < 0) || (p_End(2) < 0)
    error('penalty parameter mu (i.e., p_2) should be nonnegative')
end
if p_Init(2) > p_End(2)
    error('mu_Init should not larger than mu_End')
end

% load parameter
kappa_s_times = Option.Homotopy.kappa_s_times;
kappa_s_exp = Option.Homotopy.kappa_s_exp;

%% create record for time and log 
% evaluate the number of continuation step based on s_Init and s_End
s_Init = p_Init(1);
s_End = p_End(1);
s_test = s_Init;
continuationStepNum = 1;
while true
    if s_test == s_End
        break
    else        
        s_trial = min([kappa_s_times * s_test, s_test^kappa_s_exp], [], 2);
        s_test = max([s_trial, s_End], [], 2);
        continuationStepNum = continuationStepNum + 1;
    end
end

% log (time, param, cost, KKT error, stepSize, natRes)
Log.param               = zeros(continuationStepNum, 2); % [s, mu] 
Log.cost                = zeros(continuationStepNum, 2); % [ocp, penalty]
Log.KKT_error           = zeros(continuationStepNum, 2); % [primal, dual_scaled]
Log.stepSize_primal     = zeros(continuationStepNum, 2); % [min, average]
Log.stepSize_dual       = zeros(continuationStepNum, 2); % [min, average]
Log.VI_natural_residual = zeros(continuationStepNum, 1);
Log.iterNum             = zeros(continuationStepNum, 1);
Log.timeElapsed         = zeros(continuationStepNum, 1); % elapsed time in each continuation step

%% continuation loop (j: continuation step counter)
z_Init_j = z_Init;
p_j = p_Init; 
mu_End = p_End(2);

for j = 1 : continuationStepNum
    %% step 1: solve a NLP with given p
    % solve problem
    solution_j = self.Solver('x0', z_Init_j, 'p', p_j,...
        'lbg', [zeros(NLP.Dim.h, 1); zeros(NLP.Dim.c, 1)],...
        'ubg', [zeros(NLP.Dim.h, 1); inf*ones(NLP.Dim.c, 1)]);
    % extract solution and information
    z_Opt_j = full(solution_j.x);
    J_ocp_j = full(NLP.FuncObj.J_ocp(z_Opt_j, p_j));
    J_penalty_j = full(NLP.FuncObj.J_penalty(z_Opt_j, p_j));
    KKT_error_primal_j = self.Solver.stats.iterations.inf_pr(end);
    KKT_error_dual_j = self.Solver.stats.iterations.inf_du(end); 
    VI_nat_res_j = self.evaluateNaturalResidual(z_Opt_j);

    %% step 2: record and print information of the current continuation iterate
    Log.param(j, :) = [p_j(1), p_j(2)];
    Log.cost(j, :) = [J_ocp_j, J_penalty_j];
    Log.KKT_error(j, :) = [KKT_error_primal_j, KKT_error_dual_j];
    Log.stepSize_primal(j, :) = [min(self.Solver.stats.iterations.alpha_pr(2:end)),...
        sum(self.Solver.stats.iterations.alpha_pr(2:end))/(self.Solver.stats.iter_count)];
    Log.stepSize_dual(j, :) = [min(self.Solver.stats.iterations.alpha_du(2:end)),...
        sum(self.Solver.stats.iterations.alpha_du(2:end))/(self.Solver.stats.iter_count)];
    Log.VI_natural_residual(j) = VI_nat_res_j;
    Log.iterNum(j) = self.Solver.stats.iter_count;
    Log.timeElapsed(j) = self.Solver.stats.t_wall_total; % self.Solver.t_proc_total;
    if mod(j, 10) == 1
        disp('---------------------------------------------------------------------------------------------------------------------------------')
        headMsg = ' step  |   param(s/mu)   | cost(ocp/penalty) | KKT(primal/dual)| alpha_p(min/ave)| alpha_d(min/ave)| nat_res | iterNum | time(s) ';
        disp(headMsg)
    end
    prevIterMsg = [' ',...
        num2str(j,'%10.2d'), '/', num2str(continuationStepNum,'%10.2d'),' | ',...
        num2str(Log.param(j, 1), '%10.1e'), ' ', num2str(Log.param(j, 2), '%10.1e'), ' | ',...
        num2str(Log.cost(j, 1), '%10.2e'), ' ', num2str(Log.cost(j, 2), '%10.2e'),' | ',...
        num2str(Log.KKT_error(j, 1), '%10.1e'), ' ', num2str(Log.KKT_error(j, 2), '%10.1e'),' | ',...
        num2str(Log.stepSize_primal(j, 1), '%10.1e'), ' ', num2str(Log.stepSize_primal(j, 2), '%10.1e'), ' | ',...
        num2str(Log.stepSize_dual(j, 1), '%10.1e'), ' ' , num2str(Log.stepSize_dual(j, 2), '%10.1e'),' | ',...
        num2str(Log.VI_natural_residual(j), '%10.1e'),' |   ',...
        num2str(Log.iterNum(j), '%10.4d'),'  | ',...
        num2str(Log.timeElapsed(j), '%10.4f')];
    disp(prevIterMsg)
    
    %% step 3: check ternimation based on the current homotopy iterate
    if strcmp(self.Solver.stats.return_status, 'Solve_Succeeded') && (j == continuationStepNum)
        % IPOPT at the final homotopy iteration finds the optimal solution
        exitFlag = true;
    elseif ~strcmp(self.Solver.stats.return_status, 'Solve_Succeeded')
        % IPOPT at this homotopy iteration fails to find the optimal solution
        exitFlag = true;
    else
        % IPOPT at this homotopy iteration (not the final) finds the optimal solution, prepare for next homotopy iteration
        exitFlag = false;
        % update initial guess
        z_Init_j = z_Opt_j;
        % update relaxation parameter
        s_j = p_j(1);
        s_trial = min([kappa_s_times .* s_j, s_j.^kappa_s_exp]);
        s_j = max([s_trial, s_End]);
        % update penalty parameter
        mu_j = p_j(2);
        mu_trial = max([1/s_j, mu_j]);
        mu_j = min([mu_trial, mu_End]);
        % update parameter vector
        p_j(1) = s_j;
        p_j(2) = mu_j;
        
    end
    
    %% step 4: check exitFlag and return optimal solution
    if exitFlag
        % return the current homotopy iterate as the optimal solution
        z_Opt = z_Opt_j;
        % create Info
        Info.continuationStepNum = continuationStepNum;
        Info.terminalMsg = self.Solver.stats.return_status;        
        Info.cost.ocp = J_ocp_j;
        Info.cost.penalty = J_penalty_j;
        Info.KKT_error.primal = KKT_error_primal_j;
        Info.KKT_error.dual = KKT_error_dual_j;
        Info.VI_natural_residual = VI_nat_res_j;
        Info.time = sum(Log.timeElapsed);
        Info.iterNum = sum(Log.iterNum);
        Info.Log = Log;
        % display homotopy terminal result and then break        
        disp('*--------------------------------------------- Solution Information ----------------------------------------------*')
        disp(['1. Terminal Status: ', Info.terminalMsg]) 
        disp('2. Continuation Step Message')
        disp(['- TimeElapsed: ................................. ', num2str(Info.time,'%10.3f'), ' s'])
        disp(['- Continuation Step: ........................... ', num2str(Info.continuationStepNum)])        
        disp(['- Time Per Continuation Step: .................. ', num2str(Info.time / Info.continuationStepNum,'%10.2f'), ' s/Step'])
        disp(['- Iterations: .................................. ', num2str(Info.iterNum)])
        disp(['- Time Per Iteration: .......................... ', num2str(1000 * Info.time / Info.iterNum,'%10.2f'), ' ms/Iter'])
        disp('3. Solution Message')
        disp(['- Cost(ocp): ................................... ', num2str(Info.cost.ocp,'%10.3e'), '; '])
        disp(['- Cost(penalty): ............................... ', num2str(Info.cost.penalty,'%10.3e'), '; '])
        disp(['- KKT(primal): ................................. ', num2str(Info.KKT_error.primal,'%10.3e'), '; '])
        disp(['- KKT(dual): ................................... ', num2str(Info.KKT_error.dual,'%10.3e')  '; '])
        disp(['- equilibrium constraint(natural residual): .... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])
        break
    end
    
end

end

