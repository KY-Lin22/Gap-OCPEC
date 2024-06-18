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

% load parameter
kappa_s_times = self.Option.Homotopy.kappa_s_times;
VI_nat_res_tol = self.Option.Homotopy.VI_nat_res_tol;
sigma_Init = self.Option.Homotopy.sigma_Init;
sigma_End = self.Option.Homotopy.sigma_End;
kappa_sigma_times = self.Option.Homotopy.kappa_sigma_times;
kappa_sigma_exp = self.Option.Homotopy.kappa_sigma_exp;

% formulate parameter vector
p_Init = [s_Init; sigma_Init];
p_End = [s_End; sigma_End];

%% create record for time and log 
% evaluate the max number of continuation step
p_test = p_Init;
continuationStepMaxNum = 0;
while true
    if all(p_test == p_End)
        break
    else
        p_trial = [kappa_s_times*p_test(1);...
            min([kappa_sigma_times*p_test(2), p_test(2)^kappa_sigma_exp])];
        p_test = max([p_trial, p_End], [], 2);
        continuationStepMaxNum = continuationStepMaxNum + 1;
    end
end
% log for each step's information
Log.param       = zeros(continuationStepMaxNum + 1, 2); % [s, sigma]
Log.cost        = zeros(continuationStepMaxNum + 1, 1);
Log.KKT_error   = zeros(continuationStepMaxNum + 1, 6); % [primal, dual, dual_scaled, complementary, complementary_scaled, total]
Log.VI_nat_res  = zeros(continuationStepMaxNum + 1, 1);
Log.timeElapsed = zeros(continuationStepMaxNum + 1, 1); 

%% continuation loop (j: continuation step counter)
z_Init_j = z_Init;
p_j = p_Init; 
j = 0;
while true
    %% Step 1: solve subproblem and extract solution & information
    if j == 0
        % stage 1: solve the first parameterized NLP by non-interior-point method
        [z_Opt_j, Info_Stage_1] = self.non_interior_point_method(z_Init_j, p_j);
        gamma_h_j = Info_Stage_1.gamma_h;
        gamma_c_j = Info_Stage_1.gamma_c;
        cost_j = Info_Stage_1.cost;
        KKT_error_primal_j               = Info_Stage_1.KKT_error_primal;
        KKT_error_dual_j                 = Info_Stage_1.KKT_error_dual;
        KKT_error_dual_scaled_j          = Info_Stage_1.KKT_error_dual_scaled;
        KKT_error_complementary_j        = Info_Stage_1.KKT_error_complementary;
        KKT_error_complementary_scaled_j = Info_Stage_1.KKT_error_complementary_scaled;
        KKT_error_total_j                = Info_Stage_1.KKT_error_total;
        VI_nat_res_j = Info_Stage_1.VI_natural_residual;
        time_j = Info_Stage_1.Time.total;
    else
        % stage 2: solve the subsequent parameterized NLP by C/GMRES method

    end

    %% step 2: record and print information of the current continuation iterate
    % record
    Log.param(j + 1, :) = p_j';
    Log.cost(j + 1, :) = cost_j;
    Log.KKT_error(j + 1, :) = [KKT_error_primal_j, KKT_error_dual_j, KKT_error_dual_scaled_j,...
        KKT_error_complementary_j, KKT_error_complementary_scaled_j, KKT_error_total_j];
    Log.VI_nat_res(j + 1, :) = VI_nat_res_j;
    Log.timeElapsed(j + 1, :) = time_j;

    % print
    if mod(j, 10) ==  0
        disp('---------------------------------------------------------------------------------------------------------------------------------------------')
        headMsg = ' StepNum |    s     |   sigma  |   cost   |  KKT(P)  |  KKT(D)  | KKT(D,s) |  KKT(C)  | KKT(C,s) |  KKT(T)  | VI_nat_res |  time(s) |';
        disp(headMsg)
    end
    continuation_Step_Msg = ['  ',...
        num2str(j,'%10.2d'), '/', num2str(continuationStepMaxNum,'%10.2d'),'  | ',...
        num2str(Log.param(j + 1, 1), '%10.2e'),' | ',...
        num2str(Log.param(j + 1, 2), '%10.2e'),' | ',...
        num2str(Log.cost(j + 1),'%10.2e'), ' | ',...
        num2str(Log.KKT_error(j + 1, 1), '%10.2e'), ' | ',...
        num2str(Log.KKT_error(j + 1, 2), '%10.2e'), ' | ',...
        num2str(Log.KKT_error(j + 1, 3), '%10.2e'), ' | ',...
        num2str(Log.KKT_error(j + 1, 4), '%10.2e'), ' | ',...
        num2str(Log.KKT_error(j + 1, 5), '%10.2e'), ' | ',...
        num2str(Log.KKT_error(j + 1, 6), '%10.2e'), ' |  ',...
        num2str(Log.VI_nat_res(j + 1), '%10.2e'),'  |  ',...
        num2str(Log.timeElapsed(j + 1), '%10.4f'),'  | '];
    disp(continuation_Step_Msg)

    %% step 3: check ternimation based on the current homotopy iterate

    %% step 4: check exitFlag and return optimal solution

end



end