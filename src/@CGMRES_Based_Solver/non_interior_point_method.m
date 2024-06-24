function [z_Opt, Info] = non_interior_point_method(self, z_Init, p)
%Stage 1: solve NLP with given z_Init and p by non-interior-point method
% NLP has the form:
%  min  J(z),
%  s.t. h(z) = 0,
%       c(z, s) >= 0,
% where: z is the variable,
%        s is the parameter,
%        J is the cost, and h, c are the constraints
% Syntax:
%          [z_Opt, Info] = non_interior_point_method(self, z_Init, p)
%          [z_Opt, Info] = self.non_interior_point_method(z_Init, p)
% Argument:
%          z_Init: double, NLP.Dim.z X 1, initial guess
%          p: double, given parameter, p = [s; sigma]
% Output:
%          z_Opt: double, NLP.Dim.z X 1, optimal solution found by solver
%          Info: struct, record the iteration information
%% iteration routine (Y: previous iterate Y_{k-1}, Y_k: current iterate Y_{k}) 
% time record
Time = struct('KKTEval', 0, 'searchDirection', 0, 'lineSearch', 0, 'else', 0, 'total', 0);
% counter, beta and iterate
k = 1;
beta = self.Option.NIP.LineSearch.beta_Init;
Y = [z_Init; ones(self.NLP.Dim.h, 1); ones(self.NLP.Dim.c, 1)];
Y_Node = cumsum([self.NLP.Dim.z, self.NLP.Dim.h, self.NLP.Dim.c]);

% iteration loop
while true
    %% step 0: check iteration counter
    if k > self.Option.NIP.maxIterNum
        % failure case 1: exceed the maximum number of iteration
        terminalStatus = 0;
        terminalMsg = ['- Solver fails: ', 'because the maximum number of iteration exceeded'];
        break
    end
    timeStart_total = tic;

    %% step 1: KKT residual, matrix, and error of previous iterate z
    timeStart_KKTEval = tic;
    % KKT residual and matrix
    KKT_residual = full(self.FuncObj.KKT_residual(Y, p));
    KKT_matrix = sparse(self.FuncObj.KKT_matrix(Y, p));    
    % KKT error
    KKT_error = norm(KKT_residual, inf); 
    timeElasped_KKTEval = toc(timeStart_KKTEval);

    %% step 2: search direction evaluation based on previous iterate z
    timeStart_SearchDirection = tic;
    % solve sparse linear system using mldivide
    dY = KKT_matrix\(-KKT_residual);
    % dYNorm (L_inf norm)
    dYNorm = norm(dY, inf);    
    timeElasped_searchDirection = toc(timeStart_SearchDirection);

    %% step 3: check whether we can terminate successfully based on the previous iterate z
    if KKT_error < self.Option.NIP.tol.KKT_error
        % Success case 1: the KKT error satisfies tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because the KKT error satisfies tolerance'];
        break
    elseif dYNorm < self.Option.NIP.tol.dYNorm
        % Success case 2: the norm of search direction satisfies tolerance
        terminalStatus = 1;
        terminalMsg = ['- Solver succeeds: ', 'because the norm of search direction satisfies tolerance'];   
        break       
    end

    %% step 4: merit line search
    [Y_k, Info_LineSearch] = self.LineSearch_Merit(beta, Y, p, dY);
    % check status
    if Info_LineSearch.status == 0
        % failure case 2: line search fails
        terminalStatus = -1;
        terminalMsg = ['- Solver fails: ', 'because merit line search reaches the min stepsize'];        
        break
    else
        % line search quantities
        beta_k   = Info_LineSearch.beta;
        stepSize = Info_LineSearch.stepSize;
        merit    = Info_LineSearch.merit;
    end
    timeElasped_lineSearch = Info_LineSearch.time;

    %% step 5: record and print information of this iteration k
    timeElasped_total = toc(timeStart_total);
    Time.KKTEval = Time.KKTEval + timeElasped_KKTEval;
    Time.searchDirection = Time.searchDirection + timeElasped_searchDirection;
    Time.lineSearch = Time.lineSearch + timeElasped_lineSearch;
    Time.total = Time.total + timeElasped_total;
    % some other quantities
    J = full(self.FuncObj.J(Y(1 : Y_Node(1), 1)));
    LAG_grad_T = KKT_residual(            1 : Y_Node(1), 1);
    h          = KKT_residual(Y_Node(1) + 1 : Y_Node(2), 1);
    PSI        = KKT_residual(Y_Node(2) + 1 : Y_Node(3), 1);  
    if self.Option.NIP.printLevel == 2
        % head
        if mod(k, 10) == 1
            disp('----------------------------------------------------------------------------------------------------------------------')
            headMsg = ' Iter |   cost   | LAG_grad |    h     |    PSI   |  dYNorm  |   beta   | stepsize |  merit   | merit(t) | time(ms) |';
            disp(headMsg)
        end
        % previous iterate message
        prevIterMsg = ['  ',...
            num2str(k,'%10.3d'),' | ',...
            num2str(J,'%10.2e'), ' | ',...
            num2str(norm(LAG_grad_T, inf), '%10.2e'), ' | ',...
            num2str(norm(h, inf), '%10.2e'), ' | ',...
            num2str(norm(PSI, inf), '%10.2e'), ' | ',...
            num2str(dYNorm,'%10.2e'), ' | ',...
            num2str(beta_k,'%10.2e'), ' | ',...
            num2str(stepSize,'%10.2e'), ' | ', ...
            num2str(merit(1),'%10.2e'), ' | ', num2str(merit(2),'%10.2e'), ' | ',...
            num2str(1000 * timeElasped_total,'%10.2e'), ' | '];
        disp(prevIterMsg)
    end

    %% step 6: prepare next iteration
    k = k + 1;
    beta = beta_k;
    Y = Y_k;

end

%% return optimal solution and create information
% extract primal and dual variable
z       = Y(            1 : Y_Node(1), 1);
gamma_h = Y(Y_Node(1) + 1 : Y_Node(2), 1);
gamma_c = Y(Y_Node(2) + 1 : Y_Node(3), 1);
% return previous iterate as solution
z_Opt = z;
% create Info (basic: time, iterNum, terminalStatus)
Time.else = Time.total - Time.searchDirection - Time.lineSearch - Time.KKTEval;
Info.Time = Time;
Info.iterNum = k - 1;
Info.terminalStatus = terminalStatus;
Info.terminalMsg = terminalMsg;
% create Info (corresponds to the solution z_Opt: dual variable, cost, KKT, natural residual)
Info.gamma_h = gamma_h;
Info.gamma_c = gamma_c;
Info.cost = J;
Info.KKT_error = KKT_error;
Info.VI_natural_residual = self.evaluate_natural_residual(z);
% display termination and solution message, then break rountie
if (self.Option.NIP.printLevel == 1) || (self.Option.NIP.printLevel == 2)
    disp('*--------------------------------------------- Solution Information ----------------------------------------------*')
    disp('1. Terminal Status')
    disp(Info.terminalMsg)
    disp('2. Iteration Process Message')
    disp(['- Iterations: ................... ', num2str(Info.iterNum)])
    disp(['- TimeElapsed: .................. ', num2str(Info.Time.total,'%10.3f'), 's'])
    disp(['- AverageTime: .................. ', num2str(1000 * Info.Time.total /Info.iterNum,'%10.2f'), ' ms/Iter'])
    disp('3. Solution Message')
    disp(['- Cost: ......................... ', num2str(Info.cost,'%10.3e'), '; '])
    disp(['- KKT error: .................... ', num2str(Info.KKT_error, '%10.3e'), '; '])
    disp(['- VI natural residual: .......... ', num2str(Info.VI_natural_residual,'%10.3e'), '; '])
end

end