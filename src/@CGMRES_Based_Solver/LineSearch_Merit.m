function [Y_k, Info] = LineSearch_Merit(self, beta, Y, p, dY)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here
timeStart = tic;

Y_Node = cumsum([self.NLP.Dim.z, self.NLP.Dim.h, self.NLP.Dim.c]);
% load parameter
stepSize_min = self.Option.NIP.LineSearch.stepSize_Min;
stepSize_decayRate = self.Option.NIP.LineSearch.stepSize_DecayRate;
rho = self.Option.NIP.LineSearch.rho;
nu_D = self.Option.NIP.LineSearch.nu_D;
timeStep = self.OCPEC.timeStep;

%% some quantities at current iterate Y
% directional derivative of cost 
J_grad_times_dz = full(self.FuncObj.J_grad_times_dz(Y, dY));
% constraint violation M (L1 norm scaled by time step)
% - L1 norm follows IPOPT, and also the cost is the sum of stage cost
% - as a constraint measure, it need to be scaled by time step to consistent with the cost that has been scaled
M = timeStep * norm(full(self.FuncObj.M(Y, p)), 1);
% penalty parameter
beta_Trial = J_grad_times_dz/((1 - rho) * M);
if beta >= beta_Trial
    beta_k = beta;
else
    beta_k = beta_Trial + 1;
end
% merit and its directional derivative
z  = Y(1 : Y_Node(1), 1);
merit = full(self.FuncObj.J(z)) + beta_k * M;
merit_DD = J_grad_times_dz - beta_k * M;

%% backtracking line search
has_found_new_iterate = false;
stepSize_init = 1;

while ~has_found_new_iterate
     %% Step 1: estimate trial stepsize, iterate, and merit
     % step size
     stepSize_trial = max([stepSize_init, stepSize_min]);
     % iterate
     Y_trial = Y + stepSize_trial * dY;     
     % merit
     z_trial = Y_trial(1 : Y_Node(1), 1);
     M_trial = timeStep * norm(full(self.FuncObj.M(Y_trial, p)), 1);
     merit_trial = full(self.FuncObj.J(z_trial)) + beta_k * M_trial;

     %% Step 2: check sufficient decrease condition
     if merit_trial <= merit + stepSize_trial * nu_D * merit_DD
         has_found_new_iterate = true;
         status = 1;
     end

     %% Step 3: checking min stepsize
    if ~has_found_new_iterate
        if stepSize_trial == stepSize_min
            % linesearch fails on the min stepsize, break backtracking linesearch procedure
            status = 0;
            break
        else
            % estimate a smaller stepsize
            stepSize_init = stepSize_decayRate * stepSize_init;
        end
    end

end

%% organize output
timeElapsed = toc(timeStart);

Info.status = status;
Info.time = timeElapsed;
switch status
    case 0
        % fail, return the previous one
        Y_k = Y;           
    case 1
        % success, return the new iterate
        Y_k = Y_trial;
        Info.beta = beta_k;
        Info.stepSize = stepSize_trial;
        Info.merit = [merit, merit_trial];
end

end