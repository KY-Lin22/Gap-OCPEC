function [Y_dot, Info] = solve_differential_equation(self, Y, p, p_dot, Y_dot_Init)
%UNTITLED26 Summary of this function goes here
%   Detailed explanation goes here
timeStart = tic;

switch self.Option.Continuation.differential_equation_solve
    case 'FDGMRES'
        % TODO
        [Y_dot, rho] = self.FDGMRES_method(Y, p, p_dot, Y_dot_Init);
    case 'direct'       
        Y_dot = full(self.FuncObj.Y_dot(Y, p, p_dot));
        rho = 0;      
    otherwise
        error('specified method is not supported')
end
% info
Info.time = toc(timeStart);
Info.GMRES_res = rho;
end