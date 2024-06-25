function [Y_dot, Info] = solve_differential_equation(self, Y, p, p_dot, Y_dot_Init, epsilon)
%UNTITLED26 Summary of this function goes here
%   Detailed explanation goes here
timeStart = tic;
switch self.Option.Continuation.differential_equation_solve
    case 'FDGMRES'
        [Y_dot, rho] = self.FDGMRES_method(Y, p, p_dot, Y_dot_Init, epsilon);
    case 'direct'       
        % KKT residual, KKT matrix, and sensitivity matrix
        KKT_residual = full(self.FuncObj.KKT_residual(Y, p));
        KKT_matrix = sparse(self.FuncObj.KKT_matrix(Y, p));
        sensitivity_matrix = sparse(self.FuncObj.sensitivity_matrix(Y, p));
        % Y_dot and rho
        Y_dot = KKT_matrix\(-epsilon * KKT_residual - sensitivity_matrix * p_dot);
        rho = 0;      
    otherwise
        error('specified method is not supported')
end
% info
Info.time = toc(timeStart);
Info.GMRES_res = rho;
end