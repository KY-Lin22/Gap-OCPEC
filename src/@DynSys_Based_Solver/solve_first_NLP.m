function [Y, Info] = solve_first_NLP(self, z_Init, p)
%UNTITLED25 Summary of this function goes here
%   Detailed explanation goes here
disp('---------------------------------------------------------------------------------------------------')
% load relaxation parameter
s = p(1);
% solve
solution = self.FuncObj.IPOPT_Solver('x0', z_Init, 'p', s,...
    'lbg', [zeros(self.NLP.Dim.h, 1); zeros(self.NLP.Dim.c, 1)],...
    'ubg', [zeros(self.NLP.Dim.h, 1); inf*ones(self.NLP.Dim.c, 1)]);
% extract information
z = full(solution.x);
gamma_h = full(solution.lam_g(1 : self.NLP.Dim.h, 1));
gamma_c = full(solution.lam_g(self.NLP.Dim.h + 1 : end, 1));
solver_status = (strcmp(self.FuncObj.IPOPT_Solver.stats.return_status, 'Solve_Succeeded')) ...
    || (strcmp(self.FuncObj.IPOPT_Solver.stats.return_status, 'Solved_To_Acceptable_Level'))...
    || (strcmp(self.FuncObj.IPOPT_Solver.stats.return_status, 'Feasible_Point_Found'));
terminal_status = solver_status;
terminal_msg = self.FuncObj.IPOPT_Solver.stats.return_status;
time = self.FuncObj.IPOPT_Solver.stats.t_wall_total; % t_proc_total;
% assemble Y
Y = [z; gamma_h; gamma_c];
% info
Info.terminal_status = terminal_status;
Info.terminal_msg = terminal_msg;
Info.time = time;
end