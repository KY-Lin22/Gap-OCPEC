function [Y, Info] = solve_first_NLP(self, z_Init, p)
%UNTITLED25 Summary of this function goes here
%   Detailed explanation goes here
disp('---------------------------------------------------------------------------------------------------')
switch self.Option.Continuation.first_NLP_solve
    case 'non_interior_point'
        [z, Info_NIP] = self.non_interior_point_method(z_Init, p);
        gamma_h = Info_NIP.gamma_h;
        gamma_c = Info_NIP.gamma_c;
        terminal_status = Info_NIP.terminal_status;
        terminal_msg = Info_NIP.terminal_msg;
        time = Info_NIP.Time.total;
    case 'IPOPT'
        s = p(1);
        solution = self.FuncObj.IPOPT_Solver('x0', z_Init, 'p', s,...
            'lbg', [zeros(self.NLP.Dim.h, 1); zeros(self.NLP.Dim.c, 1)],...
            'ubg', [zeros(self.NLP.Dim.h, 1); inf*ones(self.NLP.Dim.c, 1)]);
        z = full(solution.x);
        gamma_h = full(solution.lam_g(1 : self.NLP.Dim.h, 1));
        gamma_c = full(solution.lam_g(self.NLP.Dim.h + 1 : end, 1));
        solver_status = (strcmp(self.FuncObj.IPOPT_Solver.stats.return_status, 'Solve_Succeeded')) ...
        || (strcmp(self.FuncObj.IPOPT_Solver.stats.return_status, 'Solved_To_Acceptable_Level'))...
        || (strcmp(self.FuncObj.IPOPT_Solver.stats.return_status, 'Feasible_Point_Found'));
        terminal_status = solver_status;
        terminal_msg = self.FuncObj.IPOPT_Solver.stats.return_status;
        time = self.FuncObj.IPOPT_Solver.stats.t_wall_total; % self.Solver.t_proc_total; 
    otherwise
        error('specified method is not supported')
end
% assemble Y
Y = [z; gamma_h; gamma_c];
% info
Info.terminal_status = terminal_status;
Info.terminal_msg = terminal_msg;
Info.time = time;
end