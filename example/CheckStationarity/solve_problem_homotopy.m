function output = solve_problem_homotopy(solver, Prob, x_0, s_0)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% tol
stationary_tol = 1e-3; % default 1e-3
max_vio_tol = 1e-6; % default 1e-6
s_tol = 1e-8; % default 1e-8
s_damping = 0.1;
%
x_Init_j = x_0;
s_j = s_0;
j = 1;
while true
    % solve
    solution = solver('x0', x_Init_j, 'p', s_j,...
        'lbg', zeros(size(Prob.g, 1), 1), 'ubg', inf*ones(size(Prob.g, 1), 1));
    % extract solution
    x_Opt_j = full(solution.x);
    J_Opt_j = full(solution.f);
    max_vio_j = abs(min([x_Opt_j(1), x_Opt_j(2)]));
    comp_vio_j = abs(x_Opt_j(1) * x_Opt_j(2));
    % print
    if mod(j, 10) == 1
        disp('---------------------------------------------------------------------------------------------------------------------------------------------')
        headMsg = ' Homotopy |    s     |  lambda  |   eta    |   cost   |  max_vio  |  max_comp | iterNum |  time(ms) |';
        disp(headMsg)
    end
    prevIterMsg = ['    ',...
        num2str(j,'%10.3d'), '   | ',...
        num2str(s_j, '%10.2e'),' | ',...
        num2str(x_Opt_j(1), '%10.2e'),' | ',...
        num2str(x_Opt_j(2), '%10.2e'),' | ',...
        num2str(J_Opt_j, '%10.2e'),' | ',...
        num2str(max_vio_j, '%10.2e'),'  | ',...
        num2str(comp_vio_j, '%10.2e'),'  |   ',...
        num2str(solver.stats.iter_count, '%10.3d'),'   |  ',...
        num2str(1000*solver.stats.t_wall_total, '%10.4f'),'  | '];
    disp(prevIterMsg)    
    
    % check termination based on the current homotopy iterate
    if strcmp(solver.stats.return_status, 'Solve_Succeeded') && ( (max_vio_j <= max_vio_tol) || (s_j <= s_tol) )
        % IPOPT at the final homotopy iteration finds the optimal solution
        output.x_Opt = x_Opt_j;
        % determine stationary point
        if norm(output.x_Opt - [1;0], 2) <= stationary_tol
            output.x_Opt_type = '*';
            output.x_Opt_color = 'r';
        elseif norm(output.x_Opt - [0;1], 2) <= stationary_tol
            output.x_Opt_type = '+';
            output.x_Opt_color = 'b';
        elseif norm(output.x_Opt - [0;0], 2) <= stationary_tol
            output.x_Opt_type = 'o';
            output.x_Opt_color = 'g';
        else
            output.x_Opt_type = ' ';
            output.x_Opt_color = ' ';
        end              
        break
    elseif ~strcmp(solver.stats.return_status, 'Solve_Succeeded')
        % IPOPT at this homotopy iteration fails to find the optimal solution
        output.x_Opt = [];
        output.x_Opt_type = ' ';
        output.x_Opt_color = '';
        break
    else
        % IPOPT at this homotopy iteration (not the final) finds the optimal solution, prepare for next homotopy iteration
        x_Init_j = x_Opt_j;
        s_j = s_damping * s_j;
        j = j + 1;
    end
end

end

