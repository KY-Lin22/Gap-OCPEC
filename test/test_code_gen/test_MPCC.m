clear all
clc
[~, ~, ~] = rmdir('func_gen', 's');

%%
import casadi.*
% variable
lambda = MX.sym('lambda', 1, 1);
eta = MX.sym('eta', 1, 1);
x = [lambda, eta];
% relaxation parameter
s = MX.sym('s', 1, 1);
% cost
J = (lambda - 1)^2 + (eta - 1)^2;
% symbolic expression for primal gap function
param_c = 1;
stationary_point_c = lambda - (1/param_c) * eta;
omega_c = max(0, stationary_point_c);
p_c = 0.5*lambda^2 - 0.5*omega_c^2 + lambda*(omega_c - lambda);
phi_c = eta' * (lambda - omega_c) + param_c * p_c;
% function object for primal gap function
gap_func_name = 'phi_c_func';
phi_c_func = Function(gap_func_name, {lambda, eta}, {phi_c});
gap_func_implementation = 'symbolic';
switch gap_func_implementation
    case 'symbolic'
        % direct use symbolic expression
        phi_func = phi_c_func;
    case 'codegen_fd'   
        % check funcgen folder
        if ~isfolder('./func_gen')
            mkdir('./func_gen')
            addpath(genpath('./func_gen'))
        end      
        cd './func_gen' % change current folder to next level to generate code             
        % Check if the MEX file exists in the current folder
        current_file_list = dir(); % get a list of files in the current folder
        code_filename = [gap_func_name '_' gap_func_implementation];
        is_mex_files_exists = any(strcmp({current_file_list.name}, [code_filename '.mexw64']));
        if ~is_mex_files_exists
            % code gen
            codeGenerator = CodeGenerator(code_filename, struct('mex',true));
            codeGenerator.add(phi_c_func);
            codeGenerator.generate;
            mex([code_filename '.c'], '-largeArrayDims')
        end
        % load back as external function
        phi_c_func_import = Importer([code_filename '.mexw64'],'dll');
        phi_func = external(gap_func_name, phi_c_func_import, struct('enable_fd', true));
        disp(phi_func)   
        % change current folder to main level
        cd '..'
    case 'codegen_jac'
        % NOTE: only available for 'limit_memory' hessian approximation !!!
        % custom jacobian for primal gap function (currently use automatic symbolic one, later use the explicit one)
        phi_c_jac_lambda = jacobian(phi_c, lambda);
        phi_c_jac_eta = jacobian(phi_c, eta);     
        % function object for custom jacobian of primal gap function 
        nominal_out = MX(1,1); % a empty symbolic variable for 1 x 1 phi_c and not make use of it 
        jac_phi_c_func = Function(['jac' '_' gap_func_name],...% naming convention: names as jac_fname
            {lambda, eta, nominal_out}, {phi_c_jac_lambda, phi_c_jac_eta});
        % check funcgen folder
        if ~isfolder('./func_gen')
            mkdir('./func_gen')
            addpath(genpath('./func_gen'))
        end      
        cd './func_gen' % change current folder to next level to generate code
        % Check if the MEX file exists in the current folder
        current_file_list = dir(); % get a list of files in the current folder
        code_filename = [gap_func_name '_' gap_func_implementation];
        is_mex_files_exists = any(strcmp({current_file_list.name}, [code_filename '.mexw64']));
        if ~is_mex_files_exists
            % code gen
            code_filename = [gap_func_name '_' gap_func_implementation];
            codeGenerator = CodeGenerator(code_filename, struct('mex',true));
            codeGenerator.add(phi_c_func);
            codeGenerator.add(jac_phi_c_func);
            codeGenerator.generate;
            mex([code_filename '.c'], '-largeArrayDims')
        end
        % load back as external function
        phi_c_func_import = Importer([code_filename '.mexw64'], 'dll');
        phi_func = external(gap_func_name, phi_c_func_import);
        disp(phi_func)    
        % change current folder to main level
        cd '..'        
end
% constraint
g = [lambda;...
    s - phi_func(lambda, eta)];
lbg = [0; 0];
ubg = [inf; inf];

% solver
Prob = struct('x', reshape(x, 2, 1), 'f', J, 'g', g, 'p', s);
Option = struct();
% Option.ipopt.print_level = 0;
% Option.ipopt.tol = 1e-8; % default 1e-8
% Option.ipopt.max_iter = 3000; % default 3000
Option.ipopt.hessian_approximation = 'exact'; % 'exact', 'limited-memory'
solver = casadi.nlpsol('solver', 'ipopt', Prob, Option);

%% solve
x_0 = [0; 0.5];
solution = solver('x0', x_0, 'lbg', lbg, 'ubg', ubg, 'p', 1e-8);
% extract solution
disp('optimal solution: ')
disp(full(solution.x))
disp('inequality constraint residual: ')
disp(full(solution.g))
disp('gap function residual: ')
disp(full(phi_c_func(solution.x(1), solution.x(2))));