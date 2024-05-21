function phi_c_func = create_phi_c_func(self, OCPEC, param_c)
% return phi_c_func which is a function object with output phi_c and input lambda, eta

import casadi.*

%% create function object phi_c_func: init
% create strongly convex function d and its derivative
[d_func, d_grad, d_hessian] = self.create_strongly_convex_func(OCPEC);
% init variable (function object used to generate C code can be SX functions, see: //groups.google.com/g/casadi-users/c/bFTIg2tYgVs/m/KMcCRLMYCQAJ)
lambda = SX.sym('lambda', OCPEC.Dim.lambda, 1);
eta = SX.sym('eta', OCPEC.Dim.lambda, 1); 
% init omega_solver
omega_solver = self.create_omega_solver(OCPEC, param_c);
% compute omega_c expression (pure symbolic or cvx/qp solver)
omega_c = omega_solver(lambda, eta); % - later also wrap omega_solver as a external function with output omega_c and input lambda, eta?
% compute phi_c expression (generally a mixture of symbolic and cvx/qp solver)
p_c = d_func(lambda) - d_func(omega_c) + d_grad(lambda) * (omega_c - lambda);
phi_c = eta' * (lambda - omega_c) + param_c * p_c;
% create function object for phi_c expression 
phi_c_func_init = Function(self.gap_func_name, {lambda, eta}, {phi_c}, {'lambda', 'eta'}, {'phi_c'});

%% create function object phi_c_func: final
if strcmp(OCPEC.VISetType, 'finitely_representable') || strcmp(OCPEC.VISetType, 'polyhedral')
    %% phi_c_func is a mixture of symbolic and cvx/qp solver
    switch self.gap_func_implementation
        case 'symbolic'
            error(['gap_func_implementation (symbolic) does not support VISetType: finitely_representable or polyhedral, ', ...
                'because in these cases, omega_solver is a function object that defines a convex or QP optimization solver, ', ...
                'which is non-symbolic and can not be directly embedded in CasADi symbolic graph to compute derivative'])
        case 'codegen_fd'
            % TODO (may be achieved by a dedicated wrapper using FD or the solver can provide sensitivity, i.e., variable omega w.r.t parameter lambda, eta)
            % casadi blog: https://web.casadi.org/blog/nlp_sens/
        case 'codegen_jac'
            % TODO (may be achieved by a dedicated wrapper using FD or the solver can provide sensitivity, i.e., variable omega w.r.t parameter lambda, eta)
            % casadi blog: https://web.casadi.org/blog/nlp_sens/
    end

elseif strcmp(OCPEC.VISetType, 'box_constraint') || strcmp(OCPEC.VISetType, 'nonnegative_orthant')
    %% phi_c_func is a highly nonlinear symbolic expression
    switch self.gap_func_implementation
        case 'symbolic'
            % Note: directly return phi_c_func_init as output function object phi_c_func
            phi_c_func = phi_c_func_init;
        case 'codegen_fd'
            % Note: regard phi_c_func as a black box function of lambda, eta, the jacobian is evaluated by CasADi option: finite differences

            % check func_gen folder
            if ~isfolder('./func_gen')
                mkdir('./func_gen')
                addpath(genpath('./func_gen'))
            end
            cd './func_gen' % change current folder to func_gen folder to generate code
            % Check if the MEX file exists in the func_gen folder
            current_file_list = dir(); % get a list of files in the current folder
            param_c_string = [num2str(fix(param_c)) '_' num2str(mod(floor(param_c*10^1), 10)) num2str(mod(floor(param_c*10^2), 10))];
            code_filename = [self.gap_func_name '_' self.gap_func_implementation '_' 'c' param_c_string '_' self.strongly_convex_func '_' 'Dim' num2str(OCPEC.Dim.lambda)];
            is_mex_files_exists = any(strcmp({current_file_list.name}, [code_filename '.mexw64']));
            if ~is_mex_files_exists
                % generate code file
                codeGenerator = CodeGenerator(code_filename, struct('mex',true));
                codeGenerator.add(phi_c_func_init);
                codeGenerator.generate;
                disp('generating C code for gap function successfully')
                mex([code_filename '.c'], '-largeArrayDims')
            end        
            % load mex code back as a CasADi external function
            phi_c_func_import = Importer([code_filename '.mexw64'], 'dll');
            phi_c_func = external(self.gap_func_name, phi_c_func_import, struct('enable_fd', true));
            % change current folder back to main level
            cd '..'
        case 'codegen_jac'
            % Note: regard phi_c_func as a black box function of lambda, eta, the jacobian is evaluated by CasADi option: custom jacobian

            % compute custom jacobian of phi_c expression
            % phi_c_jac_lambda = jacobian(phi_c, lambda);
            % phi_c_jac_eta = jacobian(phi_c, eta);
            phi_c_jac_lambda = eta' - param_c * (lambda - omega_c)' * d_hessian(lambda);  
            phi_c_jac_eta = (lambda - omega_c)';
            % create function object for custom jacobian of phi_c expression
            nominal_out = SX(1,1); % a empty symbolic variable for 1 x 1 phi_c and not make use of it
            jac_phi_c_func_init = Function(['jac' '_' self.gap_func_name],...% naming convention: names as jac_fname
                {lambda, eta, nominal_out}, {phi_c_jac_lambda, phi_c_jac_eta});
            % check func_gen folder
            if ~isfolder('./func_gen')
                mkdir('./func_gen')
                addpath(genpath('./func_gen'))
            end
            cd './func_gen' % change current folder to func_gen folder to generate code 
            % Check if the MEX file exists in the func_gen folder
            current_file_list = dir(); % get a list of files in the current folder
            param_c_string = [num2str(fix(param_c)) '_' num2str(mod(floor(param_c*10^1), 10)) num2str(mod(floor(param_c*10^2), 10))];
            code_filename = [self.gap_func_name '_' self.gap_func_implementation '_' 'c' param_c_string '_' self.strongly_convex_func '_' 'Dim' num2str(OCPEC.Dim.lambda)];
            is_mex_files_exists = any(strcmp({current_file_list.name}, [code_filename '.mexw64']));
            if ~is_mex_files_exists
                % generate code file
                codeGenerator = CodeGenerator(code_filename, struct('mex',true));
                codeGenerator.add(phi_c_func_init);
                codeGenerator.add(jac_phi_c_func_init);
                codeGenerator.generate;
                disp('generating C code for gap function successfully')
                mex([code_filename '.c'], '-largeArrayDims')
            end
            % load mex code back as a CasADi external function
            phi_c_func_import = Importer([code_filename '.mexw64'], 'dll');
            phi_c_func = external(self.gap_func_name, phi_c_func_import);
            % change current folder back to main level
            cd '..'
    end

end

end