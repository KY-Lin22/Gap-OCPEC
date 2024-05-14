clear all
clc

%% custom jacobian (function: single output f, multiple input x) 
x = casadi.SX.sym('x', 2, 1);
y = casadi.SX.sym('y', 2, 1);
f = x + y;

dfdx = jacobian(f, x);
dfdy = jacobian(f, y);

f_func = casadi.Function('f', {x, y}, {f});
f_func.print_dimensions

nominal_out = casadi.SX(2,1); % for 2 x 1 f
f_jac_func = casadi.Function('jac_f', {x, y, nominal_out}, {dfdx, dfdy});
f_jac_func.print_dimensions

codeGenerator = casadi.CodeGenerator('codegen_f',struct('mex',true));
codeGenerator.add(f_func);
codeGenerator.add(f_jac_func);
codeGenerator.generate;

mex codegen_f.c -largeArrayDims

%$
f_Importer = casadi.Importer('codegen_f.mexw64', 'dll');
f_Importer.has_function('f')      % 1
f_Importer.has_function('jac_f')  % 1

f_blackbox = casadi.external('f', f_Importer);
f_blackbox.print_dimensions

x_mx = casadi.MX.sym('x',2,1);
y_mx = casadi.MX.sym('y',2,1);
f_blackbox_mx = f_blackbox(x_mx, y_mx);
f_blackbox_jac_dydx_mx = jacobian(f_blackbox_mx, x_mx);
f_blackbox_jac_dydy_mx = jacobian(f_blackbox_mx, y_mx);
