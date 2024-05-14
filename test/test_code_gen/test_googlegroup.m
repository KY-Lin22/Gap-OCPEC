%% finite difference (function: single output f, single input x)
% reference://groups.google.com/g/casadi-users/c/bFTIg2tYgVs/m/KMcCRLMYCQAJ
clear all
clc

x = casadi.SX.sym('x',2);
f = x;
for i=1:10
  f = sin(f)*transpose(cos(f));
end
df = jacobian(f,x);

fun = casadi.Function('fun',{x},{f});
jacfun = casadi.Function('jac_fun',{x},{df});

fun.print_dimensions

codeGenerator = casadi.CodeGenerator('codegen_fun',struct('mex',true));
codeGenerator.add(fun);
codeGenerator.add(jacfun);
codeGenerator.generate;

mex codegen_fun.c -largeArrayDims
funImporter = casadi.Importer('codegen_fun.mexw64', 'dll');

funImporter.has_function('fun')      % 1
funImporter.has_function('jac_fun')  % 1

importedFun = casadi.external('fun', funImporter, struct('enable_fd', true));
importedFun.print_dimensions % Input 0 ("i0"): 2x1, Output 0 ("o0"): 2x2

importedJac = casadi.external('jac_fun', funImporter, struct('enable_fd', true));
importedJac.print_dimensions % Input 0 ("i0"): 2x1, Output 0 ("o0"): 4x2

x_mx = casadi.MX.sym('x',2,1);
mxFunMEX = importedFun(x_mx);

jacobian(mxFunMEX, x_mx)

%% custom jacobian (function: single output f, single input x)
% reference://groups.google.com/g/casadi-users/c/bFTIg2tYgVs/m/KMcCRLMYCQAJ
clear all
clc

x = casadi.SX.sym('x',2); 
% methods(x)
f = x;
for i=1:10
  f = sin(f)*transpose(cos(f));
end
% methods(f)
df = jacobian(f,x); % 2x2 f will be first reshaped as 4 x 1 f in a column manner
                    % and then compute the jacobian 4 X 2 df, i.e., 4 x 1 f w.r.t. 2 x 1 x

fun = casadi.Function('fun',{x},{f});
% methods(fun)
nominal_out = casadi.SX(2,2); % a empty symbolic variable for 2x2 f, required by CasADi Jacobian syntax !!
                              % see also: https://groups.google.com/g/casadi-users/c/IfDEURmTXNQ/m/QIa9Gq6hBQAJ
                              %           https://github.com/casadi/casadi/issues/2490
                              %           https://groups.google.com/g/casadi-users/c/fjStOu3YZB0
jacfun = casadi.Function('jac_fun',{x,nominal_out},{df}); % naming convention: names as jac_fname so that CasADi regard it as jacobian of function fname
                                                          % section 6.3.3 in https://web.casadi.org/docs/#importing-a-function-with-external                                                          
fun.print_dimensions
jacfun.print_dimensions

codeGenerator = casadi.CodeGenerator('codegen_fun',struct('mex',true));
codeGenerator.add(fun);
codeGenerator.add(jacfun);
codeGenerator.generate;
% methods(codeGenerator)

mex codegen_fun.c -largeArrayDims
funImporter = casadi.Importer('codegen_fun.mexw64', 'dll');
% methods(funImporter)
funImporter.has_function('fun')      % 1
funImporter.has_function('jac_fun')  % 1

importedFun = casadi.external('fun', funImporter);
importedFun.print_dimensions % Input 0 ("i0"): 2x1, Output 0 ("o0"): 2x2
% methods(importedFun)

% importedJac = casadi.external('jac_fun', funImporter);
% importedJac.print_dimensions % Input 0 ("i0"): 2x1, Output 0 ("o0"): 4x2  

x_Mx = casadi.MX.sym('x',2,1); % current only available MX, if want to use SX, need to define eval_sx in external (unknown detail) 
mxFunMEX = importedFun(x_Mx);

jacmxFunMex = jacobian(mxFunMEX, x_Mx);

