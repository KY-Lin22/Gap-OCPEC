clear all
clc

opts = struct();
param = MyParam('test', opts);

opti = casadi.Opti();

x1 = opti.variable();

dummy = opti.parameter();
opti.set_value(dummy, 1);

disp(param);
param_eval = param(dummy);

opti.minimize(x1^2)
opti.subject_to(param_eval.monitor('testje')*x1 >= 1)

opti.set_initial(x1, 4)

opti.solver('ipopt')

opti.callback(@(i) nlp_callback(param, k))
%%
sol = opti.solve();