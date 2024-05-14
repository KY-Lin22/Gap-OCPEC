clear all
clc

import casadi.*
% Use the function
f = MyCallback('f', 0.5);
res = f(2);
disp(res)

x = MX.sym('x');
y = f(x);
disp(y)
