function FuncObj = create_FuncObj(self, nlp)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
FuncObj = struct();

%% function
% cost J
FuncObj.J = Function('J', {nlp.z, nlp.p}, {nlp.J}, {'z', 'p'}, {'J'});
% constraint h and c
FuncObj.h = Function('h', {nlp.z, nlp.p}, {nlp.h}, {'z', 'p'}, {'h'});
FuncObj.c = Function('c', {nlp.z, nlp.p}, {nlp.c}, {'z', 'p'}, {'c'});

end

