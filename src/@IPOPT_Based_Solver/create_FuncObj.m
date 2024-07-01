function FuncObj = create_FuncObj(self) 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
import casadi.*
FuncObj = struct();

%% IPOPT solver for solving the parameterized NLP
NLP_Prob = struct('x', self.NLP.z, 'f', self.NLP.J, 'g', [self.NLP.h; self.NLP.c], 'p', self.NLP.s); 
IPOPT_Option = self.Option.IPOPT_Solver;
FuncObj.IPOPT_Solver = nlpsol('Solver', 'ipopt', NLP_Prob, IPOPT_Option);

end