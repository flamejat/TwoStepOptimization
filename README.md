TwoStepOptimization
===================

Energy consumption and cost minimization for integer decision variables using model predictive control with a two-step optimization algorithm.


New approach to the solution of optimization problems with discrete control variables modeled by binary integer programs (BIP). The solution of BIP is computationally demanding when the number of the BIP variables increase. Two instances that increase the number BIP variables in practical applications are the reduction of the discretization sampling time and the increase of the optimization time period. The proposed approach transforms a single BIP optimization into a linear program (LP) and $N$ feasibility BIP's, with less number of variables. The reduction of the number of variables increases the algorithm speed in providing a solution. The approach permits to solve optimization problems with longer time intervals and with a higher number of control variables, while being computationally tractable.