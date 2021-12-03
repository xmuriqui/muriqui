# Example of Mixed Integer Non-Linear Programming problem. 
# Example 2 of
# Duran and Grossmann, an outer-approximation algorithm for a class 
# of mixed-integer nonlinear programs. 
# Mathematical Programming 36 (1986) 307-339


var x3 >= 0 <= 2; 	# Continuous variables
var x5 >= 0 <= 2; 	# Continuous variables
var x9 >= 0 <= 2;	# Continuous variable
var x11 >= 0;		# Continuous variable
var x13 >= 0;		# Continuous variable
var x16 >= 0 <= 3;	# Continuous variable

var y1 binary; 		# Integer variable
var y2 binary; 		# Integer variable
var y3 binary; 		# Integer variable
var y4 binary;		# Integer variable
var y5 binary;		# Integer variable

param U;

minimize fcObj: 5*y1 + 8*y2 + 6*y3 + 10*y4 + 6*y5 - 10*x3 - 15*x5 - 15*x9 + 15*x11 + 5*x13 - 20*x16 + exp(x3) + exp(x5/1.2) - 60*log(x11 + x13 + 1) + 140;

subject to g1: -log(x11 + x13 + 1) <= 0 ;
subject to g2: -x3 - x5 - 2*x9 + x11 + 2*x16 <= 0;
subject to g3: -x3 - x5 - 0.75*x9 + x11 + 2*x16 <= 0;
subject to g4: x9 - x16 <= 0;
subject to g5: 2*x9 - x11 - 2*x16 <= 0 ;
subject to g6: -0.5*x11 + x13 <= 0;
subject to g7: 0.2*x11 - x13 <= 0;
subject to g8: exp(x3) - U*y1 <= 1;
subject to g9: exp(x5/1.2) - U*y2 <= 1;
subject to g10: 1.25*x9 - U*y3 <= 0;
subject to g11: x11 + x13 - U*y4 <= 0;
subject to g12: -2*x9 + 2*x16 - U*y5 <= 0;
subject to g13: y1 + y2 = 1;
subject to g14: y4 + y5 <= 1;


data;

param U := 10;


option solver "muriqui";

#seting algorithm to muriqui
option muriqui_alg_choice "str MRQ_OA_ALG";

#seting muriqui algorithm options, setting CPLEX as MILP solver and IPOPT as in_nlp_solver
option muriqui_options "str in_milp_solver MRQ_CPLEX         str in_nlp_solver MRQ_IPOPT ";


#seting options to milp solver (CPLEX)
option muriqui_milp_options "int CPX_PARAM_THREADS 1        dbl CPX_PARAM_TILIM 100";

#seting options to nilp solver (IPOPT)
#option muriqui_nlp_options "int max_iter 5000         str linear_solver  ma27";






solve;

