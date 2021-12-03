# Example of Mixed Integer Non-Linear Programming problem. 
# Example 1 of Duran and Grossmann, an outer-approximation algorithm for  a class
# of mixed-integer nonlinear programs. Mathematical Programming 36 (1986) 307-339


var x1 >= 0 <= 2; 	# Continuous variables
var x2 >= 0 <= 2; 	# Continuous variables
var x6 >= 0 <= 1;	# Continuous variable
var y1 binary; 		# Integer variable
var y2 binary; 		# Integer variable
var y3 binary; 		# Integer variable

param U;

minimize fcObj: 5*y1 + 6*y2 + 8*y3 + 10*x1 - 7*x6 -18 * log(x2 + 1) - 19.2*log(x1 - x2 + 1) + 10;

subject to g1: 0.8*log(x2 + 1) + 0.96*log(x1 - x2 + 1) - 0.8*x6 >= 0;
subject to g2: x2 - x1 <= 0;
subject to g3: x2 - U*y1 <= 0;
subject to g4: x1 - x2 - U*y2 <= 0;
subject to g5: log(x2 + 1) + 1.2*log(x1 - x2 + 1) - x6 - U*y3 >= -2;
subject to g6: y1 + y2 <= 1;



data;

param U := 2;

option solver "muriqui";


#in this example, we solve this problem by several different algorithms. 
#Note, you will probably need solve your models just once choosing one of the
#available algorithms.


#solving model by lp-nlp BB algorithm
options muriqui_alg_choice "str MRQ_LP_NLP_BB_OA_BASED_ALG";
#options muriqui_options "str in_milp_solver MRQ_GUROBI             str in_nlp_solver MRQ_NLP_MOSEK             dbl in_max_cpu_time 600             int in_max_iterations 1000000000";
solve;


#solving model by outer approximation algorithm
options muriqui_alg_choice "str MRQ_OA_ALG";
solve;


#solving model by extended cutting plane algorithm
options muriqui_alg_choice "str MRQ_ECP_ALG";
solve;


#solving model by extended supporting hyperplane algorithm
options muriqui_alg_choice "str MRQ_ESH_ALG";
solve;


#solving model by bonmin hybrid algorithm
options muriqui_alg_choice "str MRQ_BONMIN_HYBRID_ALG";
solve;


#solving model by branch-and-bound algorithm (hybriding with outer approximation)
options muriqui_alg_choice "str MRQ_BB_ALG";
solve;


#solving model by pure branch-and-bound algorithm 
options muriqui_alg_choice "str MRQ_BB_ALG";
options muriqui_options "int in_use_outer_app 0      int in_use_outer_app_as_heuristic  0";
solve;







