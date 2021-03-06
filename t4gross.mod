# Example of Mixed Integer Non-Linear Programming problem. 
# Example 3 of
# Duran and Grossmann, an outer-approximation algorithm for a class 
# of mixed-integer nonlinear programs. 
# Mathematical Programming 36 (1986) 307-339


var x{1..5}; 		# Continuous variables
var y{1..25} binary;	# Integer variables


param H;
param K;
param N;
param M;

param c{i in 1..25};
param w{i in 1..25, j in 1..5};
param z{i in 1..25, j in 1..5};
param d{i in 1..10, j in 1..5};
param R{i in 1..25} := min{j in 1..M}( sum{k in 1..K}( w[i,k]*(d[j,k] - z[i,k])*(d[j,k] - z[i,k]) ) ) ;

maximize fcObj:  (sum{i in 1..25}(c[i]*y[i]) -0.6*x[1]*x[1] + 0.9*x[2] + 0.5*x[3] - 0.1*x[4]*x[4] - x[5]);

subject to g1_{i in 1..N}: sum{k in 1..K}( w[i,k]*(x[k] - z[i,k])*(x[k] - z[i,k]) ) - (1 - y[i])*H <= R[i];

subject to g2: x[1] - x[2] + x[3] + x[4] + x[5] <= 10;
subject to g3: 0.6*x[1] - 0.9*x[2] - 0.5*x[3] + 0.1*x[4] + x[5] <= -0.64;
subject to g4: x[1] - x[2] + x[3] - x[4] + x[5] >= 0.69;
subject to g5: 0.157*x[1] + 0.05*x[2] <= 1.5;
subject to g6: 0.25*x[2] + 1.05*x[4] - 0.3*x[5] >= 4.5;

subject to g7:  2 <= x[1] <= 4.5;
subject to g8:  0 <= x[2] <= 8.0;
subject to g9:  3 <= x[3] <= 9.0;
subject to g10: 0 <= x[4] <= 5.0;
subject to g11: 4 <= x[5] <= 10;


data;

param H := 1000;
param K := 5;
param N := 25;
param M := 10;

param c := 
1 1 
2 0.2
3 1 
4 0.2
5 0.9
6 0.9
7 0.1
8 0.8
9 1 
10 0.4
11 1
12 0.3
13 0.1
14 0.3
15 0.5
16 0.9
17 0.8
18 0.1
19 0.9
20 1
21 1
22 1
23 0.2
24 0.7
25 0.7;

param z: 1 	2 	3 	4 	5 := 
  1	2.26 5.15 4.03 1.74 4.74
2	5.51 9.01 3.84 1.47 9.92
3	4.06 1.80 0.71 9.09 8.13
4	6.30 0.11 4.08 7.29 4.24
5	2.81 1.65 8.08 3.99 3.51
6	4.29 9.49 2.24 9.78 1.52
7	9.76 3.64 6.62 3.66 9.08
8	1.37 6.99 7.19 3.03 3.39
9	8.89 8.29 6.05 7.48 4.09
10	7.42 4.60 0.30 0.97 8.77
11	1.54 7.06 0.01 1.23 3.11
12	7.74 4.40 7.93 5.95 4.88
13	9.94 5.21 8.58 0.13 4.57
14	9.54 1.57 9.66 5.24 7.90
15	7.46 8.81 1.67 6.47 1.81
16	0.56 8.10 0.19 6.11 6.40
17	3.86 6.68 6.42 7.29 4.66
18	2.98 2.98 3.03 0.02 0.67
19	3.61 7.62 1.79 7.80 9.81
20	5.68 4.24 4.17 6.75 1.08
21	5.48 3.74 3.34 6.22 7.94
22	8.13 8.72 3.93 8.80 8.56
23	1.37 0.54 1.55 5.56 5.85
24	8.79 5.04 4.83 6.94 0.38
25	2.66 4.19 6.49 8.04 1.66 ;


param w: 1 	2 	3 	4 	5 := 
1	9.57 2.74 9.75 3.96 8.67
2	8.38 3.93 5.18 5.20 7.82
3	9.81 0.04 4.21 7.38 4.11
4	7.41 6.08 5.46 4.86 1.48
5	9.96 9.13 2.95 8.25 3.58
6	9.39 4.27 5.09 1.81 7.58
7	1.88 7.20 6.65 1.74 2.86
8	4.01 2.67 4.86 2.55 6.91
9	4.18 1.92 2.60 7.15 2.86
10	7.81 2.14 9.63 7.61 9.17
11	8.96 3.47 5.49 4.73 9.43
12	9.94 1.63 1.23 4.33 7.08
13	0.31 5.00 0.16 2.52 3.08
14	6.02 0.92 7.47 9.74 1.76
15	5.06 4.52 1.89 1.22 9.05
16	5.92 2.56 7.74 6.96 5.18
17	6.45 1.52 0.06 5.34 8.47
18	1.04 1.36 5.99 8.10 5.22
19	1.40 1.35 0.59 8.58 1.21
20	6.68 9.48 1.60 6.74 8.92
21	1.95 0.46 2.90 1.79 0.99
22	5.18 5.10 8.81 3.27 9.63
23	1.47 5.71 6.95 1.42 3.49
24	5.40 3.12 5.37 6.10 3.71
25	6.32 0.81 6.12 6.73 7.93 ;


param d: 1 2 3 4 5 :=
1	0.62 5.06 7.82 0.22 4.42
2	5.21 2.66 9.54 5.03 8.01
3	5.27 7.72 7.97 3.31 6.56
4	1.02 8.89 8.77 3.10 6.66
5	1.26 6.80 2.30 1.75 6.65
6	3.74 9.06 9.80 3.01 9.52
7	4.64 7.99 6.69 5.88 8.23
8	8.35 3.79 1.19 1.96 5.88
9	6.44 0.17 9.93 6.80 9.75
10	6.49 1.92 0.05 4.89 6.43 ;



option solver "./muriqui";
options muriqui_alg_choice "str MRQ_BB_ALG";
options muriqui_options "int in_use_outer_app 1       str in_exp_strategy MRQ_BB_ES_DEPTH     str in_nlp_solver MRQ_IPOPT        str in_milp_solver MRQ_GUROBI";

solve;

