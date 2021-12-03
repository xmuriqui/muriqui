/*
* Example of Muriqui API usage
* 
* Here, we solve the following MINLP problem:
* 
* Example of Mixed Integer Non-Linear Programming problem.
* Example 1 of Duran and Grossmann, an outer-approximation algorithm 
* for a class of mixed-integer nonlinear programs. Mathematical 
* Programming 36 (1986) 307-339
* 
* 
* var x0 >= 0 <= 2; 	# Continuous variables
* var x1 >= 0 <= 2; 	# Continuous variables
* var x2 >= 0 <= 1;	# Continuous variable
* var x3 binary; 		# Integer variable
* var x4 binary; 		# Integer variable
* var x5 binary; 		# Integer variable
* 
* param U;
* 
* minimize fcObj: 5*x3 + 6*x4 + 8*x5 + 10*x0 - 7*x2 -18 * log(x1 + 1) - 19.2*log(x0 - x1 + 1) + 10;
* 
* subject to g0: 0.8*log(x1 + 1) + 0.96*log(x0 - x1 + 1) - 0.8*x2 >= 0;
* subject to g1: x1 - x0 <= 0;
* subject to g2: x1 - U*x3 <= 0;
* subject to g3: x0 - x1 - U*x4 <= 0;
* subject to g4: log(x1 + 1) + 1.2*log(x0 - x1 + 1) - x2 - U*x5 >= -2;
* subject to g5: x3 + x4 <= 1;
* 
* data;
* 
* param U := 2;
* 
*/

#include <iostream>
#include <new>

#include "muriqui.hpp"
#include "MRQ_solvers.hpp"
#include "t1gross.hpp"


using namespace minlpproblem;
using namespace muriqui;





int main(int argc, char **argv)
{
    const int n = 6;  //number of variables
    const int m = 6;  //number of constraints
    
    const double U = 2.0;
    
    
    int ret, code = 0;
    
    MRQ_MILP_SOLVER milpSolver = MRQ_getDefaultMILPSolverCode(); //MRQ_CPLEX;
    MRQ_NLP_SOLVER nlpSolver = MRQ_getDefaultNLPSolverCode(); //MRQ_NLP_MOSEK;
    
    MyEval eval;                    //object to perform nonlinear function evaluations
    MRQ_MINLPProb prob;             //object to represent the MINLP problem
    
    MRQ_LPBBExtCutPlan lpbb;        //object to apply LP based branch-and-bound
    MRQ_LPNLPBBOuterApp	lpnlpbb;    //object to apply LP/NLP based branch-and-bound
    MRQ_OuterApp oa;                //object to apply outer approximation algorithm
    MRQ_BranchAndBound bb;          //object to apply branch-and-bound algorithm
    
    
    MRQ_GeneralSolverParams nlpsParams; //object to pass parameters to mosek solver 
    MRQ_GeneralSolverParams milpsParams; //object to pass parameters to cplex solver
    
    
    //matrix structure to represent linear constraints (only nonzeros)
    int arows[9], acols[9];
    double avals[9];
    
    //matrix structure to represent jacobian (only nonzeros)
    int jrows[7], jcols[7];
    
    //matrix structure to represent hessian of lagrangian (only nonzeros)
    int hrows[3], hcols[3];
    
    
    //first, we set the number of variables and constraints
    ret = prob.setParametersAndAllocate(n, m);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    
    /************  setting variable bounds  *************/
    //variable: x1 (index 0)
    ret = prob.setVariableLowerBound(0, 0.0);
    ret += prob.setVariableUpperBound(0, 2.0);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    //variable: x2 (index 1)
    ret =  prob.setVariableLowerBound(1, 0.0);
    ret += prob.setVariableUpperBound(1, 2.0);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    //variable: x6 (index 2)
    ret =  prob.setVariableLowerBound(2, 0.0);
    ret += prob.setVariableUpperBound(2, 1.0);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    //variable: y1 (index 3). Note this variable is binary
    ret =  prob.setVariableLowerBound(3, 0.0);
    ret += prob.setVariableUpperBound(3, 1.0);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    //variable: y2 (index 4). Note this variable is binary
    ret =  prob.setVariableLowerBound(4, 0.0);
    ret += prob.setVariableUpperBound(4, 1.0);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    //variable: y3 (index 5). Note this variable is binary
    ret =  prob.setVariableLowerBound(5, 0.0);
    ret += prob.setVariableUpperBound(5, 1.0);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    
    
    
    //setting variable types. Default type is continuous, so, we just need specify integer variables
    ret =  prob.setVariableType(3, MIP_VT_INTEGER);
    ret += prob.setVariableType(4, MIP_VT_INTEGER);
    ret += prob.setVariableType(5, MIP_VT_INTEGER);
    
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    
    //setting nonlinear evaluation object to perform evaluations
    prob.setNonLinearEvaluationObject(&eval);
    
    
    /***********  setting objective function  ***********/
    
    //MRQ_MINLPProb always minimize the problems
    
    //specifying objective has nonlinear terms. We perform the calculation in the callback object
    prob.setObjNonLinearTermFlag(true);
    
    
    
    /***********  setting constraints  ***********/
    
    //specifying constraints having nonlinear terms
    prob.setConstraintNonLinearTermFlag(0, true);
    prob.setConstraintNonLinearTermFlag(4, true);
    
    
    //setting the matrix of linear constraints. (we just need specify the nonzero elements)
    
    arows[0] = 1; acols[0] = 0; avals[0] = -1.0;   	//A[1][0] = -1.0
    arows[1] = 1; acols[1] = 1; avals[1] = 1.0;		//A[1][1] = 1.0
    arows[2] = 2; acols[2] = 1; avals[2] = 1.0;		//A[2][1] = 1.0
    arows[3] = 2; acols[3] = 3; avals[3] = -U;		//A[2][3] = -U
    arows[4] = 3; acols[4] = 0; avals[4] = 1.0;		//A[3][0] = 1.0
    arows[5] = 3; acols[5] = 1; avals[5] = -1.0;	//A[3][1] = -1.0
    arows[6] = 3; acols[6] = 4; avals[6] = -U;		//A[3][4] = -U
    arows[7] = 5; acols[7] = 3; avals[7] = 1.0;		//A[5][3] = 1.0
    arows[8] = 5; acols[8] = 4; avals[8] = 1.0;		//A[5][4] = 1.0
    
    
    ret = prob.setConstraintsLinearPart(9, arows, acols, avals);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    
    //setting constraint bounds
    ret  = prob.setConstraintLowerBound(0, 0.0);
    ret += prob.setConstraintUpperBound(1, 0.0);
    ret += prob.setConstraintUpperBound(2, 0.0);
    ret += prob.setConstraintUpperBound(3, 0.0);
    ret += prob.setConstraintLowerBound(4, -2.0);
    ret += prob.setConstraintUpperBound(5, 1.0);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    
    //we need specify the structure of jacobian matrices, i.e., positions in the constraint derivatives that could have a nonzero element.
    
    jrows[0] = 0; jcols[0] = 0;   //variable 0 has partial derivative in the contraint 0
    jrows[1] = 0; jcols[1] = 1;	//variable 1 has partial derivative in the contraint 0
    jrows[2] = 0; jcols[2] = 2;	//variable 2 has partial derivative in the contraint 0
    
    //g[2] = -( log( x[1] + 1 ) + 1.2*log( x[0] - x[1] + 1 ) - x[2] - U*y[2] + 2 );
    
    jrows[3] = 4; jcols[3] = 0;	//variable 0 has partial derivative in the contraint 4
    jrows[4] = 4; jcols[4] = 1;	//variable 1 has partial derivative in the contraint 4
    jrows[5] = 4; jcols[5] = 2;	//variable 2 has partial derivative in the contraint 4
    jrows[6] = 4; jcols[6] = 5;	//variable 5 has partial derivative in the contraint 4
    
    
    ret = prob.setJacobianStructure(7, jrows, jcols);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    
    //we need specify the structure of hessian of lagrangian, i.e, possible nonzero positions in the lagrangian hessian
    
    hrows[0] = 0; hcols[0] = 0;
    hrows[1] = 1; hcols[1] = 0;
    hrows[2] = 1; hcols[2] = 1;
    
    ret = prob.setLagrangianHessianStructure(3, hrows, hcols);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    
    //we can pritn the problem to check. Nonlinear terms will only apper as f(x) ou g(x)
    prob.print();
    
    
    //checking first order derivatives. To be sure our callbacks are calculating derivatives correctly. You just need that if you wanna check the derivatives.
    {
        bool answer = false;
        
        double x[n];
        
        //checking in the point {0.5}^n
        for(int i = 0; i < n; i++)
            x[i] = 0.5;
        
        std::cout << "Checking first order derivatives\n";
        
        prob.checkFisrtDerivatives(true, true, x, 1e-4, NULL, 1e-6, answer);
        
        if(answer)
            std::cout << "First derivatives are ok\n";
        else
            std::cout << "Inconsitency in the first derivative calculation\n";
        
        
        prob.checkSecondDerivatives(true, true, x, 1.0, NULL, 1e-4, NULL, 1e-6, answer );
        
        if(answer)
            std::cout << "Second derivatives are ok\n";
        else
            std::cout << "Inconsitency in the second derivative calculation\n";
    }
    
    
    //now, we can solve our problem using one of minlp algorithms
    
    //we can pass parameters to nlp and milp solvers
    
    //nlpsParams.storeIntegerParameter("MSK_IPAR_INTPNT_MAX_ITERATIONS", 1000); //increase number of iterations in Mosek to 1000
    //milpsParams.storeDoubleParameter("CPX_PARAM_TILIM", 40*60*60); //enforce a time limit of 4 hours to cplex
    
    
    
    //first, we use LP based branch-and-bound
    
    lpbb.in_milp_solver = milpSolver;
    lpbb.in_nlp_solver= nlpSolver; //if parameter lpbb.in_refine_final_solution_using_nlp is set to false, lpbb does not nedd a nlp solver
    
    
    
    std::cout << "Solving problem by LP based branch-and-bound\n";
    
    lpbb.run(prob, &milpsParams, &nlpsParams);
    
      std::cout << "Return code: " << lpbb.out_return_code << " status of optimization: \"" << MRQ_getStatus(lpbb.out_return_code) << "\" cpu time: " << lpbb.out_cpu_time << " clock time: " << lpbb.out_clock_time << " number of iterations: " << lpbb.out_number_of_iterations << "\n";
    
    
    //now, we use LP/NLP based branch-and-bound
    
    //setting some parameters to algorithms
    lpnlpbb.in_milp_solver = milpSolver;
    lpnlpbb.in_nlp_solver = nlpSolver;
    
    lpnlpbb.in_max_time = 600; //600 seconds
    
    std::cout << "Solving problem by LP/NLP based branch-and-bound\n";
    
    lpnlpbb.run(prob, &milpsParams, &nlpsParams);
    
    std::cout << "Return code: " << lpnlpbb.out_return_code << " status of optimization: \"" << MRQ_getStatus(lpnlpbb.out_return_code) << "\" cpu time: " << lpnlpbb.out_cpu_time << " clock time: " << lpnlpbb.out_clock_time << " number of iterations: " << lpnlpbb.out_number_of_iterations << "\n";
    
    if( lpnlpbb.out_feasible_solution ) 
    {
        std::cout << "Objective: " << lpnlpbb.out_best_obj << "\n";
        std::cout << "Solution:\n";
        for(int i = 0; i < n; i++)
            std::cout << "x[" << i << "]: " << lpnlpbb.out_best_sol[i] << "\n";
    }
    
    
    
    std::cout << "Solving problem by outer approximation\n";
    
    oa.in_milp_solver = milpSolver;
    oa.in_nlp_solver = nlpSolver;
    
    oa.run(prob, &milpsParams, &nlpsParams);
    std::cout << "Return code: " << oa.out_return_code << " status of optimization: \"" << MRQ_getStatus(oa.out_return_code) << "\" cpu time: " << oa.out_cpu_time << " clock time: " << oa.out_clock_time <<  " number of iterations: " << oa.out_number_of_iterations << "\n";
    
    if( oa.out_feasible_solution ) 
    {
        std::cout << "Objective: " << oa.out_best_obj << "\n";
        std::cout << "Solution:\n";
        for(int i = 0; i < n; i++)
            std::cout << "x[" << i << "]: " << oa.out_best_sol[i] << "\n";
    }
    
    
    
    std::cout << "Solving problem by branch-and-bound\n";
    
    bb.in_milp_solver = milpSolver;
    bb.in_nlp_solver = nlpSolver;
    
    
    if(nlpSolver == MRQ_NLP_MOSEK) //mosek cannot be used with igma
        bb.in_igma2_strategy = MRQ_BB_IHS_NO_HEURISTICS;
    
    bb.run(prob, &milpsParams, &nlpsParams);
    std::cout << "Return code: " << bb.out_return_code << " status of optimization: \"" << MRQ_getStatus(bb.out_return_code) << "\" cpu time: " << bb.out_cpu_time << " clock time: " << bb.out_clock_time << "\n";
    
    if( bb.out_feasible_solution ) 
    {
        std::cout << "Objective: " << bb.out_best_obj << "\n";
        std::cout << "Solution:\n";
        for(int i = 0; i < n; i++)
            std::cout << "x[" << i << "]: " << bb.out_best_sol[i] << "\n";
    }
    
    
    
    std::cout << "Solving problem by branch-and-bound without oa, i.e., by pure branch-and-bound\n";
    
    bb.in_use_outer_app = false;
    bb.in_use_outer_app_as_heuristic = false;
    
    bb.run(prob, &milpsParams, &nlpsParams);
    std::cout << "Return code: " << bb.out_return_code << " status of optimziation: \"" << MRQ_getStatus(bb.out_return_code) << "\" cpu time: " << bb.out_cpu_time << " clock time: " << bb.out_clock_time << "\n";
    
    if( bb.out_feasible_solution ) 
    {
        std::cout << "Objective: " << bb.out_best_obj << "\n";
        std::cout << "Solution:\n";
        for(int i = 0; i < n; i++)
            std::cout << "x[" << i << "]: " << bb.out_best_sol[i] << "\n";
    }
    
    
    

    
termination:
    
    return code;
}

