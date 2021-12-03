/*
* 
* Example of usage Muriqui API
* 
* Here, we solve the following MINLP problem:
* 
* 
* 
* 
* Example of Mixed Integer Non-Linear Programming problem. 
* Example 2 of Duran and Grossmann, an outer-approximation algorithm for a
* class of mixed-integer nonlinear programs. Mathematical Programming 36 (1986) 307-339
* 
* 
* var x0 >= 0 <= 2; 	# Continuous variables
* var x1 >= 0 <= 2; 	# Continuous variables
* var x2 >= 0 <= 2;	# Continuous variable
* var x3 >= 0;			# Continuous variable
* var x4 >= 0;			# Continuous variable
* var x5 >= 0 <= 3;	# Continuous variable
* 
* var x6 binary; 		# Integer variable
* var x7 binary; 		# Integer variable
* var x8 binary; 		# Integer variable
* var x9 binary;		# Integer variable
* var x10 binary;		# Integer variable
* 
* param U;
* 
* minimize fcObj: 5*x6 + 8*x7 + 6*x8 + 10*x9 + 6*x10 - 10*x0 - 15*x1 - 15*x2 + 15*x3 + 5*x4 - 20*x5 + exp(x0) + exp(x1/1.2) - 60*log(x3 + x4 + 1) + 140;
* 
* subject to g0: -log(x3 + x4 + 1) <= 0 ;
* subject to g1: -x0 - x1 - 2*x2 + x3 + 2*x5 <= 0;
* subject to g2: -x0 - x1 - 0.75*x2 + x3 + 2*x5 <= 0;
* subject to g3: x2 - x5 <= 0;
* subject to g4: 2*x2 - x3 - 2*x5 <= 0 ;
* subject to g5: -0.5*x3 + x4 <= 0;
* subject to g6: 0.2*x3 - x4 <= 0;
* subject to g7: exp(x0) - U*x6 <= 1;
* subject to g8: exp(x1/1.2) - U*x7 <= 1;
* subject to g9: 1.25*x2 - U*x8 <= 0;
* subject to g10: x3 + x4 - U*x9 <= 0;
* subject to g11: -2*x2 + 2*x5 - U*x10 <= 0;
* subject to g12: x6 + x7 = 1;
* subject to g13: x9 + x10 <= 1;
* 
* 
* data;
* 
* param U := 10;
* 
*/


#include <cmath>

#include <iostream>
#include <new>

#include "muriqui.hpp"
#include "MRQ_solvers.hpp"


using namespace muriqui;
using namespace minlpproblem;



//we have to define a class to perform nonlinear function evaluations

class MyEval : public MRQ_NonLinearEval
{
    
public:
    
    
    //that method is called before other functions evaluations
    virtual int initialize(const int nthreads, const int n, const int m, const int nzNLJac, const int nzNLLagHess) override;
    
    
    //to evaluate objective function nonlinear part
    virtual int eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value) override;
    
    
    //to evaluate nonlinear constraints part
    virtual int eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *constrEval, const double *x, double *values) override;
    
    
    //to evaluate nonlinear objective gradient part
    virtual int eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double* x, double* grad) override;
    
    
    //to evaluate nonlinear jacobian part
    virtual int eval_grad_nl_constrs_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *constrEval, const double *x, MIP_SparseMatrix& jacobian) override;
    
    
    //to evaluate nonlinear hessian of lagrangian (only lower triangle)
    // lagrangian:  objFactor*f(x,y) + lambda*g(x,y)
    virtual int eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double objFactor, const double *lambda, MIP_SparseMatrix& hessian) override;
    
    
    
    virtual ~MyEval(){};
};








int MyEval::initialize(const int nthreads, const int n, const int m, const int nzNLJac, const int nzNLLagHess)
{
    return 0;
}


//to evaluate objective function nonlinear part
int MyEval::eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value)
{
    //note, linear part was alredy set in the problem definition. Here, we must evaluate the strict nonlinear part
    value = exp(x[0]) + exp(x[1]/1.2) - 60*log(x[3] + x[4] + 1);
    
    return 0;
}


//to evaluate nonlinear constraints part
int MyEval::eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *constrEval, const double *x, double *values)
{
    //note, linear part was alredy set in the problem definition. Here, we must evaluate the strict nonlinear part

    if( !constrEval || constrEval[0] )
        values[0] = -log( x[3] + x[4] + 1 );
    
    if( !constrEval || constrEval[7] )
        values[7] = exp(x[0]);
    
    if( !constrEval || constrEval[8] )
        values[8] = exp(x[1]/1.2);
    
    return 0;
}


//to evaluate nonlinear objective gradient part
int MyEval::eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double* x, double* grad)
{
    //note, linear part was alredy set in the problem definition. Here, we must evaluate the strict nonlinear part
    
    grad[0] = exp(x[0]);
    grad[1] = exp(x[1]/1.2)/1.2;
    grad[2] = 0;
    grad[3] = -60/(x[3] + x[4] + 1);
    grad[4] = -60/(x[3] + x[4] + 1);
    grad[5] = 0;
    grad[6] = 0;
    grad[7] = 0;
    grad[8] = 0;
    grad[9] = 0;
    grad[10] = 0;
    
    return 0;
}


//to evaluate nonlinear jacobian part
int MyEval::eval_grad_nl_constrs_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *constrEval, const double *x, MIP_SparseMatrix& jacobian)
{
    //note, linear part was alredy set in the problem definition. Here, we must evaluate the strict nonlinear part
    
    
    if( !constrEval || constrEval[0] )
    {
        jacobian.setElement(0, 3, -1.0/(x[3] + x[4] + 1));
        jacobian.setElement(0, 4, -1.0/(x[3] + x[4] + 1) );
    }
    
    if( !constrEval || constrEval[7] )
    {
        jacobian.setElement(7, 0, exp(x[0]));
    }
    
    if( !constrEval || constrEval[8] )
    {
        jacobian.setElement(8, 1, exp( x[1]/1.2 )/1.2);
    }
    
    return 0;
}


//to evaluate nonlinear hessian of lagrangian (only lower triangle)
// lagrangian:  objFactor*f(x,y) + lambda*g(x,y)
int MyEval::eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double objFactor, const double *lambda, MIP_SparseMatrix& hessian)
{
    hessian.setAllSparseMatrix(0.0);
    
    if(objFactor != 0.0)
    {
        hessian.setElement(0, 0, objFactor* exp(x[0]));
        hessian.setElement(1, 1, objFactor* exp(x[1]/1.2)/(1.2 * 1.2) );
        hessian.setElement(3, 3, objFactor* 60/pow( x[3] + x[4] + 1, 2));
        hessian.setElement(4, 3, objFactor* 60/pow( x[3] + x[4] + 1 , 2));
        hessian.setElement(4, 4, objFactor* 60/pow( x[3] + x[4] + 1 , 2));
    }
    
    if( lambda )
    {
        hessian.addToElement(3, 3, lambda[0] * 1/pow(x[3] + x[4] + 1, 2) );
        hessian.addToElement(4, 3, lambda[0] * 1/pow(x[3] + x[4] + 1, 2) );
        hessian.addToElement(4, 4, lambda[0] * 1/pow(x[3] + x[4] + 1, 2) );
        hessian.addToElement(0, 0, lambda[7] * exp(x[0]) );
        hessian.addToElement(1, 1, lambda[8] * exp(x[1]/1.2)/(1.2*1.2) );
    }
    
    
    return 0;
}









int main(int argc, char **argv)
{
    const int nx = 6;
    const int ny = 5;
    const int n = nx + ny;  //number of variables
    const int m = 14;  //number of constraints
    
    const double U = 10;
    
    
    int ret, code = 0;
    
    MRQ_MILP_SOLVER milpSolver = MRQ_getDefaultMILPSolverCode(); //MRQ_CPLEX;
    MRQ_NLP_SOLVER nlpSolver = MRQ_getDefaultNLPSolverCode(); //MRQ_NLP_MOSEK;
    
    MyEval eval;					//object to perform nonlinear function evaluations
    MRQ_MINLPProb prob;   			//object to represent the MINLP problem
    
    
    MRQ_ContinuousRelax crelax;		//object to solve continuous relaxation
    
    MRQ_OuterApp oa;  				//object to apply outer approximation algorithm
    MRQ_BranchAndBound bb; 			//object to apply branch-and-bound algorithm
    
    
    MRQ_GeneralSolverParams nlpsParams; //object to pass parameters to mosek solver 
    MRQ_GeneralSolverParams milpsParams; //object to pass parameters to cplex solver
    
    
    double lx[n] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //variable lower bounds
    
    double ux[n] = {2.0, 2.0, 2.0, MRQ_INFINITY, MRQ_INFINITY, 3, 1, 1, 1, 1, 1}; //variable upper bounds
    
    
    //first, we set the number of variables and constraints
    ret = prob.setParametersAndAllocate(n, m);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    
    /************  setting variable bounds  *************/
    ret = prob.setVariableLowerBounds(n, lx);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    ret = prob.setVariableUpperBounds(n, ux);
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    
    //setting variable types. Default type is continuous, so, we just need specify integer variables
    for(int i = nx; i < n; i++)
    {
        ret = prob.setVariableType(i, MIP_VT_INTEGER);
        if( ret != 0 )
        {
            code = -1;
            goto termination;
        }
    }
    
    //setting nonlinear evaluation object to perform evaluations
    prob.setNonLinearEvaluationObject(&eval);
    
    
    /***********  setting objective function  ***********/
    
    //MRQ_MINLPProb always minimize the problems
    
    //specifying objective has nonlinear terms. We perform the calculation in the callback object
    prob.setObjNonLinearTermFlag(true);
    
    //here, we set the linear part of objetive function separately and let to callback only the strict nonlinear part. Off course, we could put all objective function in the callback evaluation also.
    
    {
        const double coefs[n] = {-10, -15, -15, 15, 5, -20,         5, 8, 6, 10, 6};   //linear coeficients in the objective function
        const double d = 140.0;   //constant term in the objective function
        
        //setting the objective constant term
        prob.setObjConstant(d);
        
        
        ret = prob.setObjLinearCoefficients(coefs);
        if( ret != 0 )
        {
            code = -1;
            goto termination;
        }
        
    }
    
    
    
    /***********  setting constraints  ***********/
    
    //specifying constraints having nonlinear terms
    ret = prob.setConstraintNonLinearTermFlag(0, true);
    ret += prob.setConstraintNonLinearTermFlag(7, true);
    ret += prob.setConstraintNonLinearTermFlag(8, true);
    
    if( ret != 0 )
    {
        code = -1;
        goto termination;
    }
    
    //setting the matrix of linear constraints. (we just need specify the nonzero elements). //here we specify the matrix by compresses row format
    {
        int rowStart[m+1] = {0, 0, 5, 10, 12, 15, 17, 19, 20, 21, 23, 26, 29, 31, 33};
        int acols[33];
        double avals[33];
        
        
        
        //g1: -x0 - x1 - 2*x2 + x3 + 2*x5 
        acols[0] = 0; avals[0] = -1;
        acols[1] = 1; avals[1] = -1;
        acols[2] = 2; avals[2] = -2;
        acols[3] = 3; avals[3] = 1;
        acols[4] = 5; avals[4] = 2;
        
        //g2: -x0 - x1 - 0.75*x2 + x3 + 2*x5
        acols[5] = 0; avals[5] = -1;
        acols[6] = 1; avals[6] = -1;
        acols[7] = 2; avals[7] = -0.75;
        acols[8] = 3; avals[8] = 1;
        acols[9] = 5; avals[9] = 2;
        
        //g3: x2 - x5
        acols[10] = 2; avals[10] =  1;
        acols[11] = 5; avals[11] = -1;
        
        //g4: 2*x2 - x3 - 2*x5
        acols[12] = 2; avals[12] =  2;
        acols[13] = 3; avals[13] = -1;
        acols[14] = 5; avals[14] = -2;
        
        //g5: -0.5*x3 + x4
        acols[15] = 3; avals[15] = -0.5;
        acols[16] = 4; avals[16] = 1;
        
        //g6: 0.2*x3 - x4
        acols[17] = 3; avals[17] = 0.2;
        acols[18] = 4; avals[18] = -1;
        
        //g7: exp(x0) - U*x6
        acols[19] = 6; avals[19] = -U;
        
        //g8: exp(x1/1.2) - U*x7
        acols[20] = 7; avals[20] = -U;
        
        //g9: 1.25*x2 - U*x8
        acols[21] = 2; avals[21] = 1.25;
        acols[22] = 8; avals[22] = -U;
        
        //g10: x3 + x4 - U*x9
        acols[23] = 3; avals[23] = 1;
        acols[24] = 4; avals[24] = 1;
        acols[25] = 9; avals[25] = -U;
        
        //g11: -2*x2 + 2*x5 - U*x10
        acols[26] = 2; avals[26] = -2;
        acols[27] = 5; avals[27] = 2;
        acols[28] = 10; avals[28] = -U;
        
        //g12: x6 + x7
        acols[29] = 6; avals[29] = 1;
        acols[30] = 7; avals[30] = 1;
        
        //g13: x9 + x10
        acols[31] = 9; avals[31] = 1;
        acols[32] = 10; avals[32] = 1;
        
        
        ret = prob.setConstraintsLinearPart(rowStart, acols, avals);
        
        if( ret != 0 )
        {
            code = -1;
            goto termination;
        }
    }
    
    
    //setting contsraint bounds
    {
        const double uc[m] = {0, 0, 0, 0, 0,    0, 0, 1, 1, 0,    0, 0, 1, 1};
        
        ret = prob.setConstraintUpperBounds(uc);
        if( ret != 0 )
        {
            code = -1;
            goto termination;
        }
        
        //constraint 12 is an equality constraint. So, we have to set the lower bound also
        ret = prob.setConstraintLowerBound(12, 1.0);
        if( ret != 0 )
        {
            code = -1;
            goto termination;
        }
    }
    
    
    //setting structure of jacobian. Note, linear terms in nonlinear constarints was aleardy set above;
    {
        int jrows[4], jcols[4];
        
        
        //g0: -log(x3 + x4 + 1)
        jrows[0] = 0; jcols[0] = 3;
        jrows[1] = 0; jcols[1] = 4;
        
        //g7: exp(x0) - U*x6 <= 1;
        jrows[2] = 7; jcols[2] = 0;
        
        //g8: exp(x1/1.2) - U*x7
        jrows[3] = 8; jcols[3] = 1;
        
        ret = prob.setJacobianStructure(4, jrows, jcols);
        if( ret != 0 )
        {
            code = -1;
            goto termination;
        }
    }
    
    
    //setting structure of hessian of lagrangian (only lower triangle)
    {
        int hrows[5], hcols[5];
        
        hrows[0] = 3; hcols[0] = 3;
        hrows[1] = 4; hcols[1] = 3;
        hrows[2] = 4; hcols[2] = 4;
        hrows[3] = 0; hcols[3] = 0;
        hrows[4] = 1; hcols[4] = 1;
        
        ret = prob.setLagrangianHessianStructure(5, hrows, hcols);
        if( ret != 0 )
        {
            code = -1;
            goto termination;
        }
        
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
    
    
    
    std::cout << "Solving continuous relaxation\n";
    
    crelax.in_nlp_solver = nlpSolver;
    
    crelax.run(prob, NULL, &nlpsParams);
    std::cout << "Return code: " << crelax.out_return_code << " status of optimization: \"" << MRQ_getStatus(crelax.out_return_code) << "\" cpu time: " << crelax.out_cpu_time << " clock time: " << crelax.out_clock_time << "\n";
    
    if( crelax.out_feasible_solution ) 
    {
        std::cout << "Objective: " << crelax.out_best_obj << "\n";
        std::cout << "Solution:\n";
        for(int i = 0; i < n; i++)
            std::cout << "x[" << i << "]: " << crelax.out_best_sol[i] << "\n";
    }
    
    
    
    std::cout << "Solving problem by outer approximation\n";
    
    oa.in_milp_solver = milpSolver;
    oa.in_nlp_solver = nlpSolver;
    
    oa.run(prob, &milpsParams, &nlpsParams);
    std::cout << "Return code: " << oa.out_return_code << " status of optimization: \"" << MRQ_getStatus(oa.out_return_code) << "\" cpu time: " << oa.out_cpu_time << " clock time: " << oa.out_clock_time << "\n";
    
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
    
    
    
    
termination:
    
    return code;
}


















