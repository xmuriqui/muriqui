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

#include <cmath>
#include "muriqui.hpp"


using namespace minlpproblem;
using namespace muriqui;

//we have to define a class to perform nonlinear function evaluations
class MyEval : public MRQ_NonLinearEval
{
    
public:
    
    
    //that method is called before other functions evaluations in a algorithm
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
    
    
    //that method is called after all functions evaluations in a algorithm
    virtual void finalize(const int nthreads, const int n, const int m, const int nzNLJac, const int nzNLLagHess) override;
    
    virtual ~MyEval(){};
};






int MyEval::initialize(const int nthreads, const int n, const int m, const int nzNLJac, const int nzNLLagHess)
{
    return 0;
}


//to evaluate objective function nonlinear part
int MyEval::eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value)
{
    value = 5*x[3] + 6*x[4] + 8*x[5] + 10*x[0] - 7*x[2] - 18*log(x[1]+1)               -19.2*log(x[0] -x[1] + 1) + 10;
    
    return 0;
}


//to evaluate nonlinear constraints part
int MyEval::eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *constrEval, const double *x, double *values)
{
    const double U = 2.0;
    const double *y = &x[3];
    
    //linear and quadratic constraints must not be evaluated
    
    if( !constrEval || constrEval[0] )
        values[0] = 0.8*log(x[1] + 1) + 0.96*log(x[0] - x[1] + 1) - 0.8*x[2]; //0.8*log(x2 + 1) + 0.96*log(x1 - x2 + 1) - 0.8*x6 >= 0;

    if( !constrEval || constrEval[4] )
        values[4] = log(x[1] + 1) + 1.2*log(x[0] - x[1] + 1) -x[2] - U*y[2]; //log(x2 + 1) + 1.2*log(x1 - x2 + 1) - x6 - U*y3 >= -2;
    
    return 0;
}


//to evaluate nonlinear objective gradient part
int MyEval::eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *grad)
{
    grad[0] = 10 - 19.2/(x[0] -x[1] + 1);
    grad[1] = -18/(x[1]+1) + 19.2/(x[0] -x[1] + 1);
    grad[2] = -7;
    grad[3] = 5;
    grad[4] = 6;
    grad[5] = 8;
    
    return 0;
}


//to evaluate nonlinear jacobian part
int MyEval::eval_grad_nl_constrs_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *constrEval, const double *x, MIP_SparseMatrix& jacobian)
{
    const double U = 2.0;

    if( !constrEval || constrEval[0] )
    {
        jacobian.setElement(0, 0, 0.96/(x[0] - x[1] + 1));
        
        jacobian.setElement(0, 1, ( 0.8/(x[1] + 1) - 0.96/(x[0] - x[1] + 1) ));
        
        jacobian.setElement(0, 2, -0.8);
    }
    
    if( !constrEval || constrEval[4] )
    {
        jacobian.setElement(4, 0, 1.2/( x[0] - x[1] + 1 ));
        
        jacobian.setElement(4, 1, ( 1/(x[1] + 1) - 1.2/(x[0] - x[1] + 1) ));
        
        jacobian.setElement(4, 2, -1);
        
        jacobian.setElement(4, 5, -U);
    }
    
    return 0;
}


//to evaluate nonlinear hessian of lagrangian (only lower triangle)
// lagrangian:  objFactor*f(x,y) + lambda*g(x,y)
int MyEval::eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double objFactor, const double *lambda, MIP_SparseMatrix& hessian)
{
    double aux;
    
    hessian.setAllSparseMatrix(0.0); //set all sparse matrix to zero
    
    
    //first we set the objective part
    
    if( objFactor != 0.0 )
    {
        aux = objFactor * 19.2/pow(x[0] -x[1] + 1, 2);
        
        hessian.setElement(0, 0, aux ); //setting value to position 0,0
        
        hessian.setElement(1, 0, -aux ); //setting value to position 1,0
        
        hessian.setElement(1, 1, objFactor * ( 18.0/pow( x[1] + 1, 2) + 19.2/pow(x[0] - x[1] + 1, 2) ) ); //setting value to position 1,1
    }
    
    
    //now, the constraints: if lambda is null, we must assume all lambdas are zero	
    if( lambda )
    {
        //cHess[0][0]
        aux = lambda[0]*( -0.96/pow( x[0] - x[1] + 1, 2) ) + lambda[4]*( -1.2/pow( x[0] - x[1] + 1, 2) );
        
        hessian.addToElement(0, 0, aux); //adding value to position 0,0
        
        
        //cHess[1][0]
        //aux = lambda[0]*( 0.96/pow( x[0] - x[1] + 1, 2) ) + lambda[4]*( 1.2/pow( x[0] - x[1] + 1, 2) );
        
        hessian.addToElement(1, 0, -aux); //adding value to position 1,0
        
        
        //cHess[1][1]
        aux = lambda[0]*( -0.8/pow( x[1] + 1, 2) - 0.96/pow( x[0] - x[1] + 1, 2) ) + lambda[4]*( -1/pow(x[1] + 1, 2) - 1.2/pow(x[0] - x[1] + 1, 2) );
        
        hessian.addToElement(1, 1, aux); //adding value to position 1,1
    }
    
    
    return 0;
}


void MyEval::finalize(const int nthreads, const int n, const int m, const int nzNLJac, const int nzNLLagHess)
{

}
