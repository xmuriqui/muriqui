/*
* Librarie to handle different solvers to handle optimization problems. We start from linear and quadratic problems, and we will handle famous solvers like cplex, gurobi, mosek, etc... :D
* 
* 
* 
* 
* Author: Wendel Melo
* 
* Date: 22-Jan-2015
* 
*/


#ifndef OPT_SOLVERS_HPP
#define OPT_SOLVERS_HPP


#include "OPT_basesolvers.hpp"
#include "OPT_lpsolvers.hpp"
#include "OPT_qpsolvers.hpp"
#include "OPT_nlpsolvers.hpp"





namespace optsolvers
{
    
    
    OPT_LPSolver* OPT_newLPSolver( const int solver );
    
    
    OPT_QPSolver* OPT_newQPSolver( const int solver );
    
    
    OPT_QCPSolver* OPT_newQCPSolver( const int solver );
    
    
    OPT_NLPSolver* OPT_newNLPSolver( const int solver );
    
    
    //we assume solver is a appropriate object to receive the minlp problem. For example, if prob has nonlinear functions, solver must be an OPT_NLPSolver object. This function DOES call initSolverEnv.
    int OPT_setMINLPProblem( const minlpproblem::MIP_MINLPProb& prob, optsolvers::OPT_Solver* solver, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int naddvars = 0, const int naddconstrs = 0  );
    
    
    int OPT_setQCPProblemOnCplex(const minlpproblem::MIP_MINLPProb& prob, OPT_Cplex *optCplex, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int naddvars = 0, const int naddconstrs = 0  );
    
    
    bool OPT_isSolverAvailable(int solverCode);
    
    
    int OPT_copyConstraintLinearParts( int beginOriginIndex, int endOriginIndex, int beginDestinIndex, OPT_LPSolver &origin, OPT_LPSolver &destin );
    
    
    class OPT_NLPNonObjEval : public optsolvers::OPT_NonLinearEval
    {
        bool mynewx;
        int norig, morig, nzjacorig, nzhessorig;
        minlpproblem::MIP_NonLinearEval *oeval;
        
    public:
        
        OPT_NLPNonObjEval( const int noriginal, const int moriginal, const int nzjacoriginal, const int nzhessoriginal, minlpproblem::MIP_NonLinearEval *originalEval );
        
        virtual ~OPT_NLPNonObjEval();
        
        void initialize(const int noriginal, const int moriginal, const int nzjacoriginal, const int nzhessoriginal, minlpproblem::MIP_NonLinearEval *originalEval);
        
        
        virtual int eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value) override;
        
        virtual int eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *constrEval, const double *x, double *values) override;
        
        virtual int eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *values) override;
        
        virtual int eval_grad_nl_constrs_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *constrEval, const double *x, minlpproblem::MIP_SparseMatrix& jacobian) override;
        
        virtual int eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double objFactor, const double *lambda, minlpproblem::MIP_SparseMatrix& hessian) override;
    };
    
    
    //calculator of box (lower and upper bounds) to constraints
    class OPT_ConstrsBoundsCalculator
    {
    public:
        
        /*
        * calcConstr[i] flag if method should calculate bounds to constraint i. If it is NULL, we calculate for all constraints;
        * 
        * solver is an object where will set the problem. So, You just need allocate memory for it before call this function.
        * 
        */
        static int calculate( minlpproblem::MIP_MINLPProb &prob, const double *olc, const double *ouc, OPT_LPSolver *solver, OPT_GeneralSolverParams *params, const bool calcEvenOnOriginalValues, const bool *calcConstr, double *lc, double *uc );
        
    private:
        static int setProblemBase( minlpproblem::MIP_MINLPProb &prob, const double *olc, const double *ouc, OPT_LPSolver *solver, OPT_GeneralSolverParams *params, const bool calcEvenOnOriginalValues, const bool *calcConstr);
    };
    
    
    
    class OPT_ObjCutSetter
    {
    public:
        
        int setObjCut( OPT_LPSolver *solver, const int constrIndex, const minlpproblem::MIP_MINLPProb &prob, const double zu, minlpproblem::MIP_NonLinearEval *eval ) const;
        
        int updateObjCut( OPT_LPSolver *solver, const int constrIndex, const minlpproblem::MIP_MINLPProb &prob, const double zu ) const;
    };
    
    
    //class to handle MINLPProb or solver
    class OPT_MINLPProbOrSolverUnifier
    {
        
    public:
        
        OPT_MINLPProb *prob;
        OPT_LPSolver *solver;
        
        
        OPT_MINLPProbOrSolverUnifier();
        
        void setAsMINLPProblem(OPT_MINLPProb *prob);
        
        void setAsSolver(OPT_LPSolver *solver);
        
        int setVariableBounds(int varIndex, double lb, double ub);
        
        int setConstraintBounds(int constrIndex, double lb, double ub);
        
        int setConstraintQuadMatrix( int constrIndex, int nvars, const int *rows, const int *cols, const double *values);
        
        int setConstraintLinearCoefs( int constrIndex, int nvars, const int *cols, const double *values);
        
    };
    
    
    
    /*
    * class to add a contsraint in the form:
    * 
    * sum_{i \in I} |x_i - w_i| <= b,      where I is a set of indices, w is an ordinary solution (constant) and b is a constant. x is the variable vector in the problem
    * 
    * Actually, we use the square of two-norm to implement this constraint:
    * 
    * sum_{i \in I} (w_i - x_i)² <= b²,
    * 
    * 
    */ 
    class OPT_SolutionDistanceConstraintSetter
    {
        
    public:
        
        static int setDistConstraint( OPT_MINLPProb &prob, const int constrIndex, const int nvars, const int *varIndices, const double *sol, const double maxDistance, int *auxInds = NULL, double *auxVars = NULL);
        
        static int setDistConstraint( OPT_QCPSolver &solver, const int constrIndex, const int nvars, const int *varIndices, const double *sol, const double maxDistance, int *auxInds = NULL, double *auxVars = NULL);
        
        static int setDistConstraint( OPT_MINLPProbOrSolverUnifier &unifier, const int constrIndex, const int nvars, const int *varIndices, const double *sol, const double maxDistance, int *auxInds = NULL, double *auxVars = NULL);
    };
    
    
    /*
    * class to decrease variable bounds around an input solution. Let I be a set of variable indices, and w be the input soluton. Let \alpha be a parameter in [0, 1]
    * 
    * for each variable index i in I, we update lower bound (l_i) and (u_i) as :
    * 
    *    l_i = max(  w_i - \alpha * d_i  ,  l_i )
    *    u_i = min(  w_i + \alpha * d_i  ,  u_i )
    * 
    * 	d_i can be defined as the box size, i. e.:
    * 
    *  d_i = u_i - l_i   if u_i and l_i are not infinity
    * 
    *  if u_i or l_i are infinity, we can set d_i as:
    *  
    *  d_i = |w_i| + 1.0   //we sum 1.0 because w_i can be zero. In this case, l_i and u_i would be equal
    * 
    */ 
    class OPT_SolutionDistanceBoxConstraintsSetter
    {
        
    public:
        
        static int setDistConstraint( OPT_MINLPProb &prob, const double *olx, const double *oux, const int nvars, const int *varIndices, const double *sol, const double percentageToNeighborhood );
        
        static int setDistConstraint( OPT_LPSolver &solver, const double *olx, const double *oux, const int nvars, const int *varIndices, const double *sol, const double percentageToNeighborhood );
        
        static int setDistConstraint( OPT_MINLPProbOrSolverUnifier &unifier, const double *olx, const double *oux, const int nvars, const int *varIndices, const double *sol, const double percentageToNeighborhood );
        
    };
    
    
};





























#endif
