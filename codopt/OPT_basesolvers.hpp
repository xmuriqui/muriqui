

#ifndef OPT_BASESOLVERS_HPP
#define OPT_BASESOLVERS_HPP


#include <cstdlib>
#include <cstring>

#include <ostream>
#include <string>

#include <map>



#include "OPT_config.hpp"
#include "SPM_NewSparseMatrix.hpp"
#include "SPM_DynSparseMatrix.hpp"
#include "MIP_minlpProblem.hpp"



namespace optsolvers
{
    #define OPT_PRINT_CALLBACK_ERROR_MSG OPT_PRINT_NLP_CALLBACK_FUNCTION_ERROR //to print errors on nlp evaluation erros
    #define OPT_PRINT_MAX_ITER_WARNING WAXM_PRINT_MAX_ITER_WARNING //to printing a warning message about maximum number of itertions
    #define OPT_PRINT_ERROR_RETURN_CODE_ON_SOLVE WAXM_PRINT_ERROR_RETURN_CODE_ON_SOLVE //to printing message for error code after solver call in solve method
    
    
    #define OPT_INFINITY MIP_INFINITY
    
    
    #define OPT_WORHP_PARAM_FILE_NAME "optsolverworhp.xml"
    
    
    
    
    typedef newspm::SPM_NewSparseMatrix<int, double> OPT_SparseMatrix;
    //typedef newspm::SPM_SparseRow<double> OPT_SparseRow;
    typedef newspm::SPM_NewSparseMatrix<int, unsigned int> OPT_UIntSparseMatrix;
    //typedef newspm::SPM_SparseRow<unsigned int> OPT_UIntSparseRow;
    typedef newspm::SPM_DynSparseMatrix OPT_DynSparseMatrix;
    
    typedef minlpproblem::MIP_MINLPProb OPT_MINLPProb;
    typedef minlpproblem::MIP_NonLinearEval OPT_NonLinearEval;
    
    
    
    
    
    enum OPT_OPTSENSE
    {
        OPT_MINIMIZE	=	900,
        OPT_MAXIMIZE
    };
    
    
    enum OPT_LISTSOLVERS
    {
        OPT_UNDEFINEDSOLVER	= 700,
        OPT_GLPK            = 701,
        OPT_CBC             = 702,
        OPT_CPLEX           = 703,
        OPT_GUROBI          = 704,
        OPT_XPRESS          = 705,
        OPT_MOSEK           = 706,
        OPT_KNITRO          = 707,
        OPT_IPOPT           = 708,
        OPT_ALGENCAN        = 709,
        OPT_WORHP           = 710,
        OPT_OPTIZELLE       = 711,
        
        #define OPT_NUMBER_OF_SOLVERS 12 //we use that to make sure we consider all solvers in some functions like OPT_isSolverAvailable. UPDATE THAT WHEN YOU A NEW SOLVER
        
        OPT_IQUAD           = OPT_UNDEFINEDSOLVER + OPT_NUMBER_OF_SOLVERS
    };
    
    
    enum OPT_LPSOLVERS
    {
        OPT_LP_GLPK			= OPT_GLPK,
        OPT_LP_CBC			= OPT_CBC,
        OPT_LP_CPLEX		= OPT_CPLEX,
        OPT_LP_GUROBI		= OPT_GUROBI,
        OPT_LP_XPRESS		= OPT_XPRESS,
        OPT_LP_MOSEK		= OPT_MOSEK,
        OPT_LP_KNITRO		= OPT_KNITRO,
        OPT_LP_IPOPT		= OPT_IPOPT,
        OPT_LP_ALGENCAN		= OPT_ALGENCAN,
        OPT_LP_WORHP		= OPT_WORHP,
        OPT_LP_OPTIZELLE    = OPT_OPTIZELLE,
        OPT_LP_IQUAD		= OPT_IQUAD,
        
        OPT_LPS_BEGIN		= OPT_LP_GLPK,
        OPT_LPS_END			= OPT_LP_IQUAD
    };
    
    
    enum OPT_QPSOLVERS
    {
        OPT_QP_CPLEX 		= OPT_LP_CPLEX,
        OPT_QP_GUROBI		= OPT_LP_GUROBI,
        OPT_QP_XPRESS		= OPT_LP_XPRESS,
        OPT_QP_MOSEK		= OPT_LP_MOSEK,
        OPT_QP_KNITRO		= OPT_LP_KNITRO,
        OPT_QP_IPOPT		= OPT_LP_IPOPT,
        OPT_QP_ALGENCAN		= OPT_LP_ALGENCAN,
        OPT_QP_WORHP		= OPT_LP_WORHP,
        OPT_QP_OPTIZELLE    = OPT_LP_OPTIZELLE,
        OPT_QP_IQUAD		= OPT_LP_IQUAD,
        
        OPT_QPS_BEGIN		= OPT_QP_CPLEX,
        OPT_QPS_END			= OPT_QP_IQUAD
    };
    
    
    enum OPT_QCPSOLVERS
    {
        ////cplex and gurobi are not supported because has no adequate functions to set quadratic terms in constraints...
        //OPT_QCP_CPLEX		= OPT_QP_CPLEX, 
        //OPT_QCP_GUROBI	= OPT_QP_GUROBI,
        OPT_QCP_XPRESS		= OPT_QP_XPRESS,
        OPT_QCP_MOSEK		= OPT_QP_MOSEK,
        OPT_QCP_KNITRO		= OPT_QP_KNITRO,
        OPT_QCP_IPOPT		= OPT_QP_IPOPT,
        OPT_QCP_ALGENCAN	= OPT_QP_ALGENCAN,
        OPT_QCP_WORHP		= OPT_QP_WORHP,
        OPT_QCP_OPTIZELLE   = OPT_QP_OPTIZELLE,
        OPT_QCP_IQUAD		= OPT_QP_IQUAD,
        
        OPT_QCPS_BEGIN		= OPT_QCP_XPRESS,
        OPT_QCPS_END		= OPT_QCP_IQUAD
    };
    
    
    enum OPT_NLPSOLVERS
    {
        OPT_NLP_MOSEK		= OPT_QCP_MOSEK,
        OPT_NLP_KNITRO		= OPT_QCP_KNITRO,
        OPT_NLP_IPOPT		= OPT_QCP_IPOPT,
        OPT_NLP_ALGENCAN	= OPT_QCP_ALGENCAN,
        OPT_NLP_WORHP		= OPT_QCP_WORHP,
        OPT_NLP_OPTIZELLE   = OPT_QCP_OPTIZELLE,
        OPT_NLP_IQUAD		= OPT_QCP_IQUAD,
        
        OPT_NLPS_BEGIN		= OPT_NLP_MOSEK,
        OPT_NLPS_END		= OPT_NLP_IQUAD
    };
    
    
    enum OPT_SOLVERTYPE
    {
        OPT_LP 		= 600,
        OPT_QP,
        OPT_QCP,
        OPT_NLP
    };
    
    
    inline bool OPT_isQPSolverType(OPT_SOLVERTYPE solverType)
    {
        return solverType >= OPT_QP;
    }
    
    inline bool OPT_isQCPSolverType(OPT_SOLVERTYPE solverType)
    {
        return solverType >= OPT_QCP;
    }
    
    inline bool OPT_isNLPSolverType(OPT_SOLVERTYPE solverType)
    {
        return solverType == OPT_NLP;
    }
    
    
    
    enum OPT_VARTYPE
    {
        OPT_VT_CONTINUOUS 	= 201,
        OPT_VT_INTEGER
    };
    
    
    enum OPT_RETURNCODES
    {
        OPT_FEASIBLE_SOLUTION				= 1,
        OPT_OPTIMAL_SOLUTION				= 0,
        OPT_INFEASIBLE_PROBLEM				= -100,
        OPT_UNBOUNDED_PROBLEM				= -101,
        OPT_MAX_ITERATIONS					= -102,
        OPT_MAX_TIME						= -103,
        OPT_MEMORY_ERROR					= -104,
        OPT_BAD_INPUT						= -105,
        OPT_NO_FEASIBLE_SOLUTION_FOUND		= -106,
        OPT_SOLVER_ERROR					= -107,
        OPT_SUBSOLVER_ERROR					= -108,
        OPT_INCOMPATIBLE_SOLVER				= -109,
        OPT_OPERATION_NOT_SUPPORTED			= -110,
        OPT_OPERATION_NOT_IMPLEMENTED		= -111,
        OPT_LIBRARY_NOT_AVAILABLE			= -112,
        OPT_STOP_REQUIRED_BY_USER			= -113,
        OPT_CALLBACK_FUNCTION_ERROR			= -114,
        OPT_MAX_EVALUATIONS					= -115,
        OPT_LICENSE_ERROR					= -116,
        OPT_UNDEFINED_ERROR					= -117,
        OPT_UNDEFINED 						= -118
    };
    
    
    
    class OPT_GeneralSolverParams
    {
    public:
        std::map<std::string, int64_t> intParams;
        std::map<std::string, double> dblParams;
        std::map<std::string, std::string> strParams;
        
        
        int storeIntegerParameter(const char *name, const int64_t value);
        
        int storeDoubleParameter(const char *name, const double value);
        
        int storeStringParameter(const char *name, const char *value);
        
        void desallocate();
        
        void print(std::ostream &out = std::cout);
        
        void removeIntegerParameter(const char *name);
        
        void removeDoubleParameter(const char *name);
        
        void removeStringParameter(const char *name);
        
        int storeParametersFromFile(const char *fileName, const bool printErrorMsgs = true, const bool printFileOpenError = false);
        
        
        ~OPT_GeneralSolverParams();
    };
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    class OPT_Solver
    {
        
    protected:
        
        
        int naux; //size of auxIndex and auxValues
        int maux;
        int *auxIndex, *auxIndex2;
        double *auxValues, *auxValues2;
        
        unsigned int numberOfWarningsByIterLimit;
        
        
        virtual int __addConstraints(const int nm) = 0;
        
        virtual int __addVariables(const int nn, const bool initFree) = 0;
        
        virtual int allocateAuxStructures(const int size);
        
        virtual int allocateConstrStructures(const int m);
        
        virtual int allocateVarStructures(const int n);
        
        virtual int checkConstraintIndex(const int index);
        
        virtual int checkConstraintIndices(const unsigned int size, const int *indices);
        
        virtual int checkVariableIndex(const int index);
        
        virtual int checkVariableIndices(const unsigned int size, const int *indices);
        
        
        virtual void deallocateSol();
        
        
        virtual bool getMinusLambdaOnLagran(); //return attribute to solver convention of minus signal before lambda on lagrangian definition. We use that because on the dual solution, we follow the convention of linear programming (I think linear programming solvers consider this minus signal). So, we adopt in the solution of optsolver, this minus signal in the method getDualSolution
        
        
        virtual void printDblParamErrorMsg(const int error, const char *param, const double value );
        
        
        virtual void printIntParamErrorMsg(const int error, const char *param, const int value );
        
        
        virtual void printStrParamErrorMsg(const int error, const char *param, const char *value );
        
        
        virtual void resetSol();
        
        
    public:
        
        bool feasSol;
        int retCode, origSolverRetCode;
        double objValue, dualObjValue;
        double *sol, *constr;
        double *dualSolC, *dualSolV;
        
        
        unsigned int maxNumberOfWarningsByIterLimit;
        
        
        
        OPT_Solver();
        
        virtual ~OPT_Solver();
        
        virtual int addConstraints(const int nm);
        
        
        
        virtual int addVariables(const int nn, const bool initFree = true);
        
        
        virtual void initialize();
        
        //virtual int allocateMemory( const int maxPrimal, const int maxConstr = 0 );
        
        virtual void deallocateMemory();
        
        virtual void deallocateSolverEnv() = 0;
        
        virtual double getDualObjValue();
        
        virtual void getDualSolution( double *dualConstrs, double *dualVarBounds, const bool correctSignal = true );
        
        virtual int getNumberOfConstraints(int &m) = 0;
        
        virtual int getNumberOfIterations(long unsigned int &niter) = 0;
        
        virtual int getNumberOfVars(int &n) = 0;
        
        virtual double getObjValue();
        
        virtual void getSolution( double *solution, double *constraints = NULL);
        
        virtual OPT_LISTSOLVERS getSolverCode() = 0;
        
        virtual std::string getSolverName();
        
        virtual OPT_SOLVERTYPE getSolverType() = 0;
        
        virtual int getVariableType( const int index, OPT_VARTYPE &varType ) = 0;
        
        //that function can receive estimatives of maximum number of variables, constraints and nonzeros in quadratic terms matrices, but it is not mandatory...
        virtual int initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz = 0) = 0;
        
        //if this method return true, we implement the solver class derived from OPT_MyNLPSolver, i.e.,  using MIP_MINLPProb to store coefficients
        virtual bool isMyNLPClass();
        
        virtual int removeConstraints(const int ninds, const int *indices ) = 0;
        
        //all constraints in the range [begin end] will be removed
        virtual int removeConstraintsByRange(const int begin, const int end);
        
        virtual int removeVars(const int ninds, const int *indices ) = 0;
        
        virtual int setObjCutLowerBound(const double objLBound);
        
        virtual int setObjCutUpperBound(const double objUBound);
        
        virtual int setMaxCPUTime(const double time) = 0;
        
        virtual int setMaxTime(const double time);
        
        virtual int setNumberOfThreads(const int nthreads) = 0;
        
        virtual int setOutputLevel( const int level ) = 0;
        
        virtual int setRelativeDualTol( const double tol ) = 0;
    
        virtual int setRelativeOptimalityTol( const double tol ) = 0;
        
        virtual int setRelativePrimalTol( const double tol ) = 0;
        
        virtual int setDoubleParameter(const char *param, const double value) = 0;
        
        virtual int setIntegerParameter(const char *param, const int value ) = 0;
        
        virtual int setParameters( OPT_GeneralSolverParams &params );
        
        virtual int setStringParameter(const char *param, const char *value) = 0;
        
        virtual void setThreadNumber(const unsigned int threadNumber);
        
        virtual int setVariableType( const int index, const OPT_VARTYPE varType );
        
        virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true ) = 0;
        
        virtual int solveAndGetTime(double *cpuTime, double *clockTime, const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true);
        
        //just implemented by some solvers like glpk and ipopt to enforce reoptimization... Just call this method if you did not do any change in the structure of the problem (add new non zero elements, new constraints, etc...)
        virtual int warmUp();
        
        
    };
    
    
    
    
    /* {min|max} c'x + d
        s.t:
            l_c <= Ax <= u_c
            l_v <= x <= u_v
    */
    class OPT_LPSolver : public OPT_Solver
    {
        
    protected:
        
        
        
    public:
        
        OPT_LPSolver();
        
        virtual ~OPT_LPSolver();
        
        
        int addVariablesFromProb(const minlpproblem::MIP_MINLPProb& prob, const bool setObj, const bool setVarBounds, const bool setVarType, const int naddvars);
        
        //all linear coefficients are copied, inclusive the ones in quadratic objective and constraints
        virtual int copyLPPartFrom( OPT_LPSolver &other, const bool copyObj = true, const bool copyConstrs = true, const bool copyVarBounds = true, const bool copyVarTypes = true);
        
        virtual int generateModelFile(const char *fileName);
        
        virtual int getConstraintBounds( const int index, double &lb, double &ub ) = 0;
        
        virtual int getFullConstraintLinearPart(const int constrIndex, double *values);
        
        virtual int getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value) = 0;
        
        virtual int getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values) = 0;
        
        virtual int getLinearCoefsInConstraints(int &nzs, int *rows, int *cols, double *values);
        
        virtual int getNumberOfLinearCoefsInConstraints(int &nzs);
        
        virtual int getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs) = 0;
        
        
        virtual int getObjConstant(double &ObjConstant) = 0;
        
        virtual int getObjLinearCoef( const int index, double &value ) = 0;
        
        virtual int getNumberOfIntVars(int &nI) = 0;
        
        virtual int getObjSense(OPT_OPTSENSE &sense) = 0;
        
        virtual OPT_SOLVERTYPE getSolverType();
        
        virtual int getVariableBounds(const int index, double &lb, double &ub) = 0;
        
        
        //warning: this method replaces all coefficients in a constraint.
        virtual int resetConstraintLinearPart( const int constrIndex, const int nzs, const int *cols, const double *values ) = 0;
        
        //warning: this method replaces all coefficients in a constraint by the line sourceLine in sparse matrix M.
        virtual int resetConstraintLinearPart( const int constrIndex, const int sourceLine, const OPT_SparseMatrix &M )
        {
            return resetConstraintLinearPart(constrIndex, M.getNumberOfElementsAtRow(sourceLine), M[sourceLine], M(sourceLine));
        }
        
        //virtual int resetLinearConstraintPart( const int constrIndex, spm::SPM_SparseRow< double> &row );
        
        
        //warning: this method can create a new auxiliary variable because some solvers like cplex and gurobi does not acept dual bounded constraints directly...
        virtual int setConstraintBounds( const int index, const double lb, const double ub ) = 0;
        
        
        
        //virtual int setConstraintLowerBound( const int index, const double lb ) = 0;
        
        //virtual int setConstraintUpperBound( const int index, const double ub ) = 0;
        
        
        
        //set column on linear part of constraints only...
        virtual int setLinearColumn( const int varIndex, const int nzs, const int *rows, const double *values);
        
        
        virtual int setConstraintsLinearCoefs( const int nzs, const int *rows, const int *cols, const double *values );
        
        virtual int setConstraintLinearCoefs( const int constrIndex, const int nzs, const int *cols, const double *values);
        
        virtual int setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value) = 0;
        
        
        //note: this method is just to copy strict nonlinear constraints. So, this method ignores nonlinear constraints. So, the final number of constraints will be the number of strict linear constraints in pob plus naddconstrs
        int setLinearObjAndConstraintsFrom(const minlpproblem::MIP_MINLPProb &prob, const bool setObj = true, const bool setConstrs = true, const bool setVarBounds = true, const bool setVarType = true, const int naddvars = 0, const int naddconstrs = 0 );
        
        
        virtual int setObjLinearCoef( const int index, const double value ) = 0;
        
        virtual int setObjLinearCoefs( const int nzs, const int *cols, const double *values );
        
        virtual int setObjLinearPart( const int n, const double *values );
        
        
        virtual void setObjConstant(const double value) = 0;
        
        virtual int setObjSense( const OPT_OPTSENSE sense ) = 0;
        
        
        //we prefer do not declare this method like virtual because we assume if it is called from a LP pointer, user only want LP part of the problem...
        //warning: this method overwrites any method called before to set the problem...
        int setProblemFrom(const minlpproblem::MIP_MINLPProb &prob, const bool setObj = true, const bool setConstrs = true, const bool setVarBounds = true, const bool setVarType = true, const int naddvars = 0, const int naddconstrs = 0 );
        
        
        
        //set first n variable bounds
        virtual int setnVariablesBounds( const int n, const double *lb, const double *ub );
        
        
        virtual int setVariableBounds( const int index, const double lb, const double ub ) = 0;
        
        
        virtual int setVariablesBounds( const int ninds, const int *inds, const double *lb, const double *ub );
        
        
        //virtual int setVariableLowerBound( const int index, const double lb ) = 0;
        
        //virtual int setVariableUpperBound( const int index, const double ub ) = 0;
        
        
    };
    
    
    
    /* {min|max} 0.5 x'Qx + c'x + d
        s.t:
            l_c <= Ax <= u_c
            l_v <= x <= u_v  
    */
    class OPT_QPSolver : public OPT_LPSolver
    {
        
    protected:
        
    public:
        
        OPT_QPSolver();
        
        virtual ~OPT_QPSolver();
        
        //copy QP and LP part
        virtual int copyQPPartFrom( OPT_QPSolver &other, const bool copyObj = true, const bool copyConstrs = true, const bool copyVarBounds = true, const bool copyVarTypes = true);
        
        
        virtual int getNumberOfQuadObjTerms(int &nzs) = 0;
        
        virtual int getObjQuadTerm( const int row, const int col, double &value) = 0;
        
        virtual int getObjQuadPart( int &nzs, int *rows, int *cols, double *values ) = 0;
        
        virtual OPT_SOLVERTYPE getSolverType();
        
        //we prefer do not declare this method like virtual because we assume if it is called from a LP pointer, user only want LP part of the problem...
        //warning: this method overwrites any method called before...
        int setProblemFrom(const minlpproblem::MIP_MINLPProb &prob, const bool setObj = true, const bool setConstrs = true, const bool setVarBounds = true, const bool setVarType = true, const int naddvars = 0, const int naddconstrs = 0 );
        
        //just lower triangle...
        virtual int setObjQuadCoef( const int row, const int col, const double value ) = 0;
        
        //set by triple sparse format
        virtual int setObjQuadMatrix( const int nzs, const int *rows, const int *cols, const double *values );
        
        //set by compressed row format
        virtual int setObjQuadMatrix(const int *rowStart, const int *cols, const double *values);
        
        virtual int setObjQuadMatrix(const OPT_SparseMatrix &M)
        {
            int ret;
            
            if(M.getNumberOfElements() > 0)
                ret = setObjQuadMatrix( (const int *) M.offset, M[0], M(0));
            else
                ret = setObjQuadMatrix(0, NULL, NULL, NULL);
            
            return ret;
        }
        
    };
    
    
    /* {min|max} 0.5 x'Qx + c'x + d
        s.t:
            l_c_i <=  0.5 x'Q_i x + a_i x <= u_c_i
            l_x <= x <= u_x 
    */
    class OPT_QCPSolver: public OPT_QPSolver
    {
        
    protected:
        
    public:
        
        OPT_QCPSolver();
        
        virtual ~OPT_QCPSolver();
        
        //Copy QCP, QP and LP part
        virtual int copyQCPPartFrom( OPT_QCPSolver &other, const bool copyObj = true, const bool copyConstrs = true, const bool copyVarBounds = true, const bool copyVarTypes = true);
        
        
        virtual int getNumberOfConstraintQuadTerms(const int index, int &nzs) = 0; 
        
        
        //we prefer do not declare this method like virtual because we assume if it is called from a LP pointer, user only want LP part of the problem...
        //warning: this method overwrites any method called before...
        int setProblemFrom(const minlpproblem::MIP_MINLPProb &prob, const bool setObj = true, const bool setConstrs = true, const bool setVarBounds = true, const bool setVarType = true, const int naddvars = 0, const int naddconstrs = 0 );
        
        
        virtual int getConstraintQuadMatrix( const int index, int &nzs, int *rows, int *cols, double *values ) = 0;
        
        virtual OPT_SOLVERTYPE getSolverType() override;
        
        //only lower triangle... warning: that method overwrite all quadratic matrix in a constraint...
        virtual int setConstraintQuadMatrix( const int index, const int nzs, const int *rows, const int *cols, const double *values ) = 0;
        
        //set lower triangle by compressed row format
        virtual int setConstraintQuadMatrix( const int index, const int *qrowStart, const int *qcols, const double *qvalues );
        
        virtual int setConstraintQuadMatrix( const int index, const OPT_SparseMatrix &M) {
            return setConstraintQuadMatrix(index, (const int*) M.offset, M[0], M(0) );
        }
        
    };
    
    
    /*
    * {min|max} f(x) + 0.5 x'Qx + c'x + d
    * s.t:
    * 
    * 		l_c_i <= g(x) + 0.5x'Q_i x + a_i x <= u_c_i
    * 		l_x <= x <= u_x
    */
    class OPT_NLPSolver : public OPT_QCPSolver
    {
    protected:
        
        bool useSPMtoJacAndHess; //flag to use sparse matrices stored in this classes to Jacobian (J) and Hessian (H).
        
        unsigned int threadNumber; //thread number to be passed by nlEval object. By default, is 0
        
        bool *auxCEval;
        
        OPT_NonLinearEval *nlEval; //we will padronize the callback evaluations using our OPT_NonLinearEval object.
        
        OPT_SparseMatrix J; //Jacobian structure
        
        OPT_SparseMatrix lagH; //Lagrangian Hessian
        
        
        virtual int allocateConstrStructures(const int m);
        
        virtual int allocateVarStructures(const int n);
        
        //virtual int allocateSPMatrices();
        
        virtual int __removeConstraints (const int ninds, const int* indices ) = 0;
        
        
    public:
        
        double in_nl_obj_factor; //factor to nonlinear objective terms
        /*absolute and relative feasibility tolerances. Just to calculate feasSol flag in some nlp solvers.*/
        double in_absolute_feas_tol;
        double in_relative_feas_tol;
        
        
        OPT_NLPSolver();
        
        virtual ~OPT_NLPSolver();
        
        //copy NLP, QCP, QP and LP part
        virtual int copyNLPartFrom( OPT_NLPSolver &other, const bool copyObj = true, const bool copyConstrs = true, const bool copyVarBounds = true, const bool copyVarTypes = true, const bool shareEvaluationObject = true);
        
        
        virtual void deallocateMemory();
        
        virtual void initialize();
        
        virtual int getConstrNLFlag(const int index, bool &flag) = 0;
        
        
        virtual int getLagrangianHessianStructure(int &nzs, int* rows, int* cols);
        
        virtual int getJacobianStructure(int &nzs, int* rows, int* cols);
        
        virtual int getNonLinearEvalObjectPointer( OPT_NonLinearEval* &nlEval );
        
        virtual int getNumberOfNLConstraints(int& mnl) = 0;
        
        virtual int getNumberOfNonZerosInLagrangianHessian(int &nzs);
        
        virtual int getNumberOfNonZerosInJacobian(int &nzs);
        
        virtual int getObjNLFlag(bool &flag) = 0;
        
        virtual OPT_SOLVERTYPE getSolverType();
        
        virtual int initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz) override
        {
            in_nl_obj_factor = 1.0;
            return 0;
        }
        
        virtual int removeConstraints(const int ninds, const int* indices );
        
        virtual int setConstrNLFlag(const int index, const bool flag) = 0;
        
        virtual int setInitialSolution(const double* x, const double* dualConstrs, const double* dualVars) = 0;
        
        virtual int setJacobianStructure(const int nzs, const int* rows, const int* cols);
        
        //set by compressed row format
        virtual int setJacobianStructure(const int* rowStart, const int* cols);
        
        //warning: that method store the jacobian of a  nonlinear constraint. If that line have already been defined, it will be overwritten
        virtual int setJacobianRowStructure(const int row, const int nzs, const int* cols);
        
        virtual int setLagrangianHessianStructure(const int nzs, const int* rows, const int* cols);
        
        //set by compressed row format (you must set just the lower triangle)
        virtual int setLagrangianHessianStructure(const int *rowStart, const int* cols);
        
        virtual int setLagrangianHessianRowStructure( const int row, const int nzs, const int* cols);
        
        virtual int setNonLinearEvalObject( OPT_NonLinearEval *nlEval );
        
        
        virtual int setObjNLFlag(const bool flag) = 0;
        
        //we prefer do not declare this method like virtual because we assume if it is called from a LP pointer, user only want LP part of the problem...
        //warning: this method overwrites any method called before...
        virtual int setProblemFrom(const minlpproblem::MIP_MINLPProb &prob, const bool setObj = true, const bool setConstrs = true, const bool setVarBounds = true, const bool setVarType = true, const int naddvars = 0, const int naddconstrs = 0  );
        
        virtual void setThreadNumber(const unsigned int threadNumber);
    };
    
    
    
    //that class works like NLPSolver, but here we store the problem definition in our OPT_MINLPProb object. When solve method is called, we set the coefficients. Note that is necessary in solvers like Ipopt
    class OPT_MyNLPSolver : public OPT_NLPSolver
    {
    public:
        
        OPT_MINLPProb prob; //some solvers like Ipopt has no structure to store coeficients. We will use the OPT_MINLPProb to store coefficients to pass to Ipopt. To solvers that can store those coefficients like Mosek, we will not use this object.
        
    protected:
        
        
        double *xInit;
        double *zInit; //dual variables of variable bounds [z_l; z_u]
        double *lambdaInit; //dual variable of constraints...
        
        
        bool constrsBndsChg; //changes in constraint bounds
        bool nmChg; //changes in number of variables or constraints...
        bool genQuadConstrChg; //changes in some quadratic constr matrix in constraints
        bool genConstrChg; //changes in some constraint
        bool genHessChg; //changes 
        
        bool varTypeChg;
        
        bool *quadObjChg; //changes in row i in quadratic matrix of obj function
        bool *constrChg; //changes in constraint j
        //bool **quadConstrChg; //changes in constraint j, row i of quadratic matrix...
        bool *rowQuadConstrChg; //changes in row i of some quadratic constr matrix in constraints
        bool *hessChg; //changes in row i in lagrangian hessian...
        
        
        
        unsigned int *jacRowStartIndex; //general index (considering triple sparse format)
        int **jacCols; //jacCols[i] has the columns for line i in jacobian
        
        unsigned int *hessRowStartIndex;
        int **hessCols; 
        
        
        //virtual void __deallocateSolverEnv() = 0;
        
        virtual int allocateConstrStructures(const int m);
        
        virtual int allocateVarStructures(const int n);
        
        
        // __ methods from OPT_NLPSolver __
        
        
        virtual int __removeConstraints(const int ninds, const int* indices );
        
        
    public:
        
        
        OPT_MyNLPSolver();
        
        virtual ~OPT_MyNLPSolver();
        
        
        // __methods from Solver __
        
        virtual void deallocateMemory() override;
        
        virtual void deallocateSolverEnv() override;
        
        virtual void getDualSolution( double *dualConstrs, double *dualVarBounds, const bool correctSignal = true ) override;
        
        virtual int getVariableType( const int index, OPT_VARTYPE &varType ) override;
        
        virtual void initialize() override;
        
        
        //virtual void initialize();
        
        virtual int setVariableType( const int index, const OPT_VARTYPE varType ) override;
        
        
        
        // __ methods from LPSolver __
        
        
        
        
        virtual int __addConstraints(const int nm) override;
        
        virtual int __addVariables(const int nn, const bool initFree = true) override;
        
        virtual int generateModelFile(const char* fileName) override;
        
        virtual int getConstraintBounds( const int index, double &lb, double &ub ) override;
        
        virtual int getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value) override;
        
        virtual int getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values) override;
        
        virtual int getObjConstant(double &objConstant) override;
        
        virtual int getObjLinearCoef( const int index, double &value ) override;
        
        virtual int getNumberOfConstraints(int &m) override;
        
        virtual int getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs) override;
        
        virtual int getNumberOfIntVars(int &n) override;
        
        virtual int getNumberOfVars(int &n) override;
        
        virtual int getObjSense(OPT_OPTSENSE &sense) override;
        
        virtual int getVariableBounds(const int index, double &lb, double &ub) override;
        
        virtual int removeVars(const int ninds, const int *indices ) override;
        
        
        //that method is a disaster for our Sparse Matrix implementations (oriented by rows...)
        virtual int setLinearColumn( const int varIndex, const int nzs, const int* rows, const double* values) override;
        
        
        virtual int resetConstraintLinearPart( const int constrIndex, const int nzs, const int *cols, const double *values ) override;
        
        virtual int setConstraintBounds( const int index, const double lb, const double ub ) override;
        
        //virtual int setConstraintLowerBound( const int index, const double lb );
        
        //virtual int setConstraintUpperBound( const int index, const double ub );
        
        
        virtual int setConstraintsLinearCoefs( const int nzs, const int* rows, const int* cols, const double* values ) override;
        
        virtual int setConstraintLinearCoefs( const int constrIndex, const int nzs, const int *cols, const double *values) override;
        
        virtual int setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value) override;
        
        
        virtual int setObjLinearCoef( const int index, const double value ) override;
        
        virtual int setObjLinearCoefs( const int nzs, const int* cols, const double* values ) override;
        
        virtual int setObjLinearPart( const int n, const double *values ) override;
        
        
        virtual void setObjConstant(const double value) override;
        
        virtual int setObjSense( const OPT_OPTSENSE sense ) override;
        
        //set first n variable bounds
        virtual int setnVariablesBounds( const int n, const double *lb, const double *ub ) override;
        
        
        virtual int setVariableBounds( const int index, const double lb, const double ub ) override;
        
        
        virtual int setVariablesBounds( const int ninds, const int *inds, const double *lb, const double *ub ) override;
        
        
        //__ methods from QPSolver __
        
        
        virtual int getNumberOfQuadObjTerms(int &nzs) override;
        
        virtual int getObjQuadTerm( const int row, const int col, double &value) override;
        
        virtual int getObjQuadPart( int &nzs, int *rows, int *cols, double *values ) override;
        
        virtual int setObjQuadCoef( const int row, const int col, const double value ) override;
        
        virtual int setObjQuadMatrix( const int nzs, const int *rows, const int *cols, const double *values ) override;
        
        //set by compressed row format
        virtual int setObjQuadMatrix(const int *rowStart, const int *cols, const double *values) override;
        
        // __ methods from QCPSolver __
        
        virtual int getNumberOfConstraintQuadTerms( const int index, int &nzs) override; 
        
        
        virtual int getConstraintQuadMatrix( const int index, int &nzs, int *rows, int *cols, double *values ) override;
        
        virtual int setConstraintQuadMatrix( const int index, const int nzs, const int *qrows, const int *qcols, const double *qvalues )  override;
        
        virtual int setConstraintQuadMatrix( const int index, const int *qrowStart, const int *qcols, const double *qvalues ) override;
        
        // __ methods from OPT_NLPSolver __
        
        virtual int getLagrangianHessianStructure(int &nzs, int* rows, int* cols) override;
        
        virtual int getJacobianStructure(int &nzs, int* rows, int* cols) override;
        
        virtual int getConstrNLFlag(const int index, bool &flag) override;
        
        virtual int getNonLinearEvalObjectPointer( OPT_NonLinearEval* &nlEval ) override;
        
        virtual int getNumberOfNLConstraints(int& mnl) override;
        
        virtual int getNumberOfNonZerosInLagrangianHessian(int &nzs) override;
        
        virtual int getNumberOfNonZerosInJacobian(int &nzs) override;
        
        virtual int getObjNLFlag(bool &flag) override;
        
        //if this method return true, we implement the solver class derived from OPT_MyNLPSolver, i.e.,  using MIP_MINLPProb to store coefficients
        virtual bool isMyNLPClass() override;
        
        virtual int setConstrNLFlag(const int index, const bool flag) override;
        
        virtual int setInitialSolution(const double* x, const double* dualConstrs, const double* dualVars) override;
        
        virtual int setJacobianStructure(const int nzs, const int* rows, const int* cols) override;
        
        //set by compressed row format
        virtual int setJacobianStructure(const int* rowStart, const int* cols) override;
        
        //warning: that method store the jacobian of a  nonlinear constraint. If that line have already been defined, it will be overwritten
        virtual int setJacobianRowStructure(const int row, const int nzs, const int* cols) override;
        
        
        virtual int setLagrangianHessianStructure( const int nzs, const int* rows, const int* cols) override;
        
        //set by compressed row format (you must set just the lower triangle)
        virtual int setLagrangianHessianStructure(const int *rowStart, const int* cols) override;
        
        virtual int setLagrangianHessianRowStructure( const int row, const int nzs, const int* cols) override;
        
        virtual int setNonLinearEvalObject( OPT_NonLinearEval *nlEval ) override;
        
        
        virtual int setObjNLFlag(const bool flag) override;
        
        
        //we prefer do not declare this method like virtual because we assume if it is called from a LP pointer, user only want LP part of the problem...
        //int setProblemFrom( minlpproblem::MIP_MINLPProb &prob );
        
        
    protected:
        
        int allocateAuxDerivativeIndexStructures( );
        
        void deallocateAuxDerivativeIndexStructures();
        
        #if 0
            int setHessIndexRow( const int rowIndex, int &nzs );
            int setJacIndexRow( const int rowIndex, int& nzs );
        #endif
        
        
        unsigned int setHessIndexRow( const int i, const int mquad, const int *quadIndex, int* &hessColsi , unsigned int &nzl); 
        
        int setJacIndexRow( const int i, int* &jacColsi , unsigned int &nzl );
        
        int setFullJacIndex(unsigned int* jacRowStart, int** jacCols);
        
        unsigned int setFullHessIndex(const int mquad, const int* quadIndex, unsigned int* hessRowStart, int** hessCols);
    };
    
    
    
    
    
    
    
    
    
    inline std::string OPT_getSolverName(int solver)
    {
        switch(solver)
        {
            case OPT_GLPK:
                return "glpk";
                
            case OPT_CBC:
                return "cbc";
                
            case OPT_CPLEX:
                return "cplex";
                
            case OPT_GUROBI:
                return "gurobi";
                
            case OPT_XPRESS:
                return "xpress";
                
            case OPT_MOSEK:
                return "mosek";
                
            case OPT_KNITRO:
                return "knitro";
                
            case OPT_IPOPT:
                return "ipopt";
                
            case OPT_WORHP:
                return "worhp";
                
            case OPT_ALGENCAN:
                return "algencan";
            
            case OPT_OPTIZELLE:
                return "optizelle";
                
            case OPT_IQUAD:
                return "iquad";
                
            default:
                return "undefined";
        }
    }
    
    
    
    
    class OPT_ParameterSetter
    {
    public:
        
        virtual ~OPT_ParameterSetter(){}
        
        virtual int setDoubleParameter(const char *name, const double value) = 0;
        
        virtual int setIntegerParameter(const char *name, const long int value) = 0;
        
        virtual int setStringParameter(const char *name, const char *value) = 0;
    }; //we turn this class purely virtual to take advantage in another packages like Muriqui and Iquad
    
    
    
    class OPT_GeneralSolverParamsParameterSetter : public OPT_ParameterSetter
    {
    public:
        
        OPT_GeneralSolverParams *gsparams;
        
        OPT_GeneralSolverParamsParameterSetter(OPT_GeneralSolverParams *gsparams)
        {
            this->gsparams = gsparams;
        }
        
        
        virtual ~OPT_GeneralSolverParamsParameterSetter(){}
        
        
        virtual int setDoubleParameter(const char *name, const double value)
        {
            return gsparams->storeDoubleParameter(name, value);
        }
        
        virtual int setIntegerParameter(const char *name, const long int value)
        {
            return gsparams->storeIntegerParameter(name, value);
        }
        
        virtual int setStringParameter(const char *name, const char *value)
        {
            return gsparams->storeStringParameter(name, value);
        }
    };
    
    
    int OPT_readParametersWithTypeFromFile(const char *fileName, const bool printErrorMsgs, const bool printFileOpenError, OPT_ParameterSetter &psetter);
    
    
    //return true if string has only white spaces, tabulations and \n
    inline bool OPT_isEmptyString(const char *s)
    {
        for( ; *s != '\0'; s++)
        {
            if( *s != ' ' && *s != '\t' && *s != '\v' && *s != '\n' && *s != '\f' && *s != '\r' )
                return false;
        }
        
        return true;
    }
    
}


#endif
