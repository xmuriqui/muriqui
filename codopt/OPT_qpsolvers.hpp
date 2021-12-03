

#ifndef OPT_QPSOLVERS_HPP
#define OPT_QPSOLVERS_HPP



#if OPT_HAVE_CPLEX
extern "C" {
    #include "cplex.h"
}
#endif

#if OPT_HAVE_GUROBI
extern "C" {
    #include "gurobi_c.h"
}
#endif

#if OPT_HAVE_XPRESS
extern "C" {
    #include "xprs.h"
}
#endif




#include "OPT_basesolvers.hpp"


namespace optsolvers{
    
    //we just set set cplex as QP solver becuase interface to QCP is too poor and it is not so compatible with optsolver. To avoid loss eficiency in other solvers changing our API, we prefer just set CPLEX as a QP solver
    
    class OPT_Cplex : public OPT_QPSolver
    {
    protected:
        
        double objConstant;
        
        
        // __methods from Solver __
        
        
        virtual int __addConstraints(const int nm)  override;
        
        virtual int __addVariables(const int nn, const bool initFree = true)  override;
        
        
        
    public:
        
        
        #if OPT_HAVE_CPLEX
            CPXENVptr     env;
            CPXLPptr      prob;
        #endif
        
        
        OPT_Cplex();
        
        virtual ~OPT_Cplex();
        
        // __methods from Solver __
        
        
        virtual void deallocateSolverEnv()  override;
        
        virtual int getNumberOfConstraints(int &m) override;
        
        virtual int getNumberOfVars(int &n) override;
        
        virtual int getNumberOfIterations(long unsigned int& niter) override;
        
        virtual OPT_LISTSOLVERS getSolverCode() override;
        
        virtual int getVariableType( const int index, OPT_VARTYPE &varType ) override;
        
        virtual void initialize() override;
        
        virtual int initSolverEnv(const int maxConstrs = 0, const int maxVars = 0, const int maxQuadNz = 0) override;
        
        virtual int removeConstraints(const int ninds, const int* indices ) override;
        
        //all constraints in the range [begin end] will be removed
        virtual int removeConstraintsByRange(const int begin, const int end) override;
        
        virtual int removeVars(const int ninds, const int *indices ) override;
        
        virtual int setObjCutLowerBound(const double objLBound) override;
        
        virtual int setObjCutUpperBound(const double objUBound) override;
        
        virtual int setMaxCPUTime(const double time) override;
        
        virtual int setMaxTime(const double time) override;
        
        virtual int setNumberOfThreads(const int nthreads) override;
        
        virtual int setOutputLevel( const int level ) override;
        
        virtual int setRelativeDualTol( const double tol ) override;
    
        virtual int setRelativeOptimalityTol( const double tol ) override;
        
        virtual int setRelativePrimalTol( const double tol ) override;
        
        virtual int setDoubleParameter(const char *param, const double value) override;
        
        virtual int setIntegerParameter(const char *param, const int value ) override;
        
        virtual int setStringParameter(const char *param, const char *value) override;
        
        virtual int setVariableType( const int index, const OPT_VARTYPE varType ) override;
        
        virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true) override;
        
        
        
        // __ methods from LPSolver __
        
        
        
        
        virtual int generateModelFile(const char* fileName) override;
        
        virtual int getConstraintBounds( const int index, double &lb, double &ub ) override;
        
        virtual int getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value) override;
        
        virtual int getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values) override;
        
        virtual int getLinearCoefsInConstraints(int &nzs, int *rows, int *cols, double *values) override;
        
        virtual int getObjConstant(double &objConstant) override;
        
        virtual int getObjLinearCoef( const int index, double &value ) override;
        
        virtual int getNumberOfLinearCoefsInConstraints(int &nzs) override;
        
        virtual int getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs) override;
        
        virtual int getNumberOfIntVars(int &n) override;
        
        virtual int getObjSense(OPT_OPTSENSE &sense) override;
        
        virtual int getVariableBounds(const int index, double &lb, double &ub) override;
        
        virtual int setLinearColumn( const int varIndex, const int nzs, const int* rows, const double* values) override;
        
        
        virtual int resetConstraintLinearPart( const int constrIndex, const int nzs, const int* cols, const double* values ) override;
        
        virtual int setConstraintBounds( const int index, const double lb, const double ub ) override;
        
        
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
        
        //virtual int setVariableLowerBound( const int index, const double lb );
        
        //virtual int setVariableUpperBound( const int index, const double ub );
        
        
        
        //__ methods from QPSolver __
        
        
        virtual int getNumberOfQuadObjTerms(int &nzs) override;
        
        virtual int getObjQuadTerm( const int row, const int col, double &value) override;
        
        virtual int getObjQuadPart( int &nzs, int *rows, int *cols, double *values ) override;
        
        virtual int setObjQuadCoef( const int row, const int col, const double value ) override;
        
        virtual int setObjQuadMatrix( const int nzs, const int *rows, const int *cols, const double *values ) override;
        
        //set by compressed row format
        virtual int setObjQuadMatrix(const int *rowStart, const int *cols, const double *values) override;
        
        // __ methods from QCPSolver __
        
        
        //int setQuadConstraint( const int index, const int nzs, const int *qrows, const int *qcols, const double *qvalues );
        
        
        // own methods
        virtual int cloneFrom( OPT_LPSolver &other);
        
    
    protected:
        
        //internal methods...
        
        char constraintSense( const double lb, const double ub );
        
        int deleteQuadObjMatrix();
        
        int getSolution(const bool integer, const bool storeSol, const bool storeConstrs, const bool storeDualSol);
        
        
        friend int OPT_setQCPProblemOnCplex(const minlpproblem::MIP_MINLPProb& prob, OPT_Cplex *optCplex, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int naddvars, const int naddconstrs);
        
        //friend int OPT_setQCPProblemOnCplex( minlpproblem::MIP_MINLPProb& prob, OPT_Cplex *optCplex, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType );
        
    };
    
    
    
    
    
    class OPT_Gurobi : public OPT_QPSolver
    {
    protected:
        
        int *indAuxVarConstr; //index of auxiliary variable used to set contsraints... (I hate Gurobi)
        
        
        // __ methods from LPSolver __
        
        
        
        virtual int __addConstraints(const int nm) override;
        
        virtual int __addVariables(const int nn, const bool initFree = true) override;
        
        virtual int getMyNumberOfVars(int &n);
        
    public:
        
        #if OPT_HAVE_GUROBI
            GRBenv   *env;
            GRBmodel *prob;
        #endif
        
        
        
        OPT_Gurobi();
        
        virtual ~OPT_Gurobi();
        
        // __methods from Solver __
        
        virtual void deallocateMemory() override;
        
        virtual void deallocateSolverEnv() override;
        
        virtual int getNumberOfIterations(long unsigned int &niter) override;
        
        virtual OPT_LISTSOLVERS getSolverCode() override;
        
        virtual int getVariableType( const int index, OPT_VARTYPE &varType ) override;
        
        virtual void initialize() override;
        
        virtual int initSolverEnv(const int maxConstrs = 0, const int maxVars = 0, const int maxQuadNz = 0) override;
        
        virtual int removeConstraints(const int ninds, const int *indices ) override;
        
        virtual int removeVars(const int ninds, const int *indices ) override;
        
        virtual int setObjCutLowerBound(const double objLBound) override;
        
        virtual int setObjCutUpperBound(const double objUBound) override;
        
        virtual int setMaxCPUTime(const double time) override;
        
        virtual int setNumberOfThreads(const int nthreads) override;
        
        virtual int setOutputLevel( const int level ) override;
        
        virtual int setRelativeDualTol( const double tol ) override;
    
        virtual int setRelativeOptimalityTol( const double tol ) override;
        
        virtual int setRelativePrimalTol( const double tol ) override;
        
        virtual int setDoubleParameter(const char *param, const double value) override;
        
        virtual int setIntegerParameter(const char *param, const int value ) override;
        
        virtual int setStringParameter(const char *param, const char *value) override;
        
        virtual int setVariableType( const int index, const OPT_VARTYPE varType ) override;
        
        virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true) override;
        
        
        
        // __ methods from LPSolver __
        
        
        virtual int generateModelFile(const char* fileName) override;
        
        virtual int getConstraintBounds( const int index, double &lb, double &ub ) override;
        
        virtual int getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value) override;
        
        virtual int getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values) override;
        
        virtual int getNumberOfConstraints(int &m) override;
        
        virtual int getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs) override;
        
        virtual int getNumberOfIntVars(int &nI) override;
        
        virtual int getNumberOfVars(int &n) override;
        
        virtual int getObjConstant(double &objConstant) override;
        
        virtual int getObjLinearCoef( const int index, double &value ) override;
        
        virtual int getObjSense(OPT_OPTSENSE &sense) override;
        
        virtual int getVariableBounds(const int index, double &lb, double &ub) override;
        
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
        
        virtual int setObjQuadMatrix(const int *rowStart, const int *cols, const double *values) override;
        
        // __ methods from QCPSolver __
        
        
        //virtual int setQuadConstraint( const int index, const int nzs, const int *qrows, const int *qcols, const double *qvalues );
        
        
        
        
    protected:
        
        //internal methods...
        
        
        
    };
    
    
    
    
    
    
    class OPT_Xpress : public OPT_QCPSolver
    {
    protected:
        
        
        // __ methods from LPSolver __
        
        virtual int __addConstraints(const int nm);
        
        virtual int __addVariables(const int nn, const bool initFree = true);
        
        
        
    public:
        
        #if OPT_HAVE_XPRESS
            XPRSprob prob;
        #endif
        
        
        OPT_Xpress();
        
        virtual ~OPT_Xpress();
        
        // __methods from Solver __
        
        virtual void deallocateSolverEnv() override;
        
        virtual int getNumberOfIterations(long unsigned int &niter) override;
        
        virtual OPT_LISTSOLVERS getSolverCode() override;
        
        virtual int getVariableType( const int index, OPT_VARTYPE &varType ) override;
        
        virtual void initialize() override;
        
        virtual int initSolverEnv(const int maxConstrs = 0, const int maxVars = 0, const int maxQuadNz = 0) override;
        
        virtual int removeVars(const int ninds, const int *indices) override;
        
        virtual int removeConstraints(const int ninds, const int *indices) override;
        
        virtual int setObjCutLowerBound(const double objLBound) override;
        
        virtual int setObjCutUpperBound(const double objUBound) override;
        
        virtual int setMaxCPUTime(const double time) override;
        
        virtual int setNumberOfThreads(const int nthreads) override;
        
        virtual int setOutputLevel( const int level ) override;
        
        virtual int setRelativeDualTol(const double tol) override;
    
        virtual int setRelativeOptimalityTol( const double tol ) override;
        
        virtual int setRelativePrimalTol(const double tol) override;
        
        virtual int setDoubleParameter(const char *param, const double value) override;
        
        virtual int setIntegerParameter(const char *param, const int value) override;
        
        virtual int setStringParameter(const char *param, const char *value) override;
        
        virtual int setVariableType( const int index, const OPT_VARTYPE varType ) override;
        
        virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true) override;
        
        
        
        // __ methods from LPSolver __
        
        
        virtual int generateModelFile(const char* fileName) override;
        
        virtual int getConstraintBounds( const int index, double &lb, double &ub) override;
        
        virtual int getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value) override;
        
        virtual int getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values) override;
        
        virtual int getNumberOfConstraints(int &m) override;
        
        virtual int getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs) override;
        
        virtual int getNumberOfIntVars(int &n) override;
        
        virtual int getNumberOfVars(int &n) override;
        
        virtual int getObjConstant(double &objConstant) override;
        
        virtual int getObjLinearCoef( const int index, double &value) override;
        
        virtual int getObjSense(OPT_OPTSENSE &sense) override;
        
        virtual int getVariableBounds(const int index, double &lb, double &ub) override;
        
        virtual int setLinearColumn( const int varIndex, const int nzs, const int* rows, const double* values) override;
        
        
        virtual int resetConstraintLinearPart(const int constrIndex, const int nzs, const int *cols, const double *values) override;
        
        virtual int setConstraintBounds(const int index, const double lb, const double ub) override;
        
        
        virtual int setConstraintsLinearCoefs( const int nzs, const int* rows, const int* cols, const double* values ) override;
        
        virtual int setConstraintLinearCoefs( const int constrIndex, const int nzs, const int *cols, const double *values) override;
        
        virtual int setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value) override;
        
        
        virtual int setObjLinearCoef( const int index, const double value) override;
        
        virtual int setObjLinearCoefs( const int nzs, const int* cols, const double* values ) override;
        
        virtual int setObjLinearPart( const int n, const double *values ) override;
        
        
        virtual void setObjConstant(const double value) override;
        
        virtual int setObjSense( const OPT_OPTSENSE sense ) override;
        
        
        virtual int setVariableBounds( const int index, const double lb, const double ub ) override;
        
        
        
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
        
        virtual int setConstraintQuadMatrix( const int index, const int nzs, const int *qrows, const int *qcols, const double *qvalues ) override;
        
        //virtual int setQuadConstraintMatrix( const int index, const int *qrowStart, const int *qcols, const double *qvalues ) override;
        
        
    protected:
        
        //internal methods...
        
        char boundFlag(const double lb, const double ub);
        
    };
    
    
    
    
    
}



#endif
