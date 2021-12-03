

#ifndef OPT_LPSOLVERS_HPP
#define OPT_LPSOLVERS_HPP



#if OPT_HAVE_GLPK
extern "C"{
    #include "glpk.h"
}
#endif


#if OPT_HAVE_CBC_OR_OSI
    #include "OsiSolverInterface.hpp"
#endif

#if OPT_HAVE_CBC
    #include "CbcModel.hpp"
    #include "CbcCutGenerator.hpp"
#endif




#include "OPT_basesolvers.hpp"


namespace optsolvers {
    
    class OPT_Glpk : public OPT_LPSolver
    {
    protected:
        
        
        //int *auxIndex2;
        //double *auxValues2;
            
        
        //	__methods from Solver __
            
        virtual int __addConstraints(const int nm)  override;
        
        virtual int __addVariables(const int nn, const bool initFree = true)  override;
        
        
    public:
        
        #if OPT_HAVE_GLPK
            glp_prob  *prob;
            
            glp_smcp  *simParam;
            glp_iocp  *intParam;
        #endif
        
        
        OPT_Glpk();
        
        virtual ~OPT_Glpk();
        
        // __methods from Solver __
        
        //virtual int allocateAuxStructures(const int n);
        
        virtual void deallocateMemory() override;
        
        virtual void deallocateSolverEnv() override;
        
        virtual int getNumberOfIterations(long unsigned int &niter) override;
        
        virtual OPT_LISTSOLVERS getSolverCode() override;
        
        virtual int getVariableType( const int index, OPT_VARTYPE &varType ) override;
        
        virtual void initialize() override;
        
        virtual int initSolverEnv(const int maxConstrs = 0, const int maxVars = 0, const int maxQuadNz = 0) override;
        
        virtual int setObjCutLowerBound(const double objLBound) override;
        
        virtual int setObjCutUpperBound(const double objUBound) override;
        
        virtual int setMaxCPUTime(const double time)override;
        
        virtual int setNumberOfThreads(const int nthreads) override;
        
        virtual int setOutputLevel( const int level ) override;
        
        virtual int setRelativeDualTol( const double tol ) override;
    
        virtual int setRelativeOptimalityTol( const double tol ) override;
        
        virtual int setRelativePrimalTol( const double tol ) override;
        
        virtual int setDoubleParameter(const char *param, const double value) override;
        
        virtual int setIntegerParameter(const char *param, const int value ) override;
        
        virtual int setStringParameter(const char *param, const char *value) override;
        
        virtual int setVariableType( const int index, const OPT_VARTYPE varType ) override;
        
        virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true)  override;
        
        virtual int warmUp() override;
        
        
        // __ methods from LPSolver __
        
        
        
        //virtual int addConstraints(const int nm);
        
        //virtual int addVariables(const int nn, const bool initFree);
        
        virtual int generateModelFile(const char* fileName) override;
        
        virtual int getConstraintBounds( const int index, double &lb, double &ub ) override;
        
        virtual int getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value) override;
        
        virtual int getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values) override;
        
        virtual int getNumberOfConstraints(int &m) override;
        
        virtual int getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs) override;
        
        virtual int getNumberOfIntVars(int &n) override;
        
        virtual int getNumberOfVars(int &n) override;
        
        virtual int getObjConstant(double &objConstant) override;
        
        virtual int getObjLinearCoef( const int index, double &value ) override;
        
        virtual int getObjSense(OPT_OPTSENSE &sense) override;
        
        virtual int getVariableBounds(const int index, double &lb, double &ub) override;
        
        virtual int removeVars(const int ninds, const int *indices ) override;
        
        virtual int removeConstraints(const int ninds, const int *indices ) override;
        
        
        virtual int setLinearColumn( const int varIndex, const int nzs, const int* rows, const double* values) override;
        
        
        virtual int resetConstraintLinearPart( const int constrIndex, const int nzs, const int *cols, const double *values ) override;
        
        virtual int setConstraintBounds( const int index, const double lb, const double ub ) override;
        
        
        virtual int setConstraintsLinearCoefs( const int nzs, const int* rows, const int* cols, const double* values ) override;
        
        virtual int setConstraintLinearCoefs( const int constrIndex, const int nzs, const int *cols, const double *values) override;
        
        virtual int setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value) override;
        
        
        virtual int setObjLinearCoef( const int index, const double value ) override;
        
        virtual int setObjLinearCoefs( const int nzs, const int* cols, const double* values ) override;
        
        virtual int setObjLinearPart( const int n, const double *values ) override;
        
        
        virtual void setObjConstant(const double value) override;
        
        virtual int setObjSense( const OPT_OPTSENSE sense ) override;
        
        
        virtual int setVariableBounds( const int index, const double lb, const double ub ) override;
        
        
    protected:
        
        int boundFlag(const double lb, const double ub);
        
    };
    
    
    #if OPT_HAVE_CBC_OR_OSI
        int setRetCodeAndFeasSolByOSISolverTermination( OsiSolverInterface *solver, int &retcode, bool &feasSol );
    #endif
    
    
    /*class to handle a generic solver from OSI interface*/
    class OPT_OpenSolverInterface: public OPT_LPSolver
    {
    protected:
        
        double objConstant;
        
        //	__methods from Solver __
            
        virtual int __addConstraints(const int nm)  override;
        
        virtual int __addVariables(const int nn, const bool initFree = true)  override;
        
        OPT_DynSparseMatrix A; //to store linear coefficients. Unfortunatelly, OSI does not allow change coefficients in the constraint matrix. So, we keep our own storage of this coeficients to allow modifications in efficient way.
        
        double *lc, *uc; //to store bounds os constraints
        
        int indexOfFirstChangedConstraint; 	
        
        //to update the osi solver with the constraints stored in our matrices
        virtual int updateSolver(); 
        
        
    public:
        
        #if OPT_HAVE_CBC_OR_OSI
            OsiSolverInterface* solver;
        #endif
        
        
        OPT_OpenSolverInterface();
        
        virtual ~OPT_OpenSolverInterface();
        
        // __ methods from OPT_Solver __
        
        virtual void deallocateMemory() override;
        
        virtual void deallocateSolverEnv()  override;
        
        virtual int getNumberOfConstraints(int &m) override;
        
        virtual int getNumberOfIterations(long unsigned int& niter) override;
        
        virtual int getNumberOfVars(int &n) override;
        
        virtual OPT_LISTSOLVERS getSolverCode() override;
        
        virtual int getVariableType( const int index, OPT_VARTYPE &varType ) override;
        
        virtual void initialize() override;
        
        //virtual int initSolverEnv(const int maxConstrs = 0, const int maxVars = 0, const int maxQuadNz = 0) override;
        
        virtual int removeConstraints(const int ninds, const int* indices ) override;
        
        
        virtual int removeVars(const int ninds, const int *indices ) override;
        
        virtual int setObjCutLowerBound(const double objLBound) override;
        
        virtual int setObjCutUpperBound(const double objUBound) override;
        
        virtual int setMaxCPUTime(const double time) override;
        
        virtual int setMaxTime(const double time) override;
        
        //virtual int setNumberOfThreads(const int nthreads) override;
        
        virtual int setOutputLevel( const int level ) override;
        
        virtual int setRelativeDualTol( const double tol ) override;
    
        virtual int setRelativeOptimalityTol( const double tol ) override;
        
        virtual int setRelativePrimalTol( const double tol ) override;
        
        virtual int setDoubleParameter(const char *param, const double value) override;
        
        virtual int setIntegerParameter(const char *param, const int value ) override;
        
        virtual int setStringParameter(const char *param, const char *value) override;
        
        virtual int setVariableType( const int index, const OPT_VARTYPE varType ) override;
        
        virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true) override;
        
        
        
        
        // __ methods from OPT_LPSolver __
        
        
        virtual int generateModelFile(const char *fileName) override;
        
        virtual int getConstraintBounds( const int index, double &lb, double &ub ) override;
        
        virtual int getFullConstraintLinearPart(const int constrIndex, double *values) override;
        
        virtual int getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value) override;
        
        virtual int getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values) override;
        
        virtual int getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs) override;
        
        virtual int getObjConstant(double &ObjConstant) override;
        
        virtual int getObjLinearCoef( const int index, double &value ) override;
        
        virtual int getNumberOfIntVars(int &nI) override;
        
        virtual int getObjSense(OPT_OPTSENSE &sense) override;
        
        virtual int getVariableBounds(const int index, double &lb, double &ub)  override;
        
        
        //warning: this method replaces all coefficients in a constraint.
        virtual int resetConstraintLinearPart( const int constrIndex, const int nzs, const int *cols, const double *values )  override;
        
        
        //warning: this method can create a new auxiliary variable because some solvers like cplex and gurobi does not acept dual bounded constraints directly...
        virtual int setConstraintBounds( const int index, const double lb, const double ub ) override;
        
        
        virtual int setConstraintsLinearCoefs( const int nzs, const int *rows, const int *cols, const double *values ) override;
        
        virtual int setConstraintLinearCoefs( const int constrIndex, const int nzs, const int *cols, const double *values)  override;
        
        virtual int setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value)  override;
        
        
        
        virtual int setObjLinearCoef( const int index, const double value )  override;
        
        //virtual int setObjLinearCoefs( const int nzs, const int *cols, const double *values ) override;
        
        //virtual int setObjLinearPart( const int n, const double *values )  override;
        
        
        virtual void setObjConstant(const double value)  override;
        
        virtual int setObjSense( const OPT_OPTSENSE sense )  override;
        
        
        //set first n variable bounds
        //virtual int setnVariablesBounds( const int n, const double *lb, const double *ub )  override;
        
        
        virtual int setVariableBounds( const int index, const double lb, const double ub )  override;
        
        
        //virtual int setVariablesBounds( const int ninds, const int *inds, const double *lb, const double *ub ) override;
    };
    
    
    
    class OPT_Cbc : public OPT_OpenSolverInterface
    {
        
    protected:
        
        //	__methods from Solver __
        
        
    public:
        
        #if OPT_HAVE_CBC
            CbcModel *model;
        #endif
        
        
        OPT_Cbc();
        
        virtual ~OPT_Cbc();
        
        // __ methods from OPT_Solver __
        
        virtual void deallocateMemory() override;
        
        virtual void deallocateSolverEnv()  override;
        
        //virtual int getNumberOfConstraints(int &m) override;
        
        virtual int getNumberOfIterations(long unsigned int& niter) override;
        
        //virtual int getNumberOfVars(int &n) override;
        
        virtual OPT_LISTSOLVERS getSolverCode() override;
        
        //virtual int getVariableType( const int index, OPT_VARTYPE &varType ) override;
        
        virtual void initialize() override;
        
        virtual int initSolverEnv(const int maxConstrs = 0, const int maxVars = 0, const int maxQuadNz = 0) override;
        
        //virtual int removeConstraints(const int ninds, const int* indices ) override;
        
        
        //virtual int removeVars(const int ninds, const int *indices ) override;
        
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
        
        //virtual int setVariableType( const int index, const OPT_VARTYPE varType ) override;
        
        virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true) override;
        
        
        
    };
    
}



#endif
