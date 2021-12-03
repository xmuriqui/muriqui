

#ifndef OPT_NLPSOLVERS_HPP
#define OPT_NLPSOLVERS_HPP

#include <vector>
#include <string>
#include <unordered_set>


#if OPT_HAVE_MOSEK
	#include "mosek.h"
#endif


#if OPT_HAVE_IPOPT
	#include "IpIpoptApplication.hpp"
	#include "IpTNLP.hpp"
#endif


#if OPT_HAVE_WORHP
	#include "worhp.h"
#endif


#if OPT_HAVE_KNITRO
	#include "knitro.h"
#endif


#if OPT_HAVE_IQUAD
	//#include "iquad.hpp"
#endif



#include "OPT_basesolvers.hpp"



namespace iquad
{
	class IQD_BranchAndBound;
}



namespace optsolvers{
	
	#if OPT_HAVE_MOSEK
		
		int MSKAPI OPT_mosekStruc(void    *nlhandle, MSKintt  *numgrdobjnz, MSKidxt  *grdobjsub, MSKidxt  i, int      *convali, MSKidxt  *grdconinz, MSKidxt  *grdconisub, MSKintt  yo, MSKintt  numycnz, MSKCONST MSKidxt  *ycsub, MSKintt  maxnumhesnz, MSKintt  *numhesnz, MSKidxt  *hessubi, MSKidxt  *hessubj);
		
		
		
		int MSKAPI OPT_mosekEval(void      *nlhandle, MSKCONST double    *xx, double    yo, MSKCONST double    *yc, double    *objval, MSKintt   *numgrdobjnz, MSKidxt   *grdobjsub, double    *grdobjval, MSKintt   numi, MSKCONST MSKidxt   *subi, double    *conval, MSKCONST MSKintt   *grdconptrb, MSKCONST MSKintt   *grdconptre, MSKCONST MSKidxt   *grdconsub, double    *grdconval, double    *grdlag, MSKintt   maxnumhesnz, MSKintt   *numhesnz, MSKidxt   *hessubi, MSKidxt   *hessubj, double    *hesval);
		
		
	#endif
	
	
	
	class OPT_Mosek : public OPT_NLPSolver
	{
	protected:
		
		bool nlObj; //flag to nonlinear objective function
		bool hasNLConstrs;
		bool *nlConstr; //flags to nonlinear constraints...
		double nlObjFactor; //we just use that due to obj factor  
		
		
		// __ methods from LPSolver __
		
		
		virtual int __addConstraints(const int nm) override;
		
		virtual int __addVariables(const int nn, const bool initFree = true)  override;
		
		
		
		// __ methods from OPT_NLPSolver __
		
		virtual int __removeConstraints(const int ninds, const int* indices)  override;
		
		
	public:
		
		#if OPT_HAVE_MOSEK
			MSKenv_t     env;
		    MSKtask_t    task;
		#endif
		
		
		OPT_Mosek();
		
		virtual ~OPT_Mosek();
		
		// __methods from Solver __
		
		virtual void deallocateSolverEnv() override;
		
		//virtual void getDualSolution( double *dualConstrs, double *dualVarBounds, const bool correctSignal = true );
		
		virtual int getNumberOfIterations(unsigned int &niter) override;
		
		virtual OPT_LISTSOLVERS getSolverCode() override;
		
		virtual int getVariableType( const int index, OPT_VARTYPE &varType ) override;
		
		virtual void initialize() override;
		
		virtual int initSolverEnv(const int maxConstrs = 0, const int maxVars = 0, const int maxQuadNz = 0) override;
		
		virtual int removeVars(const int ninds, const int *indices) override;
		
		virtual int setLowerObjCut(const double objLBound) override;
		
		virtual int setUpperObjCut(const double objUBound) override;
		
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
		
		virtual int getConstraintBounds( const int index, double &lb, double &ub) override;
		
		virtual int getLinearConstraintCoef( const int constrIndex, const int varIndex, double &value) override;
		
		virtual int getLinearConstraintPart(const int constrIndex, int &nzs, int *cols, double *values) override;
		
		virtual int getLinearObjCoef( const int index, double &value ) override;
		
		virtual int getNumberOfConstraints(int &m) override;
		
		virtual int getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs) override;
		
		virtual int getNumberOfIntVars(int &n) override;
		
		virtual int getNumberOfVars(int &n) override;
		
		virtual int getObjSense(OPT_OPTSENSE &sense) override;
		
		virtual int getVariableBounds(const int index, double &lb, double &ub) override;
		
		virtual int setLinearColumn( const int varIndex, const int nzs, const int* rows, const double* values) override;
		
		
		virtual int resetLinearConstraintPart( const int constrIndex, const int nzs, const int *cols, const double *values ) override;
		
		virtual int setConstraintBounds( const int index, const double lb, const double ub ) override;
		
		//virtual int setConstraintLowerBound( const int index, const double lb );
		
		//virtual int setConstraintUpperBound( const int index, const double ub );
		
		
		virtual int setLinearConstraintsCoefs( const int nzs, const int* rows, const int* cols, const double* values ) override;
		
		virtual int setLinearConstraintCoefs( const int constrIndex, const int nzs, const int *cols, const double *values) override;
		
		virtual int setLinearConstraintCoef( const int constrIndex, const int varIndex, const double value) override;
		
		
		virtual int setLinearObjCoef( const int index, const double value ) override;
		
		virtual int setLinearObjCoefs( const int nzs, const int* cols, const double* values ) override;
		
		virtual int setLinearObjFunction( const int n, const double *values ) override;
		
		
		virtual void setObjConstant(const double value) override;
		
		virtual int setObjSense( const OPT_OPTSENSE sense ) override;
		
		//set first n variable bounds
		virtual int setnVariablesBounds(const int n, const double *lb, const double *ub) override;
		
		virtual int setVariableBounds(const int index, const double lb, const double ub) override;
		
		virtual int setVariablesBounds(const int ninds, const int *inds, const double *lb, const double *ub) override;
		
		//virtual int setVariableLowerBound( const int index, const double lb );
		
		//virtual int setVariableUpperBound( const int index, const double ub );
		
		
		
		//__ methods from QPSolver __
		
		
		virtual int getNumberOfQuadObjTerms(int &nzs) override;
		
		virtual int getQuadObjTerm(const int row, const int col, double &value)  override;
		
		virtual int getQuadObjPart(int &nzs, int *rows, int *cols, double *values) override;
		
		virtual int setQuadObjCoef(const int row, const int col, const double value) override;
		
		virtual int setQuadObjMatrix(const int nzs, const int *rows, const int *cols, const double *values) override;
		
		//set by compressed row format
		//virtual int setQuadObjMatrix(const int *rowStart, const int *cols, const double *values) override;
		
		// __ methods from QCPSolver __
		
		virtual int getNumberOfConstraintQuadTerms( const int index, int &nzs) override;
		
		
		virtual int getQuadConstraintMatrix(const int index, int &nzs, int *rows, int *cols, double *values) override;
		
		virtual int setQuadConstraintMatrix( const int index, const int nzs, const int *qrows, const int *qcols, const double *qvalues ) override;
		
		
		// __ methods from OPT_NLPSolver __
		
		
		virtual int getConstrNLFlag(const int index, bool &flag) override;
		
		virtual int getNumberOfNLConstraints(int& mnl) override;
		
		virtual int setConstrNLFlag(const int index, const bool flag) override;
		
		virtual int setInitialSolution(const double* x, const double* dualConstrs, const double* dualVars) override;
		
		virtual int setProblemFrom(const minlpproblem::MIP_MINLPProb &prob, const bool setObj = true, const bool setConstrs = true, const bool setVarBounds = true, const bool setVarType = true, const int naddvars = 0, const int naddconstrs = 0) override;
		
		virtual int setObjNLFlag(const bool flag) override;
		
		
	protected:
		
		//internal methods...
		
		
		virtual int allocateConstrStructures(const int m) override;
		
		void updatehasNLConstr();
		
		
		#if OPT_HAVE_MOSEK
			MSKboundkeye boundKeye( const double lb, const double ub );
			
			int getSolution(const MSKsoltypee solType, const bool storeSol, const bool storeConstrs, const bool storeDualSol);
			
		
		
		
			friend int MSKAPI OPT_mosekStruc(void    *nlhandle, MSKintt  *numgrdobjnz, MSKidxt  *grdobjsub, MSKidxt  i, int      *convali, MSKidxt  *grdconinz, MSKidxt  *grdconisub, MSKintt  yo, MSKintt  numycnz, MSKCONST MSKidxt  *ycsub, MSKintt  maxnumhesnz, MSKintt  *numhesnz, MSKidxt  *hessubi, MSKidxt  *hessubj);
		
		
			friend int MSKAPI OPT_mosekEval(void      *nlhandle, MSKCONST double    *xx, double    yo, MSKCONST double    *yc, double    *objval, MSKintt   *numgrdobjnz, MSKidxt   *grdobjsub, double    *grdobjval, MSKintt   numi, MSKCONST MSKidxt   *subi, double    *conval, MSKCONST MSKintt   *grdconptrb, MSKCONST MSKintt   *grdconptre, MSKCONST MSKidxt   *grdconsub, double    *grdconval, double    *grdlag, MSKintt   maxnumhesnz, MSKintt   *numhesnz, MSKidxt   *hessubi, MSKidxt   *hessubj, double    *hesval);
		
		#endif
	};
	
	
	class OPT_MyProblemToIpopt;
	
	#if OPT_HAVE_IPOPT
	
		class OPT_IpoptIntermediateCallback
		{
		public:
			virtual bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value, Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq) = 0;
		};
	
	#endif
	
	
	class OPT_Ipopt : public OPT_MyNLPSolver
	{
		
	protected:
		
		bool jac_c_const;
		bool jac_d_const;
		bool hessian_const;
		
		bool enforceReoptimization;
		
		int mquad;
		unsigned int numberOfIterations;
		
		
		//int *quadIndex;
		//int *nqcons; //nqcons[i] has the number of matrices in QC having some elemnt in row i
		//int *indqcons; //indqconsi has the last constraint having some element in row i of QC
		
		
		/*unsigned int *jacIndexSh; //matrix for getting indexes for Ipopt jacobian... 
		unsigned int *hessIndexSh; ////matrix for getting indexes for Ipopt hessian...
		unsigned int *jacIndexBase; //m+1-array with base for indices in Ipopt Jacobian for each constraint
		unsigned int *hessIndexBase; //n+1-array with base for indices in Ipopt Hessian*/
		
		
		/*unsigned int **Jindex;
		unsigned int **Aindex;
		unsigned int ***QCindex;
		unsigned int **QCrowindex;
		unsigned int **lagHindex; */
		
		
		
		virtual void __deallocateSolverEnv() override;
		
		//virtual int allocateVarStructures(const int size);
		
		//virtual int allocateConstrStructures(const int m);
		
		virtual bool getMinusLambdaOnLagran() override; //return attribute to solver convention of minus signal before lambda on lagrangian definition. We use that because on the dual solution, we follow the convention of linear programming (I think linear programming solvers consider this minus signal). So, we adopt in the solution of optsolver, this minus signal in the method getDualSolution
		
		
	public:
		
		
		#if OPT_HAVE_IPOPT
			Ipopt::SmartPtr<OPT_MyProblemToIpopt> mynlp;
			Ipopt::SmartPtr<Ipopt::IpoptApplication> app;
			OPT_IpoptIntermediateCallback *ipoptInterCallback;
		#endif
		
		OPT_Ipopt();
		
		virtual ~OPT_Ipopt();
		
		// __methods from Solver __
		
		virtual void deallocateMemory() override;
		
		virtual int getNumberOfIterations(unsigned int &niter) override;
		
		virtual OPT_LISTSOLVERS getSolverCode() override;
		
		virtual int getVariableType( const int index, OPT_VARTYPE &varType) override;
		
		virtual void initialize() override;
		
		virtual int initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz = 0) override;
		
		virtual int setMaxCPUTime(const double time) override;
		
		virtual int setNumberOfThreads(const int nthreads) override;
		
		virtual int setOutputLevel(const int level) override;
		
		virtual int setRelativeDualTol( const double tol ) override;
		
		virtual int setRelativeOptimalityTol(const double tol) override;
		
		virtual int setRelativePrimalTol(const double tol) override;
		
		virtual int setDoubleParameter(const char *param, const double value) override;
		
		virtual int setIntegerParameter(const char *param, const int value ) override;
		
		virtual int setStringParameter(const char *param, const char *value) override;
		
		virtual int setVariableType( const int index, const OPT_VARTYPE varType ) override;
		
		virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true) override;
		
		virtual int warmUp() override;
		
		
		// __methods from NLPSolver __
		
		
		//virtual int setInitialSolution(const double* x, const double* dualConstrs, const double* dualVars);
		
		#if OPT_HAVE_IPOPT
			void setIntermediateCallbackPointer( OPT_IpoptIntermediateCallback *intermediateCallback);
		#endif
		
		
	/* protected:
		
		int allocateIpoptIndexStructures();
		
		void desallocateIpoptIndexStructures();
		
		int setHessIndexRow( const int rowIndex, int &nzs );
		
		int setJacIndexRow( const int rowIndex, int& nzs ); */
		
		
		
		friend class OPT_MyProblemToIpopt;
		
	};
	
	
	
	class OPT_Worhp: public OPT_MyNLPSolver
	{
		//int nzJac;
		//int nzHess;
		OPT_UIntSparseMatrix uJac; //unsigned int sparse matrix having indexes for jacobian on worph (we need to sort by column). Position (i, j) has the index of (i,j) in worph data structure
		OPT_UIntSparseMatrix uHess; //sparse matrix having indexes for hessian on worph (we need to sort by column). Position (i, j) has the index of (i,j) in worph data structure
		
		
		virtual void __deallocateSolverEnv() override;
		
		virtual bool getMinusLambdaOnLagran() override; //return attribute to solver convention of minus signal before lambda on lagrangian definition. We use that because on the dual solution, we follow the convention of linear programming (I think linear programming solvers consider this minus signal). So, we adopt in the solution of optsolver, this minus signal in the method getDualSolution
		
	public:
		
		#if OPT_HAVE_WORHP
			//OptVar    opt;
			Workspace wsp; //we set it here to get information about solution. Annyway, we have to set initialised atributte as false in solve method
			Params    par;
			//Control   cnt;
		#endif
		
		OPT_Worhp();
    		
		virtual ~OPT_Worhp();
		
		// __methods from Solver __
		
		virtual void deallocateMemory() override;
		
		virtual int getNumberOfIterations(unsigned int &niter) override;
		
		virtual OPT_LISTSOLVERS getSolverCode() override;
		
		virtual int getVariableType( const int index, OPT_VARTYPE &varType ) override;
		
		virtual void initialize() override;
		
		virtual int initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz = 0) override;
		
		virtual int setMaxCPUTime(const double time) override;
		
		virtual int setNumberOfThreads(const int nthreads) override;
		
		virtual int setOutputLevel(const int level) override;
		
		virtual int setRelativeDualTol( const double tol ) override;
		
		virtual int setRelativeOptimalityTol( const double tol ) override;
		
		virtual int setRelativePrimalTol( const double tol ) override;
		
		virtual int setDoubleParameter(const char *param, const double value) override;
		
		virtual int setIntegerParameter(const char *param, const int value ) override;
		
		virtual int setStringParameter(const char *param, const char *value) override;
		
		virtual int setVariableType( const int index, const OPT_VARTYPE varType ) override;
		
		virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true) override;
		
		virtual int warmUp() override;
		
		
		// __methods from NLPSolver __
		
		
		//virtual int setInitialSolution(const double* x, const double* dualConstrs, const double* dualVars);
		
		
	protected:
		
		int allocateAuxDerivativeIndexStructures( );
		
		void desallocateAuxDerivativeIndexStructures();
		
		
		int setuJacStructRow( const int rowIndex );
		
		int setuHessStructRow( const int rowIndex );
		
		void setuJacIndices();
		
		void setuHessIndices();
		
		//void setSPMonHessianValues( const int n, minlpproblem::MIP_SparseMatrix &M, const double factor, double *values );
	};
	
	#if OPT_HAVE_KNITRO
	
		int OPT_knitroObjConstrF (const int evalRequestCode, const int n, const int m, const int nnzJ, const int nnzH, const double * const  x, const double * const  lambda, double * const  obj, double * const  c, double * const  objGrad, double * const  jac, double * const  hessian, double * const  hessVector, void *userParams);
		
		int OPT_knitroGrads (const int evalRequestCode, const int n, const int m, const int nnzJ, const int nnzH, const double * const  x, const double * const  lambda, double * const  obj, double * const  c, double * const  objGrad, double * const  jac, double * const  hessian, double * const  hessVector, void *userParams);
		
		int OPT_knitroHess(const int evalRequestCode, const int n, const int m, const int nnzJ, const int nnzH, const double * const  x, const double * const  lambda, double * const  obj, double * const  c, double * const  objGrad, double * const  jac, double * const  hessian, double * const  hessVector, void *userParams);
	#endif
	
	
	class OPT_Knitro : public OPT_MyNLPSolver
	{
		
	protected:
		
		//bool jac_c_const;
		//bool jac_d_const;
		//bool hessian_const;
		
		//bool enforceReoptimization;
		
		//int mquad;
		//int *quadIndex;
		
		/*unsigned int *jacIndexSh; //matrix for getting indexes for Ipopt jacobian... 
		unsigned int *hessIndexSh; ////matrix for getting indexes for Ipopt hessian...
		unsigned int *jacIndexBase; //m+1-array with base for indices in Ipopt Jacobian for each constraint
		unsigned int *hessIndexBase; //n+1-array with base for indices in Ipopt Hessian*/
		
		
		/*unsigned int **Jindex;
		unsigned int **Aindex;
		unsigned int ***QCindex;
		unsigned int **QCrowindex;
		unsigned int **lagHindex; */
		
		
		//bool *allCEvalTrue;
		
		bool problemInitialized;
		
		bool consBoundChg;
		bool newInitPoint;
		bool paramChg;
		int mquad;
		int *quadIndex;
		int *nqcons; //nqcons[i] has the number of matrices in QC having some elemnt in row i
		int *indqcons; //indqconsi has the last constraint having some element in row i of QC
		double *auxLambda;
		
		
		virtual void __deallocateSolverEnv() override;
		
		//virtual int allocateVarStructures(const int size);
		
		//virtual int allocateConstrStructures(const int m);
		
		virtual bool getMinusLambdaOnLagran() override; //return attribute to solver convention of minus signal before lambda on lagrangian definition. We use that because on the dual solution, we follow the convention of linear programming (I think linear programming solvers consider this minus signal). So, we adopt in the solution of optsolver, this minus signal in the method getDualSolution
		
		
	public:
		
		
		void  *knitroContext; //variable to knitro set problem
		int objFnType; //user parameter to mixed integer optimization. It must be set as KTR_FNTYPE_UNCERTAIN, KTR_FNTYPE_CONVEX or KTR_FNTYPE_NONCONVEX for objective function
		int *cFnType; //user parameter to mixed integer optimization. It must be set as KTR_FNTYPE_UNCERTAIN, KTR_FNTYPE_CONVEX or KTR_FNTYPE_NONCONVEX for each constraint
		
		
		OPT_Knitro();
		
		virtual ~OPT_Knitro();
		
		// __methods from Solver __
		
		virtual int allocateConstrStructures(const int m) override;
		
		virtual void deallocateMemory() override;
		
		virtual int getNumberOfIterations(unsigned int &niter) override;
		
		virtual OPT_LISTSOLVERS getSolverCode() override;
		
		//virtual int getVariableType( const int index, OPT_VARTYPE &varType );
		
		virtual void initialize() override;
		
		virtual int initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz = 0) override;
		
		virtual int setConstraintBounds( const int index, const double lb, const double ub ) override;
		
		virtual int setInitialSolution(const double* x, const double* dualConstrs, const double* dualVars) override;
		
		virtual int setLowerObjCut(const double objLBound) override;
		
		virtual int setUpperObjCut(const double objUBound) override;
		
		virtual int setMaxCPUTime(const double time) override;
		
		virtual int setMaxTime(const double time) override;
		
		virtual int setNumberOfThreads(const int nthreads) override;
		
		virtual int setOutputLevel(const int level) override;
		
		virtual int setRelativeDualTol( const double tol ) override;
		
		virtual int setRelativeOptimalityTol( const double tol ) override;
		
		virtual int setRelativePrimalTol( const double tol ) override;
		
		virtual int setDoubleParameter(const char *param, const double value) override;
		
		virtual int setIntegerParameter(const char *param, const int value ) override;
		
		virtual int setStringParameter(const char *param, const char *value) override;
		
		//virtual int setVariableType( const int index, const OPT_VARTYPE varType );
		
		virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true) override;
		
		//virtual int warmUp();
		
		
		// __ methods from NLPSolver __
		
		
		//virtual int setInitialSolution(const double* x, const double* dualConstrs, const double* dualVars);
		
		
		// __ my methods __
		
		void setObjFnType(const int value);
		
		void setCFnType(const int *values);
		
		
		
		
		friend int OPT_knitroObjConstrF (const int             evalRequestCode, const int n, const int             m, const int nnzJ, const int             nnzH, const double * const  x, const double * const  lambda, double * const  obj, double * const  c, double * const  objGrad, double * const  jac, double * const  hessian, double * const  hessVector, void *userParams);
		
		friend int OPT_knitroGrads (const int evalRequestCode, const int n, const int m, const int nnzJ, const int nnzH, const double * const  x, const double * const  lambda, double * const  obj, double * const  c, double * const  objGrad, double * const  jac, double * const  hessian, double * const  hessVector, void *userParams);
		
		friend int OPT_knitroHess(const int evalRequestCode, const int n, const int m, const int nnzJ, const int nnzH, const double * const  x, const double * const  lambda, double * const  obj, double * const  c, double * const  objGrad, double * const  jac, double * const  hessian, double * const  hessVector, void *userParams);
		
	};
	
	
	
	
	class OPT_Algencan : public OPT_MyNLPSolver
	{
	protected:
		
		
		
		virtual void __deallocateSolverEnv() override;
		
		
	public:
		
		bool checkder;
		
		std::string outputfnm;
		std::string specfnm;
		
		double epsfeas; //Feasibility tolerance for the sup-norm of the constraints. (Ignored in the unconstrained and bound-constrained cases.)
		double epsopt; //Optimality tolerance for the sup-norm of the projected gradient of the Augmented Lagrangian in the constrained case and the sup-norm of the projected gradient of the objective function in the unconstrained and the bound-constrained cases.
		std::vector<std::string> params;
		
		OPT_Algencan();
		
		virtual ~OPT_Algencan();
		
		void clearParameterList();
		
		virtual int getNumberOfIterations(unsigned int &niter) override;
		
		virtual OPT_LISTSOLVERS getSolverCode() override;
		
		//that function can receive estimatives of maximum number of variables, constraints and nonzeros in quadratic terms matrices, but it is not mandatory...
		virtual int initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz = 0) override;
		
		virtual int setMaxCPUTime(const double time) override;
		
		virtual int setNumberOfThreads(const int nthreads) override;
		
		virtual int setOutputLevel( const int level ) override;
		
		virtual int setRelativeDualTol( const double tol ) override;
	
		virtual int setRelativeOptimalityTol( const double tol ) override;
		
		virtual int setRelativePrimalTol( const double tol ) override;
		
		virtual int setDoubleParameter(const char *param, const double value) override;
		
		virtual int setIntegerParameter(const char *param, const int value ) override;
		
		
		virtual int setStringParameter(const char *param, const char *value) override;
		
		virtual int setVariableType(const int index, const OPT_VARTYPE varType) override;
		
		virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true ) override;
		
		
		
		//class to store 
		class OPT_MyData
		{
		public:
			
			unsigned int thnumber;
			int mquad;
			
			bool *constrEval;
			int *quadIndex;
			double *auxVars;
			double *auxConstr;
			OPT_MINLPProb *prob;
			std::unordered_set<int> *nzRowsLagH;
			std::unordered_set<int> *colsNzRowHess;
			
			
			OPT_MyData()
			{
				thnumber = -1;
				mquad = 0;
				
				constrEval = NULL;
				quadIndex = NULL;
				auxVars = NULL;
				auxConstr = NULL;
				prob = NULL;
				nzRowsLagH = NULL;
				colsNzRowHess = NULL;
			}
		};
		
		
		
		
		
	protected:
		
		OPT_Algencan::OPT_MyData data;
		
		void buildJacIndex();
		
		
		
		
		int inform;
		double cnorm, nlpsupn, snorm;
		
		double efacc, eoacc; //Feasibility and optimality levels, below which a newton-based acceleration proccess is lauched
		double efstin, eostin; //Parameters efstain and eostain are related to feasibility and optimality tolerances, respectively, and are used by algencan to stop the execution declaring that an infeasible stationary point of the sum of the squared infeasibilities was found
	};
	
	
	
	
	
	class OPT_Iquad : public OPT_MyNLPSolver
	{
		
	protected:
		
		//methods from OPT_MyNLPSolver
		virtual void __deallocateSolverEnv() override;
		
		
		virtual bool getMinusLambdaOnLagran() override; //return attribute to solver convention of minus signal before lambda on lagrangian definition. We use that because on the dual solution, we follow the convention of linear programming (I think linear programming solvers consider this minus signal). So, we adopt in the solution of optsolver, this minus signal in the method getDualSolution
		
		
	public:
		
		#if OPT_HAVE_IQUAD
			iquad::IQD_BranchAndBound *bb;
		#endif
		
		//methods from OPT_solver
		
		OPT_Iquad();
		
		virtual ~OPT_Iquad();
		
		virtual int getNumberOfIterations(unsigned int &niter) override;
		
		virtual OPT_LISTSOLVERS getSolverCode() override;
		
		virtual void initialize() override;
		
		virtual int initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz = 0) override;
		
		virtual int setLowerObjCut(const double objLBound) override;
		
		virtual int setUpperObjCut(const double objUBound) override;
		
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
		
		
		virtual int solve(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true ) override;
		
		int solveWParams(const bool resetSol = true, const bool storeSol = true, const bool storeConstrs = true, const bool storeDualSol = true, OPT_GeneralSolverParams *subSolverParams = NULL, OPT_GeneralSolverParams *sdpParams = NULL ) ;
	};
	
	
	
}



#endif
