#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_solvers.hpp"


//class to abstract the proccess of set a lazy constraint
class MRQ_AbstractCallbackSolver
{
public:
	
	virtual void initializeSolverData(void **solverParams) = 0;
	
	virtual int getx(const int milpSolver, void **solverParams, const int n, double *x) = 0;
	
	virtual int setLazyConstraint(const int milpSolver, void **solverParams, const int nz, int *cols, double *vals, const double lb, const double ub) = 0;
};


class MRQ_CplexCallbackSolver : public MRQ_AbstractCallbackSolver
{
public:
	
	CPXCENVptr env;
	void *cbdata;
	int wherefrom;
	int *useraction_p;
	
	virtual void initializeSolverData(void **solverParams) override
	{
		#if OPT_HAVE_CPLEX
			env = *((CPXCENVptr*)  solverParams[0]);
			cbdata = solverParams[1];
			wherefrom = *((int*) solverParams[2]);
			useraction_p = (int*) solverParams[3];
		#endif
	}
	
	
	virtual int getx(const int milpSolver, void **solverParams, const int n, double *x) override
	#if OPT_HAVE_CPLEX
	{
		int r = CPXgetcallbacknodex(env, cbdata, wherefrom, x, 0, n);
		
		if(r != 0)
		{
			#if MRQ_DEBUG_MODE
				MRQ_PRINTERRORNUMBER(r);
			#endif
		}
		return r;
	}
	#else
	{
		return MRQ_LIBRARY_NOT_AVAILABLE;
	}
	#endif
	
	
	virtual int setLazyConstraint(const int milpSolver, void **solverParams, const int nz, int *cols, double *vals, const double lb, const double ub) override
	#if OPT_HAVE_CPLEX
	{
		char sense = MRQ_cplexConstraintSense(lb, ub);
		double rhs =  sense == 'L' ? ub : lb;
		
		
		int r = CPXcutcallbackadd(env, cbdata, wherefrom, nz, rhs, sense, cols, vals, CPX_USECUT_PURGE);
		
		if(r != 0)
		{
			#if MRQ_DEBUG_MODE
				MRQ_PRINTERRORNUMBER(r);
			#endif
		}
		else
			*useraction_p = CPX_CALLBACK_SET; //Tell CPLEX that cuts have been created
		
		return r;
	}
	#else
	{
		return MRQ_LIBRARY_NOT_AVAILABLE;
	}
	#endif
};



class MRQ_AbastractCallbackSolver
{
	
	static inline int getxOnCplex( void **solverParams)
	#if OPT_HAVE_CPLEX
	{
		
	}
	#else
	{
		return MRQ_LIBRARY_NOT_AVAILABLE;
	}
	#endif
	
	static int getx(const int milpSolver, void **solverParams, const int n, double *x)
	{
		int r;
		
		if(milpSolver == MRQ_CPLEX)
		{
		}
		else
		{
			r = MRQ_BAD_PARAMETER_VALUES;
		}
		
		return r;
	}
	
	#if OPT_HAVE_CPLEX
		static inline void getCplexParams(void **solverParams, CPXCENVptr &env, void* &cbdata, int &wherefrom, int* &useraction_p)
		{
			env = *((CPXCENVptr*)  solverParams[0]);
			cbdata = solverParams[1];
			wherefrom = *((int*) solverParams[2]);
			useraction_p = (int*) solverParams[3];
		}
	#endif
	
	
	static inline int setLazyConstraintOnCplex(void **solverParams, const int nz, int *cols, double *vals, const double lb, const double ub)
	#if OPT_HAVE_CPLEX
	{
		CPXCENVptr env = *((CPXCENVptr*)  solverParams[0]);
		void *cbdata = solverParams[1];
		int wherefrom = *((int*) solverParams[2]);
		int *useraction_p = (int*) solverParams[3];
		
		
		char sense = MRQ_cplexConstraintSense(lb, ub);
		double rhs =  sense == 'L' ? ub : lb;
		
		
		int r = CPXcutcallbackadd(env, cbdata, wherefrom, nz, rhs, sense, cols, vals, CPX_USECUT_PURGE);
		
		if(r != 0)
		{
			#if MRQ_DEBUG_MODE
				MRQ_PRINTERRORNUMBER(r);
			#endif
		}
		else
			*useraction_p = CPX_CALLBACK_SET; //Tell CPLEX that cuts have been created
		
		return r;
	}
	#else
	{
		return MRQ_LIBRARY_NOT_AVAILABLE;
	}
	#endif
	
	
	static int setLazyConstraint(const int milpSolver, void **solverParams, const int nz, int *cols, double *vals, const double lb, const double ub)
	{
		int r;
		
		if(milpSolver == MRQ_CPLEX)
		{
			r = setLazyConstraintOnCplex(solverParams, nz, cols, vals, lb, ub);
		}
		else
		{
			r = MRQ_BAD_PARAMETER_VALUES;
		}
		
		return r;
	}
	
	
};
