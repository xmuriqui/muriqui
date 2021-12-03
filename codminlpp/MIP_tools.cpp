

#include <math.h>
#include <cstdlib>


#include "MIP_minlpProblem.hpp"
#include "MIP_tools.hpp"



using namespace minlpproblem;






#if 0
	//this classe try get a feasible rouded solution for a integer problem loking for linear constraints
	class MIP_SmartRouding
	{
		int ntargetConstr; //number of target constraints, i.e., constraints that we consider to do the rouding
		int *targetConstr; //index of traget constraints
		
	public:
		
		MIP_SmartRouding();
		
		~MIP_SmartRouding();
		
		int buildTargetConstrArray(MIP_MINLPProb &prob);
		
		void desallocate();
		
		int smartRounding(MIP_MINLPProb &prob, const double *lx, const double *ux, const double *sol, double *rsol);
		
	};
	


MIP_SmartRouding::MIP_SmartRouding()
{
	ntargetConstr = -1;
	targetConstr = NULL;
}


MIP_SmartRouding::~MIP_SmartRouding()
{
	desallocate();
}


void MIP_SmartRouding::desallocate()
{
	MIP_secFree( targetConstr );
	ntargetConstr = -1;
}


int MIP_SmartRouding::buildTargetConstrArray( MIP_MINLPProb &prob )
{
	const int m = prob.m;
	const bool *nlConstr = prob.nlConstr;
	const int *xtype = prob.xtype;
	const double *lc = prob.lc;
	const double *uc = prob.uc;
	MIP_SparseMatrix *QC = prob.QC;
	MIP_SparseMatrix &A = prob.A;
	
	
	
	ntargetConstr = 0;
	
	if( targetConstr )
		free(targetConstr);
	
	targetConstr = (int*) malloc( m *sizeof(int) );
	if( !targetConstr )
	{
		#if MIP_DEBUG_MODE
			MIP_PRINTMEMERROR;
		#endif
		
		return MIP_MEMORY_ERROR; 
	}
	
	
	for(int i = 0; i < m; i++)
	{
		if( nlConstr[i] || QC[i].getNumberOfElements() > 0 )
			continue;
		
		
		if( A[i].getNumberOfElements() == 0 || (lc[i] <= -MIP_INFINITY && uc[i] >= MIP_INFINITY) )
			continue;
		
		//linear constraint. Now, we check if this constraint only have integer variables
		bool onlyIntVars = true;
		MIP_SparseRow &a = A[i];
		const unsigned int nel = a.getNumberOfElements();
		
		for( unsigned int j = 0; j < nel; j++ )
		{
			if( !MIP_isIntegerType( xtype[ a[j].getColumn() ] ) )
			{
				onlyIntVars = false;
				break;
			}
		}
		
		if( onlyIntVars )
		{
			targetConstr[ntargetConstr] = i;
			ntargetConstr++;
		}
	}
	
	
	
	return 0;
}



int MIP_SmartRouding::smartRounding(MIP_MINLPProb &prob, const double *lx, const double *ux, const double *sol, double *rsol)
{
	const int n = prob.n;
	const int m = prob.m;
	const double *mylx = lx ? lx : prob.lx;
	const double *myux = ux ? ux : prob.ux;
	const double *lc = prob.lc;
	const double *uc = prob.uc;
	
	MIP_SparseMatrix &A = prob.A;
	
	const int ntargetConstr = this->ntargetConstr;
	const int *targetConstr = this->targetConstr;
	
	
	//we use INFINITY like a falg
	//MIP_setAllArray<double>( n, rsol, INFINITY );
	MIP_copyArray(n, sol, rsol);
	
	
	for( int i = 0; i < n; i++ )
	{
		if( rsol[i] < mylx[i] )
			rsol[i] = mylx[i];
		else if( rsol[i] > myux[i] )
			rsol[i] = myux[i];
	}
	
	
}


#endif






MIP_Mutex::MIP_Mutex()
{
	initialize();
}


void MIP_Mutex::initialize()
{
	#if MIP_CPP_MULTITHREADING
		
	#else	
		#if MIP_OMP_MULTITHREADING
			omp_init_lock(&mutex);
		#endif
	#endif
}


int MIP_Mutex::lock( )
{
	
	#if MIP_CPP_MULTITHREADING
		//if( nthreads > 1 )
			mymutex.lock();
	#else
		#if MIP_OMP_MULTITHREADING
			//if( nthreads > 1 )
				omp_set_lock(&mutex);
		#endif
	#endif
	
	return 0;
}


int MIP_Mutex::tryLock(  )
{
	#if MIP_CPP_MULTITHREADING
		//if( nthreads > 1 )
			return mymutex.try_lock() == true ? 0 : MIP_UNDEFINED_ERROR;
	#else
		#if MIP_OMP_MULTITHREADING
			//if( nthreads > 1 )
				return omp_test_lock(&mutex) == 1 ? 0 : MIP_UNDEFINED_ERROR;
		#endif
	#endif
	
	return 0;
}


int MIP_Mutex::unlock( )
{
	#if MIP_CPP_MULTITHREADING
		//if( nthreads > 1 )
			mymutex.unlock();
	#else
		#if MIP_OMP_MULTITHREADING
			//if( nthreads > 1 )
				omp_unset_lock(&mutex);
		#endif
	#endif
	
	return 0;
}


void MIP_Mutex::destroy()
{
	#if MIP_OMP_MULTITHREADING
		omp_destroy_lock(&mutex);
	#endif
}


MIP_Mutex::~MIP_Mutex()
{
	destroy();
}










































