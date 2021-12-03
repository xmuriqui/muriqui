
#include <new>
#include "MRQ_dataStructures.hpp"
#include "MRQ_tools.hpp"

using namespace std;

using namespace muriqui;




MRQ_HistorySolution::MRQ_HistorySolution()
{
	nvars = iter = -1;
	time = cputime = -1.0;
	sol = NULL;
}


/*MRQ_HistorySolution::MRQ_HistorySolution(const int n, const long unsigned int iter, const double time, const double cputime, const double *sol, const double objF)
{
	int i;
	
	i = allocateSol(n);
	if(i != 0)
	{
		return;
	}
	
	this->nvars = n;
	this->iter = iter;
	this->time = time;
	this->cputime = cputime;
	
	for(i = 0; i < n; i++)
		this->sol[i] = sol[i];
	
	this->objF = objF;
	
}*/


MRQ_HistorySolution::~MRQ_HistorySolution()
{
	freeSolution();
}


int MRQ_HistorySolution::allocateSol(const int n)
{
	sol = (double *) malloc( n * sizeof(double) );
	if( !sol )
	{
		return MRQ_MEMORY_ERROR;
	}
	
	return 0;
}


void MRQ_HistorySolution::freeSolution()
{
	MRQ_secFree(sol);
}


int MRQ_HistorySolution::getnvars()
{
	return nvars;
}

int MRQ_HistorySolution::getiter()
{
	return iter;
}

double MRQ_HistorySolution::gettime()
{
	return time;
}

double MRQ_HistorySolution::getcputime()
{
	return cputime;
}


double MRQ_HistorySolution::getobjvalue()
{
	return objF;
}



int MRQ_HistorySolution::getsolution(double *solution)
{
	if( sol == NULL )
	{
	    return MRQ_VALUE_ERROR;
	}
	else
	{
		MRQ_copyArray(nvars, sol, solution);
		return 0;
	}
}


int MRQ_HistorySolution::set(const int n, const long unsigned int iter, const double time, const double cputime, const double *sol, const double objF)
{
	const int ret = allocateSol(n);
	if(ret != 0)
	{
		#if MRQ_DEBUG_MODE
			MRQ_PRINTERRORNUMBER(ret);
		#endif
		
		return ret;
	}
	
	this->nvars = n;
	this->iter = iter;
	this->time = time;
	this->cputime = cputime;
	
	//for(int i = 0; i < n; i++)
		//this->sol[i] = sol[i];
	
	MRQ_copyArray(n, sol, this->sol);
	
	
	this->objF = objF;
	
	return 0;
}






MRQ_SolutionHistory::MRQ_SolutionHistory()
{
	nsols = 0;
	hsols = NULL;
}


MRQ_SolutionHistory::~MRQ_SolutionHistory()
{
	desallocate();
}


void MRQ_SolutionHistory::desallocate()
{
	for(unsigned int i = 0; i < nsols; i++)
	{
		if( hsols[i] )
			delete hsols[i];
	}
	
	MRQ_secFree(hsols);
	/*if(hsols)
	{
		free(hsols);
		hsols = NULL;
	} */
	
	nsols = 0;
}


unsigned int MRQ_SolutionHistory::getnsols() const
{
	return nsols;
}


int MRQ_SolutionHistory::addSolution(const int n, const long unsigned int iter, const double time, const double cputime, const double* sol, const double objF)
{
	MRQ_HistorySolution **haux;
	
	haux = (MRQ_HistorySolution **) realloc( hsols, (nsols + 1) * sizeof(MRQ_HistorySolution *) );
	
	if(haux)
	{
		hsols = haux;
		hsols[nsols] = new (std::nothrow) MRQ_HistorySolution();
		
		if(!hsols[nsols])
		{
			#if MRQ_DEBUG_MODE
				MRQ_PRINTMEMERROR;
			#endif
			
			return MRQ_MEMORY_ERROR;
		}
		
		const int r = hsols[nsols]->set(n, iter, time, cputime, sol, objF);
		
		if(r != 0)
		{
			#if MRQ_DEBUG_MODE
				MRQ_PRINTERRORNUMBER(r);
			#endif
			
			return r;
		}
		
		
		nsols++;
		return 0;
	}
	else
	{
		#if MRQ_DEBUG_MODE
			MRQ_PRINTMEMERROR;
		#endif
		
		return MRQ_MEMORY_ERROR;
	}
}


//it is only a pointer, not a copy. Do not free. If you want a copy of solution use the method getsol of MRQ_HistorySolution pointer
MRQ_HistorySolution * MRQ_SolutionHistory::getHistSolPointer(const unsigned int index)
{
	return hsols[index];
}












