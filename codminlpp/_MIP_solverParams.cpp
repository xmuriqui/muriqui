
#include <new>
#include "MIP_minlpProblem.hpp"



using namespace std;

using namespace minlpproblem;




int MIP_GeneralSolverParams::addIntegerParameter(const char *name, const int value)
{
    MIP_SolverParam<int> *p;
    
    p = new (nothrow) MIP_SolverParam<int>;
    if(!p)
		return MIP_MEMORY_ERROR;
    p->setParam(name, value);
    
    intParams.push_back(p);
    return 0;
}

int MIP_GeneralSolverParams::addDoubleParameter(const char* name, const double value)
{
    MIP_SolverParam<double> *p;
    
    p = new (nothrow) MIP_SolverParam<double>;
    if(!p)
		return MIP_MEMORY_ERROR;
    p->setParam(name, value);
    
    dblParams.push_back(p);
	
    return 0;
}

int MIP_GeneralSolverParams::addStringParameter(const char *name, const char *value)
{
    MIP_SolverParam<char*> *p;
    
    p = new (nothrow) MIP_SolverParam<char*>;
    if(!p)
		return MIP_MEMORY_ERROR;
    p->setParam(name, value);
    
    strParams.push_back(p);
	
    return 0;
}

void MIP_GeneralSolverParams::desallocate()
{
	list< MIP_SolverParam<int>* >::iterator itInt;
    list< MIP_SolverParam<double>* >::iterator itDbl;
    list< MIP_SolverParam<char*>* >::iterator itStr;
	
	for( itInt = intParams.begin(); itInt != intParams.end(); itInt++ )
		delete (*itInt);
	
	for( itDbl = dblParams.begin(); itDbl != dblParams.end(); itDbl++ )
	    delete (*itDbl);
	
	for( itStr = strParams.begin(); itStr != strParams.end(); itStr++ )
		delete (*itStr);
	
	
    intParams.clear();
    dblParams.clear();
    strParams.clear();
}

MIP_GeneralSolverParams::~MIP_GeneralSolverParams()
{
    desallocate();
}