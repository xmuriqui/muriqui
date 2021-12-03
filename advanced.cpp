
#include <cstdlib>

#include <new>

#include "muriqui.hpp"
#include "MRQ_advanced.hpp"


using namespace std;
using namespace muriqui;



int MRQ_SolutionStorer::addSolution(const unsigned int n, const double objValue, const double *solution)
{
    double *sol;
    
    MRQ_malloc(sol, n);
    
    MRQ_IFMEMERRORRETURN(!sol);
    
    try
    {
        sols.push_back( std::pair<double, double*>(objValue, sol) );
    }
    catch(bad_alloc& ba)
    {
        return MRQ_MEMORY_ERROR;
    }
    
    MRQ_copyArray(n, solution, sol);
    
    return 0;
}


void MRQ_SolutionStorer::deallocate()
{
    for(auto p: sols)
    {
        free(p.second);
    }
    
    sols.clear();
}


MRQ_SolutionStorer::~MRQ_SolutionStorer()
{
    deallocate();
}



bool MRQ_UserCallbacks::tryUpdateBestSolution(const int threadNumber, const int n, double* solution, const double fsolution, const long unsigned int iter)
{
    return alg->tryUpdateBestSolution(threadNumber, n, solution, fsolution, iter, clockStart, timeStart, alg->in_store_history_solutions);
}


long unsigned int MRQ_UserCallbacks::BB_getOpenNodes(long unsigned int numberOfNodes, MRQ_NewBBNode* &nodes, bool sortNodesByBound)
{
    branchAndBound::BBL_Node *bblNodep;
    auto ret = bblCallBacks->getNodes(numberOfNodes, bblNodep, true);
    
    nodes = (MRQ_NewBBNode*) bblNodep;
    
    return ret;
}


double MRQ_UserCallbacks::getLowerBound() const
{
    return alg->zl;
}

double MRQ_UserCallbacks::getUpperBound() const
{
    return alg->zu;
}


int MRQ_UserCallbacks:: BB_getNumberOfOpenNodes(long unsigned int &nnodes) const
{
    if( bblCallBacks )
    {
        nnodes = bblCallBacks->getNumberOfOpenNodes();
        return 0;
    }
    else
    {
        nnodes = ULONG_MAX;
        return  MRQ_BAD_DEFINITIONS;
    }
}


int MRQ_UserCallbacks::getBestSolutionCopy(const int n, double *solution, double &fsolution)
{
    return alg->getBestSolutionCopy(n, solution, fsolution);
}



int MRQ_UserCallbacks::addGlobalQuadraticCut(const int nza, const int *acols, const double *avalues, const int nzQ, const int *Qrows, const int* Qcols, const double *Qvalues, const double lb, const double ub)
{
    //if(cutGen)
        //return cutGen->setQuadraticCut(nza, acols, avalues, nzQ, Qrows, Qcols, Qvalues, lb, ub);
    //else
        return cutGenerator->setQuadraticCut(nza, acols, avalues, nzQ, Qrows, Qcols, Qvalues, lb, ub);
}


void MRQ_UserCallbacks::BB_printOpenNodesList()
{
    if( bblCallBacks )
    {
        bblCallBacks->printOpenNodesList();
    }
    //else if(alg->out_algorithm == MRQ_BB_ALG)
        //return ((MRQ_BranchAndBound *) alg)->nodes->print();
}


MRQ_Preprocessor* MRQ_UserCallbacks::BB_getPreprocessorPointer( unsigned int threadNumber)
{
    if( bblCallBacks )
        return bblCallBacks->getPreprocessorPointer(threadNumber);
    else
        return nullptr;
}


minlpproblem:: MIP_ConstraintsByColumnsStorager * MRQ_UserCallbacks::BB_getConstraintsByColumnsStoragerPointer()
{
    if( bblCallBacks )
        return bblCallBacks->ccstorager;
    else
        return nullptr;
}



int muriqui::MRQ_insideRun( MRQ_Algorithm *alg, MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, const int thnumber,  const double insideSolverMaxTime, double *nlx, double *nux )
{
    return alg->insideRun( prob, milpSolverParams, nlpSolverParams, thnumber, insideSolverMaxTime, nlx, nux );
}




























