
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <climits>

#include <iostream>
#include <new>


#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"


#if OPT_HAVE_CBC_OR_OSI
    #include "CbcModel.hpp"
    #include "OsiSolverInterface.hpp"
    #include "CoinModel.hpp"
#endif



using namespace std;
using namespace newspm;
using namespace optsolvers;



OPT_OpenSolverInterface:: OPT_OpenSolverInterface(): OPT_LPSolver()
{
    initialize();
}


OPT_OpenSolverInterface::~OPT_OpenSolverInterface()
{
    deallocateMemory();
    deallocateSolverEnv();
}


// __methods from Solver __
int OPT_OpenSolverInterface::__addConstraints(const int nm)
#if OPT_HAVE_CBC_OR_OSI
{
    int r, oldm;
    const double lb = -OPT_INFINITY, ub = OPT_INFINITY;
    
    r = getNumberOfConstraints(oldm);
    OPT_IFERRORRETURN(r, r);
    
    
    for( int i = 0; i < nm; i++)
    {
        solver->addRow(0, NULL, NULL, lb, ub);
    }
    
    //adding constraint in our structure
    {
        int r = A.addRows(nm);
        OPT_IFERRORRETURN(r, OPT_MEMORY_ERROR);
        
        r = OPT_realloc(lc, oldm + nm);
        OPT_IFERRORRETURN(r, r);
        
        r = OPT_realloc(uc, oldm + nm);
        OPT_IFERRORRETURN(r, r);
        
        OPT_setAllArray(nm, &lc[oldm], -OPT_INFINITY);
        OPT_setAllArray(nm, &uc[oldm], OPT_INFINITY);
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::__addVariables(const int nn, const bool initFree)
#if OPT_HAVE_CBC_OR_OSI
{
    const double c = 0.0;
    const double lb = -OPT_INFINITY, ub = OPT_INFINITY;
    int r;
    
    
    for( int i = 0; i < nn; i++ )
    {
        solver->addCol(0, NULL, NULL, lb, ub, c);
    }
    
    r = A.addCols(nn);
    OPT_IFERRORRETURN(r, OPT_MEMORY_ERROR);
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


void OPT_OpenSolverInterface::deallocateMemory()
{
    #if OPT_HAVE_CBC_OR_OSI
        //OPT_secDelete(solver);
    #endif
    OPT_secFree(lc);
    OPT_secFree(uc);
    A.desallocateMemory();
    OPT_LPSolver::deallocateMemory();
}

void OPT_OpenSolverInterface::deallocateSolverEnv()
{
    A.desallocateMemory();
}


int OPT_OpenSolverInterface::getNumberOfConstraints(int &m)
#if OPT_HAVE_CBC_OR_OSI
{
    m = A.getNumberOfRows(); //solver->getNumRows();
    
    #if OPT_DEBUG_MODE
        assert( m == solver->getNumRows());
    #endif
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::getNumberOfIterations(long unsigned int& niter)
#if OPT_HAVE_CBC_OR_OSI
{
    niter = solver->getIterationCount();
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::getNumberOfVars(int &n)
#if OPT_HAVE_CBC_OR_OSI
{
    n =solver->getNumCols();
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




OPT_LISTSOLVERS OPT_OpenSolverInterface::getSolverCode()
{
    return OPT_CBC;
}


int OPT_OpenSolverInterface::getVariableType( const int index, OPT_VARTYPE &varType )
#if OPT_HAVE_CBC_OR_OSI
{
    #if OPT_DEBUG_MODE
    {
        int r = checkVariableIndex(index);
        OPT_IFERRORRETURN(r, r);
    }
    #endif
    
    varType = solver->isInteger(index) ? OPT_VT_INTEGER : OPT_VT_CONTINUOUS;
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



void OPT_OpenSolverInterface::initialize()
{
    OPT_LPSolver::initialize();
    lc = NULL;
    uc = NULL;
    indexOfFirstChangedConstraint = INT_MAX;
    objConstant = 0.0;
#if OPT_HAVE_CBC_OR_OSI
    solver = NULL;
#endif
}


/*int OPT_OpenSolverInterface::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
#if OPT_HAVE_CBC_OR_OSI
{
    
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif */



int OPT_OpenSolverInterface::removeConstraints(const int ninds, const int* indices )
#if OPT_HAVE_CBC_OR_OSI
{
    int r;
    
    #if OPT_DEBUG_MODE
    {
        r = checkConstraintIndices(ninds, indices);
        OPT_IFERRORRETURN(r, r);
    }
    #endif
    
    solver->deleteRows(ninds, indices);
    
    r = A.removeRows(ninds, (SPM_dynIndex*) indices);
    OPT_IFERRORRETURN(r, OPT_BAD_INPUT);
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::removeVars(const int ninds, const int *indices )
#if OPT_HAVE_CBC_OR_OSI
{
    int r;
    
    #if OPT_DEBUG_MODE
    {
        r = checkVariableIndices(ninds, indices);
        OPT_IFERRORRETURN(r, r);
    }
    #endif
    
    solver->deleteCols(ninds, indices);
    
    r = A.removeColumns(ninds, (SPM_dynIndex*) indices);
    OPT_IFERRORRETURN(r, OPT_BAD_INPUT);
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::setObjCutLowerBound(const double objLBound)
#if OPT_HAVE_CBC_OR_OSI
{
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::setObjCutUpperBound(const double objUBound)
#if OPT_HAVE_CBC_OR_OSI
{
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::setMaxCPUTime(const double time)
#if OPT_HAVE_CBC_OR_OSI
{
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::setMaxTime(const double time)
#if OPT_HAVE_CBC_OR_OSI
{
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


/*int OPT_OpenSolverInterface::setNumberOfThreads(const int nthreads)
#if OPT_HAVE_CBC_OR_OSI
{
    model.setNumberThreads(nthreads);
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif */


int OPT_OpenSolverInterface::setOutputLevel( const int level )
#if OPT_HAVE_CBC_OR_OSI
{
    solver->messageHandler()->setLogLevel(level); //turn-off messages in Clp solver
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::setRelativeDualTol( const double tol )
#if OPT_HAVE_CBC_OR_OSI
{
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_CBC_OR_OSI
{
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::setRelativePrimalTol( const double tol )
#if OPT_HAVE_CBC_OR_OSI
{
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_CBC_OR_OSI
{
    bool r;
    int key = OsiLastDblParam;
    
    //reading the parameter value in the string (it is not the ideal, but by now, we just do it)
    sscanf(param, "%d", &key);
    
    r = solver->setDblParam((OsiDblParam) key, value);
    if(!r)
    {
        printIntParamErrorMsg( !r, param, value );
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_CBC_OR_OSI
{
    bool r;
    int key = OsiLastIntParam;
    
    //reading the parameter value in the string (it is not the ideal, but by now, we just do it)
    sscanf(param, "%d", &key);
    
    r = solver->setIntParam((OsiIntParam) key, value);
    if(!r)
    {
        printIntParamErrorMsg( !r, param, value );
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_CBC_OR_OSI	
{
    bool r;
    int key = OsiLastStrParam;
    
    //reading the parameter value in the string (it is not the ideal, but by now, we just do it)
    sscanf(param, "%d", &key);
    
    r = solver->setStrParam((OsiStrParam) key, value);
    if(!r)
    {
        printStrParamErrorMsg( !r, param, value );
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::setVariableType( const int index, const OPT_VARTYPE varType )
#if OPT_HAVE_CBC_OR_OSI	
{
    #if OPT_DEBUG_MODE
    {
        int r = checkVariableIndex(index);
        OPT_IFERRORRETURN(r, r);
    }
    #endif
    
    
    if( varType == OPT_VT_INTEGER )
    {
        solver->setInteger(index);
    }
    else if( varType == OPT_VT_CONTINUOUS )
    {
        solver->setContinuous(index);
    }
    else
    {
        //you have to do code to treat a new type of variable
        assert(false);
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_OpenSolverInterface::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_CBC_OR_OSI
{
    int r;
    int n, m, nI;
    const double *sprimalSol, *sdualSolV;
    const double *sconstr, *sdualSolC;
    
    
    
    if(resetSol)
    {
        this->resetSol();
    }
    else
    {
        feasSol = false;
    }
    
    r = updateSolver();
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    r = getNumberOfIntVars(nI);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    r = getNumberOfConstraints(m);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    r = getNumberOfVars(n);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    
    if(nI == 0)
    {
        solver->initialSolve();
    }
    else if(nI > 0)
    {
        solver->branchAndBound();
    }
    else
    {
        assert(false);
    }
    
    
    if( solver->isAbandoned() )
    {
        retCode = OPT_UNDEFINED_ERROR;
    }
    else if( solver->isProvenOptimal() )
    {
        retCode = OPT_OPTIMAL_SOLUTION;
        feasSol = true;
    }
    else if( solver->isProvenPrimalInfeasible() )
    {
        retCode = OPT_INFEASIBLE_PROBLEM;
    }
    else if( solver->isProvenDualInfeasible() )
    {
        retCode = OPT_UNBOUNDED_PROBLEM;
    }
    else if( solver->isPrimalObjectiveLimitReached() )
    {
        retCode = OPT_OPTIMAL_SOLUTION; //I am not sure about this retcode in this case...
        feasSol = true;
    }
    else if( solver->isDualObjectiveLimitReached() )
    {
        retCode = OPT_OPTIMAL_SOLUTION; //I am not sure about this retcode in this case...
    }
    else if( solver->isIterationLimitReached() )
    {
        retCode = OPT_MAX_ITERATIONS;
    }
    
    objValue = solver->getObjValue() + objConstant;
    
    
    sprimalSol = solver->getColSolution();
    sdualSolV = solver->getReducedCost();
    sconstr = solver->getRowActivity();
    sdualSolC = solver->getRowPrice();
    
    
    if( sprimalSol )
        OPT_copyArray(n, sprimalSol, sol);
    
    if( sdualSolV )
        OPT_copyArray(n, sdualSolV, dualSolV);
    
    if( sdualSolC )
        OPT_copyArray(m, sdualSolC, dualSolC);
    
    if(sconstr)
    {
        OPT_copyArray(m, sconstr, constr);
        
        //checking if solution is feasible
        if(feasSol == false)
        {
            const double sinf = solver->getInfinity();
            const double absTol = 1e-6;
            const double relTol = 1e-6;
            
            feasSol = true;
            for(int i = 0; i < m; i++)
            {
                if( lc[i] > -sinf )
                {
                    const double lci = lc[i] - (OPT_abs( lc[i]*relTol ) +  absTol);
                    
                    if( !(constr[i] >= lci) ) //we use not due to NAN 
                    {
                        feasSol = false;
                        break;
                    }
                }
                
                if( uc[i] < sinf )
                {
                    const double uci = uc[i] + (OPT_abs( uc[i]*relTol ) + absTol);
                    
                    if( !(constr[i] <= uci) ) //we use not due to NAN
                    {
                        feasSol = false;
                        break;
                    }
                }
            }
        }
    }
    
    
    
    
    
    
termination:
    
    return retCode;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_OpenSolverInterface::generateModelFile(const char *fileName)
#if OPT_HAVE_CBC_OR_OSI
{
    int r = updateSolver();
    OPT_IFERRORRETURN(r, r);
    
    solver->writeLp(fileName);
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


// __ methods from OPT_LPSolver __


int OPT_OpenSolverInterface::getConstraintBounds( const int index, double &lb, double &ub ) 
#if OPT_HAVE_CBC_OR_OSI
{
    int r;
    
    r = checkConstraintIndex(index);
    OPT_IFERRORRETURN(r, r);
    
    const double *lc = solver->getRowLower();
    const double *uc = solver->getRowUpper();
    const double solverInf = solver->getInfinity();
    
    lb = lc[index];
    ub = uc[index];
    
    #if OPT_DEBUG_MODE
        assert(lc[index] == this->lc[index]);
        assert(uc[index] == this->uc[index]);
    #endif
    
    if( lb <= -solverInf )
        lb = -OPT_INFINITY;
    
    if( ub >= solverInf )
        ub = OPT_INFINITY;
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::getFullConstraintLinearPart(const int constrIndex, double *values) 
#if OPT_HAVE_CBC_OR_OSI
{
    int r = A.getFullRow(constrIndex, values, false, 1.0);
    OPT_IFERRORRETURN(r, OPT_BAD_INPUT);
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value) 
#if OPT_HAVE_CBC_OR_OSI
{
    int r = A.getCoefficient(constrIndex, varIndex, value);
    OPT_IFERRORRETURN(r, OPT_BAD_INPUT);
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values) 
#if OPT_HAVE_CBC_OR_OSI
{
    int r;
    SPM_dynIndex mynzs;
    
    r = A.getCoefficientsInARow(constrIndex, mynzs, (SPM_dynIndex*) cols, values);
    OPT_IFERRORRETURN(r, OPT_BAD_INPUT);
    
    nzs = mynzs;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs) 
#if OPT_HAVE_CBC_OR_OSI
{
    unsigned int unzs;
    
    int r = A.getNumberOfElementsInARow(constrIndex, unzs);
    OPT_IFERRORRETURN(r, OPT_BAD_INPUT);
    
    nzs = unzs;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::getObjConstant(double &ObjConstant) 
#if OPT_HAVE_CBC_OR_OSI
{
    ObjConstant = this->objConstant;
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::getObjLinearCoef( const int index, double &value ) 
#if OPT_HAVE_CBC_OR_OSI
{
    int r = checkVariableIndex(index);
    OPT_IFERRORRETURN(r, r);
    
    const double *c = solver->getObjCoefficients();
    
    value = c[index];
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::getNumberOfIntVars(int &nI) 
#if OPT_HAVE_CBC_OR_OSI
{
    nI = solver->getNumIntegers();
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::getObjSense(OPT_OPTSENSE &sense) 
#if OPT_HAVE_CBC_OR_OSI
{
    sense = solver->getObjSense() > 0 ? OPT_MINIMIZE : OPT_MAXIMIZE;
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::getVariableBounds(const int index, double &lb, double &ub)  
#if OPT_HAVE_CBC_OR_OSI
{
    int r = checkVariableIndex(index);
    OPT_IFERRORRETURN(r, r);
    
    lb = solver->getColLower()[index];
    ub = solver->getColUpper()[index];
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




//warning: this method replaces all coefficients in a constraint.
int OPT_OpenSolverInterface::resetConstraintLinearPart( const int constrIndex, const int nzs, const int *cols, const double *values )  
#if OPT_HAVE_CBC_OR_OSI
{
    int r;
    
    r = checkVariableIndices(nzs, cols);
    OPT_IFERRORRETURN(r, r);
    
    r = checkConstraintIndex(constrIndex);
    OPT_IFERRORRETURN(r, r);
    
    r = A.setCoefficientsInARow(constrIndex, nzs, (const SPM_dynIndex*) cols, values);
    OPT_IFERRORRETURN(r, OPT_MEMORY_ERROR);
    
    if(constrIndex < indexOfFirstChangedConstraint)
        indexOfFirstChangedConstraint = constrIndex;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




//warning: this method can create a new auxiliary variable because some solvers like cplex and gurobi does not acept dual bounded constraints directly...
int OPT_OpenSolverInterface::setConstraintBounds( const int index, const double lb, const double ub ) 
#if OPT_HAVE_CBC_OR_OSI
{
    int r;
    const double solverInf = solver->getInfinity();
    
    r = checkConstraintIndex(index);
    OPT_IFERRORRETURN(r, r);
    
    const double mylb = lb <= -OPT_INFINITY ? -solverInf : lb;
    const double myub = ub >=  OPT_INFINITY ? solverInf : ub;
    
    solver->setRowBounds(index, mylb, myub);
    
    lc[index] = mylb;
    uc[index] = myub;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_OpenSolverInterface::setConstraintsLinearCoefs( const int nzs, const int *rows, const int *cols, const double *values ) 
#if OPT_HAVE_CBC_OR_OSI
{
    int r, firstRowChanged = INT_MAX;
    //const unsigned int oldnz = A.getNumberOfElements();
    
    
    if(nzs == 0)
        return 0;
    
    r = A.setCoefficients(nzs, (const SPM_dynIndex*) rows, (const SPM_dynIndex*) cols, values);
    OPT_IFERRORRETURN(r, OPT_UNDEFINED_ERROR);
    
    for(int k = 0; k < nzs; k++)
    {
        if( rows[k] < firstRowChanged )
            firstRowChanged = rows[k];
    }
    
    #if 0
    else if(nzs > 0)
    {
        //so, we have to merge the structures and values
        int m, n;
        //bool *usedRows = (bool*) auxIndex;
        double *fullRow = auxValues;
        std::set<int> usedRows_; //std::set is an ordered set
        int lastRow = rows[0] + 1; //just to initialize lastRow with a value different of rows[0]
        
        r = getNumberOfConstraints(m);
        OPT_IFERRORRETURN(r, r);
        
        r = getNumberOfVars(n);
        OPT_IFERRORRETURN(r, r);
        
        
        //here, we are sure nzs is greater than zero due to the first if in this function
        for(int k = 0; k < nzs; k++)
        {
            const int row = rows[k];
            
            if( row != lastRow) //to try avoid reinsert several times the same row index in usedRows. At least if the elemente are ordered b rows, it will be more efficient
            {
                usedRows_.insert( row );
                lastRow = row;
            }
        }
        
        
        for(std::set<int>::iterator it = usedRows_.begin(); it != usedRows_.end() ; ++it )
        {
            const int row = *it;
            
            r = A.getFullRow(row, fullRow, false, false, 1.0);
            OPT_IFERRORRETURN(r, OPT_UNDEFINED_ERROR);
            
            for(int k = 0; k < nzs; k++)
            {
                if( rows[k] == row )
                    fullRow[ cols[k] ] = values[k];
            }
            
            r = A.setRowStructureAndValues(row, fullRow, n);
            OPT_IFERRORRETURN(r, OPT_MEMORY_ERROR);
        }
        
        firstRowChanged = *(usedRows_.begin()); //the set is ordered. So, the lower row index is the first in the set
        
    }
    #endif
    
    
    if( firstRowChanged < indexOfFirstChangedConstraint )
        indexOfFirstChangedConstraint = firstRowChanged;
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::setConstraintLinearCoefs( const int constrIndex, const int nzs, const int *cols, const double *values)  
#if OPT_HAVE_CBC_OR_OSI
{
    int r;
    //unsigned int oldnzs;
    
    
    //if(oldnzs == 0)
    {
        r = A.setCoefficientsInARow(constrIndex, nzs, (const SPM_dynIndex*) cols, values);
        OPT_IFERRORRETURN(r, OPT_UNDEFINED_ERROR);
    }
    #if 0
    else
    { //so, we have to accumulate old coefficients with the new ones
        int n;
        double *fullRow = auxValues;
        
        r = getNumberOfVars(n);
        OPT_IFERRORRETURN(r, r);
        
        r = A.getFullRow(constrIndex, fullRow, false, false, 1.0);
        OPT_IFERRORRETURN(r, OPT_MEMORY_ERROR);
        
        for(int i = 0; i < nzs; i++)
            fullRow[ cols[i] ] = values[i];
        
        r = A.setRowStructureAndValues(constrIndex, fullRow, n);
        OPT_IFERRORRETURN(r, OPT_MEMORY_ERROR);
    }
    #endif
    
    
    if( constrIndex < indexOfFirstChangedConstraint )
        indexOfFirstChangedConstraint = constrIndex;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value)  
#if OPT_HAVE_CBC_OR_OSI
{
    return setConstraintLinearCoefs(constrIndex, 1, &varIndex, &value);
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_OpenSolverInterface::setObjLinearCoef( const int index, const double value )  
#if OPT_HAVE_CBC_OR_OSI
{
    int r = checkVariableIndex(index);
    OPT_IFERRORRETURN(r, r);
    
    solver->setObjCoeff(index, value);
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




void OPT_OpenSolverInterface::setObjConstant(const double value)  
{
    objConstant= value;
}



int OPT_OpenSolverInterface::setObjSense( const OPT_OPTSENSE sense )  
#if OPT_HAVE_CBC_OR_OSI
{
    solver->setObjSense( sense == OPT_MAXIMIZE ? -1 : 1 );
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_OpenSolverInterface::setVariableBounds( const int index, const double lb, const double ub )  
#if OPT_HAVE_CBC_OR_OSI
{
    const double solverInf = solver->getInfinity();
    
    int r;
    double mylb, myub;
    
    r = checkVariableIndex(index);
    OPT_IFERRORRETURN(r, r);
    
    
    mylb = lb <= -OPT_INFINITY ? -solverInf : lb;
    myub = ub >= OPT_INFINITY ? solverInf : ub;
    
    solver->setColBounds(index, mylb, myub);
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_OpenSolverInterface::updateSolver()
#if OPT_HAVE_CBC_OR_OSI
{
    int r, m;
    SPM_dynIndex *cols = (SPM_dynIndex *) auxIndex;
    SPM_dynValue *values = auxValues;
    
    r = getNumberOfConstraints(m);
    OPT_IFERRORRETURN(r, r);
    
    if( indexOfFirstChangedConstraint < m )
    { //so, user updated matrix of constraints. By now, our unique choice is delete constraints and reinsert by our structures
        
        int *rowsToReadd = auxIndex;
        const int myIndexOfFirstChangedConstraint = indexOfFirstChangedConstraint; //just to pute the value in a ocal variable. Maybe, in this way compiler will be sure to use vectorization with #prama ivdeep since it can be sure this variable will no change.
        const int nrowsToReadd = m-myIndexOfFirstChangedConstraint;
        
        
        #pragma ivdep
        #pragma GCC ivdep
        for( int i = myIndexOfFirstChangedConstraint; i < m; i++ )
        {
            rowsToReadd[ i - myIndexOfFirstChangedConstraint ] = i;
        }
        
        solver->deleteRows( nrowsToReadd, rowsToReadd );
        
        //now, we re-add the rows. By miss of information, we will rows, onde by one because I am not sure about parameters of addRows method
        for(int i = myIndexOfFirstChangedConstraint; i < m; i++)
        {
            SPM_dynIndex rnzs;
            
            r = A.getCoefficientsInARow(i, rnzs, cols, values);
            OPT_IFERRORRETURN(r, OPT_UNDEFINED);
            
            solver->addRow( rnzs, (int *) cols, values, lc[i], uc[i] );
        }
    }
    
    
    indexOfFirstChangedConstraint = INT_MAX;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif









