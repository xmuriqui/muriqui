/* this is a copy of cbc.cpp In this version, we create and ser a solver object apart form cbcModel and we clone this solver sobject everiytime OPT_Cbc::solve() is called. So, several cones has to be done if we have to solve a sequence of MILP problems changing parameters
*/


#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <iostream>
#include <new>


#include "SPM_DynSparseMatrix.hpp"

#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"


#if OPT_HAVE_CBC_OR_OSI
    #include "OsiSolverInterface.hpp"
    #include "OsiClpSolverInterface.hpp"
#endif


#if OPT_HAVE_CBC
    #include "CbcModel.hpp"
    #include "CoinModel.hpp"
    
    #include "CbcHeuristicLocal.hpp"
    
    // Cuts
    #include "CglGomory.hpp"
    #include "CglProbing.hpp"
    #include "CglKnapsackCover.hpp"
    #include "CglRedSplit.hpp"
    #include "CglClique.hpp"
    #include "CglFlowCover.hpp"
    #include "CglMixedIntegerRounding2.hpp"
    
    // Heuristics
    #include "CbcHeuristic.hpp"
#endif





using namespace std;
using namespace optsolvers;





OPT_Cbc:: OPT_Cbc(): OPT_OpenSolverInterface()
{
    initialize();
}


OPT_Cbc::~OPT_Cbc()
{
    deallocateMemory();
    deallocateSolverEnv();
}


// __methods from Solver __
/*int OPT_Cbc::__addConstraints(const int nm)
{
    OsiSolverInterface* const  solver = model->solver();
    
    for( int i = 0; i < nm; i++)
    {
        solver->addRow(0, NULL, NULL, -OPT_INFINITY, OPT_INFINITY);
    }
    
    return 0;
}


int OPT_Cbc::__addVariables(const int nn, const bool initFree)
{
    OsiSolverInterface* const solver = model->solver();
    
    const double c = 0.0;
    const double lb = -OPT_INFINITY, ub = OPT_INFINITY;
    
    
    for( int i = 0; i < nn; i++ )
    {
        solver->addCol(0, NULL, NULL, lb, ub, c);
    }
    
    return 0;
} */



void OPT_Cbc::deallocateMemory()
{
#if OPT_HAVE_CBC
    OPT_secDelete(solver);
    OPT_secDelete(model);
#endif
    OPT_OpenSolverInterface::deallocateMemory();
}



void OPT_Cbc::deallocateSolverEnv()
{
    OPT_OpenSolverInterface::deallocateSolverEnv();
}


/*int OPT_Cbc::getNumberOfConstraints(int &m)
#if OPT_HAVE_CBC
{
    m = model->getNumRows();
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif */


int OPT_Cbc::getNumberOfIterations(long unsigned int& niter)
#if OPT_HAVE_CBC
{
    
    niter = model->getIterationCount();
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


/*int OPT_Cbc::getNumberOfVars(int &n)
#if OPT_HAVE_CBC
{
    n = model->getNumCols();
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif */




OPT_LISTSOLVERS OPT_Cbc::getSolverCode()
{
    return OPT_CBC;
}


void OPT_Cbc::initialize()
{
    model = NULL;
    OPT_OpenSolverInterface::initialize();
#if OPT_HAVE_CBC
#endif
}


int OPT_Cbc::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
#if OPT_HAVE_CBC
{
    desallocateSolverEnv();
    
    
    model = new (std::nothrow) CbcModel();
    solver = new (std::nothrow) OsiClpSolverInterface();
    OPT_IFMEMERRORRETURN(!model || !solver);
    
    
    setOutputLevel(0);
    
    

    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cbc::setObjCutLowerBound(const double objLBound)
#if OPT_HAVE_CBC
{
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cbc::setObjCutUpperBound(const double objUBound)
#if OPT_HAVE_CBC
{
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cbc::setMaxCPUTime(const double time)
#if OPT_HAVE_CBC
{
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cbc::setMaxTime(const double time)
#if OPT_HAVE_CBC
{
    bool r = model->setDblParam(CbcModel::CbcMaximumSeconds, time);
    if(!r)
    {
        printIntParamErrorMsg( !r, "CbcMaximumSeconds", time );
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cbc::setNumberOfThreads(const int nthreads)
#if OPT_HAVE_CBC
{
    
    model->setNumberThreads(nthreads);
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cbc::setOutputLevel( const int level )
#if OPT_HAVE_CBC
{
    int r;
    

    
    r = model->setPrintingMode( level );
    OPT_IFERRORRETURN(r, OPT_BAD_INPUT);
    
    //r = model->setIntParam( CbcModel::CbcPrinting, level);
    //OPT_IFERRORRETURN(r, OPT_BAD_INPUT);
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cbc::setRelativeDualTol( const double tol )
#if OPT_HAVE_CBC
{

    
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cbc::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_CBC
{

    
    bool r = model->setDblParam(CbcModel::CbcAllowableFractionGap, tol);
    if(!r)
    {
        printIntParamErrorMsg( !r, "CbcAllowableFractionGap", tol);
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cbc::setRelativePrimalTol( const double tol )
#if OPT_HAVE_CBC
{

    
    OPT_OPERATIONNOTIMPLEMENTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cbc::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_CBC
{
    bool r;
    int key = CbcModel::CbcLastDblParam;
    

    
    //reading the parameter value in the string (it is not the ideal, but by now, we just do it)
    sscanf(param, "%d", &key);
    
    r = model->setDblParam((CbcModel::CbcDblParam) key, value);
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


int OPT_Cbc::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_CBC
{
    bool r;
    int key = CbcModel::CbcLastIntParam;
    
    
    //reading the parameter value in the string (it is not the ideal, but by now, we just do it)
    sscanf(param, "%d", &key);
    
    r = model->setIntParam((CbcModel::CbcIntParam) key, value);
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


int OPT_Cbc::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_CBC
{
    OPT_OPERATIONNOTSUPORTEDRET(getSolverCode());
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


/*int OPT_Cbc::setVariableType( const int index, const OPT_VARTYPE varType )
#if OPT_HAVE_CBC
{
    OsiSolverInterface* const  solver = model->solver();
    
    
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
#endif */


int OPT_Cbc::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_CBC_OR_OSI
{
    int r;
    int n, m;
    double *lxcopy = auxValues, *uxcopy = auxValues2; //we must copy var bounds because method Cbc::branchandbound change them... :(  
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
    
    
    //std::cout << "numberCutGenerators(): " << model->numberCutGenerators() << "\n";
    
    {
        //unfortunatelly, that seems the unique way to make cbc works fine. It is necessary alway clone the solver because CbcModel class change the solver object... :(
        OsiSolverInterface *solverClone =  solver->clone();
        model->assignSolver( solverClone );
        
        model->resetModel();
    }
    
    
    
    
    #if 0
    /*generating the cut generators
    * I take this portion of code from allCuts.cpp in examples cbc folder
    */
    {
        CglProbing generator1;
        generator1.setUsingObjective(true);
        generator1.setMaxPass(1);
        generator1.setMaxPassRoot(5);
        // Number of unsatisfied variables to look at
        generator1.setMaxProbe(10);
        generator1.setMaxProbeRoot(1000);
        // How far to follow the consequences
        generator1.setMaxLook(50);
        generator1.setMaxLookRoot(500);
        // Only look at rows with fewer than this number of elements
        generator1.setMaxElements(200);
        generator1.setRowCuts(3);

        CglGomory generator2;
        // try larger limit
        generator2.setLimit(300);

        CglKnapsackCover generator3;

        CglRedSplit generator4;
        // try larger limit
        generator4.setLimit(200);

        CglClique generator5;
        generator5.setStarCliqueReport(false);
        generator5.setRowCliqueReport(false);
        
        CglMixedIntegerRounding2 mixedGen;
        CglFlowCover flowGen;
        
        
        // Add in generators
        // Experiment with -1 and -99 etc
        // This is just for one particular model
        model->addCutGenerator(&generator1,-1,"Probing");
        //model.addCutGenerator(&generator2,-1,"Gomory");
        model->addCutGenerator(&generator2,1,"Gomory");
        model->addCutGenerator(&generator3,-1,"Knapsack");
        // model.addCutGenerator(&generator4,-1,"RedSplit");
        //model.addCutGenerator(&generator5,-1,"Clique");
        model->addCutGenerator(&generator5,1,"Clique");
        model->addCutGenerator(&flowGen,-1,"FlowCover");
        model->addCutGenerator(&mixedGen,-1,"MixedIntegerRounding");
    }
    
    
    /*setting up heuristics
    * I take this portion of code from example5.cpp in examples cbc folder
    */ 
    {
        // Allow rounding heuristic

        CbcRounding heuristic1(*model);
        model->addHeuristic(&heuristic1);

        // And local search when new solution found

        CbcHeuristicLocal heuristic2(*model);
        model->addHeuristic(&heuristic2);
    }
    #endif
    
    
    
    model->branchAndBound();
    
    //std::cout << "depois de model->branchAndBound() - model->solver(): " << model->solver() << "\n";
    
    /*solver = model->solver();//model->branchAndBound can generate a new solver object, so we update and pray to a new solver object not be generated in another object
    
    //we have to restore the variable bounds
    solver->setColLower(lxcopy);
    solver->setColUpper(uxcopy); */
    
    if( model->isProvenOptimal() )
    {
        retCode = OPT_OPTIMAL_SOLUTION;
        feasSol = true;
    }
    else if( model->isProvenInfeasible() )
    {
        retCode = OPT_INFEASIBLE_PROBLEM;
    }
    else if( model->isContinuousUnbounded() )
    {
        retCode = OPT_UNBOUNDED_PROBLEM;
    }
    else if( model->isProvenDualInfeasible()  )
    {
        retCode = OPT_UNBOUNDED_PROBLEM;
    }
    else if( model->isNodeLimitReached() )
    {
        retCode = OPT_MAX_ITERATIONS;
    }
    else if( model->isSecondsLimitReached() )
    {
        retCode = OPT_MAX_TIME;
    }
    else if( model->isAbandoned() )
    {
        retCode = OPT_UNDEFINED_ERROR;
    }
    else if( model->isSolutionLimitReached() )
    {
        retCode = OPT_UNDEFINED_ERROR;
    }
    else
    {
        retCode = OPT_UNDEFINED_ERROR;
    }
    
    
    r = getNumberOfConstraints(m);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    r = getNumberOfVars(n);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    
    
    objValue = model->getObjValue() + objConstant;
    dualObjValue = model->getBestPossibleObjValue() + objConstant;
    
    sprimalSol = model->getColSolution();
    sdualSolV = model->getReducedCost();
    sconstr = model->getRowActivity();
    sdualSolC = model->getRowPrice();
    
    
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
            const double sinf = model->getInfinity();
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












