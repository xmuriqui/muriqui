

#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"

#if OPT_HAVE_IQUAD
#include "iquad.hpp"
using namespace iquad;
#endif


using namespace optsolvers;





OPT_Iquad::OPT_Iquad()
{
    initialize();
}


OPT_Iquad::~OPT_Iquad()
{
    deallocateSolverEnv();
}



void OPT_Iquad::deallocateSolverEnv()
{
#if OPT_HAVE_IQUAD
    OPT_secDelete(bb);
#endif
    OPT_MyNLPSolver::deallocateSolverEnv();
}


bool OPT_Iquad::getMinusLambdaOnLagran()
{
    return false;
}


OPT_LISTSOLVERS OPT_Iquad::getSolverCode()
{
    return optsolvers::OPT_IQUAD;
}


void OPT_Iquad::initialize()
{
    OPT_MyNLPSolver::initialize();
    
    #if OPT_HAVE_IQUAD
        bb = NULL;
    #endif
}


int OPT_Iquad::getNumberOfIterations(long unsigned int& niter)
#if OPT_HAVE_IQUAD
{
    niter = bb->out_number_of_iterations;
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Iquad::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
#if OPT_HAVE_IQUAD
{
    
    /*const int r = OPT_MyNLPSolver::initSolverEnv(maxConstrs, maxVars, maxQuadNz);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return r;
    } */
    
    
    __desallocateSolverEnv();
    
    bb = new (std::nothrow) IQD_BranchAndBound;
    
    if( !bb )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        
        return OPT_MEMORY_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Iquad::setObjCutLowerBound(const double objLBound)
#if OPT_HAVE_IQUAD
{
    bb->in_lower_bound = objLBound;
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Iquad::setObjCutUpperBound(const double objUBound)
#if OPT_HAVE_IQUAD
{
    bb->in_upper_bound = objUBound;
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Iquad::setMaxCPUTime(const double time)
#if OPT_HAVE_IQUAD
{
    bb->in_max_cpu_time = time;
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Iquad::setMaxTime(const double time)
#if OPT_HAVE_IQUAD
{
    bb->in_max_time = time;
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Iquad::setNumberOfThreads(const int nthreads)
#if OPT_HAVE_IQUAD
{
    bb->in_number_of_threads = nthreads;
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Iquad::setOutputLevel( const int level )
#if OPT_HAVE_IQUAD
{
    bb->in_print_level = level;
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Iquad::setRelativeDualTol( const double tol )
#if OPT_HAVE_IQUAD
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Iquad::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_IQUAD
{
    bb->in_relative_convergence_tol = tol;
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Iquad::setRelativePrimalTol( const double tol )
#if OPT_HAVE_IQUAD
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Iquad::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_IQUAD
{
    const int r = bb->setDoubleParameter(param, value);
    
    if( r != 0 )
    {
        printDblParamErrorMsg(r, param, value );
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Iquad::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_IQUAD
{
    const int r = bb->setIntegerParameter(param, value);
    
    if( r != 0 )
    {
        printIntParamErrorMsg(r, param, value );
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Iquad::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_IQUAD
{
    const int r = bb->setStringParameter(param, value);
    
    if( r != 0 )
    {
        printStrParamErrorMsg(r, param, value );
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Iquad::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_IQUAD
{
    return solveWParams(resetSol, storeSol, storeConstrs, storeDualSol, NULL, NULL);
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Iquad::solveWParams(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol, OPT_GeneralSolverParams *subSolverParams, OPT_GeneralSolverParams *sdpParams)
#if OPT_HAVE_IQUAD
{
    const int n = prob.n;
    
    
    if( resetSol )
        this->resetSol();
    else
    {
        retCode = OPT_UNDEFINED;
        origSolverRetCode = INT_MAX;
        feasSol = false;
    }
    
    
    
    origSolverRetCode = bb->run( prob, subSolverParams, sdpParams );

    
    switch( origSolverRetCode )
    {
        case IQD_OPTIMAL_SOLUTION:
            retCode = OPT_OPTIMAL_SOLUTION;
            break;
        
        case IQD_INFEASIBLE_PROBLEM:
            retCode = OPT_INFEASIBLE_PROBLEM;
            break;
            
        case IQD_MAX_TIME_STOP:
            retCode = OPT_MAX_TIME;
            break;
            
        case IQD_MAX_ITERATIONS_STOP:
            retCode = OPT_MAX_ITERATIONS;
            #if OPT_PRINT_MAX_ITER_WARNING
                if( numberOfWarningsByIterLimit < maxNumberOfWarningsByIterLimit )
                {
                    std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Iquad solving!\n";
                    numberOfWarningsByIterLimit++;
                    
                    if( numberOfWarningsByIterLimit == maxNumberOfWarningsByIterLimit )
                        std::cerr << OPT_PREPRINT "Warning: Maximum number of warnings by maximum iteration achieved! Stopping these warnings.\n";
                }
            #endif
            
            break;
            
        case IQD_UNBOUNDED_PROBLEM:
            retCode = OPT_UNBOUNDED_PROBLEM;
            break;
            
        case IQD_MEMORY_ERROR:
            retCode = OPT_MEMORY_ERROR;
            break;
            
        case IQD_BAD_DEFINITIONS:
            retCode = OPT_BAD_INPUT;
            break;
            
        case IQD_CALLBACK_FUNCTION_ERROR:
            retCode = OPT_CALLBACK_FUNCTION_ERROR;
            break;
            
        case IQD_LIBRARY_NOT_AVAILABLE:
            retCode = OPT_LIBRARY_NOT_AVAILABLE;
            break;
            
        case IQD_SDP_SOLVING_ERROR:
        case IQD_QCP_SOLVER_ERROR:
        case IQD_NLP_SOLVER_ERROR:
            retCode = OPT_SUBSOLVER_ERROR;
            break;
            
        default:
            retCode = OPT_UNDEFINED_ERROR;
    }
    
    
    feasSol = bb->out_feasible_sol;
    objValue = bb->out_best_obj;
    dualObjValue = bb->out_lower_bound;
    
    if(storeSol)
        OPT_copyArray(n, bb->out_best_sol, sol);
    
    
    if( storeConstrs && feasSol )
    {
        if( nmChg )
        {
            OPT_setAllArray(prob.m, auxCEval, true);
            nmChg = false;
        }
        
        prob.constraintsEval( threadNumber, true, auxCEval, bb->out_best_sol, constr );
    }
    
    if( prob.objFactor < 0 )
    {
        objValue = -objValue;
        dualObjValue = -dualObjValue;
    }
    
    
    return retCode;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



