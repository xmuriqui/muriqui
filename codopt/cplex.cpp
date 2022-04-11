
#include <cstdlib>
#include <cassert>
#include <climits>

#include <iostream>
#include <list>

#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"




using namespace std;
using namespace optsolvers;



OPT_Cplex::OPT_Cplex():OPT_QPSolver()
{
    initialize();
}





OPT_Cplex::~OPT_Cplex()
{
    //desallocateMemory();
    deallocateSolverEnv();
}





// __methods from Solver __


int  OPT_Cplex::cloneFrom( OPT_LPSolver &other)
#if OPT_HAVE_CPLEX
{
    int status;
    CPXLPptr aux;
    OPT_Cplex &cother = (OPT_Cplex&) other;
    
    if( other.getSolverCode() != OPT_CPLEX )
        return OPT_LPSolver::copyLPPartFrom( other);
    
    if(!env)
    {
        env = CPXopenCPLEX(&status);
        if( !env )
        {
            char errmsg[CPXMESSAGEBUFSIZE];
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(status);
            #endif
            
            CPXgeterrorstring (env, status, errmsg);
            std::cerr << OPT_PREPRINT "Cplex error msg: " << errmsg << "\n";
            
            return OPT_SOLVER_ERROR;
        }
    }
    
    
    aux = CPXcloneprob (cother.env, cother.prob, &status);
    if(!aux)
    {
        char errmsg[CPXMESSAGEBUFSIZE];
        OPT_PRINTERRORMSGP("Error to clone a CPLEX problem: ", status);
        
        CPXgeterrorstring (env, status, errmsg);
        std::cerr << OPT_PREPRINT "Cplex error msg: " << errmsg << "\n";
        
        return OPT_SOLVER_ERROR;
    }
    
    if(prob)
        CPXfreeprob(env, &prob);
    
    prob = aux;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


void OPT_Cplex::deallocateSolverEnv()
{
#if OPT_HAVE_CPLEX
    if (prob)
    {
        CPXfreeprob(env, &prob);
        prob = NULL;
    }

    if (env)
    {
        CPXcloseCPLEX(&env);
        env = NULL;
    }
#endif
    deallocateMemory();
}


int OPT_Cplex::getNumberOfIterations(long unsigned int& niter)
#if OPT_HAVE_CPLEX
{
    if( CPXgetnumint(env, prob) > 0 )
    {
        niter = CPXgetnodecnt(env, prob);
    }
    else
    {
        niter = OPT_max( CPXgetitcnt (env, prob), CPXgetbaritcnt (env, prob) );
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif

OPT_LISTSOLVERS OPT_Cplex::getSolverCode()
{
    return OPT_CPLEX;
}


int OPT_Cplex::getVariableType( const int index, OPT_VARTYPE &varType )
#if OPT_HAVE_CPLEX
{
    
    if( CPXgetnumint(env, prob) > 0 || CPXgetnumbin(env, prob) > 0 )  //cplex return an error if we call CPXgetctype in a non integer problem. So, we have to perform this test
    {
        char vt;
        
        const int r = CPXgetctype(env, prob, &vt, index, index);
        OPT_IFERRORRETURN(r, OPT_SOLVER_ERROR);
        
        varType = vt == CPX_INTEGER || vt == CPX_BINARY ? OPT_VT_INTEGER : OPT_VT_CONTINUOUS; 
    }
    else
    {
        varType = OPT_VT_CONTINUOUS;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


void OPT_Cplex::initialize()
{
    OPT_QPSolver::initialize();
    
#if OPT_HAVE_CPLEX
    origSolverRetCode = INT_MAX;
    env = NULL;
    prob = NULL;
#endif
    
    objConstant = 0.0;
}



int OPT_Cplex::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
#if OPT_HAVE_CPLEX
{
    int val = 0;
    
    env = CPXopenCPLEX(&val);
    
    if( !env )
    {
        char errmsg[CPXMESSAGEBUFSIZE];
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(val);
        #endif
        
        CPXgeterrorstring (env, val, errmsg);
        std::cerr << OPT_PREPRINT "Cplex error msg: " << errmsg << "\n";
        
        return OPT_SOLVER_ERROR;
    }
    
    
    prob = CPXcreateprob( env, &val, "qpCplexOptSolvers" );
    if( !prob )
    {
        char errmsg[CPXMESSAGEBUFSIZE];
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(val);
        #endif
        
        CPXgeterrorstring (env, val, errmsg);
        std::cerr << OPT_PREPRINT "Cplex error msg: " << errmsg << "\n";
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::setObjCutLowerBound(const double objLBound)
#if OPT_HAVE_CPLEX
{
    //that is a MIP parameter, but it is better than nothing...
    int r = CPXsetdblparam( env, CPX_PARAM_CUTLO, objLBound );
    
    //try it some day:
    //r = CPXsetdblparam( env, CPX_PARAM_OBJLLIM, objLBound );
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::setObjCutUpperBound(const double objUBound)
#if OPT_HAVE_CPLEX
{
    //that is a MIP parameter, but it is better than nothing...
    int r = CPXsetdblparam( env, CPX_PARAM_CUTUP, objUBound );
    
    //try it some day:
    //r = CPXsetdblparam( env, CPX_PARAM_OBJULIM, objLBound );
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::setMaxCPUTime(const double time)
#if OPT_HAVE_CPLEX
{
    int r;
    
    //set cplex to use cpu time instead of real time...
    r = CPXsetintparam(env, CPX_PARAM_CLOCKTYPE, 1);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            OPT_PRINTERRORNUMBER(r);
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    r = CPXsetdblparam(env, CPX_PARAM_TILIM, time);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            OPT_PRINTERRORNUMBER(r);
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::setMaxTime(const double time)
#if OPT_HAVE_CPLEX
{
    int r;
    
    //set cplex to use real time
    r = CPXsetintparam(env, CPX_PARAM_CLOCKTYPE, 2);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            OPT_PRINTERRORNUMBER(r);
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    r = CPXsetdblparam(env, CPX_PARAM_TILIM, time);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            OPT_PRINTERRORNUMBER(r);
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::setNumberOfThreads(const int nthreads)
#if OPT_HAVE_CPLEX
{
    int r = CPXsetintparam( env, CPX_PARAM_THREADS, nthreads );
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            OPT_PRINTERRORNUMBER(r);
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::setOutputLevel( const int level )
#if OPT_HAVE_CPLEX
{
    
    int r = CPXsetintparam( env, CPX_PARAM_SCRIND, level > 0 ? CPX_ON : CPX_OFF );
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            OPT_PRINTERRORNUMBER(r);
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::setRelativeDualTol( const double tol )
#if OPT_HAVE_CPLEX
{
    return OPT_OPERATION_NOT_IMPLEMENTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cplex::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_CPLEX
{
#if OPT_HAVE_CPLEX
    int r = CPXsetdblparam(env, CPX_PARAM_BAREPCOMP, tol );
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
    }
    
    r = CPXsetdblparam(env, CPX_PARAM_BARQCPEPCOMP, tol );
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
    }
    
    r = CPXsetdblparam(env, CPX_PARAM_EPGAP, tol );
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
    }
    
    return 0;
    
#endif
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cplex::setRelativePrimalTol( const double tol )
#if OPT_HAVE_CPLEX
{
#if OPT_HAVE_CPLEX
    return OPT_OPERATION_NOT_IMPLEMENTED;
#endif
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_CPLEX
{
    int aux, r;
    
    r = CPXgetparamnum( env, param, &aux );
    if(r == 0)
    {
        r = CPXsetdblparam( env, aux, value );
    }
    
    if(r != 0)
    {
        printDblParamErrorMsg( r, param, value );
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_CPLEX
{
    int aux, r;
    
    r = CPXgetparamnum( env, param, &aux );
    if(r == 0)
    {
        r = CPXsetintparam( env, aux, value );
    }
    
    
    if(r != 0)
    {
        printIntParamErrorMsg( r, param, value );
        
        return OPT_BAD_INPUT;
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_CPLEX
{
    int aux, r;
    
    
    r = CPXgetparamnum( env, param, &aux );
    if(r == 0)
    {
        r = CPXsetstrparam( env, aux, value );
    }
    
    
    if(r != 0)
    {
        printStrParamErrorMsg( r, param, value );
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



/*int OPT_Cplex::setParameters( OPT_GeneralSolverParams &params )
#if OPT_HAVE_CPLEX
{
    int r, aux, code = 0;
    char solver[] = "cplex";
    
    list< OPT_SolverParam<int>* >::iterator itInt;
    list< OPT_SolverParam<double>* >::iterator itDbl;
    list< OPT_SolverParam<char*>* >::iterator itStr;
    
    
    for( itInt = params.intParams.begin(); itInt != params.intParams.end(); itInt++ )
    {
        r = CPXgetparamnum(env, (*itInt)->name, &aux);
        if(r != 0)
        {
            printIntParamErrorMsg( r, solver, (*itInt)->name, (*itInt)->value );
            
            code = OPT_BAD_INPUT;
        }
        else
        {
            r = CPXsetintparam(env, aux, (*itInt)->value);
            if(r != 0)
            {
                printIntParamErrorMsg( r, solver, (*itInt)->name, (*itInt)->value );
                //getchar();
                
                code = OPT_BAD_INPUT;
            }
        }
    }
    
    
    
    for( itDbl = params.dblParams.begin(); itDbl != params.dblParams.end(); itDbl++ )
    {
        r = CPXgetparamnum(env, (*itDbl)->name, &aux);
        if(r != 0)
        {
            printDblParamErrorMsg( r, solver, (*itDbl)->name, (*itDbl)->value );
            
            code = OPT_BAD_INPUT;
        }
        else
        {
            r = CPXsetdblparam(env, aux, else
    {
    }(*itDbl)->value);
            if(r != 0)
            {
                printDblParamErrorMsg( r, solver, (*itDbl)->name, (*itDbl)->value );
                //getchar();
                
                code = OPT_BAD_INPUT;
            }
        }
    }
    
    
    
    for( itStr = params.strParams.begin(); itStr != params.strParams.end(); itStr++ )
    {
        r = CPXgetparamnum(env, (*itStr)->name, &aux);
        if(r != 0)
        {
            printStrParamErrorMsg( r, solver, (*itStr)->name, (*itStr)->value );
            
            code = OPT_BAD_INPUT;
        }
        else
        {
            r = CPXsetstrparam(env, aux, (*itStr)->value);
            if(r != 0)
            {
                printStrParamErrorMsg( r, solver, (*itStr)->name, (*itStr)->value );
                
                code = OPT_BAD_INPUT;
            }
        }
    }



    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif
*/



int OPT_Cplex::setVariableType( const int index, const OPT_VARTYPE varType )
#if OPT_HAVE_CPLEX
{
    const char type = varType == OPT_VT_INTEGER ? 'I' : 'C';
    
    int r = CPXchgctype( env, prob, 1, &index, &type );
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    //cplex is ridiculous. Even if you set all variables like continuous. It changes the problem type to MILP pnly by the fact of you called CPXchgctype
    
    if( CPXgetnumint(env, prob) == 0 )
    {
        const int cptype = CPXgetprobtype(env, prob);
        int newptype;
        
        switch( cptype )
        {
            case -1:
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                    
                return OPT_SOLVER_ERROR;
                
            case CPXPROB_MILP:
            case CPXPROB_FIXEDMILP:
            case CPXPROB_LP: //that case is crazy, but...
                newptype = CPXPROB_LP;
                break;
                
            case CPXPROB_MIQP:
            case CPXPROB_FIXEDMIQP:
            case CPXPROB_QP: //that case is crazy, but...
                newptype = CPXPROB_QP;
                break;
                
            default:
                newptype = CPXPROB_QCP;
        }
        
        
        r = CPXchgprobtype(env, prob, newptype);
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            return OPT_SOLVER_ERROR;
        }
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::getSolution(const bool integer, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_CPLEX
{
    int r, status;
    
    double *psol, *pconstr;
    double *pdualc, *pdualv; 
    
    
    if( integer || !storeDualSol )
    {
        pdualc = NULL;
        pdualv = NULL;
    }
    else
    {
        pdualc = dualSolC;
        pdualv = dualSolV;
    }
    
    
    psol = storeSol ? sol : NULL;
    pconstr = storeConstrs ? constr : NULL;
    
    
    feasSol = false;
    
    r = CPXsolution( env, prob, &status, &objValue, psol, pdualc, pconstr, pdualv );
    
    if( r == 0)
    {
        int pfeasind_p;
        
        r = CPXsolninfo(env, prob, NULL, NULL, &pfeasind_p, NULL);
        
        if(pfeasind_p)
            feasSol = true;
    }
    else
    {
        /* #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif */
        
        if(r != CPXERR_NO_SOLN)
        {
            cerr << OPT_PREPRINT << "Error " << r << " to obtain solution at CPLEX" << OPT_GETFILELINE << endl;
            
            return OPT_SOLVER_ERROR;
        }
        
    }
    
    
    if( storeConstrs )
    {
        int m;
        getNumberOfConstraints(m);
        
        //constr has the slack of constraints...
        r = CPXgetrhs( env, prob, auxValues, 0, m-1 );
        
        for( int i = 0; i < m; i++ )
            constr[i] = auxValues[i] - constr[i];
    }
    
    
    objValue += objConstant;
    dualObjValue = objValue;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_CPLEX
{
    int r, v, status;
    int m = 0;
    int paramSolutionTArgetCode;
    
    getNumberOfConstraints(m);
    
    if(resetSol)
    {
        origSolverRetCode = INT_MAX;
        this->resetSol();
    }
    else
    {
        feasSol = false;
    }
    
    //checking if we are using cplex to nonconvex target...
    v = 0;
    
    #ifdef CPX_PARAM_SOLUTIONTARGET
        paramSolutionTArgetCode = CPX_PARAM_SOLUTIONTARGET;
    #else
        paramSolutionTArgetCode = CPXPARAM_OptimalityTarget;
    #endif
    
    CPXgetintparam(env, paramSolutionTArgetCode, &v);
    
    
    if( v == 3 || CPXgetnumint(env, prob) > 0 )
    {
        int optstatus = CPXmipopt(env, prob);
        
        if(optstatus != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(optstatus);
                
                r = CPXgetsubstat(env, prob);
                std::cerr << OPT_PREPRINT " cplex opt substatus: " << r << "\n";
            #endif
                
            retCode = OPT_SOLVER_ERROR;
            //goto termination; we do not go to termination and decide about retCode using CPXgetstat
        }
        
        
        status = CPXgetstat(env, prob);
        
        r = getSolution(true, storeSol, storeConstrs, storeDualSol);
        
        /*get dual bound on integer problems. If the problem is infeasible,
         */
        {
            double myDualObjV = -INFINITY;
            
            int r2 = CPXgetbestobjval( env, prob, &myDualObjV );
            if(r2)
                OPT_PRINTERRORMSGP("Error at CPXgetbestobjval: ", r2);
            
            if( -OPT_INFINITY < myDualObjV && myDualObjV < OPT_INFINITY )
                dualObjValue = myDualObjV + objConstant;
        }
        
        origSolverRetCode = status;
        
        switch( status )
        {
        case CPXMIP_OPTIMAL:
        case CPXMIP_FEASIBLE_RELAXED_INF:
        case CPXMIP_OPTIMAL_TOL:
        case CPXMIP_OPTIMAL_POPULATED:
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                retCode = r;
                goto termination;
            }
            
            
            retCode = OPT_OPTIMAL_SOLUTION;
            
            break;
            
        case CPXMIP_INFEASIBLE:
        case CPXMIP_INForUNBD:
        case CPXMIP_ABORT_INFEAS:
            retCode = OPT_INFEASIBLE_PROBLEM;
            break;

        case CPXMIP_FAIL_INFEAS_NO_TREE:
            retCode = OPT_MEMORY_ERROR;
            break;

        case CPXMIP_TIME_LIM_FEAS:
        case CPXMIP_TIME_LIM_INFEAS:
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
                #endif
                
                retCode = r;
                goto termination;
            }
            
            retCode = OPT_MAX_TIME;
            break;
            
            
        case CPXMIP_UNBOUNDED:
            retCode = OPT_UNBOUNDED_PROBLEM;
            break;
            
        case CPXMIP_NODE_LIM_FEAS:
        case CPXMIP_NODE_LIM_INFEAS:
            
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                retCode = r;
                goto termination;
            }
            
            #if OPT_PRINT_MAX_ITER_WARNING
                if( numberOfWarningsByIterLimit < maxNumberOfWarningsByIterLimit )
                {
                    std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Cplex solving!\n";
                    numberOfWarningsByIterLimit++;
                    
                    if( numberOfWarningsByIterLimit == maxNumberOfWarningsByIterLimit )
                        std::cerr << OPT_PREPRINT "Warning: Maximum number of warnings by maximum iteration achieved! Stopping these warnings.\n";
                }
            #endif
            
            
            retCode = OPT_MAX_ITERATIONS;
            break;
            
            
        default:
            
            retCode = OPT_UNDEFINED_ERROR;
        }
        
        
    }
    else
    {
        const int cptype = CPXgetprobtype(env, prob);
        //std::cout << "cplex problem type: " << cptype << std::endl;
        
        if( cptype == CPXPROB_LP || cptype == CPXPROB_FIXEDMILP )
        {
            //I think it is faster call this function if the problem is linear. So, cplex can choose the best optimizer available for the problem... (we allow user change it by set a cplex parameter also...)
            r = CPXlpopt(env, prob);
        }
        else if( cptype == CPXPROB_QP || cptype == CPXPROB_FIXEDMIQP )
        {
            r = CPXqpopt(env, prob); //qpopt by default call baropt. However, user can change this cplex parameter. So, it is better call CPXqpopt here instead of CPXbaropt to attend a posible user choice...
        }
        else
        {
            r = CPXbaropt(env, prob);
        }
        
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            retCode = OPT_SOLVER_ERROR;
            goto termination;
        }
        
        
        status = CPXgetstat(env, prob);
        
        r = getSolution(false, storeSol, storeConstrs, storeDualSol);
        
        origSolverRetCode = status;
        
        switch(status)
        {
        case CPX_STAT_OPTIMAL:
            
            
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                retCode = r;
                goto termination;
            }
            
            retCode = OPT_OPTIMAL_SOLUTION;
            break;
            
            
        case CPX_STAT_INFEASIBLE:
            retCode = OPT_INFEASIBLE_PROBLEM;
            break;
            
            
        case CPX_STAT_INForUNBD:
        case CPX_STAT_UNBOUNDED:
            retCode = OPT_UNBOUNDED_PROBLEM;
            break;
            
            
        case CPX_STAT_ABORT_TIME_LIM:
            
            retCode = OPT_MAX_TIME;
            break;
            
        case CPX_STAT_ABORT_IT_LIM:
            
            #if OPT_PRINT_MAX_ITER_WARNING
                if( numberOfWarningsByIterLimit < maxNumberOfWarningsByIterLimit )
                {
                    std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Cplex solving!\n";
                    numberOfWarningsByIterLimit++;
                    
                    if( numberOfWarningsByIterLimit == maxNumberOfWarningsByIterLimit )
                        std::cerr << OPT_PREPRINT "Warning: Maximum number of warnings by maximum iteration achieved! Stopping these warnings.\n";
                }
            #endif
            
            
            retCode = OPT_MAX_ITERATIONS;
            break;
        
        default:
            
            retCode = OPT_UNDEFINED_ERROR;
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






// __ methods from LPSolver __



int OPT_Cplex::__addConstraints(const int nm)
#if OPT_HAVE_CPLEX
{
    //r = CPXaddrows( env, prob, 0, nm, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL );
    
    int r = CPXnewrows( env, prob, nm, NULL, NULL, NULL, NULL );
    OPT_IFERRORRETURN(r, OPT_SOLVER_ERROR);
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::__addVariables(const int nn, const bool initFree)
#if OPT_HAVE_CPLEX
{
    int r;
    double *lb;
    
    if( initFree )
    {
        lb = auxValues;
        
        for(int i = 0; i < nn; i++)
            lb[i] = -CPX_INFBOUND;
    }
    else
    {
        lb  = NULL;
    }
    
    
    r = CPXnewcols(env, prob, nn, NULL, lb, NULL, NULL, NULL );
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif





int OPT_Cplex::generateModelFile(const char* fileName)
#if OPT_HAVE_CPLEX
{
    int indNewVar = -1;
    int r, code;
    
    
    //if we have a constant, we have to add an artificial variable. To avoid confusion about indices (if new variables are added after), we delete the new variable after write the file...
    
    if( objConstant != 0.0 )
    {
        double obj = 1.0, lb, ub;
        
        
        lb = ub = objConstant;
        r = CPXnewcols(env, prob, 1, &obj, &lb, &ub, NULL, NULL);
        
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            code = OPT_SOLVER_ERROR;
            goto termination;
        }
        
        indNewVar = CPXgetnumcols(env, prob) -1;
    }
    
    
    
    
    r = CPXwriteprob(env, prob, fileName, "LP");
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        code = OPT_SOLVER_ERROR;
        goto termination;
    }
    
    
    
    code = 0;
    
termination:
    
    if( indNewVar >= 0 )
    {
        r = CPXdelcols(env, prob, indNewVar, indNewVar);
        
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            return OPT_SOLVER_ERROR;
        }
    }
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::getConstraintBounds( const int index, double &lb, double &ub )
#if OPT_HAVE_CPLEX
{
    int r;
    char sense;
    
    r = CPXgetsense(env, prob, &sense, index, index);
    OPT_IFERRORRETURN(r, OPT_SOLVER_ERROR);
    
    r = CPXgetrhs(env, prob, &lb, index, index);
    OPT_IFERRORRETURN(r, OPT_SOLVER_ERROR);
    
    
    if( sense == 'L' )
    {
        ub = lb;
        lb = -OPT_INFINITY;
    }
    else if( sense == 'E' )
    {
        ub = lb;
    }
    else if( sense == 'G' )
    {
        ub = OPT_INFINITY;
    }
    else //if( sense == 'R' ) //ranged constraint...
    {
        double range;
        r = CPXgetrngval(env, prob, &range, index, index);
        OPT_IFERRORRETURN(r, OPT_SOLVER_ERROR);
        
        
        //second by cplex manual the range can be positive or negative. 
        if( range >= 0 )
        {
            ub = lb + range;
        }
        else
        {
            ub = lb;
            lb = ub + range;
        }
         
    }
    
    
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value)
#if OPT_HAVE_CPLEX
{
    int r = CPXgetcoef(env, prob, constrIndex, varIndex, &value);
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values)
#if OPT_HAVE_CPLEX
{
    int matbeg[2];
    int dummy;
    
    const int r = CPXgetrows( env, prob, &nzs, matbeg, cols, values, naux, &dummy, constrIndex, constrIndex );
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::getLinearCoefsInConstraints(int &nzs, int *rows, int *cols, double *values)
#if OPT_HAVE_CPLEX
{
    //here we assume arrays have enough space
    
    int m, rcode;
    int *rmatbeg = NULL;
    
    nzs = 0;
    
    m = CPXgetnumrows(env, prob);
    
    if(m > 0)
    {
        const int rmatspace = INT_MAX;
        int r, k, dummy;
        
        OPT_malloc(rmatbeg, m+1);
        OPT_IFMEMERRORGOTOLABEL(!rmatbeg, rcode, termination);
        
        r = CPXgetrows(env, prob, &nzs, rmatbeg, cols, values, rmatspace, &dummy, 0, m-1);
        OPT_IFERRORRETURN(r, OPT_SOLVER_ERROR);
        
        k = 0;
        
        for(int i = 0; i < m; i++)
        {
            const int size = (i < m - 1 ? rmatbeg[i+1] : nzs )  -  rmatbeg[i];
            OPT_setAllArray( size,  &rows[k], i);
            k += size;
        }
        
        #if OPT_DEBUG_MODE
            assert(k == nzs);
        #endif
    }
    
    rcode = 0;
    
termination:

    if(rmatbeg) free(rmatbeg);
    
    return rcode;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::getObjConstant(double &objConstant) 
#if OPT_HAVE_CPLEX
{
    objConstant = this->objConstant;
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cplex::getObjLinearCoef( const int index, double &value )
#if OPT_HAVE_CPLEX
{
    int r = CPXgetobj( env, prob, &value, index, index );
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::getNumberOfConstraints(int &m)
#if OPT_HAVE_CPLEX
{
    m = CPXgetnumrows(env, prob) + CPXgetnumqconstrs(env, prob);
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex:: getNumberOfLinearCoefsInConstraints( int& nzs)
#if OPT_HAVE_CPLEX
{
    int nzcnt_p;
    int m;
    nzs = 0;
    
    m = CPXgetnumrows(env, prob);
    
    if(m > 0)
    {
        int *rmatbeg = auxIndex;
        
        int r = CPXgetrows(env, prob, &nzcnt_p, rmatbeg, NULL, NULL, 0, &nzs, 0, m-1);
        nzs = -nzs;
        
        if(r != 0 && r != CPXERR_NEGATIVE_SURPLUS)
            OPT_IFERRORRETURN(r, OPT_SOLVER_ERROR);
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif

int OPT_Cplex::getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs)
#if OPT_HAVE_CPLEX
{
    int nzcnt_p;
    int rmatbeg[2];
    int m = CPXgetnumrows(env, prob);
    
    nzs = 0;
    
    if( constrIndex < 0 || constrIndex >= m )
        return OPT_BAD_INPUT;
    
    
    CPXgetrows(env, prob, &nzcnt_p, rmatbeg, NULL, NULL, 0, &nzs, constrIndex, constrIndex); // CPLEX will return an error because we put NULL arrays. Anyway, it will return for us the value in nzs having the negative value of non zero elements...
    
    nzs = -nzs;
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::getNumberOfIntVars(int &nI)
#if OPT_HAVE_CPLEX
{
    nI = CPXgetnumint( env, prob );
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif





int OPT_Cplex::getNumberOfVars(int &n)
#if OPT_HAVE_CPLEX
{
    n = CPXgetnumcols( env, prob );
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::getObjSense(OPT_OPTSENSE &sense)
#if OPT_HAVE_CPLEX
{
    int v = CPXgetobjsen(env, prob);
    
    if( v == CPX_MIN )
        sense = optsolvers::OPT_MINIMIZE;
    else if ( v == CPX_MAX )
        sense = optsolvers::OPT_MAXIMIZE;
    else
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(v);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cplex::getVariableBounds(const int index, double &lb, double &ub)
#if OPT_HAVE_CPLEX
{
    int r = CPXgetub( env, prob, &ub, index, index );
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    if(ub >= CPX_INFBOUND)
        ub = OPT_INFINITY;
    
    
    r = CPXgetlb( env, prob, &lb, index, index );
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    if( lb <= -CPX_INFBOUND )
        lb = -OPT_INFINITY;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cplex::removeConstraints(const int ninds, const int *indices )
#if OPT_HAVE_CPLEX
{
    //cplex differentiates linear and quadratic constraints. One more good reason to do not provide suport to QCP problems on cplex
    
    //cplex requires an array having at least the size of number of constraint... :/
    
    int r;
    int m = CPXgetnumrows(env, prob);
    int *delstat = auxIndex;
    
    r = getNumberOfConstraints(m);
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return OPT_SOLVER_ERROR;
    }
    
    
    //only constraints marked with 1 in delstat will be removed
    OPT_setAllArray<int>(m, delstat, 0);
    
    for(int i = 0; i < ninds; i++)
        delstat[ indices[i] ] = 1;
    
    r = CPXdelsetrows(env, prob, delstat);
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return OPT_SOLVER_ERROR;
    }
    
    
    #if 0
    {
    //this version is wrong, beacuse delete by interval. But, after the first removing, other constraints indices can have changed, and so, we will remove the wrong constraints.
    int r, stran = 0;
    
    #if OPT_DEBUG_MODE
        int m, m0 = 0;
        getNumberOfConstraints(m0);
    #endif
    
    
    
    for(int i = 0; i < ninds-1; i++ )
    {
        if( indices[i] != indices[i+1] - 1 )
        {
            r = CPXdelrows(env, prob, indices[stran], indices[i] );
            
            if(r != 0)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                    
                return OPT_SOLVER_ERROR;
            }
            
            stran = i+1;
        }
        
    }
    
    
    r = CPXdelrows(env, prob, indices[stran], indices[ninds-1] );
            
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    #if OPT_DEBUG_MODE
        getNumberOfConstraints(m);
        
        assert( m + ninds == m0 );
    #endif
    }
    #endif
    
    return 0;
    
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



//all constraints in the range [begin end] will be removed
int OPT_Cplex::removeConstraintsByRange(const int begin, const int end)
#if OPT_HAVE_CPLEX
{
    const int r = CPXdelrows(env, prob, begin, end);
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cplex::removeVars(const int ninds, const int *indices )
#if OPT_HAVE_CPLEX
{
    int r, stran = 0;
    
    
    for(int i = 0; i < ninds-1; i++ )
    {
        if( indices[i] != indices[i+1] - 1 )
        {
            r = CPXdelcols(env, prob, indices[stran], indices[i] );
            
            if(r != 0)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                return OPT_SOLVER_ERROR;
            }
            
            stran = i+1;
        }
    }
    
    
    r = CPXdelcols(env, prob, indices[stran], indices[ninds-1] );
            
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::setLinearColumn( const int varIndex, const int nzs, const int* rows, const double* values)
#if OPT_HAVE_CPLEX
{
    int r;
    int size;
    
    
    if( naux > 0 )
    {
        for(int i = 0; i < naux; i++)
            auxIndex[i] = varIndex;
        
        
        int times = nzs/naux;
        
        if( times*naux < nzs )
            times++;
        
        size = naux;
        
        for( int i = 0; i < times; i++)
        {
            if( i == times - 1 ) //last time
                size = nzs - i*naux;
            
            
            r = CPXchgcoeflist( env, prob, size, &rows[i*naux] , auxIndex, &values[i*naux] );
            
            if(r != 0)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                    
                return OPT_SOLVER_ERROR;
            }
        }
        
    }
    else
    {
        for(int i = 0; i < nzs; i++)
        {
            r = CPXchgcoef(env, prob, rows[i], varIndex, values[i] );
            
            if(r != 0)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                    
                return OPT_SOLVER_ERROR;
            }
        }
    }
    
    
    
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




char OPT_Cplex::constraintSense( const double lb, const double ub )
{
    char value;
    
    if(lb > -OPT_INFINITY)
    {
        if(ub < OPT_INFINITY)
        {
            if(lb == ub)
                value = 'E';
            else
                value = 'R';
        }
        else
            value = 'G';
    }
    else
        value = 'L';
    
    
    return value;
}



int OPT_Cplex::resetConstraintLinearPart( const int constrIndex, const int nzs, const int *cols, const double *values )
#if OPT_HAVE_CPLEX
{
    int nzcnt, r;
    int nz, code ;
    
    
    //getting the number of nonzeros in the row...
    CPXgetrows( env, prob, &nzcnt, &r, NULL, NULL, 0, &nz, constrIndex, constrIndex );
    
    nz = -nz; //nz has a negative value of nonzeros in this constraints..
    
    if( nz > 0 )
    {
        int n;
        
        getNumberOfVars(n);
        
        
        
        //we have to replace predefined coefficinents, so we set the n columns in this constraints...
        
        OPT_setAllArray(n, auxIndex, constrIndex);
        OPT_setAllArray(n, auxValues, 0.0);
        
        
        for(int i = 0; i < n; i++)
            auxIndex2[i] = i;
        
        for(int i = 0; i < nzs; i++)
            auxValues[ cols[i] ] = values[i];
        
        
        r = CPXchgcoeflist(env, prob, n, auxIndex, auxIndex2, auxValues);
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            code = OPT_SOLVER_ERROR;
            goto termination;
        }
        
    }
    else
    {
    
        for(int i = 0; i < nzs; i++ )
            auxIndex[i] = constrIndex;
        
        
        r = CPXchgcoeflist(env, prob, nzs, auxIndex, cols, values);
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            code = OPT_SOLVER_ERROR;
            goto termination;
        }
    
    }
    
    
    
    //we are avoid call our API, but this process of setting constraints bound is little dark to us. So, we do this call at least until we are sure it works fine...
    //return setConstraintBounds(constrIndex, lb, ub);
    
    
    code = 0;
    
termination:
    
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::setConstraintBounds( const int index, const double lb, const double ub )
#if OPT_HAVE_CPLEX
{
    const char sense = constraintSense(lb, ub);
    int r;
    double rhs;
    
    
    
    r = CPXchgsense(env, prob, 1, &index, &sense);
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    rhs = sense == 'L' ? ub : lb;
    
    r = CPXchgrhs(env, prob, 1, &index, &rhs);
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    
    if( sense == 'R' )
    {
        double rngval = ub -lb;// second explanations in CPXnewrows (cplex manual), ranged constraints are between rhs and rhs + rngval.
        
        r = CPXchgrngval(env, prob, 1, &index, &rngval);
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            return OPT_SOLVER_ERROR;
        }
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif





int OPT_Cplex::setConstraintsLinearCoefs( const int nzs, const int* rows, const int* cols, const double* values )
#if OPT_HAVE_CPLEX
{
    int r;
    
    r = CPXchgcoeflist( env, prob, nzs, rows, cols, values );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::setConstraintLinearCoefs( const int constrIndex, const int nzs, const int *cols, const double *values)
#if OPT_HAVE_CPLEX
{
    int r;
    
    for(int i = 0; i < nzs; i++)
        auxIndex[i] = constrIndex;
    
    r = CPXchgcoeflist(env, prob, nzs, auxIndex, cols, values);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value)
#if OPT_HAVE_CPLEX
{
    int r = CPXchgcoef( env, prob, constrIndex, varIndex, value );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif





int OPT_Cplex::setObjLinearCoef( const int index, const double value )
#if OPT_HAVE_CPLEX
{
    int r = CPXchgobj( env, prob, 1, &index, &value );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::setObjLinearCoefs( const int nzs, const int* cols, const double* values )
#if OPT_HAVE_CPLEX
{
    int r = CPXchgobj( env, prob, nzs, cols, values );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::setObjLinearPart( const int n, const double *values )
#if OPT_HAVE_CPLEX
{
    const int on = CPXgetnumcols( env, prob );
    int r;
    
    if( on < n )
        return OPT_BAD_INPUT;
    
    for(int i = 0; i < on; i++)
        auxIndex[i] = i;
    
    
    r = CPXchgobj( env, prob, n, auxIndex, values );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    if(on != n)
    {
        for(int i = n; i < on; i++)
            auxValues[i] = 0.0;
        
        r = CPXchgobj(env, prob, on-n, &auxIndex[n], &auxValues[n] );
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            return OPT_SOLVER_ERROR;
        }
        
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif





void OPT_Cplex::setObjConstant(const double value)
{
    objConstant = value;
}




int OPT_Cplex::setObjSense( const OPT_OPTSENSE sense )
{
#if OPT_HAVE_CPLEX
    CPXchgobjsen( env, prob, sense == OPT_MINIMIZE ? CPX_MIN : CPX_MAX );
    
    return 0;
    
#else
    OPT_LIBNOTAVAILABLERET(getSolverCode());
#endif
}



int OPT_Cplex::setnVariablesBounds( const int n, const double *lb, const double *ub )
#if OPT_HAVE_CPLEX
{
    int r;
    char *lu = (char *) auxIndex2;
    
    
    for(int i = 0; i < n; i++)
        auxIndex[i] = i;
    
    for(int i = 0; i < n; i++)
        lu[i] = 'L';
    
    
    r = CPXchgbds(env, prob, n, auxIndex, lu, lb);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    for(int i = 0; i < n; i++)
        lu[i] = 'U';
    
    r = CPXchgbds(env, prob, n, auxIndex, lu, ub);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::setVariableBounds( const int index, const double lb, const double ub )
#if OPT_HAVE_CPLEX
{
    int r;
    
    const double mylb = lb <= -OPT_INFINITY ? -CPX_INFBOUND : lb ;
    
    const double myub = ub >= OPT_INFINITY ? CPX_INFBOUND : ub;
    
    
    r = CPXchgbds(env, prob, 1, &index, "L", &mylb);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    r = CPXchgbds(env, prob, 1, &index, "U", &myub);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::setVariablesBounds( const int ninds, const int *inds, const double *lb, const double *ub )
#if OPT_HAVE_CPLEX
{
    int r;
    char *lu = (char *) auxIndex2;
    
    
    for(int i = 0; i < ninds; i++)
        lu[i] = 'L';
    
    
    r = CPXchgbds(env, prob, ninds, inds, lu, lb);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    for(int i = 0; i < ninds; i++)
        lu[i] = 'U';
    
    r = CPXchgbds(env, prob, ninds, inds, lu, ub);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



/*
int OPT_Cplex::setVariableLowerBound( const int index, const double lb )
#if OPT_HAVE_CPLEX
{
    int r;
    
    const double mylb = lb <= -OPT_INFINITY ? -CPX_INFBOUND : lb ;
    
    r = CPXchgbds(env, prob, 1, &index, "L", &mylb);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Cplex::setVariableUpperBound( const int index, const double ub )
#if OPT_HAVE_CPLEX
{
    int r;
    
    const double myub = ub >= OPT_INFINITY ? CPX_INFBOUND : ub;
    
    r = CPXchgbds(env, prob, 1, &index, "U", &myub);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif
*/





//__ methods from QPSolver __


int OPT_Cplex::getNumberOfQuadObjTerms(int &nzs)
#if OPT_HAVE_CPLEX
{
    nzs = CPXgetnumqpnz( env, prob );
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Cplex::getObjQuadTerm( const int row, const int col, double &value)
#if OPT_HAVE_CPLEX
{
    int r = CPXgetqpcoef( env, prob, row, col, &value );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::getObjQuadPart( int &nzs, int *rows, int *cols, double *values )
#if OPT_HAVE_CPLEX
{
    const int n = CPXgetnumcols(env, prob);
    int code, dummy, r;
    int mynzs, *myrows = NULL, *mycols = NULL;
    double *myvalues= NULL;
    
    //I am not sure if qmatbeg should be n or  n+1 size...
    int *qmatbeg = NULL;
    
    
    
    //cplex DOES return full matrix instead of lower triangle! So, we have to fix the indexes! I HATE CPLEX!
    
    //getting the number of nonzero in full matrix Q
    r = CPXgetquad( env, prob, &dummy, NULL, NULL, NULL, 0, &mynzs, 0, n-1 ); //cplex will return a CPXERR_NEGATIVE_SURPLUS
    if( r != 0 && r != CPXERR_NEGATIVE_SURPLUS )
        OPT_PRINTERRORMSG("Warning: CPLEX did not retunr the expected return code CPXERR_NEGATIVE_SURPLUS. Check the manual!");
    
    
    mynzs = -mynzs;
    
    //cout << "mynzs: " << mynzs<< endl;
    
    
    qmatbeg = (int *) malloc( (n+1) * sizeof(int) );
    myrows = (int *) malloc( mynzs * sizeof(int) );
    mycols = (int *) malloc( mynzs * sizeof(int) );
    myvalues = (double *) malloc( mynzs * sizeof(double) );
    
    if( !qmatbeg || !myrows || !mycols || !mycols)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        
        code = OPT_MEMORY_ERROR;
        goto termination;
    }
    
    
    
    r = CPXgetquad(env, prob, &nzs, qmatbeg, myrows, myvalues,  mynzs, &dummy, 0, n-1 );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        code = OPT_SOLVER_ERROR;
        goto termination;
    }
    
    
    for(int i = 0; i < n; i++)
    {
        int start = qmatbeg[i];
        int end = i < n-1 ? qmatbeg[i+1] : nzs;
        
        for(int j = start; j < end; j++)
            mycols[j] = i;
    }
    
    //now, we discrad the upper triangle...
    
    for(int i = mynzs-1; i >= 0; i--)
    {
        if( mycols[i] > myrows[i] )
        {
            //upper triangle...
            for( int j = i+1; j < mynzs; j++ )
            {
                mycols[j-1] = mycols[j];
                myrows[j-1] = myrows[j];
                myvalues[j-1] = myvalues[j];
            }
            
            mynzs--;
        }
    }
    
    nzs = mynzs;
    
    
    for(int i = 0; i < nzs; i++)
    {
        rows[i] = myrows[i];
        cols[i] = mycols[i];
        values[i] = myvalues[i];
    }
    
    
    
    code = 0;
    
termination:
    
    
    if( qmatbeg )	free(qmatbeg);
    if( myrows )	free( myrows );
    if( mycols )	free( mycols );
    if( myvalues )	free( myvalues );
    
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::setObjQuadCoef( const int row, const int col, const double value )
#if OPT_HAVE_CPLEX
{
    int r = CPXchgqpcoef( env, prob, row, col, value );
    
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::deleteQuadObjMatrix()
#if OPT_HAVE_CPLEX
{
    unsigned int nzq = CPXgetnumqpnz(env, prob);
    
    if(nzq > 0)
    {//deleting old coefficients
        const int n = CPXgetnumcols(env, prob);
        int* rowStart = auxIndex;
        
        OPT_setAllArray(n, rowStart, 0);
        
        //we do not use CPXcopyquad directly to set the matrix because that function needs receive the complete matrix, not only the lower triangle (so stupid!).
        const int r = CPXcopyquad(env, prob, rowStart, rowStart, rowStart, auxValues); //you cannot pass a NULL pointer to cplex. Even if the matrix size is zero
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            return OPT_SOLVER_ERROR;
        }
        
        assert(CPXgetnumqpnz(env, prob) == 0);
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::setObjQuadMatrix( const int nzs, const int* rows, const int* cols, const double* values )
#if OPT_HAVE_CPLEX
{
    
    
    #if 0
    const int n = CPXgetnumcols(env, prob);
    bool rowsInOrder=true;
    int r, code;
    int* mycols = NULL;
    const int* pcols;
    double *myvalues = NULL;
    const double *pvalues;
    
    int *rowSize = auxIndex;
    int *rowStart = auxIndex2;
    
    //now, this method replaces the coeficients already stored
    
    if(n < 1)
        return 0;
    
    
    OPT_setAllArray(n, rowSize, 0);
    
    for(int i = 0; i < nzs; i++)
        rowSize[ rows[i] ]++;
    
    
    
    for(int i = 1; i < nzs; i++)
    {
        if( rows[i] < rows[i-1] )
        {
            rowsInOrder = false;
            break;
        }
    }
    
    
    rowStart[0] = 0;
    for(int i = 1; i < n; i++)
        rowStart[i] = rowStart[i-1] + rowSize[i-1];
    
    #if OPT_DEBUG_MODE
        assert( rowStart[n-1] + rowSize[n-1] == nzs );
    #endif
    
    
    if(rowsInOrder)
    {
        pcols = cols;
        pvalues = values;
    }
    else
    {
        mycols = (int *) malloc( nzs * sizeof(int) );
        myvalues = (double*) malloc( nzs*sizeof(double) );
        
        if( !mycols || !myvalues )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            code = OPT_MEMORY_ERROR;
            goto termination;
        }
        
        pcols = mycols;
        pvalues = myvalues;
        
        #if OPT_DEBUG_MODE
            OPT_setAllArray<double>(nzs, myvalues, (double) NAN);
        #endif
        
        OPT_setAllArray(n, rowSize, 0);
        
        for(int i = 0; i < nzs; i++)
        {
            const int row = rows[i];
            
            const int ind = rowStart[row] + rowSize[row];
            
            mycols[ind] = cols[i];
            myvalues[ind] = values[i];
            
            rowSize[row]++;
        }
        
        #if OPT_DEBUG_MODE
            for(int i = 0; i < nzs; i++)
                assert( !::isnan(pvalues[i]) );
        #endif
    }
    
    
    
    
    r = CPXcopyquad(env, prob, rowStart, rowSize, pcols, pvalues);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = r;
        goto termination;
    }
    
    
    
    code = 0;
    
termination:
    
    if(mycols)		free(mycols);
    if(myvalues)	free(myvalues);
    
    return code;
    #endif
    
    
    int r = deleteQuadObjMatrix();
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    for(int i = 0; i < nzs; i++)
    {
        r = CPXchgqpcoef( env, prob, rows[i], cols[i], values[i] );
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            return OPT_SOLVER_ERROR;
        }
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Cplex::setObjQuadMatrix(const int *rowStart, const int *cols, const double *values)
#if OPT_HAVE_CPLEX
{
    const int n = CPXgetnumcols(env, prob);
    
    int r = deleteQuadObjMatrix();
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    //const unsigned int nzs = rowStart[n]-rowStart[0];
    
    for(int i = 0; i < n; i++)
    {
        const int rnzs = rowStart[i+1]-rowStart[i];
        
        const int *rcols = &cols[ rowStart[i] ];
        const double *rvalues = &values[ rowStart[i] ];
        
        for(int j = 0; j < rnzs; j++)
        {
            r = CPXchgqpcoef( env, prob, i, rcols[j], rvalues[j] );
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                    
                return OPT_SOLVER_ERROR;
            }
        }
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


