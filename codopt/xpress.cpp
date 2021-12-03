

#include <cstdlib>

#include <cstring>
#include <cassert>

#include <iostream>

#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"


using namespace optsolvers;

using namespace std;



inline int OPT_getIntFromStr(const char *s)
{
    return atoi(s);
}





OPT_Xpress::OPT_Xpress():OPT_QCPSolver()
{
    initialize();
}



OPT_Xpress::~OPT_Xpress()
{
    deallocateSolverEnv();
}





// __methods from Solver __

void OPT_Xpress:: deallocateSolverEnv()
{
#if OPT_HAVE_XPRESS
    XPRSdestroyprob(prob);
#endif
    deallocateMemory();
}


int OPT_Xpress::getNumberOfIterations(long unsigned int& niter)
#if OPT_HAVE_XPRESS
{
    int alg, myniter = 0, r = 0;
    int nI;
    
    niter = 0;
    
    getNumberOfIntVars(nI);
    
    
    if(nI == 0)
    {
        XPRSgetintattrib(prob, XPRS_ALGORITHM, &alg);
        
        switch(alg)
        {
        case 2: //DUal Simplex
        case 3: //Primal simplex
        case 5: //Network simplex
            r = XPRSgetintattrib(prob, XPRS_SIMPLEXITER, &myniter);
            break;
        case 4: //Newton Barrier
            r = XPRSgetintattrib(prob, XPRS_BARITER, &myniter);
            break;
        case 1: //No LP optimization yet. 
            break;
        default:
            OPT_PRINTERRORMSG("Strange algorithm option");
        }
    }
    else
    {
        r = XPRSgetintattrib(prob, XPRS_NODES, &myniter);
    }
    
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    niter = myniter;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


OPT_LISTSOLVERS OPT_Xpress::getSolverCode()
{
    return OPT_XPRESS;
}


int OPT_Xpress::getVariableType( const int index, OPT_VARTYPE &varType )
#if OPT_HAVE_XPRESS
{
    char vt;
    
    const int r = XPRSgetcoltype(prob, &vt, index, index );
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    varType = vt == 'I' ? OPT_VT_INTEGER : OPT_VT_CONTINUOUS;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


void OPT_Xpress:: initialize()
{
#if OPT_HAVE_XPRESS
    origSolverRetCode = XPRS_LP_UNSTARTED;
    prob = NULL;
#endif
}



int OPT_Xpress::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
#if OPT_HAVE_XPRESS
{
    //that is ridiculous, but xpress has a global function to initialize the library. Xpress manual guarantees several threads can made call to XPRSinit and the initialization will be done just one time... Let us see... 
    
    int r;
    
    r = XPRSinit(NULL); 
    
    if( r )
    {
        char banner[256];
        
        XPRSgetbanner(banner);
        
        cerr << OPT_PREPRINT << "Error xpress banner: " << banner << endl;
    }
    
    r = XPRScreateprob(&prob);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    //just put it to avoi xpress print information abut empty problem being incorporated...
    r = XPRSsetintcontrol( prob, XPRS_OUTPUTLOG, XPRS_OUTPUTLOG_ERRORS_AND_WARNINGS );
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    //that is ridiculous, but if we do not read a problem form file, we have to load the matrices. Since int his moment we have no matrices to load, we are obligated to create a empty problem having no variables nor constraints...
    r = XPRSloadlp(prob, "qcpXpressOptSolvers", 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    
    
    
    
    
    
    r = XPRSsetintcontrol( prob, XPRS_IFCHECKCONVEXITY, 0 );
    
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




int OPT_Xpress::setObjCutLowerBound(const double objLBound)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSsetdblcontrol( prob, XPRS_MIPABSCUTOFF, objLBound );
    
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




int OPT_Xpress::setObjCutUpperBound(const double objUBound)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSsetdblcontrol( prob, XPRS_MIPABSCUTOFF, objUBound );
    
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




int OPT_Xpress::setMaxCPUTime(const double time)
#if OPT_HAVE_XPRESS
{
    //const int r = XPRSsetdblcontrol(prob, XPRS_MAXTIME, -time); //we need put negative values to stop after time seconds...
    
    const int r = XPRSsetintcontrol(prob, XPRS_MAXTIME, -time); //we need put negative values to stop after time seconds...
    
    
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




int OPT_Xpress::setNumberOfThreads(const int nthreads)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSsetintcontrol(prob, XPRS_THREADS, nthreads == 0 ? -1 : nthreads);
    
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




int OPT_Xpress::setOutputLevel( const int level )
#if OPT_HAVE_XPRESS
{
    int l;
    
    if( level <= 0)
        l = XPRS_OUTPUTLOG_NO_OUTPUT;
    else if( level >= 4 )
        l = XPRS_OUTPUTLOG_FULL_OUTPUT;
    else if( level == 3 )
        l = XPRS_OUTPUTLOG_ERRORS_AND_WARNINGS;
    else
        l = XPRS_OUTPUTLOG_ERRORS;
    
    
    const int r = XPRSsetintcontrol( prob, XPRS_OUTPUTLOG, l ); //XPRS_MIPLOG
    
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



int OPT_Xpress::setRelativeDualTol( const double tol )
#if OPT_HAVE_XPRESS
{
    int r;
    
    r = XPRSsetdblcontrol(prob, XPRS_BARDUALSTOP, tol);
    
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


int OPT_Xpress::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_XPRESS
{
    int r;
    
    r = XPRSsetdblcontrol(prob, XPRS_BARGAPSTOP, tol);
    
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



int OPT_Xpress::setRelativePrimalTol( const double tol )
#if OPT_HAVE_XPRESS
{
    return OPT_OPERATION_NOT_IMPLEMENTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Xpress::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_XPRESS
{
    //const int pn = OPT_getIntFromStr( param );
    int r, iHeaderId, iTypeInfo;
    
    //that is not the ideal, but since xpress works setting number to parameters, we will retrieve parameter number from string "param".
    
    r = XPRSgetcontrolinfo(prob, param, &iHeaderId, &iTypeInfo);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        printDblParamErrorMsg(r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    
    r = XPRSsetdblcontrol(prob, iHeaderId, value);
    
    if ( r != 0 )
    {
        printDblParamErrorMsg(r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Xpress::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_XPRESS
{
    //const int pn = OPT_getIntFromStr( param );
    int r, iHeaderId, iTypeInfo;
    
    
    r = XPRSgetcontrolinfo(prob, param, &iHeaderId, &iTypeInfo);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        printIntParamErrorMsg(r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    
    r = XPRSsetintcontrol(prob, iHeaderId, value);
    
    if ( r != 0 )
    {
        printIntParamErrorMsg(r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Xpress::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_XPRESS
{
    //const int pn = OPT_getIntFromStr( param );
    int r, iHeaderId, iTypeInfo;
    
    
    r = XPRSgetcontrolinfo(prob, param, &iHeaderId, &iTypeInfo);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        printStrParamErrorMsg(r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    
    r = XPRSsetstrcontrol(prob, iHeaderId, value);
    
    if ( r != 0 )
    {
        printStrParamErrorMsg(r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Xpress::setVariableType( const int index, const OPT_VARTYPE varType )
#if OPT_HAVE_XPRESS
{
    const int r = XPRSchgcoltype(prob, 1, &index, varType == OPT_VT_INTEGER ? "I" : "C" );
    
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




int OPT_Xpress::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_XPRESS
{
    int r, r2, r3, status, nI = 0;
    double *psol = storeSol ? sol : NULL;
    double *pconstr = storeConstrs ? constr : NULL;
    
    if(resetSol)
    {
        origSolverRetCode = XPRS_LP_UNSTARTED;
        this->resetSol();
    }
    else
    {
        feasSol = false;
    }
    
    getNumberOfIntVars(nI);
    
    
    if( nI > 0)
    {
        r = XPRSmipoptimize(prob, "");
        
        r2 = XPRSgetmipsol(prob, psol, pconstr);
        
        r3 = XPRSgetdblattrib(prob, XPRS_MIPOBJVAL, &objValue);
        
        
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            retCode = OPT_SOLVER_ERROR;
            goto termination;
        }
        
        
        if(r2 != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            retCode = OPT_SOLVER_ERROR;
            goto termination;
        }
        else
        {
            feasSol = true;
        }
        
        #if OPT_DEBUG_MODE
            if(r3 != 0)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
            }
        #endif
        
        
        
        
        
        r = XPRSgetintattrib(prob, XPRS_MIPSTATUS, &status);
        
        #if OPT_DEBUG_MODE
            if(r != 0)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                    
                retCode = OPT_SOLVER_ERROR;
                goto termination;
            }
        #endif
        
        origSolverRetCode = status;
        
        switch(status)
        {
            case XPRS_MIP_OPTIMAL:
                feasSol = true;
                retCode = OPT_OPTIMAL_SOLUTION;
                break;
            
            case XPRS_MIP_INFEAS:
                retCode = OPT_INFEASIBLE_PROBLEM;
                break;
                
            case XPRS_MIP_UNBOUNDED:
                retCode = OPT_UNBOUNDED_PROBLEM;
                break;
                
            case XPRS_MIP_SOLUTION:
                feasSol = true;
                retCode = OPT_FEASIBLE_SOLUTION;
                break;
                
            default:
                retCode = OPT_UNDEFINED_ERROR;
        }
        
    }
    else
    {
        double *pdualSolC = storeDualSol ? dualSolC : NULL;
        double *pdualSolV = storeDualSol ? dualSolV : NULL;
        
        
        r = XPRSlpoptimize(prob, "");
        
        
        r2 = XPRSgetlpsol(prob, psol, pconstr, pdualSolC, pdualSolV);
        
        r3 = XPRSgetdblattrib(prob, XPRS_LPOBJVAL, &objValue);
        
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            retCode = OPT_SOLVER_ERROR;
            goto termination;
        }
        
        if(r2 != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            retCode = OPT_SOLVER_ERROR;
            goto termination;
        }
        else
        {
            feasSol = true;
        }
        
        
        #if OPT_DEBUG_MODE
            if(r3 != 0)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
            }
        #endif
        
        
        
        
        
        
        r = XPRSgetintattrib(prob, XPRS_LPSTATUS, &status);
        
        #if OPT_DEBUG_MODE
            if(r != 0)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                    
                retCode = OPT_SOLVER_ERROR;
                goto termination;
            }
        #endif
        
        origSolverRetCode = status;
        
        switch(status)
        {
            case XPRS_LP_OPTIMAL:
                feasSol = true;
                retCode = OPT_OPTIMAL_SOLUTION;
                break;
                
            case XPRS_LP_INFEAS:
            case XPRS_LP_CUTOFF:
                retCode = OPT_INFEASIBLE_PROBLEM;
                break;
                
            case XPRS_LP_UNBOUNDED:
                retCode = OPT_UNBOUNDED_PROBLEM;
                break;
            
            case XPRS_LP_NONCONVEX:
                retCode = OPT_INCOMPATIBLE_SOLVER;
                break;
            
            default:
                retCode = OPT_UNDEFINED_ERROR;
        }
        
        
        
        
        
    }
    
    
    
    if( storeConstrs )
    {
        int m;
        char bf;
        double rhs;
        
        //I am not sure if this part is correct because XPRESS is not available to me to test...
        
        getNumberOfConstraints(m);
        
        for(int i = 0; i < m; i++)
        {
            r = XPRSgetrowtype( prob, &bf, i, i);
            #if OPT_DEBUG_MODE
                if( r != 0 )
                {
                    OPT_PRINTERRORNUMBER(r);
                    return OPT_SOLVER_ERROR;
                }
            #endif
            
            r = XPRSgetrhs(prob, &rhs, i, i);
            #if OPT_DEBUG_MODE
                if( r != 0 )
                {
                    OPT_PRINTERRORNUMBER(r);
                    return OPT_SOLVER_ERROR;
                }
            #endif
            
            if( bf == 'G' )
                constr[i] = rhs + constr[i];
            else if( bf != 'N' )
                constr[i] = rhs - constr[i];
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



int OPT_Xpress::__addConstraints(const int nm)
#if OPT_HAVE_XPRESS
{
    char *qrtype;
    int r, code;
    double* const rhs = auxValues;
    
    qrtype =  (char *) ( (void *) auxIndex );
    
    
    OPT_setAllArray(nm, qrtype, 'N');
    OPT_setAllArray(nm, rhs, 0.0);
    
    
    r = XPRSaddrows(prob, nm, 0, qrtype, rhs, NULL, NULL, NULL, NULL);
    
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
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Xpress::__addVariables(const int nn, const bool initFree)
#if OPT_HAVE_XPRESS
{
    int r, code;
    double *lb, *ub, *bounds = NULL;
    double *c = auxValues;
    
    
    
    
    OPT_setAllArray(nn, c, 0.0);
    
    if( initFree )
    {
        bounds = (double *) malloc( 2*nn * sizeof(double) );
    
        if(!bounds)
        {
            code = OPT_MEMORY_ERROR;
            goto termination;
        }
        
        lb = bounds;
        ub = &bounds[nn];
        
        OPT_setAllArray(nn, lb, XPRS_MINUSINFINITY);
        OPT_setAllArray(nn, ub, XPRS_PLUSINFINITY );
    }
    else
    {
        ub = lb = c; // we take advantage c and set variables fixing all to 0.0
    }
    
    
    r = XPRSaddcols(prob, nn, 0, c, NULL, NULL, NULL, lb, ub);
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return OPT_SOLVER_ERROR;
    }
    
    
    //XPRSchgbounds(prob, nn);
    
    code = 0;
termination:
    
    if(bounds)	free(bounds);
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif





int OPT_Xpress::generateModelFile(const char* fileName)
#if OPT_HAVE_XPRESS
{
    XPRSsetprobname(prob, "qcpXpressOptSolvers");
    
    const int r = XPRSwriteprob(prob, fileName, "pl");
    
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



int OPT_Xpress::getConstraintBounds( const int index, double &lb, double &ub )
#if OPT_HAVE_XPRESS
{
    int r;
    char bf;
    
    
    r = XPRSgetrhs( prob, &lb, index, index );
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    r = XPRSgetrowtype( prob, &bf, index, index );
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    if( bf == 'L' )
    {
        ub = lb;
        lb = -OPT_INFINITY;
    }
    else if( bf == 'E' )
    {
        ub = lb;
    }
    else if( bf == 'G' )
    {
        ub = OPT_INFINITY;
    }
    else if( bf == 'R' )
    {
        //ranged constraint
        
        r = XPRSgetrhsrange(prob, &ub, index, index);
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
        
        ub = ub + lb;
    }
    else //if( bf == 'N' )
    {
        lb = -OPT_INFINITY;
        ub = OPT_INFINITY;
    }
    
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Xpress::getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSgetcoef(prob, constrIndex, varIndex, &value);
    
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



int OPT_Xpress::getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values)
#if OPT_HAVE_XPRESS
{
    int mstart[2];
    
    const int r = XPRSgetrows( prob, mstart, cols, values, naux, &nzs, constrIndex, constrIndex );
    
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



int OPT_Xpress::getObjLinearCoef( const int index, double &value )
#if OPT_HAVE_XPRESS
{
    const int r = XPRSgetobj(prob, &value, index, index);
    
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




int OPT_Xpress::getNumberOfConstraints(int &m)
#if OPT_HAVE_XPRESS
{
    /*int r, iHeaderId, iTypeInfo;
    
    r = XPRSgetattribinfo( prob, "ROWS", &iHeaderId, &iTypeInfo );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif
        
        return OPT_SOLVER_ERROR;
    } */
    
    
    
    const int r = XPRSgetintattrib( prob, XPRS_ORIGINALROWS, &m );
    
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



int OPT_Xpress::getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSgetrows( prob, NULL, NULL, NULL, 0, &nzs, constrIndex, constrIndex );
    
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



int OPT_Xpress::getNumberOfIntVars(int &n)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSgetintattrib( prob, XPRS_ORIGINALMIPENTS, &n );
    
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




int OPT_Xpress::getNumberOfVars(int &n)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSgetintattrib(prob, XPRS_ORIGINALCOLS, &n );
    
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



int  OPT_Xpress:: getObjConstant(double &objConstant)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSgetobj(prob, &objConstant, -1, -1); //-1 index is the constant term
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
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



int OPT_Xpress::getObjSense(OPT_OPTSENSE &sense)
#if OPT_HAVE_XPRESS
{
    double s;
    const int r = XPRSgetdblattrib( prob, XPRS_OBJSENSE, &s );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    sense = s == XPRS_OBJ_MINIMIZE ? OPT_MINIMIZE : OPT_MAXIMIZE;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Xpress::getVariableBounds(const int index, double &lb, double &ub)
#if OPT_HAVE_XPRESS
{
    int r;
    
    r = XPRSgetlb(prob, &lb, index, index);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    r = XPRSgetub(prob, &ub, index, index);
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



int OPT_Xpress::removeConstraints(const int ninds, const int *indices )
#if OPT_HAVE_XPRESS
{
    const int r = XPRSdelrows( prob, ninds, indices );
    
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



int OPT_Xpress::removeVars(const int ninds, const int *indices)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSdelcols( prob, ninds, indices );
    
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



int OPT_Xpress::setLinearColumn( const int varIndex, const int nzs, const int* rows, const double* values)
#if OPT_HAVE_XPRESS
{
    int r;
    
    OPT_setAllArray(nzs, auxIndex, varIndex);
    
    r = XPRSchgmcoef(prob, nzs, rows, auxIndex, values);
    
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





int OPT_Xpress::resetConstraintLinearPart( const int constrIndex, const int nzs, const int *cols, const double *values )
#if OPT_HAVE_XPRESS
{
    int code, r;
    
    
    
    int nz = XPRSgetrows(prob, NULL, NULL, NULL, 0, &r, constrIndex, constrIndex);
    
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
        
        
        
        r = XPRSchgmcoef(prob, n, auxIndex, auxIndex2, auxValues);
        
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
        OPT_setAllArray(nzs, auxIndex, constrIndex);
        
        r = XPRSchgmcoef(prob, nzs, auxIndex, cols, values);
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            code = OPT_SOLVER_ERROR;
            goto termination;
        }
    }
    
    
    code = 0;
    
termination:
    
    return code; //setConstraintBounds(constrIndex, lb, ub);
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



char OPT_Xpress::boundFlag(const double lb, const double ub)
{
    char qrtype;
    
    if( lb > -OPT_INFINITY )
    {
        if( ub < OPT_INFINITY )
        {
            if( lb == ub )
                qrtype = 'E';
            else
                qrtype = 'R';
        }
        else
            qrtype = 'G';
    }
    else
    {
        if( ub < OPT_INFINITY )
            qrtype = 'L';
        else
            qrtype = 'N';
    }
    
    return qrtype;
}



int OPT_Xpress::setConstraintBounds( const int index, const double lb, const double ub )
#if OPT_HAVE_XPRESS
{
    int r;
    const char bf = boundFlag(lb, ub);
    
    
    r = XPRSchgrowtype(prob, 1, &index, &bf);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    r = XPRSchgrhs(prob, 1, &index, bf == 'G' ? &lb: &ub);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    if( bf == 'R' )
    {
        const double range = ub - lb;
        
        r = XPRSchgrhsrange(prob, 1, &index, &range);
        
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





int OPT_Xpress::setConstraintsLinearCoefs( const int nzs, const int* rows, const int* cols, const double* values )
#if OPT_HAVE_XPRESS
{
    const int r = XPRSchgmcoef( prob, nzs, rows, cols, values );
    
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




int OPT_Xpress::setConstraintLinearCoefs( const int constrIndex, const int nzs, const int *cols, const double *values)
#if OPT_HAVE_XPRESS
{
    int r;
    
    
    OPT_setAllArray( nzs, auxIndex, constrIndex );
    
    r = XPRSchgmcoef( prob, nzs, auxIndex, cols, values );
    
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




int OPT_Xpress::setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSchgcoef(prob, constrIndex, varIndex, value);
    
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





int OPT_Xpress::setObjLinearCoef( const int index, const double value )
#if OPT_HAVE_XPRESS
{
    const int r = XPRSchgobj(prob, 1, &index, &value);
    
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




int OPT_Xpress::setObjLinearCoefs( const int nzs, const int* cols, const double* values )
#if OPT_HAVE_XPRESS
{
    const int r = XPRSchgobj(prob, nzs, cols, values);
    
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




int OPT_Xpress::setObjLinearPart( const int n, const double *values )
#if OPT_HAVE_XPRESS
{
    int r, on;
    
    getNumberOfVars(on);
    
    if(n > on)
        return OPT_BAD_INPUT;
    
    
    
    for(int i = 0; i < on; i++)
        auxIndex[i] = i;
    
    r = XPRSchgobj(prob, n, auxIndex, values);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    if( n != on )
    {
        for(int i = n; i < on; i++)
            auxValues[i] = 0.0;
        
        
        r = XPRSchgobj(prob, on-n, &auxIndex[n], &auxValues[n] );
        
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





void OPT_Xpress:: setObjConstant(const double value)
{
#if OPT_HAVE_XPRESS
    
    const int cInd = -1;
    double nvalue = -value; //I do not know because, but we have to reverse the signal...
    //cout << "obj const value: " << value << endl;
    
    const int r = XPRSchgobj(prob, 1, &cInd, &nvalue);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
    #endif
    
#endif
}




int OPT_Xpress::setObjSense( const OPT_OPTSENSE sense )
#if OPT_HAVE_XPRESS
{
    const int r = XPRSchgobjsense(prob, sense == OPT_MINIMIZE ? XPRS_OBJ_MINIMIZE : XPRS_OBJ_MAXIMIZE);
    
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





int OPT_Xpress::setVariableBounds( const int index, const double lb, const double ub )
#if OPT_HAVE_XPRESS
{
    int code = 0;
    int r;
    double v;
    
    v = ub < OPT_INFINITY ? ub : XPRS_PLUSINFINITY;
    
    r = XPRSchgbounds(prob, 1, &index, "U", &v );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR;
    }
    
    
    
    v = lb > -OPT_INFINITY ? lb : XPRS_MINUSINFINITY;
    
    r = XPRSchgbounds(prob, 1, &index, "L", &v );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR;
    }
    
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif






//__ methods from QPSolver __



int OPT_Xpress::getNumberOfQuadObjTerms(int &nzs)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSgetqrowqmatrixtriplets( prob, -1, &nzs, NULL, NULL, NULL );
    
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



int OPT_Xpress::getObjQuadPart( int &nzs, int *rows, int *cols, double *values )
#if OPT_HAVE_XPRESS
{
    //so crazy: xpress manual says that function returns lower triangular part, but, actually, ir returns upper triangular part. So, we hve to change row and col indexes...
    const int r = XPRSgetqrowqmatrixtriplets( prob, -1, &nzs, cols, rows, values );
    
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



int OPT_Xpress::getObjQuadTerm( const int row, const int col, double &value)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSgetqobj(prob, row, col, &value);
    
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




int OPT_Xpress::setObjQuadCoef( const int row, const int col, const double value )
#if OPT_HAVE_XPRESS
{
    const int r = XPRSchgqobj(prob, row, col, value);
    
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




int OPT_Xpress::setObjQuadMatrix( const int nzs, const int* rows, const int* cols, const double* values )
#if OPT_HAVE_XPRESS
{
    const int r = XPRSchgmqobj(prob, nzs, rows, cols, values);
    
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


int OPT_Xpress::setObjQuadMatrix(const int *rowStart, const int *cols, const double *values)
#if OPT_HAVE_XPRESS
{
    int n, code = 0;
    int *myrows = auxIndex;
    
    
    getNumberOfVars(n);
    
    //deleting the current q obj matrix
    XPRSdelqmatrix(prob, -1);
    
    
    for(int i = 0; i < n; i++)
    {
        const int rowSize = rowStart[i+1] - rowStart[i];
        
        if(rowSize > 0)
        {
            const int ind = rowStart[i];
            
            OPT_setAllArray(rowSize, myrows, i);
            
            const int r = XPRSchgmqobj(prob, rowSize, myrows, &cols[ind], &values[ind]);
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                code = OPT_SOLVER_ERROR;
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


// __ methods from QCPSolver __



int OPT_Xpress::getNumberOfConstraintQuadTerms( const int index, int &nzs)
#if OPT_HAVE_XPRESS
{
    const int r = XPRSgetqrowqmatrixtriplets( prob, index, &nzs, NULL, NULL, NULL );
    
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



int OPT_Xpress::getConstraintQuadMatrix( const int index, int &nzs, int *rows, int *cols, double *values )
#if OPT_HAVE_XPRESS
{
    //althoug manual says xpress returns the lower triangle, actually it returns the upper triangle. So, we have to change row and column indexes
    
    const int r = XPRSgetqrowqmatrixtriplets( prob, index, &nzs, cols, rows, values );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    //we have to multiply by 2 because xpress do not divide this matrix by 2 like other solvers, note we had to multiply entrues by 0.5 in setQuadConstraintPart...
    
    for(int i = 0; i < nzs; i++)
        values[i] *= 2.0;
    
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Xpress::setConstraintQuadMatrix( const int index, const int nzs, const int *qrows, const int *qcols, const double *qvalues )
#if OPT_HAVE_XPRESS
{
    //bool isq = false;
    int r;
    //int nq = 0;
    int code = 0;
    double *hqvalues = NULL; //we have to multiply values by 0.5
        
    
    #if 0
    r = XPRSgetqrows( prob, &nq, auxIndex );
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    for( int i = 0; i < nq; i++ )
    {
        if( auxIndex[i] == index  )
        {
            isq = true;
            break;
        }
    }
    
    
    
    if(isq)
    {
        
        
        // we have to change coefficinents 
        for(int i = 0; i < nzs; i++)
        {
            //we have to multiply by 0.5 becuase xpress do not divide by 2 like other solvers...
            
            r = XPRSchgqrowcoeff(prob, index, qrows[i], qcols[i], 0.5*qvalues[i] );
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                code = OPT_SOLVER_ERROR;
            }
        }
        
    }
    #endif
    
    //deleting the current q matrix in the constraint
    XPRSdelqmatrix(prob, index);
    
    
    {
        //constraint has no quadratic terms yet. So, we can add all matrix in a once time...
        
        if( nzs <= naux || nzs <= maux )
        {
            hqvalues = auxValues;
        }
        else
        {
            hqvalues = (double *) malloc( nzs * sizeof(double) );
            
            if(!hqvalues)
            {
                code = OPT_MEMORY_ERROR;
                goto termination;
            }
        }
        
        for(int i = 0; i < nzs; i++)
            hqvalues[i] = 0.5 * qvalues[i];
        
        
        r = XPRSaddqmatrix(prob, index, nzs, qrows, qcols, hqvalues);
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            code = OPT_SOLVER_ERROR;
            goto termination;
        }
    }
    
    
termination:
    
    if( hqvalues && hqvalues != auxValues )
        free(hqvalues);
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



/*int OPT_Xpress::setQuadConstraintMatrix( const int index, const int *qrowStart, const int *qcols, const double *qvalues )
#if OPT_HAVE_XPRESS
{
    int n, code = 0;
    int *myrows = auxIndex;
    getNumberOfVars(n);
    
    return OPT_OPERATION_NOT_IMPLEMENTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif */












