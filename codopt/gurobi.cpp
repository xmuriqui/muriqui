

#include <cstdlib>
#include <cassert>
#include <climits>
#include <math.h>

#include <iostream>
#include <list>

#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"


//we set flag bellow while gurobi add auxiliary variables to range constraints... (I HATE GUROBI!)
#define OPT_GUROBI_USE_AUX_VARS 1



using namespace std;
using namespace optsolvers;




OPT_Gurobi::OPT_Gurobi():OPT_QPSolver()
{
    initialize();
}




        
OPT_Gurobi::~OPT_Gurobi()
{
    deallocateSolverEnv();
}





// __methods from Solver __



void OPT_Gurobi::deallocateMemory()
{
    OPT_LPSolver::deallocateMemory();
    
    OPT_secFree(indAuxVarConstr);
}


void OPT_Gurobi::deallocateSolverEnv()
{
#if OPT_HAVE_GUROBI
    if( prob )
    {
        GRBfreemodel( prob );
        prob = NULL;
    }
    
    if( env )
    {
        GRBfreeenv( env );
        env = NULL;
    }
#endif

    deallocateMemory();
}


int OPT_Gurobi::getNumberOfIterations(long unsigned int &niter)
#if OPT_HAVE_GUROBI
{
    int nI;
    
    getNumberOfIntVars(nI);
    
    if(nI == 0)
    {
        double value1 = 0, value2 = 0;
        
        GRBgetdblattr(prob, "IterCount", &value1);
        GRBgetdblattr(prob, "BarIterCount", &value2);
        
        niter = OPT_max(value1, value2);
    }
    else
    {
        double value = 0;
        
        GRBgetdblattr(prob, "NodeCount", &value);
        niter = value;
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


OPT_LISTSOLVERS OPT_Gurobi::getSolverCode()
{
    return OPT_GUROBI;
}


int OPT_Gurobi::getVariableType( const int index, OPT_VARTYPE &varType )
#if OPT_HAVE_GUROBI
{
    char vt;
    
    const int r = GRBgetcharattrelement( prob, GRB_CHAR_ATTR_VTYPE, index, &vt );
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    varType = vt == GRB_INTEGER ? OPT_VT_INTEGER : OPT_VT_CONTINUOUS;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


void OPT_Gurobi::initialize()
{
    OPT_QPSolver::initialize();
    
#if OPT_HAVE_GUROBI
    origSolverRetCode = INT_MAX;
    env = NULL;
    prob = NULL;
#endif
    
    indAuxVarConstr = NULL;
}





int OPT_Gurobi::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
#if OPT_HAVE_GUROBI
{
    int r;
    
    r = GRBloadenv(&env, NULL);
    if(r)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        if( r == GRB_ERROR_NO_LICENSE )
        {
            OPT_PRINTERRORMSG("Gurobi failed to find a valid license.");
            return OPT_LICENSE_ERROR;
        }
        
        return OPT_SOLVER_ERROR;
    }
    
    
    
    r = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
    #if OPT_DEBUG_MODE
        if(r)
        {
            OPT_PRINTERRORNUMBER(r);
        }
    #endif
    
    //setNumberOfThreads(1);
    
    
    r = GRBnewmodel(env, &prob, "qcpGurobiOptSolvers", 0, NULL, NULL, NULL, NULL, NULL);
    if(r)
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




int OPT_Gurobi::setObjCutLowerBound(const double objLBound)
#if OPT_HAVE_GUROBI
{
    OPT_OPTSENSE sense;
    int r;
    
    getObjSense( sense );
    
    if( sense != OPT_MAXIMIZE )
        return OPT_BAD_INPUT;
    
    r = GRBsetdblparam(env, GRB_DBL_PAR_CUTOFF, objLBound);
    if(r)
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




int OPT_Gurobi::setObjCutUpperBound(const double objUBound)
#if OPT_HAVE_GUROBI
{
    OPT_OPTSENSE sense;
    int r;
    
    getObjSense( sense );
    
    if( sense != OPT_MINIMIZE )
        return OPT_BAD_INPUT;
    
    r = GRBsetdblparam(env, GRB_DBL_PAR_CUTOFF, objUBound);
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



int OPT_Gurobi::setMaxCPUTime(const double time)
#if OPT_HAVE_GUROBI
{
    //we should not get the original environment, sice parameter changing does not afect it. We have to use GRBgetenv
    GRBenv *env = GRBgetenv(prob);
    const int r = GRBsetdblparam(env, GRB_DBL_PAR_TIMELIMIT, time);
    
    #if OPT_DEBUG_MODE
        if(r != 0)
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



int OPT_Gurobi::setNumberOfThreads(const int nthreads)
#if OPT_HAVE_GUROBI
{
    //we should not get the original environment, sice parameter changing does not afect it. We have to use GRBgetenv
    GRBenv *env = GRBgetenv(prob);
    int r = GRBsetintparam(env, GRB_INT_PAR_THREADS, nthreads);
    
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



int OPT_Gurobi::setOutputLevel( const int level )
#if OPT_HAVE_GUROBI
{
    //we should not get the original environment, sice parameter changing does not afect it. We have to use GRBgetenv
    GRBenv *env = GRBgetenv(prob);
    int r = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, level > 0 ? 1 : 0);
    
    #if OPT_DEBUG_MODE
        if(r != 0)
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



int OPT_Gurobi::setRelativeDualTol( const double tol )
#if OPT_HAVE_GUROBI
{
    return OPT_OPERATION_NOT_IMPLEMENTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Gurobi::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_GUROBI
{
    //we should not get the original environment, sice parameter changing does not afect it. We have to use GRBgetenv
    GRBenv *env = GRBgetenv(prob);
    
    int r = GRBsetdblparam(env, GRB_DBL_PAR_BARCONVTOL, tol);
    #if OPT_DEBUG_MODE
        if(r != 0)
        {
            OPT_PRINTERRORNUMBER(r);
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    r = GRBsetdblparam(env, GRB_DBL_PAR_BARQCPCONVTOL, tol);
    #if OPT_DEBUG_MODE
        if(r != 0)
        {
            OPT_PRINTERRORNUMBER(r);
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    r = GRBsetdblparam(env, GRB_DBL_PAR_MIPGAP, tol);
    #if OPT_DEBUG_MODE
        if(r != 0)
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



int OPT_Gurobi::setRelativePrimalTol( const double tol )
#if OPT_HAVE_GUROBI
{
    return OPT_OPERATION_NOT_IMPLEMENTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Gurobi::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_GUROBI
{
    //we should not get the original environment, since parameter changing does not afect it. We have to use GRBgetenv
    GRBenv *env = GRBgetenv(prob);
    const int r = GRBsetdblparam(env, param, value);
    
    
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



int OPT_Gurobi::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_GUROBI
{
    //we should not get the original environment, sice parameter changing does not afect it. We have to use GRBgetenv
    GRBenv *env = GRBgetenv(prob);
    const int r = GRBsetintparam(env, param, value);
    
    
    if( r != 0 )
    {
        printIntParamErrorMsg(r,  param, value );
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Gurobi:: setStringParameter(const char *param, const char *value)
#if OPT_HAVE_GUROBI
{
    //we should not get the original environment, sice parameter changing does not afect it. We have to use GRBgetenv
    GRBenv *env = GRBgetenv(prob);
    const int r = GRBsetstrparam(env, param, value);
    
    
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




/*
int OPT_Gurobi::setParameters( OPT_GeneralSolverParams &params )
#if OPT_HAVE_GUROBI
{
    int r, code = 0;
    char solver[] = "gurobi";
    
    std::list< OPT_SolverParam<int>* >::iterator itInt;
    std::list< OPT_SolverParam<double>* >::iterator itDbl;
    std::list< OPT_SolverParam<char*>* >::iterator itStr;
    
    
    for(itInt = params.intParams.begin(); itInt != params.intParams.end(); itInt++)
    {
        r = GRBsetintparam(env, (*itInt)->name, (*itInt)->value);
        if( r != 0 )
        {
            printIntParamErrorMsg( r, solver, (*itInt)->name, (*itInt)->value );
            
            code = OPT_BAD_INPUT;
        }
    }
    
    
    for(itDbl = params.dblParams.begin(); itDbl != params.dblParams.end(); itDbl++)
    {
        r = GRBsetdblparam(env, (*itDbl)->name, (*itDbl)->value );

        if(r != 0)
        {
            printDblParamErrorMsg( r, solver, (*itDbl)->name, (*itDbl)->value );
            
            code = OPT_BAD_INPUT;
        }
    }
    
    
    for(itStr = params.strParams.begin(); itStr != params.strParams.end(); itStr++)
    {
        r = GRBsetstrparam(env, (*itStr)->name, (*itStr)->value );

        if(r != 0)
        {
            printStrParamErrorMsg( r, solver, (*itStr)->name, (*itStr)->value );
            
            code = OPT_BAD_INPUT;
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



int OPT_Gurobi::setVariableType( const int index, const OPT_VARTYPE varType )
#if OPT_HAVE_GUROBI
{
    const int r = GRBsetcharattrelement( prob, GRB_CHAR_ATTR_VTYPE, index,  varType == OPT_VT_INTEGER ? GRB_INTEGER : GRB_CONTINUOUS ) ;
    
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




int OPT_Gurobi::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_GUROBI
{
    int r, status;
    int n, m;
    
    
    if(resetSol)
    {
        origSolverRetCode = INT_MAX;
        this->resetSol();
    }
    else
    {
        feasSol = false;
    }
    
    r = GRBoptimize(prob);
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        retCode = OPT_SOLVER_ERROR;
        goto termination;
    }
    
    
    r = GRBgetintattr(prob, GRB_INT_ATTR_STATUS, &status);
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        retCode = OPT_SOLVER_ERROR;
        goto termination;
    }
    
    getMyNumberOfVars( n );
    
    
    GRBgetintattr(prob, GRB_INT_ATTR_NUMCONSTRS, &m);
    
    
    if( storeSol )
    {
        GRBgetdblattrarray(prob, GRB_DBL_ATTR_X, 0, n, sol); // we use naux becuase the toal number of variables contains the artifical variables introduced by range constraints (I hate Gurobi)
    }
    
    if( storeConstrs )
    {
        //constraint values are in the arificial variables
        for(int i = 0; i < m; i++)
        {
            GRBgetdblattrarray(prob, GRB_DBL_ATTR_X, indAuxVarConstr[i], 1, &constr[i] );
            
            constr[i] = -constr[i];
        }
    }
    
    
    
    GRBgetdblattr(prob, GRB_DBL_ATTR_OBJVAL, &objValue);
    
    //dualObjValue = objValue;
    
    if( storeDualSol )
    {
        GRBgetdblattrarray(prob, GRB_DBL_ATTR_PI, 0, m, dualSolC);
        
        GRBgetdblattrarray(prob, GRB_DBL_ATTR_RC, 0, n, dualSolV);
    }
    
    origSolverRetCode = status;
    
    
    switch( status )
    {
    case GRB_OPTIMAL:
    case GRB_SUBOPTIMAL:
        
        feasSol = true;
        retCode = OPT_OPTIMAL_SOLUTION;
        break;
        
    case GRB_INFEASIBLE:
        retCode = OPT_INFEASIBLE_PROBLEM;
        break;
        
    case GRB_INF_OR_UNBD:
        retCode = OPT_UNDEFINED_ERROR;
        break;
        
    case GRB_UNBOUNDED:
        
        retCode = OPT_UNBOUNDED_PROBLEM;
        break;
        
    case GRB_ITERATION_LIMIT:
        
        #if OPT_PRINT_MAX_ITER_WARNING
            if( numberOfWarningsByIterLimit < maxNumberOfWarningsByIterLimit )
            {
                std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Gurobi solving!\n";
                numberOfWarningsByIterLimit++;
                
                if( numberOfWarningsByIterLimit == maxNumberOfWarningsByIterLimit )
                    std::cerr << OPT_PREPRINT "Warning: Maximum number of warnings by maximum iteration achieved! Stopping these warnings.\n";
            }
        #endif
        
        retCode = OPT_MAX_ITERATIONS;
        break;
        
    case GRB_TIME_LIMIT:
        
        retCode = OPT_MAX_TIME;
        break;
        
    default:
        retCode = OPT_UNDEFINED_ERROR;
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



int OPT_Gurobi::__addConstraints(const int nm)
#if OPT_HAVE_GUROBI
{
    //char *sense = NULL;
    int r, code, n;
    int *auxp;
    double *lb = NULL, *ub;
    
    auxp = (int *) realloc( indAuxVarConstr, maux * sizeof(int) );
    //sense = (char *) malloc( nm * sizeof(char) );
    lb = (double *) malloc( 2*nm * sizeof(double) );
    
    if( !auxp || !lb )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        
        code = OPT_MEMORY_ERROR;
        goto termination;
    }
    
    indAuxVarConstr = auxp;
    ub = &lb[nm];
    
    //for(int i = 0; i < nm; i++)
        //sense[i] = GRB_LESS_EQUAL;
    
    OPT_setAllArray<double>( nm, lb, 0.0 );
    OPT_setAllArray<double>( nm, ub, 0.0 );
    
    
    getNumberOfVars(n); //here, we need the number of vars including auxiliary variables generated by constraints...
    
    //r = GRBaddconstrs( prob, nm, 0, NULL, NULL, NULL, sense, NULL, NULL );
    r = GRBaddrangeconstrs( prob, nm, 0, NULL, NULL, NULL, lb, ub, NULL );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR;
        goto termination;
    }
    
    GRBupdatemodel(prob);
    
    
    
    
    //setting the indices of auxiliaries variables used to guroby in this ranged constraints...
    auxp = &indAuxVarConstr[ maux - nm ];
    for(int i = 0; i < nm; i++)
        auxp[i] = n + i;
    
    
    code = 0;
    
termination:
    
    //if( sense ) free(sense);
    
    if(lb)		free(lb);
    
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Gurobi::__addVariables(const int nn, const bool initFree)
#if OPT_HAVE_GUROBI
{
    int r, m;
    double *lb;
    
    if( initFree )
    {
        lb = auxValues;
        for(int i = 0; i < nn; i++)
            lb[i] = -GRB_INFINITY;
    }
    else
    {
        lb = NULL;
    }
    
    #if OPT_DEBUG_MODE
        r = getNumberOfConstraints(m);
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
        
        if( m > 0 )
        {
            std::cout << OPT_PREPRINT << "warning: Adding variables in gurobi after constraints. Atention to new variables indices.\n";
        }
    #endif
    
    r = GRBaddvars( prob, nn, 0, NULL, NULL, NULL, NULL, lb, NULL, NULL, NULL );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    GRBupdatemodel(prob);
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif





int OPT_Gurobi::generateModelFile(const char* fileName)
#if OPT_HAVE_GUROBI
{
    int r;
    
    GRBupdatemodel(prob);
    
    r = GRBwrite(prob, fileName );
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



int OPT_Gurobi::getConstraintBounds( const int index, double &lb, double &ub )
#if OPT_HAVE_GUROBI
{
    int code = 0, r;
    
    
    //note we have to set -ub as lower bound because coefficient of auxiliary variable in the constraint is 1.0 instead of -1.0
    r = GRBgetdblattrelement( prob, GRB_DBL_ATTR_LB, indAuxVarConstr[index], &ub );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR; 
    }
    
    ub = -ub;
    
    if(ub >= GRB_INFINITY)
        ub = OPT_INFINITY;
    
    
    //note we have to set -lb as upper bound because coefficient of auxiliary variable in the constraint is 1.0 instead of -1.0
    r = GRBgetdblattrelement( prob, GRB_DBL_ATTR_UB, indAuxVarConstr[index], &lb );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR; 
    }
    
    lb = -lb;
    if( lb <= -GRB_INFINITY )
        lb = -OPT_INFINITY;
    
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Gurobi::getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value)
#if OPT_HAVE_GUROBI
{
    GRBupdatemodel(prob);
    
    int r = GRBgetcoeff(prob, constrIndex, varIndex, &value);
    
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


int OPT_Gurobi::getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values)
#if OPT_HAVE_GUROBI
{
    int cbeg[2];
    
    const int r = GRBgetconstrs( prob, &nzs, cbeg, cols, values, constrIndex, 1 );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    //we have to disconsider the auxiliary variable added by range constraints... aff!
    
    //we run by the end because is most probable the searched column be the last...
    for(int i = nzs-1; i >= 0; i-- )
    {
        if( cols[i] == indAuxVarConstr[constrIndex] )
        {
            for(int j = i+1; j < nzs; j++)
            {
                cols[j-1] = cols[j];
                values[j-1] = values[j];
            }
            
            nzs--; 
            
            break;
        }
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Gurobi::getObjLinearCoef( const int index, double &value )
#if OPT_HAVE_GUROBI
{
    GRBupdatemodel(prob);
    
    int r = GRBgetdblattrelement(prob, "Obj", index, &value);
    
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




int OPT_Gurobi::getNumberOfConstraints(int &m)
#if OPT_HAVE_GUROBI
{
    GRBupdatemodel(prob);
    
    int r = GRBgetintattr(prob, GRB_INT_ATTR_NUMCONSTRS, &m);
    
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



int OPT_Gurobi::getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs)
#if OPT_HAVE_GUROBI
{
    const int r = GRBgetconstrs( prob, &nzs, NULL, NULL, NULL, constrIndex, 1 );
    
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



int OPT_Gurobi::getNumberOfIntVars(int &nI)
#if OPT_HAVE_GUROBI
{
    GRBupdatemodel(prob);
    
    int r = GRBgetintattr(prob, "NumIntVars", &nI);
    
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



int OPT_Gurobi::getMyNumberOfVars(int &n)
#if OPT_HAVE_GUROBI
{
    int r = getNumberOfVars(n);
    
    #if OPT_GUROBI_USE_AUX_VARS
        int m;
        r += getNumberOfConstraints(m);
        n = n - m; //we have to disconsider the auxiliary (I HATE GUROBI)
    #endif
    
    return r == 0 ? 0 : OPT_SOLVER_ERROR;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Gurobi::getNumberOfVars(int &n)
#if OPT_HAVE_GUROBI
{
    GRBupdatemodel(prob);
    
    int r = GRBgetintattr(prob, "NumVars", &n);
    
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



int OPT_Gurobi::getObjConstant(double &objConstant)
#if OPT_HAVE_GUROBI
{
    int r = GRBgetdblattr( prob, "ObjCon", &objConstant );
    
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


int OPT_Gurobi::getObjSense(OPT_OPTSENSE &sense)
#if OPT_HAVE_GUROBI
{
    GRBupdatemodel(prob);
    
    int s;
    int r = GRBgetintattr(prob, "ModelSense", &s);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    sense = s == GRB_MINIMIZE? OPT_MINIMIZE : OPT_MAXIMIZE;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Gurobi::getVariableBounds(const int index, double &lb, double &ub)
#if OPT_HAVE_GUROBI
{
    
    int r = GRBgetdblattrelement(prob, GRB_DBL_ATTR_LB, index, &lb);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    r = GRBgetdblattrelement(prob, GRB_DBL_ATTR_UB, index, &ub);
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



int OPT_Gurobi::removeConstraints(const int ninds, const int *indices )
#if OPT_HAVE_GUROBI
{
    GRBupdatemodel(prob);
    
    int *aindices = OPT_getPointerToConstArray( indices);
    
    const int r = GRBdelconstrs( prob, ninds, aindices );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    GRBupdatemodel(prob);
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Gurobi::removeVars(const int ninds, const int *indices )
#if OPT_HAVE_GUROBI
{
    GRBupdatemodel(prob);
    
    int *aindices = OPT_getPointerToConstArray( indices);
    
    const int r = GRBdelvars( prob, ninds, aindices );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    GRBupdatemodel(prob);
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Gurobi::setLinearColumn( const int varIndex, const int nzs, const int* rows, const double* values)
#if OPT_HAVE_GUROBI
{
    int r, code = 0;
    int *arows = OPT_getPointerToConstArray(rows);
    double *avalues = OPT_getPointerToConstArray(values);
    
    
    if( naux > 0 )
    {
        for(int i = 0; i < naux; i++)
            auxIndex[i] = varIndex;
        
        int times, size;
        
        times = nzs/naux;
        
        if( times*naux < nzs )
            times++;
        
        size = naux;
        
        for(int i = 0; i < times; i++)
        {
            
            if( i == times - 1 )
                size = nzs - i*naux;
            
            r = GRBchgcoeffs( prob, size, &arows[i*naux], auxIndex, &avalues[i*naux]);
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                code = OPT_SOLVER_ERROR; // we do not return because we let the remainder coeffs be set
            }
        }
        
    }
    else
    {
        int vindex;
        
        for(int i = 0; i < nzs; i++)
        {
            r = GRBchgcoeffs( prob, 1, &arows[i], &vindex, &avalues[i] );
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                code = OPT_SOLVER_ERROR; // we do not return because we let the remainder coeffs be set
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




int OPT_Gurobi::resetConstraintLinearPart( const int constrIndex, const int nzs, const int *cols, const double *values )
#if OPT_HAVE_GUROBI
{
    int r, n, code = 0;
    //int *auxIndex2 = NULL;
    //int *acols = OPT_getPointerToConstArray(cols );
    //double *avalues = OPT_getPointerToConstArray(values);
    
    
    getMyNumberOfVars(n);
    
    
    
    OPT_setAllArray(n, auxIndex, constrIndex);
    OPT_setAllArray(n, auxValues, 0.0);
    
    for(int i = 0; i < n; i++)
        auxIndex2[i] = i;
    
    for(int i = 0; i < nzs; i++)
        auxValues[ cols[i] ] = values[i];
    
    
    
    r = GRBchgcoeffs( prob, n, auxIndex, auxIndex2, auxValues);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR; 
    }
    
    
    /*r = setConstraintBounds(constrIndex, lb, ub); //that is a complicated process. We are avoid call interface methods, but I think here isa good idea...
    if( r != 0 )
        code = r; */
    
    
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Gurobi::setConstraintBounds( const int index, const double lb, const double ub )
#if OPT_HAVE_GUROBI
{
    int r, code = 0;
    
    
    
    //note we have to set -ub as lower bound because coefficient of auxiliary variable in the constraint is 1.0 instead of -1.0
    r = GRBsetdblattrelement(prob, GRB_DBL_ATTR_LB, indAuxVarConstr[index], ub < OPT_INFINITY ? -ub : -GRB_INFINITY );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR; 
    } 
    
    
    /*r = GRBsetdblattrelement(prob, GRB_DBL_ATTR_RHS, index, ub); //change rhs of the range constraint to ub. Nw, we let variable range from 0 to ub + lb (Note that auxiliary variable has coef 1.0 int this equality constraint)
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR; 
    } */
    
    
    
    //note we have to set -lb as upper bound because coefficient of auxiliary variable in the constraint is 1.0 instead of -1.0
    r = GRBsetdblattrelement(prob, GRB_DBL_ATTR_UB, indAuxVarConstr[index], lb > -GRB_INFINITY ? -lb : GRB_INFEASIBLE );
    
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






int OPT_Gurobi::setConstraintsLinearCoefs( const int nzs, const int* rows, const int* cols, const double* values )
#if OPT_HAVE_GUROBI
{
    int *arows = OPT_getPointerToConstArray(rows);
    int *acols = OPT_getPointerToConstArray(cols);
    double *avalues = OPT_getPointerToConstArray(values);
    
    
    int r = GRBchgcoeffs( prob, nzs, arows, acols, avalues);
    
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




int OPT_Gurobi::setConstraintLinearCoefs( const int constrIndex, const int nzs, const int *cols, const double *values)
#if OPT_HAVE_GUROBI
{
    int r;
    int *acols = OPT_getPointerToConstArray(cols);
    double *avalues = OPT_getPointerToConstArray(values);
    
    
    
    for(int i = 0; i < nzs; i++)
        auxIndex[i] = constrIndex;
    
    r = GRBchgcoeffs( prob, nzs, auxIndex, acols, avalues);
    
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




int OPT_Gurobi::setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value)
#if OPT_HAVE_GUROBI
{
    int r;
    
    int rindex = constrIndex;
    int vindex = varIndex;
    double val = value;
    
    
    r = GRBchgcoeffs(prob, 1, &rindex, &vindex, &val);
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





int OPT_Gurobi::setObjLinearCoef( const int index, const double value )
#if OPT_HAVE_GUROBI
{
    int r = GRBsetdblattrelement(prob, "Obj", index, value);
    
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




int OPT_Gurobi::setObjLinearCoefs( const int nzs, const int* cols, const double* values )
#if OPT_HAVE_GUROBI
{
    //I hate cplex and gurobi... input array arguments should be const, but they are not. So, I am obligated to doing this kind of terrible thing...
    
    int *acols = OPT_getPointerToConstArray(cols);
    double *avalues = OPT_getPointerToConstArray(values);
    
    
    int r = GRBsetdblattrlist(prob, "Obj", nzs, acols, avalues);
    
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




int OPT_Gurobi::setObjLinearPart( const int n, const double *values )
#if OPT_HAVE_GUROBI
{
    int on;
    double *avalues = OPT_getPointerToConstArray(values);
    
    getMyNumberOfVars(on);
    
    
    int r = GRBsetdblattrarray(prob, "Obj", 0, n, avalues);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    if( on != n )
    {
        const int size = on - n;
        
        OPT_setAllArray(size, auxValues, 0.0);
        
        r = GRBsetdblattrarray(prob, "Obj", n, size, auxValues);
        
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





void OPT_Gurobi::setObjConstant(const double value)
{
    #if OPT_HAVE_GUROBI
    int r = GRBsetdblattr( prob, "ObjCon", value );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
    }
    
    #endif
}




int OPT_Gurobi::setObjSense( const OPT_OPTSENSE sense )
#if OPT_HAVE_GUROBI
{
    int r = GRBsetintattr(prob, "ModelSense", sense == OPT_MINIMIZE ? GRB_MINIMIZE : GRB_MAXIMIZE);
    
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



int OPT_Gurobi::setnVariablesBounds( const int n, const double *lb, const double *ub )
#if OPT_HAVE_GUROBI
{
    int code = 0, r;
    
    
    r = GRBgetdblattrarray( prob, "LB", 0, n, OPT_getPointerToConstArray(lb) );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR; 
    }
    
    
    r = GRBgetdblattrarray( prob, "UB", 0, n, OPT_getPointerToConstArray(ub) );
    
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



int OPT_Gurobi::setVariableBounds( const int index, const double lb, const double ub )
#if OPT_HAVE_GUROBI
{
    int code= 0;
    
    int r = GRBsetdblattrelement(prob, "LB", index, lb);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR; 
    }
    
    
    r = GRBsetdblattrelement(prob, "UB", index, ub);
    
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



int OPT_Gurobi::setVariablesBounds( const int ninds, const int *inds, const double *lb, const double *ub )
#if OPT_HAVE_GUROBI
{
    int r, code = 0;
    
    r = GRBsetdblattrlist(prob, "LB", ninds, OPT_getPointerToConstArray(inds), OPT_getPointerToConstArray(lb) );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR; 
    }
    
    r = GRBsetdblattrlist(prob, "UB", ninds, OPT_getPointerToConstArray(inds), OPT_getPointerToConstArray(ub) );
    
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


int OPT_Gurobi::getNumberOfQuadObjTerms(int &nzs)
#if OPT_HAVE_GUROBI
{
    const int r = GRBgetintattr( prob, "NumQNZs", &nzs );
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


int OPT_Gurobi::getObjQuadPart( int &nzs, int *rows, int *cols, double *values )
#if OPT_HAVE_GUROBI
{
    //Gurobi return the upper triangle... so, we have to change indexes...
    const int r = GRBgetq( prob, &nzs, cols, rows, values );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    //we have to multiply diag termbs by 2 because Gurobi does not adopt the 0.5 before Q
    
    for(int i = 0; i < nzs; i++)
    {
        if( cols[i] == rows[i] )
            values[i] *= 2.0;
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Gurobi::getObjQuadTerm( const int row, const int col, double &value)
#if OPT_HAVE_GUROBI
{
    int r, code, nzq;
    
    int *qrow = NULL, *qcol = NULL;
    double *qvalue = NULL; //that is ridiculous
    
    
    value = NAN;
    
    
    GRBupdatemodel(prob);
    
    r = GRBgetintattr( prob, "NumQNZs", &nzq );
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR; 
        goto termination;
    }
    
    
    if( nzq > naux && nzq > maux )
    {
        qrow = (int *) malloc( nzq * sizeof(int) );
        qcol = (int *) malloc( nzq * sizeof(int) );
        qvalue = (double *) malloc( nzq * sizeof(double) );
        
        if( !qrow || !qcol || !qvalue )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            
            code = OPT_MEMORY_ERROR;
            goto termination;
        }
    }
    else
    {
        qrow = auxIndex;
        qcol = auxIndex2;
        qvalue = auxValues;
    }
    
    
    r = GRBgetq( prob, &nzq, qrow, qcol, qvalue );
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR; 
        goto termination;
    }
    
    
    for(int i = 0; i < nzq; i++ )
    {
        if( qrow[i] == row && qcol[i] == col )
        {
            value = qvalue[i];
            code = 0;
            goto termination;
        }
    }
    
    
    value = 0.0;
    
    code = 0;
    
termination:
    
    
    if( nzq > naux && nzq > maux )
    {	
        if(qrow)	free(qrow);
        if(qcol)	free(qcol);
        if(qvalue)	free(qvalue);
    }
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Gurobi::setObjQuadCoef( const int row, const int col, const double value )
#if OPT_HAVE_GUROBI
{
    int nzq = 0;
    int r, code;
    
    int myrow = row;
    int mycol = col;
    double myvalue = myrow == mycol ? 0.5*value : value; //Gurobi has a different way of handle this coeffcients... You have to put on values the real coefficients that you want on qp terms, i.e, you do not have to multiply diagonal terms by 2.
    
    int *qrows = NULL, *qcols;
    double *qvalues = NULL;
    
    
    if( myrow < mycol )
        OPT_swap(myrow, mycol);
    
    
    
    GRBupdatemodel(prob);
    
    
    GRBgetintattr( prob, "NumQNZs", &nzq );
    /*GRBgetintattr( prob, "NumVars", &n );
    
    #if OPT_GUROBI_USE_AUX_VARS
        int m;
        getNumberOfConstraints(m);
        n = n - m;
    #endif */
    
    
    if( nzq > 0 )
    {
        qrows = (int *) malloc(2*nzq * sizeof(int));
        qvalues = (double *) malloc( nzq * sizeof(double) );
        
        if( !qrows || !qvalues )
        {
            code = OPT_MEMORY_ERROR;
            goto termination;
        }
        
        qcols = &qrows[nzq];
        
        r = GRBgetq( prob, &nzq, qrows, qcols, qvalues );
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            code = OPT_SOLVER_ERROR; 
            goto termination;
        }
        
        
        for(int i = 0; i < nzq; i++)
        {
            if( qrows[i] < qcols[i] )
                OPT_swap( qrows[i], qcols[i] );
            
            if(qrows[i] == myrow && qcols[i] == mycol)
            {
                myvalue = myvalue -qvalues[i]; //gurobi add this term to therm already in the matrix, so, we have to do it to final value be the wanted value...
                break;
            }
        }
        
    }
        
        
    r = GRBaddqpterms( prob, 1, &myrow, &mycol, &myvalue );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR;
        goto termination;
    }
    
    
    
    
    code = 0;
    
termination:
    
    if( qrows )		free(qrows);
    if( qvalues )	free(qvalues);
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Gurobi::setObjQuadMatrix( const int nzs, const int* rows, const int* cols, const double* values )
#if OPT_HAVE_GUROBI
{
    int code, r, nzq = 0;
    int n;
    int *arows = OPT_getPointerToConstArray(rows );
    int *acols = OPT_getPointerToConstArray(cols );
    double *myvalues = NULL;
    
    
    int *pqrow = NULL;
    //int *pqcol;
    double *pqvalue = NULL;
    double *q = NULL; 
    
    
    //NumQNZs 
    
    GRBupdatemodel(prob);
    
    GRBgetintattr( prob, "NumQNZs", &nzq );
    GRBgetintattr( prob, "NumVars", &n );
    
    #if OPT_GUROBI_USE_AUX_VARS
        int m;
        getNumberOfConstraints(m);
        
        n -= m;
    #endif
    
    
    
    
    if( nzs > naux && nzs > maux )
    {
        myvalues = (double *) malloc( nzs * sizeof(double) );
        
        if( !myvalues )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            
            code = OPT_MEMORY_ERROR;
            goto termination;
        }
        
    }
    else
    {
        myvalues = auxValues;
    }
    
    
    
    for(int i = 0; i < nzs; i++)
    {
        if( arows[i] == acols[i] )
            myvalues[i] = 0.5*values[i];
        else
            myvalues[i] = values[i];
    }
    
    
    //Gurobi has no function to change the coefficients... just add...
    
    if( nzq > 0 )
    {
        #if 0
        pqrow = (int *) malloc( 2*nzq * sizeof(int) );
        pqvalue = (double *) malloc( nzq * sizeof(double) );
        q = (double *) calloc( (n*(n+1))/2, sizeof(double) );
        
        
        if( !pqrow || !pqvalue || !q )
        {
            code = OPT_MEMORY_ERROR;
            goto termination;
        }
        
        pqcol = &pqrow[nzq];
        
        
        r = GRBgetq( prob, &nzq, pqrow, pqcol, pqvalue );
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            code = OPT_SOLVER_ERROR; 
            goto termination;
        }
        
        
        //checking if we have to mix the current coefiicients and the new coefficients...
        for(int i = 0; i < nzq; i++)
        {
            if( pqrow[i] < pqcol[i] )
                OPT_swap( pqrow[i], pqcol[i] );
            
            
            int r = pqrow[i];
            int c = pqcol[i];
            
            q[ (r*(r+1))/2 + c ] = pqvalue[i];
        }
        
        /*cout << "q: " << endl;
        for( int i = 0, k = 0; i < n; i++ )
        {
            for( int j = 0; j <= i; j++, k++ )
                cout << q[k] << " ";
            cout << endl;
        } */
            
        
        for(int i = 0; i < nzs; i++)
        {
            int r = rows[i];
            int c = cols[i];
            
            if( r < c )
                OPT_swap(r, c);
            
            double v = q[ (r*(r+1))/2 + c ];
            
            if( v  ) //v != 0.0
            {
                myvalues[i] -=  v; //gurobi add this term to therm already in the matrix, so, we have to do it to final value be the wanted value...
                
                //cout << "Cai! i: " << i << " r: " << r << " c: " << c << endl;
            }
        }
        
        
        //we free the memory here to try avoid memory fragmentantion when gurobi have to allocate new structures in the call GRBaddqpterms below.
        free(pqrow); 	pqrow = NULL;
        free(pqvalue);	pqvalue = NULL;
        free(q);		q = NULL;
        
        #endif
        
        //gurobi just has functions to add values to quadratics coefficients. So, we delete all quadratic coefiicients
        r = GRBdelq(prob);
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            code = OPT_SOLVER_ERROR;
            goto termination;
        }
        
    }
    
    
    r = GRBaddqpterms(prob, nzs, arows, acols, myvalues);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR;
        goto termination;
    }
    
    //generateModelFile("optgurobi.lp");
    //cout << "Gerei optgurobi.lp em optsolvers" << endl;
    //getchar();
    
    
    
    code = 0;
    
    
termination:
    
    if(nzs > naux && nzs > maux )
    {
        if(myvalues)	free(myvalues);
    }
    
    if(pqrow)		free(pqrow);
    if(pqvalue)		free(pqvalue);
    
    if(q)			free(q);
    
    return code;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Gurobi::setObjQuadMatrix(const int *rowStart, const int *cols, const double *values)
#if OPT_HAVE_GUROBI
{
    int n, nzq = 0;
    int code = 0;
    int *myrows = auxIndex;
    double *myvalues = auxValues;
    
    
    GRBupdatemodel(prob);
    
    GRBgetintattr( prob, "NumQNZs", &nzq );
    GRBgetintattr( prob, "NumVars", &n );
    
    #if OPT_GUROBI_USE_AUX_VARS
        int m;
        getNumberOfConstraints(m);
        
        n -= m;
    #endif
    
    
    if( nzq > 0 )
    {
        //gurobi just has functions to add values to quadratics coefficients. So, we delete all quadratic coefiicients
        const int r = GRBdelq(prob);
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            code = OPT_SOLVER_ERROR;
            //we let thigns go on anyway and do not return now...
        }
    }
    
    //I am not sure if it is better set coeffcients one by one, or row by row... :(
    //by now, we are performing row by row
    for(int i = 0; i < n; i++)
    {
        const int rowSize = rowStart[i+1]-rowStart[i];
        const int ind = rowStart[i];
        int* rcols = OPT_getPointerToConstArray( &cols[ind] ) ;// &cols[ind];
        
        
        if( rowSize > 0 )
        {
            #if OPT_DEBUG_MODE
                assert( rowSize <= n );
            #endif
            
            OPT_setAllArray(rowSize, myrows, i);
            OPT_copyArray(rowSize, &values[ind], myvalues);
            
            for(int j = 0; j < rowSize; j++)
            {
                if( rcols[j] == i )
                {
                    myvalues[j] *= 0.5;
                    break; //we assume there is no repeated positions
                }
            }
            
            const int r = GRBaddqpterms(prob, rowSize, myrows, rcols, myvalues);
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                code = OPT_SOLVER_ERROR;
                //goto termination;
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


/*int OPT_Gurobi::setQuadConstraint( const int index, const int nzs, const int *qrows, const int *qcols, const double *qvalues )
#if OPT_HAVE_GUROBI
{
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif
*/



    
