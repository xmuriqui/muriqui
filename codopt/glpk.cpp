

#include <math.h>
#include <cstdlib>
#include <cassert>

#include <iostream>

#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"


#if OPT_HAVE_GLPK
#include "glpk.h"
#endif



using namespace std;
using namespace optsolvers;





OPT_Glpk:: OPT_Glpk():OPT_LPSolver()
{
    initialize();
}


OPT_Glpk::~OPT_Glpk()
{
    //desallocateMemory();
    deallocateSolverEnv();
}




// __methods from Solver __


void OPT_Glpk::deallocateMemory()
{
    OPT_LPSolver::deallocateMemory();

    //OPT_secFree(auxIndex2);
    //OPT_secFree(auxValues2);
}



void OPT_Glpk::deallocateSolverEnv()
{
#if OPT_HAVE_GLPK
    glp_delete_prob(prob);
    prob = NULL;

    OPT_secFree( simParam );
    OPT_secFree( intParam );
#endif
    
    deallocateMemory();
}


int OPT_Glpk::getNumberOfIterations(long unsigned int& niter)
#if OPT_HAVE_GLPK
{
    niter = glp_get_it_cnt(prob);
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


OPT_LISTSOLVERS OPT_Glpk::getSolverCode()
{
    return OPT_GLPK;
}


int OPT_Glpk::getVariableType( const int index, OPT_VARTYPE &varType )
#if OPT_HAVE_GLPK
{
    const int r = glp_get_col_kind(prob, index + 1);
    
    varType = r == GLP_IV? OPT_VT_INTEGER : OPT_VT_CONTINUOUS;
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



void OPT_Glpk::initialize()
{
    OPT_LPSolver::initialize();

#if OPT_HAVE_GLPK
    origSolverRetCode = GLP_UNDEF;
    prob = NULL;
    simParam = NULL;
    intParam = NULL;
#endif

    //auxIndex2 = NULL;
    //auxValues2 = NULL;

    //we put that because glpk counts indices starting from 1, and we need one more position in our arrays. So, we initialize naux and maux as 1 instead of 0;
    naux = 1;
    maux = 1;
}


int OPT_Glpk::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
#if OPT_HAVE_GLPK
{
    int code;

    prob = glp_create_prob();
    simParam = (glp_smcp *) malloc( sizeof(glp_smcp) );
    intParam = (glp_iocp *) malloc( sizeof(glp_iocp) );

    if( !prob || !simParam || !intParam )
    {
        code = OPT_MEMORY_ERROR;
        goto termination;
    }


    glp_init_smcp(simParam);
    glp_init_iocp(intParam);


    //enabling pre-processing
    simParam->presolve = GLP_ON;
    intParam->presolve = GLP_ON;


    //enabling feasibility heuristics for integer programming
    intParam->fp_heur = GLP_ON;
    intParam->ps_heur = GLP_ON;


    //enabling cuts in the B&B
    intParam->gmi_cuts = GLP_ON;
    intParam->mir_cuts = GLP_ON;
    intParam->cov_cuts = GLP_ON;
    intParam->clq_cuts = GLP_ON;

    setOutputLevel(1);

    code = 0;

termination:

    if( code != 0 )
        deallocateSolverEnv();



    return code;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setObjCutLowerBound(const double objLBound)
#if OPT_HAVE_GLPK
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setObjCutUpperBound(const double objUBound)
#if OPT_HAVE_GLPK
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setMaxCPUTime(const double time)
#if OPT_HAVE_GLPK
{
    const double mtime = time * 1000; //miliseconds
    simParam->tm_lim = mtime;
    intParam->tm_lim = mtime;

    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Glpk::setNumberOfThreads(const int nthreads)
#if OPT_HAVE_GLPK
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setOutputLevel( const int level )
#if OPT_HAVE_GLPK
{
    glp_term_out( level > 1 ? GLP_ON : GLP_OFF );


    if( level <= 0 )
        simParam->msg_lev = GLP_MSG_OFF;
    else if (level <= 1)
        simParam->msg_lev = GLP_MSG_ERR;
    else if (level <= 9)
        simParam->msg_lev = GLP_MSG_ON;
    else
        simParam->msg_lev = GLP_MSG_ALL;


    intParam->msg_lev = simParam->msg_lev;


    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setRelativeDualTol( const double tol )
#if OPT_HAVE_GLPK
{
    #if OPT_HAVE_GLPK
        intParam->mip_gap = tol;
    #endif
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Glpk::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_GLPK
{
    #if OPT_HAVE_GLPK
        return OPT_OPERATION_NOT_IMPLEMENTED;
    #endif
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Glpk::setRelativePrimalTol( const double tol )
#if OPT_HAVE_GLPK
{
    #if OPT_HAVE_GLPK
        return OPT_OPERATION_NOT_IMPLEMENTED;
    #endif
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Glpk::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_GLPK
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_GLPK
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_GLPK
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setVariableType( const int index, const OPT_VARTYPE varType )
#if OPT_HAVE_GLPK
{
    glp_set_col_kind( prob, index + 1, varType == OPT_VT_INTEGER ? GLP_IV : GLP_CV );

    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_GLPK
{
    int r, status;
    const int n  =  glp_get_num_cols(prob);
    const int nI = glp_get_num_int(prob);
    const int m = glp_get_num_rows(prob);
    
    
    
    if(resetSol)
    {
        origSolverRetCode = GLP_UNDEF;
        this->resetSol();
    }
    else
    {
        feasSol = false;
    }

    if( nI > 0 )
    {
        r = glp_intopt(prob, intParam);
        
        //std::cout << "r: " << r << "\n";

        status = glp_mip_status(prob);

        objValue = glp_mip_obj_val(prob);
        
        if( storeSol )
        {
            for(int i = 0; i < n; i++)
                sol[i] = glp_mip_col_val(prob, i+1);
        }
        
        if( storeConstrs )
        {
            for(int i = 0; i < m; i++)
                constr[i] = glp_mip_row_val(prob, i+1);
        }

    }
    else
    {
        r = glp_simplex(prob, simParam);

        status = glp_get_status(prob);

        objValue = glp_get_obj_val(prob);
        
        if( storeSol )
        {
            for(int i = 0; i < n; i++)
                sol[i] = glp_get_col_prim( prob, i+1 );
        }
        
        
        if( storeConstrs )
        {
            for(int i = 0; i < m; i++)
                constr[i] = glp_get_row_prim(prob, i+1);
        }
        
        
        if( storeDualSol )
        {
            for(int i = 0; i < n; i++)
                dualSolV[i] = glp_get_col_dual(prob, i+1);
            
            for(int i = 0; i < m; i++)
                dualSolC[i] = glp_get_row_dual(prob, i+1);
        }
    }

    origSolverRetCode = r;

    switch(r)
    {
    case 0:
        break;

    case GLP_EITLIM:
        
        #if OPT_PRINT_MAX_ITER_WARNING
            if( numberOfWarningsByIterLimit < maxNumberOfWarningsByIterLimit )
            {
                std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Glpk solving!\n";
                numberOfWarningsByIterLimit++;
                
                if( numberOfWarningsByIterLimit == maxNumberOfWarningsByIterLimit )
                    std::cerr << OPT_PREPRINT "Warning: Maximum number of warnings by maximum iteration achieved! Stopping these warnings.\n";
            }
        #endif
        
        retCode = OPT_MAX_ITERATIONS;
        goto termination;
        break;

    case GLP_ETMLIM:
        retCode = OPT_MAX_TIME;
        goto termination;
        break;

    case GLP_ENOPFS:
        retCode = OPT_INFEASIBLE_PROBLEM;
        goto termination;
        break;

    default:
        retCode = OPT_UNDEFINED_ERROR;
        goto termination;
        break;
    }
    
    origSolverRetCode = status; //so function return 0, we store in origSolverRetCode the status. Otherwise, we let origSolverRetCode with r
    
    switch( status )
    {
    case GLP_OPT:

        feasSol = true;
        retCode = OPT_OPTIMAL_SOLUTION;
        break;

    case GLP_FEAS:

        feasSol = true;
        retCode = OPT_UNDEFINED_ERROR;
        break;

    case GLP_NOFEAS:

        retCode = OPT_INFEASIBLE_PROBLEM;
        break;

    case GLP_UNBND:

        feasSol = true; //I am sue about that, but I supposethe solution is feasible...
        retCode = OPT_UNBOUNDED_PROBLEM;
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



int OPT_Glpk::warmUp()
#if OPT_HAVE_GLPK
{
    return glp_warm_up(prob) == 0 ? 0 : OPT_SOLVER_ERROR;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


// __ methods from LPSolver __



int OPT_Glpk::__addConstraints(const int nm)
#if OPT_HAVE_GLPK
{
    glp_add_rows( prob, nm );

    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::__addVariables(const int nn, const bool initFree )
#if OPT_HAVE_GLPK
{
    glp_add_cols(prob, nn);

    if( initFree )
    {
        const int n = glp_get_num_cols(prob);

        for(int i = n; i > n-nn; i-- )
            glp_set_col_bnds( prob, i, GLP_FR, 0.0, 0.0 );
    }

    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


/*
int OPT_Glpk::allocateAuxStructures(const int n)
{
    int *auxi;
    double *auxd;


    auxd = (double *) realloc( auxValues2, n * sizeof(double) );
    if( !auxd )
        return OPT_MEMORY_ERROR;

    auxValues2 = auxd;


    auxi = (int *) realloc( auxIndex2, n * sizeof(int) );
    if( !auxi )
        return OPT_MEMORY_ERROR;

    auxIndex2 = auxi;


    return OPT_LPSolver::allocateAuxStructures(n);
}
*/





int OPT_Glpk::generateModelFile(const char* fileName)
#if OPT_HAVE_GLPK
{
    glp_set_prob_name(prob, "lpGlpkOptSolvers");

    int r = glp_write_lp( prob, NULL, fileName );

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



int OPT_Glpk::getConstraintBounds( const int index, double &lb, double &ub )
#if OPT_HAVE_GLPK
{
    lb = OPT_max(glp_get_row_lb(prob, index + 1), -OPT_INFINITY);
    ub = OPT_min(glp_get_row_ub(prob, index + 1), OPT_INFINITY);
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value)
#if OPT_HAVE_GLPK
{
    value = 0.0;

    const int c = varIndex + 1;

    int nz = glp_get_mat_row(prob, constrIndex+1, auxIndex, auxValues );


    for( int i = 1; i <= nz; i++ )
    {
        if( auxIndex[i] == c )
        {
            value = auxValues[i];
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



int OPT_Glpk::getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values)
#if OPT_HAVE_GLPK
{
    nzs = glp_get_mat_row(prob, constrIndex+1, auxIndex, auxValues);
    
    for(int i = 0; i < nzs; i++)
        cols[i] = auxIndex[i+1] - 1;
    
    for(int i = 0; i < nzs; i++)
        values[i] = auxValues[i+1];
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Glpk::getObjConstant(double &objConstant)
#if OPT_HAVE_GLPK
{
    objConstant = glp_get_obj_coef(prob, 0);
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Glpk::getObjLinearCoef( const int index, double &value )
#if OPT_HAVE_GLPK
{
    value = glp_get_obj_coef(prob, index + 1 );
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::getNumberOfConstraints(int &m)
#if OPT_HAVE_GLPK
{
    m = glp_get_num_rows(prob);
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif

        
int OPT_Glpk::getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs)
#if OPT_HAVE_GLPK
{
    nzs = glp_get_mat_row(prob, constrIndex, NULL, NULL);
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Glpk::getNumberOfIntVars(int &n)
#if OPT_HAVE_GLPK
{
    n = glp_get_num_int(prob);
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::getNumberOfVars(int &n)
#if OPT_HAVE_GLPK
{
    n = glp_get_num_cols(prob);
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::getObjSense(OPT_OPTSENSE &sense)
#if OPT_HAVE_GLPK
{
    sense = glp_get_obj_dir(prob) == GLP_MIN ? OPT_MINIMIZE : OPT_MAXIMIZE;

    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::getVariableBounds(const int index, double &lb, double &ub)
#if OPT_HAVE_GLPK
{
    lb = OPT_max(glp_get_col_lb(prob, index+1), -OPT_INFINITY);
    
    ub = OPT_min(glp_get_col_ub(prob, index+1), OPT_INFINITY);
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::removeConstraints(const int ninds, const int *indices )
#if OPT_HAVE_GLPK
{
    int code;
    int *rows = NULL;


    //if( ninds < naux ) //note we need one more index to glpk, so do not use <=
    {
        rows = auxIndex;
    }
    /*else
    {
        rows = (int *) malloc( (ninds + 1) * sizeof(int) );

        if( !rows )
        {
            code = OPT_MEMORY_ERROR;
            goto termination;
        }
    } */

    OPT_copyShiftArray( ninds, 1, indices, &rows[1] );

    glp_del_rows( prob, ninds, rows );

    code = 0;

//termination:


    //if( rows != auxIndex && rows != NULL )
    //free(rows);

    return code;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::removeVars(const int ninds, const int *indices)
#if OPT_HAVE_GLPK
{
    OPT_copyShiftArray(ninds, 1, indices, &auxIndex[1] );

    glp_del_cols(prob, ninds, auxIndex);

    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setLinearColumn( const int varIndex, const int nzs, const int* rows, const double* values)
#if OPT_HAVE_GLPK
{
    const int m = glp_get_num_rows(prob);
    int nz, code;



    //int *auxIndex3 = NULL;
    //double *auxValues3 = NULL;


    //if( m < naux ) //note we need one more index to glpk, so do not use <=
    {
        //auxIndex = this->auxIndex;
        //auxValues = this->auxValues;

        //auxIndex3 = this->auxIndex3;
        //auxValues3 = this->auxValues3;
    }
    /* else
    {
        auxIndex = (int *) malloc( 2*(m + 1) * sizeof(int) );
        auxVals = (double *) malloc( 2*(m + 1) * sizeof(double) );

        if( !auxIndex || !auxVals )
        {
            code = OPT_MEMORY_ERROR;
            goto termination;
        }

        auxIndex2 = &auxIndex[m + 1];
        auxVals2 = &auxVals[m + 1];
    } */


    nz = glp_get_mat_col( prob, varIndex+1, auxIndex, auxValues );

    if( nz == 0 )
    {
        OPT_copyArray(nzs, values, &auxValues[1]);

        OPT_copyShiftArray(nzs, 1, rows, &auxIndex[1]);

        glp_set_mat_col(prob, varIndex+1, nzs, auxIndex, auxValues);
    }
    else
    {
        OPT_setAllArray(m, &auxValues2[1], 0.0);

        for(int i = 1; i <= nz; i++)
            auxValues2[ auxIndex[i] ] = auxValues[i];

        for(int i = 0; i < nzs; i++)
            auxValues2[ rows[i]+1 ] = values[i];

        for(int i = 1; i <= m; i++)
            auxIndex2[i] = i;


        glp_set_mat_col(prob, varIndex+1, m, auxIndex2, auxValues2);
    }



    code = 0;

    /* termination:


        if( m >= naux )
        {
            if( auxIndex )	free(auxIndex);
            if( auxVals )	free(auxVals);
        } */

    return code;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::boundFlag( const double lb, const double ub )
#if OPT_HAVE_GLPK
{
    int f;

    if( lb > -OPT_INFINITY )
    {
        if( ub < OPT_INFINITY )
            f = lb == ub ? GLP_FX : GLP_DB;
        else
            f = GLP_LO;
    }
    else
    {
        f = ub < OPT_INFINITY ? GLP_UP : GLP_FR;
    }


    return f;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::resetConstraintLinearPart( const int constrIndex, const int nzs, const int *cols, const double *values )
#if OPT_HAVE_GLPK
{
    //int r = setLinearConstraintCoefs( constrIndex, nzs, cols, values );
    //if( r != 0 )
        //return r;
    
    OPT_copyShiftArray( nzs, 1, cols, &auxIndex[1] );
    
    OPT_copyArray( nzs, values, &auxValues[1]  ); //we have to start form index 1 :(
    
    glp_set_mat_row(prob, constrIndex+1, nzs, auxIndex, auxValues);
    

    return 0; //setConstraintBounds(constrIndex, lb, ub);
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setConstraintBounds( const int index, const double lb, const double ub )
#if OPT_HAVE_GLPK
{
    glp_set_row_bnds( prob, index+1, boundFlag(lb, ub), lb, ub );

    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



//that method overload all prealocated elements...
int OPT_Glpk::setConstraintsLinearCoefs( const int nzs, const int* rows, const int* cols, const double* values )
#if OPT_HAVE_GLPK
{
    const int n = glp_get_num_cols(prob);
    const int m = glp_get_num_rows(prob);


    for(int i = 1; i <= n; i++)
        auxIndex[i] = i;


    for(int i = 0; i < m; i++)
    {
        //getting the actual coefficients...
        int nz = glp_get_mat_row(prob, i+1, auxIndex2, auxValues2);

        if( nz == 0 )
        {
            nz = 0;
            for(int k = 0; k < nzs; k++)
            {
                if( rows[k] == i )
                {
                    nz++;
                    auxIndex2[nz] = cols[k]+1;
                    auxValues2[nz] = values[k];
                }
            }


            glp_set_mat_row(prob, i+1, nz, auxIndex2, auxValues2);
        }
        else
        {
            OPT_setAllArray( n, &auxValues[1], 0.0 );

            for(int j = 1; j <= nz; j++)
                auxValues[ auxIndex2[j] ] = auxValues2[j];


            for(int k = 0; k < nzs; k++)
            {
                if( rows[k] == i )
                    auxValues[cols[k]+1] = values[k];
            }


            glp_set_mat_row(prob, i+1, n, auxIndex, auxValues);
        }
    }


    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setConstraintLinearCoefs( const int constrIndex, const int nzs, const int *cols, const double *values)
#if OPT_HAVE_GLPK
{
    int pnz;


    //unfortunatelly glpk replaces the entire row when we call glp_set_mat_row. So, we will get the current values in the respective row to do not lost them...

    pnz = glp_get_mat_row( prob, constrIndex+1, auxIndex2, auxValues2 );


    if( pnz == 0 )
    {
        OPT_copyShiftArray( nzs, 1, cols, &auxIndex[1] );

        OPT_copyArray( nzs, values, &auxValues[1]  ); //we have to star form index 1 :(

        glp_set_mat_row(prob, constrIndex+1, nzs, auxIndex, auxValues );
    }
    else
    {
        const int n = glp_get_num_cols(prob);

#if OPT_DEBUG_MODE
        assert(naux > n);
#endif


        OPT_setAllArray( n, &auxValues[1], 0.0 );

        for(int i = 1; i <= pnz; i++)
            auxValues[ auxIndex2[i] ] = auxValues2[i];

        //setting now the new values. Note new values can overwrite the olders.
        for(int i = 0; i < nzs; i++)
            auxValues[ cols[i]+1 ] = values[i];


        for(int i = 1; i <= n; i++)
            auxIndex[i] = i;


        //note, this functions ignores zeros... :)
        glp_set_mat_row(prob, constrIndex+1, n, auxIndex, auxValues );
    }




    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value)
#if OPT_HAVE_GLPK
{
    //too many work to reply all here. In this case, it is better call our API function...

    return setConstraintsLinearCoefs(1, &constrIndex, &varIndex, &value);
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Glpk::setObjLinearCoef( const int index, const double value )
#if OPT_HAVE_GLPK
{
    glp_set_obj_coef(prob, index+1, value );

    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setObjLinearCoefs( const int nzs, const int* cols, const double* values )
#if OPT_HAVE_GLPK
{
    for(int i = 0; i < nzs; i++)
    {
        glp_set_obj_coef(prob, cols[i] + 1, values[i] );
    }

    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Glpk::setObjLinearPart( const int n, const double *values )
#if OPT_HAVE_GLPK
{
    const int on = glp_get_num_cols(prob);
    
    for(int i = 1; i <= n; i++)
    {
        glp_set_obj_coef( prob, i, values[i-1] );
    }
    
    for(int i = n+1; i <= on; i++)
    {
        glp_set_obj_coef( prob, i, 0.0 );
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




void OPT_Glpk::setObjConstant(const double value)
{
    #if OPT_HAVE_GLPK
        glp_set_obj_coef(prob, 0, value);
    #endif
}




int OPT_Glpk::setObjSense( const OPT_OPTSENSE sense )
#if OPT_HAVE_GLPK
{
    glp_set_obj_dir( prob, sense == OPT_MINIMIZE ? GLP_MIN : GLP_MAX);

    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Glpk::setVariableBounds( const int index, const double lb, const double ub )
#if OPT_HAVE_GLPK
{
    glp_set_col_bnds( prob, index+1, boundFlag(lb, ub), lb, ub );

    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif








