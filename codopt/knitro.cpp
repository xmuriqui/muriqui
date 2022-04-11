

#include <cmath>
#include <iostream>


#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"


#if OPT_HAVE_KNITRO
    #include "knitro.h"
#endif



using namespace optsolvers;
//using namespace newspm;
using namespace minlpproblem;




#if OPT_HAVE_KNITRO

int optsolvers::OPT_knitroObjConstrF (const int             evalRequestCode, 
                    const int             n,
                    const int             m,
                    const int             nnzJ,
                    const int             nnzH,
                    const double * const  x,
                    const double * const  lambda,
                        double * const  obj,
                        double * const  c,
                        double * const  objGrad,
                        double * const  jac,
                        double * const  hessian,
                        double * const  hessVector,
                        void   *        userParams)
{
    OPT_Knitro *myknitro = (OPT_Knitro *) userParams;
    OPT_MINLPProb &prob = myknitro->prob;
    
    const bool hasFreeConstrs = myknitro->hasFreeConstrs;
    const unsigned int thnumber = myknitro->threadNumber;
    const bool *auxCEval = myknitro->auxCEval;
    double *g = myknitro->auxValues;
    double *pconstr = hasFreeConstrs ? g : c;
    
    bool newx = true;
    int r;
    
    
    /*#if OPT_DEBUG_MODE
        if(evalRequestCode != KTR_RC_EVALFC)
        {
            assert(false);
        }
    #endif*/
    
    
    
    r = prob.objEval(thnumber, newx, x, *obj, myknitro->in_nl_obj_factor);
    if( r != 0 )
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            std::cerr <<  OPT_PREPRINT "Callback function error " << r << " at objective evaluation"  OPT_GETFILELINE "\n";
        #endif
            
        return r;
    }
    
    
    if( prob.hasNlObj )
        newx = false;
    
    
    {
        
        
        r = prob.constraintsEval(thnumber, newx, auxCEval, x, pconstr);
        if( r != 0 )
        {
            #if OPT_PRINT_CALLBACK_ERROR_MSG
                std::cerr <<  OPT_PREPRINT "Callback function error " << r << " at nonlinear constraints evaluation"  OPT_GETFILELINE "\n";
            #endif
                
            return r;
        }
        
        
        if(hasFreeConstrs)
        {
            const int om = prob.m;
            int k = 0;
            
            for(int i = 0; i < om; i++)
            {
                if( auxCEval[i] )
                {
                    c[k] = pconstr[i]; 
                    k++;
                }
            }
            #if OPT_DEBUG_MODE
                assert(k == m);
            #endif
        }
        
    }
    
    
    
    return 0;
}




int  optsolvers::OPT_knitroGrads (const int             evalRequestCode,
                    const int             n,
                    const int             m,
                    const int             nnzJ,
                    const int             nnzH,
                    const double * const  x,
                    const double * const  lambda,
                            double * const  obj,
                            double * const  c,
                            double * const  objGrad,
                            double * const  jac,
                            double * const  hessian,
                            double * const  hessVector,
                            void   *        userParams)
{
    OPT_Knitro *myknitro = (OPT_Knitro *) userParams;
    OPT_MINLPProb &prob = myknitro->prob;
    
    const unsigned int thnumber = myknitro->threadNumber;
    const bool *auxCEval = myknitro->auxCEval;
    bool newx = true;
    int r;
    
    
    /*#if OPT_DEBUG_MODE
        if (evalRequestCode != KTR_RC_EVALGA)
        {
            assert(false);
        }
    #endif*/
    
    
    {
        //const bool *nlConstr = prob.nlConstr;
        
        MIP_SparseMatrix &J = prob.J;
        
        
        int nzs;
        const bool hasFreeConstrs = myknitro->hasFreeConstrs;
        const int mquad = myknitro->mquad;
        //const int *quadIndex = myknitro->quadIndex;
        const OPT_SparseMatrix *singleJ = OPT_getSingleJacobianPointer(prob, mquad);
        
        const int *sizeColsNzRowJac = myknitro->sizeColsNzRowJac;
        int **colsNzRowJac = myknitro->colsNzRowJac;
        
        
        //OPT_setAllArray( nnzJ, jac, 0.0 );
        
        if( prob.hasNlConstrs )
        {
            int r = prob.nlJacobianEval( thnumber, newx, auxCEval, x, J);
            
            newx = false;
            
            if( r != 0 )
            {
                #if OPT_PRINT_CALLBACK_ERROR_MSG
                    std::cerr <<  OPT_PREPRINT  "Callback function error " << r << " at nonlinear jaobian evaluation"  OPT_GETFILELINE "\n";
                #endif
                return r;
            }
            
        }
        
        
        if(singleJ && !hasFreeConstrs)
        {
            nzs = singleJ->getValues(jac);
        }
        else
        {
            nzs = OPT_getValuesFromCompleteJacobian( prob, J, sizeColsNzRowJac, colsNzRowJac, auxCEval, m, x, objGrad, jac);
        }
        
        #if OPT_DEBUG_MODE
            assert(nzs == nnzJ);
        #endif
        
    }
    
    
    //objective function. We eval after constraints because we use objGrad as temp array...
    r = prob.objGradEval( thnumber, newx, x, objGrad, myknitro->in_nl_obj_factor );
    if( r != 0 )
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            std::cerr <<  OPT_PREPRINT "Callback function error " << r << " at nonlinear constraints evaluation"  OPT_GETFILELINE "\n";
        #endif
            
        return r;
    }
    
    
    return 0;
}




int  optsolvers::OPT_knitroHess (const int             evalRequestCode,
                    const int             n,
                    const int             m,
                    const int             nnzJ,
                    const int             nnzH,
                    const double * const  x,
                    const double * const  lambda,
                            double * const  obj,
                            double * const  c,
                            double * const  objGrad,
                            double * const  jac,
                            double * const  hessian,
                            double * const  hessVector,
                            void   *        userParams)
{
    OPT_Knitro *myknitro = (OPT_Knitro *) userParams;
    OPT_MINLPProb &prob = myknitro->prob;
    
    const unsigned int thnumber = myknitro->threadNumber;
    const int mquad = myknitro->mquad;
    //const bool hasFreeNLConstrs = myknitro->hasFreeNLConstrs;
    
    const bool *auxCEval = myknitro->auxCEval;
    const int *quadIndex = myknitro->quadIndex;
    
    const double objFactor = evalRequestCode == KTR_RC_EVALH_NO_F ? 0.0 : 1.0  ; //prob->objFactor is already considered in prob->nlpHessianEval...
    
    OPT_SparseMatrix &lagH = prob.lagH;
    
    const int nNzRowsLagH = myknitro->nNzRowsLagH;
    const int *nzRowsLagH = myknitro->nzRowsLagH;
    const int *sizeColsNzLagH = myknitro->sizeColsNzLagH;
    int **colsNzRowLagH = myknitro->colsNzRowLagH;
    
    double *auxValues = myknitro->auxValues;
    double *auxValues2= myknitro->auxValues2;
    
    
    //delet this some day...
    #if OPT_DEBUG_MODE
        assert( evalRequestCode != KTR_RC_EVALHV_NO_F && evalRequestCode != KTR_RC_EVALHV );
    #endif
    
    int nzs;
    
    
    int r = OPT_evalCompleteLagrangianHessian(thnumber, true, x, prob, lagH, mquad, quadIndex, objFactor,  myknitro->in_nl_obj_factor, lambda, auxCEval, m, nNzRowsLagH, nzRowsLagH, sizeColsNzLagH, colsNzRowLagH, auxValues, auxValues2, nzs, hessian );
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return 0;
    }
    
    
    return 0;
}

#endif







OPT_Knitro::OPT_Knitro():OPT_MyNLPSolver()
{
    initialize();
}


OPT_Knitro::~OPT_Knitro()
{
    deallocateSolverEnv();
    deallocateMemory();
}



// __methods from Solver __



void OPT_Knitro::deallocateSolverEnv()
{
    #if OPT_HAVE_KNITRO
        KTR_context** kc = (KTR_context**) &knitroContext;
        
        if(*kc)
        {
            KTR_free(kc);
            *kc = NULL;
        }
        
        problemInitialized = false;
    #endif
    
    desallocateMyAuxDerivativeIndexStructures();
    
    OPT_MyNLPSolver::deallocateSolverEnv();
    
    //printf("__desallocateSolverEnv");
    //OPT_getchar();
}


bool OPT_Knitro::getMinusLambdaOnLagran()
{
    return false;
}


/*int OPT_Knitro::allocateConstrStructures(const int m)
{
    
    
    if( m > maux )
    {
    
        bool *auxb = (bool *) realloc( allCEvalTrue, m*sizeof(bool) ); 
        
        if( !auxb )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            return OPT_MEMORY_ERROR;
        }
        
        allCEvalTrue = auxb;
        
        OPT_setAllArray( m-maux, &(allCEvalTrue[maux]), true );
        
    }
    
    
    //we are calling allocateConstrStructures at end because we use variable maux 
    const int r = OPT_MyNLPSolver::allocateConstrStructures(m);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    
    return 0;
}
*/


int OPT_Knitro::allocateConstrStructures(const int m)
{
    if( m > maux )
    {
        /*int *pi = (int *) realloc( cFnType, m*sizeof(int) ); 
        
        if( !pi )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            return OPT_MEMORY_ERROR;
        }
        
        cFnType = pi; */
        
        int r = OPT_realloc( cFnType, m );
        OPT_IFERRORRETURN(r, r);
        
        
        #if OPT_HAVE_KNITRO
            OPT_setAllArray( m-maux, &(cFnType[maux]), KTR_FNTYPE_UNCERTAIN );
        #endif
    }
    
    
    //we are calling allocateConstrStructures at end because we use variable maux 
    const int r = OPT_MyNLPSolver::allocateConstrStructures(m);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    
    return 0;
}


void OPT_Knitro::deallocateMemory()
{
    OPT_MyNLPSolver::deallocateMemory();
    
    OPT_secFree(auxLambda);
    OPT_secFree(quadIndex);
    OPT_secFree(nqcons);
    OPT_secFree(cFnType);
    //OPT_secFree(allCEvalTrue);
    
    desallocateMyAuxDerivativeIndexStructures();
}


void OPT_Knitro::desallocateMyAuxDerivativeIndexStructures()
{
    if(colsNzRowJac)
    {
        
        for(int i = 0; i < nRowsJac; i++)
        {
            if(colsNzRowJac[i])
                free(colsNzRowJac[i]);
        }
        
        free(colsNzRowJac);
        colsNzRowJac = NULL;
    }
    
    OPT_secFree(sizeColsNzRowJac);
    
    if(colsNzRowLagH)
    {
        for(int i = 0; i < nNzRowsLagH; i++)
        {
            if(colsNzRowLagH[i])
                free(colsNzRowLagH[i]);
        }
        free(colsNzRowLagH);
        colsNzRowLagH = NULL;
        nRowsJac = 0;
    }
    
    OPT_secFree(sizeColsNzLagH);
    OPT_secFree(nzRowsLagH);
    
    nNzRowsLagH = 0;
}


void OPT_Knitro::initialize()
{
    OPT_MyNLPSolver::initialize();
    
    problemInitialized = false;
    
    consBoundChg = false;
    newInitPoint = false;
    paramChg = false;
    
    mquad = -1;
    quadIndex = NULL;
    
    nqcons = NULL;
    indqcons = NULL;
    
    auxLambda = NULL;
    
#if OPT_HAVE_KNITRO
    objFnType = KTR_FNTYPE_UNCERTAIN;
#endif
    cFnType = NULL;
    
    //allCEvalTrue = NULL;
    knitroContext = NULL;
    
    newm = -1;
    
    nRowsJac = 0;
    sizeColsNzRowJac = NULL;
    colsNzRowJac = NULL;
    nNzRowsLagH = 0;
    nzRowsLagH = NULL;
    sizeColsNzLagH = NULL;
    colsNzRowLagH = NULL;
}



int OPT_Knitro::getNumberOfIterations( long unsigned int &niter)
#if OPT_HAVE_KNITRO
{
    int nI;
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    getNumberOfIntVars(nI);
    
    
    niter = nI == 0 ? KTR_get_number_iters(kc) : KTR_get_mip_num_nodes(kc);
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


OPT_LISTSOLVERS OPT_Knitro::getSolverCode()
{
    return optsolvers::OPT_KNITRO;
}



int OPT_Knitro::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
#if OPT_HAVE_KNITRO
{
    int r;
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    
    
    //we call it here because maybe user need call this method because we lost knitro parameters when we have to reset KTR_init_problem in OPT_Knitro::solve()
    deallocateSolverEnv();
    
    
    kc = KTR_new();
    if( kc == NULL )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORMSG("Error at knitro KTR_new. Try to check your knitro license.");
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    //setting
    KTR_set_func_callback (kc, &OPT_knitroObjConstrF);
    KTR_set_grad_callback (kc, &OPT_knitroGrads);
    KTR_set_hess_callback (kc, &OPT_knitroHess);
    
    
    //we does not allow knitro perform callback evaluations in parallel by now... That is because maybe minlpProb callbacks can do not hope evaluations in parallel are made, or maybe it is prepared for a smaller number of threads...
    
    r = KTR_set_int_param (kc, KTR_PARAM_PAR_CONCURRENT_EVALS, KTR_PAR_CONCURRENT_EVALS_NO);
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
    #endif
    
    
    
    //---- SPECIFY THAT THE USER IS ABLE TO PROVIDE EVALUATIONS
    // *---- OF THE HESSIAN MATRIX WITHOUT THE OBJECTIVE COMPONENT.
    // *---- TURNED OFF BY DEFAULT BUT SHOULD BE ENABLED IF POSSIBLE.
    r = KTR_set_int_param (kc, KTR_PARAM_HESSIAN_NO_F, KTR_HESSIAN_NO_F_ALLOW);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
    #endif
    
    //to knitro honor bounds...
    r = KTR_set_int_param (kc, KTR_PARAM_HONORBNDS, KTR_HONORBNDS_ALWAYS);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
    #endif
    
    
    r = KTR_set_int_param(kc, KTR_PARAM_OUTLEV, 0);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
    #endif
    
    
    paramChg = true;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::setConstraintBounds( const int index, const double lb, const double ub )
{
    const int r = OPT_MyNLPSolver::setConstraintBounds( index, lb, ub);
    
    if( r == 0 )
    {
        consBoundChg = true;
        return 0;
    }
    
    #if OPT_DEBUG_MODE
        OPT_PRINTERRORNUMBER(r);
    #endif
    
    return r;
}



int OPT_Knitro::setInitialSolution(const double* x, const double* dualConstrs, const double* dualVars)
{
    const int r = OPT_MyNLPSolver::setInitialSolution(x, dualConstrs, dualVars);
    
    if( r == 0 )
    {
        newInitPoint = true;
        return 0;
    }
    
    #if OPT_DEBUG_MODE
        OPT_PRINTERRORNUMBER(r);
    #endif
    
    return r;
}



int OPT_Knitro::setMaxCPUTime(const double time)
#if OPT_HAVE_KNITRO
{
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_double_param(kc, KTR_PARAM_MAXTIMECPU, time);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    paramChg = true;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::setObjCutLowerBound(const double objLBound)
#if OPT_HAVE_KNITRO
{
    //knitro just have a parameter KTR_PARAM_FSTOPVAL to stop if a feasible solution is better than a threshould. Unfortunatelly, that is not enough to us...
    paramChg = true;
    
    return OPT_OPERATION_NOT_IMPLEMENTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Knitro::setObjCutUpperBound(const double objUBound)
#if OPT_HAVE_KNITRO
{
    //knitro just have a parameter KTR_PARAM_FSTOPVAL to stop if a feasible solution is better than a threshould. Unfortunatelly, that is not enough to us...
    paramChg = true;
    
    return OPT_OPERATION_NOT_IMPLEMENTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::setMaxTime(const double time)
#if OPT_HAVE_KNITRO
{
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_double_param(kc, KTR_PARAM_MAXTIMEREAL, time);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    paramChg = true;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::setNumberOfThreads(const int nthreads)
#if OPT_HAVE_KNITRO
{
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_int_param(kc, KTR_PARAM_PAR_NUMTHREADS, nthreads);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    paramChg = true;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::setOutputLevel(const int level)
#if OPT_HAVE_KNITRO
{
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_int_param(kc, KTR_PARAM_OUTLEV, level);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    paramChg = true;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::setRelativeDualTol( const double tol )
#if OPT_HAVE_KNITRO
{
    paramChg = true;
    return OPT_OPERATION_NOT_IMPLEMENTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_KNITRO
{
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_double_param(kc, KTR_PARAM_OPTTOL, tol);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    paramChg = true;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::setRelativePrimalTol( const double tol )
#if OPT_HAVE_KNITRO
{
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_double_param(kc, KTR_PARAM_FEASTOL, tol);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    paramChg = true;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Knitro::_setParameters(const OPT_GeneralSolverParams &params)
#if OPT_HAVE_KNITRO
{
    int r, code = 0;


    for(auto &pair : params.intParams)
    {
        r = _setIntegerParameter(pair.first.c_str(),  pair.second);
        if ( r != 0 )
        {
            code = OPT_BAD_INPUT;
        }
    }

    for(auto &pair : params.dblParams)
    {
        r = _setDoubleParameter(pair.first.c_str(),  pair.second);
        if ( r != 0 )
        {
            code = OPT_BAD_INPUT;
        }
    }

    for(auto &pair : params.strParams)
    {
        r = _setStringParameter(pair.first.c_str(),  pair.second.c_str());

        if ( r != 0 )
        {
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


int OPT_Knitro::_setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_KNITRO
{
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_double_param_by_name(kc,  param, value);
    
    if( r != 0 )
    {
        //#if OPT_DEBUG_MODE
            //OPT_PRINTERRORNUMBER(r);
        //#endif
        
        printDblParamErrorMsg(!r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    paramChg = true;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::_setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_KNITRO
{
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_int_param_by_name(kc,  param, value);
    
    if( r != 0 )
    {
        //#if OPT_DEBUG_MODE
            //OPT_PRINTERRORNUMBER(r);
        //#endif
        
        printIntParamErrorMsg(!r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    paramChg = true;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::_setStringParameter(const char *param, const char *value)
#if OPT_HAVE_KNITRO
{
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_char_param_by_name(kc,  param, value);
    
    if( r != 0 )
    {
        //#if OPT_DEBUG_MODE
            //OPT_PRINTERRORNUMBER(r);
        //#endif
        
        printStrParamErrorMsg(r, param, value);
        
        return OPT_BAD_INPUT;
    }
        
    
    paramChg = true;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_KNITRO
{
    const int r = params.storeDoubleParameter(param, value);
        
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
    }
    
    paramChg = true;
    
    return r;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_KNITRO
{
    const int r = params.storeIntegerParameter(param, value);
        
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
    }
    
    paramChg = true;
    
    return r;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Knitro::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_KNITRO
{
    const int r = params.storeStringParameter(param, value);
        
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
    }
    
    paramChg = true;
    
    return r;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif






int OPT_Knitro::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_KNITRO
{
    const int n = prob.n;
    const int om = prob.m;
    const int nI = prob.getNumberOfIntegerVars();
    
    int r;
    
    
    const double *lx = prob.lx, *ux = prob.ux;
    const double *olc = prob.lc, *ouc = prob.uc;
    const MIP_SparseMatrix *QC = prob.QC;
    
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    
    int *jacIndexCons = NULL, *jacIndexVars;
    int *myhessRows = NULL, *myhessCols;
    double *mylambdaInit = NULL;
    
    
    
    
    //std::cout << "1 nmChg: " << nmChg << " genConstrChg: " << genConstrChg << " genHessChg: " << genHessChg << " consBoundChg: " << consBoundChg << "\n";
    
    
    if(resetSol)
    {
        origSolverRetCode = KTR_RC_NULL_POINTER; //I think we have to put some knitro code here...
        this->resetSol();
    }
    else
    {
        feasSol = false;
    }
    
    retCode = OPT_UNDEFINED;
    
    
    //mylambdaInit = (double*) calloc( (om+n), sizeof(double) );
    OPT_calloc(mylambdaInit, om+n);
    OPT_IFMEMERRORGOTOLABEL(!mylambdaInit, retCode, termination);
                    
    
    
    
    if(nmChg || genConstrChg || genHessChg || consBoundChg || varTypeChg)
    {
        int nnzJac;
        int nnzLagH;
        int objType;
        
        int *cType = auxIndex;
        int *xType = auxIndex2;
        
        double *lc = auxValues;//prob.lc;
        double *uc = auxValues2;//prob.uc;
        
        const bool *nlConstr = prob.nlConstr;
        double *plambdaInit = NULL;
        
        
        
        OPT_secFree(quadIndex);
        
        
        {
            newm = 0;
            mquad = 0;
            hasFreeConstrs = false;
            hasFreeNLConstrs = false;
            
            for(int i = 0; i < om; i++)
            {
                if(olc[i] > -MIP_INFINITY || ouc[i] < MIP_INFINITY)
                {
                    auxCEval[i] = true;
                    newm++;
                    
                    if(QC[i].getNumberOfElements() > 0)
                        mquad++;
                }
                else //free constraint
                {
                    hasFreeConstrs = true;
                    auxCEval[i] = false;
                    
                    if( nlConstr[i] )
                        hasFreeNLConstrs = true;
                }
            }
        }
        
        
        OPT_malloc(quadIndex, mquad);
        OPT_IFMEMERRORGOTOLABEL( !quadIndex, retCode, termination );
        
        
        {
            int k = 0;
            
            for(int i = 0; i < om; i++)
            {
                if(auxCEval[i] && QC[i].getNumberOfElements() > 0)
                {
                    quadIndex[k] = i;
                    k++;
                }
            }
            
            #if OPT_DEBUG_MODE
                assert(k == mquad);
            #endif
        }
        
        
        
        if( problemInitialized )
        {
            initSolverEnv(-1, -1, -1); //initSolverEnv already free memory from knitro
        }
        
        
        desallocateMyAuxDerivativeIndexStructures();
        
        //calculating number of nonzeros in jacobian
        {
            int nzs;
            const OPT_SparseMatrix *singleJ;
            
            int r = OPT_getSingleJacPointerOrCalcJacIndices(prob, mquad, hasFreeConstrs, auxCEval, singleJ, nnzJac, sizeColsNzRowJac, colsNzRowJac);
            
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            
            nRowsJac = prob.m; //we have to save the current number of constraints because user can add or delete constraints, and so, we would have a trouble to desallocate colsNzRowJac
            
            //jacIndexCons = (int *) malloc( nnzJac*2* sizeof(int) );
            
            OPT_malloc( jacIndexCons, nnzJac*2 );
            OPT_IFMEMERRORGOTOLABEL(!jacIndexCons, retCode, termination);
            
            
            jacIndexVars = &jacIndexCons[nnzJac];
            
            
            if(singleJ)
            {
                nzs = singleJ->getStructure(jacIndexCons, jacIndexVars);
            }
            else
            {
                int k = 0;
                
                nzs = 0;
                for(int i = 0; i < om; i++)
                {
                    if( !auxCEval[i] )
                        continue;
                    
                    const int ncols = sizeColsNzRowJac[i];
                    
                    OPT_setAllArray(ncols, &jacIndexCons[nzs], k);
                    OPT_copyArray(ncols, colsNzRowJac[i], &jacIndexVars[nzs]);
                    
                    nzs += ncols;
                    k++;
                }
                
                #if OPT_DEBUG_MODE
                    assert(k == newm);
                #endif
            }
            
            #if OPT_DEBUG_MODE
                assert(nzs == nnzJac);
            #endif
            
        }
        
        
        
        
        //calculating number of nonzeros in hessian
        {
            int nzs;
            const OPT_SparseMatrix *singleH;
            
            int r = OPT_getSingleLagHPointerOrCalcLagHIndices(prob, mquad, quadIndex, auxCEval, singleH, nnzLagH, nNzRowsLagH, nzRowsLagH, sizeColsNzLagH, colsNzRowLagH);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            
            //myhessRows = (int *) malloc( 2*nnzLagH*sizeof(int) );
            OPT_malloc( myhessRows, 2*nnzLagH );
            OPT_IFMEMERRORGOTOLABEL(!myhessRows, retCode, termination);
            
            
            myhessCols = &myhessRows[nnzLagH];
            
            
            if(singleH)
            {
                nzs = singleH->getStructure(myhessRows, myhessCols);
            }
            else
            {
                nzs = 0;
                
                for(int i = 0; i < nNzRowsLagH; i++)
                {
                    const int ncols = sizeColsNzLagH[i];
                    
                    OPT_setAllArray(ncols, &myhessRows[nzs], nzRowsLagH[i]);
                    
                    OPT_copyArray(ncols, colsNzRowLagH[i], &myhessCols[nzs] );
                    
                    nzs += ncols;
                }
            }
            
            #if OPT_DEBUG_MODE
                assert(nzs == nnzLagH);
            #endif
        }
        
        
        
        
        if( prob.hasNlObj )
        {
            objType = KTR_OBJTYPE_GENERAL;
        }
        else
        {
            //objType = prob.getNumberOfObjQuadTerms() == 0 ? KTR_OBJTYPE_LINEAR : KTR_OBJTYPE_QUADRATIC ;
            
            if( prob.getNumberOfObjQuadTerms() > 0 )
            {
                objType = KTR_OBJTYPE_QUADRATIC;
            }
            else
            {
                objType = prob.hasLinCoefObj() ? KTR_OBJTYPE_LINEAR : KTR_OBJTYPE_CONSTANT ;
            }
            
        }
        
        
        {
            int k = 0;
            for(int i = 0; i < om; i++)
            {
                if(auxCEval[i])
                {
                    if( nlConstr[i] )
                        cType[k] = KTR_CONTYPE_GENERAL;
                    else
                    {
                        cType[k] = QC[i].getNumberOfElements() == 0 ? KTR_CONTYPE_LINEAR : KTR_CONTYPE_QUADRATIC ;
                    }
                    k++;
                }
            }
            
            #if OPT_DEBUG_MODE
                assert(k == newm);
            #endif
        }
        
        
            
            
        if( om == newm )
        {
            OPT_copyArray(om, olc, lc);
            OPT_copyArray(om, ouc, uc);
        }
        else
        {
            int k = 0;
            
            for(int i = 0; i < om; i++)
            {
                if( auxCEval[i])
                {
                    lc[k] = olc[i];
                    uc[k] = ouc[i];
                    k++;
                }
            }
            
            #if OPT_DEBUG_MODE
                assert(k == newm);
            #endif
        }
        
        
        
        if(newm > 0)
        {
            if( !std::isnan(lambdaInit[0]) )
            {
                if( om == newm )
                {
                    OPT_copyArray(om, (const double*) lambdaInit, mylambdaInit);
                }
                else
                {
                    int k = 0;
                    
                    for(int i = 0; i < om; i++)
                    {
                        if(auxCEval[i])
                        {
                            mylambdaInit[k] = lambdaInit[i];
                            k++;
                        }
                    }
                    
                    #if OPT_DEBUG_MODE
                        assert(k == newm);
                    #endif
                }
                    
                    
                plambdaInit = mylambdaInit;
            }
        }
        
        
        if( n > 0  )
        {
            
            if( std::isnan(zInit[0]) )
            {
                if( plambdaInit ) //if plambdaInit is not defined, we let in this way
                    OPT_setAllArray(n, &plambdaInit[newm], 0.0);
            }
            else
            {
                if(!plambdaInit)
                {
                    OPT_setAllArray(newm, mylambdaInit, 0.0);
                    plambdaInit = mylambdaInit;
                }
                
                OPT_copyArray(n, zInit, &plambdaInit[newm]);
            }
        }
        
        
        _setParameters(params);
        
        
        if( nI == 0 )
        {
            //unfortunatelly, hessian must be represented by upper triangle. That is not good to us, but I hope solve this problem just changing the myhessRows and myhessCols parameters order
            
            //std::cout << "kc: " << kc << " n: " << n << " m: " << m << " objType: " << objType << " nnzJ: " << nnzJ << " nnzH: " << nnzH;
            
            
            r = KTR_init_problem(kc, n, KTR_OBJGOAL_MINIMIZE, objType, lx, ux, newm, cType, lc, uc, nnzJac, jacIndexVars, jacIndexCons, nnzLagH, myhessCols, myhessRows, std::isnan(xInit[0]) ? NULL : xInit, plambdaInit );
            
        }
        else
        {
            int *myxtype = prob.xtype;
            
            for(int i = 0; i < n; i++)
                xType[i] = MIP_isIntegerType(myxtype[i]) ? KTR_VARTYPE_INTEGER : KTR_VARTYPE_CONTINUOUS ;
            
            r = KTR_mip_init_problem(kc, n, KTR_OBJGOAL_MINIMIZE, objType, objFnType, xType, lx, ux, newm, cType, cFnType, lc, uc, nnzJac, jacIndexVars, jacIndexCons, nnzLagH, myhessCols, myhessRows,  std::isnan(xInit[0]) ? NULL : xInit, mylambdaInit );
        }
        
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            origSolverRetCode = r;
            retCode = OPT_SOLVER_ERROR;
            goto termination;
        }
        
        problemInitialized = true;
        
        
        nmChg = false;
        genConstrChg = false;
        genHessChg = false;
        consBoundChg = false;
        varTypeChg = false;
    }
    else
    {
        r = KTR_chgvarbnds(kc, lx, ux);
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            retCode = OPT_SOLVER_ERROR;
            goto termination;
        }
        
        
        if( paramChg || ( newInitPoint && !std::isnan(xInit[0]) ) )
        {
            double *px, *plambda;
            
            if( newInitPoint && !std::isnan(xInit[0]) )
            {
                
                if( hasFreeConstrs )
                {
                    int k = 0;
                    for(int i = 0; i < om; i++)
                    {
                        if(auxCEval[i])
                        {
                            mylambdaInit[k] = lambdaInit[i];
                            k++;
                        }
                    }
                    
                    #if OPT_DEBUG_MODE
                        assert(k == newm);
                    #endif
                }
                else
                {
                    OPT_copyArray( newm, lambdaInit, mylambdaInit );
                }
                
                OPT_copyArray( n, zInit, &(mylambdaInit[newm]) );
                
                px = xInit;
                plambda = mylambdaInit;
            }
            else
            {
                px = NULL;
                plambda = NULL;
            }
            
            
            r = KTR_restart( kc, px, plambda );
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                //we do not finish the function if we fail here...
            }
            
            newInitPoint = false;
        }
        
        
        
    }
    
    
    
    
    if( nI == 0 )
    {
        origSolverRetCode = KTR_solve( kc, storeSol ? sol : auxValues, mylambdaInit, 0, &objValue, NULL, NULL, NULL, NULL, NULL, this ); //knitro does not return constraint values in c input argument
    }
    else
    {
        origSolverRetCode = KTR_mip_solve( kc, storeSol ? sol : auxValues, auxLambda, 0, &objValue, NULL, NULL, NULL, NULL, NULL, this );
    }
    
    switch(origSolverRetCode)
    {
        case KTR_RC_OPTIMAL:
        case KTR_RC_NEAR_OPT:
            
            retCode = OPT_OPTIMAL_SOLUTION;
            feasSol = true;
            break;
            
        case KTR_RC_FEAS_XTOL:
        case KTR_RC_FEAS_NO_IMPROVE:
        case KTR_RC_FEAS_FTOL:
            
            retCode = OPT_FEASIBLE_SOLUTION;
            feasSol = true;
            break;
            
        case KTR_RC_INFEASIBLE:
        case KTR_RC_INFEAS_XTOL:
        case KTR_RC_INFEAS_NO_IMPROVE:
        case KTR_RC_INFEAS_MULTISTART:
        case KTR_RC_INFEAS_VAR_BOUNDS:
        case KTR_RC_INFEAS_CON_BOUNDS:
            
            retCode = OPT_INFEASIBLE_PROBLEM;
            break;
            
        case KTR_RC_ITER_LIMIT_FEAS:
            //retCode = OPT_MAX_ITERATIONS;
            feasSol = true;
            //do not break here
        case KTR_RC_ITER_LIMIT_INFEAS:
            retCode = OPT_MAX_ITERATIONS;
            #if OPT_PRINT_MAX_ITER_WARNING
                if( numberOfWarningsByIterLimit < maxNumberOfWarningsByIterLimit )
                {
                    std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Knitro solving!\n";
                    numberOfWarningsByIterLimit++;
                    
                    if( numberOfWarningsByIterLimit == maxNumberOfWarningsByIterLimit )
                        std::cerr << OPT_PREPRINT "Warning: Maximum number of warnings by maximum iteration achieved! Stopping these warnings.\n";
                }
            #endif
            break;
        
        case KTR_RC_TIME_LIMIT_FEAS:
            feasSol = true;
            //do not break here
        case KTR_RC_TIME_LIMIT_INFEAS:
            retCode = OPT_MAX_TIME;
            break;
        
        case KTR_RC_UNBOUNDED:
            retCode = OPT_UNBOUNDED_PROBLEM;
            feasSol = true;
            break;
        
        case KTR_RC_FEVAL_LIMIT_FEAS:
        case KTR_RC_MIP_TERM_FEAS :
        case KTR_RC_MIP_SOLVE_LIMIT_FEAS:
        case KTR_RC_MIP_NODE_LIMIT_FEAS:
        case KTR_RC_MIP_EXH_FEAS:
            feasSol = true;
            //do not break here
        default:
            
            #if OPT_DEBUG_MODE
                std::cout << OPT_PREPRINT "Status " << origSolverRetCode << " on Knitro termination code at OPT_Knitro::solve\n";
                //getchar();
            #endif
            retCode = OPT_UNDEFINED_ERROR;
    }
    
    
    if( storeConstrs && feasSol )
    {
        r = prob.constraintsEval( threadNumber, true, auxCEval,  storeSol ? sol : auxValues, constr);
        
        if( r != 0 )
        {
            #if OPT_PRINT_CALLBACK_ERROR_MSG
                std::cerr << OPT_PREPRINT << "Callback function error " << r << OPT_GETFILELINE << "\n";
            #endif
        }
    }
    
    
    if( storeDualSol )
    {
        int k = 0;
        for(int i = 0; i < om; i++)
        {
            if( auxCEval[i] )
            {
                dualSolC[i] = mylambdaInit[k];
                k++;
            }
            else
            {
                dualSolC[i] = 0.0;
            }
        }
        
        #if OPT_DEBUG_MODE
            assert(k == newm);
        #endif
        
        //OPT_copyArray( newm, mylambdaInit, dualSolC );
        OPT_copyArray( n, &(mylambdaInit[newm]), dualSolV );
    }
    
    
    if( prob.objFactor < 0 )
    {
        objValue = -objValue;
    }
    
    
    
    
termination:
    
    if( jacIndexCons )	free(jacIndexCons);
    if( myhessRows )	free(myhessRows);
    if( mylambdaInit )	free(mylambdaInit);
    
    return retCode;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




void OPT_Knitro::setObjFnType(const int value)
{
    objFnType = value;
}


void OPT_Knitro::setCFnType(const int *values)
{
    OPT_copyArray(prob.m, values, cFnType );
}





















































