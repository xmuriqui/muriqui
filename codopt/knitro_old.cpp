

#include <math.h>
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
    
    const unsigned int thnumber = myknitro->threadNumber;
    const bool *auxCEval = myknitro->auxCEval;
    bool newx = true;
    int r;
    
    
    #if OPT_DEBUG_MODE
        if(evalRequestCode != KTR_RC_EVALFC)
        {
            assert(false);
        }
    #endif
    
    
    
    r = prob.objEval(thnumber, newx, x, *obj);
    if( r != 0 )
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return r;
    }
    
    
    if( prob.hasNlObj )
        newx = false;
    
    
    r = prob.constraintsEval(thnumber, newx, auxCEval, x, c);
    if( r != 0 )
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return r;
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
    
    
    #if OPT_DEBUG_MODE
        if (evalRequestCode != KTR_RC_EVALGA)
        {
            assert(false);
        }
    #endif
    
    
    //objective function. We eval after constraints because we use objGrad as temp array...
    r = prob.objGradEval( thnumber, newx, x, objGrad );
    if( r != 0 )
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            OPT_PRINTERRORNUMBER(r);
        #endif
            
        return r;
    }
    
    if( prob.hasNlObj )
        newx = false;
    
    
    
    {
        const bool *nlConstr = prob.nlConstr;
        
        int r;
        double *grad = myknitro->auxValues2;
        
        MIP_SparseMatrix &J = prob.J;
        const MIP_SparseMatrix &A = prob.A;
        const MIP_SparseMatrix *QC = prob.QC;
        
        unsigned int *jacRowStartIndex = myknitro->jacRowStartIndex ;
        int **jacCols = myknitro->jacCols;
        
        
        
        //OPT_setAllArray( nnzJ, jac, 0.0 );
        
        if( prob.hasNlConstrs )
        {
            r = prob.nlJacobianEval( thnumber, newx, auxCEval, x, J);
            
            if( r != 0 )
            {
                #if OPT_PRINT_CALLBACK_ERROR_MSG
                    std::cerr << OPT_PREPRINT << "Callback function error " << r << std::endl;
                #endif
                return r;
            }
        }
        
        for( int i = 0; i < m; i++ )
        {
            const unsigned int nz = jacRowStartIndex[i+1] - jacRowStartIndex[i];
            const int *jacCol = jacCols[i];
            
            
            double *v = &jac[ jacRowStartIndex[i] ];
            
            if( nlConstr[i] && QC[i].getNumberOfElements() == 0 && A.getNumberOfElementsAtRow(i) == 0 )
            {
                int *p = NULL;
                #if OPT_DEBUG_MODE
                    int rnz;
                    p = &rnz;
                #endif
                
                J.getRowValues(i, v, p);
                
                #if OPT_DEBUG_MODE
                    assert( nz == rnz );
                #endif
            }
            else if( QC[i].getNumberOfElements() == 0 && (!nlConstr[i] || J.getNumberOfElementsAtRow(i) == 0) )
            {
                int *p = NULL;
                #if OPT_DEBUG_MODE
                    int rnz;
                    p = &rnz;
                #endif
                
                A.getRowValues(i, v, p);
                
                #if OPT_DEBUG_MODE
                    assert( nz == rnz );
                #endif
            }
            else
            {
                OPT_setAllArray(nz, v, 0.0);
                MIP_constrCompleteGrad( prob, J,  i, x, grad, true);
                
                unsigned int k = 0;
                for( int j = 0; j < n; j++ )
                {
                    //columns are ordered in jacCols
                    if( grad[j] != 0.0 )
                    {
                        while( (int) jacCol[k] != j )
                            k++; //I hope jacCols has the non zeros index correctly in order. Otherwise, we will have serious problems here...
                        
                        v[k] = grad[j];
                        k++; //we can walk one position in k
                    }
                }
                
                #if OPT_DEBUG_MODE
                    assert( (int) k <= n );
                    assert( k <= nz );
                #endif
                
            }
            
        }
        
        
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
    const int *quadIndex = myknitro->quadIndex;
    
    const unsigned int *hessRowStartIndex = myknitro->hessRowStartIndex;
    const int *nqcons = myknitro->nqcons;
    const int *indqcons = myknitro->indqcons;
    
    bool newx = true;
    double *grad = myknitro->auxValues2;
    
    
    MIP_SparseMatrix &lagH = prob.lagH;
    const MIP_SparseMatrix &Q = prob.Q;
    const MIP_SparseMatrix *QC = prob.QC;
    
    
    const double objf = evalRequestCode == KTR_RC_EVALH_NO_F ? 0.0 : 1.0 ; //prob->objFactor is already considered in prob->nlpHessianEval...
    const bool objQ = Q.getNumberOfElements() > 0 && objf != 0.0;
    const bool constrQ = mquad > 0;
    const bool evalLag = ( prob.hasNlObj && objf != 0.0 ) || prob.hasNlConstrs;
    
    
    //delet this some day...
    #if OPT_DEBUG_MODE
        assert( evalRequestCode != KTR_RC_EVALHV_NO_F && evalRequestCode != KTR_RC_EVALHV );
    #endif
    
    
    if( evalLag )
    {
        //prob->objFactor is already considered in prob->nlpHessianEval...
        const int r = prob.nlpHessianEval(thnumber, newx, x, prob.hasNlObj ? objf : 0.0, lambda, lagH );
        
        newx = false;
            
        if( r != 0 )
        {
            #if OPT_PRINT_CALLBACK_ERROR_MSG
                std::cerr << OPT_PREPRINT << "Error " << r << "in prob.nlpHessianEval" << OPT_GETFILELINE << std::endl;
            #endif
            return r;
        }
    }
    
    
    for( int i = 0; i < n; i++ )
    {
        const unsigned int nz = hessRowStartIndex[i+1] - hessRowStartIndex[i];
        double *v = &( hessian[ hessRowStartIndex[i] ] );
        
        if( Q.getNumberOfElementsAtRow(i) == 0 && nqcons[i] == 0 )
        {
            //std::cout << "i: " << i << "\n";
            //prob.writeMIQCPModelInAMPLFile("knitro.mod");
            //OPT_getchar();
            
            if( evalLag )
            {
                int *p = NULL;
                #if OPT_DEBUG_MODE
                    int rnz;
                    p = &rnz;
                #endif
                lagH.getRowValues(i, v, p);
                    
                #if OPT_DEBUG_MODE
                    assert( nz == rnz );
                #endif
            }
            else
            {
                OPT_setAllArray(nz, v, 0.0);
            }
        }
        else if( lagH.getNumberOfElementsAtRow(i) == 0 && nqcons[i] == 0 )
        {
            
            const double f = objf*prob.objFactor;
            
            if( f == 0.0 )
            {
                OPT_setAllArray(nz, v, 0.0);
            }
            else
            {
                /*MIP_SparseRow &qrow = Q[i];
                const unsigned int nel = qrow.getNumberOfElements();
                unsigned int k;
                
                for(k = 0; k < nel; k++)
                    v[k] = qrow[k].getValue() * f; */
                
                
                const double* rvalues = Q.operator()(i);
                const unsigned int nel = Q.getNumberOfElementsAtRow(i);
                for(unsigned int k = 0; k < nel; k++)
                    v[k] = rvalues[k] * f;
                
                
                #if OPT_DEBUG_MODE
                    assert( nz == nel );
                #endif
            }
        }
        else if( Q.getNumberOfElementsAtRow(i) == 0 && lagH.getNumberOfElementsAtRow(i) == 0 && nqcons[i] == 1 )
        {
            const double f = lambda[ indqcons[i] ];
            
            if( f == 0.0 )
            {
                OPT_setAllArray(nz, v, 0.0);
            }
            else
            {
                /*MIP_SparseRow &qcrow = QC[indqcons[i]][i];
                const unsigned int nel = qcrow.getNumberOfElements();
                unsigned int k;
                
                for(k = 0; k < nel; k++)
                    v[k] = qcrow[k].getValue() * f; */
                
                
                const double* rvalues = QC[indqcons[i]](i);
                const unsigned int nel = QC[indqcons[i]].getNumberOfElementsAtRow(i);
                
                for(unsigned int k = 0; k < nel; k++)
                    v[k] = rvalues[k] * f;
                
                #if OPT_DEBUG_MODE
                    assert( nz == nel );
                #endif
            }
        }
        else
        {
            OPT_setAllArray( nz, v, 0.0 );
                
            OPT_setAllArray( i+1, grad, 0.0 ); //we just use lower triangle
            
            if( evalLag )
                lagH.accumulateRowInArray(i, grad, 1.0); //lagH.getFullRowAccumulation( i, grad );
            
            if( objQ )
                Q.accumulateRowInArray(i, grad,  objf*prob.objFactor); //OPT_accumulateRowInArray(Q, i, objf*prob.objFactor, grad);
            
            if( constrQ )
            {
                
                for(int j = 0; j < mquad; j++)
                {
                    const int ind = quadIndex[j];
                    
                    if( lambda[ind] != 0.0 )
                        QC[ind].accumulateRowInArray(i, grad, lambda[ind]); //OPT_accumulateRowInArray(QC[ind], i, lambda[ind], grad);
                }
            }
            
            const int *hessCol = myknitro->hessCols[i];
            
            
            unsigned int w = 0;
            for( int k = 0; k <= i; k++ ) //just lower triangle
            {
                //columns are ordered in hessCols
                if( grad[k] != 0.0 )
                {
                    while( (int) hessCol[w] != k )
                        w++; //I hope jhessCols has the non zeros index correctly in order. Otherwise, we will have serious problems here...
                    
                    v[w] = grad[k];
                    w++; //we can walk one position in w
                }
            }
            
            #if OPT_DEBUG_MODE
                assert( w <= nz );
            #endif
        }
        
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
    __desallocateSolverEnv();
    deallocateMemory();
}



// __methods from Solver __



void OPT_Knitro::__deallocateSolverEnv()
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
        
        cFnType = pi;*/
        
        int r = OPT_realloc(cFnType, m);
        OPT_IFERRORRETURN(r, r);
        
        #if OPT_HAVE_KNITRO
            if( m > maux )
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
}



int OPT_Knitro::getNumberOfIterations( unsigned int &niter)
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
return OPT_LIBRARY_NOT_AVAILABLE;
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
    __desallocateSolverEnv();
    
    //printf("Entrei em initSolverEnv!");
    //OPT_getchar();
    
    kc = KTR_new();
    if( kc == NULL )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORMSG("Error at knitro KTR_new");
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
return OPT_LIBRARY_NOT_AVAILABLE;
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
return OPT_LIBRARY_NOT_AVAILABLE;
}
#endif



int OPT_Knitro::setLowerObjCut(const double objLBound)
#if OPT_HAVE_KNITRO
{
    //knitro just have a parameter KTR_PARAM_FSTOPVAL to stop if a feasible solution is better than a threshould. Unfortunatelly, that is not enough to us...
    paramChg = true;
    
    return OPT_OPERATION_NOT_IMPLEMENTED;
}
#else
{
return OPT_LIBRARY_NOT_AVAILABLE;
}
#endif


int OPT_Knitro::setUpperObjCut(const double objUBound)
#if OPT_HAVE_KNITRO
{
    //knitro just have a parameter KTR_PARAM_FSTOPVAL to stop if a feasible solution is better than a threshould. Unfortunatelly, that is not enough to us...
    paramChg = true;
    
    return OPT_OPERATION_NOT_IMPLEMENTED;
}
#else
{
return OPT_LIBRARY_NOT_AVAILABLE;
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
return OPT_LIBRARY_NOT_AVAILABLE;
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
return OPT_LIBRARY_NOT_AVAILABLE;
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
return OPT_LIBRARY_NOT_AVAILABLE;
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
return OPT_LIBRARY_NOT_AVAILABLE;
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
return OPT_LIBRARY_NOT_AVAILABLE;
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
return OPT_LIBRARY_NOT_AVAILABLE;
}
#endif



int OPT_Knitro::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_KNITRO
{
    const char solver[] = "knitro";
    
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_double_param_by_name(kc,  param, value);
    
    if( r != 0 )
    {
        //#if OPT_DEBUG_MODE
            //OPT_PRINTERRORNUMBER(r);
        //#endif
        
        printDblParamErrorMsg(!r, solver, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    paramChg = true;
    
    return 0;
}
#else
{
return OPT_LIBRARY_NOT_AVAILABLE;
}
#endif



int OPT_Knitro::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_KNITRO
{
    const char solver[] = "knitro";
    
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_int_param_by_name(kc,  param, value);
    
    if( r != 0 )
    {
        //#if OPT_DEBUG_MODE
            //OPT_PRINTERRORNUMBER(r);
        //#endif
        
        printIntParamErrorMsg(!r, solver, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    paramChg = true;
    
    return 0;
}
#else
{
return OPT_LIBRARY_NOT_AVAILABLE;
}
#endif



int OPT_Knitro::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_KNITRO
{
    const char solver[] = "knitro";
    
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    const int r = KTR_set_char_param_by_name(kc,  param, value);
    
    if( r != 0 )
    {
        //#if OPT_DEBUG_MODE
            //OPT_PRINTERRORNUMBER(r);
        //#endif
        
        printStrParamErrorMsg(r, solver, param, value);
        
        return OPT_BAD_INPUT;
    }
        
    
    paramChg = true;
    
    return 0;
}
#else
{
return OPT_LIBRARY_NOT_AVAILABLE;
}
#endif





int OPT_Knitro::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_KNITRO
{
    const int n = prob.n;
    const int m = prob.m;
    const int nI = prob.getNumberOfIntegerVars();
    
    int r;
    
    const double *lx = prob.lx;
    const double *ux = prob.ux;
    
    KTR_context* &kc = (KTR_context* &) knitroContext;
    
    
    bool resetProblem = false;
    int *jacIndexCons = NULL, *jacIndexVars;
    int *myhessRows = NULL, *myhessCols;
    
    
    
    
    //std::cout << "1 nmChg: " << nmChg << " genConstrChg: " << genConstrChg << " genHessChg: " << genHessChg << " consBoundChg: " << consBoundChg << "\n";
    
    
    if(resetSol)
    {
        origSolverRetCode = KTR_RC_NULL_POINTER; //I think we have to put spme knitro code here...
        this->resetSol();
    }
    
    retCode = OPT_UNDEFINED;
    feasSol = false;
    
    
    
    if( nmChg || genConstrChg || genHessChg || consBoundChg )
        resetProblem = true;
    
    
    
    if( genQuadConstrChg || mquad < 0 )
    {
        mquad = prob.getNumberOfQuadMatricesInConstrs();
        
        if( quadIndex )
            free( quadIndex );
        
        if( mquad > 0 )
        {
            quadIndex = (int*) malloc( mquad * sizeof(int) );
            
            if( !quadIndex )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                retCode = OPT_MEMORY_ERROR;
                goto termination;
            }
            
            prob.getQuadMatrixConstraintInds( quadIndex );
        }
        else
            quadIndex = NULL;
        
    }
    
    
    
    
    if( genQuadConstrChg || nmChg )
    {
        //MIP_SparseMatrix *QC = prob.QC;
        
        int *pi = (int *) realloc( nqcons,  2*n * sizeof(int) );
        if( !pi )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            retCode = OPT_MEMORY_ERROR;
            goto termination;
        }
        
        nqcons = pi;
        indqcons = &(nqcons[n]);
        
        /*for(int i = 0; i < n; i++)
        {
            nqcons[i] = 0;
            for( int j = 0; j < mquad; j++ )
            {
                if( QC[ quadIndex[j] ][i].getNumberOfElements() > 0 )
                {
                    nqcons[i]++;
                    indqcons[i] = quadIndex[j];
                }
            }
        } */
        
        
        OPT_calcNqconstQaudIndex(n, mquad, quadIndex, prob.QC, nqcons, indqcons);
        
        
        genQuadConstrChg = false;
    }
    
    
    
    
    if( nmChg )
    {
        int r;
        
        
        OPT_setAllArray(m, auxCEval, true);
        
        desallocateAuxDerivativeIndexStructures();
        
        
        double *pd = (double *) realloc( auxLambda, (n+m)*sizeof(double) );
        
        if( !pd )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            retCode = OPT_MEMORY_ERROR;
            goto termination;
        }
        
        auxLambda = pd;
        
        
        
        r = allocateAuxDerivativeIndexStructures();
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            retCode = r;
            goto termination;
        }
        
        setFullJacIndex( jacRowStartIndex, jacCols );
        setFullHessIndex( mquad, quadIndex, hessRowStartIndex, hessCols );
        
        genConstrChg = false;
        genHessChg = false;
        nmChg = false;
        
        OPT_setAllArray( m, constrChg, false );	
        OPT_setAllArray( n, quadObjChg, false );
        OPT_setAllArray( n, rowQuadConstrChg, false);
        OPT_setAllArray( n, hessChg, false );
        
        /*for(int i = 0; i < n; i++)
        {
            const int nzl = hessRowStartIndex[i+1] - hessRowStartIndex[i];
            
            std::cout << "hessRowStartIndex["<<i<<"]: " << hessRowStartIndex[i] << " nzl: " << nzl << "\n";
            
            for(int j = 0; j < nzl; j++)
            {
                std::cout << "hessCol["<<i<<"]["<<j<<"]: " << hessCols[i][j] << "\t";
            }
            
            std::cout << "\n";
        }
        
        std::cout << "hessRowStartIndex["<<n<<"]: " << hessRowStartIndex[n] << "\n"; */
        
    }
    
    
    
    
    if( resetProblem )
    {
        bool haslambdaInit = true;
        int r;
        int objType;
        int *cType = auxIndex;
        int *xType = auxIndex2;
        
        double *lc = auxValues;//prob.lc;
        double *uc = auxValues2;//prob.uc;
        
        const bool *nlConstr = prob.nlConstr;
        MIP_SparseMatrix *QC = prob.QC;
        
        
        
        
        
        if( genConstrChg )
        {
            for(int i = 0; i < m; i++)
            {
                if( constrChg[i] )
                {
                    unsigned int nzl;
                    int dif;
                    
                    r = setJacIndexRow(i, jacCols[i], nzl);
                    
                    if( r != 0 )
                    {
                        #if OPT_DEBUG_MODE
                            OPT_PRINTERRORNUMBER(r);
                        #endif
                        retCode = OPT_UNDEFINED_ERROR;
                        goto termination;
                    }
                    
                    dif = nzl - ( (int) jacRowStartIndex[i+1] - jacRowStartIndex[i] ); //diference between current shift and new shift to jacIndexBase
                    
                    if( dif != 0 ) 
                    {
                        //we have to shift lines start indices
                        for( int j = i+1; j <= m; j++ )
                            jacRowStartIndex[j] += dif;
                    }
                    
                    constrChg[i] = false;
                    
                    //reoptimize = false;
                }
            }
            
            genConstrChg = false;
        }
        
        
        //std::cout << "genHessChg: " << genHessChg << "\n";
        
        if( genHessChg )
        {
            for( int i = 0; i < n; i++ )
            {
                if( hessChg[i] || rowQuadConstrChg[i] || quadObjChg[i] )
                {
                    unsigned int nzl;
                    int dif;
                    
                    r = setHessIndexRow( i, mquad, quadIndex, hessCols[i], nzl );
                    if( r != 0 )
                    {
                        #if OPT_DEBUG_MODE
                            OPT_PRINTERRORNUMBER(r);
                        #endif
                        retCode = OPT_UNDEFINED_ERROR;
                        goto termination;
                    }
                    
                    
                    dif = nzl - ( (int) hessRowStartIndex[i+1] - hessRowStartIndex[i] ); //diference between current shift and new shift to hessIndexBase
                    
                    if( dif != 0 )
                    {
                        for(int j = i+1; j <= n; j++)
                            hessRowStartIndex[j] += dif;
                    }
                    
                    hessChg[i] = false;
                    rowQuadConstrChg[i] = false;
                    quadObjChg[i] = false;
                    
                }
                
            }
            
            genHessChg = false;
            
        }
        
        const int nnzJ = jacRowStartIndex[m];
        const int nnzH = hessRowStartIndex[n];
        
        
        
        
        if( prob.hasNlObj )
        {
            objType = KTR_OBJTYPE_GENERAL;
        }
        else
        {
            objType = prob.getNumberOfObjQuadTerms() == 0 ? KTR_OBJTYPE_LINEAR : KTR_OBJTYPE_QUADRATIC ;
        }
        
        
        for(int i = 0; i < m; i++)
        {
            
            if( nlConstr[i] )
                cType[i] = KTR_CONTYPE_GENERAL;
            else
            {
                cType[i] = QC[i].getNumberOfElements() == 0 ? KTR_CONTYPE_LINEAR : KTR_CONTYPE_QUADRATIC ;
            }
        }
        
        if( nnzJ > 0 )
        {
            unsigned int anz = 0;
            
            
            jacIndexCons = (int *) malloc( 2*nnzJ * sizeof(int) );
            
            if( !jacIndexCons )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                retCode = OPT_MEMORY_ERROR;
                goto termination;
            }
            
            jacIndexVars = &jacIndexCons[nnzJ];
            
            
            for(int i = 0; i < m; i++)
            {
                const unsigned int nz = jacRowStartIndex[i+1] - jacRowStartIndex[i];
                
                OPT_setAllArray<int>(nz, &(jacIndexCons[anz]), i);
                
                OPT_copyArray(nz, jacCols[i], &(jacIndexVars[anz]) );
                
                anz += nz;
            }
        }
        
        if( nnzH > 0 )
        {
            unsigned int anz = 0;
            
            myhessRows = (int *) malloc( 2* nnzH * sizeof(int) );
            
            if( !myhessRows )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                retCode = OPT_MEMORY_ERROR;
                goto termination;
            }
            
            
            myhessCols = &myhessRows[nnzH];
            
            for(int i = 0; i < n; i++)
            {
                const unsigned int nz = hessRowStartIndex[i+1] - hessRowStartIndex[i];
                
                OPT_setAllArray( nz, &(myhessRows[anz]), i);
                
                OPT_copyArray( nz, hessCols[i], &(myhessCols[anz])  );
                
                anz += nz;
            }
            
            
            /*std::cout << "opa nnzH: " << nnzH << "\n";
            for(int i = 0; i < nnzH; i++)
            {
                std::cout << "hrows["<<i<<"]: " << myhessRows[i] << " hcols["<<i<<"]: " << myhessCols[i] << "\n";
            }*/
        }
        
        
        if( m > 0 )
        {
            if( isnan(lambdaInit[0]) )
                haslambdaInit = false;
            else if( isnan(zInit[0]) )
                haslambdaInit = false;
            else
            {
                
                OPT_copyArray( m, lambdaInit, auxLambda );
                OPT_copyArray( n, zInit, &(auxLambda[m]) );
            }
        }
        
        /*std::cout << "objType: " << objType << std::endl;
        
        for(int i = 0; i < m; i++)
            std::cout << "cType["<<i<<"]: " << cType[i] << std::endl; */
        
        
        {
            const double *olc = prob.lc, *ouc = prob.uc;
            
            for(int i = 0; i < m; i++)
            {
                if( olc[i] > -KTR_INFBOUND  ||  ouc[i] < KTR_INFBOUND )
                {
                    lc[i] = olc[i];
                    uc[i] = ouc[i];
                }
                else
                {
                    lc[i] = -KTR_INFBOUND/10; //KNITRO does not allow free constraints...
                    uc[i] = KTR_INFBOUND;
                }
            }
        }
        
        //std::cout << "problemInitialized: " << problemInitialized << "\n";
        
        if( problemInitialized )
        {
            OPT_PRINTERRORMSG(".......................................................\nWarning: problem structure or constraint bounds were changed. It is necessary to optsolvers to reinitialize knitro data structures. \nUsers parameters will be lost. To remediate it, call:\n\tinitSolverEnv \n\tso, set again your knitro parameters and call solve.\n......................................................\n");
            
            //KTR_free(&kc);
            initSolverEnv(-1, -1, -1); //initSolverEnv already free memory from knitro
            
            //setIntegerParameter( "derivcheck", 3 );
        }
        
        
        //prob.print();
        
        
        if( nI == 0 )
        {
            //unfortunatelly, hessian must be represented by upper triangle. That is not good to us, but I hope solve this problem just changing the myhessRows and myhessCols parameters order
            
            //std::cout << "kc: " << kc << " n: " << n << " m: " << m << " objType: " << objType << " nnzJ: " << nnzJ << " nnzH: " << nnzH;
            
            
            
            r = KTR_init_problem(kc, n, KTR_OBJGOAL_MINIMIZE, objType, lx, ux, m, cType, lc, uc, nnzJ, jacIndexVars, jacIndexCons, nnzH, myhessCols, myhessRows, newInitPoint ? xInit : NULL, haslambdaInit ? auxLambda : NULL );
            
        }
        else
        {
            int *myxtype = prob.xtype;
            
            #pragma ivdep
            #pragma gcc ivdep
            for(int i = 0; i < n; i++)
                xType[i] = MIP_isIntegerType(myxtype[i]) ? KTR_VARTYPE_INTEGER : KTR_VARTYPE_CONTINUOUS ;
            
            r = KTR_mip_init_problem(kc, n, KTR_OBJGOAL_MINIMIZE, objType, objFnType, xType, lx, ux, m, cType, cFnType, lc, uc, nnzJ, jacIndexVars, jacIndexCons, nnzH, myhessCols, myhessRows,  newInitPoint ? xInit : NULL, haslambdaInit ? auxLambda : NULL );
        }
        
        problemInitialized = true;
        
        newInitPoint = false;
        consBoundChg = false;
        paramChg = false;
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            /*for(int i = 0; i < m; i++)
                std::cout << "i: " << i << " lc: " << lc[i] << " uc: " << uc[i] << "\n"; */
            
            if( r == KTR_RC_ILLEGAL_CALL )
            {
                OPT_PRINTERRORMSG("License error on knitro. Probably size problem is grater than your license allows.");
                retCode = OPT_LICENSE_ERROR;
            }
            else
                retCode = OPT_SOLVER_ERROR;
            goto termination;
        }
        
        
        if( jacIndexCons )	
        {
            free(jacIndexCons);
            jacIndexCons = NULL;
        }
        if( myhessRows )
        {
            free(myhessRows);
            myhessRows = NULL;
        }
        
        
    }
    else
    {
        r = KTR_chgvarbnds( kc, lx, ux );
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            retCode = OPT_SOLVER_ERROR;
            goto termination;
        }
        
        
        if( paramChg || ( newInitPoint && !isnan(xInit[0]) ) )
        {
            double *px, *plambda;
            
            if( newInitPoint && !isnan(xInit[0]) )
            {
                
                OPT_copyArray( m, lambdaInit, auxLambda );
                OPT_copyArray( n, zInit, &(auxLambda[m]) );
                
                px = xInit;
                plambda = auxLambda;
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
    
    
    
    //falta capturar a solução... :(
    
    
    if( nI == 0 )
    {
        origSolverRetCode = KTR_solve( kc, storeSol ? sol : auxValues, auxLambda, 0, &objValue, NULL, NULL, NULL, NULL, NULL, this ); //knitro does not return constraint values in c input argument
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
                std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Knitro solving!\n";
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
        OPT_copyArray( m, auxLambda, dualSolC );
        OPT_copyArray( n, &(auxLambda[m]), dualSolV );
    }
    
    
    if( prob.objFactor < 0 )
    {
        objValue = -objValue;
    }
    
    
    
    
termination:
    
    if( jacIndexCons )	free(jacIndexCons);
    if( myhessRows )	free(myhessRows);
    
    
    return retCode;
}
#else
{
return OPT_LIBRARY_NOT_AVAILABLE;
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





















































