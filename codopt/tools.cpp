
#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"




using namespace minlpproblem;
using namespace optsolvers;






#if OPT_SUPER_THREAD_DEBUG_MODE
    std::map< std::thread::id, std::ostream* >  optsolvers::OPT_thsOut;
    OPT_Mutex SEMAPH_thsOut;
    
    void optsolvers::OPT_createFileToThread( const std::thread::id &tid )
    {
        std::stringstream fname;
        fname << "opt_thread_output_" << tid << ".dbg";
        
        
        SEMAPH_thsOut.lock( );
        
        if( OPT_thsOut.count( tid ) == 0 )
        {
            std::ofstream *fout = new std::ofstream;
            
            fout->open( fname.str() );
            
            OPT_thsOut[ tid ] = fout;
        }
        
        SEMAPH_thsOut.unlock( );
    }
    
    
    void optsolvers::OPT_closeFileToThread( const std::thread::id &tid )
    {
        
        SEMAPH_thsOut.lock( );
        
        if( OPT_thsOut.count( tid ) > 0 )
        {
            auto p = OPT_thsOut[tid];
            OPT_thsOut.erase(tid); //removing the file from map. In this, it will be created again in the next execution of some function being debuged
            delete p; //it will call close method (I hope)
        }
        
        SEMAPH_thsOut.unlock( );
    }
    
    
#endif





OPT_Mutex::OPT_Mutex()
{
    initialize();
}


void OPT_Mutex::initialize()
{
    #if OPT_CPP_MULTITHREADING
        
    #else	
        #if OPT_OMP_MULTITHREADING
            omp_init_lock(&mutex);
        #endif
    #endif
}


int OPT_Mutex::lock( )
{
    
    #if OPT_CPP_MULTITHREADING
        //if( nthreads > 1 )
            mymutex.lock();
    #else
        #if OPT_OMP_MULTITHREADING
            //if( nthreads > 1 )
                omp_set_lock(&mutex);
        #endif
    #endif
    
    return 0;
}


int OPT_Mutex::tryLock(  )
{
    #if OPT_CPP_MULTITHREADING
        //if( nthreads > 1 )
            return mymutex.try_lock() == true ? 0 : OPT_UNDEFINED_ERROR;
    #else
        #if OPT_OMP_MULTITHREADING
            //if( nthreads > 1 )
                return omp_test_lock(&mutex) == 1 ? 0 : OPT_UNDEFINED_ERROR;
        #endif
    #endif
    
    return 0;
}


int OPT_Mutex::unlock( )
{
    #if OPT_CPP_MULTITHREADING
        //if( nthreads > 1 )
            mymutex.unlock();
    #else
        #if OPT_OMP_MULTITHREADING
            //if( nthreads > 1 )
                omp_unset_lock(&mutex);
        #endif
    #endif
    
    return 0;
}


void OPT_Mutex::destroy()
{
    #if OPT_OMP_MULTITHREADING
        omp_destroy_lock(&mutex);
    #endif
}


OPT_Mutex::~OPT_Mutex()
{
    destroy();
}






int optsolvers::OPT_getSingleJacPointerOrCalcJacIndices( const OPT_MINLPProb &prob, const int mquad, const bool hasFreeConstrs, const bool *constrsEval, const OPT_SparseMatrix* &singleJ, int &nnz_jac_g, int* &sizeColsNzRowJac, int** &colsNzRowJac)
{
    nnz_jac_g = 0;
    
    singleJ = OPT_getSingleJacobianPointer(prob, mquad);
    
    
    if(singleJ && !hasFreeConstrs)
    {
        nnz_jac_g = singleJ->getNumberOfElements();
        
        //if we have only one matrix composing Jacobian, we do not need allocate colsNzRowJac
    }
    else
    {
        int ret;
        const int om = prob.m;
        
        singleJ = NULL;
        
        //std::cout << "om: " << om << "\n";
        //OPT_getchar();
        
        colsNzRowJac = (int**) calloc( om, sizeof(int*) );
        sizeColsNzRowJac = (int*) malloc( om* sizeof(int) );
        OPT_IFMEMERRORRETURN(!colsNzRowJac || !sizeColsNzRowJac);
        
        
        ret = OPT_fillColsNzRowJac(prob, constrsEval, colsNzRowJac, sizeColsNzRowJac);
        OPT_IFERRORRETURN(ret, ret);
        
        
        nnz_jac_g = 0;
        for(int i = 0; i < om; i++)
            nnz_jac_g += sizeColsNzRowJac[i];
        
    }
    
    return 0;
}





int optsolvers::OPT_getSingleLagHPointerOrCalcLagHIndices( const OPT_MINLPProb &prob, const int mquad, const int *quadIndex, const bool *constrsEval, const OPT_SparseMatrix* &singleLagH, int &nnz_h_lag, int &nNzRowsLagH, int* &nzRowsLagH, int* &sizeColsNzLagH, int** &colsNzRowLagH)
{
    nnz_h_lag = 0;
    
    //testing if we have a hessian composed by more than one sparse matriz
    
    singleLagH = OPT_getSingleLagHPointer(prob, mquad, quadIndex);
    
    
    if(singleLagH)
    {
        nnz_h_lag = singleLagH->getNumberOfElements();
    }
    else
    {
        int ret;
        std::unordered_set<int> setnzRowsLagH;
        
        
        ret = OPT_fillNzRowsLagH(prob, constrsEval, mquad, quadIndex, setnzRowsLagH);
        OPT_IFERRORRETURN(ret, ret);
        
        
        nNzRowsLagH = setnzRowsLagH.size();
        
        OPT_malloc(nzRowsLagH, nNzRowsLagH); //nzRowsLagH = (int *) malloc( nNzRowsLagH * sizeof(int) );
        OPT_IFMEMERRORRETURN(!nzRowsLagH);
        
        
        {
            int i = 0;
            for(auto it: setnzRowsLagH)
            {
                nzRowsLagH[i] = it;
                i++;
            }
        }
        
        
        OPT_malloc(sizeColsNzLagH, nNzRowsLagH); //sizeColsNzLagH = (int*) malloc( nNzRowsLagH * sizeof(int) );
        OPT_calloc(colsNzRowLagH, nNzRowsLagH); //colsNzRowLagH = (int**) calloc( nNzRowsLagH , sizeof(int*) );
        OPT_IFMEMERRORRETURN(!sizeColsNzLagH || !colsNzRowLagH);
        
        
        ret = OPT_fillColsNzRowHess(prob, constrsEval, mquad, quadIndex, nNzRowsLagH , nzRowsLagH, colsNzRowLagH, sizeColsNzLagH);
        OPT_IFERRORRETURN(ret, ret);
        
        
        nnz_h_lag = 0;
        
        for(int k = 0; k < nNzRowsLagH; k++)
            nnz_h_lag += sizeColsNzLagH[k];
    }
    
    
    return 0;
}




int optsolvers::OPT_getValuesFromCompleteJacobian( const OPT_MINLPProb &prob, const OPT_SparseMatrix &J, const int *sizeColsNzRowJac, int **colsNzRowJac, const bool *auxCEval, const int newm, const double *x, double *auxVars, double *values)
{
    //const int n = prob.n;
    const int origm = prob.m;
    const bool *nlConstr = prob.nlConstr;
    const OPT_SparseMatrix &A = prob.A;
    const OPT_SparseMatrix *QC = prob.QC;
    
    
    int nzs = 0;
    
    //OPT_setAllArray(n, auxVars, 0.0);
    
    
    for(int i = 0; i < origm; i++)
    {
        if(!auxCEval[i])
            continue;
        
        bool set = false;
        
        
        if( QC[i].getNumberOfElements() == 0 )
        {
            int nzrow;
            
            if( !nlConstr[i] )
            {
                A.getRowValues(i, &values[nzs], &nzrow);
                
                nzs += nzrow;
                set = true;
                
            }
            else if( A.getNumberOfElementsAtRow(i) == 0 )
            {
                J.getRowValues(i, &values[nzs], &nzrow);
                nzs += nzrow;
                set = true;
            }
        }
        
        
        
        if(!set)
        {
            //we will eval the complete gradient. So, we initialize the nonzero cols as zero
            
            
            const int nminds = sizeColsNzRowJac[i];
            const int *pcols = colsNzRowJac[i];
            
            for(int j = 0; j < nminds; j++)
            {
                auxVars[ pcols[j] ] = 0.0; 
            }
            
            
            MIP_constrCompleteGrad(prob, J, i, x, auxVars, false);
            
            
            //std::cout << "i: " << i << "\n";
            
            for(int j = 0; j < nminds; j++)
            {
                //const int ind = pcols[j];
                
                //std::cout << "col " << j << ": " << ind << "\n";
                
                values[nzs] = auxVars[ pcols[j] ];
                //auxVars[ind] = 0.0; //initializing again this position
                nzs++;
            }
            
        }
        
    }
    
    
    return nzs;
}



int optsolvers::OPT_getValuesFromCompleteLagHessian( const OPT_MINLPProb& prob, const OPT_SparseMatrix& lagH, const int mquad, const int *quadIndex, const double obj_factor, const double *lambda, const int nNzRowsLagH, const int *nzRowsLagH, const int *sizeColsNzLagH, int **colsNzRowLagH, double *auxVars, double *values)
{
    //const int n = prob.n;
    double *rowHvalues = auxVars;
    
    
    //OPT_setAllArray(n, rowHvalues, 0.0);
        
    int nzs = 0;
    
    for(int i = 0; i < nNzRowsLagH; i++)
    {
        const int row = nzRowsLagH[i];
        const int ncols = sizeColsNzLagH[i];
        const int *pcols = colsNzRowLagH[i];
        
        //we will eval the complete gradient. So, we initialize the nonzero cols as zero
        for(int j = 0; j < ncols; j++)
            rowHvalues[ pcols[j] ] = 0.0;
        
        
        //prob.objFactor is already considered inside this function
        MIP_completeLagHessianRow(prob, mquad, quadIndex, lagH, obj_factor, lambda, row, rowHvalues, false); //do not use mylambda: we put zero inpositions having nlConstr[i] === false. We also are sure constraints in quadIndex are not redundante.
        
        
        for(int j = 0; j < ncols; j++)
        {
            const int col = pcols[j];
            values[nzs] = rowHvalues[col];
            nzs++;
            //rowHvalues[col] = 0.0;
        }
    }
    
    /*#if OPT_DEBUG_MODE //TODO: remove it some day
        for(int i = 0; i < n; i++)
        {
            assert(rowHvalues[i] == 0.0);
        }
    #endif */
    
    return nzs;
}






int optsolvers::OPT_evalCompleteLagrangianHessian( const int thnumber, const bool newx, const double *x, const OPT_MINLPProb& prob, OPT_SparseMatrix& lagH, const int mquad, const int *quadIndex, const double obj_factor, const double aditional_nl_obj_factor, const double *lambda, const bool* auxCEval, const int newm, const int nNzRowsLagH, const int *nzRowsLagH, const int *sizeColsNzLagH, int **colsNzRowLagH, double *auxConstrsVars, double *auxConstrs, int &sizeValues, double *values)
{
    const double realObjF = prob.objFactor * obj_factor;
    
    const bool evalLag = prob.hasNlObj || prob.hasNlConstrs;
    const int om = prob.m;
    const bool hasFreeConstrs = newm != om;
    
    //const bool *nlConstr = prob.nlConstr;
    const OPT_SparseMatrix &Q = prob.Q;
    const OPT_SparseMatrix *QC = prob.QC;
    
    const int hasH = prob.hasNlObj || prob.hasNlConstrs ? 1 : 0;
    
    const int hasQ = Q.getNumberOfElements() > 0 ? 1 : 0;
    
    int nzs;
    
    
    double *mylambda = auxConstrs;
    const double *pmylambda;
    const double *plambda;
    
    sizeValues = 0;
    
    
    if( hasFreeConstrs )
    {
        int k = 0;
        
        
        //OPT_setAllArray(om, mylambda, 0.0);
            
        for(int i = 0; i < om; i++)
        {
            if(auxCEval[i])
            {
                mylambda[i] = lambda[k];
                k++;
            }
            else
            {
                mylambda[i] = 0.0;
            }
        }
        
        #if OPT_DEBUG_MODE
            assert(k == newm);
        #endif
            
        pmylambda = mylambda;
    }
    else
    {
        //OPT_copyArray(om, lambda, mylambda);
        pmylambda = lambda;
    }
    
    plambda = pmylambda;
    
    
    if(evalLag)
    {
        #if 0
        
        if(hasFreeNLConstrs)
        {
            double *mylambda2 = auxConstrsVars;
            
            OPT_copyArray(om, plambda, mylambda2);
            
            for(int i = 0; i < om; i++)
            {
                if(!nlConstr[i] || !auxCEval[i])
                    mylambda2[i] = 0.0; //we only eval nonlinear (nonredundant) constraints
            }
            
            plambda = mylambda2;
        }
        
        #endif
        
        /* we do not more set lambda on zero at quadratic constraints. We assume is explicit user evaluation callbacks should eval only strict nonlinear terms.
        if( mquad > 0 )
        {
            double *mylambda2 = auxConstrsVars;
            
            OPT_copyArray(om, plambda, mylambda2);
            
            for(int i = 0; i < mquad; i++)
            {
                mylambda2[ quadIndex[i] ] = 0.0; //we do not want evaluate quadratics in 
            }
            
            plambda = mylambda2;
        }*/
        
        
        if(!prob.hasNlConstrs && realObjF == 0.0)
        {
            lagH.setAllSparseMatrix(0.0);
        }
        else
        {
            //prob->objFactor is already considered in prob->nlpHessianEval...
            
            int r = prob.nlpHessianEval(thnumber, newx, x, prob.hasNlObj ? obj_factor * aditional_nl_obj_factor : 0.0, plambda, lagH);
            
            if( r != 0 )
            {
                #if OPT_PRINT_CALLBACK_ERROR_MSG
                    std::cerr << OPT_PREPRINT << "Error " << r << "in prob->constraintsEval" << OPT_GETFILELINE << "\n";
                #endif
                return r;
            }
            
        }
    }
        
    if(hasQ + mquad + hasH == 1)
    {
        if(hasH)
            nzs = lagH.getValues(values);
        
        if(hasQ)
        {
            nzs =Q.getNumberOfElements();
            if(realObjF == 0.0)
            {
                OPT_setAllArray(nzs, values, 0.0);
            }
            else
            {
                Q.getValues(values);
                if(realObjF != 1.0)
                    OPT_multiplyAllArray(nzs, realObjF, values);
            }
        }
        
        if(mquad)
        {
            const int ind = quadIndex[0];
            nzs = QC[ind].getNumberOfElements();
            
            if(pmylambda[ind] == 0.0)
            {
                OPT_setAllArray(nzs, values, 0.0);
            }
            else
            {
                QC[ind].getValues(values);
                if(pmylambda[ind] != 1.0)
                    OPT_multiplyAllArray(nzs, pmylambda[ind], values);
            }
        }
    }
    else
    {
        nzs = OPT_getValuesFromCompleteLagHessian(prob, lagH, mquad, quadIndex, obj_factor, pmylambda, nNzRowsLagH, nzRowsLagH, sizeColsNzLagH, colsNzRowLagH, auxConstrsVars, values); //do not use plambda: we put zero in positions having nlConstr[i] === false. We also are sure constarints in quadIndex are not redundante.
        
        #if 0
        double *rowHvalues = nlp->auxValues2;
        
        OPT_setAllArray(n, rowHvalues, 0.0);
        
        
        nzs = 0;
        
        for(int i = 0; i < nNzRowsLagH; i++)
        {
            const int row = nzRowsLagH[i];
            const int ncols = sizeColsNzLagH[i];
            const int *pcols = colsNzRowLagH[i];
            
            //prob.objFactor is already considered inside this function
            MIP_constrCompleteLagHessianRow(*prob, mquad, quadIndex, lagH, obj_factor, lambda, row, x, rowHvalues, false);
            
            
            for(int j = 0; j < ncols; j++)
            {
                const int col = pcols[j];
                
                
                values[nzs] = rowHvalues[col];
                nzs++;
                rowHvalues[col] = 0.0;
            }
        }
        
        //TODO: remove it some day
            for(int i = 0; i < n; i++)
            {
                assert(rowHvalues[i] == 0.0);
            }
        
        #endif
    }
    
    
    
    
    sizeValues = nzs;
    
    return 0;
}








