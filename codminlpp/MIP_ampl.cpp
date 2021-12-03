
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>


#include "MIP_ampl.hpp"
#include "MIP_tools.hpp"

#include <new>
#include <iostream>


//using namespace std;
using namespace minlpproblem;









MIP_NLEvalAmlp::MIP_NLEvalAmlp(const bool quadObj, const bool quadConstrs, char* stub)
{
    //quadObj = quadConstrs = true;
    this->quadObj = quadObj;
    this->quadConstrs = quadConstrs;
    
    #if MIP_HAVE_ASL
        aslv = NULL;
    #endif
    this->stub = stub;
    //ctrAux = NULL;
    auxVals = NULL;
    
    #if MIP_AMPL_DEBUG_MODE
        checkJacInd = false;
        checkHessInd = false;
    #endif
    
    nEvalErrors = 0;
}


MIP_NLEvalAmlp::~MIP_NLEvalAmlp()
{
    desallocate();
}


void MIP_NLEvalAmlp::desallocate()
{
    #if MIP_HAVE_ASL
        
        if(aslv)
        {
            for(unsigned int i = 0; i < nthreads; i++)
            {
                if(aslv[i])
                    ASL_free(&aslv[i]);
            }
            
            free(aslv);
            aslv = NULL;
        }
        if(auxVals)
        {
            for(unsigned int i = 0; i < nthreads; i++)
            {
                if(auxVals[i])
                    free(auxVals[i]);
            }
            
            free(auxVals); 
            auxVals = NULL;
        }
        
        /*if(ctrAux)
        {
            if(ctrAux[0])
                free(ctrAux[0]);
            
            free(ctrAux);
            ctrAux = NULL;
        } */
    
    #endif
}




int MIP_NLEvalAmlp::initialize(const int nthreads, const int n, const int m, const int nzJac, const int nzLagHess)
#if MIP_HAVE_ASL
{
    unsigned int j;
    int aux, aux2;
    FILE *nl = NULL;
    
    //printf("Entrei na callback de inicializacao do AMPL\n");
    
    #if MIP_AMPL_DEBUG_MODE
        checkJacInd = true;
        checkHessInd = true;
    #endif
    
    nEvalErrors = 0;
    
    desallocate();
    
    MIP_calloc(aslv, nthreads); //aslv = (ASL **) malloc(nthreads *sizeof(ASL *));
    if(!aslv)
    {
        #if MIP_MEMORY_ERROR
            MIP_PRINTMEMERROR;
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    
    for(decltype(this->nthreads) i = 0; i < (unsigned) nthreads; i++)
    {
        aslv[i] = ASL_alloc(ASL_read_pfgh);
        
        nl = jac0dim_ASL(aslv[i], stub, (fint)strlen(stub));
        
        //making only the non-linear constraints be evaluated at conval and jacval functions
        aslv[i]->i.n_conjac_[0] = 0;
        aslv[i]->i.n_conjac_[1] = aslv[i]->i.nlc_;
        aslv[i]->i.congrd_mode = 2;
        
        pfgh_read_ASL(aslv[i], nl, 0);
        
        aux = quadObj ? -1 : 0;
        aux2= quadConstrs ? 0 : 1;
        
        aslv[i]->p.Sphset(aslv[i], 0, aux, 1, aux2, 1); //to eval hessian...
    }
    
    
    MIP_calloc(auxVals, nthreads); //auxVals = (double **) malloc( nthreads * sizeof(double *) );
    //ctrAux = (double **) malloc( nthreads * sizeof(double *) );
    if( !auxVals ) //|| !ctrAux )
    {
        #if MIP_MEMORY_ERROR
            MIP_PRINTMEMERROR;
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    
    {
        unsigned int aux = quadConstrs ? 0 : aslv[0]->i.nzc_ ;
        j = MIP_max<unsigned int>( aux, nzLagHess );
    }
    
    //MIP_malloc(auxVals[0], j*nthreads);
    /*if(!auxVals[0]) // || !ctrAux[0])
    {
        #if MIP_MEMORY_ERROR
            MIP_PRINTMEMERROR;
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    
    for(i = 1; i < nthreads; i++)
    {
        auxVals[i] = &auxVals[0][ i*j ];
        //ctrAux[i] = &ctrAux[0][ i*aslv[0]->i.n_con_ ];
    } */
    
    
    for(unsigned int i = 0; i < (unsigned) nthreads; i++)
    {
        MIP_malloc( auxVals[i], j );
        MIP_IFMEMERRORRETURN(!auxVals[i]);
    }
    
    
    /*for(i = 0; i < nthreads; i++)
    {
        //to linear constraints in hessian evaluation
        for(j = m; j < aslv[0]->i.n_con_; j++) 
            ctrAux[i][j] = 0.0;
    } */
    
    
    this->nthreads = nthreads;
    
    
    return 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif


int MIP_NLEvalAmlp::eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value)
#if MIP_HAVE_ASL
{
    fint nerror = 1;
    ASL *asl = aslv[threadnumber];
    void *p = &x;
    double *myx = *( (double **) p); //countourning the problem of vector x be constraint
    
    if( newx ) //maybe it is better put it inside if( nlobj ), but maybe it is a good idea do ampl knows about that solution anyway...
        asl->p.Xknown(asl, myx, 0); 
    
    
    
    value = objval(0, myx, &nerror);
    
    if(nerror != 0)
    {
        nEvalErrors++;
        
        #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
            if( nEvalErrors <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                std::cerr << MIP_PREPRINT  "Error " << nerror << " at AMPL objective function evaluation" << MIP_GETFILELINE << "\n";
        #endif
        
        return nerror;
    }
    
    
    
    return 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif


int MIP_NLEvalAmlp::eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *values)
#if MIP_HAVE_ASL
{
    fint nerror = 1;
    ASL *asl = aslv[threadnumber];
    void *p = &x;
    double *myx = *( (double **) p);
    
    
    if(newx)
        asl->p.Xknown(asl, myx, 0);
    
    
    objgrd(0, myx, values, &nerror);
    
    if( nerror != 0 )
    {
        nEvalErrors++;
        #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
            if( nEvalErrors <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                std::cerr << MIP_PREPRINT "Error " << nerror << " at AMPL objective function gradient evaluation" << MIP_GETFILELINE << "\n";
        #endif
        
        return nerror;
    }
    
    
    
    return 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif


int MIP_NLEvalAmlp::eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *ctrEval, const double *x, double *values)
#if MIP_HAVE_ASL
{
    ASL *asl = aslv[threadnumber];
    const int mnl = asl->i.nlc_;
    bool evalAll = true;
    bool evalError = false;
    fint nerror = 1; //we should put a nonegative value to avoid AMPL abort if some evaluation error occour
    void *p = &x;
    double *myx = *( (double **) p); //countourning the problem of vector x be constant
    
    
    if(newx)
        asl->p.Xknown(asl, myx, 0);
    
    if(ctrEval)
    {
        for( int i = 0; i < mnl; i++ )
        {
            if( !ctrEval[i] )
            {
                evalAll = false;
                break;
            }
        }
    }
    //if user do not pass ctrEval, we consider it want evaluate all constraints
    
    //Unfourtunatelly, AMPL only does not allow separate quadratic constraint from general nonlinear ones. So, we only take advantage of quadratic constraints constraints if all of them are quadratic...
    
    if(evalAll)
    {
        nerror = 1;
        asl->p.Conval(asl, myx, values, &nerror);
        
        if(nerror)
        {
            nEvalErrors++;
            
            #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
                if( nEvalErrors <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                    std::cerr << MIP_PREPRINT "Error " << nerror << " at AMPL constraint " " evaluation" << MIP_GETFILELINE << "\n";
            #endif
            
            evalError = true;
        }
    }
    else
    {
        for(int i = 0; i < mnl; i++)
        {
            if(ctrEval[i])
            {
                nerror = 1;
                
                //values[i] = conival(i, myx, &nerror);
                values[i] = asl->p.Conival(asl, i, myx, &nerror);
                
                if(nerror)
                {
                    nEvalErrors++;
                    
                    #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
                        if( nEvalErrors <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                            std::cerr << MIP_PREPRINT "Error " << nerror << " at AMPL constraint " << i << " evaluation" << MIP_GETFILELINE << "\n";
                    #endif
                    
                    evalError = true;
                }
            }
        }
    }
    
    
    /*if(evalError)
    {
        double s = 0.0;
        for(i = 0; i < n; i++)
        {
            s += x[i];
            printf("x[%d]: %0.2f ", i, x[i]);
        }
        printf("\ns: %f\n", s);
    } */
    
    
    return evalError ? -1 : 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif


int MIP_NLEvalAmlp::eval_grad_nl_constrs_part( const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *ctrEval, const double *x, MIP_SparseMatrix& J)
#if MIP_HAVE_ASL
{
    bool evalAll = true;
    bool evalError = false;
    fint nerror = 1; //we should put a nonegative value to avoid AMPL abort if some evaluation error occour
    //double aux;
    
    
    ASL *asl = aslv[threadnumber];
    const int mnl = asl->i.nlc_;
    cgrad *cg;
    //double* const auxVals = J.getRowValuesPointer(0);
    double* const auxVals = this->auxVals[threadnumber];
    void *p = &x;
    double *myx = *( (double **) p); //countourning the problem of vector x be constant
    
    
    if( ctrEval )
    {
        for( int i = 0; i < mnl; i++ )
        {
            if( !ctrEval[i] )
            {
                evalAll = false;
                break;
            }
        }
    }
    
    
    if(newx)
        asl->p.Xknown(asl, myx, 0); //xknown(myx);
    
    
    //evalAll = true;
    
    
    if( evalAll )
    {
        nerror = 1;
        //std::cout << "Avaliando jac completo" << std::endl;
        asl->p.Jacval( asl, myx, auxVals, &nerror );
        
        if(nerror)
        {
            nEvalErrors++;
            #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
                if( nEvalErrors <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                    std::cerr << MIP_PREPRINT "Error " << nerror << " at AMPL jacobian constraint " " evaluation" << MIP_GETFILELINE << "\n";
            #endif
                
            evalError = true;
        }
    }
    else
    {
        for(int i = 0; i < mnl; i++)
        {
            if(ctrEval[i])
            {
                nerror = 1;
                //congrd(i, myx, auxVals, &nerror);
                asl->p.Congrd(asl, i, myx, auxVals, &nerror);
                
                if(nerror)
                {
                    nEvalErrors++;
                    #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
                        if( nEvalErrors <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                            std::cerr << MIP_PREPRINT "Error " << nerror << " at AMPL jacobian constraint " << i << " evaluation" << MIP_GETFILELINE << "\n";
                    #endif
                        
                    evalError = true;
                }
            }
        }
        
    }
    
    const int nnzrows = J.nNzRowIndices;
    const unsigned int *nzrow = J.nzRowIndex;
    
    for(int w = 0; w < nnzrows; w++)
    {
        const int i = nzrow[w];
        
        if( i < mnl && (evalAll || ctrEval[i]) )
        {
            //cg = Cgrad[i];
            int k = 0;
            
            //J.setRowValues(i, &auxVals[Cgrad[i]->goff] );
            
            double *rvalues = J.getRowValuesPointer(i);
            
            for(cg = Cgrad[i]; cg; cg = cg->next, k++)
            {
                rvalues[k] = auxVals[cg->goff];
            }
            
            
            
            #if MIP_AMPL_DEBUG_MODE
            
            assert(k == J.getNumberOfElementsAtRow(i));
            
            if( checkJacInd )
            {
                int *rcols = J.getRowColsPointer(i);
                
                int k = 0;
                for(cg = Cgrad[i]; cg; cg = cg->next, k++)
                    assert(rcols[k] == cg->varno);
                
                checkJacInd = false;
            }
            #endif
        }
    }
    
    //for(int w = 0; w < aslv[0]->i.nzc_ ; w++)
        //std::cout << "auxVals["<<w<<"]: " << auxVals[w] << "\n";
    //J.printSparseMatrix();
    
    #if 0
    for(i = 0; i < mnl; i++)
    {
        if(ctrEval[i])
        {
            MIP_SparseRow &jrow = J[i];
            
            for(k = 0, cg = Cgrad[i]; cg; cg = cg->next, k++)
            {
                //aux = auxVals[cg->goff];
                
                /*if( isinf(aux)  )
                {
                    if( aux > 0 )
                        aux = MIP_INFINITY; // +Infinity
                    else 
                        aux = -MIP_INFINITY; // -Infinity
                } */
                
                //J.setElementNoCheckSymmetry(i, cg->varno, auxVals[cg->goff]);
                
                jrow[k].setValue( auxVals[cg->goff] );
                
                #if MIP_AMPL_DEBUG_MODE
                    assert( (int) jrow[k].getColumn() == cg->varno );
                #endif
            }
            
        }
    }
    #endif
    
    /*#if MIP_AMPL_DEBUG_MODE
        if(checkJacInd)
        {
            unsigned int k = 0;
            int* const jacCols = J.getRowColsPointer(0);
            
            for(int i = 0; i < mnl; i++)
            {
                for( cg = Cgrad[i]; cg; cg = cg->next )
                {
                    assert( jacCols[k]  == cg->varno );
                    k++;
                }
            }
            //checkJacInd = false; TODO: uncomment it
        }
    #endif */
    
    
    return evalError ? -1 : 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif


int MIP_NLEvalAmlp::eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double obj_factor, const double *lambda, MIP_SparseMatrix& hessian)
#if MIP_HAVE_ASL
{
    fint nerror; //we should put a nonegative value to avoid AMPL abort if some evaluation error occour
    void *p = &x;
    double *myx = *( (double **) p);
    ASL *asl = aslv[threadnumber];
    double* const auxVals =  hessian.getRowValuesPointer(0) ;//this->auxVals[threadnumber];
    //double* const ctrAux = this->ctrAux[threadnumber];
    double objFactor = quadObj ? 0 : obj_factor;
    void *p2 = &lambda;
    double *mylambda = *( (double **) p2);
    
    
    
    //hessian.setAllSparseMatrix(NAN);
    
    
    if(newx)
    {
        asl->p.Xknown(asl, myx, 0);
        
        //It is ridiculous, but for guarantee that we are evaluating in the desired point, we need to make a function evaluation because the hessian eval does not get the curretn point as a parameter. That function assumes we want evaluate in the last point evaluated in the objective and constraints functions...
        
        //if(obj_factor != 0.0)
        if( !quadObj )
        {
            nerror = 1;
            asl->p.Objval(asl, 0, myx, &nerror);
        }
        
        if( !quadConstrs && m > 0 )
        {
            nerror = 1;
            asl->p.Conival(asl, 0, myx, &nerror);
        }
    }
    
    
    //we eval the hessian only in the general NLP constraints... do not use lambda directly on sphes because ASL acess positions referent to linear cosntraints too (ridiculous...)
    //for(i = 0; i < m; i++)
        //ctrAux[i] = lambda[i];
    
    
    nerror = 1;
    asl->p.Sphes(asl, 0, auxVals, quadObj ? - 1 : 0, &objFactor,  quadConstrs ?  NULL:mylambda);
    
    
    #if MIP_AMPL_DEBUG_MODE
        if(checkHessInd)
        {
            fint *hrownos = sputinfo->hrownos;
            fint *hcolstarts = sputinfo->hcolstarts;
            
            unsigned int nEls = hcolstarts[n]; //size 
            auto hessCols = hessian.getRowColsPointer(0);
            
            /*for(unsigned int i = 0; i < nEls; i++)
            {
                //std::cout << "hessCols["<<i<<"]: " << hessCols[i] << "  hrownos["<<i<<"]: " << hrownos[i] <<  "\n";
                assert( hessCols[i] == hrownos[i] );
            } */
            
            
            for(unsigned int i = 0; i < n; i++)
            {
                assert(hcolstarts[i+1] - hcolstarts[i] == hessian.getNumberOfElementsAtRow(i) );
                
                for(unsigned int j = hcolstarts[i]; j < hcolstarts[i+1] ; j++ )
                {
                    assert( hessCols[j] == sputinfo->hrownos[j] );
                }
            }
            
            
            //MIP_getchar();
            
            checkHessInd = false; //TODO: uncomment it
        }
    #endif
    
    #if 0
    #if MIP_AMPL_DEBUG_MODE
        static bool check = true;
    #endif
        
    for(i = 0; i < n; i++)
    {
        /*aux = hcolstarts[i+1];
        
        //for(j = sputinfo->hcolstarts[i]; j < sputinfo->hcolstarts[i+1]; j++)
        for(j = hcolstarts[i]; j < aux; j++)
        {
            //hessian.setElement(i, sputinfo->hrownos[j], auxVals[j]);
            
            //AMPL evals the upper triangle and sparse matrix store the lower. So, we pass i as col as hrownos[j] as row... 
            
            hessian.setElementNoCheckSymmetry(i, hrownos[j], auxVals[j]);
        } */
        
        
        #if MIP_AMPL_DEBUG_MODE
            if( check )
            {
                MIP_SparseRow &row = hessian[i];
                
                for( int j = hcolstarts[i], w = 0; j < hcolstarts[i+1]; j++, w++)
                {
                    assert( hrownos[j] == row[w].getColumn() ); //if we pass by here, we are sure we can set directly all sparse matrix line using auxVals
                }
            }
        #endif
        
        
        hessian[i].setElementsByOrder( hcolstarts[i+1] -hcolstarts[i], &auxVals[ hcolstarts[i] ] );
    }
    
    
    #if MIP_AMPL_DEBUG_MODE
        check = false;
    #endif
    #endif
    
    /*printf("objFactor: %f m: %d quadObj: %d\n", objFactor, n_con, (int) quadObj);
    
    for(int i = 0; i < m; i++)
        printf("l[%d]: %0.6f \t", i, mylambda[i]);
    
    printf("\n");
    hessian.printSparseMatrix();
    getchar(); */
    
    return 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif



void MIP_NLEvalAmlp::finalize(const int nthreads, const int n, const int mnl, const int nzNLJac, const int nzNLLagHess)
{
    desallocate();
}







MIP_ReadAmplModel::MIP_ReadAmplModel()
{
    #if MIP_HAVE_ASL
        asl = NULL;
    #endif
}

MIP_ReadAmplModel::~MIP_ReadAmplModel()
{
    desallocate();
}


void MIP_ReadAmplModel::desallocate()
{
    #if MIP_HAVE_ASL
        if(asl)
        {
            /*if(asl->i.Urhsx_)	free(asl->i.Urhsx_);
            if(asl->i.Uvx_) 	free(asl->i.Uvx_);
            if(asl->i.havex0_) 	free(asl->i.havex0_);
            if(asl->i.X0_) 	free(asl->i.X0_); */
            
            ASL_free(&asl);
            asl = NULL;
        }
    #endif
}


bool MIP_ReadAmplModel::isMaximizationProblem()
{
    #if MIP_HAVE_ASL
        return asl->i.objtype_[0];
    #else
        return false;
    #endif
}




int MIP_ReadAmplModel::readProblem(char *stub, MIP_MINLPProb &prob, MIP_NonLinearEval* &eval)
#if MIP_HAVE_ASL
{
    bool quadObj = false, quadConstrs;
    int ret, i, j, k, r;
    int code;
    
    int n, m;
    fint maxVals, aux; //we use fint to keep compatibility on ASL for windows 
    int *iRow = NULL, *jCol = NULL;
    fint *nzqp = NULL, **rowqp = NULL, **colqp = NULL; //we use fint to keep compatibility on ASL for windows
    double **delsqp = NULL, *vals = NULL;
    //double *lc = NULL, *uc = NULL, *l = NULL, *u = NULL;
    
    //MIP_NLEvalAmlp *eval = NULL;
    
    FILE *nl = NULL;
    cgrad *cglc;	//for non-zero values of linear constraints
    ograd *og;		//for linear part of objective function...
    
    
    
    //first, we check if we have a quadratic problem:
    asl = ASL_alloc(ASL_read_fg); //to read the quadratic parts we need ASL_read_fg
    nl  = jac0dim_ASL(asl, stub, (fint)strlen(stub));
    
    if( n_obj > 1)
    {
        MIP_PRINTERRORMSG("Error: multiple objective functions defined!");
        code = MIP_BAD_DEFINITIONS;

        goto desallocate_memory;
    }
    
    
    asl->i.Urhsx_ = (real *)  M1alloc( asl->i.n_con_ *sizeof(real) );
    asl->i.Uvx_ = (real *) M1alloc( asl->i.n_var_ * sizeof(real) );
    asl->i.havex0_ = (char *) M1alloc( asl->i.n_var_ * sizeof(char) );
    asl->i.X0_ = (real *) M1alloc( asl->i.n_var_ * sizeof(real) );
    
    if( !asl->i.Urhsx_ || !asl->i.Uvx_ || !asl->i.havex0_ || !asl->i.X0_ )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        
        code = MIP_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    qp_read(nl, 0);
    
    //printf("nvar: %d  nbv: %d niv: %d\n", n_var, nbv, niv);
    
    
    n = asl->i.n_var_;
    m = asl->i.n_con_;
    
    
    prob.setParametersAndAllocate(n, m);
    
    MIP_malloc(nzqp, (1 + asl->i.nlc_) ); //nzqp  = (int *) malloc( (1 + asl->i.nlc_) * sizeof(int) );
    MIP_calloc(rowqp, (1 + asl->i.nlc_) ); //rowqp = (int **) malloc( (1 + asl->i.nlc_)  * sizeof(int *) );
    MIP_calloc(colqp, (1 + asl->i.nlc_) ); //colqp = (int **) malloc( (1 + asl->i.nlc_) * sizeof(int *) );
    MIP_calloc(delsqp, (1+ asl->i.nlc_) ); //delsqp= (double **) malloc( (1+ asl->i.nlc_) * sizeof(double *) );
    
    if( !nzqp || !rowqp || !colqp || !delsqp )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        
        code = MIP_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    maxVals = 0;
    
    if(asl->i.n_obj_ > 0)
    {
        //reading quadratic terms of objective functions
        nzqp[0] = nqpcheck(0, &rowqp[0], &colqp[0], &delsqp[0]);
        
        if( nzqp[0] >= 0 )
        {
            quadObj = true;
            maxVals = nzqp[0];
        }
        /*else
        {
            quadObj = false;
        }*/
    }
    
    
    //check if all constraints are quadratic. Unfourtunatelly, AMPL only does not allow separate quadratic constraint from general nonlinear ones. So, we only take advantage of quadratic constraints constraints if all of them are quadratic...
    
    quadConstrs = true;
    
    for(i = 1; i <= asl->i.nlc_; i++)
    {
        //printf("checando restricao i: %d asl->i.nlc_: %d\n", i, asl->i.nlc_);
        //index 0 is to objective function...
        nzqp[i] = nqpcheck(-i, &rowqp[i], &colqp[i], &delsqp[i]);
        
        if( nzqp[i] < 0 )
        {
            //that constraint is not quadratic... :(
            quadConstrs = false;
            break;
        }
        
        maxVals = MIP_max(maxVals, nzqp[i]);
    }
    
    
    //counting the nonzero of linear part of nonlinear constraints constraints...
    aux = 0;
    for(k = 1; k <= asl->i.nlc_; k++)
    {
        for( cglc = asl->i.Cgrad_[k-1]; cglc; cglc = cglc->next )
        {
            if(cglc->coef != 0.0)
            {
                aux++;
            }
        }
    }
    maxVals = MIP_max(maxVals, aux);
    
    
    //counting the nonzero of linear constraints...
    aux = 0;
    for(int i = asl->i.nlc_; i < asl->i.n_con_ ; i++ )
    {
        for( cglc = asl->i.Cgrad_[i]; cglc; cglc = cglc->next )
            aux++;
    }
    
    maxVals = MIP_max(maxVals, aux);
    
    MIP_malloc(iRow, maxVals); //iRow = (int *) malloc( maxVals * sizeof(int) );
    MIP_malloc(jCol, maxVals); //jCol = (int *) malloc( maxVals * sizeof(int) );
    MIP_malloc(vals, maxVals); //vals = (double *) malloc(maxVals * sizeof(double));
    
    if(!iRow || !jCol || !vals)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        
        code = MIP_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    
    //setting linear and quadratic part of objective function
    if( quadObj )
    {
        //setting the constant for objective function
        prob.setObjConstant( objconst_ASL(asl, 0) );
        
        //setting the linear part of objective function
        for(og = Ograd[0]; og; og = og->next )
            prob.setObjLinearCoefficient(og->varno, og->coef);
        
        //setting the quadratic part of objective function
        if( nzqp[0] > 0 )
        {
            aux = 0;
            for(i = 0; i < n; i++)
            {
                for(j = colqp[0][i]; j < colqp[0][i+1]; j++)
                {
                    //It return the full matrix and we only want the lower triangle...
                    
                    if(rowqp[0][j] >= i)
                    {
                        iRow[aux] = rowqp[0][j];
                        jCol[aux] = i;
                        vals[aux] = delsqp[0][j]; //if the problem is a maximization, we set objFactor
                        aux++;
                    }
                    
                }
            }
            
            int r = prob.setObjQuadCoefsMatrix( aux, iRow, jCol, vals);
            
            if(r != 0)
            {
                #if MIP_DEBUG_MODE
                    MIP_PRINTERRORNUMBER(r);
                #endif
                code = r;
                goto desallocate_memory;
            }
        }
        
    }
    
    
    if(quadConstrs )
    {
        
        //seting the linear part of the quadratic constraint
        int nzs = 0;
        for(k = 1; k <= asl->i.nlc_; k++)
        {
            for( cglc = asl->i.Cgrad_[k-1]; cglc; cglc = cglc->next )
            {
                if(cglc->coef != 0.0)
                {
                    iRow[nzs] = k-1;//newIndC[k-1];
                    jCol[nzs] = cglc->varno;
                    vals[nzs] = cglc->coef;
                    
                    nzs++;
                }
            }
        }
        
        if( nzs > 0 )
        {
            //const int r = prob.setNewConstraintsLinearElements(nzs, iRow, jCol, vals);
            const int r = prob.setConstraintsLinearPart(nzs, iRow, jCol, vals);
            
            if(r != 0)
            {
                //printf("Error in seting linear coefficients of quadratic constraints: %d\n", aux2);
                MIP_PRINTERRORNUMBER(r);
                MIP_getchar();
                
                code = r;
                goto desallocate_memory;
            }
        }
        
        
        for(k = 1; k <= asl->i.nlc_; k++) //all nonlinear constraints are quadratic...
        {
            nzs = 0;
            for(i = 0; i < n; i++)
            {
                
                //std::cout << "k: " << k << " i: " << i << " asl->i.nlc_: " << asl->i.nlc_ << " nzqp[k]: " << nzqp[k] <<  std::endl;
                
                if( nzqp[k] == 0 )
                {
                    MIP_PRINTERRORMSG("Warning: An inconsistency was found in AMPL Solver Library (ASL): Index of nonlinear constraint being considered linear by ASL nqpcheck function...");
                    continue; //that is really strange. nzqp should not be 0, since it means the constarint is linear. However, we have a teste problem (nesting_6EMA.nl from comp_geo) where i.nlc_ is 545, i.e., we have 545 nonlinear constraints and constraint 474 has nzqp = 0! I hate ASL! I think it is a ASL bug...
                }
                    
                
                for(j = colqp[k][i]; j < colqp[k][i+1]; j++)
                {
                    if(rowqp[k][j] >= i) //we want only the lower triangle...
                    {
                        iRow[nzs] = rowqp[k][j];
                        jCol[nzs] = i;
                        vals[nzs] = delsqp[k][j];
                        nzs++;
                    }
                }
            }
            
            i = prob.setConstraintQuadCoefsMatrix(k-1, nzs, iRow, jCol, vals);
            if(i != 0)
            {
                #if MIP_DEBUG_MODE
                    MIP_PRINTERRORNUMBER(i); //printf("Error in seting quadratic matrix of constraint %d. return: %d\n", k-1, i);
                    MIP_getchar();
                #endif
                
                code = i;
                goto desallocate_memory;
            }
        }
    }
    
    
    //seting now the linear constraints...
    if( prob.getNumberOfLinearCoefs() == 0 )
    {
        int nzs = 0;
        for(int i = asl->i.nlc_; i < asl->i.n_con_ ; i++)
        {
            for( cglc = asl->i.Cgrad_[i]; cglc; cglc = cglc->next, nzs++ )
            {
                iRow[nzs] = i; 
                jCol[nzs] = cglc->varno; 
                vals[nzs] = cglc->coef;
            }
        }
        
        if( nzs > 0 )
        {
            //const int r = prob.setNewConstraintsLinearElements(nzs, iRow, jCol, vals);
            const int r = prob.setConstraintsLinearPart(nzs, iRow, jCol, vals);
            if(r != 0)
            {
                #if MIP_DEBUG_MODE
                    MIP_PRINTERRORNUMBER(r);// printf("Error in seting coefficients of linear constraints: %d\n", aux2);
                    MIP_getchar();
                #endif
                
                code = r;
                goto desallocate_memory;
            }
        }
    }
    else
    {
        //We already have linear coeficients (problem from quadratic constraints). So, we set each new constraint one by one...
        
        
        for(int i = asl->i.nlc_; i < asl->i.n_con_ ; i++)
        {
            int nzs = 0;
            for( cglc = asl->i.Cgrad_[i]; cglc; cglc = cglc->next, nzs++ )
            {
                jCol[nzs] = cglc->varno; 
                vals[nzs] = cglc->coef;
            }
            
            if( nzs > 0 )
            {
                const int r = prob.setConstraintLinearPart(i, nzs, jCol, vals);
                
                if(r != 0)
                {
                    #if MIP_DEBUG_MODE
                        MIP_PRINTERRORNUMBER(r); // printf("Error in seting coefficients of linear constraints: %d\n", aux2);
                        MIP_getchar();
                    #endif
                    
                    code = r;
                    goto desallocate_memory;
                }
            }
        }
        
    }
    
    
    free(rowqp);
    free(colqp);
    free(delsqp);
    
    rowqp = colqp = NULL;
    delsqp = NULL;
    
    free(iRow);
    free(jCol);
    free(vals);
    
    iRow = jCol = NULL;
    vals = NULL;
    
    
    
    //setting rhs..
    for(int i = 0; i < asl->i.nlc_ ; i++ )
    {
        if( asl->i.LUrhs_[i] > negInfinity )
        {
            prob.setConstraintLowerBound( i, asl->i.LUrhs_[i] );
            
            //if( asl->i.Urhsx_[i] < Infinity )
                //fprintf(stderr, "Warning: nonconvexity detected\n");
        }
        if( asl->i.Urhsx_[i] < Infinity )
            prob.setConstraintUpperBound(i, asl->i.Urhsx_[i] );
    }
    
    
    for(int i = asl->i.nlc_; i < asl->i.n_con_; i++ )
    {
        if( asl->i.LUrhs_[i] > negInfinity )
            prob.setConstraintLowerBound(i, asl->i.LUrhs_[i] );
        
        if( asl->i.Urhsx_[i] < Infinity )
            prob.setConstraintUpperBound(i, asl->i.Urhsx_[i] );
    }
    
    
    for(int i = 0; i < n; i++)
    {
        if( asl->i.LUv_[i] > negInfinity )
            prob.setVariableLowerBound(i, asl->i.LUv_[i] );
        
        if( asl->i.Uvx_[i] < Infinity )
            prob.setVariableUpperBound(i, asl->i.Uvx_[i] );
    }
    
    
    
    
    
    //setting variable types
    
    //integer in an objective and in a nonlinear constraint
    for(i = nlvb - nlvbi ; i < nlvb; i++)
    {
        prob.setVariableType(i, MIP_VT_INTEGER);
    }
    
    //integer just in nonlinear constraints
    for(i = nlvc - nlvci; i < nlvc; i++)
    {
        prob.setVariableType(i, MIP_VT_INTEGER);
    }
    
    //integer just in objectives
    for(i = nlvo - nlvoi; i < nlvo ; i++)
    {
        prob.setVariableType(i, MIP_VT_INTEGER);
    }
    
    //binary and other integer (linear)
    for(i = n_var - nbv - niv; i < n_var; i++)
    {
        prob.setVariableType(i, MIP_VT_INTEGER);
    }
    
    
    
    
    //setting the initial solution
    for(int i = 0; i < n; i++)
    {
        if( asl->i.havex0_[i] )
        {
            //prob.setInitialSolution(i, asl->i.X0_[i]);
        }
    }
    
    if ( asl->i.objtype_[0])
    {
        //maximization problem
        prob.setObjFactor(-1.0);
    }
    
    
    //Now, we address the general non-linearity in the problem
    //restoring the pointers in ASL
    qp_opify();
    
    
    /*if(asl->i.Urhsx_)	free(asl->i.Urhsx_);
    if(asl->i.Uvx_) 	free(asl->i.Uvx_);
    if(asl->i.havex0_) 	free(asl->i.havex0_);
    if(asl->i.X0_) 		free(asl->i.X0_);
    
    asl->i.Urhsx_ = NULL;
    asl->i.Uvx_ = NULL;
    asl->i.havex0_ = NULL;
    asl->i.X0_ = NULL; */
    
    
    
    if( !quadObj || !quadConstrs )
    {
        ASL_free(&asl);
        
        eval = new (std::nothrow) MIP_NLEvalAmlp(quadObj, quadConstrs, stub);
        if(!eval)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTMEMERROR;
            #endif
            code = MIP_MEMORY_ERROR;
            goto desallocate_memory;
        }
        
        prob.setNonLinearEvaluationObject(eval);
        
        asl  = ASL_alloc(ASL_read_pfgh);
        
        nl  = jac0dim_ASL(asl, stub, (fint)strlen(stub));
    
        asl->i.Urhsx_ = (real *)  M1alloc( asl->i.n_con_ *sizeof(real) );
        asl->i.Uvx_ = (real *) M1alloc( asl->i.n_var_ * sizeof(real) );
        asl->i.havex0_ = (char *) M1alloc( asl->i.n_var_ * sizeof(char) );
        asl->i.X0_ = (real *) M1alloc( asl->i.n_var_ * sizeof(real) );
        
        if( !asl->i.Urhsx_ || !asl->i.Uvx_ || !asl->i.havex0_ || !asl->i.X0_)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTMEMERROR;
            #endif
            code = MIP_MEMORY_ERROR;
            goto desallocate_memory;
        }
        
        //making only the non-linear constraints be evaluated at conval and jacval functions
        asl->i.n_conjac_[0] = 0;
        asl->i.n_conjac_[1] = asl->i.nlc_;
        asl->i.congrd_mode = 2;
    
        pfgh_read_ASL(asl, nl, 0);
        
        maxVals = 0;
        
        if( !quadConstrs )
        {
            //counting non-zeros of gradients of nonlinear constraints... 
            aux = 0;
            for(i = 0; i < asl->i.nlc_; i++)
            {
                for(cglc = Cgrad[i]; cglc; cglc = cglc->next)
                    aux++;
            }
        
            maxVals = aux;
        }
        
        {
            //counting nonzeros of hessian of constraints...
            fint aux = quadObj ? -1 : 0; //we only consider objective function if it is nonlinear...
            
            r= quadConstrs ? 0 : 1; //we only consider constraints if it is nonlinear
            
            aux = sphsetup(aux, 1, r, 1);
            //std::cout << "aux: " << aux << "\n";
            //MIP_getchar();
            
            maxVals = MIP_max(aux, maxVals);
            
            if(maxVals > 0)
            {
                MIP_malloc(iRow, maxVals); //iRow = (int *) malloc( maxVals * sizeof(int) );
                MIP_malloc(jCol, maxVals); //jCol = (int *) malloc( maxVals * sizeof(int) );
                
                if( !iRow || !jCol )
                {
                    #if MIP_DEBUG_MODE
                        MIP_PRINTMEMERROR;
                    #endif
                    code = MIP_MEMORY_ERROR;
                    goto desallocate_memory;
                }
            }
        }
        
        //seting the structure of hessian of lagrangian  
        r = 0;
        for(i = 0; i < n; i++)
        {
            for( j = sputinfo->hcolstarts[i]; j < sputinfo->hcolstarts[i+1] ; j++ )
            {
                ////the hessian is stored using only the upper triangle in AMPL. (That is manual says, but it is not true (actually is lower triangle)... anyway sparse matrix is capable to invert indexes if col is greater than row)
                
                iRow[r] = i;
                jCol[r] = sputinfo->hrownos[j];
                
                #if MIP_AMPL_DEBUG_MODE
                    assert( iRow[r] >= jCol[r] );
                #endif
                
                
                r++;
            }
        }
        
        
        /*#if MIP_AMPL_DEBUG_MODE
            if(aux != aux2)
            {
                printf("Error in hessian of constraints coefficients counting!\n");
                getchar();
            }
        #endif */
        
        ret = prob.setLagrangianHessianStructure(r, iRow, jCol);
        if(ret != 0)
        {
            MIP_PRINTERRORNUMBER(ret);
            //printf("Error at seting structure of hessian of constraints. return: %d\n", aux);
            MIP_getchar();
            
            code = ret;
            goto desallocate_memory;
        }
        
        
        if( !quadConstrs )
        {
            aux = 0;
            for(i = 0; i < asl->i.nlc_; i++)
            {
                for(cglc = Cgrad[i]; cglc; cglc = cglc->next)
                {
                    iRow[aux] = i;//newIndC[i] -mrqc -mlin;
                    jCol[aux] = cglc->varno;
                    aux++;
                }
            }
            
            r = prob.setJacobianStructure(aux, iRow, jCol);
            if(r != 0)
            {
                MIP_PRINTERRORNUMBER(r);
                //printf("Error at seting jacobian structure of constraints. return: %d\n", aux2);
                
                MIP_getchar();
            
                code = r;
                goto desallocate_memory;
            }
        }
        
        if( !quadObj )
            prob.setObjNonLinearTermFlag(true);
    }
    
    
    
    #if MIP_AMPL_DEBUG_MODE
        
        if( prob.getNumberOfVars() <= 50 )
            prob.print();
        /*printf("----------------------------------------------------------------------");
        
        printf("c:\n");
        for(int i = 0; i < n; i++)
            printf("%f ", prob.c[i]);
        printf("\n");
        
        printf("Q:\n");
        prob.Q.printSparseMatrix();
        
        
        printf("A:\n");
        prob.A.printSparseMatrix();
        
        printf("lc:\n");
        for(int i = 0; i < m; i++)
        printf("%f ", prob.lc[i]);
        printf("\n");
        
        printf("uc:\n");
        for(int i = 0; i < m; i++)
        printf("%f ", prob.uc[i]);
        printf("\n");
        
        
        printf("lx:\n");
        for(int i = 0; i < n; i++)
        printf("%f ", prob.lx[i]);
        printf("\n");
        
        printf("ux:\n");
        for(int i = 0; i < n; i++)
        printf("%f ", prob.ux[i]);
        printf("\n");
        
        
        printf("Jacobian of constraints:\n");
        prob.J.printSparseMatrix();
        
        
        printf("Hessian of lagrangian:\n");
        prob.lagH.printSparseMatrix();
        
        
        printf("Integer variables:\n");
        for(int i = 0; i < n; i++)
        {
        if(prob.xtype[i] == MIP_VT_INTEGER)
            printf("%d ", i);
        }
        printf("\n"); */
    #endif
    
    
    
    code = 0;
    
desallocate_memory:
    
    
    #if MIP_SAVE_OUTPUT_FILE
        if(outFile)
            fclose(outFile);
    #endif
    
    //if(myEval) delete myEval;
    
    if(nzqp)	free(nzqp);
    if(rowqp)	free(rowqp); // the vectors was allocated in ALS
    if(colqp)	free(colqp); // the vectors was allocated in ALS
    
    if(iRow)	free(iRow);
    if(jCol)	free(jCol);
    
    if(delsqp)	free(delsqp); // the vectors was allocated in ALS
    
    if(vals)	free(vals);
    
    
    
    return code;
    
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif




void MIP_ReadAmplModel::putInfoToAmpl(const char* msg, const int return_code, double* primalSol, double* dualSol)
{
#if MIP_HAVE_ASL
    solve_result_num = return_code;
    
    write_sol_ASL(asl, msg, primalSol, dualSol, NULL);
#endif
}


void MIP_ReadAmplModel::readParameters(char** argv, Option_Info* opinfo)
{
#if MIP_HAVE_ASL
    getopts_ASL(asl, argv, opinfo);
#endif
}




/*#if MIP_HAVE_ASL
extern "C" { //these definitions are necessary to run ASL in multithread environment. We must define semaphores to ASL.
    
    MIP_Mutex *WAXM_MIP_amplMutexes;
    MIP_Mutex WAXM_MIP_amplMutex;
    
    
void ACQUIRE_DTOA_LOCK(unsigned int n)
{
    WAXM_MIP_amplMutex.lock();
}


void FREE_DTOA_LOCK(unsigned int n)
{
    WAXM_MIP_amplMutex.unlock();
}


int dtoa_get_threadno(void)
{
    return 1;
}


void init_dtoa_locks(void)
{
    
}


//void set_max_dtoa_threads(unsigned int);
    
};
#endif */











