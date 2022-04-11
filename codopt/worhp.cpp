

#include <cmath>
#include <cstdlib>
#include <climits>
#include <cassert>
#include <iostream>
#include <new>


#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"


//using namespace std;
using namespace optsolvers;
using namespace newspm;
using namespace minlpproblem;





OPT_Worhp::OPT_Worhp():OPT_MyNLPSolver()
{
    initialize();
}


            
OPT_Worhp::~OPT_Worhp()
{
    deallocateMemory();
}



void OPT_Worhp::deallocateSolverEnv()
{
/* #if OPT_HAVE_WORHP
    WorhpFree(&opt, &wsp, &par, &cnt);
#endif */
    OPT_MyNLPSolver::deallocateSolverEnv();
}



// __methods from Solver __
void OPT_Worhp::deallocateMemory()
{
    OPT_MyNLPSolver::deallocateMemory();
}



bool OPT_Worhp::getMinusLambdaOnLagran()
{
    return false;
}


int OPT_Worhp::getNumberOfIterations(long unsigned int &niter)
#if OPT_HAVE_WORHP
{
    niter = wsp.MajorIter;
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


OPT_LISTSOLVERS OPT_Worhp::getSolverCode()
{
    return optsolvers::OPT_WORHP;
}



int OPT_Worhp::getVariableType( const int index, OPT_VARTYPE &varType )
{
    varType = optsolvers::OPT_VT_CONTINUOUS;
    
    return 0;
}



void OPT_Worhp::initialize()
{
/* #if OPT_HAVE_WORHP
    opt.initialised = false;
    wsp.initialised = false;
    par.initialised = false;
    cnt.initialised = false;
#endif */
    
    //nzJac = 0;
    //nzHess = 0;
    
    OPT_MyNLPSolver::initialize();
}



int OPT_Worhp::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
#if OPT_HAVE_WORHP
{
    int status;
    
    
    // Properly zeros everything, or else the following routine could get confused.
    WorhpPreInit(&opt, &wsp, &par, &cnt);
    
    
    /* Make sure so the initialisers know they have to work! */
    //opt.initialised = false;
    //wsp.initialised = false;
    par.initialised = false;
    //cnt.initialised = false;
    
    
    /*
    * Parameter initialisation routine that must be called 
    * when using ReadParamsNoInit instead of ReadParams.
    */ 
    InitParams(&status, &par);
    
    
    
    //ReadParamsNoInit(&status, (char*)OPT_WORHP_PARAM_FILE_NAME, &par);
    
    #if OPT_DEBUG_MODE
        assert(par.initialised == true);
    #endif
    
    //printf("par.NLPprint: %d\n", par.NLPprint);
        
    par.MatrixCC = false;
    par.Timeout = INFINITY;
    par.LogLevel = 0;
    par.LogResult = 0;
    par.NLPprint = -1;
    par.MaxIter = 1000;
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Worhp::setMaxCPUTime(const double time)
#if OPT_HAVE_WORHP
{
    par.Timeout = time;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif





int OPT_Worhp::setNumberOfThreads(const int nthreads)
#if OPT_HAVE_WORHP
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Worhp::setOutputLevel(const int level)
#if OPT_HAVE_WORHP
{
    par.NLPprint = level-1;
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Worhp::setRelativeDualTol(const double tol)
#if OPT_HAVE_WORHP
{
    return OPT_OPERATION_NOT_IMPLEMENTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Worhp::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_WORHP
{
    par.AcceptTolOpti = tol;
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Worhp::setRelativePrimalTol( const double tol )
#if OPT_HAVE_WORHP
{
    par.AcceptTolFeas = tol;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Worhp::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_WORHP
{
    bool r;
    
    r = WorhpSetDoubleParam(&par, param, value);
    
    if( !r )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        printDblParamErrorMsg(!r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Worhp::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_WORHP
{
    bool r;
    
    r = WorhpSetIntParam(&par, param, value);
    
    if( !r )
    {
        r = WorhpSetBoolParam(&par, param, value);
        if( !r )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            printIntParamErrorMsg(!r, param, value);
            
            return OPT_BAD_INPUT;
        }
    }
    
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Worhp::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_WORHP
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Worhp::setVariableType( const int index, const OPT_VARTYPE varType )
#if OPT_HAVE_WORHP
{
    //anyway we set var type in MIP_MINLPProb...
    OPT_MyNLPSolver::setVariableType(index, varType);
    
    
    //worhp is a continuous solver...
    if( varType == OPT_VT_INTEGER )
        return OPT_OPERATION_NOT_SUPPORTED;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Worhp::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_WORHP
{
    const bool evalLag = prob.hasNlObj || prob.hasNlConstrs;
    const bool constrQ = prob.getNumberOfQuadMatricesInConstrs() > 0;
    
    bool newx;
    //bool reinit = false;
    bool reSetJacIndex = constrsBndsChg; //worhp needs indices in column order. So, if bounds changes, maybe some contsraints turn free and we have to reset the column order.
    bool reSetAllJac = false;
    bool reSetHessIndex = false, reSetAllHess = false;
    bool *constrEval = (bool *) auxIndex ;
    int n, m;
    int r, w; //, nzJac = 0, nzHess = 0;
    int newm = 0;
    unsigned int nzJac;
    double *X, *G, *L, *U, *Mu;
    double *auxConstr = auxValues;
    
    unsigned int mquad = 0;
    int *quadIndex = nullptr;
    
    const bool *nlConstr = prob.nlConstr; 
    const double *lx = prob.lx, *ux = prob.ux;
    const double *lc = prob.lc, *uc = prob.uc;
    
    //const bool *nlConstr = prob.nlConstr;
    //const int mq = prob.getNumberOfQuadMatricesInConstrs();
    const MIP_SparseMatrix &Q = prob.Q;
    const MIP_SparseMatrix &A = prob.A;
    const MIP_SparseMatrix *QC= prob.QC;
    MIP_SparseMatrix &J = prob.J;
    MIP_SparseMatrix &lagH = prob.lagH;
    
    //OptVar    opt;
    //Workspace wsp;
    //Control   cnt;
    //Params    par2;
    
    Params    dummyPar; //that is ridicuous, but we need pass a parameter structure to worhpFree... so we just create this for it
    
    
    opt.initialised = false;
    wsp.initialised = false;
    cnt.initialised = false;
    
    
    InitParams(&w, &dummyPar);
    
    
    getNumberOfVars(n);
    getNumberOfConstraints(m);
    
    if( resetSol )
    {
        this->resetSol();
    }
    else
    {
        feasSol = false;
    }
    
    
    if( nmChg )
    {
        deallocateAuxDerivativeIndexStructures();
        allocateAuxDerivativeIndexStructures();
        
        nmChg = false;
        //reinit = true;
        
        reSetAllJac = true;
        reSetAllHess = true;
        
        //reSetHessIndex = true; //we need set it becuase we have to put diagonal positions in hessian even if they are zero
    }
    
    
    for( int i = 0; i < m; i++ )
    {
        if( reSetAllJac || constrChg[i] )
        {
            r = setuJacStructRow(i);
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                retCode = OPT_MEMORY_ERROR;
                goto termination;
            }
            
            constrChg[i] = false;
            
            reSetJacIndex = true;
        }
    }
    
    
    for( int i = 0; i < n; i++ )
    {
        if( reSetAllHess || hessChg[i] || rowQuadConstrChg[i] || quadObjChg[i] )
        {
            r = setuHessStructRow(i);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            hessChg[i] = false;
            rowQuadConstrChg[i] = false;
            quadObjChg[i] = false;
            
            reSetHessIndex = true;
        }
    }
    
    
    
    if( reSetJacIndex )
    {
        int r = setuJacIndices();
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        //std::cout << "uJac: " << std::endl;
        //uJac.printSparseMatrix();
        
        constrsBndsChg = false;
    }
    
    
    if( reSetHessIndex )
    {
        //first we check if all lines really has the diagonal position (we could have a row empty in hessian. Note if row is not empty, so, the diagonal position is over there due to setuHessStructRow)
        
        for(int i = 0; i < n; i++)
        {
            const unsigned int nel = uHess.getNumberOfElementsAtRow(i);
            
            if( nel == 0 )
            {
                r = uHess.setRowStructure(i, 1, &i);
                OPT_IFMEMERRORGOTOLABEL(r, retCode, termination);
                
            }
        }
        
        r = setuHessIndices();
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        //std::cout << "uHess: " << std::endl;
        //uHess.printSparseMatrix();
        
        //OPT_getchar();
    }
    
    
    
    
    
    
    //if( reinit )
    {
        mat_int *row, *col;
        
        /* Make sure so the initialisers know they have to work! */
        opt.initialised = false;
        wsp.initialised = false;
        dummyPar.initialised = false;
        cnt.initialised = false;
        
        
        for(int i = 0; i < m; i++)
        {
            if( lc[i] > -OPT_INFINITY || uc[i] < OPT_INFINITY )
            {
                newm++;
                constrEval[i] = true;
                
                if( QC[i].getNumberOfElements() > 0 )
                    mquad++;
                
            }
            else
            {
                constrEval[i] = false;
            }
        }
        
        if( mquad > 0 )
        {
            OPT_malloc(quadIndex, mquad);
            OPT_IFMEMERRORGOTOLABEL( !quadIndex, retCode, termination );
            
            int k = 0;
            for(int i = 0; i < m; i++)
            {
                if( constrEval[i] && QC[i].getNumberOfElements() > 0 )
                {
                    quadIndex[k] = i;
                    k++;
                }
            }
            
            #if OPT_DEBUG_MODE
                assert( k == mquad );
            #endif
        }
        
        
        //Specify number of variables and constraints.
        opt.n = n;
        opt.m = newm;
        
        
        //Specify nonzeros of derivative matrices.
        wsp.DF.nnz = n;
        wsp.HM.nnz = uHess.getNumberOfElements();
        
        //to calculate number of nonzeros in jacobian, we should discard free constraints
        nzJac = 0;
        for(int i = 0; i < m; i++)
        {
            if( constrEval[i] )
                nzJac += uJac.getNumberOfElementsAtRow(i);
        }
        
        wsp.DG.nnz = nzJac;
        
        
        
        
        if( newm > 0 )
        {
            if( std::isnan( lambdaInit[0] ) )
                par.InitialLMest = true;
        }
        
        
        // Data structure initialisation. We have to do it every time we call worhp...
        
        WorhpInit(&opt, &wsp, &par, &cnt);
        
        if( cnt.status != FirstCall )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(cnt.status);
            #endif
            
            retCode = OPT_SOLVER_ERROR;
            goto termination;
        }
        
        //setting objective type
        if( prob.objFactor == 0 )
        {
            opt.FType = WORHP_CONSTANT;
        }
        else
        {
            if( prob.hasNlObj )
            {
                opt.FType = WORHP_NONLINEAR;
            }
            else
            {
                if( Q.getNumberOfElements() > 0 )
                {
                    opt.FType = WORHP_QUADRATIC;
                }
                else
                {
                    bool hasObjLinCoef = prob.hasLinCoefObj();
                    /*const double *c = prob.c;
                    
                    for(int i = 0;i < n; i++)
                    {
                        if(c[i] != 0.0)
                        {
                            hasObjLinCoef = true;
                            break;
                        }
                    } */
                    
                    opt.FType = hasObjLinCoef ? WORHP_LINEAR : WORHP_CONSTANT ;
                }
            }
        }
        
        
        par.MaxIter = 1000; //that is reaaly strange, but worhp only accepts change this number here...
        
        //WorhpInit(&opt, &wsp, &par, &cnt);
        
        //Specify matrix structures in CS. We have to use Fortran indexing
        
        //Gradient of f
        row = wsp.DF.row;
        for( int i = 0; i < n; i++ )
            row[i] = i + 1; //fortran indices
        
        
        //Gradient of g
        
        if( uJac.getNumberOfElements() > 0 )
        {
            row = wsp.DG.row;
            col = wsp.DG.col;
            
            w = 1; //we start w from 1instead of 0 because worhp adotps fortran indixing
            for(int i = 0; i < m; i++)
            {
                if( lc[i] <= -OPT_INFINITY && uc[i] >= OPT_INFINITY )
                    continue;
                
                const unsigned int nel = uJac.getNumberOfElementsAtRow(i);
                const int* const rcols = uJac[i];
                const unsigned int* const rvalues = uJac(i);
                
                for( unsigned int j = 0; j < nel; j++ )
                {
                    const unsigned int ind = rvalues[j];
                    
                    row[ind] = w; //fortran index (w already starts from 1) :(
                    col[ind] = rcols[j] + 1; //fortran index
                }
                
                w++;
            }
            
            #if OPT_DEBUG_MODE
                assert(w == newm + 1);
            #endif
            
            
            /*std::cout << "\njac: " << std::endl;
            uJac.printSparseMatrix();
            for(unsigned int j = 0; j < uJac.getNumberOfElements(); j++)
                std::cout << "row["<<j<<"]: " << row[j] << " col["<<j<<"]: " << col[j] << std::endl; */
        }
        
        
        if( uHess.getNumberOfElements() > 0 )
        {
            row = wsp.HM.row;
            col = wsp.HM.col;
            
            for(int i = 0; i < n; i++)
            {
                
                const int *hcols = uHess[i];
                const unsigned int *hvalues = uHess(i);
                
                const unsigned int nel = uHess.getNumberOfElementsAtRow(i);
                for(unsigned int j = 0; j < nel; j++)
                {
                    const unsigned int ind = hvalues[j];
                    
                    col[ind] = hcols[j] + 1; //fortran index :(
                    row[ind] = i+1; //fortran index :(
                }
                
            }
            
            /*std::cout << "\nhess: \n";
            uHess.printSparseMatrix();
            for(unsigned int j = 0; j < uHess.getNumberOfElements(); j++)
                std::cout << "row["<<j<<"]: " << row[j] << " col["<<j<<"]: " << col[j] << std::endl; */
        }
        
    }
    
    
    //Abbreviate notation
    X = opt.X;
    G = opt.G;
    Mu= opt.Mu;
    
    //setting variable and constraint bounds. We are sorry, but we have to set all every time solve is called. I do not want put flags to sinalize bounds changing...
    
    
    L = opt.GL;
    U = opt.GU;
    w = 0;
    for(int i = 0; i < m; i++)
    {
        double lb = lc[i], ub = uc[i];
        //r = getConstraintBounds(i, lb, ub);
        
        if( !constrEval[i] )
        {
            continue; //free constraint
        }
        
        
        L[w] = lb <= -OPT_INFINITY ? -par.Infty : lb;
        U[w] = ub >=  OPT_INFINITY ?  par.Infty : ub;
        
        if( nlConstr[i] )
        {
            opt.GType[w] = WORHP_NONLINEAR;
        }
        else
        {
            if( QC[i].getNumberOfElements() > 0 )
                opt.GType[w] = WORHP_QUADRATIC;
            else
            {
                opt.GType[w] = A.getNumberOfElementsAtRow(i) > 0 ? WORHP_LINEAR : WORHP_CONSTANT ;
            }
        }
        
        w++;
    }
    
    #if OPT_DEBUG_MODE
        assert( w == opt.m );
    #endif
    
    
    L = opt.XL;
    U = opt.XU;
    for(int i = 0; i < n; i++)
    {
        L[i] = lx[i] <= -OPT_INFINITY ? -par.Infty : lx[i];
        U[i] = ux[i] >=  OPT_INFINITY ?  par.Infty : ux[i];
    }
    
    
    if( n > 0 )
    {
        double *X = opt.X;
        
        if( std::isnan( xInit[0] ) )
        {
            //setting initial solution
            for(int i = 0; i < n; i++)
            {
                X[i] = L[i] > -par.Infty ? L[i] : (U[i] < par.Infty ? U[i] : 0) ;
            }
        }
        else
        {
            OPT_copyArray( n, xInit, X );
        }
        
        if( std::isnan( zInit[0] ) )
            OPT_setAllArray(n, opt.Lambda, 0.0);
        else
            OPT_copyArray( n, zInit, opt.Lambda );
        
    }
    
    if( m > 0 )
    {
        double *Mu = opt.Mu;
        
        if( std::isnan( lambdaInit[0] ) )
        {
            OPT_setAllArray( opt.m, opt.Mu, 0.0);
        }
        else
        {
            
            if( m == opt.m )
            {
                //no free constraints
                OPT_copyArray( m, lambdaInit, Mu );
            }
            else
            {
                w = 0;
                for(int i = 0; i < m; i++)
                {
                    if( lc[i] > -OPT_INFINITY || uc[i] < OPT_INFINITY )
                    {
                        Mu[w] = lambdaInit[i];
                        w++;
                    }
                }
                
                #if OPT_DEBUG_MODE
                    assert(w == opt.m);
                #endif
            }
        }
    }
    
    
   
    
    
    
    //nzJac = wsp.DG.nnz; // in nonfirst calls, nzJac can be set as 0 
    //nzHess = wsp.HM.nnz; // in nonfirst calls, nzHess can be set as 0 
    
    
    //for(int i = 0; i < n; i++)
        //cout << "opt.X[" << i << "]: " << opt.X[i] << endl;
    
    
    
    while( cnt.status < TerminateSuccess && cnt.status > TerminateError )
    {
        if( GetUserAction(&cnt, callWorhp) )
        {
            Worhp(&opt, &wsp, &par, &cnt);
            //No DoneUserAction!
        }
        
        
        if( GetUserAction(&cnt, iterOutput) )
        {
            IterationOutput(&opt, &wsp, &par, &cnt);
            DoneUserAction(&cnt, iterOutput);
        }
        
        if( opt.newX)
            newx = true;
        
        if( GetUserAction(&cnt, evalF) )
        {
            if( wsp.ScaleObj != 0.0 )
            {
                r = prob.objEval(threadNumber, true, X, opt.F, in_nl_obj_factor);
                
                if( prob.hasNlObj )
                    newx = false;
                
                if( r != 0 )
                {
                    #if OPT_DEBUG_MODE
                        std::cerr <<  OPT_PREPRINT "Callback function error " << r << OPT_GETFILELINE << "\n";
                    #endif
                    
                    //retCode = OPT_CALLBACK_FUNCTION_ERROR;
                    //goto termination;
                    
                    opt.F = NAN;
                }
            }
            
            opt.F *= wsp.ScaleObj;
            
            DoneUserAction(&cnt, evalF);
        }
        
        
        if( GetUserAction(&cnt, evalG) )
        {
            double *pG = newm == m ? G : auxConstr ;
            
            
            r = prob.constraintsEval(threadNumber, newx, constrEval, X, pG);
            
            if( prob.hasNlConstrs )
                newx = false;
            
            if(r == 0)
            {
                if( newm != m )
                {
                    w = 0;
                    for(int i = 0; i < m; i++)
                    {
                        if( constrEval[i] )
                        {
                            G[w] = pG[i];
                            w++;
                        }
                    }
                    
                    #if OPT_DEBUG_MODE
                        assert(pG != G);
                        assert(w == newm);
                    #endif
                }
            }
            else //( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    std::cerr <<  OPT_PREPRINT << "Callback function error " << r << OPT_GETFILELINE << "\n";
                #endif
                
                //retCode = OPT_CALLBACK_FUNCTION_ERROR;
                //goto termination;
                OPT_setAllArray<double>(newm, G, NAN);
            }
            
            DoneUserAction(&cnt, evalG);
        }
        
        
        if( GetUserAction(&cnt, evalDF) )
        {
            const double f = wsp.ScaleObj;
            double *grad = wsp.DF.val;
            
            r = prob.objGradEval( threadNumber, newx, X, grad, in_nl_obj_factor );
            
            if( prob.hasNlObj )
                newx = false;
            
            if( r == 0 )
            {
                if( f != 1.0 )
                {
                    OPT_multiplyAllArray(n, f, grad);
                    //for( int i = 0; i < n; i++ )
                        //grad[i] *= f;
                }
            }
            else
            {
                #if OPT_DEBUG_MODE
                    std::cerr << OPT_PREPRINT "Callback function error " << r << OPT_GETFILELINE << "\n";
                #endif
                
                //retCode = OPT_CALLBACK_FUNCTION_ERROR;
                //goto termination;
                OPT_setAllArray<double>( n, grad, NAN );
            }
            
            
            DoneUserAction(&cnt, evalDF);
        }
        
        
        if( GetUserAction(&cnt, evalDG) )
        {
            double *values = wsp.DG.val;
            double *grad = auxValues;
            int r = 0;
            
            //OPT_setAllArray<double>( nzJac, values, 0.0 );
            
            if( prob.hasNLConstraints() )
            {
                r = prob.nlJacobianEval( threadNumber, newx, constrEval, X, J);
                
                newx = false;
                
                if( r != 0 )
                {
                    #if OPT_DEBUG_MODE
                        std::cerr <<  OPT_PREPRINT << "Callback function error " << r << OPT_GETFILELINE << "\n";
                    #endif
                    
                    //retCode = OPT_CALLBACK_FUNCTION_ERROR;
                    //goto termination;
                        
                    OPT_setAllArray<double>( nzJac, values, NAN );
                }
            }
            
            if( r == 0 )
            {
                for( int i = 0; i < m; i++ )
                {
                    if( !constrEval[i] )
                        continue;
                    
                    const auto rnel = uJac.getNumberOfElementsAtRow(i);
                    const int* jcols = uJac[i];
                    const unsigned int* jvals = uJac(i);
                    
                    if(  !( nlConstr[i] && QC[i].getNumberOfElements() == 0 && A.getNumberOfElementsAtRow(i) == 0 )   ) //if constraint has only the nonlinear term, we do not need initialize with zeros because we do not accumulate several structures in gradient
                    {
                        for( unsigned int p = 0; p < rnel; p++ )
                            grad[ jcols[p] ] = 0.0;
                    }
                    
                    MIP_constrCompleteGrad( prob, J, i, X, grad, false );
                    
                    /*unsigned int k = 0;
                    
                    for( unsigned int j = 0; j < (unsigned int) n; j++ )
                    {
                        //columns are ordered in uJac
                        if( grad[j] != 0.0 )
                        {
                            while( jcols[k] != j ) //( jrow[k].getColumn() != j )
                                k++; //I hope uJac has the non zeros index correctly in order. Otherwise, we will have serious problems here...
                            
                            values[ jvals[k] ] = grad[j]; //values[ jrow[k].getValue() ] = grad[j];
                            
                            k++; //we can walk one position in k
                        }
                    }
                    
                    #if OPT_DEBUG_MODE
                        //assert( k <= jrow.getNumberOfElements() );
                        assert( k <= uJac.getNumberOfElementsAtRow(i) );
                    #endif */
                        
                    for( unsigned int p = 0; p < rnel; p++ )
                    {
                        values[ jvals[p] ] = grad[ jcols[p] ];  //jvals has the index of position in the array values
                    }
                }
            }
            
            
            DoneUserAction(&cnt, evalDG);
        }
        
        
        
        if( GetUserAction(&cnt, evalHM) )
        {
            const double obj_factor = wsp.ScaleObj;
            const double objf = obj_factor* prob.objFactor;
            const bool objQ = Q.getNumberOfElements() > 0 && objf != 0.0;
            int r = 0;
            double *values = wsp.HM.val;
            double *line = auxValues;
            double *oMu = opt.Mu;
            double *pMu;
            
            
            if( newm < m )
                OPT_setAllArray(uHess.getNumberOfElements(), values, 0.0);  //we have free constraints
            
            
                
            if( newm != m )
            {
                w = 0;
                for(int j = 0; j < m; j++)
                {
                    if( constrEval[j] )
                    {
                        auxConstr[j] = oMu[w];
                        w++;
                    }
                    else
                    {
                        auxConstr[j] = 0.0;
                    }
                }
                
                pMu = auxConstr;
                
                #if OPT_DEBUG_MODE
                    assert( w == newm );
                #endif
            }
            else
            {
                pMu = oMu;
            }
            
            
            
            if( evalLag )
            {
                
                
                //prob.objFactor is already considered in prob->nlpHessianEval...
                r = prob.nlpHessianEval(threadNumber, newx, X, prob.hasNlObj ? obj_factor*in_nl_obj_factor : 0.0, pMu, lagH);
                
                newx = false;
                
                if( r != 0 )
                {
                    #if OPT_DEBUG_MODE
                        std::cerr <<  OPT_PREPRINT "Callback function error " << r << OPT_GETFILELINE << "\n";
                    #endif
                    
                    //retCode = OPT_CALLBACK_FUNCTION_ERROR;
                    //goto termination;
                    
                    OPT_setAllArray<double>(wsp.HM.nnz , values, NAN );
                }
                
            }
            
            
            if(r == 0)
            {
                for( auto rowit = uHess.beginRowIndex() ; *rowit < n; ++rowit )
                {
                    int i = *rowit;
                    
                    const unsigned int ncols = uHess.getNumberOfElementsAtRow(i);
                    const int* hcol = uHess[i];
                    const unsigned int* hvalue = uHess(i);
                    
                    for(int j = 0; j < ncols; j++)
                        line[ hcol[j] ] = 0.0;
                    
                    //prob.objFactor is already considered inside this function
                    MIP_completeLagHessianRow(prob, mquad, quadIndex, lagH, obj_factor, pMu, i, line, false);
                    
                    for(int j = 0; j < ncols; j++)
                    {
                        values[ hvalue[j] ] = line[ hcol[j] ];
                    }
                }
                
                
                /* for(int  i = 0; i < n; i++)
                {
                    OPT_setAllArray( i+1, line, 0.0 ); //we just use lower triangle
                    
                    if( evalLag )
                        lagH.getFullRow(i, line, true, false, 1); //lagH.getFullRowAccumulation( i, grad );
                    
                    if( objQ )
                        Q.accumulateRowInArray(i, line, objf); //OPT_accumulateRowInArray( Q, i, objf, grad );
                    
                    if( constrQ )
                    {
                        
                        w = 0;
                        for(int j = 0; j < m; j++)
                        {
                            if( !constrEval[j] )
                                continue;
                            
                            if( QC[j].getNumberOfElements() > 0 && oMu[w] != 0 )
                                QC[j].accumulateRowInArray(i, line, oMu[w]); 
                            
                            w++;
                        }
                        
                        #if OPT_DEBUG_MODE
                            assert( w == newm );
                        #endif
                    }
                    
                    
                    
                    const int* hcol = uHess[i];
                    const unsigned int* hvalue = uHess(i);
                    
                    unsigned int w = 0;
                    for(unsigned int k = 0; k <= (unsigned int) i; k++) //just lower triangle
                    {
                        //columns are ordered in uHess
                        if( line[k] != 0.0 )
                        {
                            while( hcol[w] != k ) //while( hrow[w].getColumn() != k )
                                w++; //I hope uHess has the non zeros index correctly in order. Otherwise, we will have serious problems here...
                            
                            values[ hvalue[w] ] = line[k]; //values[ hrow[w].getValue() ] = grad[k];
                            w++; //we can walk one position in w
                        }
                    }
                    
                    #if OPT_DEBUG_MODE
                        assert( w <=  uHess.getNumberOfElementsAtRow(i) );
                    #endif
                } */
            }
            
            
            DoneUserAction(&cnt, evalHM);
        }
        
        
        
        
        #if OPT_DEBUG_MODE
            assert( !GetUserAction(&cnt, fidif) ); // I think we do not need call function to finite difference...
        #endif
        
    }
    
    

    
    
    
    
    if( storeSol )
    {
        OPT_copyArray(n, opt.X, sol);
    }
    
    if( storeConstrs )
    {
        
        if( opt.m == m ) //no free constraints
        {
            OPT_copyArray(m, G, constr);
        }
        else
        {
            bool *ceval = (bool *) auxIndex2;
            double *auxConstr = auxValues2;
            
            X = opt.X;
            
            for(int i = 0; i < m; i++)
                ceval[i] = !constrEval[i];
            
            int r = prob.constraintsEval( threadNumber, true, ceval, X, auxConstr );
            
            if(r != 0)
                OPT_PRINTERRORMSG("warning: evaluation error on redundant constraints");
            
            w = 0;
            for(int i = 0; i < m; i++)
            {
                if( constrEval[i] )
                {
                    constr[i] = G[w];
                    w++;
                }
                else
                {
                    constr[i] = auxConstr[i];
                }
            }
            
            #if OPT_DEBUG_MODE
                assert( w == opt.m );
            #endif
        }
        
        
    }
    
    if( storeDualSol )
    {
        double *Mu = opt.Mu;
        
        if( m == opt.m )
        {
            OPT_copyArray(m, Mu, dualSolC);
        }
        else
        {
            w = 0;
            for(int i = 0; i <m; i++)
            {
                if( constrEval[i] )
                {
                    dualSolC[i] = Mu[w];
                    w++;
                }
                else
                    dualSolC[i] = 0.0;
            }
            
            #if OPT_DEBUG_MODE
                assert( w == opt.m );
            #endif
        }
        
        //so, worhp only define lambda as n dimensional array.
        OPT_copyArray(n, opt.Lambda, dualSolV);
    }
    
    //objValue = opt.F;
    
    //worhp is a really crazy solver. We are having problems about final objective value. So, we perform a function evaluation to try get the correct value 
    r = prob.objEval(threadNumber, true, sol, objValue, in_nl_obj_factor);
    if( r != 0 )
    {
        objValue = NAN;
        
        #if OPT_DEBUG_MODE
            std::cerr << OPT_PREPRINT  "Error " << r << OPT_GETFILELINE << "\n";
        #endif
        
        retCode = OPT_CALLBACK_FUNCTION_ERROR;
        goto termination;
    }
    
    
    
    if( prob.objFactor < 0 )
    {
        objValue = -objValue;
        dualObjValue = -dualObjValue;
    }
    
    
    
    origSolverRetCode = cnt.status;
    
    
    switch( cnt.status )
    {
        case OptimalSolution:
        case LowPassFilterOptimal:
        case AcceptableSolution:
        case OptimalSolutionConstantF:
        case NotDiffable:
        case AcceptableSolutionSKKT:
            retCode = OPT_OPTIMAL_SOLUTION;
            feasSol = true;
            break;
            
        case GlobalInfeas:
        case LocalInfeas:
        case LocalInfeasOptimal:
            retCode = OPT_INFEASIBLE_PROBLEM;
            break;
            
        case MaxIter:
            #if OPT_PRINT_MAX_ITER_WARNING
                if( numberOfWarningsByIterLimit < maxNumberOfWarningsByIterLimit )
                {
                    std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Worhp solving!\n";
                    numberOfWarningsByIterLimit++;
                    
                    if( numberOfWarningsByIterLimit == maxNumberOfWarningsByIterLimit )
                        std::cerr << OPT_PREPRINT "Warning: Maximum number of warnings by maximum iteration achieved! Stopping these warnings.\n";
                }
            #endif
            retCode = OPT_MAX_ITERATIONS;
            break;
            
        case Timeout:
            retCode = OPT_MAX_TIME;
            break;
            
        case Unbounded:
            retCode = OPT_UNBOUNDED_PROBLEM;
            break;
            
        case FeasibleSolution:
        //case StationaryPointFound:
        case SearchDirectionSmall:
        case SearchDirectionZero:
            retCode = OPT_FEASIBLE_SOLUTION;
            feasSol = true;
            break;
            
        case LicenseError:
            OPT_PRINTERRORMSG("License error at worhp solver. Check your worhp license!");
            retCode = OPT_LICENSE_ERROR;
            break;
        
        default:
            retCode = OPT_UNDEFINED;
    }
    
    if(newm == m)
    {
        //for(int i = 0; i < newm; i++)
            //printf("G[%d]: %f\n", i, G[i]);
        
        #if OPT_DEBUG_MODE
            //printf("status: %d\n", cnt.status);
            feasSol = prob.isConstrValuesFeasible(1.0e-5, 1.0e-5, G);
            if( retCode == OPT_OPTIMAL_SOLUTION || retCode == OPT_FEASIBLE_SOLUTION  )
            {
                assert( feasSol );
            }
        #endif
    }
    
termination:
    
    WorhpFree(&opt, &wsp, &dummyPar, &cnt);
    
    if(quadIndex) free(quadIndex);
    
    
    return retCode;
    
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Worhp::warmUp()
#if OPT_HAVE_WORHP
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



// __methods from NLPSolver __




// __ methods from MyNLPSolver __


int OPT_Worhp::allocateAuxDerivativeIndexStructures( )
{
    const int &n = prob.n;
    const int &m = prob.m;
    
    int r;
    
    r = uJac.initialize(m, n, false);
    r +=uHess.initialize(n, n, true);
    
    //r = uJac.allocateSparseRows(m);
    //r+= uHess.allocateSparseRows(n);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        
        return OPT_MEMORY_ERROR;
    }
    
    
    return 0;
}



void OPT_Worhp::deallocateAuxDerivativeIndexStructures()
{
    uJac.desallocateMemory();
    uHess.desallocateMemory();
}




int OPT_Worhp::setuJacStructRow( const int rowIndex )
{
    const unsigned int n = prob.n;
    bool *cols = nullptr;
    const minlpproblem::MIP_SparseMatrix *QC = prob.QC;
    
    int code, r;
    
    
    OPT_malloc(cols, n);
    OPT_IFMEMERRORGOTOLABEL(!cols, code, termination);
    
    
    OPT_setAllArray(n, cols, false);
    
    
    prob.J.getRowStructure(rowIndex, cols, NULL, true ); //prob.J[rowIndex].getStructure( cols );
    prob.A.getRowStructure(rowIndex, cols, NULL, true ); //prob.A[rowIndex].getStructure( cols );
    
    if( QC[rowIndex].getNumberOfElements() > 0 )
    {
        auto &rowQC = QC[rowIndex];
        
        //const unsigned int nqrows = QC[rowIndex].getNumberOfRows();
        
        
        for( auto rowit = rowQC.beginRowIndex() ; *rowit < n ; ++rowit )
        {
            const unsigned int j = *rowit;
            
            int nel;
            rowQC.getRowStructure(j, cols, &nel, true); // unsigned int nel = QC[rowIndex][j].getStructure( cols ); //for each nondiagonal position in QC[i], we have two positions to set in the gradient
            
            #if OPT_DEBUG_MODE
                assert( nel > 0 );
            #endif
            cols[j] = true; //if there is some column in line j, both the column and j must be set as true
        }
        
        
        /*for( unsigned int j = 0; j < nqrows; j++ )
        {
            int nel;
            QC[rowIndex].getRowStructure(j, cols, &nel, true); // unsigned int nel = QC[rowIndex][j].getStructure( cols ); //for each nondiagonal position in QC[i], we have two positions to set in the gradient
            
            if( nel > 0 )
                cols[j] = true; //if there is some column in line j, both the column and j must be set as true
        }*/
        
    }
    
    
    r = uJac.setRowStructure(rowIndex, n, cols);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    
    code = 0;
    
termination:

    if(cols)   free(cols);
    
    return code;
}


int OPT_Worhp::setuHessStructRow( const int rowIndex )
{
    const int &n = prob.n;
    const int &m = prob.m;
    bool *cols = nullptr;
    minlpproblem::MIP_SparseMatrix *QC = prob.QC;
    
    int code, r;
    
    
    
    OPT_malloc(cols, n);
    OPT_IFMEMERRORGOTOLABEL(!cols, code, termination);
    
    
    OPT_setAllArray(n, cols, false);
    
    
    prob.lagH.getRowStructure(rowIndex, cols, NULL, true); //prob.lagH[rowIndex].getStructure(cols);
    
    prob.Q.getRowStructure(rowIndex, cols, NULL, true); //prob.Q[rowIndex].getStructure(cols);
    
    for(int i = 0; i < m; i++)
    { //maybe we should filter the free constraints, but we let them in the sctruture also. In a B&B procedure, maybe constarints turns free and undo that from some iteration to another.
        
        if( QC[i].getNumberOfElements() > 0 )
        {
            QC[i].getRowStructure(rowIndex, cols, NULL, true ); //QC[i][rowIndex].getStructure( cols );
        }
    }
    
    //that is ridiculous, but worphs needs we set diagonal posiition, even if is zero
    
    cols[rowIndex] = true;
    
    
    r = uHess.setRowStructure(rowIndex, n, cols);
    OPT_IFERRORGOTOLABEL(r, code, OPT_MEMORY_ERROR, termination);
    
    code = 0;
    
termination:

    if(cols)   free(cols);
    
    return code;
}




int OPT_Worhp::setuJacIndices()
{
    const unsigned int n = uJac.getNumberOfColumns();
    const unsigned int m = uJac.getNumberOfRows();
    
    const double *lc = prob.lc;
    const double *uc = prob.uc;
    
    
    int code;
    unsigned int *auxIndex = nullptr;
    
    //worhp needs we set indexes in column order (I HATE IT!!!)
    
    
    #if OPT_DEBUG_MODE
        assert(  m == (unsigned) prob.m );
        uJac.setAllSparseMatrix(UINT_MAX); //if we miss to set some position, we will have a memory error when try use the respective index...
    #endif
    
    
    OPT_malloc(auxIndex, n+1);
    OPT_IFMEMERRORGOTOLABEL( !auxIndex, code, termination );
    
    //we use auxIndex[i] to store the next position to column i
    
    //uJac.countRowsEachColumn(&auxIndex[1]);
    
    //counting the number of rows for each columns
    OPT_setAllArray<unsigned int>(n+1, auxIndex, 0.0);
    
    {
        //auxIndex[i] should have the number of elements in col i-1
        auto pauxIndex = &auxIndex[1];
        
        for( auto rowit = uJac.beginRowIndex(); *rowit < m ; ++rowit )
        {
            auto row = *rowit;
            
            if( lc[row] > -OPT_INFINITY || uc[row] < OPT_INFINITY ) //we should not consider the free constraint
            {
                const unsigned int rnel = uJac.getNumberOfElementsAtRow(row);
                
                const auto cols = uJac.getRowColsPointer( row );
                
                for(unsigned int j = 0; j < rnel; j++)
                    pauxIndex[ cols[j] ]++;
            }
        }
    }
    
    
    
    for(unsigned int i = 2; i <= n; i++)
        auxIndex[i] += auxIndex[i-1]; //we have to accumulate the sum of itens in columns
    
    
    /*cout << "auxIndex: " << endl;
    for(int i = 0; i <= n; i++)
        cout << "auxIndex[" << i << "]: " << auxIndex[i] << endl;
    
    cout << "prob.J: " << endl;
    prob.J.printSparseMatrix();
    
    cout << "prob.A: " << endl;
    prob.A.printSparseMatrix(); */
    
    
    for( auto rowit = uJac.beginRowIndex(); *rowit < m ; ++rowit )
    {
        auto i = *rowit;
        
        if( lc[i] > -OPT_INFINITY || uc[i] < OPT_INFINITY ) //we should not consider the free constraint
        {
            const unsigned int nel = uJac.getNumberOfElementsAtRow(i);
            const auto *rcols = uJac[i];
            unsigned int *rvalues = uJac(i);
            
            for( unsigned int j = 0; j < nel; j++ )
            {
                const unsigned int col =  rcols[j];//row[j].getColumn();
                
                //auxIndex[col] is initialized having the number of elements until the next guy in column col considering column ordering.
                
                rvalues[j] = auxIndex[col]; //row[j].setValue(  auxIndex[col] ); //do not forget: that sparse matrix is about unsigned integer, not double
                
                auxIndex[col]++;
            }
        }
    }
    
    
    code = 0;
termination:
    
    if(auxIndex) free(auxIndex);
    
    return code;
}



int OPT_Worhp::setuHessIndices()
{
    const unsigned int n = uHess.getNumberOfColumns();
    const unsigned int nnondiag = uHess.getNumberOfElements() - n; // all diagonal positions should be here
    
    int code;
    unsigned int *auxIndex = nullptr;
    
    #if OPT_DEBUG_MODE
        assert( uHess.getNumberOfElements() >= n );
    #endif
    
    
    //worhp needs we set indexes in column order (I HATE IT!!!). For hessians, nondiagonal elements should came first. Diagonal elements should be in the last positions
    
    #if OPT_DEBUG_MODE
        uHess.setAllSparseMatrix(UINT_MAX); //if we miss to set some position, we will have a memory error when try use the respective index...
    #endif
    
        
    OPT_malloc(auxIndex, n+1);
    OPT_IFMEMERRORGOTOLABEL( !auxIndex, code, termination );    
    
    
    //we use auxIndex[i] to store the next position to column i
    auxIndex[0] = 0;
    uHess.countRowsEachColumn( &auxIndex[1] );
    
    
    //unfortunatelly, diagonal positions should be the last indices. So, we have to uncount diagonal positions
    for(unsigned int i = 1; i <= n; i++)
        auxIndex[i] += auxIndex[i-1] - 1; //we have to accumulate the sum of itens in columns //note we start from the 1 to count the indices
    
    
    //for(unsigned int i = 1; i <= n; i++)
        //auxIndex[i] -= i; //note we start from the 1 to count the indices
    
    
    for(unsigned int i = 0; i < n; i++)
    {
        //OPT_UIntSparseRow &row = uHess[i];
        //const unsigned int nel = row.getNumberOfElements();
        
        const unsigned int nel = uHess.getNumberOfElementsAtRow(i);
        
        #if OPT_DEBUG_MODE 
            assert( nel > 0 ); //at least the adiagonal should be int he row...
        #endif
        
        //the elements in uHess are sorted, here we assume the adiagonal element is the last in the row...
        
        const int* const rcols = uHess[i];
        unsigned int* const rvalues = uHess(i);
        
        for( unsigned int j = 0; j < nel-1; j++)
        {
            const unsigned int col = rcols[j]; //row[j].getColumn();
            
            rvalues[j] = auxIndex[col]; // row[j].setValue( auxIndex[col] ); //do not forget: that sparse matrix is about unsigned integer, not double
            
            auxIndex[col]++;
        }
        
        #if OPT_DEBUG_MODE
            
            assert( (unsigned int) rcols[nel-1] == i ); //assert( row[nel-1].getColumn() == i ); //the last element in the row must be the diagonal position...
        #endif
        
        rvalues[nel-1] = nnondiag + i; //row[nel-1].setValue( nnondiag + i );
        
    }
    
    
    code = 0;
termination:
    
    if(auxIndex) free(auxIndex);
    
    return code;
}
































