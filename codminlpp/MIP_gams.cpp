/*
* Interface to GAMS (General Algebraic Modeling System) 
* 
* GAMS has no reference guide about how to use its library. Our references are the examples:
* 
* https://forum.gamsworld.org/viewtopic.php?t=8206
* 
* Examples from GAMS -CoinOR:
* 
* https://projects.coin-or.org/GAMSlinks/browser/trunk/GAMSlinks/src/SolverInterfaces
* 
* 
* Reference to GAMS users learn abou setting solver option: 
* 
* 	https://www.gams.com/latest/docs/UG_SolverUsage.html
* 
* 
* Author: Wendel Melo
* 
* 30-January-2018
*
*/ 


#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>

#include <new>

#include "MIP_gams.hpp"
#include "MIP_tools.hpp"




using namespace minlpproblem;



#if MIP_HAVE_GAMS
static inline int MIP_initializeGamsStructures(const char* stub, struct gmoRec* &gmo, struct gevRec* &gev)
{
    int r, retCode;
    char buffer[MIP_GAMS_BUFFER_SIZE] = "\0";
    r = gmoCreate(&gmo, buffer, MIP_GAMS_BUFFER_SIZE); // initialize GMO  libraries
    if(!gmo)
    {
        MIP_PRINTERRORMSGP("Error to greate GAMS object check if GAMS directory is in your system librray path. Original error msg: ", buffer);
        retCode = MIP_UNDEFINED_ERROR;
        goto termination;
    }
    
    
    r = gevCreate(&gev, buffer, MIP_GAMS_BUFFER_SIZE); // initialize GEV  libraries
    MIP_IFERRORGOTOLABELANDPRINTBUFFER(!r || !gev, retCode, MIP_UNDEFINED_ERROR, buffer, termination); //this gams function returns 1 on success. Aff...
    
    // load control file
    r = gevInitEnvironmentLegacy(gev, stub);
    if(r)
        MIP_PRINTERRORMSGP("Could not load control file ", stub);
    MIP_IFERRORGOTOLABEL(r, retCode, MIP_UNDEFINED_ERROR, termination);
    
    
    r = gmoRegisterEnvironment(gmo, gev, buffer);
    MIP_IFERRORGOTOLABELANDPRINTBUFFER(r, retCode, MIP_UNDEFINED_ERROR, buffer, termination);
    
    r = gmoLoadDataLegacy(gmo, buffer);
    MIP_IFERRORGOTOLABELANDPRINTBUFFER(r, retCode, MIP_UNDEFINED_ERROR, buffer, termination);
    
    /* reformulate objective variable out of model, if possible */
    gmoObjStyleSet(gmo, gmoObjType_Fun);
    
    
    {
        int do2dir = 0;
        int dohess = 1;
        r = gmoHessLoad(gmo, 0, &do2dir, &dohess); //I have no idea about do2dir and dohess parameters
        if( !dohess )
        {
            MIP_IFERRORGOTOLABEL(r, retCode, MIP_UNDEFINED_ERROR, termination);
        }
    }
    
    
    retCode = 0;
    
termination:
    
    return retCode;
}


static inline void MIP_desallocateGamsStructures( struct gmoRec* &gmo, struct gevRec* &gev)
{
    if(gmo)
    {
        gmoFree(&gmo);
        gmo = NULL;
    }
    if(gev)
    {
        gevFree(&gev);
        gev = NULL;
    }
}
#endif






MIP_NLEvalGams::MIP_NLEvalGams(const char *stub, MIP_MINLPProb *prob):MIP_NonLinearEval()
{
    this->stub = stub;
    this->prob = prob;
    #if MIP_HAVE_GAMS
        gmos_ = NULL;
        gevs_ = NULL;
    #endif
    auxValues_ = NULL;
    nthreads = 0;
}


MIP_NLEvalGams::~MIP_NLEvalGams()
{
    desallocate();
}


void MIP_NLEvalGams::desallocate()
{
    #if MIP_HAVE_GAMS
        for(decltype(nthreads) k = 0; k < nthreads; k++)
            MIP_desallocateGamsStructures(gmos_[k], gevs_[k]);
        
        MIP_secFree(gmos_);
        MIP_secFree(gevs_);
        
        if(auxValues_)
        {
            for(decltype(nthreads) k = 0; k < nthreads; k++)
            {
                if(auxValues_[k])
                    free(auxValues_[k]);
            }
            free(auxValues_);
            auxValues_ = NULL;
        }
    #endif
    nthreads = 0;
}


int MIP_NLEvalGams::initialize(const int nthreads, const int n, const int m, const int nzJac, const int nzLagHess)
#if MIP_HAVE_GAMS
{
    int r;
    this->nthreads = nthreads;
    
    
    nEvalErros = 0;
    
    //setting gams objects to perform evaluations
    MIP_calloc(gmos_, nthreads);
    MIP_calloc(gevs_, nthreads);
    MIP_calloc(auxValues_, nthreads);
    
    MIP_IFMEMERRORRETURN(!gmos_ || !gevs_ || !auxValues_);
    
    for(int k = 0; k < nthreads; k++)
    {
        r = MIP_initializeGamsStructures(stub, gmos_[k], gevs_[k]);
        MIP_IFERRORRETURN(r, r);
        
        MIP_malloc(auxValues_[k], MIP_max(n,m) );
        MIP_IFMEMERRORRETURN(!auxValues_[k]);
    }
    
    
    return 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif


int MIP_NLEvalGams::eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value)
#if MIP_HAVE_GAMS
{
    struct gmoRec* &gmo = gmos_[threadnumber];
    int r, nerror;
    
    if( newx )
        gmoEvalNewPoint(gmo, x);
    
    
    #if MIP_DEBUG_MODE
        assert( prob->hasObjNLTerm() );
    #endif
    
    r = gmoEvalFuncObj(gmo, x, &value, &nerror);
    
    if( r != 0 )
    {
        #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
            nEvalErros++;
            if( nEvalErros <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
            {
                std::cerr << MIP_PREPRINT "Error " << r << " at GAMS objective function evaluation" << MIP_GETFILELINE << "\n";
            }
        #endif
        return r;
    }
    
    if(nerror != 0)
    {
        #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
            nEvalErros++;
            if( nEvalErros <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                std::cerr << MIP_PREPRINT "Error " << nerror << " at GAMS objective function evaluation" << MIP_GETFILELINE << "\n";
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


int MIP_NLEvalGams::eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *values)
#if MIP_HAVE_GAMS
{
    struct gmoRec* &gmo = gmos_[threadnumber];
    int r, nerror;
    double val, gx;
    
    if( newx )
        gmoEvalNewPoint(gmo, x);
    
    
    #if MIP_DEBUG_MODE
        assert( prob->hasObjNLTerm() );
    #endif
    
    
    MIP_setAllArray<double>(n, values, 0.0); //gams coin-or example perform that, so I am just copying...
    
    r = gmoEvalGradObj(gmo, x, &val, values, &gx, &nerror); //I have no idea what gx and val are (maybe val is the objective value...)
    if( r != 0 )
    {
        #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
            nEvalErros++;
            if( nEvalErros <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                std::cerr << MIP_PREPRINT "Error " << r << " at GAMS objective function gradient evaluation" << MIP_GETFILELINE << "\n";
        #endif
        return r;
    }
    
    if(nerror != 0)
    {
        #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
            nEvalErros++;
            if( nEvalErros <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                std::cerr << MIP_PREPRINT "Error " << nerror << " at GAMS objective function gradient evaluation" << MIP_GETFILELINE << "\n";
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


int MIP_NLEvalGams::eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *ctrEval, const double *x, double *values)
#if MIP_HAVE_GAMS
{
    struct gmoRec* &gmo = gmos_[threadnumber];
    int r, nerror;
    bool *nlConstr = prob->nlConstr;
    
    
    if( newx )
        gmoEvalNewPoint(gmo, x);
    
    
    for(int i = 0; i < m; i++)
    {
        if(nlConstr[i] && (!ctrEval || ctrEval[i]) )
        {
            r = gmoEvalFunc(gmo, i, x, &values[i], &nerror);
            
            if( r != 0 )
            {
                #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
                    nEvalErros++;
                    if( nEvalErros <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                        std::cerr << MIP_PREPRINT "Error " << r << " at GAMS constraint " << i << " evaluation" << MIP_GETFILELINE << "\n";
                #endif
                return r;
            }
            
            if(nerror != 0)
            {
                #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
                    nEvalErros++;
                    if( nEvalErros <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                        std::cerr << MIP_PREPRINT "Error " << nerror << " at GAMS constraint " << i << " evaluation" << MIP_GETFILELINE << "\n";
                #endif
                return nerror;
            }
        }
    }
    
    return 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif


int MIP_NLEvalGams::eval_grad_nl_constrs_part( const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *ctrEval, const double *x, MIP_SparseMatrix& J)
#if MIP_HAVE_GAMS
{
    struct gmoRec* &gmo = gmos_[threadnumber];
    double* values = auxValues_[threadnumber];
    bool *nlConstr = prob->nlConstr;
    
    int r, nerror;
    double val, gx;
    
    if( newx )
        gmoEvalNewPoint(gmo, x);
    
    
    for(int i = 0; i < m; i++)
    {
        if(nlConstr[i] && (!ctrEval || ctrEval[i]) )
        {
            r = gmoEvalGrad(gmo, i, x, &val, values, &gx, &nerror); //gmoEvalGrad evaluates the fullgradient
            
            if( r != 0 )
            {
                #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
                    nEvalErros++;
                    if( nEvalErros <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                        std::cerr << MIP_PREPRINT "Error " << r << " at GAMS gradient of constraint " << i << " evaluation" << MIP_GETFILELINE << "\n";
                #endif
                return r;
            }
            
            if(nerror != 0)
            {
                #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
                    nEvalErros++;
                    if( nEvalErros <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                        std::cerr << MIP_PREPRINT "Error " << nerror << " at GAMS gradient of constraint " << i << " evaluation" << MIP_GETFILELINE << "\n";
                #endif
                return nerror;
            }
            
            
            J.setFullRowValues(i, values);
        }
    }
    
    
    return 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif


int MIP_NLEvalGams::eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double obj_factor, const double *lambda, MIP_SparseMatrix& hessian)
#if MIP_HAVE_GAMS
{
    struct gmoRec* &gmo = gmos_[threadnumber];
    double *myLambda = auxValues_[threadnumber];
    
    bool *nlConstr = prob->nlConstr;
    
    int r, nerror;
    double myObjF = prob->hasNlObj ? obj_factor : 0.0;
    
    if( newx )
        gmoEvalNewPoint(gmo, x);
    
    if( lambda )
    {
        for(int k = 0; k < m; k++)
            myLambda[k] = nlConstr[k] ? lambda[k] : 0.0;
    }
    else //if we do not have lambda, we must assume all lambda are zeros
    {
        MIP_setAllArray<double>( m, myLambda, 0.0 );
    }
    
    // coment from GamsNLP.cpp: // for GAMS, lambda would need to be multiplied by -1, we do this via the constraint weight
    r = gmoHessLagValue(gmo, x, myLambda, hessian.getRowValuesPointer(0), myObjF, -1.0, &nerror );
    
    if( r != 0 )
    {
        #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
            nEvalErros++;
            if( nEvalErros <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                std::cerr << MIP_PREPRINT "Error " << r << " at GAMS hessian of lagrangian evaluation" << MIP_GETFILELINE << "\n";
        #endif
        return r;
    }
    
    if(nerror != 0)
    {
        #if MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR
            nEvalErros++;
            if( nEvalErros <= MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS )
                std::cerr << MIP_PREPRINT "Error " << nerror << " at GAMS hessian of lagrangian evaluation" << MIP_GETFILELINE << "\n";
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


void MIP_NLEvalGams::finalize(const int nthreads, const int n, const int mnl, const int nzNLJac, const int nzNLLagHess)
{
    desallocate();
}








MIP_GamsModelReader::MIP_GamsModelReader()
{
    #if MIP_HAVE_GAMS
        gmo = NULL;
        gev = NULL;
        //opt = NULL;
    #endif
}


MIP_GamsModelReader::~MIP_GamsModelReader()
{
    desallocate();
}


void MIP_GamsModelReader::desallocate()
{
    #if MIP_HAVE_GAMS
        MIP_desallocateGamsStructures(gmo, gev);
    #endif
}



int MIP_GamsModelReader::getDblOption(const char *optionName, double &value)
#if MIP_HAVE_GAMS
{
    value = gevGetDblOpt(gev,optionName);
    return 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif



//to gen an integer option (parameter)
int MIP_GamsModelReader::getIntOption(const char *optionName, int &value)
#if MIP_HAVE_GAMS
{
    value = gevGetIntOpt(gev, optionName);
    return 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif



int MIP_GamsModelReader::getStrOption(const char *optionName, std::string &value)
#if MIP_HAVE_GAMS
{
    char myvalue[10000];
    
    gevGetStrOpt(gev, optionName, myvalue);
    value = myvalue;
    
    return 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif



int MIP_GamsModelReader::getInitialSolution(double *solution)
#if MIP_HAVE_GAMS
{
    
    //setting initial values of variables
    if( gmo )
    {
        double *xInit = solution;
        
        int r = gmoGetVarL(gmo, xInit);
        MIP_IFERRORRETURN(r, MIP_UNDEFINED_ERROR);
        
        return 0;
    }
    
    return MIP_UNDEFINED_ERROR;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif


bool MIP_GamsModelReader::isMaximizationProblem()
#if MIP_HAVE_GAMS
{
    return gmoSense(gmo) == gmoObj_Max ; //maximization problem
}
#else
{
    return false;
}
#endif


int MIP_GamsModelReader::getOptionFileName(std::string &fileName)
#if MIP_HAVE_GAMS
{
    //trying reading options
    if( gmoOptFile(gmo) > 0 )
    {
        char optFileName[10000];
        gmoNameOptFile(gmo, optFileName);
        
        fileName = optFileName;
        
        return 0;
    }
    else
    {
        return MIP_NOT_APPLICABLE;
    }
    
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif




int MIP_GamsModelReader::readProblem(const char* stub, MIP_MINLPProb& prob, MIP_NonLinearEval* &eval)
#if MIP_HAVE_GAMS
{
    int n, m, nzJac, nzHess, nzObjQ, maxNz;
    int r, retCode;
    
    //char buffer[MIP_GAMS_BUFFER_SIZE] = "\0";
    int *rows = NULL, *cols = NULL;
    double *values = NULL;
    
    
    r = MIP_initializeGamsStructures(stub, gmo, gev);
    MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    
    gmoUseQSet(gmo, 1); //to active detection of quadratic objective and constraints
    
    
    n = gmoN(gmo); //number of variables
    m = gmoM(gmo); //number of constraints
    
    //std::cout << MIP_PREPRINT "n: " << n << " m: " << m << "\n";
    
    
    r = prob.addVariables(n);
    MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    r = prob.addConstraints(m);
    MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    nzObjQ = gmoObjQNZ(gmo);
    nzJac = gmoNZ(gmo);
    nzHess = gmoHessLagNz(gmo);
    
    //std::cout << MIP_PREPRINT "nzObjQ: " << nzObjQ << " nzJac: " << nzJac << " nzHess: " << nzHess << "\n";
    
    maxNz = MIP_max(n,m);
    
    for(int i = 0; i < m; i++)
    {
        const int rnzQ = gmoGetRowQNZOne(gmo,i);
        if( rnzQ > maxNz )
            maxNz = rnzQ;
    }
    if( nzObjQ > maxNz )
        maxNz = nzObjQ;
    if( nzJac > maxNz )
        maxNz = nzJac;
    if( nzHess > maxNz )
        maxNz = nzHess;
    
    
    MIP_malloc(values, maxNz );
    MIP_malloc(rows, maxNz);
    MIP_malloc(cols, maxNz);
    
    MIP_IFMEMERRORGOTOLABEL( !values, retCode, termination );
    
    {//setting variables bounds
        double *lx = values, *ux = values;
        
        gmoGetVarLower(gmo, lx);
        r = prob.setVariableLowerBounds(n, lx);
        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        gmoGetVarUpper(gmo, ux);
        r = prob.setVariableUpperBounds(n, ux);
        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
    }
    
    
    //setting variables types:  gmovar_I, gmovar_B
    for(int i = 0; i < n; i++)
    {
        auto vt = gmoGetVarTypeOne(gmo, i);
        MIP_VARTYPE myvt;
        
        if( vt == gmovar_I || vt == gmovar_B )
            myvt = MIP_VT_INTEGER;
        else
            myvt = MIP_VT_CONTINUOUS;
        
        r = prob.setVariableType(i, myvt);
        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
    }
    
    
    
    {//setting constraints bounds
        double *rhs = values;
        double lci, uci;
        
        gmoGetRhs(gmo, rhs);
        
        for(int i = 0; i < m; i++)
        {
            switch( gmoGetEquTypeOne(gmo, i) )
            {
                case gmoequ_E:	//equality
                    lci = uci = rhs[i];
                    break;
                case gmoequ_G: // greater-equal
                    lci = rhs[i];
                    uci = MIP_INFINITY;
                    break;
                case gmoequ_L:	// lower-equal
                    lci = -MIP_INFINITY;
                    uci = rhs[i];
                    break;
                case gmoequ_N:
                    lci = -MIP_INFINITY;
                    uci = MIP_INFINITY;
                    break;
                default:
                    MIP_PRINTERRORMSGP("Unsupported equation type ", i);
            }
            
            r = prob.setConstraintLowerBound(i, lci);
            MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            r = prob.setConstraintUpperBound(i, uci);
            MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
        }
    }
    
    
    //setting objective
    {
        if( gmoGetObjOrder(gmo) == gmoorder_L || gmoGetObjOrder(gmo) == gmoorder_Q ) //check if objective is linear os quadratic
        {
            double *c = values;
            
            r = gmoGetObjVector( gmo, c, NULL );
            MIP_IFERRORGOTOLABEL(r, retCode, MIP_UNDEFINED_ERROR, termination);
            
            r = prob.setObjLinearCoefficients(c);
            MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            prob.setObjConstant( gmoObjConst(gmo) );
            
            if( gmoGetObjOrder(gmo) == gmoorder_Q ) //check if objective is quadratic
            {
                r = gmoGetObjQ(gmo, cols, rows, values);
                MIP_IFERRORGOTOLABEL(r, retCode, MIP_UNDEFINED_ERROR, termination); //we invert cols and rows to get the lower triangle, 
                
                for(int i = 0; i < nzObjQ; i++)
                    std::cout << "obj Q - 1 - row: " << rows[i] << " col: " << cols[i] << " values: " << values[i] << "\n";
                
                //gams already assumes qudratic Q is multiplied by 0.5
                r = prob.setObjQuadCoefsMatrix( nzObjQ, rows, cols, values );
                MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
            }
        }
        else
        {
            prob.setObjNonLinearTermFlag(true);
        }
        
        if( gmoSense(gmo) == gmoObj_Max ) //maximization problem
            prob.setObjFactor(-1.0);
        
    }
    
    
    
    //gmoGetVarPrior  //to get priorities, i think...
    
    
    /***********
    * setting constraints linear terms and jacobian
    * We cannot use gmoGetMatrixRow if we have used gmoUseQSet(gmo, 1). So, first we try detect quadratic constraints, and, after, we will use gmoUseQSet(gmo, 0) to detect the general constraints.
    ***********/ 
    {
        int *nlFlag = rows;
        int rlnz, rnlnz; //linear terms in a row, nonlinear termos in a row
        
        
        for(int i = 0; i < m; i++)
        {
            const auto constrType = gmoGetEquOrderOne(gmo, i);
            
            r = gmoGetRowSparse(gmo, i, cols, values, nlFlag, &rlnz, &rnlnz); //cols and values will have the row jacobian structure. I think when in mode gmoUseQSet(gmo, 1), rnlnz always will be zero
            MIP_IFERRORGOTOLABEL(r, retCode, MIP_UNDEFINED_ERROR, termination);
            
            std::cout << "constr: " << i << " rlnz: " << rlnz << " rnlnz: " << rnlnz << "\n";
            
            if( constrType == gmoorder_L || constrType == gmoorder_Q ) //linear constraint
            {
                //looking the coenne example, I have the ideia linear terms appears first in cols and values, so, lets check
                #if MIP_DEBUG_MODE
                    for(int k = 0; k < rlnz; k++)
                        assert( nlFlag[k] == 0 );
                #endif
                
                r = prob.setConstraintLinearPart(i, rlnz, cols, values );
                MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                
                if( constrType == gmoorder_Q )
                {
                    int rnzQ = gmoGetRowQNZOne(gmo,i);
                    
                    r = gmoGetRowQ(gmo, i, cols, rows, values); //we invert rows and cols to get the lower triangle
                    MIP_IFERRORGOTOLABEL(r, retCode, MIP_UNDEFINED_ERROR, termination);
                    
                    for(int k = 0; k < rnzQ; k++)
                        std::cout << "\tconstr " << i << " - Q - " << k << " - row: " << rows[k]  << " col: " << cols[k] << " value: " << values[k] << "\n";
                    
                    r = prob.setConstraintQuadCoefsMatrix(i, rnzQ, rows, cols, values); 
                    MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                }
                
            }
            else if( constrType == gmoorder_NL ) //nonlinear constraint
            {
                r = prob.setConstraintNonLinearTermFlag(i, true);
                MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                
                //r = prob.setJacobianRowStructure(i, nzrow, cols);
                //MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
            }
            else
            {
                MIP_PRINTERRORMSGP("Invalid value for constraint type: ", constrType);
                retCode = MIP_UNDEFINED_ERROR;
                goto termination;
            }
        }
    }
    
    
    
    if( prob.hasNLConstraints() )
    {
        gmoUseQSet(gmo, 0); //so, we desative the detection of quadratic to get the general nonlinear constraints
        
        
        int *rowStart = rows;
        
        r = gmoGetMatrixRow(gmo, rowStart, cols, values, NULL); //function gmoGetMatrixRow gets the jacobian
        MIP_IFERRORGOTOLABEL(r, retCode, MIP_UNDEFINED_ERROR, termination);
        
        rowStart[m] = nzJac; //just to turn thigs easier below
        
        
        for(int i = 0; i < m; i++)
        {
            bool nlConstr;
            
            int r = prob.getConstraintsNonLinearTermFlag(i, nlConstr);
            MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            if( nlConstr ) //note, now quadratic constraints are considered nonlinear by gams, so do not use gmoGetEquOrderOne(gmo, i) here
            {
                #if MIP_DEBUG_MODE
                    assert( gmoGetEquOrderOne(gmo, i) == gmoorder_NL );
                #endif
                const int nzrow = rowStart[i+1] - rowStart[i];
                
                
                for(int k = 0; k < nzrow; k++)
                    std::cout << "\tconstr " << i << " - Jac - " << k  << " col: " << cols[rowStart[i] + k] << "\n";
                
                r = prob.setJacobianRowStructure(i, nzrow, &cols[rowStart[i]]);
                MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
            }
        }
    }
    
    
    //set lagrangian hessian structure
    {
        r = gmoHessLagStruct(gmo, rows, cols);
        MIP_IFERRORGOTOLABEL(r, retCode, MIP_UNDEFINED_ERROR, termination);
        
        //for(int i = 0; i < nzHess; i++)
            //std::cout << "rows["<<i<<"]: " << rows[i] << " cols["<<i<<"]: " << cols[i] <<  std::endl;
        
        r = prob.setLagrangianHessianStructure(nzHess, rows, cols);
        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
    }
    
    
    
    if( prob.hasObjNLTerm() || prob.hasNLConstraints() )
    {
        eval = new (std::nothrow) MIP_NLEvalGams(stub, &prob);
        
        MIP_IFMEMERRORGOTOLABEL(!eval, retCode, termination);
        prob.setNonLinearEvaluationObject(eval);
    }
    
    
    #if 0
    //testando avaliação de função. Testado a função gmoEvalFunc avalia toda a restricao, incluindo termos lineares
    {
        int nerror;
        double g1, g2, x[n];
        
        MIP_setAllArray<double>(n, x, 10.0);
        
        r = gmoEvalFunc(gmo, 2, x, &g1, &nerror);
        MIP_IFERRORGOTOLABEL(r, retCode, MIP_UNDEFINED_ERROR, termination);
        
        r = gmoEvalFuncNL(gmo, 2, x, &g2, &nerror);
        
        std::cout << "Avaliei restricao 2. gmoEvalFunc: " << g1 << " gmoEvalFuncNL: " << g2 << "\n";
    }
    #endif
    
    
    prob.print();
    
    
    #if 0
    //trying reading options
    if( gmoOptFile(gmo) > 0 )
    {
        char optFileName[10000];
        gmoNameOptFile(gmo, optFileName);
        
        std::cout << "optFileName: " << optFileName << "\n";
    }
    #endif
    
    /*std::cout << "iterlim: " << gevGetIntOpt(gev, "iterlim") << "\n";
    std::cout << "tolproj: " << gevGetDblOpt(gev, "tolproj") << "\n";
    std::cout << "nlp: " << gevGetDblOpt(gev, "NLP") << "\n"; */
    
    
    
    retCode = 0;
    
termination:
    
    if(values)	free(values);
    if(rows)	free(rows);
    if(cols)	free(cols);
    
    return retCode;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif


int MIP_GamsModelReader::setSolution( const int modelStatus, const int solverStatus, const double *sol)
#if MIP_HAVE_GAMS
{
    int r;
    
    gmoModelStatSet(gmo, modelStatus);
    
    gmoSolveStatSet(gmo, solverStatus);
    
    if( sol )
    {
        r = gmoSetSolutionPrimal(gmo, sol); // to set dual variables or consraint values, use some other function like gmoSetSolution, gmoSetSolution2, gmoSetSolution8
        MIP_IFERRORRETURN(r, MIP_UNDEFINED_ERROR);
    }
    
    
    //this function should make library send solution to Gams
    r = gmoUnloadSolutionLegacy(gmo);
    MIP_IFERRORRETURN(r, MIP_UNDEFINED_ERROR);
    
    return 0;
}
#else
{
    return MIP_LIBRARY_NOT_AVAILABLE;
}
#endif



