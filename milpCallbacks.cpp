

#include <cstdlib>
#include <cassert>
#include <new>
#include "MRQ_solvers.hpp"
#include "MRQ_milpCallbacks.hpp"

using namespace muriqui;





#if OPT_HAVE_GLPK
void muriqui::MRQ_bboaGlpkCallback(glp_tree *T, void *info)
{
    MRQ_getchar();
    
    if(glp_ios_reason(T) != GLP_IROWGEN)
        return; //we just run callback if we are in a MIP incombent solution. Unfortunatelly, glpk only accepts one callback to call to everything. Node, glpk manual says we should not modify the problem if the status is GLP_IBINGO
    
    const unsigned int nparams = -1; //size of array params
    void* params[] = {T};
    unsigned int thnumber = 0;
    MRQ_MILPCallbackData &data = ((MRQ_MILPCallbackData*) info)[thnumber];
    
    
    const int n = data.prob->n;
    const int nI = data.nI;
    const int *intVars = data.intVars;
    double *sol = data.auxVars;
    MRQ_LinearApproxAlgorithm *alg = data.alg;
    
    data.callbackSolver->initializeSolverData(nparams, params);
    data.callbackSolver->getNodeSolution(n, sol);
    
    MRQ_getchar();
    
    //we just run our callback if we have a integer solution
    
    if(MRQ_isIntegerSol(nI, intVars, sol, alg->in_integer_tol) == false)
        return;
    
    MRQ_getchar();
    
    int r = alg->solverCallbackLazyConstraints(data);
    
    if(r != 0)
    {
        if(alg->in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        
        glp_ios_terminate(T);
    }
    
}
#endif



#if OPT_HAVE_GUROBI
int __stdcall muriqui::MRQ_bboaGurobiCallback( GRBmodel *model, void *cbdata, int where, void *usrdata)
{
    if(where != GRB_CB_MIPSOL)
        return 0; //we just run callback if we are in a MIP incombent solution. Unfortunatelly, gurobi only accepts one callback to call to everything (I hate Gurobi)
    
    MRQ_MILPCallbackData *data = (MRQ_MILPCallbackData*) usrdata;
    
    const unsigned int nparams = 3;
    void* params[] = {model, cbdata, &where};
    
    unsigned int thnumber = 0; //gurobi only call this callback from thread 0. So, we allocate memory for all threads but only use for threa 0. We make this to keep compatibility with others solvers...
    
    data[thnumber].callbackSolver->initializeSolverData(nparams, params);
    
    return data[thnumber].alg->solverCallbackLazyConstraints( data[thnumber]);
}
#endif


#if OPT_HAVE_CPLEX

int CPXPUBLIC muriqui::MRQ_labbLazyConstraintCplexCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
    MRQ_MILPCallbackData *data = (MRQ_MILPCallbackData*) cbhandle;
    
    unsigned int thnumber; // I am not sure about the correct type to thread number. Unfortunatelly, CPXgetcallbackinfo receives a void pointer.
    
    const unsigned int nparams = 4; //size of array params
    void* params[] = {&env, cbdata, &wherefrom, useraction_p};
    
    
    int r = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &thnumber);
    if(r != 0)
    {
        if(data->alg->in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        
        return MRQ_MILP_SOLVER_ERROR;
    }
    
    //std::cout << "thnumber: " << thnumber << "\n";
    
    
    data[thnumber].callbackSolver->initializeSolverData(nparams, params);
    
    r = data[thnumber].alg->solverCallbackLazyConstraints( data[thnumber]);
    
    if(r != 0 && r != MRQ_LAZY_MILP_BB_SOLUTION_LOWER_THAN_ZL)
    {
        if(data->alg->in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
    }
    
    return r;
}


int CPXPUBLIC muriqui::MRQ_labbCplexBeforeSolveCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
    MRQ_MILPCallbackData *data = (MRQ_MILPCallbackData*) cbhandle;
    
    unsigned int thnumber; // I am not sure about the correct type to thread number. Unfortunatelly, CPXgetcallbackinfo receives a void pointer.
    //double myzl;
    MRQ_LinearApproxAlgorithm *alg = data[0].alg;
    
    const unsigned int nparams = 4; //size of array params
    void* params[] = {&env, cbdata, &wherefrom, useraction_p};
    
    
    int r = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &thnumber);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    //std::cout << "MRQ_labbCplexBeforeSolveCallback - thnumber: " << thnumber << "\n";
    
    /*r = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_REMAINING, &myzl);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    //std::cout << "myzl: " << myzl << "\n";
    
    if( myzl > alg->zl ) //we should protect zl updating with a semaphore. Howerver, that would be so expensive and would not bring us many benefits. So, we perform a non-protected updating. It is enough for ur porpouses.
        alg->zl = myzl; */
    
    
    data[thnumber].callbackSolver->initializeSolverData(nparams, params);
    
    
    r = alg->solverCallbackBeforeSolve( data[thnumber] );
    
    if(r != 0 )
    {
        if(data->alg->in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
    }
    
    
    return r;
}


int CPXPUBLIC muriqui::MRQ_labbCplexBranchingCallback(CPXCENVptr env, void *cbdata, int          wherefrom, void *cbhandle, int brtype, int sos, int nodecnt, int bdcnt, const int *nodebeg, const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p)
{
    
    /* nodecnt is the number of new nodes being generates. If it equals to zero, the current node being "branched will be already fathomed"*/
    if( nodecnt == 0)
        return 0;
    
    {
        MRQ_MILPCallbackData *data = (MRQ_MILPCallbackData*) cbhandle;
        
        unsigned int thnumber; // I am not sure about the correct type to thread number. Unfortunatelly, CPXgetcallbackinfo receives a void pointer.
        //double myzl;
        MRQ_LinearApproxAlgorithm *alg = data[0].alg;
        
        const unsigned int nparams = 4; //size of array params
        void* params[] = {&env, cbdata, &wherefrom, useraction_p};
        
        const unsigned int nbranchParams = 9;
        void* branchParams[] = {&brtype, &sos, &nodecnt, &bdcnt, &nodebeg, &indices, &lu, &bd, &nodeest};
        
        
        
        int r = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &thnumber);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        
        data[thnumber].callbackSolver->initializeSolverData(nparams, params);
        
        
        r = data[thnumber].callbackSolver->initializeExtraBranchingSolverParameters( nbranchParams, branchParams );
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        
        r = alg->solverCallbackBranching( data[thnumber] );
        
        if(r != 0 )
        {
            if(data->alg->in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
        }
        
        
        return r;
    }
}


void CPXPUBLIC muriqui:: MRQ_labbCplexDeleteNodeCallback( CPXCENVptr env, int wherefrom, void *cbhandle, int seqnum, void *handle)
{
    //handle is the pointer pointer to the user private data that was assigned to the node when it was created with one of the callback branching routines.
    if(handle)
    {
        MRQ_SizeIndexSol *p = (MRQ_SizeIndexSol *) handle;
        delete p;
    }
}


#endif





MRQ_MILPSolverCallbackInterface::MRQ_MILPSolverCallbackInterface(unsigned int nthreads, MRQ_Mutex *SEMAPH_addConstr)
{
    this->nthreads = nthreads;
    this->SEMAPH_addConstr = SEMAPH_addConstr;
    prunedNode = false; 
}


int MRQ_MILPSolverCallbackInterface::getNodeObjValue(double &objValue)
{
    return MRQ_NONIMPLEMENTED_ERROR;
}


int MRQ_MILPSolverCallbackInterface::setCut(const int nz, const int *cols, const double *vals, const double lb, const double ub)
{
    return MRQ_NONIMPLEMENTED_ERROR;
}


int MRQ_MILPSolverCallbackInterface::setSolution(const int n, const double *sol)
{
    return MRQ_NONIMPLEMENTED_ERROR;
}


int MRQ_MILPSolverCallbackInterface::getVarsLowerBoundsOnNode( const int n, double *lbx)
{
    return MRQ_NONIMPLEMENTED_ERROR;
}


int MRQ_MILPSolverCallbackInterface::getVarsUpperBoundsOnNode( const int n, double *ubx)
{
    return MRQ_NONIMPLEMENTED_ERROR;
}


int MRQ_MILPSolverCallbackInterface::getNumberOfIterations(long int &niters)
{
    return MRQ_NONIMPLEMENTED_ERROR;
}


void MRQ_MILPSolverCallbackInterface::printConstraint(const int nz, const int *cols, const double *vals, const double lb, const double ub)
{
    if(lb > -MIP_INFINITY && ub < MIP_INFINITY)
        printf("%f <=  ", lb);
    
    for(int i = 0; i < nz; i++)
    {
        if(vals[i] != 0.0)
        {
            printf("%+fx%d  ", vals[i], cols[i]);
        }
    }
    
    double rhs;
    if(lb > -MIP_INFINITY)
    {
        if(ub < MIP_INFINITY)
            rhs = ub;
        else
            rhs = lb;
    }
    else
        rhs = ub;
    
    printf("<= %f\n", rhs);
}


#if OPT_HAVE_GUROBI
inline char MRQ_GurobiCallbackSolver::constraintSense(const double lb, const double ub, double &rhs)
{
    char sense;
    
    if(lb > -MIP_INFINITY)
    {
        rhs = lb;
        if(ub < MIP_INFINITY)
            sense = lb == ub ? GRB_EQUAL : 'R';
        else
            sense = GRB_GREATER_EQUAL;
    }
    else
    {
        sense = GRB_LESS_EQUAL;
        rhs = ub;
    }
    
    return sense;
}
#endif



void MRQ_GurobiCallbackSolver::initializeSolverData( const unsigned int nparams, void **solverParams)
{
#if OPT_HAVE_GUROBI
    
    #if OPT_DEBUG_MODE
        assert(nparams == 3);
    #endif
    
    model = (GRBmodel *) solverParams[0];
    cbdata = solverParams[1];
    where = *((int*) solverParams[2]);
#endif
}


int MRQ_GurobiCallbackSolver::getNumberOfIterations( long int &niters)
#if OPT_HAVE_GUROBI
{
    double value;
    int what;
    
    if( where == GRB_CB_MIPSOL )
        what = GRB_CB_MIPSOL_NODCNT;
    else if( where == GRB_CB_MIPNODE)
        what = GRB_CB_MIPNODE_NODCNT;
    else 
        what = GRB_CB_MIP_NODCNT;
    
    
    int r = GRBcbget(cbdata, where, what, &value);
    if(r != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
    }
    
    niters = value;
    
    return r;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_GurobiCallbackSolver::getNodeObjValue(double &objValue)
#if OPT_HAVE_GUROBI
{
    int r = GRBcbget(cbdata, where, GRB_CB_MIPSOL_OBJ, &objValue);
    if(r != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
    }
    return r;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_GurobiCallbackSolver::getNodeSolution(const int n, double *sol)
#if OPT_HAVE_GUROBI
{
    // I believe Gurobi does not require Mutex...
    int r = GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, sol);
    if(r != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
    }
    return r;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_GurobiCallbackSolver::setLazyConstraint( const int nz, const int *cols, const double *vals, const double lb, const double ub)
#if OPT_HAVE_GUROBI
{
    double rhs;
    char sense = constraintSense(lb, ub, rhs);
    
    if(sense == 'R')
    {
        //we have a double bounded constraint. So, we add two constraints. In this case, rhs is lb.
        int r = GRBcblazy(cbdata, nz, cols, vals, GRB_LESS_EQUAL, ub);
        if(r != 0)
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
        
        sense = GRB_GREATER_EQUAL;
    }
    
    int r = GRBcblazy(cbdata, nz, cols, vals, sense, rhs);
    if(r != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
    }
    
    return r;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_GurobiCallbackSolver::setSolution( const int n, const double *sol)
#if OPT_HAVE_GUROBI
{
    double val;
    
    int r = GRBcbsolution(cbdata, sol, &val);
    if(r != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
    }
    
    return r;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif




#if OPT_HAVE_GLPK
inline int MRQ_GlpkCallbackSolver::constraintSense( const double lb, const double ub)
{
    int sense;
    
    if( lb > -MIP_INFINITY )
    {
        if(ub < MIP_INFINITY)
            sense = GLP_DB;
        else
            sense = GLP_LO;
    }
    else
    {
        if(ub < MIP_INFINITY)
            sense = GLP_UP;
        else
            sense = GLP_FR;
    }
    
    return sense;
}
#endif



int MRQ_GlpkCallbackSolver::allocateMyIndices(const unsigned int size)
{
    myindices = (int*) malloc((size+1) * sizeof(int));
    if(!myindices)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTMEMERROR;
        #endif
        
        return MRQ_MEMORY_ERROR;
    }
    
    return 0;
}


MRQ_GlpkCallbackSolver::MRQ_GlpkCallbackSolver(unsigned int nthreads, MRQ_Mutex *SEMAPH_addConstr) : MRQ_MILPSolverCallbackInterface(nthreads, SEMAPH_addConstr)
{
    myindices = NULL;
}

MRQ_GlpkCallbackSolver::~MRQ_GlpkCallbackSolver()
{
    desallocate();
}


void MRQ_GlpkCallbackSolver::desallocate()
{
    MRQ_secFree(myindices);
}


void MRQ_GlpkCallbackSolver::initializeSolverData(const unsigned int nparams, void **solverParams)
{
    #if OPT_HAVE_GLPK
        T = (glp_tree*) solverParams[0];
        prob = glp_ios_get_prob(T);
    #endif
}


int MRQ_GlpkCallbackSolver::getNodeSolution( const int n, double *sol)
#if OPT_HAVE_GLPK
{
    for(int i = 1; i <= n; i++) //glpk counts indices from 1, not 0 (aff)
        sol[i-1] = glp_get_col_prim(prob, i);
    
    if(myindices == NULL)
    {
        //we take advantage to allocate here because here we know the correct number of variables
        const int r = allocateMyIndices(n);
        if(r != 0)
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
    }
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_GlpkCallbackSolver::setLazyConstraint(const int nz, const int *cols, const double *vals, const double lb, const double ub)
#if OPT_HAVE_GLPK
{
    #if MRQ_DEBUG_MODE
        assert(myindices); //getSolution should have been called before to allocate myindices
    #endif
    
    int r = glp_add_rows(prob, 1);
    
    if(r != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
        
        return MRQ_MILP_SOLVER_ERROR;
    }
    
    
    
    const int rindex = glp_get_num_rows(prob);
    const int sense = constraintSense(lb, ub);
    
    #pragma GCC ivdep
    #pragma ivdep
    for(int i = 0; i < nz; i++)
        myindices[i+1] = cols[i]+1; //glpk counts indices from 1, not from zero. Fist position in arrays are ignored also.. aff
    
    //nrows has the index of new row because glpk counts indices from 1, not from 0
    
    glp_set_mat_row(prob, rindex, nz, myindices, vals);
    
    glp_set_row_bnds(prob, rindex, sense, lb, ub);
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif





inline char  MRQ_CplexCallbackSolver::constraintSense( const double lb, const double ub)
{
    char value;
    
    if(lb > -MIP_INFINITY)
    {
        if(ub < MIP_INFINITY)
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



void MRQ_CplexCallbackSolver::initializeSolverData( const unsigned int nparams, void **solverParams)
{
    prunedNode = false;
    
#if OPT_HAVE_CPLEX
    
    #if OPT_DEBUG_MODE
        assert(nparams == 4);
    #endif
    
    env = *((CPXCENVptr*)  solverParams[0]);
    cbdata = solverParams[1];
    wherefrom = *((int*) solverParams[2]);
    useraction_p = (int*) solverParams[3];
    
    *useraction_p = CPX_CALLBACK_DEFAULT;
    
    //std::cout << "\n wherefrom: " << wherefrom << " ";
    /*{
        double nodeobj, ub;
        
        int r = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &nodeobj);
        if(r != 0)
        {
            MRQ_PRINTERRORNUMBER(r);
            
            //retCode = MRQ_MILP_SOLVER_ERROR;
            //goto termination;
        }
        //std::cout << "nodeobj: " << nodeobj << " ";
        
        
        r = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &ub);
        if(r != 0)
        {
            MRQ_PRINTERRORNUMBER(r);
            
            //retCode = MRQ_MILP_SOLVER_ERROR;
            //goto termination;
        }
        //std::cout << "ub: " << ub << " ";
    }*/
    
    /*{
        CPXINT nodeseq;
        long int iter;
        
        CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &nodeseq);
        
        std::cout << "cplex node sequence: " << nodeseq << " iter: " << iter << "\n";
    }*/
    
#endif
}



int MRQ_CplexCallbackSolver:: initializeExtraBranchingSolverParameters(const unsigned int nParams, void **solverParams)
#if OPT_HAVE_CPLEX
{
    #if MRQ_DEBUG_MODE
        assert(nParams == 9);
        assert(wherefrom == CPX_CALLBACK_MIP_BRANCH);
    #endif
    
    brtype = *((int*) solverParams[0]); 
    sos = *((int*) solverParams[1]);
    nodecnt = *((int*)solverParams[2]);
    bdcnt = *((int*) solverParams[3]);
    nodebeg = *((const int**) solverParams[4]);
    indices = *((const int**) solverParams[5]);
    lu = *((const char**)  solverParams[6]);
    bd = *((const double**) solverParams[7]);
    nodeest = *((const double**) solverParams[8] );
    
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::getNumberOfIterations(long int &niters)
#if OPT_HAVE_CPLEX
{
    int r = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT_LONG, &niters);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::getNodeObjValue(double &objValue) 
#if OPT_HAVE_CPLEX
{
    int r = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objValue);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::getBestSolutionObjValue( double &objValue)
#if OPT_HAVE_CPLEX
{
    int r = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &objValue);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::getNodeSolution(const int n, double *sol)
#if OPT_HAVE_CPLEX
{
    //printf("env: %p cbdata: %p wherefrom: %d\n", env, cbdata, wherefrom);
    
    int r = CPXgetcallbacknodex(env, cbdata, wherefrom, sol, 0, n-1);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver:: getNodeBranchVarIndex(int &nIndices, int *indices)
#if OPT_HAVE_CPLEX
{
    const int nodeIndex = 0; //Index of current node in cplex.
    CPXINT cplexIndex;
    
    *indices = -1;
    nIndices = 0;
    
    int r = CPXgetcallbacknodeinfo(env,cbdata, wherefrom, nodeIndex, CPX_CALLBACK_INFO_NODE_VAR,                                 &cplexIndex);
    
    //printf("r: %d  cplexIndex: %d  indices[0]: %d\n", r, cplexIndex, indices[0]);
    
    if(r == -1 || cplexIndex == -1)
    { //no variable branching. In this case, cplex did another kind of branchig like SOS, constraint, etc...
        return MRQ_UNDEFINED_ERROR;
    }
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    *indices = cplexIndex;
    nIndices = 1;
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::setLazyConstraint(const int nz, const int *cols, const double *vals, const double lb, const double ub)
#if OPT_HAVE_CPLEX
{
    const int purgeable = CPX_USECUT_PURGE;//CPX_USECUT_FORCE;//CPX_USECUT_PURGE;
    
    char sense = constraintSense(lb, ub);
    int r3 = 0, r1 = 0;
    double rhs =  sense == 'L' ? ub : lb;
    
    *useraction_p = CPX_CALLBACK_SET; //Tell CPLEX that cuts have been created
    
    //std::cout << "lb: " << lb << " ub: " << ub << " sense: " << sense << "\n";
    
    //std::cout << "wherefrom: " << wherefrom << "\n";
    
    //I am not sure if cplex requires Mutex to add lazy constraints. Anyway, we set...
    #if MRQ_SET_MUTEX_TO_ADD_LAZY_CONSTRAINTS
    SEMAPH_addConstr->lock(nthreads);
    #endif
    {
        if(sense == 'R')
        {
            //MRQ_PRINTERRORMSG("Unfortunatelly, cplex cannot handle double bound range constraints. So, you cannot provide a model having a double bounded nonlinear constraint.");
            
            //return MRQ_MILP_SOLVER_ERROR;
            
            //cplex cannot handle double bounded constraints here. So, we add two constraints; Note, in this case RHS is equal to lb
            
            //std::cout << "Adding constraint:\n";
            //printConstraint(nz, cols, vals, lb, ub);
            
            r1 = CPXcutcallbackadd(env, cbdata, wherefrom, nz, ub, 'L', cols, vals, purgeable);
            
            if(r1 != 0)
            {
                #if MRQ_DEBUG_MODE
                    MRQ_PRINTERRORNUMBER(r1);
                #endif
                *useraction_p = CPX_CALLBACK_FAIL;
            }
            
            sense = 'G';
        }
        
        //std::cout << "Adding constraint:\n";
        //printConstraint(nz, cols, vals, lb, ub);
        
        r3 = CPXcutcallbackadd(env, cbdata, wherefrom, nz, rhs, sense, cols, vals, purgeable);
        
    }
    #if MRQ_SET_MUTEX_TO_ADD_LAZY_CONSTRAINTS
    SEMAPH_addConstr->unlock(nthreads);
    #endif
    
    
    if(r1 != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r1);
        #endif
        *useraction_p = CPX_CALLBACK_FAIL;
        
        return r1;
    }
    
    
    if(r3 != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r3);
        #endif
        *useraction_p = CPX_CALLBACK_FAIL;
        
        return r3;
    }
    
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::getVarBoundsOnNode(const int index, double &lb, double &ub)
#if OPT_HAVE_CPLEX
{
    int r = CPXgetcallbacknodelb(env, cbdata, wherefrom, &lb, index, index);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    r = CPXgetcallbacknodeub(env, cbdata, wherefrom, &ub, index, index);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver:: getUserPointerAtCurrentNode(void **p)
#if OPT_HAVE_CPLEX
{
    const int currentNodeIndex = 0;
    
    *p = NULL;
    
    //the unique way that I could find is use the function to reset the user pointer. Remember after execute that, node will store NULL POINTER
    
    int r = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, currentNodeIndex, CPX_CALLBACK_INFO_NODE_USERHANDLE, p);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::getVarsLowerBoundsOnNode(const int n, double *lbx)
#if OPT_HAVE_CPLEX
{
    int r = CPXgetcallbacknodelb(env, cbdata, wherefrom, lbx, 0, n-1);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return r;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::getVarsUpperBoundsOnNode( const int n, double *ubx)
#if OPT_HAVE_CPLEX
{
    int r = CPXgetcallbacknodeub(env, cbdata, wherefrom, ubx, 0, n-1);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::getCurrentNodeLowerBound(double &nlb)
#if OPT_HAVE_CPLEX
{
    if( wherefrom == CPX_CALLBACK_MIP_BRANCH )
    {
        //I am not sure if CPXgetcallbacknodeobjval and CPXgetcallbacknodeinfo provide the same value. Since in the branching callbak in admipex1.c they have used CPXgetcallbacknodeobjval, I will use it here.
        int r = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &nlb);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    }
    else
    {
        const int thisnode = 0; //value 0 indicates we are refering to current node
        
        int r = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, thisnode, CPX_CALLBACK_INFO_NODE_OBJVAL, &nlb);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    }
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::getMIPDualBound(double &dualBound)
#if OPT_HAVE_CPLEX
{
    int r = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_REMAINING, &dualBound);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif



int MRQ_CplexCallbackSolver:: getNumberOfVariables(int &nvars)
#if OPT_HAVE_CPLEX
{
    int r;
    CPXCLPptr lp;
    
    r = CPXgetcallbacklp (env, cbdata, wherefrom, &lp);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    nvars = CPXgetnumcols (env, lp);
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver:: getIntegerVariableIndices(int &nIndices, int *indices)
#if OPT_HAVE_CPLEX
{
    int r, nvars, nI;
    CPXCLPptr lp;
    
    nIndices = 0;
    
    r = CPXgetcallbacklp (env, cbdata, wherefrom, &lp);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    nvars = CPXgetnumcols(env, lp);
    
    nI = CPXgetnumint(env, lp);
    
    if( nI > 0 )  //cplex return an error if we call CPXgetctype in a non integer problem. So, we have to perform this test
    {
        for(int i = 0; i < nvars; i++)
        {
            char vt;
        
            const int r = CPXgetctype(env, lp, &vt, i, i);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            if( vt == 'I')
            {
                indices[ nIndices ] = i;
                nIndices++;
            }
        }
    }
    
    #if MRQ_DEBUG_MODE
        assert( nI == nIndices);
    #endif
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver:: getPseudoCosts(int n, double *dowpcosts, double *uppcosts)
#if OPT_HAVE_CPLEX
{
    int r = CPXgetcallbackpseudocosts(env, cbdata, wherefrom, uppcosts, dowpcosts, 0, n-1);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::getVarPseudoCosts (const int index, double &dowpcost, double &uppcost)
#if OPT_HAVE_CPLEX
{
    int r = CPXgetcallbackpseudocosts(env, cbdata, wherefrom, &uppcost, &dowpcost, index, index);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver:: getNumberOfOpenNodes(long unsigned int &nnodes)
#if OPT_HAVE_CPLEX
{
    int r =  CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODES_LEFT_LONG, &nnodes);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    #if MRQ_DEBUG_MODE
        assert(nnodes > 0); //if cplex count current node like an open node, the number of nodes should be positive
    #endif
    
    nnodes--; //we subtract one here because the current node is counted by cplex like an open node
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver::pruneCurrentNode()
#if OPT_HAVE_CPLEX
{ 
    //to prune a node, we will turn the node infeasible changing the  first variable bounds 
    int r;
    const int firstVarIndex = 0;
    double lbFirstVar, ubFisrtVar;
    const double newlb = 1.0, newub = 0.0; //bounds to turn problem infeasible
    CPXLPptr nodelp;
    
    
    if( wherefrom == CPX_CALLBACK_MIP_BRANCH)
    {
        //since we are in a branching callback, we just do not generate any new node
        
        //Set useraction to indicate a user-specified branch
        *useraction_p = CPX_CALLBACK_SET;
        prunedNode = true;
    }
    else if( wherefrom == CPX_CALLBACK_MIP_SOLVE )
    {
        //Get pointer to LP subproblem
        r = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        
        //saving the original bounds of first variable
        r = CPXgetub( env, nodelp, &ubFisrtVar, firstVarIndex, firstVarIndex );
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        r = CPXgetlb( env, nodelp, &lbFirstVar, firstVarIndex, firstVarIndex );
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        //printf("limites originais da primeira variavel - l: %f u: %f\n", lbFirstVar, ubFisrtVar);
        
        //changing the bounds of first var to turn current node infeasible.
        r = CPXchgbds(env, nodelp, 1, &firstVarIndex, "L", &newlb);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        r = CPXchgbds(env, nodelp, 1, &firstVarIndex, "U", &newub);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        //solving the problem
        r = CPXprimopt (env, nodelp); //could be CPXdualopt (env, nodelp);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        *useraction_p = CPX_CALLBACK_SET; //telling cplex we already solve the node
        prunedNode = true;
        
        #if MRQ_DEBUG_MODE
            r = CPXgetstat(env, nodelp);
            assert(r == CPX_STAT_INFEASIBLE);
        #endif
        
        /*//restoring original bounds of the node
        r = CPXchgbds(env, nodelp, 1, &firstVarIndex, "L", &lbFirstVar);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        r = CPXchgbds(env, nodelp, 1, &firstVarIndex, "U", &ubFisrtVar);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);*/
    }
    else
    {
        MRQ_PRINTERRORMSGP("Non implemented action for wherefrom: ", wherefrom);
        return MRQ_NONIMPLEMENTED_ERROR;
    }
    
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


int MRQ_CplexCallbackSolver:: tryGenerateChildNodesSavingTheParentSolOnBranchVars(const double *sol)
#if OPT_HAVE_CPLEX
{
    //we only will save the parent solution if the branching type is from variable (cplex can perform other kinds of branching like SOS, Constraints)
    
    
    if( brtype == CPX_TYPE_VAR )
    {
        //int r;
        int nodeSeqNumber; //nodeSeqNumber will be filled by cplex when a node will be created.
        
        /*printf("nodecnt: %d\n", nodecnt);
        for(int i = 0; i < nodecnt; i++)
        {
            printf("nodebeg[%d]: %d nodeest[%d]: %lf\t", i, nodebeg[i], i, nodeest[i]);
        }
        printf("\n");
        for(int i = 0; i < bdcnt; i++ )
        {
            printf("i: %d - indices: %d  lu: %c  bd: %lf\n", i, indices[i], lu[i], bd[i]);
        }
        printf("\n"); */
        
        
        for(int i = 0; i < nodecnt; i++)
        {
            const unsigned int indStart = nodebeg[i];
            
            //const unsigned int indEnd = (i < nodecnt-1 ? nodebeg[i+1] : bdcnt ) - 1;
            
            const unsigned int sizeBounds = (i < nodecnt-1 ? nodebeg[i+1] : bdcnt) - indStart; //size of bounds
            
            
            //generating out array to store the parent solution
            
                MRQ_SizeIndexSol *spsol = new (std::nothrow) MRQ_SizeIndexSol;
                MRQ_IFMEMERRORRETURN( !spsol );
                
                int r3 = spsol->allocate(sizeBounds);
                MRQ_IFERRORRETURN(r3, r3);
            
                
                MRQ_IndexSol * psol = spsol->indexSol;
                
                for(unsigned int j = 0; j < sizeBounds; j++)
                {
                    unsigned int ind = indices[indStart + j];
                    psol[j].index = ind;
                    psol[j].sol = sol[ind];
                }
                
            
            
            
            int r2 = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, sizeBounds , &indices[indStart], &lu[indStart], &bd[indStart], nodeest[i], spsol, &nodeSeqNumber );
            
            if(r2)
            {
                MRQ_PRINTERRORNUMBER(r2); //do not return here because it would dangerous since previous node where already created.
                MRQ_getchar();
            }
        }
        
        *useraction_p = CPX_CALLBACK_SET;
        
    }
    else
    {
        return MRQ_BAD_PARAMETER_VALUES;
    }
    
    return 0;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


#if 0
MRQ_PointToStore::MRQ_PointToStore()
{
    id = -1;
    point = NULL;
}


void MRQ_PointToStore::desallocate()
{
    MRQ_secFree(point);
}


int MRQ_PointToStore::setPoint(long unsigned int id, unsigned int n, const double *point)
{
    this->id = id;
    
    MRQ_malloc(this->point, n);
    if( !this->point )
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTMEMERROR;
        #endif
        
        return MRQ_MEMORY_ERROR;
    }
    
    MRQ_copyArray(n, point, this->point);
    
    return 0;
}


MRQ_PointToStore::~MRQ_PointToStore()
{
    desallocate();
}

#endif


MRQ_PointsStoreKeeper::MRQ_PointsStoreKeeper(unsigned int nthreads)
{
    this->nthreads = nthreads;
}


MRQ_PointsStoreKeeper::~MRQ_PointsStoreKeeper()
{
    clear();
}


void MRQ_PointsStoreKeeper::clear()
{
    for( auto it = points.rbegin(); it != points.rend(); ++it  )
    {
        free( *it );
    }
    
    points.clear();
}


int MRQ_PointsStoreKeeper::insertPoint(unsigned int n, const double *point)
{
    int r, retCode = 0;
    double *p = NULL;
    
    MRQ_malloc(p, n);
    
    if(!p)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTMEMERROR;
        #endif
        
        retCode = MRQ_MEMORY_ERROR;
        goto termination;
    }
    
    MRQ_copyArray(n, point, p);
    
    r = 0;
    SEMAPH_access.lock(nthreads);
    {
        try
        {
            points.push_back( p );
        }
        catch(std::bad_alloc e)
        {
            r = MRQ_MEMORY_ERROR;
        }
    }
    SEMAPH_access.unlock(nthreads);
    
    if(r != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTMEMERROR;
        #endif
        
        return r;
    }
    
    
termination:
    
    if( retCode != 0 )
    {
        if(p)		free(p);
    }
    
    return retCode;
}


//this method get a solver object and put linearizations over the points stored from startId + 1. Only 
int MRQ_PointsStoreKeeper::updateLinearizationOnPoints( MRQ_MasterMILPProb &masterMilp, bool linearizeObj, double epsToActConstr, unsigned int *constrLinearsSaved, bool quadsInMaster, int constrLinStrategy, const bool *constrEval, int objLinStrategy, double zu, MRQ_LAAPointsStoring *laps)
{
    if( !points.empty() )
    {
        int r;
        const bool hasNlConstr = masterMilp.minlpProbhasNLConstraints();
        
        SEMAPH_access.lock(nthreads);
        {
            
            for( auto  it = points.rbegin();  it != points.rend(); ++it  )
            {
                double *point = *it;
                
                r = masterMilp.addLinearizedNLConstraintsByStrategy( epsToActConstr, constrLinearsSaved, true, point, quadsInMaster, constrLinStrategy, constrEval );
                
                if(r != 0)
                {
                    #if MRQ_DEBUG_MODE
                        MRQ_PRINTERRORNUMBER(r);
                    #endif
                    
                    break;
                }
                
                if( linearizeObj )
                {
                    r = masterMilp.addLinearizedObjFunction( !hasNlConstr, point, quadsInMaster, objLinStrategy, zu, laps );
                    
                    if(r != 0)
                    {
                        #if MRQ_DEBUG_MODE
                            MRQ_PRINTERRORNUMBER(r);
                        #endif
                        
                        break;
                    }
                }
                
            }
            
            clear();
        }
        SEMAPH_access.unlock(nthreads);
        
        
        //r = masterMilp.addLinearizedNLConstraintsByStrategy();
        if(r != 0)
        {
            return MRQ_MILP_SOLVER_ERROR;
        }
    }
    
    return 0;
}











MRQ_MILPCallbackData::MRQ_MILPCallbackData()
{
    initialize();
}


MRQ_MILPCallbackData::~MRQ_MILPCallbackData()
{
    desallocate();
}


int MRQ_MILPCallbackData::allocateBase(const unsigned int thnumber, const unsigned int nthreads, MRQ_MINLPProb* prob, const int milpSolver, const bool allocMaster, const MRQ_NLP_SOLVER *nlpSolver, MRQ_Mutex *SEMAPH_addConstr)
{
    const int n = prob->n;
    const int m = prob->m;
    
    int r, sizeAuxVars = (n+1);
    
    if(milpSolver == MRQ_GUROBI)
        sizeAuxVars += m; //Unfortunatelly, Gurobi generates auxiliary variables for double ranged constraints (I hate Gurobi).
    
    this->prob = prob;
    this->thnumber = thnumber;
    
    MRQ_malloc(auxConstrEval2, m);
    MRQ_malloc(auxVars, 3*sizeAuxVars);
    MRQ_malloc(auxConstr, 2*m); 
    MRQ_IFMEMERRORRETURN(!auxConstrEval2 || !auxVars || !auxConstr);
    
    auxConstr2 = &auxConstr[m];
    auxVars2 = &auxVars[sizeAuxVars];
    auxVars3 = &auxVars[2*sizeAuxVars];
    
    if(nlpSolver)
    {
        nlp = optsolvers::OPT_newLPSolver(*nlpSolver);
        nlpFeas = new (std::nothrow) MRQ_NLPFeasProb;
        MRQ_IFMEMERRORRETURN(!nlp || !nlpFeas);
    }
    
    if(allocMaster)
    {
        relaxMaster = new (std::nothrow) MRQ_MasterMILPProb;
        MRQ_IFMEMERRORRETURN(!relaxMaster);
    }
    
    r = gradEval.initialize(thnumber, prob);
    MRQ_IFERRORRETURN(r, MRQ_MEMORY_ERROR);
    
    callbackSolver = MRQ_newCallbackSolver(milpSolver, nthreads, SEMAPH_addConstr);
    MRQ_IFMEMERRORRETURN(!callbackSolver)
    
    return 0;
}



void MRQ_MILPCallbackData::desallocate()
{
    MRQ_secFree(auxConstrEval2);
    MRQ_secFree(auxVars);
    MRQ_secFree(auxConstr);
    
    MRQ_secDelete(nlp);
    
    gradEval.desallocate();
    
    MRQ_secDelete(nlpFeas);
    MRQ_secDelete(relaxMaster);
    MRQ_secDelete(callbackSolver);
}


void MRQ_MILPCallbackData::initialize()
{
    lastiter = -1;
    lastLinPoint = -1;
    
    intVars = NULL;
    constrEval = NULL;
    auxConstrEval2 = NULL;
    indices = NULL;
    auxVars = NULL;
    auxVars2= NULL;
    auxVars3= NULL;
    auxConstr = NULL;
    auxConstr2 = NULL;
    refSol = NULL;
    plc = NULL;
    puc = NULL;
    
    prob = NULL;
    relaxMaster = NULL;
    nlpFeas = NULL;
    
    nlp = NULL;
    
    alg = NULL;
    
    laps = NULL;
    
    callbackSolver = NULL;
    allPointsKeepers = NULL;
    
    out_sol_lower_than_zl = false;
    out_number_of_nlp_probs_solved = 0;
    out_cpu_time_of_nlp_solving = 0.0;
    out_clock_time_of_nlp_solving = 0.0;
}







//linearising objective 
int muriqui::MRQ_addLazyConstraintObjectiveLinearizationOnSolution( const int thnumber, MRQ_MILPSolverCallbackInterface *callbackSolver, MRQ_MINLPProb& prob, const bool incQuadsInMaster, const bool newx, const double *psol, const double *objSol, const int *indices, double *auxVars )
{
    const int n = prob.n;
    int r, retCode = 0;
    double rhs;
    
    
    r = MRQ_calculateObjLinearizedConstraint(prob, thnumber, newx, psol, incQuadsInMaster, objSol, NULL, auxVars, rhs);
    MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    //SEMAPH_addLazy.lock(nthreads);
    {
        r = callbackSolver->setLazyConstraint(n+1, indices, auxVars, -INFINITY, rhs);
    }
    //SEMAPH_addLazy.unlock(nthreads);
    MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_MILP_SOLVER_ERROR, termination);
    
    //std::cout << "linearizei funcao objetivo\n";
termination:
    
    return retCode;
}


//linearising constraints
int muriqui::MRQ_addLazyConstraintConstraintLinearizationOnSolution( const unsigned int thnumber, MRQ_MILPSolverCallbackInterface *callbackSolver, MRQ_MINLPProb &prob, MRQ_GradientsEvaluation &gradEval, const bool incQuadsInMaster, const bool *constrEval, bool newx, const double *psol, const double *pconstr, const double *masterConstr, const int *indices, const double *plc, const double *puc, double *auxVars, bool *auxConstrEval2, const int in_constr_linearisation_strategy, const double in_eps_to_active_constr_to_linearisation, unsigned int *out_number_of_constr_linears_saved)
{
    const int n = prob.n;
    const int m = prob.m;
    const bool *nlConstr = prob.nlConstr;
    MRQ_SparseMatrix *QC = prob.QC;
    double *grad = auxVars;
    
    int r, retCode = 0;
    
    
    //TODO: remove linearizations strategy by active and infeas in master solution...
    MRQ_calculateConstraintsToBeLinearizedByStrategy(prob, in_constr_linearisation_strategy, in_eps_to_active_constr_to_linearisation, constrEval, pconstr, masterConstr, plc, puc, out_number_of_constr_linears_saved, auxConstrEval2);
    
    if(prob.hasNlConstrs)
    {
        r = gradEval.evaluateJacobian(newx, auxConstrEval2, psol);
        MRQ_IFCALLBACKERRORGOTOLABEL(r, retCode, termination);
    }
    
    
    
    for(int i = 0; i < m; i++)
    {
        if(auxConstrEval2[i] == false)
            continue;
        
        if( nlConstr[i] == false && ( QC[i].getNumberOfElements() == 0u || incQuadsInMaster ) )
            continue;
        
        
        gradEval.constraintCompleteGradient(i, psol, grad);
        
        const double deltaRHS = MRQ_calculateDeltaRHSToConstraintLinearisation(pconstr[i], n, grad, psol);
        
        const double lc = plc[i] > -MIP_INFINITY ? plc[i] + deltaRHS : -INFINITY;
        const double uc = puc[i] < MIP_INFINITY ? puc[i] + deltaRHS : INFINITY;
        
        //const char sense = MRQ_cplexConstraintSense(plc[i], puc[i]);
        
        #if 0
        {
            const int ind = 10;
            
            if( i == ind  &&  !(  plc[ind] <= pconstr[ind] && pconstr[ind] <= puc[ind]   )     )
            {
                printf("rest: %f <= ", lc);
                for(int j = 0; j < n; j++)
                    printf("  %f", grad[j] );
                printf("  <= %f\n", uc);
                
                double v = 0.0;
                for(int j = 0; j < n; j++)
                    v = v  + grad[j]*psol[j] ;
                
                printf("valor da restrição na solução: %f dif: %f\n", v, uc - v);
                
                MRQ_getchar();
            }
        }
        #endif
        
        //SEMAPH_addLazy.lock(nthreads);
        {
            r = callbackSolver->setLazyConstraint(n, indices, grad, lc, uc);
        }
        //SEMAPH_addLazy.unlock(nthreads);
        MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_MILP_SOLVER_ERROR, termination);
        
    }
    
    
    
    
termination:
    
    return retCode;
}



