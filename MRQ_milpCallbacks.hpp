

#ifndef __MRQ_MILP_CALLBACKS_HPP__
#define __MRQ_MILP_CALLBACKS_HPP__

#include <iostream>
#include <list>

#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_solvers.hpp"

#define MRQ_SET_MUTEX_TO_ADD_LAZY_CONSTRAINTS 1


namespace muriqui
{


inline bool MRQ_isLazyConstraintsAvaliable(const int milpSolver)
{
    return milpSolver == MRQ_CPLEX || milpSolver == MRQ_GUROBI;
}


//class to store index and solution for a variable node. We create it to use when perform branchings in pseudo prunning, where it is necessary store optimal values of variables.
class MRQ_IndexSol
{
public:
    unsigned int index;
    double sol;
};


class MRQ_SizeIndexSol
{
public:
    unsigned short int size;
    MRQ_IndexSol *indexSol;
    
    ~MRQ_SizeIndexSol()
    {
        if(indexSol)    free(indexSol);
    }
    
    int allocate(unsigned short int size)
    {
        this->size = size;
        MRQ_malloc(indexSol, size);
        MRQ_IFMEMERRORRETURN(!indexSol);
        return 0;
    }
};



//class to abstract the proccess of set a lazy constraint
class MRQ_MILPSolverCallbackInterface
{
public:
    
    unsigned int nthreads;
    MRQ_Mutex *SEMAPH_addConstr;
    bool prunedNode;
    
    
    MRQ_MILPSolverCallbackInterface(unsigned int nthreads, MRQ_Mutex *SEMAPH_addConstr = NULL);
    
    virtual ~MRQ_MILPSolverCallbackInterface(){}
    
    virtual void initializeSolverData(const unsigned int nParams, void **solverParams) = 0;
    
    virtual int getNodeObjValue(double &objValue);
    
    virtual int getNodeSolution(const int n, double *sol) = 0;
    
    virtual int setLazyConstraint(const int nz, const int *cols, const double *vals, const double lb, const double ub) = 0;
    
    
    virtual int setCut(const int nz, const int *cols, const double *vals, const double lb, const double ub);
    
    virtual int setSolution(const int n, const double *sol);
    
    virtual int getVarsLowerBoundsOnNode(const int n, double *lbx);
    
    virtual int getVarsUpperBoundsOnNode(const int n, double *ubx);
    
    virtual int getNumberOfIterations(long int &niters);
    
    
    virtual int getUserPointerAtCurrentNode(void **p)
    {
        return MRQ_NONIMPLEMENTED_ERROR; 
    }
    
    virtual int initializeExtraBranchingSolverParameters(const unsigned int nParams, void **solverParams)
    {
        return MRQ_NONIMPLEMENTED_ERROR; 
    }
    
    
    //index of branching variable in the current node
    virtual int getNodeBranchVarIndex(int &nIndices, int *indices)
    {
        return MRQ_NONIMPLEMENTED_ERROR; 
    }
    
    
    virtual int getVarBoundsOnNode(const int index, double &lb, double &ub)
    {
        return MRQ_NONIMPLEMENTED_ERROR; 
    }
    
    //to get the number of variables in subproblems (pressolve )
    virtual int getNumberOfVariables(int &nvars)
    {
        return MRQ_NONIMPLEMENTED_ERROR; 
    }
    
    virtual int getPseudoCosts(int n, double *lowpcosts, double *uppcosts)
    {
        return MRQ_NONIMPLEMENTED_ERROR;
    }
    
    virtual int getVarPseudoCosts(const int index, double &dowpcost, double &uppcost)
    {
        return MRQ_NONIMPLEMENTED_ERROR;
    }
    
    virtual int getBestSolutionObjValue(double &objValue)
    {
        return MRQ_NONIMPLEMENTED_ERROR;
    }
    
    /*methods to muriqui server callbacks in any nodes*/
    
    virtual int getCurrentNodeLowerBound(double &nlb)
    {
        return MRQ_NONIMPLEMENTED_ERROR;
    }
    
    virtual int getMIPDualBound(double &dualBound)
    {
        return MRQ_NONIMPLEMENTED_ERROR;
    }
    
    virtual int getNumberOfOpenNodes(long unsigned int &nnodes)
    {
        return MRQ_NONIMPLEMENTED_ERROR;
    }
    
    virtual int pruneCurrentNode()
    {
        return MRQ_NONIMPLEMENTED_ERROR;
    }
    
    virtual int getIntegerVariableIndices(int &nIndices, int *indices)
    {
        return MRQ_NONIMPLEMENTED_ERROR;
    }
    
    //to be used in the branching callbacks
    virtual int tryGenerateChildNodesSavingTheParentSolOnBranchVars(const double *sol)
    {
        return MRQ_NONIMPLEMENTED_ERROR;
    }
    
    
    
protected:
    
    void printConstraint(const int nz, const int *cols, const double *vals, const double lb, const double ub);
    
};



class MRQ_GurobiCallbackSolver : public MRQ_MILPSolverCallbackInterface
{
    #if OPT_HAVE_GUROBI
    static inline char constraintSense(const double lb, const double ub, double &rhs);
    #endif
    
public:
    
    #if OPT_HAVE_GUROBI
        GRBmodel *model;
        void *cbdata;
        int where;
    #endif
    
    MRQ_GurobiCallbackSolver(unsigned int nthreads, MRQ_Mutex *SEMAPH_addConstr = NULL):MRQ_MILPSolverCallbackInterface(nthreads, SEMAPH_addConstr){}
    
    virtual ~MRQ_GurobiCallbackSolver(){}
    
    virtual void initializeSolverData(const unsigned int nParams, void **solverParams) override;
    
    virtual int getNumberOfIterations(long int &niters) override;
    
    virtual int getNodeObjValue(double &objValue) override;
    
    virtual int getNodeSolution(const int n, double *sol) override;
    
    virtual int setLazyConstraint(const int nz, const int *cols, const double *vals, const double lb, const double ub) override;
    
    virtual int setSolution(const int n, const double *sol) override;
};



class MRQ_GlpkCallbackSolver : public MRQ_MILPSolverCallbackInterface
{
    #if OPT_HAVE_GLPK
    static inline int constraintSense(const double lb, const double ub);
    #endif
    
    int allocateMyIndices(const unsigned int size);
    
    
public:
    
    #if OPT_HAVE_GLPK
        glp_tree *T;
        glp_prob *prob;
    #endif
    
    int *myindices; //unfortunatelly glpk starts counts from 1. So, we have to correct the indices
    
    
    MRQ_GlpkCallbackSolver(unsigned int nthreads, MRQ_Mutex *SEMAPH_addConstr = NULL);
    
    virtual ~MRQ_GlpkCallbackSolver();
    
    
    void desallocate();
    
    
    virtual void initializeSolverData(const unsigned int nParams, void **solverParams) override;
    
    virtual int getNodeSolution(const int n, double *sol) override;
    
    
    virtual int setLazyConstraint(const int nz, const int *cols, const double *vals, const double lb, const double ub) override;
};


class MRQ_CplexCallbackSolver : public MRQ_MILPSolverCallbackInterface
{
    
    static inline char constraintSense(const double lb, const double ub);
    
public:
    
    #if OPT_HAVE_CPLEX
        CPXCENVptr env;
        void *cbdata;
        int wherefrom;
        int *useraction_p;
        
        /* extra parameters for branching callback funcions */
        
        int brtype;
        int sos;
        int nodecnt;
        int bdcnt;
        const int *nodebeg;
        const int *indices;
        const char *lu;
        const double *bd;
        const double *nodeest;
        
        /* end of extra parameters for branching callback funcions */
    #endif
    
    MRQ_CplexCallbackSolver(unsigned int nthreads, MRQ_Mutex *SEMAPH_addConstr = NULL) : MRQ_MILPSolverCallbackInterface(nthreads, SEMAPH_addConstr){}
    
    virtual ~MRQ_CplexCallbackSolver(){}
    
    virtual void initializeSolverData(const unsigned int nParams, void **solverParams) override;
    
    virtual int initializeExtraBranchingSolverParameters(const unsigned int nParams, void **solverParams) override;
    
    virtual int getNumberOfIterations(long int &niters) override;
    
    virtual int getNodeObjValue(double &objValue) override;
    
    virtual int getBestSolutionObjValue(double &objValue) override;
    
    virtual int getNodeSolution(const int n, double *sol) override;
    
    virtual int getNodeBranchVarIndex(int &nIndices, int *indices) override;
    
    virtual int getVarsLowerBoundsOnNode(const int n, double *lbx) override;
    
    virtual int getVarsUpperBoundsOnNode(const int n, double *ubx) override;
    
    virtual int setLazyConstraint(const int nz, const int *cols, const double *vals, const double lb, const double ub) override;
    
    virtual int getCurrentNodeLowerBound(double &nlb) override;
    
    virtual int getMIPDualBound(double &dualBound) override;
    
    virtual int getNumberOfOpenNodes(long unsigned int &nnodes) override;
    
    virtual int getNumberOfVariables(int &nvars) override;
    
    virtual int getPseudoCosts(int n, double *dowpcosts, double *uppcosts) override;
    
    virtual int getVarPseudoCosts(const int index, double &dowpcost, double &uppcost) override;
    
    virtual int getIntegerVariableIndices(int &nIndices, int *indices) override;
    
    virtual int getVarBoundsOnNode(const int index, double &lb, double &ub) override;
    
    int getUserPointerAtCurrentNode(void **p) override;
    
    virtual int pruneCurrentNode() override;
    
    virtual int tryGenerateChildNodesSavingTheParentSolOnBranchVars( const double *sol) override;
};



static inline MRQ_MILPSolverCallbackInterface* MRQ_newCallbackSolver(const int milpSolver, unsigned int nthreads, MRQ_Mutex *SEMAPH_addConstr = NULL)
{
    if(milpSolver == MRQ_CPLEX)
    {
        return new (std::nothrow) MRQ_CplexCallbackSolver(nthreads, SEMAPH_addConstr);
    }
    else if(milpSolver == MRQ_GUROBI)
    {
        return new (std::nothrow) MRQ_GurobiCallbackSolver(nthreads, SEMAPH_addConstr);
    }
    else if(milpSolver == MRQ_GLPK)
    {
        return new (std::nothrow) MRQ_GlpkCallbackSolver(nthreads, SEMAPH_addConstr);
    }
    else
    {
        #if MRQ_DEBUG_MODE
            std::cerr << MRQ_PREPRINT "Invalid MILP solver " << milpSolver << MRQ_GETFILELINE << "\n";
        #endif
        return NULL;
    }
}

#if 0
class MRQ_PointToStore
{
public:
    long unsigned int id;
    double *point;
    
    MRQ_PointToStore();
    
    void desallocate(); 
    
    int setPoint(long unsigned int id, unsigned int n, const double *point);
    
    ~MRQ_PointToStore();
};
#endif

//class to store points in a Multithreading environmnet. This class has a mutex to control the storing. Its objective is to be used by MRQ_BBLP_ECP, where we have several threadings generating linearization points by the callbacks functions
class MRQ_PointsStoreKeeper
{
public:
    
    unsigned int nthreads;
    MRQ_Mutex SEMAPH_access;
    
    //linked list: first element is the MRQ_PointToStore. The second element is a counter to threaas that still needs get the pointer. When this counter get 0, we remove the node
    std::list< double* > points;
    
    MRQ_PointsStoreKeeper(unsigned int nthreads = 0);
    
    ~MRQ_PointsStoreKeeper();
    
    void clear();
    
    int insertPoint(unsigned int n, const double *point);
    
    //this method get a solver object and put linearizations over the points stored from startId + 1. Only 
    int updateLinearizationOnPoints( MRQ_MasterMILPProb &masterMilp, bool linearizeObj, double epsToActConstr, unsigned int *constrLinearsSaved, bool quadsInMaster, int constrLinStrategy, const bool *constrEval, int objLinStrategy, double zu, MRQ_LAAPointsStoring *laps );
    
};


class MRQ_MILPCallbackData
{
public:
    
    bool binProblem;
    bool linearizeObj;
    bool setQuadsInMaster;
    int nI;
    unsigned int thnumber;
    
    long int lastiter;
    long unsigned int lastLinPoint;  //last linearization point
    
    double timeStart;
    clock_t clockStart;
    
    
    const int *intVars;
    const bool *constrEval;
    bool *auxConstrEval2;
    const int *indices;
    double *auxVars, *auxVars2, *auxVars3;
    double *auxConstr;
    double *auxConstr2;
    double *plc, *puc;
    
    double *refSol; //reference solution. We use thism for example, to store the interior solution used by esh algorithm
    
    MRQ_MINLPProb *prob;
    
    MRQ_LPSolver *nlp;  //we use MRQ_LPSolver instead of nlpSOlver to take advantage this data structure to BB_LP_ECP_BASED
    MRQ_NLPFeasProb *nlpFeas;
    MRQ_MasterMILPProb *relaxMaster;
    
    MRQ_LinearApproxAlgorithm *alg;
    
    MRQ_LAAPointsStoring *laps;
    MRQ_GradientsEvaluation gradEval;
    
    MRQ_MILPSolverCallbackInterface *callbackSolver;
    
    MRQ_PointsStoreKeeper *allPointsKeepers;
    
    
    double out_cpu_time_of_nlp_solving;
    double out_clock_time_of_nlp_solving;
    bool out_sol_lower_than_zl;
    long unsigned int out_number_of_nlp_probs_solved;
    
    
    MRQ_MILPCallbackData();
    
    ~MRQ_MILPCallbackData();
    
    int allocateBase(const unsigned int thnumber, const unsigned int nthreads, MRQ_MINLPProb* prob, const int milpSolver, const bool allocMaster, const MRQ_NLP_SOLVER *nlpSolver, MRQ_Mutex *SEMAPH_addConstr);
    
    void desallocate();
    
    void initialize();
    
};



#if OPT_HAVE_CPLEX

//implemented in bonmin hybrid file implementation
int CPXPUBLIC MRQ_bonminHybCplexCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);

#endif




#if OPT_HAVE_GLPK
void MRQ_bboaGlpkCallback(glp_tree *T, void *info);
#endif

#if OPT_HAVE_GUROBI
int __stdcall MRQ_bboaGurobiCallback( GRBmodel *model, void *cbdata, int where, void *usrdata);
#endif

#if OPT_HAVE_CPLEX
int CPXPUBLIC MRQ_labbLazyConstraintCplexCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);

int CPXPUBLIC MRQ_labbCplexBeforeSolveCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);

//int CPXPUBLIC MRQ_labbCplexBranchingCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);

int CPXPUBLIC MRQ_labbCplexBranchingCallback(CPXCENVptr env, void *cbdata, int          wherefrom, void *cbhandle, int brtype, int sos, int nodecnt, int bdcnt, const int *nodebeg, const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p);

void CPXPUBLIC MRQ_labbCplexDeleteNodeCallback(CPXCENVptr env, int wherefrom, void *cbhandle, int seqnum,       void *handle);
#endif



//linearising objective 
int MRQ_addLazyConstraintObjectiveLinearizationOnSolution(const int thnumber, MRQ_MILPSolverCallbackInterface *callbackSolver, MRQ_MINLPProb& prob, const bool incQuadsInMaster, const bool newx, const double *psol, const double *objSol, const int *indices, double *auxVars );



//linearising constraints
int MRQ_addLazyConstraintConstraintLinearizationOnSolution(const unsigned int thnumber, MRQ_MILPSolverCallbackInterface *callbackSolver, MRQ_MINLPProb &prob, MRQ_GradientsEvaluation &gradEval, const bool incQuadsInMaster, const bool *constrEval, bool newx, const double *psol, const double *pconstr, const double *masterConstr, const int *indices, const double *plc, const double *puc, double *auxVars, bool *auxConstrEval2, const int in_constr_linearisation_strategy, const double in_eps_to_active_constr_to_linearisation, unsigned int *out_number_of_constr_linears_saved);




}



#endif
