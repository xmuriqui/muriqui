/*
* That file implements the Integrality Gap Minimization Algorithm, version 2, by Wendel Melo, Marcia Fampa and Fernanda Raupp
* 
* Author: Wendel Melo
* 
* Date: 10-December-2016
* 
* */


#include "BBL_branchAndBound.hpp"
#include "BBL_tools.hpp"
#include "muriqui.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_bb.hpp"



typedef branchAndBound::BBL_Mutex MRQ_Mutex;



namespace muriqui 
{


#if OPT_HAVE_IPOPT

class MRQ_IpoptIntermediateCallback : public optsolvers::OPT_IpoptIntermediateCallback
{
    double lastObj;
    double lastInf_pr;
    
    
public:
    
    //parameters
    
    int minimumItersToabortIfNoprogress; //if we reached this number of iterations and have no progreess in objective and feasibility, we abort the procedure
    
    //if ipopt reaches the minimumItersToabortIfNoprogress and improvment about objective is lower than minImprovToObj and improvment about infeasibility is lower than minImprovToInfeas, we abort ipopt solving
    double minImprovToObj;
    double minImprovToInfeas;
    
    double integer_tol;
    double absFeasTol;
    double relFeasTol;
    
    int nI, thnumber;
    const int *intVars;
    optsolvers::OPT_Ipopt *ipopt;
    
    bool *constrEval;
    double *sol;
    double *constrs;
    
    const MRQ_MINLPProb *prob;
    
    
    MRQ_IpoptIntermediateCallback();
    
    ~MRQ_IpoptIntermediateCallback();
    
    void initialize( const int thnumber = 0, const MRQ_MINLPProb* prob = NULL, const double integerTol = 0.0, const double absFeasTol = 0.0, const double relFeasTol = 0.0, const int nI = 0, const int* intVars = 0);
    
    int allocate(const int n, const int m);
    
    void desallocate();
    
    /* That callback is to be used in the Integrality gap minimization process with Ipopt. We just check if we can stop since we have a integer feasible solution
    *
    */ 
    virtual bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value, Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq);
    
};


/* In this intermediate callback for ipopt, we just stop when we have a integer feasible solution, even if this integer solution is not feasible. We hope correct this by the local search solving the nlp fixing on integer vars.
* 
*/
class MRQ_IpoptIntermediateCallback2 : public MRQ_IpoptIntermediateCallback{
    
public:
    
    virtual bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value, Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq);
    
};



/* In this intermediate callback for ipopt, we just stop when we have a integer feasible solution, even if this integer solution is not feasible. We hope correct this by the local search solving the nlp fixing on integer vars.
* 
*/
class MRQ_IpoptIntermediateCallback3 : public MRQ_IpoptIntermediateCallback{
    
public:
    
    bool nlpSolved;
    double zu;
    MRQ_NLPSolver *nlp;
    
    MRQ_IpoptIntermediateCallback3();
    
    virtual bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value, Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq);
    
};



#endif




class MRQ_IGMA1BBCallbacks : public branchAndBound::BBL_UserCallbacks
{
    
    int allocateThreadStructures(const unsigned int nthreads, double* lx, double* ux);
    
    void updateSumGapObj( const double *sol, const double *nlx, const double *nux );
    
public:
    
    unsigned int nthreads;
    unsigned int nnonimprovs;
    branchAndBound::BBL_BRANCH_STRATEGY origBBLBranchStrategy;
    MRQ_BB_CONSTRAINT_BRANCH_STRATEGY constrBranchStrat;
    int nI;
    int *intVars;
    //unsigned int nimprovs; //number of impovments for the best sol
    
    double abseps;
    double releps;
    
    clock_t clockStart;
    double timeStart;
    
    
    MRQ_IGMA1 *igma1;
    MRQ_MINLPProb *prob;
    MRQ_GeneralSolverParams *gapMinSolverParams;
    MRQ_GeneralSolverParams *nlpSolverParams;
    
    double *oplc; //we do not save a pointer to puc. We assume puc is allocated m positions after plc
    
    
    double *sumGapObj; //to objective in minimzation gaps if MRQ_IGMA2_GMOS_SAME_WEIGHT is used
    unsigned int *nGapObj; //to objective in minimzation gaps if MRQ_IGMA2_GMOS_SAME_WEIGHT is used
    
    MRQ_BinSumConstrsInds binSumConstrs;
    
    
    //thread structures
    
    unsigned int *prunesByObjCutActive;
    unsigned int *nimprovments;
    //unsigned int *nnonimprovments; //when we do not adopt objective cut in the heuristic mode
    unsigned long int *nsubiters;
    double **tplc;
    MRQ_GapMinProb *gapmins;
    MRQ_NLPSolver **nlps;
    MRQ_Preprocessor *preprocessors;
    MRQ_Random *randoms;
    MRQ_BinSumConstrsChooser *constrChoosers;
    #if OPT_HAVE_IPOPT
        MRQ_IpoptIntermediateCallback3 *ipoptInterCall;
    #endif
    
    MRQ_Mutex SEMAPH_sumGapObj;
    
    
    MRQ_IGMA1BBCallbacks( MRQ_MINLPProb *prob, MRQ_IGMA1 *igma2, MRQ_GeneralSolverParams *gapMinParams, MRQ_GeneralSolverParams *nlpParams);
    
    ~MRQ_IGMA1BBCallbacks();
    
    
    void desallocate();
    
    
    void initialize(MRQ_MINLPProb *prob, MRQ_IGMA1 *igma2, MRQ_GeneralSolverParams* gapMinSolverParams, MRQ_GeneralSolverParams *nlpParams);
    
    
    
    
    
    
    virtual int beforeAll(const unsigned int numberOfThreads, double *lx, double *ux) override;
    
    
    virtual int beforeSolvingRelaxation( const unsigned int threadNumber, branchAndBound::BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, bool &pruneNode) override;
    
    
    virtual int solveSubProblem(const unsigned int threadNumber, branchAndBound::BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, branchAndBound::BBL_RETURN_CODES &retCode, double &objValue, double &dualObjValue, double *sol, double *dualSol, bool &generalFeasibleSol, bool &pruneNode, double &nodeLowerBound, branchAndBound::BBL_BRANCH_STRATEGY &branchStrategy) override;
    
    
    virtual int chooseIndexToBranch(const int threadNumber, branchAndBound::BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, branchAndBound::BBL_RETURN_CODES retCode, double objValue, double dualObjValue, double *sol, double *dualSol, unsigned int &sizeIndices, unsigned int *indices, double *breakValues1, double *breakValues2, branchAndBound::BBL_Node* &nodes) override;
    
    
    virtual int generateNodes(const int thnumber, branchAndBound::BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, branchAndBound::BBL_RETURN_CODES retCode, const double objValue, double *sol, double *dualSol, branchAndBound::BBL_UserNodeGenerator &userNodeGenerator ) override;
    
    
    virtual int endOfIteration(const int threadNumber, const long unsigned int iter, const double cpuTime, const double wallTime, const double lb, const double ub, branchAndBound::BBL_Node &node, const double *nlx, const double *nux) override;
    
    
    virtual int updatingBestSolution(const int threadNumber, double* sol, double &objValue, const double ub, const long unsigned int iter) override;
    
    
    //virtual void afterAll(const long unsigned int iters, const double cpuTime, const double wallTime, const double lb, const double ub);
    
};


}























