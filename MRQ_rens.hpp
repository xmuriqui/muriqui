
#ifndef __MRQ_RENS_HPP__
#define __MRQ_RENS_HPP__

#include "MRQ_algClasses.hpp"


namespace muriqui {


class MRQ_RensExecutor
{
private:
    
    int rensNeighStrategy;
    int neighConstrIndex;
    int *auxInds;
    int nI, *intVars;
    int nC, *contVars;
    double *auxValues;
    double *nlx, *nux;
    minlpproblem::MIP_NonLinearEval *newEval; //actually, newEval is a MIP_EncapsulatedNonLinearEval;
    MRQ_MINLPProb *mysubprob, *subprob;
    MRQ_HeuristicExecutor *heurs;
    
    
    //int setRensNeighborhood(const double *inputlx, const double *inputux, const double *sol, double *nlx, double *nux);
    
public:
    
    bool in_apply_only_heuristics_on_subproblems; //apply several heuristcs
    bool in_solve_continuous_relaxation_on_subproblem; //solve continuous relaxation for the subproblem, considering the reduced feasible region. Since we solve continuous relaxation for the original problem, we coud skip this step, but maybe the new constraints can help preprocessor to reduce bounds and get a better solution
    bool in_stop_on_first_sol; //stop execution in the first feasible solution founded
    
    int in_print_level;
    MRQ_NLP_SOLVER in_nlp_solver;
    MRQ_MILP_SOLVER in_milp_solver;
    
    double in_integer_tol;
    double in_integer_neighborhood_factor; //factor to define local branching constraint (only if you do not use MRQ_SNS_ORIGINAL as in_neighborhood_strategy)
    double in_continuous_neighborhood_factor; //factor to define local branching constraint (only if you use MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD as in_neighborhood_strategy). Radius of the neighborhood will be defined line in_continuous_neighborhood_factor*nC*alpha, where nC is the number of continuous variables and alpha is the greatest absolute value in continuous initial solution
    double in_max_cpu_time;
    double in_max_time;
    
    MRQ_ALG_CODE in_algorithm_to_solve_subproblem;
    
    MRQ_Algorithm *in_alg; //if user put a pointer here, we use to aply algorithm. User has still the ownership of the pointer.
    MRQ_HeuristicExecutor *in_heurs; //if user put a pointer here, we use to aply heuristcs
    
    
    MRQ_Algorithm *out_alg;
    
    
    MRQ_RensExecutor();
    
    ~MRQ_RensExecutor();
    
    void deallocate();
    
    void deallocateOutAlg();
    
    int setSubProb(MRQ_MINLPProb &prob, int rensNeighStrategy);
    
    void resetOutput();
    
    void resetParameters();
    
    int run(unsigned int thnumber, unsigned int nThreads, double insideSolverMaxTime, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, MRQ_GeneralSolverParams* algParams, const double *lx, const double *ux, const double *sol, double upper_bound, double &outObj, double *outSol);
};



int MRQ_setSubproblemNeighborhood(const double *inputlx, const double *inputux, const double *sol, const int nI, const int *intVars, const int nC, const int *contVars, const int neighStrategy, const double in_integer_neighborhood_factor, const double in_continuous_neighborhood_factor, const double in_integer_tol, const int neighConstrIndex, int *auxInds, double *auxValues, MRQ_MINLPProb *prob, double *nlx, double *nux);



}

#endif
