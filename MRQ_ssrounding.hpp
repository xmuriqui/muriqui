/*
 * This file implements the structured stochastic rounding heurist. An experimental feasibility heuristic from our heads
 * 
 * 02, September, 2020
 * 
 * Author: Wendel Melo
 */ 

#ifndef __MRQ_SSROUNDING_HPP__
#define __MRQ_SSROUNDING_HPP__


#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"


namespace muriqui {



//class to perform the stochastic rounding looking for the problem structure
class MRQ_SSRCore
{
protected:
    inline void fixVarAt1(const int indVar, double *lx, double *ux, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *uc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars )
    {
        lx[indVar] = 1.0;
        ux[indVar] = 1.0;
        
        if(updateVarBoundsByKnapsackConstrs && binSumConstrInds && binSumConstrInds->knapsackInds)
        {
            const int indIntVar = reverseIntVars[indVar];
            
            const auto *kinds = binSumConstrInds->knapsackInds[indIntVar];
            auto nkinds = binSumConstrInds->nKnapsackInds[indIntVar];
            
            const minlpproblem::MIP_SparseMatrix &A = prob->A; 
            #if MRQ_DEBUG_MODE
                //const double *lc = prob->lc;
            #endif
            
            for( decltype(nkinds) k = 0; k < nkinds; k++)
            {
                const auto cind = kinds[k];
                
                const int* acols = A[cind];
                const double* avalues = A(cind);
                const unsigned int nel = A.getNumberOfElementsAtRow(cind);
                
                double rhs = uc[cind];
                double highestFreeVarCoef = -INFINITY;
                
                
                //we can run this inside branch and bound and preprocessor can disable this constraint
                if( uc[cind] >= MIP_INFINITY  )
                    continue;
                
                #if MRQ_DEBUG_MODE
                    //assert( lc[cind] <= -MIP_INFINITY && uc[cind] < MIP_INFINITY);
                #endif
                
                    
                for(unsigned int j = 0; j < nel; j++)
                {
                    const auto col = acols[j];
                    const auto coef = avalues[j];
                    
                    if( lx[col] == ux[col])
                    {
                        if( lx[col] != 0.0 )
                            rhs -= coef * lx[col];
                    }
                    else
                    {
                        if( coef > highestFreeVarCoef )
                            highestFreeVarCoef = coef;
                    }
                }
                
                //now, we check if some variable should be set to zero
                if( highestFreeVarCoef > rhs )
                {
                    for(unsigned int j = 0; j < nel; j++)
                    {
                        const auto col = acols[j];
                        const auto coef = avalues[j];
                        
                        if( lx[col] != ux[col] && coef > rhs)
                        {
                            ux[col] = 0.0;
                            #if MRQ_DEBUG_MODE
                                assert(lx[col] == 0.0);
                            #endif
                                
                            //std::cout << "Setting " << col << " to zero when fix " << indVar << " to 1 by constraint " << cind << "!\n";
                           // MRQ_getchar();
                        }
                    }
                }
                
            } //end of for( decltype(nkinds) k = 0; k < nkinds; k++)
            
        } //end of if(updateVarBoundsByKnapsackConstrs && binSumConstrInds)
        
    }
    
    
    int fixAllNonFixedVarsByStochRounding(const int n, const int nI, const int *intVars, const double *startSol, const double minimumProbToRound, const int print_level, MRQ_Random &random, int *auxInds, int &nVarsFixed, double *lx, double *ux, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager,  MRQ_BoundsUpdaterSolver *boundsUpdaterSolver, double *currlc, double *curruc, int contRelaxStrategy, optsolvers::OPT_LPSolver *nlpSolver);
    
    
    void fixAllNonFixedAt0(const int sizea, const int *acols, double *lx, double *ux);
    
    
    //fix first nVarsToFix non fixed variables to 1. We assume all nonfixed variable are binary. Returns the number of variables fixed
    int fixFirstVarsAt0( const int nVarsToFix, const int sizea, const int *acols, double *lx, double *ux );
    
    
    //fix first nVarsToFix non fixed variables to 1. We assume all nonfixed variable are binary. Returns the number of variables fixed (can be not possible fix nVarsToFix).
    int fixFirstVarsAt1( const int nVarsToFix, const int sizea, const int *acols, double *lx, double *ux, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *uc, const minlpproblem:: MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars );
    
    
    //Returns number of variables fixed
    int fixRandomVarsAt0(MRQ_Random &random, const int nVarsToFix, const int sizea, const int *acols, double *lx, double *ux, const int *reverseIntVars);
    
    //pick nVarsToFix nonfixed vars randomly and fix to 1. We assume all nonfixed variable are binary. Returns the number ofa variable fixed.
    int fixRandomVarsAt1(MRQ_Random &random, const int nVarsToFix, const int sizea, const int *acols, double *lx, double *ux, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *uc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars );
    
    //Returns number of variables fixed
    int fixRandomVarsAt01(MRQ_Random &random, const int nVarsToFix, const int sizea, const int *acols, double *lx, double *ux, const bool fixAt0, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *uc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars );
    
    
    
    
    //new version where we let some variables free in <= constraints
    virtual int class1_2_5_6Handle(const int sizea, const int *acols, const double *avalues, const double lc, const double uc, MRQ_Random &random, int *auxInds1, int *auxInds2, double *lx, double *ux, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *ubc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager );

    
    //new version where we let some variables free in <= constraints
    virtual int class0_3_4Handle(const int sizea, const int *acols, const double *avalues, const double lc, const double uc, MRQ_Random &random, int *auxInds, double *lx, double *ux,  const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *ubc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager );
    
    
    
    int class8Handle(const int sizea, const int *acols, const double *avalues, const double lc, const double uc, MRQ_Random &random, int *auxInds, double *lx, double *ux,  const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *ubc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars);
    
    
    //generic function to handle constraint from any classe. This function redirects  to the correct Handle based on constraint class
    inline int classXHandle(const int constrClass, const int sizea, const int *acols, const double *avalues, const double lc, const double uc, MRQ_Random &random, int *auxInds1, int *auxInds2, double *lx, double *ux, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *ubc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager )
    {
        switch(constrClass)
        {
            case 0:
            case 3:
            case 4:
                return class0_3_4Handle(sizea, acols, avalues, lc, uc, random, auxInds1, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars, preprocessor, ccstorager);
                break;
            
            case 1:
            case 2:
            case 5:
            case 6:
                return class1_2_5_6Handle(sizea, acols, avalues, lc, uc, random, auxInds1, auxInds2, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars, preprocessor, ccstorager);
                break;
            
            case 8:
                return class8Handle(sizea, acols, avalues, lc, uc, random, auxInds1, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
                break;
            
            default:
                assert(false);
                return MRQ_NONIMPLEMENTED_ERROR;
        }
    }
    

public:
    
    
    //if a preprocessor object is passed, we perform a preproessing after we threat each constraint
    int strucStochRounding(const MRQ_MINLPProb &prob, const double *startSol, const double *lx, const double *ux, const double *lc, const double *uc, const int nI, const int *intVars, const int *reverseIntVars, const bool randomOrderToThreatClasses, const bool randomOrderInsideClasses, const double minimumProbToRound, const int print_level, const minlpproblem::MIP_BinSumConstrsIndsByClass &binSumConstrInds, MRQ_Random &random, int *auxInds1, int *auxInds2, unsigned int *auxConstrInds, double *auxConstrs1, double *auxConstrs2, double *outlx, double *outux, int &nVarsFixedByStochRounding, MRQ_SSR_VARIABLES_ADDITIONAL_FIXING_STRATEGY in_vars_additional_fixing_strategy, MRQ_SSR_PREPROCESSING_STRATEGY preprocStrategy, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager, const bool preprocessAfterVarRounding, const int contRelaxStrategy, optsolvers::OPT_LPSolver &nlpSolver, MRQ_BoundsUpdaterSolver *boundsUpdaterSolver );
};


class MRQ_SSRCore2: public MRQ_SSRCore
{

protected:

    //new version where we let some variables free in <= constraints
    virtual int class1_2_5_6Handle(const int sizea, const int *acols, const double *avalues, const double lc, const double uc, MRQ_Random &random, int *auxInds1, int *auxInds2, double *lx, double *ux, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *ubc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager ) override;

    
    //new version where we let some variables free in <= constraints
    virtual int class0_3_4Handle(const int sizea, const int *acols, const double *avalues, const double lc, const double uc, MRQ_Random &random, int *auxInds, double *lx, double *ux,  const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *ubc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager ) override;

};


//core of Structured Stochastic Rounding
class MRQ_SSRoundingExecutor
{
    
protected:
    
    
    unsigned int *auxConstrInds;
    int *auxVarInds1, *auxVarInds2;
    double *auxVarValues1, *roundedLx, *roundedUx;
    double *auxConstrs1, *auxConstrs2; 
    
    MRQ_MINLPProb *subProb;
    
    MRQ_BoundsUpdaterSolver *boundsUpdaterSolver;
    
    
    int allocateAuxMemory(const unsigned int n, const unsigned int m, const unsigned int sizeAuxConstrIndex);
    
    void deallocate();
    
public:
    
    //bool in_preprocess_after_handling_constraints;
    bool in_preprocess_after_variable_rounding;
    bool in_random_order_to_threat_classes;
    
    bool in_random_order_to_threat_constraints_in_each_class;
    
    bool in_solve_minlp_as_local_search;
    bool in_stop_local_search_solving_on_first_improvment_solution;
    
    unsigned int in_max_number_of_main_iterations;
    unsigned int in_max_number_of_improvments;
    
    int in_print_level;
    MRQ_NLP_SOLVER in_nlp_solver;
    MRQ_MILP_SOLVER in_milp_solver;
    
    MRQ_SUBPROBLEM_NEIGHBORHOOD_STRATEGY in_neighborhood_strategy;
    
    MRQ_SSR_PREPROCESSING_STRATEGY in_preprocessing_point_strategy;
    MRQ_SSR_VARIABLES_ADDITIONAL_FIXING_STRATEGY in_additional_vars_fixing_strategy;
    
    MRQ_SSR_VARIABLE_BOUNDS_UPDT_STRATEGY in_rounding_var_bounds_updt_strategy;
    MRQ_SSR_CONTINUOUS_RELAXATION_STRATEGY_TO_STOCHASTIC_ROUNDING in_cont_relax_strategy_to_stoch_rounding;
    
    
    double in_absolute_feasibility_tol;
    double in_relative_feasibility_tol;
    
    double in_integer_tol;
    double in_integer_neighborhood_factor; //factor to define local branching constraint (only if you do not use MRQ_SNS_ORIGINAL as in_neighborhood_strategy)
    double in_continuous_neighborhood_factor; //factor to define local branching constraint (only if you use MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD as in_neighborhood_strategy). Radius of the neighborhood will be defined line in_continuous_neighborhood_factor*nC*alpha, where nC is the number of continuous variables and alpha is the greatest absolute value in continuous initial solution
    
    double in_relative_convergence_tol_to_local_search; //we will impose the local search only stop when it finds a solution at least in_fator_to_stop_local_search_on_first_improvment lower than the current best solution
    
    double in_max_cpu_time;
    double in_max_time;
    
    double in_min_probability_to_round;
    
    MRQ_ALG_CODE in_local_search_algorithm;
    
    MRQ_GeneralSolverParams* in_milp_solver_params;
    MRQ_GeneralSolverParams* in_nlp_solver_params;
    MRQ_GeneralSolverParams* in_alg_params;
    
    MRQ_Algorithm *in_alg; //if user put a pointer here, we use to aply algorithm. User has still the ownership of the pointer
    
    bool out_vars_fixed_by_stoch_rounding;
    
    unsigned int out_number_of_main_iterations;
    unsigned int out_number_of_improvments;
    
    double out_cpu_time_at_nlp_integer_fixed_sol;
    double out_obj_at_nlp_integer_fixed_sol;
    
    
    //MRQ_Algorithm *out_alg;
    MRQ_ALG_CODE out_local_search_alg_code;
    
    
    MRQ_SSRoundingExecutor();
    
    ~MRQ_SSRoundingExecutor();
    
    
    void resetOutput();
    
    void resetParameters();
    
    /*
     * 
     * lx and ux will not be changed. Unfortunatelly, we could not declare them as const for practical reazons (run MRQ_Algorithms by inside requires non const arrays, but they will not changed)
     * 
     * nlpSolver is a pointer to an solver object which nlp relaxation is already set.
     */
    int run(const MRQ_MINLPProb &prob, const minlpproblem::MIP_BinSumConstrsIndsByClass &binSumConstrInds, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager,  MRQ_Random &random, optsolvers::OPT_LPSolver &nlpSolver, unsigned int thnumber, unsigned int nThreads, double insideSolverMaxTime, const double *lc, const double *uc, double *lx, double *ux, const int nI, const int *intVars, const int *reverseIntVars, const int nC, const int *contVars, const double *relaxSol, int &algReturnCode, double &outObj, double *outSol);
    
};




    
}

#endif
