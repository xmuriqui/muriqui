
#include <cstring>

#include <string>

#include "MRQ_constants.hpp"
#include "MRQ_tools.hpp"

using namespace muriqui;






MRQ_RETURN_CODE muriqui::MRQ_intToReturnCode( const int returnCode)
{
    //warning:: MRQ_SUCCESS constant wont work with this function, because it will consider optimal solution
    switch( returnCode )
    {
        MRQ_CASE_OF_SWITCH(MRQ_OPTIMAL_SOLUTION);
        MRQ_CASE_OF_SWITCH(MRQ_INFEASIBLE_PROBLEM);
        MRQ_CASE_OF_SWITCH(MRQ_UNBOUNDED_PROBLEM);
        MRQ_CASE_OF_SWITCH(MRQ_ALG_NOT_APPLICABLE);
        MRQ_CASE_OF_SWITCH(MRQ_BAD_DEFINITIONS);
        MRQ_CASE_OF_SWITCH(MRQ_BAD_PARAMETER_VALUES);
        MRQ_CASE_OF_SWITCH(MRQ_INDEX_FAULT);
        MRQ_CASE_OF_SWITCH(MRQ_MEMORY_ERROR);
        MRQ_CASE_OF_SWITCH(MRQ_UNDEFINED_ERROR);
        MRQ_CASE_OF_SWITCH(MRQ_MAX_TIME_STOP);
        MRQ_CASE_OF_SWITCH(MRQ_MAX_ITERATIONS_STOP);
        MRQ_CASE_OF_SWITCH(MRQ_CALLBACK_FUNCTION_ERROR);
        MRQ_CASE_OF_SWITCH(MRQ_STOP_REQUIRED_BY_USER);
        MRQ_CASE_OF_SWITCH(MRQ_MILP_SOLVER_ERROR);
        MRQ_CASE_OF_SWITCH(MRQ_NLP_SOLVER_ERROR);
        MRQ_CASE_OF_SWITCH(MRQ_LIBRARY_NOT_AVAILABLE);
        MRQ_CASE_OF_SWITCH(MRQ_NAME_ERROR);
        MRQ_CASE_OF_SWITCH(MRQ_VALUE_ERROR);
        
        MRQ_CASE_OF_SWITCH(MRQ_MINLP_SOLVER_ERROR);
        
        MRQ_CASE_OF_SWITCH(MRQ_INITIAL_SOLUTION_ERROR);
        MRQ_CASE_OF_SWITCH(MRQ_HEURISTIC_FAIL);
        MRQ_CASE_OF_SWITCH(MRQ_LACK_OF_PROGRESS_ERROR);
        MRQ_CASE_OF_SWITCH(MRQ_NLP_NO_FEASIBLE_SOLUTION);
        MRQ_CASE_OF_SWITCH(MRQ_NONIMPLEMENTED_ERROR);
        
        
        //specially for heuristics like Feasibility Pump and Diving);
        MRQ_CASE_OF_SWITCH(MRQ_HEURISTIC_SUCCESS);
        MRQ_CASE_OF_SWITCH(MRQ_CONT_RELAX_OPTIMAL_SOLUTION);
        //specially for early branching);
        MRQ_CASE_OF_SWITCH(MRQ_NLP_EARLY_STOP);
        //specially for NLP solving
        MRQ_CASE_OF_SWITCH(MRQ_NLP_FEASIBLE_SOLUTION);
        //specially for MILP solving
        MRQ_CASE_OF_SWITCH(MRQ_MILP_FEASIBLE_SOLUTION);
        //specially for DEV);
        MRQ_CASE_OF_SWITCH(MRQ_DEV_CONVERGENCE);
        //internal code for nlp/nlp bb
        MRQ_CASE_OF_SWITCH(MRQ_LAZY_MILP_BB_SOLUTION_LOWER_THAN_ZL);
        
        default:
            return MRQ_UNDEFINED_ERROR;
    }
}


MRQ_NLP_SOLVER muriqui::MRQ_intToNLPSolver( const int nlpSolver)
{
    switch( nlpSolver )
    {
        MRQ_CASE_OF_SWITCH(MRQ_NLP_MOSEK);
        
        MRQ_CASE_OF_SWITCH(MRQ_IPOPT);
        
        MRQ_CASE_OF_SWITCH(MRQ_NLP_KNITRO);
        
        MRQ_CASE_OF_SWITCH(MRQ_ALGENCAN);
            
        MRQ_CASE_OF_SWITCH(MRQ_WORHP);
        
        default:
            std::cerr << MRQ_PREPRINT << "NLP solver " << nlpSolver << " not coded " << MRQ_GETFILELINE << "\n";
            //assert(false);
            return MRQ_UNDEFINED_NLP;
    }
}


MRQ_MILP_SOLVER muriqui::MRQ_intToMILPSolver( const int milpSolver)
{
    switch( milpSolver )
    {
        MRQ_CASE_OF_SWITCH(MRQ_GLPK);
        
        MRQ_CASE_OF_SWITCH(MRQ_CPLEX);
        
        MRQ_CASE_OF_SWITCH(MRQ_GUROBI);
            
        MRQ_CASE_OF_SWITCH(MRQ_XPRESS);
            
        MRQ_CASE_OF_SWITCH(MRQ_MILP_MOSEK);
            
        default:
            std::cerr << MRQ_PREPRINT << "Warning: MILP solver " << milpSolver << " not coded " << MRQ_GETFILELINE << "\n";
            //assert(false);
            return MRQ_UNDEFINED_MILP;
    }
}





int muriqui::MRQ_strToEnum(const char *svalue, MRQ_ALG_CODE &value)
{
    int ret = 0;
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_UNDEFINED_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_OA_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_NLP_BB_OA_BASED_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_ECP_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_IGMA0_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_IGMA1_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_IGMA2_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LOC_BRANCH_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_FP_HEUR_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_DIVE_HEUR_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_OA_FP_HEUR_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_ESH_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_CONT_RELAX_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BONMIN_HYBRID_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_RENS_HEUR_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_BB_ECP_BASED_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_BB_ESH_BASED_ALG), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_SSR_HEUR_ALG), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}







int muriqui::MRQ_strToEnum(const char *svalue, MRQ_MILP_SOLVER &value)
{
    int ret = 0;	
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_UNDEFINED_MILP), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_GLPK), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_CPLEX), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_GUROBI), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_XPRESS), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_MILP_MOSEK), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_LP_SOLVER &value)
{
    int ret = 0;	
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_UNDEFINED_LP), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_GLPK), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_CPLEX), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_GUROBI), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_XPRESS), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_MOSEK), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_IPOPT), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_WORHP), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_KNITRO), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_LP_ALGENCAN), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_NLP_SOLVER &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_UNDEFINED_NLP), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_NLP_MOSEK), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_NLP_KNITRO), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_IPOPT), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_WORHP), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_ALGENCAN), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_OPTIZELLE), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_GLOBAL_SOLVER &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_UNDEFINED_GLOBAL), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_IQUAD), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}



int muriqui::MRQ_strToEnum(const char *svalue, MRQ_BB_EXP_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_ES_DEPTH), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_ES_WIDTH), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_ES_BEST_LIMIT), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_ES_DEPTH_BEST_LIMIT), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}



int muriqui::MRQ_strToEnum(const char *svalue, MRQ_CONSTR_LIN_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_CLS_ALL_CONSTRS), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_CLS_ONLY_INFEAS_AND_ACTIVE), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_MASTER_SOL_ALSO), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}



int muriqui::MRQ_strToEnum(const char *svalue, MRQ_OBJ_LIN_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_OLS_ALL_POINTS), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_OLS_NON_OBJ_CUT_POINTS), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}



int muriqui::MRQ_strToEnum(const char *svalue, MRQ_QUAD_APP_MASTER_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_QAMS_NO_QUAD_APP), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_QAMS_ON_BEST_POINT), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_QAMS_ON_LAST_POINT), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}



int muriqui::MRQ_strToEnum(const char *svalue, MRQ_BB_BRANCH_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_BS_HIGHEST_INT_GAP), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_BS_BIN_FIRST_HIGHEST_INT_GAP), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_BS_STBRANCH_PSEUDO_COSTS), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_BS_VAR_PRIORITIES), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_BS_USER_INDEX_CHOICE), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_BS_USER_NODE_GENERATION), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_BB_UPDT_BOUND_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_UBS_ALL_VARS), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_UBS_INT_VARS), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_BB_PARENT_SOL_STORING_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_PSSS_NO_STORING), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_PSSS_ONLY_PRIMAL), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_PSSS_PRIMAL_AND_DUAL), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_BB_CONSTRAINT_BRANCH_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_CBS_NO_CONSTRAINT_BRANCH), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_CBS_STBRANCH_PSEUDO_COSTS), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_CBS_LOWEST_NUMBER_OF_VARS), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_CBS_HIGHEST_NUMBER_OF_VARS), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_BB_BOUND_LIN_UPDT_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_BLUS_NO_UPDT), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_BLUS_ONE_VAR), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_BLUS_SUBGROUP), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_BLUS_ALL_VARS), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}



/*int muriqui::MRQ_strToEnum(const char *svalue, MRQ_BB_IGMA2_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_I2S_NO_IGMA2), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_I2S_UNTIL_FIRST_FEAS_SOL), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_I2S_ALWAYS), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}*/


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_BB_INT_HEURISTICS_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_IHS_NO_HEURISTICS), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_IHS_UNTIL_FIRST_FEAS_SOL), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_IHS_ALWAYS), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}



int muriqui::MRQ_strToEnum(const char *svalue, MRQ_BB_PARENT_NODE_BOUNDS_STORAGE_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_PNBSS_SCHAR), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_PNBSS_SHORT_INT), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_PNBSS_FLOAT), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_PNBSS_DOUBLE), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_BB_PSEUDO_PRUNING_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_PPS_NO_PSEUDO_PRUNING), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_PPS_ONLY_ON_NODE_EXPLORATION), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_BB_PPS_ON_NODE_EXPLORATION_AND_BRANCHING), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}



int muriqui::MRQ_strToEnum(const char *svalue, MRQ_DIVE_SELECT_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_DIVE_SS_FRACTIONAL), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_DIVE_SS_VECTORLENGHT), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}



int muriqui::MRQ_strToEnum(const char *svalue, MRQ_IGMA0_EPS_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_IGMA0_EPS_NO_UPDATE), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_IGMA0_EPS_BIN_SEARCH), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}



int muriqui::MRQ_strToEnum(const char *svalue, MRQ_IGMA_GAP_MIN_OBJ_STRATEGY &value)
{
    int ret = 0;
    
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_IGMA_GMOS_SAME_WEIGHT), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_IGMA_GMOS_BY_GAP_AVERAGE), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}




int muriqui::MRQ_strToEnum(const char *svalue, MRQ_ESHP_INTERIOR_POINT_STRATEGY &value)
{
    int ret = 0;
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_ESHP_IPS_MOST_INTERIOR), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_ESHP_IPS_CLOSED_TO_CONT_RELAX_SOL), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_SUBPROBLEM_NEIGHBORHOOD_STRATEGY &value)
{
    int ret = 0;
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_SNS_ORIGINAL), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_SNS_EUCLIDEAN_INTEGER_NEIGHBORHOOD), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_ROUNDING_STRATEGY &value)
{
    int ret = 0;
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_RS_NO_ROUNDING), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_RS_NEAREST_INTEGER), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_RS_PROBABILISTIC), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_RS_SSR), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}



int muriqui::MRQ_strToEnum(const char *svalue, MRQ_SSR_PREPROCESSING_STRATEGY &value)
{
    int ret = 0;
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_SSRPS_NO_PREPROCESSING), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_SSRPS_AFTER_EACH_CONSTRAINT), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_SSRPS_AFTER_EACH_CONSTRAINT_CLASS), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_SSR_VARIABLES_ADDITIONAL_FIXING_STRATEGY &value)
{
    int ret = 0;
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_SSR_VAFS_PREPROCESSING), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_SSR_VARIABLE_BOUNDS_UPDT_STRATEGY &value)
{
    int ret = 0;
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_SSR_VBUS_NO_UPDATING), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_SSR_VBUS_LINEAR_AUXILIAR_PROBLEM), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_SSR_VBUS_AUXILIAR_PROBLEM), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_SSR_CONTINUOUS_RELAXATION_STRATEGY_TO_STOCHASTIC_ROUNDING &value)
{
    int ret = 0;
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_SSR_CRSSR_NO_CONTINUOUS_RELAXATION), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_SSR_CRSSR_ONLY_BEFORE_STARTING), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_SSR_CRSSR_BEFORE_ALL_FIXINGS), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_strToEnum(const char *svalue, MRQ_IGMA2_NEIGHBORHOOD_STRATEGY &value)
{
    int ret = 0;
    
    if( MRQ_setEnum( MRQ_STRATT(MRQ_IGMA2_NS_RECTANGULAR), svalue, value ) == 0);
    else if( MRQ_setEnum( MRQ_STRATT(MRQ_IGMA2_NS_SPHERIC), svalue, value ) == 0);
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}



int muriqui::MRQ_enumToStr(const MRQ_ALG_CODE value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_UNDEFINED_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_OA_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LP_NLP_BB_OA_BASED_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_ECP_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_IGMA0_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_IGMA1_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_IGMA2_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LOC_BRANCH_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_FP_HEUR_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_DIVE_HEUR_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_OA_FP_HEUR_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_CONT_RELAX_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BONMIN_HYBRID_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_RENS_HEUR_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LP_BB_ECP_BASED_ALG), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LP_BB_ESH_BASED_ALG), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}
    
int muriqui::MRQ_enumToStr(const MRQ_MILP_SOLVER value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_UNDEFINED_MILP ), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_GLPK), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_CBC), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_CPLEX), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_GUROBI), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_XPRESS), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_MILP_MOSEK), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_MILP_KNITRO), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_NLP_SOLVER value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_UNDEFINED_NLP), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_NLP_MOSEK), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_NLP_KNITRO), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_IPOPT), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_OPTIZELLE), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_WORHP), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_LP_SOLVER value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_UNDEFINED_LP), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LP_GLPK), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LP_CPLEX), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LP_GUROBI), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LP_XPRESS), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LP_MOSEK), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LP_IPOPT), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LP_WORHP), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_LP_KNITRO), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_GLOBAL_SOLVER value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_UNDEFINED_GLOBAL), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_IQUAD), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_CONSTR_LIN_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_CLS_ALL_CONSTRS), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_CLS_ONLY_INFEAS_AND_ACTIVE), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_MASTER_SOL_ALSO), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_BY_BOX_FILL), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}

int muriqui::MRQ_enumToStr(const MRQ_OBJ_LIN_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_OLS_ALL_POINTS), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_OLS_NON_OBJ_CUT_POINTS), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}

int muriqui::MRQ_enumToStr(const MRQ_QUAD_APP_MASTER_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_QAMS_NO_QUAD_APP), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_QAMS_ON_BEST_POINT), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_QAMS_ON_LAST_POINT), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_BB_EXP_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_ES_DEPTH), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_ES_WIDTH), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_ES_BEST_LIMIT), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_ES_DEPTH_BEST_LIMIT), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_BB_BRANCH_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_BS_HIGHEST_INT_GAP), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_BS_BIN_FIRST_HIGHEST_INT_GAP), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_BS_STBRANCH_PSEUDO_COSTS), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_BS_VAR_PRIORITIES), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_BS_USER_INDEX_CHOICE), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_BS_USER_NODE_GENERATION), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_BB_UPDT_BOUND_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_UBS_ALL_VARS), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_UBS_INT_VARS), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_BB_PARENT_SOL_STORING_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_PSSS_NO_STORING), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_PSSS_ONLY_PRIMAL), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_PSSS_PRIMAL_AND_DUAL), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}



int muriqui::MRQ_enumToStr(const MRQ_BB_CONSTRAINT_BRANCH_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_CBS_NO_CONSTRAINT_BRANCH), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_CBS_STBRANCH_PSEUDO_COSTS), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_CBS_LOWEST_NUMBER_OF_VARS), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_CBS_HIGHEST_NUMBER_OF_VARS), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_BB_BOUND_LIN_UPDT_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_BLUS_NO_UPDT), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_BLUS_ONE_VAR), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_BLUS_SUBGROUP), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_BLUS_ALL_VARS), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


/*int muriqui::MRQ_enumToStr(const MRQ_BB_IGMA2_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_I2S_NO_IGMA2), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_I2S_UNTIL_FIRST_FEAS_SOL), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_I2S_ALWAYS), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
} */


int muriqui::MRQ_enumToStr(const MRQ_BB_INT_HEURISTICS_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_IHS_NO_HEURISTICS), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_IHS_UNTIL_FIRST_FEAS_SOL), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_IHS_ALWAYS), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_BB_PARENT_NODE_BOUNDS_STORAGE_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_PNBSS_SCHAR), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_PNBSS_SHORT_INT), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_PNBSS_FLOAT), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_PNBSS_DOUBLE), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_BB_PSEUDO_PRUNING_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_PPS_NO_PSEUDO_PRUNING), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_PPS_ONLY_ON_NODE_EXPLORATION), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_BB_PPS_ON_NODE_EXPLORATION_AND_BRANCHING), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_DIVE_SELECT_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_DIVE_SS_FRACTIONAL), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_DIVE_SS_VECTORLENGHT), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_IGMA0_EPS_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_IGMA0_EPS_NO_UPDATE), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_IGMA0_EPS_BIN_SEARCH), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_IGMA_GAP_MIN_OBJ_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_IGMA_GMOS_SAME_WEIGHT), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_IGMA_GMOS_BY_GAP_AVERAGE), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_ESHP_INTERIOR_POINT_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_ESHP_IPS_MOST_INTERIOR), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_ESHP_IPS_CLOSED_TO_CONT_RELAX_SOL), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_SUBPROBLEM_NEIGHBORHOOD_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SNS_ORIGINAL), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SNS_EUCLIDEAN_INTEGER_NEIGHBORHOOD), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}



int muriqui::MRQ_enumToStr(const MRQ_ROUNDING_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_RS_NO_ROUNDING), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_RS_NEAREST_INTEGER), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_RS_PROBABILISTIC), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_RS_SSR), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_SSR_PREPROCESSING_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SSRPS_NO_PREPROCESSING), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SSRPS_AFTER_EACH_CONSTRAINT), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SSRPS_AFTER_EACH_CONSTRAINT_CLASS), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_SSR_VARIABLES_ADDITIONAL_FIXING_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SSR_VAFS_PREPROCESSING), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_SSR_VARIABLE_BOUNDS_UPDT_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SSR_VBUS_NO_UPDATING), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SSR_VBUS_LINEAR_AUXILIAR_PROBLEM), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SSR_VBUS_AUXILIAR_PROBLEM), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_SSR_CONTINUOUS_RELAXATION_STRATEGY_TO_STOCHASTIC_ROUNDING value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SSR_CRSSR_NO_CONTINUOUS_RELAXATION), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SSR_CRSSR_ONLY_BEFORE_STARTING), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_SSR_CRSSR_BEFORE_ALL_FIXINGS), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


int muriqui::MRQ_enumToStr(const MRQ_IGMA2_NEIGHBORHOOD_STRATEGY value, char *svalue)
{
    int ret = 0;
    
    if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_IGMA2_NS_RECTANGULAR), value, svalue) == 0 );
    else if( MRQ_setStrByEnum( MRQ_ATTSTR(MRQ_IGMA2_NS_SPHERIC), value, svalue) == 0 );
    else
        ret = MRQ_VALUE_ERROR;
    
    return ret;
}


std::string muriqui::MRQ_getAlgName(MRQ_ALG_CODE code)
{
    switch(code)
    {
    case MRQ_OA_ALG:
        return "Outer Approximation";
    case MRQ_LP_NLP_BB_OA_BASED_ALG:
        return "LP/NLP based Branch-And-Bound based on Outer Approximation";
    case MRQ_ECP_ALG:
        return "Extended Cutting Plane";
    case MRQ_BB_ALG:
        return "Nonlinear Branch-And-Bound";
    case MRQ_IGMA0_ALG:
        return "Integrality Gap Minimiziation Algorithm - version 0";
    case MRQ_IGMA1_ALG:
        return "Integrality Gap Minimiziation Algorithm - version 1";
    case MRQ_IGMA2_ALG:
        return "Integrality Gap Minimiziation Algorithm - version 2";
    case MRQ_LOC_BRANCH_ALG:
        return "Local Branching";
    case MRQ_FP_HEUR_ALG:
        return "Feasibility Pump Heuristic";
    case MRQ_DIVE_HEUR_ALG:
        return "Diving Heuristic";
    case MRQ_OA_FP_HEUR_ALG:
        return "Feasibility Pump Heuristic based on Outer Approximation";
    case MRQ_ESH_ALG:
        return "Extended Supporting Hyperplane";
    case MRQ_CONT_RELAX_ALG:
        return "Continuous Relaxing";
    case MRQ_BONMIN_HYBRID_ALG:
        return "Bonmin Hybrid Branch-and-Bound/Outer Approximation";
    case MRQ_RENS_HEUR_ALG:
        return "Relaxation Enforced Neighborhood Search Heuristic";
    case MRQ_LP_BB_ECP_BASED_ALG:
        return "LP based Branch-And-Bound based on Extended Cutting Plane";
    case MRQ_LP_BB_ESH_BASED_ALG:
        return "LP based Branch-And-Bound based on Extended Supporting Hyperplane";
    case MRQ_SSR_HEUR_ALG:
        return "Structured Stochastic Rounding Heuristic";
    default:
        return "Undefined";
    }
}

std::string muriqui::MRQ_getAlgInitials(const MRQ_ALG_CODE code)
{
    switch(code)
    {
    case MRQ_OA_ALG:
        return "OA";
    case MRQ_LP_NLP_BB_OA_BASED_ALG:
        return "LP/NLP BB";
    case MRQ_ECP_ALG:
        return "ECP";
    case MRQ_BB_ALG:
        return "NL BB";
    case MRQ_IGMA0_ALG:
        return "IGMA0";
    case MRQ_IGMA1_ALG:
        return "IGMA1";
    case MRQ_IGMA2_ALG:
        return "IGMA2";
    case MRQ_LOC_BRANCH_ALG:
        return "LB";
    case MRQ_FP_HEUR_ALG:
        return "FP";
    case MRQ_DIVE_HEUR_ALG:
        return "Dive";
    case MRQ_OA_FP_HEUR_ALG:
        return "OAFP";
    case MRQ_ESH_ALG:
        return "ESH";
    case MRQ_CONT_RELAX_ALG:
        return "CRelax";
    case MRQ_BONMIN_HYBRID_ALG:
        return "BHBB";
    case MRQ_RENS_HEUR_ALG:
        return "RENS";
    case MRQ_LP_BB_ECP_BASED_ALG:
        return "LP BB";
    case MRQ_SSR_HEUR_ALG:
        return "SSR";
    default:
        return "Undefined";
    }
}




std::string muriqui::MRQ_getStatus(const int returnCode)
{
    switch(returnCode)
    {
    case MRQ_OPTIMAL_SOLUTION:
        return "optimal solution";
    case MRQ_INFEASIBLE_PROBLEM:
        return "infeasible problem";
    case MRQ_UNBOUNDED_PROBLEM:
        return "unbounded problem";
    case MRQ_ALG_NOT_APPLICABLE:
        return "algorithm not applicable to input problem";
    case MRQ_BAD_DEFINITIONS:
        return "inconsistent definitions";
    case MRQ_BAD_PARAMETER_VALUES:
        return "inconsistent value for parameter";
    case MRQ_INDEX_FAULT:
        return "wrong value for index";
    case MRQ_MEMORY_ERROR:
        return "memory error";
    case MRQ_UNDEFINED_ERROR:
        return "undefined error";
    case MRQ_MAX_TIME_STOP:
        return "maximum time achieved";
    case MRQ_MAX_ITERATIONS_STOP:
        return "maximum number of iterations achieved";
    case MRQ_CALLBACK_FUNCTION_ERROR:
        return "user callback function returned an error";
    case MRQ_STOP_REQUIRED_BY_USER:
        return "user callback function required stop";
    case MRQ_MILP_SOLVER_ERROR:
        return "some error was gotten from milp solver";
    case MRQ_NLP_SOLVER_ERROR:
        return "some error was gotten from nlp solver";
    case MRQ_LIBRARY_NOT_AVAILABLE:
        return "code was not compiled with a required library";
    case MRQ_NAME_ERROR:
        return "name error";
    case MRQ_VALUE_ERROR:
        return "value error";
    case MRQ_MINLP_SOLVER_ERROR:
        return "some error was gotten from minlp solver";
    case MRQ_INITIAL_SOLUTION_ERROR:
        return "initial solution error";
    case MRQ_HEURISTIC_FAIL:
        return "heuristic failed to find a solution";
    case MRQ_LACK_OF_PROGRESS_ERROR:
        return "lack of progress error";
    case MRQ_NLP_NO_FEASIBLE_SOLUTION:
        return "no nlp feasible solution";
    case MRQ_NONIMPLEMENTED_ERROR:
        return "non implemented procedure";
    case MRQ_HEURISTIC_SUCCESS:
        return "heuristic found a solution";
    case MRQ_CONT_RELAX_OPTIMAL_SOLUTION:
        return "subsolver optimal solution";
    case MRQ_NLP_EARLY_STOP:
        return "nlp early stop";
    case MRQ_NLP_FEASIBLE_SOLUTION:
        return "nlp solver got a feasible solution";
    case MRQ_MILP_FEASIBLE_SOLUTION:
        return "milp solver got a feasible solution";
    case MRQ_DEV_CONVERGENCE:
        return "metaheuristic got convergence";
    case MRQ_LAZY_MILP_BB_SOLUTION_LOWER_THAN_ZL:
        return "bb-lp/nlp: nlp objective value ";
    default:
        return "###############";
    }
}





































