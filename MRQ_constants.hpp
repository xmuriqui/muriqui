/*
* MRQ_constants.hpp
*
*  Created on: 27/08/2013
*      Author: yo
*/

#ifndef MRQ_CONSTANTS_HPP_
#define MRQ_CONSTANTS_HPP_

#include <string>

#include "MRQ_config.hpp"
#include "OPT_solvers.hpp"
#include "BBL_constants.hpp"


#define MRQ_CASE_OF_SWITCH(str)   case str: return str


namespace muriqui
{
    #define MRQ_INFINITY MIP_INFINITY
    
    #define MRQ_MULTITHREADING 1
    
    #define MRQ_READ_MODELING_SYSTEM_PARAMS 1 //to set if we read parameters from systems like AMPL and GAMS
    
    
    
    
    # define MRQ_VERSION "0.7.04"
    # define MRQ_DATEV "24-Nov-2021"
    # define MRQ_AUTHOR "Wendel Melo"
    # define MRQ_AUTHOR_FILIATION "Computer Scientist, Federal University of Uberlandia, Brazil"
    # define MRQ_COLLABORATORS "Marcia Fampa (UFRJ, Brazil), Fernanda Raupp (LNCC, Brazil)"
    # define MRQ_EMAIL "wendelmelo@ufu.br"
    # define MRQ_BASE_PAPER "W Melo, M Fampa & F Raupp. An overview of MINLP algorithms and their \nimplementation in Muriqui Optimizer. Annals of Operations Research. 2018. \nDOI 10.1007/s10479-018-2872-5. https://rdcu.be/N5z0 ."
    
    //"W Melo, M Fampa & F Raupp Integrating nonlinear branch-and-bound and outer \napproximation for convex Mixed Integer Nonlinear Programming. Journal of Global \nOptimization, v. 60, p. 373-389, 2014."
    
    
    #define MRQ_MURIQUI_PARAMS_FILE "muriqui_params.opt"
    #define MRQ_MILP_SOLVER_PARAMS_FILE "muriqui_milp_params.opt"
    #define MRQ_NLP_SOLVER_PARAMS_FILE "muriqui_nlp_params.opt"
    
    
    #define MRQ_MURIQUI_ALG_CHOICE_FILE "muriqui_algorithm.opt"
    
    #define MRQ_CHAR_SEP "#"		//character separator to output file
    #define MRQ_OUT_FILE_NAME "muriqui_output.txt"
    
    
    #define MRQ_BINTOL 1.0e-1
    #define MRQ_CHECK_VAR_BOUNDS_ON_SOLVER_SOL_IN_BB 0


    enum MRQ_RETURN_CODE
    {
        MRQ_OPTIMAL_SOLUTION 		= 0,
        MRQ_INFEASIBLE_PROBLEM		= -1,
        MRQ_UNBOUNDED_PROBLEM		= -2,
        MRQ_ALG_NOT_APPLICABLE		= -3,
        MRQ_BAD_DEFINITIONS			= -4,
        MRQ_BAD_PARAMETER_VALUES	= -5,
        MRQ_INDEX_FAULT				= -6,
        MRQ_MEMORY_ERROR			= -7,
        MRQ_UNDEFINED_ERROR			= -8,
        MRQ_MAX_TIME_STOP			= -9,
        MRQ_MAX_ITERATIONS_STOP		= -10,
        MRQ_CALLBACK_FUNCTION_ERROR	= -11,
        
        MRQ_MILP_SOLVER_ERROR		= -12,
        MRQ_NLP_SOLVER_ERROR		= -13,
        MRQ_LIBRARY_NOT_AVAILABLE	= -14,
        MRQ_NAME_ERROR				= -15,
        MRQ_VALUE_ERROR				= -16,
        
        MRQ_MINLP_SOLVER_ERROR		= -17, //specially for local branching
        
        MRQ_INITIAL_SOLUTION_ERROR	= -18, //specially for local branching
        MRQ_HEURISTIC_FAIL			= -19, //specially for heuristics like Feasibility Pump and Diving:
        MRQ_LACK_OF_PROGRESS_ERROR	= -20,
        //specially for NLP solving
        MRQ_NLP_NO_FEASIBLE_SOLUTION= -21,
        MRQ_NONIMPLEMENTED_ERROR	= -22,
        
        
        //specially for heuristics like Feasibility Pump and Diving:
        MRQ_HEURISTIC_SUCCESS		= 1,
        MRQ_CONT_RELAX_OPTIMAL_SOLUTION		= 2,
        //specially for early branching:
        MRQ_NLP_EARLY_STOP			= 3,
        //specially for NLP solving
        MRQ_NLP_FEASIBLE_SOLUTION	= 4,
        //specially for MILP solving
        MRQ_MILP_FEASIBLE_SOLUTION	= 5,
        //specially for DEV:
        MRQ_DEV_CONVERGENCE			= 6,
        //internal code for nlp/nlp bb
        MRQ_LAZY_MILP_BB_SOLUTION_LOWER_THAN_ZL = 7,
        MRQ_STOP_REQUIRED_BY_USER	= 8,
        
        //code for generic success in any function. Note, is the same code of optimal solution
        MRQ_SUCCESS					= 0,
        
        
        MRQ_RETURN_CODE_BEGIN 		= MRQ_NONIMPLEMENTED_ERROR,
        MRQ_RETURN_CODE_END			= MRQ_STOP_REQUIRED_BY_USER
    };
    
    
    enum MRQ_ALG_CODE
    {
        MRQ_UNDEFINED_ALG 			    = 1000,
        MRQ_OA_ALG 					    = 1001,
        MRQ_LP_NLP_BB_OA_BASED_ALG 	    = 1002,
        MRQ_ECP_ALG  				    = 1003,
        MRQ_BB_ALG 					    = 1004,
        MRQ_IGMA0_ALG 				    = 1005,
        MRQ_IGMA1_ALG 				    = 1006,
        MRQ_IGMA2_ALG 				    = 1007,
        MRQ_LOC_BRANCH_ALG 			    = 1008,
        MRQ_FP_HEUR_ALG 			    = 1009,
        MRQ_DIVE_HEUR_ALG 			    = 1010,
        MRQ_OA_FP_HEUR_ALG 			    = 1011,
        MRQ_ESH_ALG					    = 1012,
        MRQ_CONT_RELAX_ALG 			    = 1013,
        MRQ_BONMIN_HYBRID_ALG		    = 1014,
        MRQ_RENS_HEUR_ALG			    = 1015,
        MRQ_LP_BB_ECP_BASED_ALG 	    = 1016,
        MRQ_LP_BB_ESH_BASED_ALG 	    = 1017,
        MRQ_SSR_HEUR_ALG = 1018,
        
        MRQ_ALG_CODE_BEGIN = MRQ_UNDEFINED_ALG,
        MRQ_ALG_CODE_END   = MRQ_RENS_HEUR_ALG
    };
    
    
    enum MRQ_LP_SOLVER
    {
        MRQ_UNDEFINED_LP 	= optsolvers::OPT_UNDEFINEDSOLVER,
        MRQ_LP_GLPK 		= optsolvers::OPT_LP_GLPK,
        MRQ_LP_CBC			= optsolvers::OPT_LP_CBC,
        MRQ_LP_CPLEX 		= optsolvers::OPT_LP_CPLEX,
        MRQ_LP_GUROBI 		= optsolvers::OPT_LP_GUROBI,
        MRQ_LP_XPRESS 		= optsolvers::OPT_LP_XPRESS,
        MRQ_LP_MOSEK 		= optsolvers::OPT_LP_MOSEK,
        MRQ_LP_IPOPT 		= optsolvers::OPT_LP_IPOPT,
        MRQ_LP_KNITRO 		= optsolvers::OPT_LP_KNITRO,
        MRQ_LP_OPTIZELLE    = optsolvers::OPT_LP_OPTIZELLE,
        MRQ_LP_WORHP 		= optsolvers::OPT_LP_WORHP,
        MRQ_LP_ALGENCAN 	= optsolvers::OPT_LP_ALGENCAN
    };
    
    
    enum MRQ_MILP_SOLVER
    {
        MRQ_UNDEFINED_MILP 		= optsolvers::OPT_UNDEFINEDSOLVER,
        MRQ_GLPK				= optsolvers::OPT_LP_GLPK,
        MRQ_CBC					= optsolvers::OPT_LP_CBC,
        MRQ_CPLEX				= optsolvers::OPT_LP_CPLEX,
        MRQ_GUROBI 				= optsolvers::OPT_LP_GUROBI,
        MRQ_XPRESS				= optsolvers::OPT_LP_XPRESS,
        MRQ_MILP_MOSEK 			= optsolvers::OPT_LP_MOSEK,
        MRQ_MILP_KNITRO			= optsolvers::OPT_LP_KNITRO,
        
        
        MRQ_MILP_SOLVER_BEGIN 	= MRQ_UNDEFINED_MILP,
        MRQ_MILP_SOLVER_END 	= MRQ_MILP_KNITRO
    };
    
    //nonlinear solvers codes
    enum MRQ_NLP_SOLVER
    {
        MRQ_UNDEFINED_NLP 		= optsolvers::OPT_UNDEFINEDSOLVER,
        MRQ_NLP_MOSEK 			= optsolvers::OPT_NLP_MOSEK,
        MRQ_NLP_KNITRO			= optsolvers::OPT_NLP_KNITRO,
        //MRQ_OLDIPOPT 			= optsolvers::OPT_NLP_OLDIPOPT,
        MRQ_IPOPT 				= optsolvers::OPT_NLP_IPOPT,
        MRQ_OPTIZELLE           = optsolvers::OPT_OPTIZELLE,
        MRQ_WORHP				= optsolvers::OPT_NLP_WORHP,
        MRQ_ALGENCAN 			= optsolvers::OPT_NLP_ALGENCAN,
        
        
        MRQ_NLP_SOLVER_BEGIN 	= MRQ_UNDEFINED_NLP,
        MRQ_NLP_SOLVER_END 		= MRQ_WORHP
    };
    
    
    enum MRQ_GLOBAL_SOLVER
    {
        MRQ_UNDEFINED_GLOBAL 	= optsolvers::OPT_UNDEFINEDSOLVER,
        MRQ_IQUAD				= optsolvers::OPT_NLP_IQUAD,
        
        MRQ_GLOBAL_SOLVER_BEGIN = MRQ_UNDEFINED_GLOBAL,
        MRQ_GLOBAL_SOLVER_END	= MRQ_IQUAD
    };
    
    
    //strategies for constraints linearizations in each linearization point, designed to work on linear approximation algorithms
    enum MRQ_CONSTR_LIN_STRATEGY
    {
        MRQ_CLS_ALL_CONSTRS								= 310,
        MRQ_CLS_ONLY_INFEAS_AND_ACTIVE,
        //experimental: only for Outer Approximation
        MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_MASTER_SOL_ALSO,
        //experimental: only for Outer Approximation
        MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_BY_BOX_FILL,
        
        MRQ_CONSTR_LIN_STRATEGY_BEGIN = MRQ_CLS_ALL_CONSTRS,
        MRQ_CONSTR_LIN_STRATEGY_END =	MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_BY_BOX_FILL
    };
    
    
    //strategies to liniearize objective function in the set of linearization points, designed to work on linear approximation algorithms
    enum MRQ_OBJ_LIN_STRATEGY
    {
        MRQ_OLS_ALL_POINTS			= 320,
        MRQ_OLS_NON_OBJ_CUT_POINTS	= 321,
        
        MRQ_OBJ_LIN_STRATEGY_BEGIN	= MRQ_OLS_ALL_POINTS,
        MRQ_OBJ_LIN_STRATEGY_END	= MRQ_OLS_NON_OBJ_CUT_POINTS
    };
    
    
    //strategies for quadratic approximation to master problem in linear approximation algorithm
    enum MRQ_QUAD_APP_MASTER_STRATEGY
    {
        MRQ_QAMS_NO_QUAD_APP		= 350,
        MRQ_QAMS_ON_BEST_POINT,
        MRQ_QAMS_ON_LAST_POINT,
        
        MRQ_QUAD_APP_MASTER_STRATEGY_BEGIN = MRQ_QAMS_NO_QUAD_APP,
        MRQ_QUAD_APP_MASTER_STRATEGY_END = MRQ_QAMS_ON_LAST_POINT
    };
    
    
    enum MRQ_BB_EXP_STRATEGY
    {
        MRQ_BB_ES_DEPTH			= 		701,
        MRQ_BB_ES_WIDTH			= 		702,
        MRQ_BB_ES_BEST_LIMIT	= 		703,
        MRQ_BB_ES_DEPTH_BEST_LIMIT = 	704,
        
        MRQ_BB_EXP_STRATEGY_BEGIN = 	MRQ_BB_ES_DEPTH,
        MRQ_BB_EXP_STRATEGY_END = 		MRQ_BB_ES_DEPTH_BEST_LIMIT
    };
    
    
    enum MRQ_BB_BRANCH_STRATEGY
    {
        MRQ_BB_BS_HIGHEST_INT_GAP = 			801,
        MRQ_BB_BS_BIN_FIRST_HIGHEST_INT_GAP = 	802,
        MRQ_BB_BS_STBRANCH_PSEUDO_COSTS =		803,
        MRQ_BB_BS_VAR_PRIORITIES = 804,
        MRQ_BB_BS_USER_INDEX_CHOICE = 			805,
        MRQ_BB_BS_USER_NODE_GENERATION = 		806,
        
        MRQ_BB_BRANCH_STRATEGY_BEGIN = 			MRQ_BB_BS_HIGHEST_INT_GAP,
        MRQ_BB_BRANCH_STRATEGY_END = 			MRQ_BB_BS_USER_NODE_GENERATION
    };
    
    
    enum MRQ_BB_UPDT_BOUND_STRATEGY
    {
        MRQ_BB_UBS_ALL_VARS = 810,
        MRQ_BB_UBS_INT_VARS,
        
        MRQ_BB_UPDT_BOUND_STRATEGY_BEGIN 	= MRQ_BB_UBS_ALL_VARS,
        MRQ_BB_UPDT_BOUND_STRATEGY_END 		= MRQ_BB_UBS_INT_VARS
    };
    
    
    
    
    
    //strategy to warm start on branch-and-bound (and parent solution storing in each node). If you use pseudo cost as branching strategy, you need at lest primal strategy to warm start...
    enum MRQ_BB_PARENT_SOL_STORING_STRATEGY
    {
        MRQ_BB_PSSS_NO_STORING = 901,
        MRQ_BB_PSSS_ONLY_PRIMAL,
        MRQ_BB_PSSS_PRIMAL_AND_DUAL,
        
        MRQ_BB_PARENT_SOL_STORING_STRATEGY_BEGIN = MRQ_BB_PSSS_NO_STORING,
        MRQ_BB_PARENT_SOL_STORING_STRATEGY_END = MRQ_BB_PSSS_PRIMAL_AND_DUAL
    };
    
    
    
    //startegy to choose an binary sum equality constraint to perform branching...
    enum MRQ_BB_CONSTRAINT_BRANCH_STRATEGY
    {
        MRQ_BB_CBS_NO_CONSTRAINT_BRANCH = 910,
        MRQ_BB_CBS_STBRANCH_PSEUDO_COSTS,
        MRQ_BB_CBS_LOWEST_NUMBER_OF_VARS,
        MRQ_BB_CBS_HIGHEST_NUMBER_OF_VARS,
        
        
        MRQ_BB_CONSTRAINT_BRANCH_STRATEGY_BEGIN = MRQ_BB_CBS_NO_CONSTRAINT_BRANCH,
        MRQ_BB_CONSTRAINT_BRANCH_STRATEGY_END = MRQ_BB_CBS_HIGHEST_NUMBER_OF_VARS
    };
    
    
    enum MRQ_BB_BOUND_LIN_UPDT_STRATEGY
    {
        MRQ_BB_BLUS_NO_UPDT		= 720,
        MRQ_BB_BLUS_ONE_VAR,
        MRQ_BB_BLUS_SUBGROUP,
        MRQ_BB_BLUS_ALL_VARS,
        
        MRQ_BB_BOUND_LIN_UPDT_STRATEGY_BEGIN	= MRQ_BB_BLUS_NO_UPDT,
        MRQ_BB_BOUND_LIN_UPDT_STRATEGY_END		= MRQ_BB_BLUS_ALL_VARS
    };
    
    
    /*enum MRQ_BB_IGMA2_STRATEGY
    {
        MRQ_BB_I2S_NO_IGMA2				= 730,
        MRQ_BB_I2S_UNTIL_FIRST_FEAS_SOL,
        MRQ_BB_I2S_ALWAYS,
        
        MRQ_BB_IGMA2_STRATEGY_BEGIN		= MRQ_BB_I2S_NO_IGMA2,
        MRQ_BB_IGMA2_STRATEGY_END		= MRQ_BB_I2S_ALWAYS
    }; */
    
    
    enum MRQ_BB_INT_HEURISTICS_STRATEGY
    {
        MRQ_BB_IHS_NO_HEURISTICS		= 740,
        MRQ_BB_IHS_UNTIL_FIRST_FEAS_SOL,
        MRQ_BB_IHS_ALWAYS,
        
        MRQ_BB_INT_HEURISTICS_STRATEGY_BEGIN	= MRQ_BB_IHS_NO_HEURISTICS,
        MRQ_BB_INT_HEURISTICS_STRATEGY_END 		= MRQ_BB_IHS_ALWAYS
    };
    
    
    enum MRQ_BB_PARENT_NODE_BOUNDS_STORAGE_STRATEGY
    {
        MRQ_BB_PNBSS_SCHAR 		= branchAndBound::BBL_PNBSS_SCHAR,
        MRQ_BB_PNBSS_SHORT_INT 	= branchAndBound::BBL_PNBSS_SHORT_INT,
        MRQ_BB_PNBSS_FLOAT 		= branchAndBound::BBL_PNBSS_FLOAT,
        MRQ_BB_PNBSS_DOUBLE 	= branchAndBound::BBL_PNBSS_DOUBLE,
    };
    
    
    enum MRQ_BB_PSEUDO_PRUNING_STRATEGY
    {
        MRQ_BB_PPS_NO_PSEUDO_PRUNING	= 780,
        MRQ_BB_PPS_ONLY_ON_NODE_EXPLORATION,
        MRQ_BB_PPS_ON_NODE_EXPLORATION_AND_BRANCHING,
    };
    
    
    //diving heuristic strategies for select a fractional variable to fix
    enum MRQ_DIVE_SELECT_STRATEGY
    {
        MRQ_DIVE_SS_FRACTIONAL =	820,
        MRQ_DIVE_SS_VECTORLENGHT =	821,
        
        MRQ_DIVE_SELECT_STRATEGY_BEGIN		= MRQ_DIVE_SS_FRACTIONAL,
        MRQ_DIVE_SELECT_STRATEGY_END		= MRQ_DIVE_SS_VECTORLENGHT
    };
    
    
    //strategies to update epsilon parameter to objective cut og IGMA1
    enum MRQ_IGMA0_EPS_STRATEGY
    {
        MRQ_IGMA0_EPS_NO_UPDATE 	= 830,
        MRQ_IGMA0_EPS_BIN_SEARCH	= 831,
        
        MRQ_IGMA0_EPS_STRATEGY_BEGIN	= MRQ_IGMA0_EPS_NO_UPDATE,
        MRQ_IGMA0_EPS_STRATEGY_END		= MRQ_IGMA0_EPS_BIN_SEARCH
    };
    
    
    enum MRQ_IGMA_GAP_MIN_OBJ_STRATEGY
    {
        MRQ_IGMA_GMOS_SAME_WEIGHT		= 840,
        MRQ_IGMA_GMOS_BY_GAP_AVERAGE	= 841,
        
        MRQ_IGMA_GAP_MIN_OBJ_STRATEGY_BEGIN = MRQ_IGMA_GMOS_SAME_WEIGHT,
        MRQ_IGMA_GAP_MIN_OBJ_STRATEGY_END = MRQ_IGMA_GMOS_BY_GAP_AVERAGE
    };
    
    
    enum MRQ_ESHP_INTERIOR_POINT_STRATEGY
    {
        MRQ_ESHP_IPS_MOST_INTERIOR				= 850,
        MRQ_ESHP_IPS_CLOSED_TO_CONT_RELAX_SOL	= 851,
        
        MRQ_ESHP_INTERIOR_POINT_STRATEGY_BEGIN	= MRQ_ESHP_IPS_MOST_INTERIOR,
        MRQ_ESHP_INTERIOR_POINT_STRATEGY_END	= MRQ_ESHP_IPS_CLOSED_TO_CONT_RELAX_SOL
    };
    
    
    //codes to specify target of general parameters in distributed computation 
    enum MRQ_DC_GENERAL_PARAMS_TARGET
    {
        MRQ_DC_GPT_MURIQUI		= 870,
        MRQ_DC_GPT_MILP_SOLVER,
        MRQ_DC_GPT_NLP_SOLVER,
        MRQ_DC_GPT_GLOBAL_SOLVER,
        
        MRQ_DC_GENERAL_PARAMS_TARGET_BEGIN	= MRQ_DC_GPT_MURIQUI,
        MRQ_DC_GENERAL_PARAMS_TARGET_END	= MRQ_DC_GPT_GLOBAL_SOLVER,
    };
    
    
    enum MRQ_SUBPROBLEM_NEIGHBORHOOD_STRATEGY
    {
        MRQ_SNS_ORIGINAL		= 890, //find the optimal rounding
        MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD,
        MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD,
        MRQ_SNS_EUCLIDEAN_INTEGER_NEIGHBORHOOD,
        
        MRQ_SUBPROBLEM_NEIGHBORHOOD_STRATEGY_BEGIN = MRQ_SNS_ORIGINAL,
        MRQ_SUBPROBLEM_NEIGHBORHOOD_STRATEGY_END = MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD
    };
    
    enum MRQ_ROUNDING_STRATEGY
    {
        MRQ_RS_NO_ROUNDING			=	950,
        MRQ_RS_NEAREST_INTEGER,
        MRQ_RS_PROBABILISTIC,
        MRQ_RS_SSR,
        
        MRQ_ROUNDING_STRATEGY_BEGIN	= MRQ_RS_NO_ROUNDING,
        MRQ_ROUNDING_STRATEGY_END	= MRQ_RS_SSR,
    };
    
    
    enum MRQ_SSR_PREPROCESSING_STRATEGY
    {
        MRQ_SSRPS_NO_PREPROCESSING = 960,
        MRQ_SSRPS_AFTER_EACH_CONSTRAINT,
        MRQ_SSRPS_AFTER_EACH_CONSTRAINT_CLASS,
        
        MRQ_SSR_PREPROCESSING_STRATEGY_BEGIN = MRQ_SSRPS_NO_PREPROCESSING,
        MRQ_SSR_PREPROCESSING_STRATEGY_END = MRQ_SSRPS_AFTER_EACH_CONSTRAINT_CLASS
    };
    
    
    enum MRQ_SSR_VARIABLES_ADDITIONAL_FIXING_STRATEGY
    {
        MRQ_SSR_VAFS_PREPROCESSING = 965,
        //MRQ_SSR_VAFS_LINEAR_AUXILIAR_PROBLEM,
        //MRQ_SSR_VAFS_AUXILIAR_PROBLEM,
        
        MRQ_SSR_VARIABLE_ADDITIONAL_FIXING_STRATEGY_BEGIN = MRQ_SSR_VAFS_PREPROCESSING,
        MRQ_SSR_VARIABLE_ADDITIONAL_FIXING_STRATEGY_END = MRQ_SSR_VAFS_PREPROCESSING
    };
    
    //just before each variable stochastic rounding inside SSR
    enum MRQ_SSR_VARIABLE_BOUNDS_UPDT_STRATEGY
    {
        MRQ_SSR_VBUS_NO_UPDATING = 970,
        MRQ_SSR_VBUS_LINEAR_AUXILIAR_PROBLEM,
        MRQ_SSR_VBUS_AUXILIAR_PROBLEM,
        
        MRQ_SSR_VARIABLE_BOUNDS_UPDT_STRATEGY_BEGIN = MRQ_SSR_VBUS_NO_UPDATING,
        MRQ_SSR_VARIABLE_BOUNDS_UPDT_STRATEGY_END = MRQ_SSR_VBUS_AUXILIAR_PROBLEM
    };
    
    
    enum MRQ_SSR_CONTINUOUS_RELAXATION_STRATEGY_TO_STOCHASTIC_ROUNDING
    {
        MRQ_SSR_CRSSR_NO_CONTINUOUS_RELAXATION = 980,
        MRQ_SSR_CRSSR_ONLY_BEFORE_STARTING,
        MRQ_SSR_CRSSR_BEFORE_ALL_FIXINGS
    };
    
    
    enum MRQ_IGMA2_NEIGHBORHOOD_STRATEGY
    {
        MRQ_IGMA2_NS_RECTANGULAR	= 990,
        MRQ_IGMA2_NS_SPHERIC,
        
        MRQ_IGMA2_NEIGHBORHOOD_STRATEGY_BEGIN = MRQ_IGMA2_NS_RECTANGULAR,
        MRQ_IGMA2_NEIGHBORHOOD_STRATEGY_END   = MRQ_IGMA2_NS_SPHERIC,
    };
    
    
    //to muriqui client and server
    enum MRQ_INPUT_FILE_TYPE
    {
        MRQ_IFT_AMPL_MODEL_FILE = 9000,
        MRQ_IFT_GAMS_MODEL_FILE,
    };
    
    
    
    
    MRQ_RETURN_CODE MRQ_intToReturnCode( int returnCode);
    
    MRQ_NLP_SOLVER MRQ_intToNLPSolver( int nlpSolver );
    
    MRQ_MILP_SOLVER MRQ_intToMILPSolver( int milpSolver );
    
    
    
    
    
    
    
    int MRQ_strToEnum(const char *svalue, MRQ_ALG_CODE &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_MILP_SOLVER &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_LP_SOLVER &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_NLP_SOLVER &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_GLOBAL_SOLVER &value);
    
    
    int MRQ_strToEnum(const char *svalue, MRQ_CONSTR_LIN_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_OBJ_LIN_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_QUAD_APP_MASTER_STRATEGY &value);
    
    
    int MRQ_strToEnum(const char *svalue, MRQ_BB_EXP_STRATEGY &value);
    
    
    int MRQ_strToEnum(const char *svalue, MRQ_BB_BRANCH_STRATEGY &value);
    
    
    int MRQ_strToEnum(const char *svalue, MRQ_BB_UPDT_BOUND_STRATEGY &value);
    
    
    int MRQ_strToEnum(const char *svalue, MRQ_BB_PARENT_SOL_STORING_STRATEGY &value);
    
    
    int MRQ_strToEnum(const char *svalue, MRQ_BB_CONSTRAINT_BRANCH_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_BB_BOUND_LIN_UPDT_STRATEGY &value);
    
    //int MRQ_strToEnum(const char *svalue, MRQ_BB_IGMA2_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_BB_INT_HEURISTICS_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_BB_PARENT_NODE_BOUNDS_STORAGE_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_BB_PSEUDO_PRUNING_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_DIVE_SELECT_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_IGMA0_EPS_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_IGMA_GAP_MIN_OBJ_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_ESHP_INTERIOR_POINT_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_SUBPROBLEM_NEIGHBORHOOD_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_ROUNDING_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_SSR_PREPROCESSING_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_SSR_VARIABLES_ADDITIONAL_FIXING_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_SSR_VARIABLE_BOUNDS_UPDT_STRATEGY &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_SSR_CONTINUOUS_RELAXATION_STRATEGY_TO_STOCHASTIC_ROUNDING &value);
    
    int MRQ_strToEnum(const char *svalue, MRQ_IGMA2_NEIGHBORHOOD_STRATEGY &value);
    
    
    int MRQ_enumToStr(const MRQ_ALG_CODE value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_MILP_SOLVER value, char *svalue);
    
    
    int MRQ_enumToStr(const MRQ_NLP_SOLVER value, char *svalue);
    
    
    int MRQ_enumToStr(const MRQ_LP_SOLVER value, char *svalue);
    
    
    int MRQ_enumToStr(const MRQ_GLOBAL_SOLVER value, char *svalue);
    
    
    int MRQ_enumToStr(const MRQ_CONSTR_LIN_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_OBJ_LIN_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_QUAD_APP_MASTER_STRATEGY value, char *svalue);
    
    
    int MRQ_enumToStr(const MRQ_BB_EXP_STRATEGY value, char *svalue);
    
    
    int MRQ_enumToStr(const MRQ_BB_BRANCH_STRATEGY value, char *svalue);
    
    
    int MRQ_enumToStr(const MRQ_BB_UPDT_BOUND_STRATEGY value, char *svalue);
    
    
    int MRQ_enumToStr(const MRQ_BB_PARENT_SOL_STORING_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_BB_CONSTRAINT_BRANCH_STRATEGY value, char *svalue);
    
    
    int MRQ_enumToStr(const MRQ_BB_BOUND_LIN_UPDT_STRATEGY value, char *svalue);
    
    
    //int MRQ_enumToStr(const MRQ_BB_IGMA2_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_BB_INT_HEURISTICS_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_BB_PARENT_NODE_BOUNDS_STORAGE_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_BB_PSEUDO_PRUNING_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_DIVE_SELECT_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_IGMA0_EPS_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_IGMA_GAP_MIN_OBJ_STRATEGY value, char *svalue);
    
    
    int MRQ_enumToStr(const MRQ_ESHP_INTERIOR_POINT_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_SUBPROBLEM_NEIGHBORHOOD_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_ROUNDING_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_SSR_PREPROCESSING_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_SSR_VARIABLES_ADDITIONAL_FIXING_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_SSR_VARIABLE_BOUNDS_UPDT_STRATEGY value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_SSR_CONTINUOUS_RELAXATION_STRATEGY_TO_STOCHASTIC_ROUNDING value, char *svalue);
    
    int MRQ_enumToStr(const MRQ_IGMA2_NEIGHBORHOOD_STRATEGY value, char *svalue);
    
    
    std::string MRQ_getAlgName(MRQ_ALG_CODE code);
    
    std::string MRQ_getAlgInitials(MRQ_ALG_CODE code);
    
    //this function receives a MRQ_RETURN_CODE and return a string with the status
    std::string MRQ_getStatus(int returnCode);
}






#endif /* MRQ_CONSTANTS_HPP_ */




