
#ifndef OPT_CONFIG_HPP_
#define OPT_CONFIG_HPP_

#include "../WAXM_config.h"

#define OPT_HAVE_CPP_2011	WAXM_HAVE_CPP_2011

#define OPT_CPP_MULTITHREADING	 	WAXM_CPP_MULTITHREADING
#define OPT_OMP_MULTITHREADING		WAXM_OMP_MULTITHREADING


#define OPT_HAVE_CLOCK_GETTIME 	WAXM_HAVE_CLOCK_GETTIME

#define OPT_HAVE_CPLEX 		WAXM_HAVE_CPLEX
#define OPT_HAVE_GUROBI 	WAXM_HAVE_GUROBI
#define OPT_HAVE_MOSEK		WAXM_HAVE_MOSEK
#define OPT_HAVE_GLPK		WAXM_HAVE_GLPK
#define	OPT_HAVE_OSI		WAXM_HAVE_OSI
#define OPT_HAVE_CBC		WAXM_HAVE_CBC
#define OPT_HAVE_XPRESS		WAXM_HAVE_XPRESS
#define OPT_HAVE_IPOPT		WAXM_HAVE_IPOPT
#define OPT_HAVE_WORHP		WAXM_HAVE_WORHP
#define OPT_HAVE_KNITRO		WAXM_HAVE_KNITRO
#define OPT_HAVE_ALGENCAN	WAXM_HAVE_ALGENCAN
#define OPT_HAVE_OPTIZELLE  WAXM_HAVE_OPTIZELLE
#define OPT_HAVE_IQUAD		WAXM_HAVE_IQUAD

#define OPT_HAVE_CBC_OR_OSI		OPT_HAVE_CBC || OPT_HAVE_OSI

#define OPT_DEBUG_MODE  WAXM_OPT_DEBUG_MODE

#define OPT_SUPER_THREAD_DEBUG_MODE WAXM_SUPER_THREAD_DEBUG_MODE

#define OPT_CHARAC_COMENT_ON_PARAMS_FILE  WAXM_CHARAC_COMENT_ON_PARAMS_FILE

#define OPT_PRINT_NLP_CALLBACK_FUNCTION_ERROR WAXM_PRINT_NLP_CALLBACK_FUNCTION_ERROR

#define OPT_DEFAULT_MAX_NUMBER_OF_WARNINGS_BY_ITERATION_LIMIT   WAXM_DEFAULT_MAX_NUMBER_OF_WARNINGS_BY_ITERATION_LIMIT  

#endif