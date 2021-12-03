
#ifndef MIP_CONFIG_H_
#define MIP_CONFIG_H_

#include "../WAXM_config.h"

#define MIP_HAVE_ASL WAXM_HAVE_ASL
#define MIP_HAVE_GAMS WAXM_HAVE_GAMS


#define MIP_CPP_MULTITHREADING	 	WAXM_CPP_MULTITHREADING
#define MIP_OMP_MULTITHREADING		WAXM_OMP_MULTITHREADING


//this flag is to NLP callbacks calling from MIP_MINLPProb
#define MIP_PRINT_NLP_CALLBACK_FUNCTION_ERROR WAXM_PRINT_NLP_CALLBACK_FUNCTION_ERROR

//this flag is to NLP callbacks from AMPL. Note we decide having a flag only to our NLP evaluation objects, since the callback errors was being printed twince (in NLPEval and MIP_MINLPProb)
#define MIP_PRINT_NLPEVALOBJECT_EVAL_ERROR 1


#define MIP_INFINITY  WAXM_INFINITY
    
#define MIP_DEBUG_MODE WAXM_DEBUG_MODE


//maximum number of evaluation error messages printed
#define MIP_MAX_PRINTS_NLPOBJECTEVAL_EVAL_ERRORS 100


#define MIP_AUTHOR	"Wendel Melo"
# define MIP_AUTHOR_FILIATION "Computer Scientist, Federal University of Uberlandia, Brazil"

#define MIP_PRINT_PREPROC_CUT_MSG 0



#endif
