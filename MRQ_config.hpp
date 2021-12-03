#ifndef MRQ_CONFIG_H_
#define MRQ_CONFIG_H_


#include "WAXM_config.h"

#define MRQ_HAVE_CLOCK_GETTIME WAXM_HAVE_CLOCK_GETTIME

#define MRQ_HAVE_ASL WAXM_HAVE_ASL
#define MRQ_HAVE_GAMS WAXM_HAVE_GAMS
#define MRQ_HAVE_CPP_2011 WAXM_HAVE_CPP_2011

#define MRQ_CPP_MULTITHREADING WAXM_CPP_MULTITHREADING
#define MRQ_OMP_MULTITHREADING WAXM_OMP_MULTITHREADING

#define MRQ_DEBUG_MODE WAXM_DEBUG_MODE

//#define MRQ_HAVE_CPLEX 1
//#define MRQ_HAVE_MOSEK 1
#define MRQ_HAVE_IQUAD WAXM_HAVE_IQUAD
//#define MRQ_HAVE_GUROBI 0


#define MRQ_SET_LBHEUR_ON_BB_NODE	1		//setting this as 0 will save 8 bytes per node. You can do it if you are reaching the memory bound in the branch and bound algorithm. However, if you are using hybrid BB + OA, this can damage the performance a little bit because this will change the exploration order in the best limit strategy, specially with pseudocosts. Prunes will occur in the same way because the lower bound from OA will stored in th lbbound.


#define MRQ_SAVE_OUTPUT_FILE 1

#define MRQ_CHARAC_COMENT_ON_PARAMS_FILE WAXM_CHARAC_COMENT_ON_PARAMS_FILE

#define MRQ_BB_SUPER_THREAD_DEBUG_MODE WAXM_SUPER_THREAD_DEBUG_MODE

#define MRQ_HANDLE_TERMINATE_SIGNALS_ON_SERVER 0 //WAXM_HANDLE_TERMINATE_SIGNALS_ON_SERVER

#endif
