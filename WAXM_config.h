#ifndef WAXM_CONFIG_H_
#define WAXM_CONFIG_H_

//flags and constants definitions to my projects

#define WAXM_HAVE_ASL                   1
#define WAXM_HAVE_GAMS                  1

#define WAXM_HAVE_CPLEX                 1
#define WAXM_HAVE_GUROBI                1
#define WAXM_HAVE_MOSEK                 1
#define WAXM_HAVE_GLPK                  0
#define WAXM_HAVE_CBC                   0
#define WAXM_HAVE_OSI                   0
#define WAXM_HAVE_XPRESS                0
#define WAXM_HAVE_IPOPT                 1
#define WAXM_HAVE_KNITRO                0


#define WAXM_HAVE_WORHP                 0
#define WAXM_HAVE_ALGENCAN              0
#define WAXM_HAVE_OPTIZELLE             0


#define WAXM_HAVE_BONMIN                0
#define WAXM_HAVE_IQUAD                 0
#define WAXM_HAVE_MURIQUI               0


#define WAXM_HAVE_LAPACK                1
#define WAXM_HAVE_CSDP                  0



#define WAXM_HAVE_CPP_2011		1
#define WAXM_CPP_MULTITHREADING 1


#define WAXM_HAVE_PUGIXML	1


#ifdef _MSC_VER
    //microsoft compiler
    #define WAXM_HAVE_POSIX		0
#else
    #define WAXM_HAVE_POSIX		1
#endif


#if WAXM_HAVE_POSIX
    #define WAXM_HAVE_CLOCK_GETTIME 1
#else
    #define WAXM_HAVE_CLOCK_GETTIME 0
#endif



#define WAXM_DEBUG_MODE		1

#define WAXM_SUPER_THREAD_DEBUG_MODE 	0

//to muriqui server
#define WAXM_HANDLE_TERMINATE_SIGNALS_ON_SERVER 1

//to optsolvers:
#define WAXM_OPT_DEBUG_MODE		1
#define WAXM_PRINT_MAX_ITER_WARNING 1
#define WAXM_PRINT_ERROR_RETURN_CODE_ON_SOLVE 1
#define WAXM_DEFAULT_MAX_NUMBER_OF_WARNINGS_BY_ITERATION_LIMIT 0

//to optsolver and minlpproblem
#define WAXM_PRINT_NLP_CALLBACK_FUNCTION_ERROR 0

#define WAXM_INFINITY 1e20

#define WAXM_CHARAC_COMENT_ON_PARAMS_FILE '#'

#endif
