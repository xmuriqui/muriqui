/** That file contains a interface to AMPL modeling system. To work, is necessary
* download and install the AMPL Solver Library (AMPL)
* 
* 
* 
* 
*/



#include <cstdio>
#include <cstdlib>

#include <new>
#include <vector>
#include <string>

#include "MRQ_dataStructures.hpp"
#include "MRQ_tools.hpp"
//#include "MRQ_functions.hpp"
#include "MRQ_algClasses.hpp"
#include "MRQ_ampl.hpp"



using namespace std;
using namespace muriqui;
using namespace minlpproblem;



#if !MRQ_READ_MODELING_SYSTEM_PARAMS

void inline MRQ_writeSolOnFile(FILE *file, double obj, const int n, const double *sol)
{
    fprintf( file, "%d\n%0.20f\n", obj, n );
    
    for(int i = 0; i < n; i++)
        fprintf( file, "%0.20f\n", sol[i] );
}


double inline MRQ_norm1Distance(const int n, const double *sol1, const double *sol2)
{
    double d = 0.0;
    
    for(int i = 0; i < n; i++)
    {
        double f = sol1[i] - sol2[i];
        d += MRQ_abs(f);
    }
    
    return d;
}


double inline MRQ_norm2Distance(const int n, const double *sol1, const double *sol2)
{
    double d = 0.0;
    
    for(int i = 0; i < n; i++)
    {
        double f = sol1[i] - sol2[i];
        d += f*f;
    }
    
    return sqrt(d);
}

double inline MRQ_integerDistance(const int nI, const int *intVars, const double *sol1, const double *sol2)
{
    double d = 0.0;
    
    for(int i = 0; i < nI; i++)
    {
        const int ind = intVars[i];
        d += MRQ_abs( sol1[ind] - sol2[ind] );
    }
    
    return d;
}

double inline MRQ_integer2Distance(const int nI, const int *intVars, const double *sol1, const double *sol2)
{
    double d = 0.0;
    
    for(int i = 0; i < nI; i++)
    {
        const int ind = intVars[i];
        double f = sol1[ind] - sol2[ind];
        d += f*f;
    }
    
    return sqrt(d);
}

#endif



void MRQ_ReadAmplModel::readAlgChoice(MRQ_AMPLParamData *paramData)
{
#if MRQ_HAVE_ASL
    
    char **argv = NULL;
    static keyword keywords[1] = {
    {(char *)"str", MRQ_readAMPLParams, paramData, (char *) "string  parameters at Muriqui"}
    };
    
    //only to set algorithm option (str)...
    static Option_Info opInfoAlg = {(char *)"muriqui", (char *) "muriqui", (char *)"muriqui_alg_choice", keywords, 1, 0, 0, 0, 0, 0, 0, 0 };
    
    argv = (char **) malloc( 2*sizeof(char*) );
    if(!argv)
    {
        MRQ_PRINTMEMERROR;
        goto termination;
    }
    
    argv[0] = (char *) malloc( sizeof(char) );
    if(!argv)
    {
        MRQ_PRINTMEMERROR;
        goto termination;
    }
    
    argv[1] = NULL;
    
    argv[0][0] = '\0';
    
    readParameters(argv, &opInfoAlg);
    
    
    
    
termination:
    
    if(argv)
    {
        if( argv[0] )	free(argv[0]);
        free(argv);
    }
    
#endif
}



void MRQ_ReadAmplModel::readMuriquiParameters( MRQ_AMPLParamData *paramData)
{
#if MRQ_HAVE_ASL
    
    char* argv[2];
    char argv0[] = "\0";
    char *argv1 = NULL;
    
    argv[0] = argv0;
    argv[1] = argv1;
    
    static keyword keywords[3] = {
    {(char *)"dbl", MRQ_readAMPLParams, paramData, (char *) "double/float (real) parameters at Muriqui"}, 
    {(char *)"int", MRQ_readAMPLParams, paramData, (char *) "integer parameters at Muriqui"}, 
    {(char *)"str", MRQ_readAMPLParams, paramData, (char *) "string  parameters at Muriqui"}
    };
    
    
    
    
    static Option_Info opInfo = {(char *)"muriqui", (char *) "muriqui", (char *)"muriqui_options", keywords, 3, 0, 0, 0, 0, 0, 0, 0 };
    
    static Option_Info opInfoGlobal = {(char *)"muriqui", (char *)"muriqui", (char *)"muriqui_global_options", keywords, 3, 0, 0, 0, 0, 0, 0, 0 }; //we do not use currently, but maybe in the future...
    
    static Option_Info opInfoMilp = {(char *)"muriqui", (char *)"muriqui", (char *)"muriqui_milp_options", keywords, 3, 0, 0, 0, 0, 0, 0, 0 };
    
    static Option_Info opInfoMinlp = {(char *)"muriqui", (char *)"muriqui", (char *)"muriqui_minlp_options", keywords, 3, 0, 0, 0, 0, 0, 0, 0 }; //we do not use currently, but maybe in the future...
    
    static Option_Info opInfoNlp = {(char *)"muriqui", (char *)"muriqui", (char *)"muriqui_nlp_options", keywords, 3, 0, 0, 0, 0, 0, 0, 0 };
    
    /*argv = (char **) malloc( 2*sizeof(char*) );
    MRQ_calloc(argv, 2);
    if(!argv)
        goto termination;
    
    argv[0] = (char *) malloc( sizeof(char) );
    if(!argv)
        goto termination;
        
    argv[1] = NULL;
    
    argv[0][0] = '\0';*/
    
    
    
    readParameters(argv, &opInfo);
    readParameters(argv, &opInfoMilp);
    readParameters(argv, &opInfoNlp);
    readParameters(argv, &opInfoMinlp);
    readParameters(argv, &opInfoGlobal);
    
    
    
//termination:
    
    /*if(argv)
    {
        if( argv[0] )	free(argv[0]);
        free(argv);
    }*/
    
    return;
#endif
}












int muriqui::MRQ_ampl(char *stub, const bool printAlgParameters, const bool printProblem)
#if MRQ_HAVE_ASL
{
    int code = 0, ret;
    char myMsg[200];
    double *psol, *dsol;
    
    MIP_NonLinearEval *myEval = NULL;
    MRQ_ReadAmplModel reader;
    
    
    FILE *outFile = NULL, *algFile = NULL;
    
    
    
    MRQ_MINLPProb prob;
    
    MRQ_GeneralSolverParams milpParams, nlpParams, globalParams, minlpParams;
    MRQ_Algorithm *alg = NULL;
    
    
    
    MRQ_AMPLParamData paramData;
    
    
    
    
    
    
    std::cout << MRQ_PREPRINT "Reading user model by AMPL interface:\n";
    ret = reader.readProblem(stub, prob, myEval);
    
    if(ret != 0)
    {
        #if MRQ_DEBUG_MODE
            //fprintf(stderr, "Error at reading ampl model! ret: %d\n", ret);
            MRQ_PRINTERRORNUMBER(ret);
        #endif
        
        goto termination;
        code = MRQ_UNDEFINED_ERROR;
    }
    
    if( reader.isMaximizationProblem() )
        std::cout << MRQ_PREPRINT "Maximization problem addressed as minimization\n";
    
    std::cout << MRQ_PREPRINT "Done\n";
    
    if( printProblem )
        prob.print();
    
    
    #if MRQ_READ_MODELING_SYSTEM_PARAMS
    {
        #if MRQ_SAVE_OUTPUT_FILE
            outFile = fopen( MRQ_OUT_FILE_NAME, "a" );

            if(outFile)
            {
                MRQ_writeProblemParameters(stub, prob, outFile);
            }
        #endif
        
        MRQ_GeneralSolverParams *sparam1, *sparam2;
        
        std::cout << MRQ_PREPRINT "Reading user parameters:\n";
        
        //we call gtops first time only to get the algorithm code 
        reader.readAlgChoice(&paramData);
        
        
        //we check if user define algorithm choice file
        {
            MRQ_ALG_CODE algCode = MRQ_UNDEFINED_ALG;
            int r = MRQ_tryReadAlgChoiceFile(MRQ_MURIQUI_ALG_CHOICE_FILE, 1, algCode);
            if( r == 0)
                paramData.algCode = algCode;
        }
        
        
        
        //if( paramData.algCode == MRQ_UNDEFINED_ALG )
            //paramData.algCode = MRQ_BB_ALG;
        
        
        alg = MRQ_newAlgorithm( paramData.algCode, prob.getNumberOfNLEqualityConstraints() );
        MRQ_IFMEMERRORGOTOLABEL(!alg, code, termination);
        
        
        alg->in_print_parameters_values = printAlgParameters;
        
        
        paramData.alg = alg;
        paramData.globalParams = &globalParams;
        paramData.milpParams = &milpParams;
        paramData.minlpParams = &minlpParams;
        paramData.nlpParams = &nlpParams;
        
        
        
        
        //now, we call gtops to get all remainder parameters
        reader.readMuriquiParameters(&paramData);
        
        
        if( alg->out_algorithm == MRQ_IGMA1_ALG )
        {
            sparam1 = &nlpParams;
            sparam2 = &globalParams;
        }
        else
        {
            sparam1 = &milpParams;
            sparam2 = &nlpParams;
        }
        
        
        std::cout << MRQ_PREPRINT "Done\n" ;
        
        {
            std::cout << MRQ_PREPRINT "Trying read input parameter file " MRQ_MURIQUI_PARAMS_FILE ".";
            
            int r = alg->readParametersWithTypeFromFile( MRQ_MURIQUI_PARAMS_FILE, true);
            if(r == 0)
                std::cout << " Done\n" ;
            else
                std::cout << " Failure\n" ;
        }
        
        
        {
            std::cout << MRQ_PREPRINT "Trying read milp solver input parameter file " MRQ_MILP_SOLVER_PARAMS_FILE  " .";
            
            int r = sparam1->storeParametersFromFile( MRQ_MILP_SOLVER_PARAMS_FILE, true);
            if(r == 0)
                std::cout << " Done\n" ;
            else
                std::cout << " Failure\n" ;
        }
        
        
        {
            std::cout << MRQ_PREPRINT "Trying read nlp solver input parameter file " MRQ_NLP_SOLVER_PARAMS_FILE  " .";
            
            int r = sparam2->storeParametersFromFile( MRQ_NLP_SOLVER_PARAMS_FILE, true);
            if(r == 0)
                    std::cout << " Done\n" ;
                else
                    std::cout << " Failure\n" ;
        }
        
        
        
        if( reader.isMaximizationProblem() )
        {
            //we have to reverse possible lower and upper bounds given by the user
            
            double aux = alg->in_upper_bound;
            
            if( alg->in_lower_bound > -MRQ_INFINITY )
                alg->in_upper_bound = -alg->in_lower_bound;
            
            if( aux < MRQ_INFINITY )
                alg->in_lower_bound = -aux;
        }
        
        
        
        alg->run(prob, sparam1, sparam2);
        
        printf("_________________________________________________________________________________\n");
        printf("Problem: %s Algorithm: %s (%d) \nReturn code: %d lower bound: %0.10f obj function: %0.10f \ntime: %0.2f cpu time: %0.2f iters: %ld\n", stub, alg->getAlgorithmName().c_str(), alg->out_algorithm,  alg->out_return_code, alg->out_lower_bound, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("_________________________________________________________________________________\n");
    
        #if MRQ_SAVE_OUTPUT_FILE
        if(outFile)
            fprintf(outFile, "%d%s%f%s%f%s%f%s%ld%s", alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
        
        if(outFile)	    fflush(outFile);
        #endif
        
        
        
        switch( alg->out_return_code )
        {
            case MRQ_OPTIMAL_SOLUTION:
                ret = 0;
                break;
            
            case MRQ_UNBOUNDED_PROBLEM:
                ret = 300;
                break;
                
            case MRQ_MAX_ITERATIONS_STOP:
                ret = 400;
                break;
            
            case MRQ_MAX_TIME_STOP:
                ret = 401;
                break;
            
            case MRQ_INFEASIBLE_PROBLEM:
                ret = 200;
                break;
            
            case MRQ_LIBRARY_NOT_AVAILABLE:
                ret = 501;
                break;
                
            case MRQ_CALLBACK_FUNCTION_ERROR:
                ret = 502;
                break;
            
            case MRQ_NLP_SOLVER_ERROR:
                ret = 505;
                break;
                
            default:
                ret = 500;
        }
        
        
        if( alg->out_best_obj < MRQ_INFINITY )
        {
            if( reader.isMaximizationProblem() )
                alg->out_best_obj = -alg->out_best_obj;
            
            
            if(alg->out_return_code == MRQ_OPTIMAL_SOLUTION)
                sprintf( myMsg, "An optimal solution was found! Obj function: %f. CPU Time: %f", alg->out_best_obj, alg->out_cpu_time );
            else
                sprintf( myMsg, "A feasible solution was found. Obj function: %f. CPU Time: %f", alg->out_best_obj, alg->out_cpu_time );
            
            psol= alg->out_best_sol;
        }
        else
        {
            sprintf( myMsg, "No feasible solution was found. CPU Time: %f", alg->out_cpu_time );
            
            psol = NULL;
        }
        
        dsol = NULL;
        reader.putInfoToAmpl(myMsg, ret, psol, dsol);
        
    }
    #else
    {
        int r = MRQ_testRun(prob, stub, printAlgParameters);
        if(r != 0)
            MRQ_PRINTERRORNUMBER(r);
    }
    #endif
    
    
    
    
    
    code = 0;
    
termination:
    
    #if !MRQ_READ_MODELING_SYSTEM_PARAMS
        alg = NULL;
    #endif
    
    #if MRQ_SAVE_OUTPUT_FILE
        if(outFile)
            fclose(outFile);
    #endif
    
    if( algFile )
        fclose(algFile);
    
    if(myEval) 	delete myEval;
    if(alg)		delete alg;
    
    
    return code;
}
#else
{
    return MRQ_LIBRARY_NOT_AVAILABLE;
}
#endif


























