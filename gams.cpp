
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <string>

#include "MIP_gams.hpp"

#include "muriqui.hpp"
#include "MRQ_tools.hpp"



using namespace muriqui;





int MRQ_gams(char *stub, const bool printAlgParameters, const bool printProblem)
{
    int r, retCode;
    
    FILE *outFile = NULL;
    
    MRQ_MINLPProb prob;
    MRQ_GeneralSolverParams milpParams, nlpParams;
    
    minlpproblem::MIP_NonLinearEval *eval = NULL;
    minlpproblem::MIP_GamsModelReader reader;
    
    
    
    
    
    std::cout << MRQ_PREPRINT "Reading user model by GAMS interface:\n";
    
    r =  reader.readProblem(stub, prob, eval);
    MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_UNDEFINED_ERROR, termination);
    
    
    
    if( reader.isMaximizationProblem() )
        std::cout << MRQ_PREPRINT "Maximization problem addressed as minimization\n";
    
    
    std::cout << MRQ_PREPRINT "Done\n";
    
    
    if( printProblem )
        prob.print();
    
    #if MRQ_READ_MODELING_SYSTEM_PARAMS
    {
        MRQ_Algorithm *alg = NULL;
        
        int gamsModelReturnCode = 
        #if MRQ_HAVE_GAMS
            gmoModelStat_ErrorUnknown;
        #else
            -1;
        #endif
        int gamsSolverReturnCode =
        #if MRQ_HAVE_GAMS
            gmoSolveStat_InternalErr;
        #else
            -1;
        #endif
            
        
        MRQ_ALG_CODE algCode = MRQ_UNDEFINED_ALG;
        
        
        #if MRQ_SAVE_OUTPUT_FILE
            outFile = fopen( MRQ_OUT_FILE_NAME, "a" );

            if(outFile)
            {
                MRQ_writeProblemParameters(stub, prob, outFile);
            }
        #endif
        
        
        
        
        std::cout << MRQ_PREPRINT "Reading user parameters:\n";
        
        std::cout << MRQ_PREPRINT "Trying reading algorithm choice file " MRQ_MURIQUI_ALG_CHOICE_FILE " . ";
        
        
        MRQ_tryReadAlgChoiceFile(MRQ_MURIQUI_ALG_CHOICE_FILE, 1, algCode);
        
        
        alg = MRQ_newAlgorithm(algCode, prob.getNumberOfNLEqualityConstraints() );
        MRQ_IFMEMERRORGOTOLABEL(!alg, retCode, termination);
        
        
        alg->in_print_parameters_values = printAlgParameters;
        
        
        //reading gams options
        {
            int ivalue;
            double dvalue;
            /*
            * List of gams options that we have interest:
            * 
            * iterLim:	Iteration limit of solver
            * optCA:	Absolute Optimality criterion solver default
            * optCR:	Relative Optimality criterion solver default
            * resLim:	Wall-clock time limit for solver
            * threads	Number of threads to be used by a solver
            */ 
            
            std::cout << MRQ_PREPRINT "Reading GAMS options: \n";
            
            ivalue = -1;
            r = reader.getIntOption("iterLim", ivalue);
            if(r == 0)
            {
                std::cout << "\t iterLim: " <<  ivalue << "\n";
                //our max number of iteratios is a long integer. So, if user let number of iterations negative, we do not set and use our default value ;)
                if(ivalue >= 0)
                    alg->in_max_iterations = ivalue; 
            }
            
            dvalue = 1e-6;
            r = reader.getDblOption("optCA", dvalue);
            if(r == 0)
            {
                std::cout << "\t optCA: " <<  dvalue <<"\n";
                alg->in_absolute_convergence_tol = dvalue;
            }
            
            dvalue = 1e-6;
            r = reader.getDblOption("optCR", dvalue);
            if(r == 0)
            {
                std::cout << "\t optCR: " <<  dvalue << "\n";
                alg->in_absolute_convergence_tol = dvalue;
            }
            
            dvalue = INFINITY;
            r = reader.getDblOption("resLim", dvalue);
            if(r == 0)
            {
                std::cout << "\t resLim: " <<  dvalue << "\n";
                alg->in_max_time = dvalue;
            }
            
            ivalue = 0;
            r = reader.getIntOption("threads", ivalue);
            if(r == 0)
            {
                std::cout << "\t threads: " <<  ivalue << "\n";
                alg->in_number_of_threads = ivalue;
            }
            
            std::cout << MRQ_PREPRINT "Done\n";
        }
        
        
        {
            std::string optFileName;
            
            int r = reader.getOptionFileName(optFileName);
            if( r == 0 )
            {
                
                std::cout << MRQ_PREPRINT "Trying read input parameter file " << optFileName << ". ";
                
                r = alg->readParametersWithTypeFromFile( optFileName.c_str(), true);
                
                if(r == 0)
                    std::cout << " Done\n" ;
                else
                    std::cout << " Failure\n" ;
            }
        
            if(r != 0) //do not put an else here necause que can get a failure inside if above
            {
                std::cout << MRQ_PREPRINT "Trying read input parameter file " MRQ_MURIQUI_PARAMS_FILE " .";
                
                r = alg->readParametersWithTypeFromFile( MRQ_MURIQUI_PARAMS_FILE, true);
                
                if(r == 0)
                    std::cout << " Done\n" ;
                else
                    std::cout << " Failure\n" ;
            }
        
        }
        
        
        {
            std::cout << MRQ_PREPRINT "Trying read milp solver input parameter file " MRQ_MILP_SOLVER_PARAMS_FILE  " .";
            
            r = milpParams.storeParametersFromFile( MRQ_MILP_SOLVER_PARAMS_FILE, true);
            if(r == 0)
                    std::cout << " Done\n" ;
                else
                    std::cout << " Failure\n" ;
        }
        
        
        {
            std::cout << MRQ_PREPRINT "Trying read nlp solver input parameter file " MRQ_NLP_SOLVER_PARAMS_FILE  " .";
            
            r = nlpParams.storeParametersFromFile( MRQ_NLP_SOLVER_PARAMS_FILE, true);
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
        
        
        
        
        
        alg->run(prob, &milpParams, &nlpParams);
        
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
        
        
        
        if(alg)
        {
            if( alg->out_best_obj < MRQ_INFINITY )
            {
                if( reader.isMaximizationProblem() )
                    alg->out_best_obj = -alg->out_best_obj;
            }
            
            
            #if MRQ_HAVE_GAMS
            
            /*gmoModelStat_OptimalGlobal        = 1,
            gmoModelStat_OptimalLocal         = 2,
            gmoModelStat_Unbounded            = 3,
            gmoModelStat_InfeasibleGlobal     = 4,
            gmoModelStat_InfeasibleLocal      = 5,
            gmoModelStat_InfeasibleIntermed   = 6,
            gmoModelStat_Feasible             = 7,
            gmoModelStat_Integer              = 8,
            gmoModelStat_NonIntegerIntermed   = 9,
            gmoModelStat_IntegerInfeasible    = 10,
            gmoModelStat_LicenseError         = 11,
            gmoModelStat_ErrorUnknown         = 12,
            gmoModelStat_ErrorNoSolution      = 13,
            gmoModelStat_NoSolutionReturned   = 14,
            gmoModelStat_SolvedUnique         = 15,
            gmoModelStat_Solved               = 16,
            gmoModelStat_SolvedSingular       = 17,
            gmoModelStat_UnboundedNoSolution  = 18,
            gmoModelStat_InfeasibleNoSolution = 19	*/					
            
            
            if( alg->out_return_code == MRQ_OPTIMAL_SOLUTION)
            {
                gamsModelReturnCode = gmoModelStat_OptimalLocal;
            }
            else 
            {
                if( alg->out_return_code == MRQ_UNBOUNDED_PROBLEM )
                {
                    gamsModelReturnCode = gmoModelStat_UnboundedNoSolution ;
                }
                else if( alg->out_return_code == MRQ_INFEASIBLE_PROBLEM )
                {
                    gamsModelReturnCode = gmoModelStat_InfeasibleLocal;
                }
                else
                {
                    gamsModelReturnCode = gmoModelStat_ErrorUnknown;
                }
                
                if( alg->out_feasible_solution )
                {
                    gamsModelReturnCode = gmoModelStat_Feasible; //we frgte other status and just set feasible
                }
            }
            
            
            
            /*
            *	gmoSolveStat_Normal      = 1,
                gmoSolveStat_Iteration   = 2,
                gmoSolveStat_Resource    = 3,
                gmoSolveStat_Solver      = 4,
                gmoSolveStat_EvalError   = 5,
                gmoSolveStat_Capability  = 6,
                gmoSolveStat_License     = 7,
                gmoSolveStat_User        = 8,
                gmoSolveStat_SetupErr    = 9,
                gmoSolveStat_SolverErr   = 10,
                gmoSolveStat_InternalErr = 11,
                gmoSolveStat_Skipped     = 12,
                gmoSolveStat_SystemErr   = 13  
            */ 
            
            
            if( alg->out_return_code == MRQ_OPTIMAL_SOLUTION || alg->out_return_code == MRQ_UNBOUNDED_PROBLEM || alg->out_return_code == MRQ_INFEASIBLE_PROBLEM )
            {
                gamsSolverReturnCode = gmoSolveStat_Normal;
            }
            else if( alg->out_return_code == MRQ_MAX_ITERATIONS_STOP )
            {
                gamsSolverReturnCode = gmoSolveStat_Iteration;
            }
            else if( alg->out_return_code == MRQ_MAX_TIME_STOP )
            {
                gamsSolverReturnCode = gmoSolveStat_SolverErr;
            }
            else if( alg->out_return_code == MRQ_CALLBACK_FUNCTION_ERROR )
            {
                gamsSolverReturnCode = gmoSolveStat_EvalError;
            }
            else if( alg->out_return_code == MRQ_ALG_NOT_APPLICABLE )
            {
                gamsSolverReturnCode = gmoSolveStat_SetupErr;
            }
            else if( alg->out_return_code == MRQ_BAD_PARAMETER_VALUES )
            {
                gamsSolverReturnCode = gmoSolveStat_SetupErr;
            }
            else if( alg->out_return_code == MRQ_STOP_REQUIRED_BY_USER )
            {
                gamsSolverReturnCode = gmoSolveStat_User;
            }
            else if( alg->out_return_code == MRQ_LIBRARY_NOT_AVAILABLE )
            {
                gamsSolverReturnCode = gmoSolveStat_SetupErr;
            }
            else if( alg->out_return_code == MRQ_MILP_SOLVER_ERROR || alg->out_return_code == MRQ_NLP_SOLVER_ERROR )
            {
                gamsSolverReturnCode = gmoSolveStat_SolverErr;
            }
            else if( alg->out_return_code == MRQ_MEMORY_ERROR )
            {
                gamsSolverReturnCode = gmoSolveStat_SystemErr;
            }
            else
            {
                gamsSolverReturnCode = gmoSolveStat_InternalErr;
            }
            
            #endif
            
            
            
            r = reader.setSolution(gamsModelReturnCode, gamsSolverReturnCode, alg->out_feasible_solution ? alg->out_best_sol : NULL);
            MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_UNDEFINED_ERROR, termination);
            
            //gmoSetHeadnTail(reader.gmo, gmoHiterused, alg->out_number_of_iterations);
            
            //gmoSetHeadnTail(reader.gmo, gmoHobjval, 1986.0);
            
            //std::cout << "gamsSolverReturnCode: " << gamsSolverReturnCode << " gamsModelReturnCode: " << gamsModelReturnCode << "\n";
            
            //r = reader.setSolverStatus(gamsSolverReturnCode);
            //MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_UNDEFINED_ERROR, termination);
            
            //r = reader.setModelStatus(gamsModelReturnCode);
            //MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_UNDEFINED_ERROR, termination);
            
            //MRQ_getchar();
            //gmoCompleteSolution( reader.gmo );
            
            //gmoWriteSolDone( reader.gmo, "iaia" );
            
            //gmoLoadSolutionLegacy(reader.gmo);
        }
    
        
        
        
    }
    #else
    {
        int r = MRQ_testRun(prob, stub, printAlgParameters);
        if(r != 0)
            MRQ_PRINTERRORNUMBER(r);
    }
    #endif
    
    
    
    
    retCode = 0;
    
termination:
    
    if(eval)	delete eval;
    if(outFile) fclose(outFile);
    
    return retCode;
}



