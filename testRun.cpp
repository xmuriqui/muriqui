/*
* File to test the algorithms in Muriqui Optimizer. It's purpose is making researches about MINLP algorithms.
* 
* 
* Author: Wendel Melo
* 
*/



#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include <iostream>
#include <vector>

#include "muriqui.hpp"

#include "MRQ_tools.hpp"


using namespace muriqui;



void inline MRQ_writeSolOnFile(FILE *file, double obj, const int n, const double *sol)
{
    fprintf( file, "%d\n%0.20f\n", n, obj );
    
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













int muriqui::MRQ_testRun(MRQ_MINLPProb &prob, const char *probName, const bool printAlgParameters)
{
    FILE *outFile = NULL;

    MRQ_GeneralSolverParams milpParams, nlpParams;

    MRQ_Algorithm *alg = NULL;

    MRQ_ContinuousRelax crelax;
    MRQ_ExtCutPlan ecp;
    MRQ_OuterApp oa;
    MRQ_LPNLPBBOuterApp bboa;
    MRQ_LPBBExtCutPlan bbecp;
    MRQ_LPBBExtSupHypPlane bbesh;
    MRQ_BonminHybrid bonminhyb;
    MRQ_ExtSupHypPlane eshp;
    MRQ_BranchAndBound newbb;
    MRQ_IGMA0 igma0;
    MRQ_IGMA1 igma1;
    MRQ_IGMA2 igma2;
    MRQ_Diving dive;
    MRQ_OAFeasibilityPump oafp;
    MRQ_FeasibilityPump fp;
    MRQ_RENS rens;
    MRQ_StructuredStochasticRounding ssr;

    const std::vector<MRQ_Algorithm *> algPointersVector = {&crelax, &ecp, &oa, &bboa, &bbecp, &bbesh, &bonminhyb, &eshp, &newbb, &igma0, &igma1, &igma2, &dive, &oafp, &fp, &rens, &ssr};

    
    
    #if 1 //to handle only binary problems 
    {
        if( !prob.isBinaryProblem() )
        {
            MRQ_PRINTERRORMSG("not threat this problem because it is not binary");
            return 0;
        }
        else
        {
            MRQ_PRINTERRORMSG("WARNING: BINARY FILTER IS ENABLED");
            MRQ_PRINTERRORMSG("WARNING: BINARY FILTER IS ENABLED");
            MRQ_PRINTERRORMSG("WARNING: BINARY FILTER IS ENABLED");
            MRQ_PRINTERRORMSG("WARNING: BINARY FILTER IS ENABLED");
            MRQ_PRINTERRORMSG("WARNING: BINARY FILTER IS ENABLED");
            MRQ_getchar();
        }
    }
    #endif
    


    #if MRQ_SAVE_OUTPUT_FILE
        outFile = fopen( MRQ_OUT_FILE_NAME, "a" );

        if(outFile)
        {
            MRQ_writeProblemParameters(probName, prob, outFile);
        }
        else
        {
            MRQ_PRINTERRORMSGP("Failure to open file: ", MRQ_OUT_FILE_NAME);
            goto termination;
        }
    #endif





    for(auto &algP: algPointersVector)
    {
        algP->in_print_parameters_values = printAlgParameters;
    }

    //if( prob.n <= 30 )
        //prob.print();

    //prob.print();
    //MRQ_getchar();

    //cplexParams.storeIntegerParameter("CPX_PARAM_THREADS", 1);
    //sParams.storeIntegerParameter("Threads", 1);
    //cplexParams.storeDoubleParameter("CPX_PARAM_TILIM", 4*60*60);
    //cplexParams.storeIntegerParameter("XPRS_THREADS", 1);
    //cplexParams.storeIntegerParameter("XPRS_MAXTIME", 600);

    //nlpParams.storeIntegerParameter("MSK_IPAR_INTPNT_MAX_ITERATIONS", 1000);
    //nlpParams.storeIntegerParameter("max_iter", 5000); //ipopt max number of iteration
    //nlpParams.storeStringParameter("derivative_test", "first-order");
    //nlpParams.storeStringParameter("derivative_test", "second-order");
    //nlpParams.storeIntegerParameter("print_level", 4);
    //mosekParams.storeIntegerParameter("MaxIter", 1500);
    //nlpParams.storeDoubleParameter("max_cpu_time", 120);
    //mosekParams.storeStringParameter("linear_solver", "pardiso");
    //mosekParams.storeStringParameter("linear_solver", "ma27"); //ma27 was faster than pardiso to solve continuous relaxation...

    //mosekParams.storeIntegerParameter("derivcheck", KTR_DERIVCHECK_ALL);
    //mosekParams.storeIntegerParameter("outlev", 3);
    //mosekParams.storeDoubleParameter("derivcheck_tol", 1e-4);


    #define MRQ_CHECK_CONVEX 0

    #if MRQ_CHECK_CONVEX
    {
        MRQ_ContinuousRelax crelax;
        MRQ_GeneralSolverParams mosekParams;
        
        //turn on mosek's convexity checker
        mosekParams.storeIntegerParameter("MSK_IPAR_CHECK_CONVEXITY", 1);
        mosekParams.storeIntegerParameter("MSK_IPAR_INTPNT_MAX_ITERATIONS", 10000);
        mosekParams.storeIntegerParameter("max_iter", 10000);
        
        
        crelax.in_nlp_solver = MRQ_NLP_MOSEK;
        //crelax.in_nlp_solver = MRQ_IPOPT;
        crelax.run(prob, NULL, &mosekParams);
        
        alg = &crelax;
        printf("____________________________________________________________________________\n");
        printf("Continuous Relaxation Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld solver retcode: %d\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, crelax.out_original_solver_return_code);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");
        
        
        //MRQ_getchar();
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                if(outFile)
                    fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%ld%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        if(crelax.out_nlp_feasible_solution)
        {
            char cmd[200];
            system("mkdir convex");
            
            sprintf(cmd, "mv %s convex", probName);
            system(cmd);
        }
        
        if(crelax.out_return_code == MRQ_INFEASIBLE_PROBLEM)
        {
            char cmd[200];
            system("mkdir infeasible");
            
            sprintf(cmd, "mv %s infeasible", probName);
            system(cmd);
        }
        
        
        if( crelax.out_original_solver_return_code == 1432 )
        {
            char cmd[200];
            system("mkdir eval_error");
            
            sprintf(cmd, "mv %s eval_error", probName);
            system(cmd);
        }
        
        
        int ret = crelax.out_original_solver_return_code;
        
        std::cerr << "crelax.out_return_code: " << crelax.out_return_code << " crelax.out_original_solver_return_code: " << ret << "\n";
        
        if(ret == 0)
        {
            char cmd[200];
            system("mkdir undefined");
            
            sprintf(cmd, "mv %s undefined", probName);
            system(cmd);
        }
        
        
        if(ret != 1290 && ret != 1291 && ret != 1293 && ret != 1294 && ret != 1295 && ret != 1432 && ret != 1500)
        {
            MRQ_getchar();
        }
        else
        {
            char cmd[200];
            system("mkdir nonconvex");
            
            sprintf(cmd, "mv %s nonconvex", probName);
            system(cmd);
        }
        
        
        goto termination;
    }
    #endif


    #define MRQ_CHECK_BINARY 0

    #if MRQ_CHECK_BINARY
    {
        const int nB = prob.getNumberOfBinaryVars();
        const int nI = prob.getNumberOfIntegerVars();
        
        if( nB == nI )
        {
            char cmd[1000];
            const char tdir[] = "binary_minlp";
            
            sprintf(cmd, "mkdir %s; mv %s %s", tdir, probName, tdir);
            system(cmd);
        }
        
        goto termination;
    }
    #endif


    #if 0 //checking derivative by minlpproblem
    {
        const double alx = -100, aux = 1000;
        const double *lx = prob.lx, *ux = prob.ux;
        
        bool answer;
        int r;
        double x[ prob.n ];
        
        for(int k = 0; k < 1000; k++)
        {
            MRQ_Random random(k);
            
            
            for(int i = 0; i < prob.n; i++)
            {
                x[i] = random.random( MRQ_max(alx, lx[i]), MRQ_min(aux, ux[i]) );
            }
            
            r = prob.checkFisrtDerivatives( prob.hasNlObj, prob.hasNLConstraints(), x, 1e-5, NULL, 1e-4, answer);
            
            std::cout << "First Derivative Checker. return code: " << r << " answer: " << answer << "\n";
            
            if (answer == false)
                MRQ_getchar();
        }
    }
    #endif


    #if 1 //continuous relaxation
    {
        //prob.print();
        
        
        crelax.in_number_of_threads = 1;
        crelax.in_nlp_solver = MRQ_IPOPT;
        crelax.in_set_special_nlp_solver_params = false;
        
        
        crelax.run(prob, NULL, &nlpParams);
        alg = &crelax;
        printf("____________________________________________________________________________\n");
        printf("Continuous Relaxation Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld solver: %s solver retcode: %d\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, optsolvers::OPT_getSolverName(crelax.in_nlp_solver).c_str(), crelax.out_original_solver_return_code );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n"); 
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                if(outFile)
                    fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%ld%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
            
        
        /*crelax.in_nlp_solver = MRQ_IPOPT;
        
        crelax.run(prob, NULL, &nlpParams);
        alg = &crelax;
        printf("____________________________________________________________________________\n");
        printf("Continuous Relaxation Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld solver: %s solver retcode: %d\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, optsolvers::OPT_getSolverName(crelax.in_nlp_solver).c_str(), crelax.out_original_solver_return_code );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n"); 
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                if(outFile)
                    fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%ld%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
            
        
            
        
        crelax.in_nlp_solver = MRQ_NLP_KNITRO;
        
        crelax.run(prob, NULL, &nlpParams);
        alg = &crelax;
        printf("____________________________________________________________________________\n");
        printf("Continuous Relaxation Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld solver: %s solver retcode: %d\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, optsolvers::OPT_getSolverName(crelax.in_nlp_solver).c_str(), crelax.out_original_solver_return_code );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n"); 
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                if(outFile)
                    fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%ld%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
            
        crelax.in_nlp_solver = MRQ_ALGENCAN;
        
        crelax.run(prob, NULL, &nlpParams);
        alg = &crelax;
        printf("____________________________________________________________________________\n");
        printf("Continuous Relaxation Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld solver: %s solver retcode: %d\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, optsolvers::OPT_getSolverName(crelax.in_nlp_solver).c_str(), crelax.out_original_solver_return_code );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n"); 
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                if(outFile)
                    fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%ld%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif 
        
        
        
        crelax.in_nlp_solver = MRQ_WORHP;
        
        crelax.run(prob, NULL, &nlpParams);
        alg = &crelax;
        printf("____________________________________________________________________________\n");
        printf("Continuous Relaxation Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld solver: %s solver retcode: %d\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, optsolvers::OPT_getSolverName(crelax.in_nlp_solver).c_str(), crelax.out_original_solver_return_code );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n"); 
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                if(outFile)
                    fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%ld%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
        */
        
    }
    #endif


    #if 1 //just to set up MILP libraries
    {
        MRQ_ExtCutPlan ecp;
        ecp.in_max_cpu_time = 20;
        ecp.in_max_iterations = 1;
        ecp.in_number_of_threads = 1;
        ecp.in_milp_solver = MRQ_CPLEX;
        ecp.in_refine_final_solution_using_nlp = false;
        ecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 0;
        
        ecp.run(prob, &milpParams, NULL);
        
        alg = &ecp;
        
        /*printf("____________________________________________________________________________\n");
        printf("ECP1 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %ld repetitions: %u all Py completelly solved: %d number of Py completelly solved: %ld \"nlp time\": %f nlp cpu time: %f\n", stub, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, ecp.out_number_of_milp_solver_iters, ecp.out_number_of_integer_solution_repetitions, (int) ecp.out_all_subproblems_py_completely_solved, ecp.out_number_of_subproblems_py_completely_solved, ecp.out_clock_time_of_nlp_solving, ecp.out_cpu_time_of_nlp_solving);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n"); */
        
        
        printf("\n\n");
    }
    #endif


    #if 0  //ecp
    {
        ecp.in_max_cpu_time = 48*60*60;
        //ecp.in_max_iterations = 20000;
        
        ecp.in_print_level = 5;
        ecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 0;
        ecp.in_milp_solver = MRQ_CPLEX;
        ecp.in_nlp_solver = MRQ_NLP_MOSEK;
        
        ecp.in_refine_final_solution_using_nlp = false;
        ecp.in_measure_nlp_time = true;
        //ecp.in_constr_linearisation_strategy = MRQ_CLS_ONLY_INFEAS_AND_ACTIVE;
        
        ecp.in_number_of_threads = 0;
        
        
        ecp.run(prob, &milpParams, NULL);
        
        alg = &ecp;
        
        printf("____________________________________________________________________________\n");
        printf("ECP0 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %ld repetitions: %lu all Py completelly solved: %d number of Py completelly solved: %ld \"nlp time\": %f nlp cpu time: %f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, ecp.out_number_of_milp_solver_iters, ecp.out_number_of_integer_solution_repetitions, (int) ecp.out_all_subproblems_py_completely_solved, ecp.out_number_of_subproblems_py_completely_solved, ecp.out_clock_time_of_nlp_solving, ecp.out_cpu_time_of_nlp_solving);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%f%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, ecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP, ecp.out_number_of_integer_solution_repetitions, MRQ_CHAR_SEP, (int) ecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP, ecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                ecp.out_cpu_time_of_nlp_solving,
                MRQ_CHAR_SEP,
                ecp.out_clock_time_of_nlp_solving,
                MRQ_CHAR_SEP
                );
            
            if(outFile)	    fflush(outFile);
        #endif
        
        if( ecp.out_return_code == MRQ_INFEASIBLE_PROBLEM )
        {
            goto termination;
        }
        
        
        /*ecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 1;
        ecp.in_measure_nlp_time = true;
        ecp.in_delete_intermediate_linearizations_in_each_iteration = false;
        
        ecp.run(prob, &milpParams, NULL);
        
        alg = &ecp;
        
        printf("____________________________________________________________________________\n");
        printf("ECP1 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %ld repetitions: %lu all Py completelly solved: %d number of Py completelly solved: %ld \"nlp time\": %f nlp cpu time: %f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, ecp.out_number_of_milp_solver_iters, ecp.out_number_of_integer_solution_repetitions, (int) ecp.out_all_subproblems_py_completely_solved, ecp.out_number_of_subproblems_py_completely_solved, ecp.out_clock_time_of_nlp_solving, ecp.out_cpu_time_of_nlp_solving);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%f%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, ecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP, ecp.out_number_of_integer_solution_repetitions, MRQ_CHAR_SEP, (int) ecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP, ecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                ecp.out_cpu_time_of_nlp_solving,
                MRQ_CHAR_SEP,
                ecp.out_clock_time_of_nlp_solving,
                MRQ_CHAR_SEP
                );
            
            if(outFile)	    fflush(outFile);
        #endif
        */
        
        /*
        ecp.in_constr_linearisation_strategy = MRQ_CLS_ONLY_INFEAS_AND_ACTIVE;
        ecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 0;
        
        ecp.run(prob, &milpParams, NULL);
        
        alg = &ecp;
        
        printf("____________________________________________________________________________\n");
        printf("ECP0 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %ld repetitions: %lu all Py completelly solved: %d number of Py completelly solved: %ld \"nlp time\": %f nlp cpu time: %f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, ecp.out_number_of_milp_solver_iters, ecp.out_number_of_integer_solution_repetitions, (int) ecp.out_all_subproblems_py_completely_solved, ecp.out_number_of_subproblems_py_completely_solved, ecp.out_clock_time_of_nlp_solving, ecp.out_cpu_time_of_nlp_solving);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%f%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, ecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP, ecp.out_number_of_integer_solution_repetitions, MRQ_CHAR_SEP, (int) ecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP, ecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                ecp.out_cpu_time_of_nlp_solving,
                MRQ_CHAR_SEP,
                ecp.out_clock_time_of_nlp_solving,
                MRQ_CHAR_SEP
                );
            
            if(outFile)	    fflush(outFile);
        #endif
        */
        
        
        /*ecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 1;
        
        ecp.run(prob, &milpParams, NULL);
        
        alg = &ecp;
        
        printf("____________________________________________________________________________\n");
        printf("ECP1 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %ld repetitions: %lu all Py completelly solved: %d number of Py completelly solved: %ld \"nlp time\": %f nlp cpu time: %f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, ecp.out_number_of_milp_solver_iters, ecp.out_number_of_integer_solution_repetitions, (int) ecp.out_all_subproblems_py_completely_solved, ecp.out_number_of_subproblems_py_completely_solved, ecp.out_clock_time_of_nlp_solving, ecp.out_cpu_time_of_nlp_solving);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%f%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, ecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP, ecp.out_number_of_integer_solution_repetitions, MRQ_CHAR_SEP, (int) ecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP, ecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                ecp.out_cpu_time_of_nlp_solving,
                MRQ_CHAR_SEP,
                ecp.out_clock_time_of_nlp_solving,
                MRQ_CHAR_SEP
                );
            
            if(outFile)	    fflush(outFile);
        #endif
        */
        
        
        /*ecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 5;
        
        ecp.run(prob, &milpParams, NULL);
        
        alg = &ecp;
        
        printf("____________________________________________________________________________\n");
        printf("ECP5 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %ld repetitions: %lu all Py completelly solved: %d number of Py completelly solved: %ld \"nlp time\": %f nlp cpu time: %f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, ecp.out_number_of_milp_solver_iters, ecp.out_number_of_integer_solution_repetitions, (int) ecp.out_all_subproblems_py_completely_solved, ecp.out_number_of_subproblems_py_completely_solved, ecp.out_clock_time_of_nlp_solving, ecp.out_cpu_time_of_nlp_solving);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%f%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, ecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP, ecp.out_number_of_integer_solution_repetitions, MRQ_CHAR_SEP, (int) ecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP, ecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                ecp.out_cpu_time_of_nlp_solving,
                MRQ_CHAR_SEP,
                ecp.out_clock_time_of_nlp_solving,
                MRQ_CHAR_SEP
                );
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        ecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 10;
        
        
        ecp.run(prob, &milpParams, NULL);
        
        alg = &ecp;
        
        printf("____________________________________________________________________________\n");
        printf("ECP10 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %ld repetitions: %lu all Py completelly solved: %d number of Py completelly solved: %ld \"nlp time\": %f nlp cpu time: %f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, ecp.out_number_of_milp_solver_iters, ecp.out_number_of_integer_solution_repetitions, (int) ecp.out_all_subproblems_py_completely_solved, ecp.out_number_of_subproblems_py_completely_solved, ecp.out_clock_time_of_nlp_solving, ecp.out_cpu_time_of_nlp_solving);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%f%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, ecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP, ecp.out_number_of_integer_solution_repetitions, MRQ_CHAR_SEP, (int) ecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP, ecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                ecp.out_cpu_time_of_nlp_solving,
                MRQ_CHAR_SEP,
                ecp.out_clock_time_of_nlp_solving,
                MRQ_CHAR_SEP
                );
            
            if(outFile)	    fflush(outFile);
        #endif */
    }
    #endif
    
    
    #if 0  //oa
    {
        //testes sobre oa
        //oa.in_print_level = 10;
        //oa.in_round_first_nlp_relaxation_solution = true;
        //oa.in_set_quadratics_in_master_problem = true;
        
        //oa.in_set_obj_lower_bound_on_master_problem = true;
        
        oa.in_max_cpu_time = 4*60*60;
        //oa.in_max_time = 120;
        oa.in_number_of_threads = 1;
        oa.in_measure_nlp_time = true;
        oa.in_milp_solver = MRQ_CPLEX;
        oa.in_nlp_solver = MRQ_NLP_MOSEK;
        
        
        oa.run(prob, &milpParams, &nlpParams);
        
        alg = &oa;
        printf("____________________________________________________________________________\n");
        printf("OA Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %ld \nobj lin saved: %u, %u constr lin saved: %u cpu time to box to constr: %0.2f nlp subprobs: %ld cpu time of nlp: %0.2f clock time of nlp: %0.2f iter to first feas sol: %lu iter to last feas sol: %lu cpu time to first feas sol: %0.2f clock time to first feas sol: %0.2f \n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, oa.out_number_of_milp_solver_iters, oa.out_number_of_obj_linears_saved_by_zu, oa.out_number_of_obj_linears_saved_by_infeasibility, oa.out_number_of_constr_linears_saved,  oa.out_cpu_time_on_box_to_constr_calculation,
        oa.out_number_of_nlp_probs_solved,
        oa.out_cpu_time_of_nlp_solving, oa.out_clock_time_of_nlp_solving, alg->out_number_of_iterations_to_first_feas_sol, alg->out_number_of_iterations_to_best_sol, alg->out_cpu_time_to_first_feas_sol, alg->out_clock_time_to_first_feas_sol );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]= %0.10f;\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d" MRQ_CHAR_SEP  "%d" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%ld" MRQ_CHAR_SEP "%ld" MRQ_CHAR_SEP "%u" MRQ_CHAR_SEP "%u" MRQ_CHAR_SEP "%u" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%ld" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP,  alg->out_algorithm,  alg->out_return_code, alg->out_lower_bound, alg->out_best_obj, alg->out_cpu_time, alg->out_clock_time,  alg->out_number_of_iterations,
                oa.out_number_of_milp_solver_iters,
                oa.out_number_of_obj_linears_saved_by_zu, oa.out_number_of_obj_linears_saved_by_infeasibility, oa.out_number_of_constr_linears_saved, oa.out_cpu_time_on_box_to_constr_calculation, oa.out_number_of_nlp_probs_solved, oa.out_cpu_time_of_nlp_solving, oa.out_clock_time_of_nlp_solving );
            
            if(outFile)	    fflush(outFile);
        #endif
        
        /*oa.in_constr_linearisation_strategy = MRQ_CLS_ONLY_INFEAS_AND_ACTIVE;
        oa.run(prob, &milpParams, &nlpParams);
        
        alg = &oa;
        printf("____________________________________________________________________________\n");
        printf("OA Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %ld \nobj lin saved: %u, %u constr lin saved: %u cpu time to box to constr: %0.2f nlp subprobs: %ld cpu time of nlp: %0.2f clock time of nlp: %0.2f iter to first feas sol: %lu iter to last feas sol: %lu cpu time to first feas sol: %0.2f clock time to first feas sol: %0.2f \n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, oa.out_number_of_milp_solver_iters, oa.out_number_of_obj_linears_saved_by_zu, oa.out_number_of_obj_linears_saved_by_infeasibility, oa.out_number_of_constr_linears_saved,  oa.out_cpu_time_on_box_to_constr_calculation,
        oa.out_number_of_nlp_probs_solved,
        oa.out_cpu_time_of_nlp_solving, oa.out_clock_time_of_nlp_solving, alg->out_number_of_iterations_to_first_feas_sol, alg->out_number_of_iterations_to_best_sol, alg->out_cpu_time_to_first_feas_sol, alg->out_clock_time_to_first_feas_sol );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]= %0.10f;\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d" MRQ_CHAR_SEP  "%d" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%ld" MRQ_CHAR_SEP "%ld" MRQ_CHAR_SEP "%u" MRQ_CHAR_SEP "%u" MRQ_CHAR_SEP "%u" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%ld" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP,  alg->out_algorithm,  alg->out_return_code, alg->out_lower_bound, alg->out_best_obj, alg->out_cpu_time, alg->out_clock_time,  alg->out_number_of_iterations,
                oa.out_number_of_milp_solver_iters,
                oa.out_number_of_obj_linears_saved_by_zu, oa.out_number_of_obj_linears_saved_by_infeasibility, oa.out_number_of_constr_linears_saved, oa.out_cpu_time_on_box_to_constr_calculation, oa.out_number_of_nlp_probs_solved, oa.out_cpu_time_of_nlp_solving, oa.out_clock_time_of_nlp_solving );
            
            if(outFile)	    fflush(outFile);
        #endif */
    }
    #endif


    #if 0  //bboa
    {
        bboa.in_number_of_threads = 1;
        bboa.in_max_cpu_time = 8*60*60;
        bboa.in_measure_nlp_time = true;
        bboa.in_nlp_solver = MRQ_NLP_MOSEK;
        bboa.in_milp_solver = MRQ_CPLEX;
        
        bboa.run(prob, &milpParams, &nlpParams);
        
        alg = &bboa;
        printf("____________________________________________________________________________\n");
        printf("BBOA Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %ld nlp subprobs: %ld cpu time of nlp: %0.2f clock time of nlp: %0.2f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bboa.out_number_of_milp_solver_iters, bboa.out_number_of_nlp_probs_solved, bboa.out_cpu_time_of_nlp_solving, bboa.out_clock_time_of_nlp_solving);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%ld%s" "%ld%s" "%ld%s" "%f%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP,
                bboa.out_number_of_milp_solver_iters, MRQ_CHAR_SEP,
                bboa.out_number_of_nlp_probs_solved, MRQ_CHAR_SEP, bboa.out_cpu_time_of_nlp_solving, MRQ_CHAR_SEP, bboa.out_clock_time_of_nlp_solving, MRQ_CHAR_SEP );
            
            if(outFile)	    fflush(outFile);
        #endif
        
        bboa.in_constr_linearisation_strategy = MRQ_CLS_ONLY_INFEAS_AND_ACTIVE;
        
        bboa.run(prob, &milpParams, &nlpParams);
        
        alg = &bboa;
        printf("____________________________________________________________________________\n");
        printf("BBOA Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %ld nlp subprobs: %ld cpu time of nlp: %0.2f clock time of nlp: %0.2f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bboa.out_number_of_milp_solver_iters, bboa.out_number_of_nlp_probs_solved, bboa.out_cpu_time_of_nlp_solving, bboa.out_clock_time_of_nlp_solving);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%ld%s" "%ld%s" "%ld%s" "%f%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP,
                bboa.out_number_of_milp_solver_iters, MRQ_CHAR_SEP,
                bboa.out_number_of_nlp_probs_solved, MRQ_CHAR_SEP, bboa.out_cpu_time_of_nlp_solving, MRQ_CHAR_SEP, bboa.out_clock_time_of_nlp_solving, MRQ_CHAR_SEP );
            
            if(outFile)	    fflush(outFile);
        #endif
        }
    #endif


    #if 0  //bbecp to test pseudo prune
    {
        bbecp.in_number_of_threads = 1;
        bbecp.in_max_cpu_time = 48*60*60;
        bbecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 0;
        bbecp.in_milp_solver = MRQ_CPLEX;
        bbecp.in_nlp_solver = MRQ_IPOPT;
        bbecp.in_measure_nlp_time = true;
        //bbecp.in_print_level = 5;
        bbecp.in_refine_final_solution_using_nlp = true;
        
        bbecp.in_pseudo_pruning_strategy = MRQ_BB_PPS_ON_NODE_EXPLORATION_AND_BRANCHING;
        
        
        
        
        bbecp.run(prob, &milpParams, NULL);
        
        alg = &bbecp;
        printf("____________________________________________________________________________\n");
        printf("BBECP 0 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu milp iters: %lu lp subprobs: %lu all Py completelly solved: %d number of Py completelly solved: %lu \"nlp time\": %f nlp cpu time: %f pseudo prune 1: %lu pseudo prune 2: %lu lower bound: %lf\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bbecp.out_number_of_milp_solver_iters, bbecp.out_number_of_master_relaxation_solved, (int) bbecp.out_all_subproblems_py_completely_solved, bbecp.out_number_of_subproblems_py_completely_solved, bbecp.out_clock_time_of_nlp_solving, bbecp.out_cpu_time_of_nlp_solving, bbecp.out_number_of_pseudo_prunes_on_solving, bbecp.out_number_of_pseudo_prunes_on_branching, alg->out_lower_bound);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%lu%s" "%lu%s" "%f%s" "%f%s", 
                alg->out_algorithm, MRQ_CHAR_SEP, 
                alg->out_return_code, MRQ_CHAR_SEP,
                alg->out_lower_bound, MRQ_CHAR_SEP,
                alg->out_best_obj, MRQ_CHAR_SEP,
                alg->out_cpu_time, MRQ_CHAR_SEP,
                alg->out_clock_time, MRQ_CHAR_SEP, 
                alg->out_number_of_iterations, MRQ_CHAR_SEP, bbecp.out_number_of_master_relaxation_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP,
                (int) bbecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                
                bbecp.out_number_of_pseudo_prunes_on_solving, MRQ_CHAR_SEP,
                bbecp.out_number_of_pseudo_prunes_on_branching, MRQ_CHAR_SEP,
                
                bbecp.out_cpu_time_of_nlp_solving, MRQ_CHAR_SEP,
                bbecp.out_clock_time_of_nlp_solving, MRQ_CHAR_SEP  );
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        bbecp.in_try_pseudo_pruning_before_solving_relaxations = true;
        bbecp.run(prob, &milpParams, NULL);
        
        alg = &bbecp;
        printf("____________________________________________________________________________\n");
        printf("BBECP 1 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu milp iters: %lu lp subprobs: %lu all Py completelly solved: %d number of Py completelly solved: %lu \"nlp time\": %f nlp cpu time: %f pseudo prune 1: %lu pseudo prune 2: %lu lower bound: %lf\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bbecp.out_number_of_milp_solver_iters, bbecp.out_number_of_master_relaxation_solved, (int) bbecp.out_all_subproblems_py_completely_solved, bbecp.out_number_of_subproblems_py_completely_solved, bbecp.out_clock_time_of_nlp_solving, bbecp.out_cpu_time_of_nlp_solving, bbecp.out_number_of_pseudo_prunes_on_solving, bbecp.out_number_of_pseudo_prunes_on_branching, alg->out_lower_bound);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%lu%s" "%lu%s" "%f%s" "%f%s", 
                alg->out_algorithm, MRQ_CHAR_SEP, 
                alg->out_return_code, MRQ_CHAR_SEP,
                alg->out_lower_bound, MRQ_CHAR_SEP,
                alg->out_best_obj, MRQ_CHAR_SEP,
                alg->out_cpu_time, MRQ_CHAR_SEP,
                alg->out_clock_time, MRQ_CHAR_SEP, 
                alg->out_number_of_iterations, MRQ_CHAR_SEP, bbecp.out_number_of_master_relaxation_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP,
                (int) bbecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                
                bbecp.out_number_of_pseudo_prunes_on_solving, MRQ_CHAR_SEP,
                bbecp.out_number_of_pseudo_prunes_on_branching, MRQ_CHAR_SEP,
                
                bbecp.out_cpu_time_of_nlp_solving, MRQ_CHAR_SEP,
                bbecp.out_clock_time_of_nlp_solving, MRQ_CHAR_SEP  );
            
            if(outFile)	    fflush(outFile);
        #endif
            
            
            
            
        bbecp.in_pseudo_pruning_strategy = MRQ_BB_PPS_NO_PSEUDO_PRUNING;
        bbecp.run(prob, &milpParams, NULL);
        
        alg = &bbecp;
        printf("____________________________________________________________________________\n");
        printf("BBECP 2 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu milp iters: %lu lp subprobs: %lu all Py completelly solved: %d number of Py completelly solved: %lu \"nlp time\": %f nlp cpu time: %f pseudo prune 1: %lu pseudo prune 2: %lu lower bound: %lf\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bbecp.out_number_of_milp_solver_iters, bbecp.out_number_of_master_relaxation_solved, (int) bbecp.out_all_subproblems_py_completely_solved, bbecp.out_number_of_subproblems_py_completely_solved, bbecp.out_clock_time_of_nlp_solving, bbecp.out_cpu_time_of_nlp_solving, bbecp.out_number_of_pseudo_prunes_on_solving, bbecp.out_number_of_pseudo_prunes_on_branching, alg->out_lower_bound);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%lu%s" "%lu%s" "%f%s" "%f%s", 
                alg->out_algorithm, MRQ_CHAR_SEP, 
                alg->out_return_code, MRQ_CHAR_SEP,
                alg->out_lower_bound, MRQ_CHAR_SEP,
                alg->out_best_obj, MRQ_CHAR_SEP,
                alg->out_cpu_time, MRQ_CHAR_SEP,
                alg->out_clock_time, MRQ_CHAR_SEP, 
                alg->out_number_of_iterations, MRQ_CHAR_SEP, bbecp.out_number_of_master_relaxation_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP,
                (int) bbecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                
                bbecp.out_number_of_pseudo_prunes_on_solving, MRQ_CHAR_SEP,
                bbecp.out_number_of_pseudo_prunes_on_branching, MRQ_CHAR_SEP,
                
                bbecp.out_cpu_time_of_nlp_solving, MRQ_CHAR_SEP,
                bbecp.out_clock_time_of_nlp_solving, MRQ_CHAR_SEP  );
            
            if(outFile)	    fflush(outFile);
        #endif
    }
    #endif
    
    
    #if 0  //bbecp
    {
        
        //bbecp.in_constr_linearisation_strategy = MRQ_CLS_ONLY_INFEAS_AND_ACTIVE;
        bbecp.in_max_cpu_time = 600;
        bbecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 0;
        bbecp.in_number_of_threads = 1;
        
        bbecp.run(prob, &milpParams, NULL);
        
        alg = &bbecp;
        printf("____________________________________________________________________________\n");
        printf("BBECP 0 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu milp iters: %lu lp subprobs: %lu all Py completelly solved: %d number of Py completelly solved: %lu \"nlp time\": %f nlp cpu time: %f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bbecp.out_number_of_milp_solver_iters, bbecp.out_number_of_master_relaxation_solved, (int) bbecp.out_all_subproblems_py_completely_solved, bbecp.out_number_of_subproblems_py_completely_solved, bbecp.out_clock_time_of_nlp_solving, bbecp.out_cpu_time_of_nlp_solving);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%f%s" "%f%s", 
                alg->out_algorithm, MRQ_CHAR_SEP, 
                alg->out_return_code, MRQ_CHAR_SEP,
                alg->out_lower_bound, MRQ_CHAR_SEP,
                alg->out_best_obj, MRQ_CHAR_SEP,
                alg->out_cpu_time, MRQ_CHAR_SEP,
                alg->out_clock_time, MRQ_CHAR_SEP, 
                alg->out_number_of_iterations, MRQ_CHAR_SEP, bbecp.out_number_of_master_relaxation_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP,
                (int) bbecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
            
                bbecp.out_cpu_time_of_nlp_solving,
                MRQ_CHAR_SEP,
                bbecp.out_clock_time_of_nlp_solving,
                MRQ_CHAR_SEP  );
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        /*
        bbecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 1;
        
        bbecp.run(prob, &milpParams, NULL);
        
        alg = &bbecp;
        printf("____________________________________________________________________________\n");
        printf("BBECP 1 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu milp iters: %lu lp subprobs: %lu all Py completelly solved: %d number of Py completelly solved: %lu \"nlp time\": %f nlp cpu time: %f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bbecp.out_number_of_milp_solver_iters, bbecp.out_number_of_master_relaxation_solved, (int) bbecp.out_all_subproblems_py_completely_solved, bbecp.out_number_of_subproblems_py_completely_solved, bbecp.out_clock_time_of_nlp_solving, bbecp.out_cpu_time_of_nlp_solving);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%f%s" "%f%s", 
                alg->out_algorithm, MRQ_CHAR_SEP, 
                alg->out_return_code, MRQ_CHAR_SEP,
                alg->out_lower_bound, MRQ_CHAR_SEP,
                alg->out_best_obj, MRQ_CHAR_SEP,
                alg->out_cpu_time, MRQ_CHAR_SEP,
                alg->out_clock_time, MRQ_CHAR_SEP, 
                alg->out_number_of_iterations, MRQ_CHAR_SEP, bbecp.out_number_of_master_relaxation_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP,
                (int) bbecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
            
                bbecp.out_cpu_time_of_nlp_solving,
                MRQ_CHAR_SEP,
                bbecp.out_clock_time_of_nlp_solving,
                MRQ_CHAR_SEP  );
            
            if(outFile)	    fflush(outFile);
        #endif */
        
        
        /*bbecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 5;
        
        bbecp.run(prob, &milpParams, NULL);
        
        alg = &bbecp;
        printf("____________________________________________________________________________\n");
        printf("BBECP 5 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu milp iters: %lu lp subprobs: %lu all Py completelly solved: %d number of Py completelly solved: %lu \"nlp time\": %f nlp cpu time: %f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bbecp.out_number_of_milp_solver_iters, bbecp.out_number_of_master_relaxation_solved, (int) bbecp.out_all_subproblems_py_completely_solved, bbecp.out_number_of_subproblems_py_completely_solved, bbecp.out_clock_time_of_nlp_solving, bbecp.out_cpu_time_of_nlp_solving);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%f%s" "%f%s", 
                alg->out_algorithm, MRQ_CHAR_SEP, 
                alg->out_return_code, MRQ_CHAR_SEP,
                alg->out_lower_bound, MRQ_CHAR_SEP,
                alg->out_best_obj, MRQ_CHAR_SEP,
                alg->out_cpu_time, MRQ_CHAR_SEP,
                alg->out_clock_time, MRQ_CHAR_SEP, 
                alg->out_number_of_iterations, MRQ_CHAR_SEP, bbecp.out_number_of_master_relaxation_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP,
                (int) bbecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
            
                bbecp.out_cpu_time_of_nlp_solving,
                MRQ_CHAR_SEP,
                bbecp.out_clock_time_of_nlp_solving,
                MRQ_CHAR_SEP  );
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        
        
        bbecp.in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 10;
        
        bbecp.run(prob, &milpParams, NULL);
        
        alg = &bbecp;
        printf("____________________________________________________________________________\n");
        printf("BBECP 10 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu milp iters: %lu lp subprobs: %lu all Py completelly solved: %d number of Py completelly solved: %lu \"nlp time\": %f nlp cpu time: %f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bbecp.out_number_of_milp_solver_iters, bbecp.out_number_of_master_relaxation_solved, (int) bbecp.out_all_subproblems_py_completely_solved, bbecp.out_number_of_subproblems_py_completely_solved, bbecp.out_clock_time_of_nlp_solving, bbecp.out_cpu_time_of_nlp_solving);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%d%s" "%lu%s" "%f%s" "%f%s", 
                alg->out_algorithm, MRQ_CHAR_SEP, 
                alg->out_return_code, MRQ_CHAR_SEP,
                alg->out_lower_bound, MRQ_CHAR_SEP,
                alg->out_best_obj, MRQ_CHAR_SEP,
                alg->out_cpu_time, MRQ_CHAR_SEP,
                alg->out_clock_time, MRQ_CHAR_SEP, 
                alg->out_number_of_iterations, MRQ_CHAR_SEP, bbecp.out_number_of_master_relaxation_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP,
                (int) bbecp.out_all_subproblems_py_completely_solved, MRQ_CHAR_SEP,
                bbecp.out_number_of_subproblems_py_completely_solved, MRQ_CHAR_SEP,
            
                bbecp.out_cpu_time_of_nlp_solving,
                MRQ_CHAR_SEP,
                bbecp.out_clock_time_of_nlp_solving,
                MRQ_CHAR_SEP  );
            
            if(outFile)	    fflush(outFile);
        #endif */
        }
    #endif 


    #if 0  //bbesh
    {
        bbesh.in_number_of_threads = 1;
        bbesh.in_max_cpu_time = 4*60*60;
        bbesh.in_fix_int_vars_to_try_improve_cut = false;
        bbesh.in_print_level = 4;
        bbesh.in_milp_solver = MRQ_CPLEX;
        bbesh.in_nlp_solver = MRQ_NLP_MOSEK;
        //bbecp.in_absolute_feasibility_tol = 5e-3;
        
        //cplexParams.storeDoubleParameter( "FeasibilityTol", 1e-8 );
        //cplexParams.storeIntegerParameter( "OutputFlag", 1 );
        
        bbesh.run(prob, &milpParams, NULL);
        
        alg = &bbesh;
        printf("____________________________________________________________________________\n");
        printf("BBESH Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu milp iters: %lu lp subprobs: %lu\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bbesh.out_number_of_milp_solver_iters, bbesh.out_number_of_master_relaxation_solved);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" , alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, bbesh.out_number_of_master_relaxation_solved, MRQ_CHAR_SEP
                    );
            
            if(outFile)	    fflush(outFile);
        #endif
        
        bbesh.in_max_lp_subiters = 50;
        bbesh.run(prob, &milpParams, NULL);
        
        alg = &bbesh;
        printf("____________________________________________________________________________\n");
        printf("BBESH Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu milp iters: %lu lp subprobs: %lu\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bbesh.out_number_of_milp_solver_iters, bbesh.out_number_of_master_relaxation_solved);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" , alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, bbesh.out_number_of_master_relaxation_solved, MRQ_CHAR_SEP
                    );
            
            if(outFile)	    fflush(outFile);
        #endif
    }
    #endif 


    #if 0  //eshp
    {
        eshp.in_fix_int_vars_to_try_improve_cut = false;
        eshp.in_print_level = 3;
        eshp.in_printing_frequency = 20;
        eshp.in_max_cpu_time = 4*60*60;
        //eshp.in_eps_to_line_search = 0.0001;
        
        eshp.in_interior_point_strategy = MRQ_ESHP_IPS_MOST_INTERIOR;
        eshp.in_constr_linearisation_strategy = MRQ_CLS_ALL_CONSTRS;
        //eshp.in_eps_to_enforce_interior_sol_on_cont_relax_sol = 0.2;
        
        //eshp.in_absolute_tol_to_check_previous_sol = 1e-3;
        //eshp.in_relative_tol_to_check_previous_sol = 1e-3;
        //eshp.in_try_solve_interior_problem_if_cont_relax_fail = false;
        
        eshp.in_number_of_threads = 1;
        eshp.in_measure_nlp_time = true;
        eshp.in_nlp_solver = MRQ_NLP_MOSEK;
        eshp.in_milp_solver = MRQ_CPLEX;
        
        //eshp.in_integer_tol = 1e-6;
        
        eshp.run(prob, &milpParams, &nlpParams);
        
        alg = &eshp;
        printf("____________________________________________________________________________\n");
        printf("ESHP Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu milp iters: %lu nlp subprobs: %lu cpu time of nlp: %0.2f clock time of nlp: %0.2f \n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, eshp.out_number_of_milp_solver_iters, eshp.out_number_of_nlp_probs_solved,
        eshp.out_cpu_time_of_nlp_solving, eshp.out_clock_time_of_nlp_solving);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%f%s" "%u%s" "%lu%s" "%f%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP,
                eshp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP,
                eshp.out_obj_opt_at_continuous_relax, MRQ_CHAR_SEP, eshp.out_number_of_lp_iterations, MRQ_CHAR_SEP, eshp.out_number_of_nlp_probs_solved, MRQ_CHAR_SEP,
                eshp.out_cpu_time_of_nlp_solving, MRQ_CHAR_SEP, eshp.out_clock_time_of_nlp_solving, MRQ_CHAR_SEP
                );
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        /*eshp.in_fix_int_vars_to_try_improve_cut = true;
        eshp.run(prob, &cplexParams, &mosekParams);
        
        alg = &eshp;
        printf("____________________________________________________________________________\n");
        printf("ESHP Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu milp iters: %lu nlp subprobs: %lu cpu time of nlp: %0.2f clock time of nlp: %0.2f \n", stub, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, eshp.out_number_of_milp_solver_iters, eshp.out_number_of_nlp_probs_solved,
        eshp.out_cpu_time_of_nlp_solving, eshp.out_clock_time_of_nlp_solving);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%f%s" "%u%s" "%lu%s" "%f%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP,
                eshp.out_number_of_milp_solver_iters, MRQ_CHAR_SEP,
                eshp.out_obj_opt_at_continuous_relax, MRQ_CHAR_SEP, eshp.out_number_of_lp_iterations, MRQ_CHAR_SEP, eshp.out_number_of_nlp_probs_solved, MRQ_CHAR_SEP,
                eshp.out_cpu_time_of_nlp_solving, MRQ_CHAR_SEP, eshp.out_clock_time_of_nlp_solving, MRQ_CHAR_SEP
                );
            
            if(outFile)	    fflush(outFile);
        #endif */
        
        
    }
    #endif


    #if 0  //bonminhyb
    {
        bonminhyb.in_outer_app_max_cpu_time = 30.0;
        bonminhyb.in_max_cpu_time = 4*60*60;
        bonminhyb.in_number_of_threads = 1;
        bonminhyb.in_measure_nlp_time = true;
        bonminhyb.in_nlp_solver = MRQ_NLP_MOSEK;
        bonminhyb.in_milp_solver = MRQ_CPLEX;
        
        bonminhyb.run(prob, &milpParams, &nlpParams);
        
        alg = &bonminhyb;
        printf("____________________________________________________________________________\n");
        printf("Bonmin hybrid Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld milp iters: %lu nlp subprobs: %lu cpu time of nlp: %0.2f clock time of nlp: %0.2f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, bonminhyb.out_number_of_milp_solver_iters, bonminhyb.out_number_of_nlp_probs_solved, bonminhyb.out_cpu_time_of_nlp_solving, bonminhyb.out_clock_time_of_nlp_solving);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%f%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, bonminhyb.out_number_of_milp_solver_iters, MRQ_CHAR_SEP, bonminhyb.out_number_of_nlp_probs_solved, MRQ_CHAR_SEP, bonminhyb.out_cpu_time_of_nlp_solving, MRQ_CHAR_SEP, bonminhyb.out_clock_time_of_nlp_solving, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
    }
    #endif

    
    #if 0 //to test branch-bound with ssr
    {
        newbb.in_number_of_threads = 1;
        newbb.in_print_level = 5;
        newbb.in_printing_frequency = 1;
        
        newbb.in_rounding_heuristic_strategy = MRQ_RS_SSR;
        newbb.in_rounding_heuristic_call_iter_frequence = 100;
        
        newbb.in_use_outer_app = false;
        newbb.in_use_outer_app_as_heuristic = false;
        
        newbb.in_set_special_nlp_solver_params = false;
        
        newbb.run(prob, &milpParams, &nlpParams);
        
        alg = &newbb;
        printf("____________________________________________________________________________\n");
        printf("New BB - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu lower bound: %0.2f iter to first feas sol: %lu iter to last feas sol: %lu cpu time to first feas sol: %0.2f clock time to first feas sol: %0.2f \n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, alg->out_lower_bound, alg->out_number_of_iterations_to_first_feas_sol, alg->out_number_of_iterations_to_best_sol, alg->out_cpu_time_to_first_feas_sol, alg->out_clock_time_to_first_feas_sol );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        std::cout << "out_prune_counter:: bound: " << newbb.out_prune_counter.bound << " infeas: " << newbb.out_prune_counter.infeas << " opt: " << newbb.out_prune_counter.opt << " user: " << newbb.out_prune_counter.user << " bad lower bounds: " << newbb.out_nlp_failure_in_some_relaxation << " before solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning1 << " after solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning2 << "\n";
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, newbb.out_number_of_iters_having_wrong_lower_bounds, MRQ_CHAR_SEP, alg->out_number_of_iterations_to_best_sol, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning1, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning2, MRQ_CHAR_SEP, newbb.out_prune_counter.bound, MRQ_CHAR_SEP, newbb.out_prune_counter.infeas, MRQ_CHAR_SEP, newbb.out_prune_counter.opt, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
    }
    #endif
    
    
    #if 0  //branch-and-bound: to test steiner
    {
        //mosekParams.storeIntegerParameter("MSK_IPAR_INTPNT_MAX_ITERATIONS", 1000);
        
        //mosekParams.storeStringParameter("linear_solver", "ma27");
        //mosekParams.storeIntegerParameter("max_iter", 1000);
        
        //MRQ_NewBB newbb;
        newbb.in_number_of_threads = 1;
        newbb.in_print_level = 5;
        //newbb.in_printing_frequency = 500000;
        //newbb.in_lower_bound = MRQ_INFINITY*0.1;
        
        newbb.in_use_outer_app = false;
        newbb.in_use_outer_app_as_heuristic = false;
        //newbb.in_outer_app_subprob_frequence = 4;
        
        //newbb.in_stop_multibranch_after_first_bound_prune = false;
        
        newbb.in_max_cpu_time = 24*60*60; //24 hours!!!
        
        
        //newbb.in_constr_branching_strategy = MRQ_BB_CBS_NO_CONSTRAINT_BRANCH;// MRQ_BB_CBS_LOWEST_NUMBER_OF_VARS;
        newbb.in_nlp_solver = MRQ_NLP_MOSEK;
        newbb.in_milp_solver = MRQ_CPLEX;
        newbb.in_igma2_gap_min_solver = MRQ_IPOPT;
        
        
        newbb.in_set_special_nlp_solver_params = false;
        
        newbb.run(prob, &milpParams, &nlpParams);
        
        alg = &newbb;
        printf("____________________________________________________________________________\n");
        printf("New BB - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu lower bound: %0.2f iter to first feas sol: %lu iter to last feas sol: %lu cpu time to first feas sol: %0.2f clock time to first feas sol: %0.2f \n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, alg->out_lower_bound, alg->out_number_of_iterations_to_first_feas_sol, alg->out_number_of_iterations_to_best_sol, alg->out_cpu_time_to_first_feas_sol, alg->out_clock_time_to_first_feas_sol );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        std::cout << "out_prune_counter:: bound: " << newbb.out_prune_counter.bound << " infeas: " << newbb.out_prune_counter.infeas << " opt: " << newbb.out_prune_counter.opt << " user: " << newbb.out_prune_counter.user << " bad lower bounds: " << newbb.out_nlp_failure_in_some_relaxation << " before solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning1 << " after solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning2 << "\n";
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, newbb.out_number_of_iters_having_wrong_lower_bounds, MRQ_CHAR_SEP, alg->out_number_of_iterations_to_best_sol, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning1, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning2, MRQ_CHAR_SEP, newbb.out_prune_counter.bound, MRQ_CHAR_SEP, newbb.out_prune_counter.infeas, MRQ_CHAR_SEP, newbb.out_prune_counter.opt, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
    }
    #endif
    
    
    #if 0 //branch-and-bound: to test binary variable fixing by dual nlp relax solution
    {
        newbb.in_number_of_threads = 1;
        newbb.in_printing_frequency = 1;
        newbb.in_fix_int_vars_from_nlp_relax_sol = true;
        
        newbb.in_max_cpu_time = 4*60*60;
        
        newbb.in_nlp_solver = MRQ_NLP_MOSEK;
        newbb.in_milp_solver = MRQ_CPLEX;
        newbb.in_igma2_gap_min_solver = MRQ_IPOPT;
        
        newbb.in_use_outer_app = false;
        newbb.in_use_outer_app_as_heuristic = false;
        
        newbb.run(prob, &milpParams, &nlpParams);
        alg = &newbb;
        printf("____________________________________________________________________________\n");
        printf("New BB - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu lower bound: %0.2f iter to first feas sol: %lu iter to last feas sol: %lu cpu time to first feas sol: %0.2f clock time to first feas sol: %0.2f \n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, alg->out_lower_bound, alg->out_number_of_iterations_to_first_feas_sol, alg->out_number_of_iterations_to_best_sol, alg->out_cpu_time_to_first_feas_sol, alg->out_clock_time_to_first_feas_sol );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        std::cout << "constraint branchings: " << newbb.out_number_of_constraint_branchings  << " out_prune_counter:: bound: " << newbb.out_prune_counter.bound << " infeas: " << newbb.out_prune_counter.infeas << " opt: " << newbb.out_prune_counter.opt << " user: " << newbb.out_prune_counter.user << " bad lower bounds: " << newbb.out_nlp_failure_in_some_relaxation << " before solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning1 << " after solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning2 << "\n";
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%u%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, newbb.out_number_of_iters_having_wrong_lower_bounds, MRQ_CHAR_SEP, alg->out_number_of_iterations_to_best_sol, MRQ_CHAR_SEP, newbb.out_number_of_constraint_branchings, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning1, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning2, MRQ_CHAR_SEP, newbb.out_prune_counter.bound, MRQ_CHAR_SEP, newbb.out_prune_counter.infeas, MRQ_CHAR_SEP, newbb.out_prune_counter.opt, MRQ_CHAR_SEP, newbb.out_number_of_strong_branching_calculations_to_pseudo_costs, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
    }
    #endif

    #if 0  //branch-and-bound: to test constraint branching
    {
        newbb.in_number_of_threads = 1;
        newbb.in_print_level = 1;
        newbb.in_printing_frequency = 20000;
        
        newbb.in_use_outer_app = false;
        newbb.in_use_outer_app_as_heuristic = false;
        
        newbb.in_max_cpu_time = 4*60*60; //24 hours!!!
        
        newbb.in_nlp_solver = MRQ_IPOPT; //MRQ_NLP_MOSEK;
        newbb.in_milp_solver = MRQ_CPLEX;
        newbb.in_igma2_gap_min_solver = MRQ_IPOPT;
        
        newbb.in_igma2_strategy = MRQ_BB_IHS_NO_HEURISTICS;
        
        newbb.in_constr_branching_strategy = MRQ_BB_CBS_STBRANCH_PSEUDO_COSTS;
        
        newbb.in_max_number_of_branchings_in_constraint_branching = -1;
        newbb.in_count_total_prunes = true;
        
        //newbb.in_call_end_of_iteration_callback = true;
        
        //prob.print();
        //MRQ_getchar();
        
        newbb.run(prob, &milpParams, &nlpParams);
        
        alg = &newbb;
        printf("____________________________________________________________________________\n");
        printf("New BB - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu lower bound: %0.2f iter to first feas sol: %lu iter to last feas sol: %lu cpu time to first feas sol: %0.2f clock time to first feas sol: %0.2f \n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, alg->out_lower_bound, alg->out_number_of_iterations_to_first_feas_sol, alg->out_number_of_iterations_to_best_sol, alg->out_cpu_time_to_first_feas_sol, alg->out_clock_time_to_first_feas_sol );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        std::cout << "constraint branchings: " << newbb.out_number_of_constraint_branchings  << " out_prune_counter:: bound: " << newbb.out_prune_counter.bound << " infeas: " << newbb.out_prune_counter.infeas << " opt: " << newbb.out_prune_counter.opt << " user: " << newbb.out_prune_counter.user << " bad lower bounds: " << newbb.out_nlp_failure_in_some_relaxation << " before solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning1 << " after solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning2 << "\n";
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%u%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, newbb.out_number_of_iters_having_wrong_lower_bounds, MRQ_CHAR_SEP, alg->out_number_of_iterations_to_best_sol, MRQ_CHAR_SEP, newbb.out_number_of_constraint_branchings, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning1, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning2, MRQ_CHAR_SEP, newbb.out_prune_counter.bound, MRQ_CHAR_SEP, newbb.out_prune_counter.infeas, MRQ_CHAR_SEP, newbb.out_prune_counter.opt, MRQ_CHAR_SEP, newbb.out_number_of_strong_branching_calculations_to_pseudo_costs, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        if( newbb.out_number_of_constraint_branchings > 0 )
        {
            newbb.in_max_number_of_branchings_in_constraint_branching = 2;
            
            newbb.run(prob, &milpParams, &nlpParams);
        
            alg = &newbb;
            printf("____________________________________________________________________________\n");
            printf("New BB - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu lower bound: %0.2f iter to first feas sol: %lu iter to last feas sol: %lu cpu time to first feas sol: %0.2f clock time to first feas sol: %0.2f \n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, alg->out_lower_bound, alg->out_number_of_iterations_to_first_feas_sol, alg->out_number_of_iterations_to_best_sol, alg->out_cpu_time_to_first_feas_sol, alg->out_clock_time_to_first_feas_sol );
            //for(int i = 0; i < prob.n; i++)
                //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
            std::cout << "constraint branchings: " << newbb.out_number_of_constraint_branchings  << " out_prune_counter:: bound: " << newbb.out_prune_counter.bound << " infeas: " << newbb.out_prune_counter.infeas << " opt: " << newbb.out_prune_counter.opt << " user: " << newbb.out_prune_counter.user << " bad lower bounds: " << newbb.out_nlp_failure_in_some_relaxation << " before solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning1 << " after solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning2 << "\n";
            printf("____________________________________________________________________________\n");
            
            #if MRQ_SAVE_OUTPUT_FILE
                if(outFile)
                    fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%u%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, newbb.out_number_of_iters_having_wrong_lower_bounds, MRQ_CHAR_SEP, alg->out_number_of_iterations_to_best_sol, MRQ_CHAR_SEP, newbb.out_number_of_constraint_branchings, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning1, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning2, MRQ_CHAR_SEP, newbb.out_prune_counter.bound, MRQ_CHAR_SEP, newbb.out_prune_counter.infeas, MRQ_CHAR_SEP, newbb.out_prune_counter.opt, MRQ_CHAR_SEP, newbb.out_number_of_strong_branching_calculations_to_pseudo_costs, MRQ_CHAR_SEP);
                
                if(outFile)	    fflush(outFile);
            #endif
            
            
            
            
            newbb.in_constr_branching_strategy = MRQ_BB_CBS_NO_CONSTRAINT_BRANCH;
            
            newbb.run(prob, &milpParams, &nlpParams);
        
            alg = &newbb;
            printf("____________________________________________________________________________\n");
            printf("New BB - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu lower bound: %0.2f iter to first feas sol: %lu iter to last feas sol: %lu cpu time to first feas sol: %0.2f clock time to first feas sol: %0.2f \n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, alg->out_lower_bound, alg->out_number_of_iterations_to_first_feas_sol, alg->out_number_of_iterations_to_best_sol, alg->out_cpu_time_to_first_feas_sol, alg->out_clock_time_to_first_feas_sol );
            //for(int i = 0; i < prob.n; i++)
                //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
            std::cout << "constraint branchings: " << newbb.out_number_of_constraint_branchings  << " out_prune_counter:: bound: " << newbb.out_prune_counter.bound << " infeas: " << newbb.out_prune_counter.infeas << " opt: " << newbb.out_prune_counter.opt << " user: " << newbb.out_prune_counter.user << " bad lower bounds: " << newbb.out_nlp_failure_in_some_relaxation << " before solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning1 << " after solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning2 << "\n";
            printf("____________________________________________________________________________\n");
            
            #if MRQ_SAVE_OUTPUT_FILE
                if(outFile)
                    fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%u%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, newbb.out_number_of_iters_having_wrong_lower_bounds, MRQ_CHAR_SEP, alg->out_number_of_iterations_to_best_sol, MRQ_CHAR_SEP, newbb.out_number_of_constraint_branchings, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning1, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning2, MRQ_CHAR_SEP, newbb.out_prune_counter.bound, MRQ_CHAR_SEP, newbb.out_prune_counter.infeas, MRQ_CHAR_SEP, newbb.out_prune_counter.opt, MRQ_CHAR_SEP, newbb.out_number_of_strong_branching_calculations_to_pseudo_costs, MRQ_CHAR_SEP);
                
                if(outFile)	    fflush(outFile);
            #endif
        }
        
    }
    #endif


    #if 0  //branch-and-bound: to test pseudo-prune
    {
        //mosekParams.storeIntegerParameter("MSK_IPAR_INTPNT_MAX_ITERATIONS", 1000);
        
        //mosekParams.storeStringParameter("linear_solver", "ma27");
        //mosekParams.storeIntegerParameter("max_iter", 1000);
        
        
        //nlpParams.storeIntegerParameter("print_level", 4);
        //nlpParams.storeIntegerParameter("max_iter", 5000); //ipopt max number of iteration
        //nlpParams.storeStringParameter("derivative_test", "first-order");
        //nlpParams.storeStringParameter("derivative_test", "second-order");
        
        if( prob.getProblemType() == minlpproblem::MIP_PT_MINLP )
        { //in quadratic problems, Mosek gets bug to set this parameter
            nlpParams.storeIntegerParameter("MSK_IPAR_INTPNT_STARTING_POINT", MSK_STARTING_POINT_SATISFY_BOUNDS);
        }
        
        
        
        //MRQ_NewBB newbb;
        newbb.in_number_of_threads = 1;
        newbb.in_print_level = 4;
        newbb.in_printing_frequency = 100;
        //newbb.in_lower_bound = MRQ_INFINITY*0.1;
        
        newbb.in_use_outer_app = false;
        newbb.in_use_outer_app_as_heuristic = false;
        //newbb.in_outer_app_subprob_frequence = 4;
        
        //newbb.in_stop_multibranch_after_first_bound_prune = false;
        
        newbb.in_max_cpu_time = 24*60*60; //24 hours!!!
        newbb.in_max_time = 1.5 * 60 * 60;
        
        
        //newbb.in_igma2_strategy = MRQ_BB_IHS_ALWAYS;
        //newbb.in_igma2_frequency = 100000;
        
        //newbb.in_call_end_of_iteration_callback = true;
        
        //newbb.in_pseudo_cost_mu = 0.8;
        
        //newbb.in_use_early_branching = true;
        
        
        //newbb.in_constr_branching_strategy = MRQ_BB_CBS_NO_CONSTRAINT_BRANCH;// MRQ_BB_CBS_LOWEST_NUMBER_OF_VARS;
        newbb.in_nlp_solver = MRQ_ALGENCAN;
        newbb.in_milp_solver = MRQ_CPLEX;
        
        
        newbb.in_int_feas_heurs_strategy = MRQ_BB_IHS_ALWAYS;
        newbb.in_feas_heuristic_frequency = 1000;
        newbb.in_feas_heuristic_max_time = 60;
        newbb.in_use_feas_heuristic_diving = false;
        newbb.in_use_feas_heuristic_oafp = false;
        newbb.in_use_feas_heuristic_fp = false;
        
        newbb.in_rounding_heuristic_call_iter_frequence = 10;
        
        newbb.in_set_special_nlp_solver_params = false;
        
        newbb.in_count_total_prunes = true;
        newbb.in_pseudo_pruning_strategy = MRQ_BB_PPS_ON_NODE_EXPLORATION_AND_BRANCHING;
        newbb.in_min_number_of_bound_prunes_per_var_before_pseudo_pruning = 10;
        
        newbb.in_alpha_to_balance_estimative_in_pseudo_pruning = 0.2;
        newbb.in_relative_upper_bound_slack_factor_for_pseudo_pruning = 0.05; //change for 0.01
        newbb.in_calculate_pseudo_cost_average_above_error_estimative = true;
        
        newbb.in_relative_convergence_tol_for_pseudo_pruning = 0.1;
        
        
        newbb.run(prob, &milpParams, &nlpParams);
        
        alg = &newbb;
        printf("____________________________________________________________________________\n");
        printf("New BB - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu lower bound: %0.2f iter to first feas sol: %lu iter to last feas sol: %lu cpu time to first feas sol: %0.2f clock time to first feas sol: %0.2f pseudo cost standard deviation estimative: %0.2f \n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, alg->out_lower_bound, alg->out_number_of_iterations_to_first_feas_sol, alg->out_number_of_iterations_to_best_sol, alg->out_cpu_time_to_first_feas_sol, alg->out_clock_time_to_first_feas_sol, newbb.out_pseudo_cost_average_above_error_estimative );
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        std::cout << "out_prune_counter:: bound: " << newbb.out_prune_counter.bound << " infeas: " << newbb.out_prune_counter.infeas << " opt: " << newbb.out_prune_counter.opt << " user: " << newbb.out_prune_counter.user << " bad lower bounds: " << newbb.out_nlp_failure_in_some_relaxation << " before solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning1 << " after solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning2 << "\n";
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, newbb.out_number_of_iters_having_wrong_lower_bounds, MRQ_CHAR_SEP, alg->out_number_of_iterations_to_best_sol, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning1, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning2, MRQ_CHAR_SEP, newbb.out_prune_counter.bound, MRQ_CHAR_SEP, newbb.out_prune_counter.infeas, MRQ_CHAR_SEP, newbb.out_prune_counter.opt, MRQ_CHAR_SEP, newbb.out_pseudo_cost_average_above_error_estimative, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
            
        #if 0
        //if( newbb.out_return_code == MRQ_OPTIMAL_SOLUTION || newbb.out_return_code == MRQ_MAX_TIME_STOP )
        {
            
            //if( newbb.out_prune_counter.pseudopruning1 + newbb.out_prune_counter.pseudopruning2 > 0 )
            {
                newbb.in_pseudo_pruning_strategy = MRQ_BB_PPS_NO_PSEUDO_PRUNING;
                
                newbb.run(prob, &milpParams, &nlpParams);
                
                alg = &newbb;
                printf("____________________________________________________________________________\n");
                printf("New BB - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu lower bound: %0.2f iter to first feas sol: %lu iter to last feas sol: %lu cpu time to first feas sol: %0.2f clock time to first feas sol: %0.2f", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, alg->out_lower_bound, alg->out_number_of_iterations_to_first_feas_sol, alg->out_number_of_iterations_to_best_sol, alg->out_cpu_time_to_first_feas_sol, alg->out_clock_time_to_first_feas_sol );
                //for(int i = 0; i < prob.n; i++)
                    //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
                std::cout << "out_prune_counter:: bound: " << newbb.out_prune_counter.bound << " infeas: " << newbb.out_prune_counter.infeas << " opt: " << newbb.out_prune_counter.opt << " user: " << newbb.out_prune_counter.user << " bad lower bounds: " << newbb.out_nlp_failure_in_some_relaxation << " before solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning1 << " after solving pseudo prunes: " << newbb.out_prune_counter.pseudopruning2 << "\n";
                printf("____________________________________________________________________________\n");
                
                #if MRQ_SAVE_OUTPUT_FILE
                    if(outFile)
                        fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%f%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, newbb.out_number_of_iters_having_wrong_lower_bounds, MRQ_CHAR_SEP, alg->out_number_of_iterations_to_best_sol, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning1, MRQ_CHAR_SEP, newbb.out_prune_counter.pseudopruning2, MRQ_CHAR_SEP, newbb.out_prune_counter.bound, MRQ_CHAR_SEP, newbb.out_prune_counter.infeas, MRQ_CHAR_SEP, newbb.out_prune_counter.opt, MRQ_CHAR_SEP, newbb.out_pseudo_cost_average_above_error_estimative, MRQ_CHAR_SEP);
                    
                    if(outFile)	    fflush(outFile);
                #endif
            }
            
        }
    
        #endif
    }
    #endif


    #if 0  //igma0
    {
        nlpParams.storeIntegerParameter("MSK_IPAR_INTPNT_MAX_ITERATIONS", 1000);
        
        igma0.in_preprocess_lin_constr = false;
        igma0.in_preprocess_obj_function = false;
        igma0.in_preprocess_quad_constrs = false;
        
        igma0.in_max_cpu_time = 4*60*60;
        
        igma0.in_use_heuristcs = true;
        igma0.in_nlp_solver = MRQ_NLP_MOSEK;
        
        //igma1.in_use_integers_vars_on_gap_min_problem = true;
        
        igma0.in_relative_obj_cut_epsilon = 2.0e-3;
        igma0.in_number_of_threads = 1;
        
        igma0.run(prob, NULL, &nlpParams);
        
        alg = &igma0;
        printf("____________________________________________________________________________\n");
        printf("IGMA 1  Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
    }
    #endif


    #if 0  //heuristics: fp, diving, oafp
    {
        fp.in_max_cpu_time = 600;
        fp.in_nlp_solver = MRQ_IPOPT;
        
        fp.in_number_of_threads = 1;
        
        fp.run(prob, NULL, &nlpParams);
        
        alg = &fp;
        printf("____________________________________________________________________________\n");
        printf("FP Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif 
        
        
        dive.in_max_cpu_time = 600;
        dive.in_nlp_solver = MRQ_IPOPT;
        dive.in_dive_selec_strategy = MRQ_DIVE_SS_FRACTIONAL;
        
        dive.in_number_of_threads = 1;
        
        dive.run(prob, NULL, &nlpParams);
        
        alg = &dive;
        printf("____________________________________________________________________________\n");
        printf("Diving Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif 
        
        
        
        oafp.in_max_cpu_time = 600;
        oafp.in_nlp_solver = MRQ_IPOPT;
        oafp.in_milp_solver = MRQ_CPLEX;
        oafp.in_number_of_threads = 1;
        
        oafp.run(prob, &milpParams, &nlpParams);
        
        alg = &oafp;
        printf("____________________________________________________________________________\n");
        printf("OAFP Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
    }
    #endif


    #if 0  //igma1
    {
        
        
        MRQ_GeneralSolverParams ipoptParams;
        
        //ipoptParams.storeIntegerParameter("max_iter", 6000);
        //ipoptParams.storeIntegerParameter("MaxIter", 1500);
        
        //ipoptParams.storeStringParameter("derivative_test", "first-order");
        //ipoptParams.storeStringParameter("derivative_test", "only-second-order");
        //ipoptParams.storeIntegerParameter("print_level", 4);
        
        //ipoptParams.storeStringParameter("evaluate_orig_obj_at_resto trial", "no");
        
        
        //ipoptParams.storeStringParameter("linear_solver", "pardiso");
        //ipoptParams.storeStringParameter("linear_solver", "mumps");
        //ipoptParams.storeStringParameter("linear_solver", "ma57");
        //ipoptParams.storeDoubleParameter("max_cpu_time", 600);
        //mosekParams.storeIntegerParameter("print_level", 6);
        //ipoptParams.storeStringParameter("evaluate_orig_obj_at_rest_trial", "no");
        
        /*ipoptParams.storeDoubleParameter("acceptable_tol", 1e-4);
        ipoptParams.storeIntegerParameter("acceptable_iter", 10);
        ipoptParams.storeDoubleParameter("tol", 1e-5);
        ipoptParams.storeDoubleParameter("compl_inf_tol", 1e-3);*/
        
        //ipoptParams.storeDoubleParameter("tol", 1e-2);
        
        //ipoptParams.storeIntegerParameter("tuner", 1);
        
        //ipoptParams.storeIntegerParameter("derivcheck", KTR_DERIVCHECK_ALL);
        //ipoptParams.storeIntegerParameter("outlev", 3);
        //ipoptParams.storeDoubleParameter("derivcheck_tol", 1e-4);
        
        
        //mosekParams.storeIntegerParameter("MSK_IPAR_CONCURRENT_NUM_OPTIMIZERS", 2);
        
        
        igma1.in_lower_bound_to_random_sol = -100;
        igma1.in_upper_bound_to_random_sol = 100;
        
        
        igma1.in_number_of_threads = 1;
        igma1.in_use_general_feas_heuristics = false;
        igma1.in_heuristic_mode = true; //set it as true get solution quality poor
        igma1.in_adopt_obj_cut = true;
        
        //igma1.in_enable_gap_min_solver_premature_stoping = true;
        igma1.in_set_special_gap_min_solver_params = true;
        igma1.in_consider_relax_infeas_if_solver_fail = true;
        
        igma1.in_printing_frequency = 100;
        igma1.in_print_level = 4;
        
        
        igma1.in_gap_min_obj_strategy = MRQ_IGMA_GMOS_BY_GAP_AVERAGE;
        
        
        igma1.in_constr_branching_strategy = MRQ_BB_CBS_NO_CONSTRAINT_BRANCH;
        
        igma1.in_max_cpu_time = 600;
        
        igma1.in_integer_tol = 1.0e-3;
        
        igma1.in_gap_min_solver = MRQ_IPOPT;
        //igma1.in_gap_min_solver = MRQ_NLP_KNITRO;
        
        
        igma1.in_nlp_solver = igma1.in_gap_min_solver;
        
        //igma1.in_min_number_of_iters_on_gap_min_premature_stop = 1000000000;
        
        
        
        igma1.run( prob, &ipoptParams, &ipoptParams );
        alg = &igma1;
        printf("____________________________________________________________________________\n");
        printf("IGMA 1 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f subiters: %lu iters: %lu prunes by obj active: %lu\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, igma1.out_number_of_inner_iterations, alg->out_number_of_iterations, igma1.out_number_of_prunes_by_obj_cut_active);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %0.12f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");
        
        /*{
            int r;
            bool feas;
            
            r = prob.isFeasibleToConstraints(0, alg->out_best_sol, true, NULL, 1e-6, 1e-6, feas);
            
            std::cout << "Testei viabilidade da solucao:: r: " << r << " feas: " << feas << "\n";
        }*/
        
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s" "%lu%s",
                alg->out_algorithm, MRQ_CHAR_SEP, 
                alg->out_return_code, MRQ_CHAR_SEP,
                alg->out_lower_bound, MRQ_CHAR_SEP,
                alg->out_best_obj, MRQ_CHAR_SEP,
                alg->out_cpu_time, MRQ_CHAR_SEP,
                alg->out_clock_time, MRQ_CHAR_SEP,
                igma1.out_number_of_inner_iterations, MRQ_CHAR_SEP,
                alg->out_number_of_iterations, MRQ_CHAR_SEP,
                igma1.out_number_of_iterations_to_best_sol, MRQ_CHAR_SEP,
                igma1.out_number_of_feas_sols, MRQ_CHAR_SEP,
                igma1.out_number_of_prunes_by_obj_cut_active, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
    }
    #endif


    #if 0  //igma2
    {
        //MRQ_GeneralSolverParams knitroParams;
        MRQ_GeneralSolverParams ipoptParams;
        
        //ipoptParams.storeStringParameter("linear_solver", "ma57");
        //ipoptParams.storeDoubleParameter("max_cpu_time", 200);
        ipoptParams.storeIntegerParameter("max_iter", 6000);
        
        //knitroParams.storeIntegerParameter("algorithm", 4);
        //knitroParams.storeIntegerParameter("act_lpsolver", 2); //KTR_ACT_LPSOLVER_CPLEX
        //knitroParams.storeIntegerParameter("blasoption", 1);
        //knitroParams.storeIntegerParameter("bar_switchrule", 1); //KTR_BAR_SWITCHRULE_NEVER
        
        //ipoptParams.storeStringParameter("expect_infeasible_problem", "no"); //igma3 does not solve infeasible gap min problems
        
        
        //igma2.in_print_level = 10;
        //igma2.in_printing_frequency = 1;
        igma2.in_number_of_threads = 1;
        igma2.in_nlp_solver =  MRQ_IPOPT;
        igma2.in_max_cpu_time = 600;
        
        igma2.in_solve_local_search_problem_even_on_non_int_sol = true;
        
        igma2.in_factor_to_max_dist_constr_on_bin_vars = 0.2;
        igma2.in_set_max_dist_constr_on_bin_vars = false;
        igma2.in_number_of_threads = 1;
        
        //igma2.in_neighborhood_strategy = MRQ_IGMA2_NS_RECTANGULAR;
        //igma2.in_percentual_to_rectangular_neighborhood = 0.025;
        
        igma2.run( prob, &ipoptParams, &ipoptParams );
        alg = &igma2;
        printf("____________________________________________________________________________\n");
        printf("IGMA2 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu feas sol on gap min: %d feas sol by igma2: %d\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, (int) igma2.out_feas_sol_on_gap_min_problem, (int) igma2.out_feas_sol_found_by_igma2_procedure);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%d%s" "%d%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, (int) igma2.out_feas_sol_on_gap_min_problem, MRQ_CHAR_SEP, (int) igma2.out_feas_sol_found_by_igma2_procedure, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
            
        /*igma2.in_nlp_solver = MRQ_ALGENCAN;
        
        igma2.run( prob, nullptr, nullptr );
        alg = &igma2;
        printf("____________________________________________________________________________\n");
        printf("IGMA2 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu feas sol on gap min: %d feas sol by igma2: %d\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, (int) igma2.out_feas_sol_on_gap_min_problem, (int) igma2.out_feas_sol_found_by_igma2_procedure);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%d%s" "%d%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, (int) igma2.out_feas_sol_on_gap_min_problem, MRQ_CHAR_SEP, (int) igma2.out_feas_sol_found_by_igma2_procedure, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
            
            
            
        igma2.in_nlp_solver = MRQ_WORHP;
        
        igma2.run( prob, nullptr, nullptr );
        alg = &igma2;
        printf("____________________________________________________________________________\n");
        printf("IGMA2 Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu feas sol on gap min: %d feas sol by igma2: %d\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, (int) igma2.out_feas_sol_on_gap_min_problem, (int) igma2.out_feas_sol_found_by_igma2_procedure);
        //for(i = 0; i < n; i++)
            //printf("x[%d]: %f\n", i, ecp.out_best_sol[i]);
        printf("____________________________________________________________________________\n");

        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%d%s" "%d%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, (int) igma2.out_feas_sol_on_gap_min_problem, MRQ_CHAR_SEP, (int) igma2.out_feas_sol_found_by_igma2_procedure, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif */
        
    }
    #endif
    
    
    
    #if 1//Structured Stochastic Rounding - testing MIP_BinSumConstrsIndsByClass
    {
        minlpproblem::MIP_BinSumConstrsIndsByClass binSumInds;
        
        //prob.print();
        
        #if 0
        {
            const int n = prob.n;
            minlpproblem::MIP_ConstraintsByColumnsStorager ccstorager;
            
            #if 1
            prob.print();
            
            int r = ccstorager. storageConstraintsByColumns( prob, true );
            
            {  
                std::cout << "constraints per variable:\n";
                
                auto *constrColsOffset = ccstorager.constrColsOffset;
                auto *constrCols = ccstorager.constrCols;
                
                for(unsigned int i = 0; i < n; i++)
                {
                    const unsigned int nzs = constrColsOffset[i+1] - constrColsOffset[i];
                    
                    auto constrIndices = &constrCols[ constrColsOffset[i] ];
                    
                    std::cout << i << ":";
                    
                    for(unsigned int j = 0; j < nzs; j++)
                    {
                        std::cout << " " << constrIndices[j];
                    }
                    
                    std::cout << "\n";
                }
                
                MRQ_getchar();
            }
            #endif
        }
        #endif
        
        int r = binSumInds.calculateIndices(prob, nullptr, nullptr, nullptr, nullptr, nullptr, true, false, true);
        if(r)
            MRQ_getchar();
        
        #if 0
            #if MRQ_SAVE_OUTPUT_FILE
                if(outFile)
                    fprintf(outFile, "%s%s", "CLASS STATICS", MRQ_CHAR_SEP );
            #endif
            
            for( unsigned int i = 0; i < binSumInds.nBinSumConstrClasses ; i++ )
            {
                
                #if MRQ_SAVE_OUTPUT_FILE
                    if(outFile)
                        fprintf(outFile, "%d%s", binSumInds.nClasses[i], MRQ_CHAR_SEP ); //saving the number of constraints in each classes
                #endif
                
                #if 0
                std::cout << "class " << i << " (" << binSumInds.nClasses[i] << " constraints) : ";
                
                
                for( unsigned int j = 0; j < binSumInds.nClasses[i]; j++ )
                {
                    std::cout << binSumInds.classes[i][j] << " ";
                }
                std::cout << "\n";
                #endif
            } 
            #if MRQ_SAVE_OUTPUT_FILE  
                if(outFile)	    fflush(outFile);
            #endif
        #endif
        
        
        
        /*if(binSumInds.nKnapsackInds)
        {
            int intVars[prob.n];
            int nI;
            
            nI = prob.getIntegerIndices(intVars); 
            
            for( unsigned int i = 0; i < prob.getNumberOfIntegerVars(); i++ )
            {
                std::cout << intVars[i] << ": ";
                
                for(unsigned j = 0; j < binSumInds.nKnapsackInds[i]; j++)
                {
                    std::cout << binSumInds.knapsackInds[i][j] << " ";
                }
                std::cout << "\n";
            }
        } */
        
        /*ssr.in_preprocess_lin_constr = false;
        ssr.in_preprocess_obj_function = false;
        ssr.in_preprocess_quad_constrs = false;*/
        
        
        
        ssr.in_print_level = 3;
        ssr.in_max_iterations = 10000;
        ssr.in_solve_minlp_as_local_search = true;
        ssr.in_preprocess_after_variable_rounding = true;
        
        //ssr.in_additional_vars_fixing_strategy = MRQ_SSRVAFS_LINEAR_AUXILIAR_PROBLEM;
        
        ssr.in_neighborhood_strategy = MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD;
        ssr.in_stop_local_search_solving_on_first_improvment_solution = false;
        ssr.in_max_number_of_improvments = 1;
        
        ssr.in_milp_solver = MRQ_CPLEX;
        ssr.in_nlp_solver = MRQ_IPOPT;
        
        ssr.in_relative_convergence_tol_to_local_search = 0.3;
        
        ssr.in_integer_neighborhood_factor = 0.3;
        
        ssr.in_max_cpu_time = 600;
        ssr.in_number_of_threads = 1;
        
        ssr.in_min_probability_to_round = 0.0;
        
        ssr.in_solve_continuous_relaxation = true;
        //ssr.in_rounding_var_bounds_updt_strategy = MRQ_SSRVBUS_AUXILIAR_PROBLEM;
        
        ssr.in_solve_continuous_relaxation = true;
        ssr.in_cont_relax_strategy_to_stoch_rounding = MRQ_SSR_CRSSR_NO_CONTINUOUS_RELAXATION;
        
        
        //prob.print();
        
        ssr.run(prob, &milpParams, &nlpParams);
        
        alg = &ssr;
        printf("____________________________________________________________________________\n");
        printf("SSR Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu obj at nlp fixing: %0.10f cpu time after nlp fixing: %0.2f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, ssr.out_obj_at_nlp_integer_fixed_sol, ssr.out_cpu_time_at_nlp_integer_fixed_sol);
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%d%s" "%lf%s" "%lf%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, (int) ssr.out_vars_fixed_by_stoch_rounding, MRQ_CHAR_SEP, ssr.out_obj_at_first_sol, MRQ_CHAR_SEP, ssr.out_cpu_time_at_nlp_integer_fixed_sol, MRQ_CHAR_SEP );
            
            if(outFile)	    fflush(outFile);
        #endif
        
        /*ssr.in_solve_continuous_relaxation = false;
        ssr.in_cont_relax_strategy_to_stoch_rounding = MRQ_SSR_CRSSR_NO_CONTINUOUS_RELAXATION;
        
        ssr.run(prob, &milpParams, &nlpParams);
        
        alg = &ssr;
        printf("____________________________________________________________________________\n");
        printf("SSR Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu obj at nlp fixing: %0.10f cpu time after nlp fixing: %0.2f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, ssr.out_obj_at_nlp_integer_fixed_sol, ssr.out_cpu_time_at_nlp_integer_fixed_sol);
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%d%s" "%lf%s" "%lf%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, (int) ssr.out_vars_fixed_by_stoch_rounding, MRQ_CHAR_SEP, ssr.out_obj_at_first_sol, MRQ_CHAR_SEP, ssr.out_cpu_time_at_nlp_integer_fixed_sol, MRQ_CHAR_SEP );
            
            if(outFile)	    fflush(outFile);
        #endif */
        
            
        /*    
        //ssr.in_random_order_to_threat_classes = true;
        
        //ssr.in_relative_convergence_tol_to_local_search = 0.5;
        
        ssr.in_integer_neighborhood_factor = 0.1; 
        
        ssr.run(prob, &milpParams, &nlpParams);
        
        alg = &ssr;
        printf("____________________________________________________________________________\n");
        printf("SSR Problem: %s Algorithm: %d Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu obj at nlp fixing: %0.10f cpu time after nlp fixing: %0.2f\n", probName, alg->out_algorithm, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations, ssr.out_obj_at_nlp_integer_fixed_sol, ssr.out_cpu_time_at_nlp_integer_fixed_sol);
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s" "%d%s" "%lf%s" "%lf%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP, (int) ssr.out_vars_fixed_by_stoch_rounding, MRQ_CHAR_SEP, ssr.out_obj_at_nlp_integer_fixed_sol, MRQ_CHAR_SEP, ssr.out_cpu_time_at_nlp_integer_fixed_sol, MRQ_CHAR_SEP );
            
            if(outFile)	    fflush(outFile);
        #endif */
        
    }
    #endif
    

    

    #if 0 //rens
    {
        MRQ_GeneralSolverParams nlpParams;
        MRQ_GeneralSolverParams minlpParams;
        
        //nlpParams.storeStringParameter("derivative_test", "first-order");
        //nlpParams.storeIntegerParameter("print_level", 4);
        
        
        rens.in_nlp_solver = MRQ_NLP_MOSEK;
        rens.in_milp_solver = MRQ_CPLEX;
        
        rens.in_number_of_threads = 1;
        rens.in_algorithm_to_solve_subproblem = MRQ_OA_ALG; //MRQ_LP_BB_ECP_BASED_ALG;
        
        rens.in_algorithm_object_to_solve_subproblem = &oa;
        
        
        rens.in_max_cpu_time= 120;
        
        rens.in_integer_neighborhood_factor = 0.40;
        rens.in_continuous_neighborhood_factor = 0.40;
        rens.in_neighborhood_strategy = MRQ_RENS_NS_ORIGINAL;
        
        rens.in_apply_only_heuristics_on_subproblem = false;
        rens.in_stop_on_first_sol = true;
        //rens.in_print_level = 10;
        
        rens.run(prob, &milpParams, &nlpParams, &minlpParams);
        
        alg = &rens;
        printf("____________________________________________________________________________\n");
        printf("RENS original - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu\n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        rens.in_integer_neighborhood_factor = 0.4;
        rens.in_neighborhood_strategy = MRQ_RENS_NS_EUCLIDEAN_INTEGER_NEIGHBORHOOD;
        rens.run(prob, &milpParams, &nlpParams, &minlpParams);
        
        alg = &rens;
        printf("____________________________________________________________________________\n");
        printf("RENS 0.4 - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu\n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        //rens.in_neighborhood_strategy = MRQ_RENS_NS_EUCLIDEAN_NEIGHBORHOOD;
        //rens.in_neighborhood_strategy = MRQ_RENS_NS_EUCLIDEAN_NEIGHBORHOOD;
        rens.in_integer_neighborhood_factor = 0.3;
        rens.run(prob, &milpParams, &nlpParams, &minlpParams);
        
        alg = &rens;
        printf("____________________________________________________________________________\n");
        printf("RENS 0.3 - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu\n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        
        //rens.in_algorithm_to_solve_subproblem =  MRQ_FP_HEUR_ALG;
        //rens.in_neighborhood_strategy = MRQ_RENS_NS_EUCLIDEAN_INTEGER_NEIGHBORHOOD;
        
        rens.in_integer_neighborhood_factor = 0.20;
        
        rens.run(prob, &nlpParams, &minlpParams);
        
        alg = &rens;
        printf("____________________________________________________________________________\n");
        printf("RENS 0.2 - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu\n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        
        rens.in_integer_neighborhood_factor = 0.1;
        
        rens.run(prob, &nlpParams, &minlpParams);
        
        alg = &rens;
        printf("____________________________________________________________________________\n");
        printf("RENS 0.1 - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu\n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif
        
        
        /*rens.in_integer_neighborhood_factor = 0.05;
        
        rens.run(prob, &nlpParams, &minlpParams);
        
        alg = &rens;
        printf("____________________________________________________________________________\n");
        printf("RENS 0.05 - Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %lu\n", probName, alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        //for(int i = 0; i < prob.n; i++)
            //printf("x[%d]: %f\n", i, alg->out_best_sol[i]);
        printf("____________________________________________________________________________\n");
        
        #if MRQ_SAVE_OUTPUT_FILE
            if(outFile)
                fprintf(outFile, "%d%s" "%d%s" "%f%s" "%f%s" "%f%s" "%f%s" "%lu%s", alg->out_algorithm, MRQ_CHAR_SEP, alg->out_return_code, MRQ_CHAR_SEP, alg->out_lower_bound, MRQ_CHAR_SEP, alg->out_best_obj, MRQ_CHAR_SEP, alg->out_cpu_time, MRQ_CHAR_SEP, alg->out_clock_time, MRQ_CHAR_SEP, alg->out_number_of_iterations, MRQ_CHAR_SEP);
            
            if(outFile)	    fflush(outFile);
        #endif */
        
    }
    #endif


    #if 0 //testing heuristic executor
    {
        MRQ_HeuristicExecutor heurs;
        double obj;
        
        heurs.in_stop_on_first_sol = false;
        
        heurs.run(prob, &milpParams, &nlpParams, INFINITY, obj);
        
        printf("____________________________________________________________________________\n");
        std::cout << "Heuristic executor - obj: " << obj << "\n" ;
        printf("____________________________________________________________________________\n");
    }
    #endif


    #if 0 //here, we print oa and continuous relaxation solution
    {
        char cRelaxSolFileName[200]; 
        char oaSolFileName[200];
        const int n = prob.n;
        
        std::cout << "stub: " << probName << "\n";
        
        sprintf(cRelaxSolFileName, "%s.crelax.solution", probName);
        
        sprintf(oaSolFileName, "%s.oa.solution", probName);
        
        if( crelax.out_return_code == MRQ_CONT_RELAX_OPTIMAL_SOLUTION
        )
        {
            FILE *file = fopen( cRelaxSolFileName, "w" );
            
            if(!file)
            {
                MRQ_PRINTERRORMSG("Error to open solution file to continuous relaxation");
            }
            else
            {
                MRQ_writeSolOnFile(file, crelax.out_best_obj, n, crelax.out_best_sol);
                
                fclose(file);
            }
        }
        
        if(oa.out_return_code == MRQ_OPTIMAL_SOLUTION)
        {
            FILE *file = fopen( oaSolFileName, "w" );
            
            if(!file)
            {
                MRQ_PRINTERRORMSG("Error to open solution file to continuous relaxation");
            }
            else
            {
                MRQ_writeSolOnFile(file, oa.out_best_obj, n, oa.out_best_sol);
                
                fclose(file);
            }
        }
        
        
        {
            int nI;
            double rsol[n];
            int intVars[n];
            double *sol1 = crelax.out_best_sol;
            double *sol2 = oa.out_best_sol;
            
            
            nI = prob.getIntegerIndices(intVars);
            
            MRQ_copyArray(n, sol1, rsol);
            
            for(int i = 0; i < nI; i++)
            {
                const int ind = intVars[i];
                rsol[ind] = round( sol1[ind] );
            }
            
            double ed1 = MRQ_norm1Distance(n, sol1, sol2);
            double ed2 = MRQ_norm2Distance(n, sol1, sol2);
            
            double eid1 = MRQ_norm1Distance(n, rsol, sol2);
            double eid2 = MRQ_norm2Distance(n, rsol, sol2);
            
            double id = MRQ_integerDistance(nI, intVars, rsol, sol2);
            double id2 = MRQ_integer2Distance(nI, intVars, rsol, sol2);
            
            #if MRQ_SAVE_OUTPUT_FILE
                if(outFile)
                    fprintf(outFile, "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP "%f" MRQ_CHAR_SEP, ed1, ed2, eid1, eid2, id, id2);
            #endif
        }
        
    }
    #endif




    termination:

    if(outFile)
            fclose(outFile);

    return 0;
}





