#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

#include <new>

#include "MRQ_solvers.hpp"

#include "MRQ_algClasses.hpp"
//#include "MRQ_nlpSolvers.hpp"
#include "MRQ_tools.hpp"


using namespace optsolvers;
using namespace muriqui;

//using namespace std;




MRQ_OAFeasibilityPump::MRQ_OAFeasibilityPump():MRQ_LinearApproxAlgorithm()//, MRQ_Heuristic()
{
	resetParameters();
	resetOutput();
	out_algorithm = MRQ_OA_FP_HEUR_ALG;
}



MRQ_OAFeasibilityPump::~MRQ_OAFeasibilityPump(){ }


void MRQ_OAFeasibilityPump::printParameters( std::ostream &out) const
{
	MRQ_LinearApproxAlgorithm::printParameters(out);
	out << "\n"
	
	MRQ_STRFFATT(in_set_norm1_on_nlp) << "\n"
	MRQ_STRFFATT(in_set_linear_obj_term_on_bin_vars_at_nlp) << "\n"
	MRQ_STRFFATT(in_solve_nlp_as_local_search_at_end) << "\n";
}


void MRQ_OAFeasibilityPump::resetParameters()
{
	MRQ_LinearApproxAlgorithm::resetParameters();
	
	//in_max_number_of_gap_equals = 3;
	in_set_norm1_on_nlp = true;
	in_set_linear_obj_term_on_bin_vars_at_nlp = true;
	in_solve_nlp_as_local_search_at_end = true;
}



int MRQ_OAFeasibilityPump::setIntegerParameter( const char *name, const long int value)
{
	int ret = MRQ_LinearApproxAlgorithm::setIntegerParameter(name, value);
	
	if(ret == 0)
		return 0;
	
	ret = 0;
	
	
	if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_norm1_on_nlp), name, value ) == 0 );
	else if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_linear_obj_term_on_bin_vars_at_nlp), name, value ) == 0 );
	else if( MRQ_setAtt<bool>( MRQ_STRATT(in_solve_nlp_as_local_search_at_end), name, value ) == 0 );
	else
		ret = MRQ_NAME_ERROR;
	
	return ret;
}



int MRQ_OAFeasibilityPump::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
	const int n = prob.n;
	const int m = prob.m;
	const int nI= prob.getNumberOfIntegerVars();
	const bool preprocess = in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_obj_function;
	
	bool flag, updtConstrBounds;
	int i, r;
	unsigned long int iter = 0;
	double timeStart;
	
	clock_t clockStart;
	
	OPT_LPSolver *master = NULL;
	MRQ_MasterMILPProb masterMilp;
	
	
	MRQ_NLPSolver *nlp = NULL;
	//MRQ_NLPSolver *nlpoafp = NULL;
	MRQ_NLPFeasPumpProb nlpfp;
	MRQ_Preprocessor preprocessor(&prob);
	
	double *lx = run_by_inside ? nlx : prob.lx;
	double *ux = run_by_inside ? nux : prob.ux;
	
	bool *auxEval = NULL;
	int *intVars = NULL;
	double *auxVars = NULL, *auxConstr = NULL;
	double *auxp1, *auxp2 = NULL; //auxp2 is not allocated
	double *plc = NULL, *puc = NULL;
	
	timeStart = MRQ_getTime();
	clockStart = clock();
	
	//zl = in_lower_bound;
	//zu = in_upper_bound;
	
	{
		auto ret = algorithmInitialization(1, (in_preprocess_lin_constr || in_preprocess_obj_function || in_preprocess_quad_constrs), milpSolverParams, nlpSolverParams, prob, lx, ux, &preprocessor, &updtConstrBounds, &plc, &puc); //that algorithm is monothread...
		MRQ_IFERRORGOTOLABEL(ret, out_return_code, ret, termination);
	}
	
	if(in_print_level > 1)
		std::cout << "\n\n" << MRQ_PREPRINT << "Starting Outer Approximation Feasibility Pump heuristic" << "\n\n";
	
	if(in_print_level > 3)
		printSubSolvers(true, true, false);
	
	
	MRQ_malloc(auxEval, m); //auxEval = (bool *) malloc(m * sizeof(bool) );
	MRQ_malloc(auxConstr, m); //auxConstr = (double *) malloc(m * sizeof(double) );
	MRQ_malloc(auxVars, n); //auxVars = (double *) malloc( n * sizeof(double) );
	MRQ_malloc(intVars, nI); //intVars = (int *) malloc( nI * sizeof(int) );
	MRQ_IFMEMERRORGOTOLABEL(!auxEval || !auxConstr || !auxVars || !intVars, out_return_code, termination);
	
	for(i = 0; i < m; i++)
		auxEval[i] = prob.nlConstr[i] || prob.QC[i].getNumberOfElements() > 0;
	
	prob.getIntegerIndices(intVars);
	
	
	if( !in_use_initial_solution || in_solve_nlp_as_local_search_at_end)
	{
		nlp = OPT_newNLPSolver( in_nlp_solver );
		MRQ_IFMEMERRORGOTOLABEL(!nlp, out_return_code, termination);
		
		r = MRQ_setNLPRelaxProb(prob, lx, ux, plc, puc, nlp, true, true, true, false, thnumber, in_set_special_nlp_solver_params, nlpSolverParams, in_number_of_threads, in_max_cpu_time, in_max_time, 0, 0);
		MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
		
		if( run_by_inside )
		{
			if( !std::isinf(insideSolverMaxTime) )
			{
				r = nlp->setMaxTime(insideSolverMaxTime );
				MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
			}
		}
		
		
		if( !in_use_initial_solution )
		{
			r = nlp->solve(false);
			
			//printf("code1: %d obj1: %f\n", nlp->retCode, nlp->objValue);
			
			//for(int i = 0; i < n; i++)
				//printf("sol1[%d]: %f \n", i, nlp->sol[i] );
			
			
			if( r == OPT_OPTIMAL_SOLUTION )
			{
				out_obj_opt_at_continuous_relax = nlp->objValue;
				zl = MRQ_max( zl,  nlp->objValue);
				
				
				if(in_print_level > 1)
					std::cout << MRQ_PREPRINT "NLP relaxation solution: " << nlp->objValue << "\n";
				
				
				if( prob.isIntegerSolution(nlp->sol, in_integer_tol) )
				{
					if( in_print_level > 1 )
						std::cout << MRQ_PREPRINT "An integer optimal solution was gotten as NLP relaxation solution\n";
					
					//MRQ_copyArray( n, nlp->sol, out_best_sol);
					
					//out_best_obj = zu = nlp->objValue;
					
					const bool flag = tryUpdateBestSolution( thnumber, n, nlp->sol, nlp->objValue, iter, clockStart, timeStart, in_store_history_solutions );
					
					#if MRQ_DEBUG_MODE
						if(!run_by_inside)
							assert(flag);
					#endif
					
					out_return_code = MRQ_OPTIMAL_SOLUTION;
					
					goto termination;
				}
				
				//MRQ_copyArray(n, nlp->sol, auxVars);
				auxp1 = nlp->sol;
				auxp2 = nlp->constr;
			}
			else if( nlp->feasSol ) //( ret == MRQ_NLP_FEASIBLE_SOLUTION )
			{
				//MRQ_copyArray(n, nlp->sol, auxVars);
				auxp1 = nlp->sol;
				auxp2 = nlp->constr;
			}
			else if( r == OPT_INFEASIBLE_PROBLEM )
			{
				out_return_code = MRQ_INFEASIBLE_PROBLEM;
				goto termination;
			}
			else
			{
				//we need a first solution...
				
				#pragma ivdep
				#pragma GCC ivdep
				for(int i = 0; i < n; i++)
					auxVars[i] = MRQ_min( ux[i], MRQ_max(0.0, lx[i]) );
				
				auxp1 = auxVars;
			}
			
		}
		
	}
	
	
	if(in_use_initial_solution)
	{
		//MRQ_copyArray(n, xInit, auxVars);
		auxp1 = xInit;
	}
	
	
	r = masterMilp.setProblemBase(thnumber, prob, in_milp_solver, false, false, true, lx, ux, nI, milpSolverParams, in_number_of_threads);
	MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
	
	master = masterMilp.master;
	
	if( run_by_inside )
	{
		if( !std::isinf(insideSolverMaxTime) )
		{
			r = master->setMaxTime( insideSolverMaxTime );
			MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
		}
	}
	
	
	r = masterMilp.addNorm1constr(n, nI, intVars, auxp1);
	MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
	
	r = masterMilp.addLinearizedNLConstraints( true, auxp1, false, auxEval, NULL, NULL, auxp2, auxp2 != NULL);
	MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
	
	
	while(true)
	{
		iter++;
		
		if( in_max_cpu_time < INFINITY || in_max_time < INFINITY )
		{
			const double rcputime = in_max_cpu_time - ( ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC);
			const double rtime = in_max_time - (MRQ_getTime() - timeStart);
			int r = 0;
			
			
			if(rcputime <= 0.0 || rtime <= 0.0)
			{
				out_return_code = MRQ_MAX_TIME_STOP;
				break;
			}
			
			if(in_max_cpu_time < INFINITY)
				r = master->setMaxCPUTime(rcputime);
			
			if(in_max_time < INFINITY)
				r += master->setMaxTime(rtime);
			
			if(r != 0)
			{
				if(in_print_level > 0)
					MRQ_PRINTERRORMSG("Error to set MILP maximum time!");
			}
		}
		
		
		
		//master->generateModelFile("oamaster.lp");
		//std::cout << "gerei oamaster.lp" << std:: endl;
		//getchar();
		
		r = master->solve(true);
		
		/*printf("iret code: %d obj value: %f\n", master->retCode, master->objValue);
		for(int w  = 0; w < prob.n; w++)
			printf("isol[%d]: %0.3f \t", w, master->sol[w]);
		printf("\n\n"); */
		
		{  //couting number of iterations
			long unsigned int masterIters;
			int r = master->getNumberOfIterations(masterIters);
			if(r == 0)
				out_number_of_milp_solver_iters += masterIters;
			else if(in_print_level > 0)
				MRQ_PRINTERRORNUMBER(r);
		}
		
		if( r != OPT_OPTIMAL_SOLUTION && !master->feasSol )
		{
			if( r == OPT_INFEASIBLE_PROBLEM )
			{
				if(in_print_level > 1)
					MRQ_PRINTERRORMSG("MILP problem is infeasible!\n");
				
				out_return_code = MRQ_INFEASIBLE_PROBLEM;
				break;
			}
			else if( r == OPT_UNBOUNDED_PROBLEM )
			{
				if(in_print_level > 1)
					MRQ_PRINTERRORMSG("MILP problem is unbounded!\n");
				
				out_return_code = MRQ_UNBOUNDED_PROBLEM;
				break;
			}
			else
			{
				if(in_print_level > 1)
					MRQ_PRINTERRORMSG("Error at MILP solving!\n");
				
				out_return_code = MRQ_MILP_SOLVER_ERROR;
				break;
			}
		}
		
		
		//we do not need evaluate linear constraints because they are already enforced in master
		r = prob.isFeasibleToConstraints(thnumber, master->sol, true, auxEval, in_absolute_feasibility_tol, in_relative_feasibility_tol, flag, auxConstr);
		
		if(flag)
		{
			double obj;
			
			r = prob.objEval(thnumber, !prob.hasNlConstrs, master->sol, obj);
			
			if(r == 0)
			{
				const bool flag = tryUpdateBestSolution( thnumber, n, master->sol, obj, iter, clockStart, timeStart, in_store_history_solutions );
				
				#if MRQ_DEBUG_MODE
					if(!run_by_inside)
						assert(flag);
				#endif
				
				//std::cout << "encontrei a melhor solucao pelo master. flag: " << flag << " obj: "<< obj << " zu: " << zu <<  "\n";
				
				out_return_code = flag ? MRQ_HEURISTIC_SUCCESS : MRQ_HEURISTIC_FAIL;
				break;
			}
		}
		
		
		//solving the nlp fp problem
		if( nlpfp.solver == NULL )
		{
			r = nlpfp.setProblemBase( in_nlp_solver, prob, lx, ux, plc, puc, nlpSolverParams, thnumber, in_set_special_nlp_solver_params, in_set_linear_obj_term_on_bin_vars_at_nlp, in_set_norm1_on_nlp, in_number_of_threads, in_max_cpu_time, in_max_time );
			
			MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
			
			if( run_by_inside )
			{
				if( !std::isinf(insideSolverMaxTime) )
				{
					r = nlpfp.solver->setMaxTime( insideSolverMaxTime );
					MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
				}
			}
		}
		
		r = nlpfp.setObjective( nI, intVars, lx, ux, master->sol);
		MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
		
		
		r = nlpfp.solver->solve(false);
		
		/*printf( "cret code: %d obj value: %f \n", nlpfp.solver->retCode, nlpfp.solver->objValue );
		for(int w = 0; w < prob.n; w++)
			printf("csol[%d]: %0.5f \t", w, nlpfp.solver->sol[w]);
		printf("\n"); */
		
		if( r == OPT_OPTIMAL_SOLUTION || nlpfp.solver->feasSol )
		{
			if( MRQ_isIntegerSol(nI, intVars, nlpfp.solver->sol, in_integer_tol) )
			{
				double obj;
				
				r = prob.objEval(thnumber, true, nlpfp.solver->sol, obj);
				
				if(r == 0)
				{
					const bool r = tryUpdateBestSolution( thnumber, n, nlpfp.solver->sol, obj, iter, clockStart, timeStart, in_store_history_solutions );
					
					out_return_code = r ? MRQ_HEURISTIC_SUCCESS : MRQ_HEURISTIC_FAIL;
					break;
				}
				//TODO: if ret != 0, probably we must declare heuristic failure because we could not find other solution, but by now, I will let in this way...
			}
			
		}
		else
		{
			if(in_print_level > 1)
				std::cerr << MRQ_PREPRINT << "Error " << r << " at NLP solving!\n";
			
			out_return_code = MRQ_NLP_SOLVER_ERROR;
			break;
		}
		
		if( in_print_level > 1 )
			std::cout << MRQ_PREPRINT << iter << " " << master->objValue << " " << nlpfp.solver->objValue << "\n";
		
		
		if( checkTerminationCriterions(thnumber, zl, zu, iter, timeStart, clockStart, out_return_code) )
		{
			break;
		}
		
		
		r = masterMilp.changeNorm1constr( nI, intVars, nlpfp.solver->sol );
		MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
		
		r = masterMilp.addLinearizedNLConstraints(true, nlpfp.solver->sol, false, auxEval, NULL, NULL, nlpfp.solver->constr, true);
		MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
	}
	
	
	
	if( out_return_code == MRQ_HEURISTIC_SUCCESS )
	{
		if( in_solve_nlp_as_local_search_at_end )
		{
			if(in_print_level > 2)
				std::cout << MRQ_PREPRINT << "Current obj function: " << out_best_obj << ". Solving nlp problem as local search.\n";
			
			if(preprocess)
			{
				bool f1, f2;
				int r;
				//lx and ux has their own space different of nlx and prob.lx
				
				for( int i = 0; i < nI; i++ )
				{
					const int ind = intVars[i];
					lx[ind] = ux[ind] = out_best_sol[ind];
				}
				
				
				preprocessor.preprocess( in_preprocess_quad_constrs, in_preprocess_obj_function, zu, lx, ux, f1, f2, plc, puc, plc, puc);
				
				r = nlp->setnVariablesBounds( n, lx, ux );
				MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
				
				if( f2 )
				{
					for(int i = 0; i < m; i++)
					{
						r = nlp->setConstraintBounds(i, plc[i], puc[i]);
						MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
					}
				}
				
			}
			else
			{
				/*std::cout << "nI: " << nI << " out_best_obj: " << out_best_obj << "\n";
				
				for(int i = 0; i < nI; i++)
				{
					const int ind = intVars[i];
					
					std::cout << "i: " << i << " intVars[i]: " << ind << " out_best_sol[" << ind << "]: " << out_best_sol[ind] << "\n"; 
				} */
				
				int r = MRQ_fixIntVarsOnSolByList( nI, intVars, out_best_sol, *nlp );
				MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
				
			}
			
			r = nlp->solve(false);
			
			if( r == OPT_OPTIMAL_SOLUTION || nlp->feasSol )
			{
				if(in_print_level > 2)
					std::cout << MRQ_PREPRINT << "New objective function: " << nlp->objValue << "\n";
				
				if( nlp->objValue < out_best_obj )
				{
					const bool flag = tryUpdateBestSolution( thnumber, n, nlp->sol, nlp->objValue, iter, clockStart, timeStart, in_store_history_solutions  );
					
					#if MRQ_DEBUG_MODE
						if(!run_by_inside)
							assert(flag);
					#endif
				}
			}
			else
			{
				if(in_print_level > 2)
					std::cout << MRQ_PREPRINT << "Failure." << std::endl;
			}
			
		}
	}

	
termination:
	
	
	if( in_print_level > 1 )
	{
		if( out_return_code == MRQ_HEURISTIC_SUCCESS )
			std::cout << MRQ_PREPRINT << "OA Feasibility Pump heuristic found a feasible solution! ";
		else
			std::cout << MRQ_PREPRINT << "OA Feasibility Pump did not find a feasible solution! ";
	}
	
	if(plc)			free(plc);
	if(auxEval)		free(auxEval);
	if(intVars)		free(intVars);
	if(auxVars)		free(auxVars);
	if(auxConstr)	free(auxConstr);
	
	
	//if(master)		delete master;
	if(nlp)		delete nlp;
	
	
	out_feasible_solution = out_best_obj < MRQ_INFINITY;
	out_number_of_iterations = iter;
	out_lower_bound = zl;
	out_upper_bound = zu;
	
	algorithmFinalization(1, prob, lx, ux);
	
	out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
	out_clock_time = MRQ_getTime() - timeStart;
	
	if(in_print_level > 1)
		std::cout << "cpu time: " << out_cpu_time << std::endl;
	
	return out_return_code;
}
