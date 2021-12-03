/*
 * That file implements the Integrality Gap Minimization Algorithm. Here solve a sequence of problems
 * 
 * min y^t (1 - y)
 * s.t:
 * s.t.: g(x,y) <= 0
 *       zl <= f(x,y) <= zu
 *       lx <= x <= ux
 *       ly <= y <= uy
 * 
 * x and y are treated as continuous variables
 * 
 * 
 * */

#include <ctime>
#include <cmath>

#include <new>

#include "MRQ_dataStructures.hpp"
#include "MRQ_algClasses.hpp"
#include "MRQ_solvers.hpp"
//#include "MRQ_nlpSolvers.hpp"
#include "MRQ_tools.hpp"

#if MRQ_HAVE_IQUAD
	#include "iquad.hpp"
#endif

//using namespace std;

using namespace optsolvers;
using namespace muriqui;




MRQ_IGMA0::MRQ_IGMA0():MRQ_Algorithm()
{
	resetParameters();
	resetOutput();
	out_algorithm = MRQ_IGMA0_ALG; 
}



MRQ_IGMA0::~MRQ_IGMA0()
{
}


int MRQ_IGMA0::checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx, const double *ux)
{
	
	return  MRQ_isBinarieProblemAtRegion(prob, lx, ux) ? 0 : -1;
}


void MRQ_IGMA0::printParameters(std::ostream &out) const
{
	char strValue[100];
	MRQ_Algorithm::printParameters(out);
	
	out << "\n"
	MRQ_STRFFATT(in_declare_relax_infeas_if_solver_fail) << "\n"
	MRQ_STRFFATT(in_use_integers_vars_on_gap_min_problem) << "\n"
	MRQ_STRFFATT(in_stop_heuristcs_on_first_solution) << "\n"
	MRQ_STRFFATT(in_use_heuristcs) << "\n"
	MRQ_STRFFATT(in_number_of_threads) << "\n";
	
	
	MRQ_enumToStr( in_eps_updt_strategy, strValue );
	out << MRQ_STR(in_eps_updt_strategy) " " << strValue <<  '\n'; 
	
	MRQ_enumToStr( in_global_solver, strValue );
	out << MRQ_STR(in_global_solver) " " << strValue <<  '\n';
	
	
	
	out << MRQ_STRFFATT(in_absolute_obj_cut_epsilon) << "\n"
	MRQ_STRFFATT(in_decrease_fator_integer_tol) << "\n"
	MRQ_STRFFATT(in_relative_obj_cut_epsilon) << "\n"
	MRQ_STRFFATT(in_max_time_on_heuristcs) << "\n"
	MRQ_STRFFATT(in_alpha) << "\n"
	;
}


void MRQ_IGMA0::resetParameters()
{
	MRQ_Algorithm::resetParameters();
	
	in_declare_relax_infeas_if_solver_fail = true;
	in_use_integers_vars_on_gap_min_problem = false;
	in_stop_heuristcs_on_first_solution = false;
	in_use_heuristcs = true;
	in_number_of_threads = 0;
	
	in_decrease_fator_integer_tol = 0.1;
	in_absolute_obj_cut_epsilon = 1.0e-3;
	in_relative_obj_cut_epsilon = 5.0e-4;
	in_max_time_on_heuristcs = 5*60;
	in_alpha = 0.5;
	
	in_global_solver = MRQ_IQUAD;
	in_eps_updt_strategy = MRQ_IGMA0_EPS_NO_UPDATE;
}



int MRQ_IGMA0::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* globalSolverParams, MRQ_GeneralSolverParams* subGlobalSolverParams)
{
	const int n = prob.n;
	const int m = prob.m;
	const int nI = prob.getNumberOfIntegerVars();
	
	bool updtConstrBounds;
	bool intSol, improved, objCutSet = false;
	int ret, nthreads;
	const int objCutIndex = m;
	unsigned long int iter = 1;
	double aux, timeStart;
	double objCut = -MRQ_INFINITY;
	double intTolerance = in_integer_tol;
	double *plc = NULL, *puc = NULL;
	
	optsolvers::OPT_ObjCutSetter ocset;
	
	MRQ_Preprocessor preprocessor(&prob);
	
	int *intVars = NULL;
	double *auxVals = NULL;
	//double *globalSol = NULL;
	
	double *lx = run_by_inside ? nlx : prob.lx;
	double *ux = run_by_inside ? nux : prob.ux;
	
	clock_t clockStart;
	
	MRQ_NLPSolver *global = NULL;
	MRQ_NLPSolver *nlp = NULL;
	
	MRQ_NoObjWithObjCutNLPEval *gapMinEval = NULL;
	
	
	timeStart = MRQ_getTime();
	clockStart = clock();
	
	
	//zl = MRQ_max(in_lower_bound, -MRQ_INFINITY);
	//zu = in_upper_bound;
	
	
	#if MRQ_MULTITHREADING
		nthreads = in_number_of_threads < 1 ? MRQ_getNumCores() : in_number_of_threads;
	#else
		nthreads = 1;
	#endif
	
	
	{
		auto ret = algorithmInitialization(nthreads, (in_preprocess_lin_constr || in_preprocess_obj_function || in_preprocess_quad_constrs), globalSolverParams, subGlobalSolverParams, prob, lx, ux, &preprocessor, &updtConstrBounds, &plc, &puc);
		if(ret != 0)
		{
			if(in_print_level > 0 && ret != MRQ_INFEASIBLE_PROBLEM )
				MRQ_PRINTERRORNUMBER(ret);
			
			out_return_code = ret;
			goto termination;
		}
	}
	
	
	if( in_print_level > 1 )
		std::cout << "\nStarting Integrality Gap Minimization Algorithm, version 0 (IGMA0)\n\n";
	
	if(in_print_level > 3)
		printSubSolvers(true, true, true);
	
	intVars = (int *) malloc( nI * sizeof(int) );
	auxVals = (double *) malloc( n * sizeof(double) );
	global = OPT_newNLPSolver( MRQ_IQUAD );
	
	if( !intVars || !auxVals || !global)
	{
		if(in_print_level > 1)
			MRQ_PRINTMEMERROR;
		out_return_code = MRQ_MEMORY_ERROR;
		goto termination;
	}
	
	
	prob.getIntegerIndices(intVars);
	
	
	nlp = OPT_newNLPSolver( in_nlp_solver );
	
	if(!nlp)
	{
		if(in_print_level > 0)
			MRQ_PRINTMEMERROR;
	
		out_return_code = MRQ_MEMORY_ERROR;
		goto termination;
	}
	
	ret = MRQ_setNLPRelaxProb( prob, lx, ux, plc, puc, nlp, true, true, true, false, thnumber, in_set_special_nlp_solver_params, NULL, in_number_of_threads, in_max_cpu_time, in_max_time, 0, 0 );
	if( ret != 0 )
	{
		if(in_print_level > 0)
			MRQ_PRINTERRORNUMBER(ret);
		out_return_code = MRQ_NLP_SOLVER_ERROR;
		goto termination;
	}
	
	if( run_by_inside )
	{
		if( !std::isinf(insideSolverMaxTime) )
			ret += nlp->setMaxTime(insideSolverMaxTime );
	}
		
	//solving continuous relaxation
	
	
	if(in_eps_updt_strategy == MRQ_IGMA0_EPS_BIN_SEARCH && zl <= -MRQ_INFINITY )
	{
		if( in_print_level > 1 )
			std::cout << "Solving NLP relaxation\n";
		
		ret = nlp->solve(false);
		
		if( ret == OPT_OPTIMAL_SOLUTION )
		{
			if(in_print_level > 1)
				std::cout << "NLP relaxation solution: " << nlp->objValue << "\n";
			
			if( nlp->objValue > zl )
				zl = nlp->objValue;
		}
		else
		{
			if(in_print_level > 0)
				std::cerr << "Failure for solving continuous relaxation\n";
		}
	}
	
	
	
	//tryng get a good feasilbe solution by means of heuristics. Note that IGMA does not look for objective function. So, the first solution gotten can be very bad, and it can be perverse to the algorithm. To remediate it, we try start from a good solution found by some heuristic
	
	if( in_use_heuristcs )
	{
		int ret;
		double obj;
		MRQ_HeuristicExecutor heurExec;
		
		if(in_print_level > 2)
			std::cout << "Running heuristics to get an initial feasible solution\n";
		
		heurExec.in_max_total_clock_time = in_max_time_on_heuristcs;
		
		heurExec.setMaxTimes( in_max_time_on_heuristcs );
		
		if( nlp->feasSol )
			heurExec.setInitialSolution( n, nlp->sol );
		
		heurExec.in_stop_on_first_sol = in_stop_heuristcs_on_first_solution;
		
		ret = heurExec.insideRun( prob, NULL, NULL, zu, obj, auxVals, false, NULL, thnumber, 1, in_max_time_on_heuristcs, lx, ux );
		
		
		
		
		//i = nlp->retCode == OPT_OPTIMAL_SOLUTION || nlp->retCode == MRQ_NLP_FEASIBLE_SOLUTION;
		//ret = getFeasibleSolution(prob, lx, ux, nlp->feasSol ? nlp->sol : NULL,  in_max_time_on_heuristcs, in_stop_heuristcs_on_first_solution, out_best_sol, obj);
		
		if(ret == MRQ_HEURISTIC_SUCCESS)
		{
			if(obj <= zu)
			{
				//out_best_obj = zu = obj;
				//MRQ_copyArray(n, auxVals, out_best_sol);
				
				
				const bool flag = tryUpdateBestSolution( thnumber, n, auxVals, obj, iter, clockStart, timeStart, in_store_history_solutions );
				
				#if MRQ_DEBUG_MODE
					assert(flag);
				#endif
				
				if( in_print_level > 1 )
					std::cout << MRQ_PREPRINT "Heuristcs found an integer feasible solution! Obj: " << out_best_obj << "\n";
			}
		}
		else if(ret == MRQ_OPTIMAL_SOLUTION)
		{
			if( obj <= zu )
			{
				if( in_print_level > 1 )
					MRQ_PRINTMSG("Heuristcs found an optimal solution!\n");
				
				out_return_code = MRQ_OPTIMAL_SOLUTION;
			}
			else
			{
				if( in_print_level > 1 )
					MRQ_PRINTERRORMSG("Heuristcs found an optimal solution worse than upper bound!\n");
				
				out_return_code = MRQ_INFEASIBLE_PROBLEM;
			}
			
			//out_best_obj = zu = obj;
			//MRQ_copyArray(n, auxVals, out_best_sol);
			
			const bool flag = tryUpdateBestSolution( thnumber, n, auxVals, obj, iter, clockStart, timeStart, in_store_history_solutions );
			
			#if MRQ_DEBUG_MODE
				assert(flag);
			#endif
			
			goto termination;
		}
		else
		{
			if( in_print_level > 1 )
				MRQ_PRINTMSG("Heuristcs had failure to find a integer feasible solution!\n");
		}
	}
	
	
	
	ret = MRQ_setNLPRelaxProb( prob, lx, ux, plc, puc, global, false, true, true, in_use_integers_vars_on_gap_min_problem, thnumber, false, globalSolverParams, nthreads, in_max_cpu_time, in_max_time, 0, 0);
	
	if( ret != 0 )
	{
		if(in_print_level > 0)
			MRQ_PRINTERRORNUMBER(ret);
		
		out_return_code = MRQ_NLP_SOLVER_ERROR;
		goto termination;
	}
	
	
	//seting objective funtion
	{
		int r;
		
		global->setObjSense(optsolvers::OPT_MINIMIZE);
		
		MRQ_setAllArray(nI, auxVals, 1.0);
		
		r = global->setObjLinearCoefs(nI, intVars, auxVals);
		
		
		MRQ_setAllArray(nI, auxVals, -2.0);
		r += global->setObjQuadMatrix(nI, intVars, intVars, auxVals);
		if( r != 0 )
		{
			if(in_print_level > 0)
				MRQ_PRINTERRORNUMBER(r);
			
			out_return_code = MRQ_NLP_SOLVER_ERROR;
			goto termination;
		}
	}
	
	
	if( run_by_inside )
	{
		if( !std::isinf(insideSolverMaxTime) )
			ret += global->setMaxTime(insideSolverMaxTime );
	}
	
	#if MRQ_HAVE_IQUAD
	if( in_global_solver == MRQ_IQUAD )
	{
		iquad::IQD_BranchAndBound *bb = ((OPT_Iquad *) global)->bb;
		
		
		bb->mrq_gap_min = true;
		
		bb->in_print_level = 100;
		bb->in_printing_frequency = 1;
		
		bb->in_number_of_threads = 1;
		bb->in_set_special_solver_parameters = false; //we really need that
		bb->in_use_heuristics_on_integer_probs = false;
		
		bb->in_split_way = iquad::IQD_SQ_DIAG_DOM;
		
		bb->in_relative_convergence_tol = 1.0e-5;
		
		bb->in_certain_of_nonconvex_Q_on_obj_f = true;
		
		//bb->in_alpha_to_integer_branch_choosing = 0.0;
		
		//bb->in_relax_prob_writing_mode = iquad::IQD_RPWM_WRITING;
		
		if( in_nlp_solver == MRQ_NLP_MOSEK )
			bb->in_nlp_solver = iquad::IQD_NLP_MOSEK;
		else if( in_nlp_solver == MRQ_IPOPT )
			bb->in_nlp_solver = iquad::IQD_NLP_IPOPT;
		else if( in_nlp_solver == MRQ_WORHP )
			bb->in_nlp_solver = iquad::IQD_NLP_WORHP;
		else
			assert(false);
	}
	#endif
	
	
	aux = (intTolerance - (intTolerance*intTolerance))*nI; //maximum value that a feasible solution can be o igm problem.
	
	//printf("zuCut: %f\n", aux);
	
	ret += global->setObjCutUpperBound( aux + 1.0e-5 ); //if lower bound pass by that value, the optimization will stop on iquad. Note we still can have problems due iquad will stop before zl achieve that value due to absolute and relative tolerances. You have to solve it some day... Anyway, to calculate that value, we consider all integer variable are in the maximum allowable gap. So, is little unprobable we have problem about this...
	
	
	
	if( prob.hasObjNLTerm() )
	{
		gapMinEval = new (std::nothrow) MRQ_NoObjWithObjCutNLPEval( prob.hasNLConstraints(),  prob.getNonLinearEvaluationObject(), prob.m, prob.objFactor );
		
		if( !gapMinEval )
		{
			MRQ_PRINTMEMERROR;
			return MRQ_MEMORY_ERROR;
		}
		
		//gapMinEval->auxVals = auxVals;
	}
	
	
	if( zu < MRQ_INFINITY )
	{
		const int re = global->addConstraints(1);
		if( re != 0)
		{
			if(in_print_level > 0)
				MRQ_PRINTERRORNUMBER(re);
			
			out_return_code = MRQ_NLP_SOLVER_ERROR;
			goto termination;
		}
		
		
		aux = MRQ_zuWithTol(zu, in_absolute_convergence_tol, in_relative_convergence_tol);
		
		//ret += global->setObjectiveCut( aux );
		
		const int r = ocset.setObjCut( global, objCutIndex, prob, aux, gapMinEval );
		
		if(r != 0)
		{
			if( in_print_level > 0 )
				MRQ_PRINTERRORNUMBER(r);
			out_return_code = MRQ_NLP_SOLVER_ERROR;
			goto termination;
		}
		
		objCutSet = true;
	}
	
	
	if(ret != 0)
	{
		if( in_print_level > 0 )
			MRQ_PRINTERRORMSG("Error at setting global optimization IGMA problem!");
		
		out_return_code = MRQ_NLP_SOLVER_ERROR;
		goto termination;
	}
	
	if( in_global_solver == MRQ_IQUAD )
		ret = ( (OPT_Iquad *) global )->solveWParams(false, true, false, false, subGlobalSolverParams, NULL);
	else
		ret = global->solve(false);
	
	
	//printf("Global return: %d\n", ret);
	
	if( ret != OPT_OPTIMAL_SOLUTION )
	{
		if( ret == OPT_INFEASIBLE_PROBLEM )
		{
			if( out_best_obj < MRQ_INFINITY )
				out_return_code = MRQ_OPTIMAL_SOLUTION;
			else
				out_return_code = MRQ_INFEASIBLE_PROBLEM;
		}
		else
			out_return_code = MRQ_NLP_SOLVER_ERROR;
		
		goto termination;
	}
	
	intSol = prob.isIntegerSolution(global->sol, intTolerance);
	if( !intSol )
	{
		MRQ_PRINTERRORMSG("Solution is not integer\n");
		if( out_best_obj < MRQ_INFINITY )
			out_return_code = MRQ_OPTIMAL_SOLUTION;
		else
			out_return_code = MRQ_INFEASIBLE_PROBLEM;
		goto termination;
	}
	
	
	if(in_print_level > 1)
		std::cout << "iter   igmp obj   zl     zu      gap\n";
	
	
	
	while(true)
	{
		
		improved = false;
		
		if( intSol )
		{
			//applying nlp resolution as local search
			
			/*for(i= 0; i < nI; i++)
			{
				k = intVars[i];
				aux = round( global->sol[k] );
				
				ret = nlp->setVarBound( k, aux, aux );
				if( ret != 0 )
				{
					out_return_code = MRQ_NLP_SOLVER_ERROR;
					goto termination;
				}
			} */
			
			if(in_print_level > 4)
				std::cout << "Applying local search. ";
			
			ret = MRQ_fixIntVarsOnSolByList(nI, intVars, global->sol, *nlp);
			
			#if MRQ_DEBUG_MODE
				assert(ret == 0);
			#endif
			
			
			ret = nlp->solve(false, true, false, false);
			
			if( nlp->feasSol ) //if( ret == MRQ_OPTIMAL_SOLUTION || ret == MRQ_NLP_FEASIBLE_SOLUTION )
			{
				//in general solution is better, but we can have only a nlp feasible solution
				
				if(in_print_level > 4)
					std::cout << "Solution : " << nlp->objValue << "\n";
				
				if( nlp->objValue < zu )
				{
					//out_best_obj = zu = nlp->objValue;
					//MRQ_copyArray(n, nlp->sol, out_best_sol );
					//improved = true;
					
					improved = tryUpdateBestSolution( thnumber, n, nlp->sol, nlp->objValue, iter, clockStart, timeStart, in_store_history_solutions );
					
					#if MRQ_DEBUG_MODE
						assert( improved == true);
					#endif
				}
				/*else if( nlp->objValue == lastObjLS)
				{//we repeat the solution in terms of objective function...
					printf("detectei ciclo!\n");
					MRQ_getchar();
				}
				
				lastObjLS = nlp->objValue; */
			}
			else if( ret == OPT_INFEASIBLE_PROBLEM )
			{
				if(in_print_level > 4)
					std::cout << "Infeasible subproblem.\n";
					
				if( in_eps_updt_strategy == MRQ_IGMA0_EPS_BIN_SEARCH )
				{
					if( objCut > zl )
						zl = objCut;
				}
			}
			else
			{
				if(in_print_level > 4)
					std::cout << "Error.\n";
				
				
				if( in_eps_updt_strategy == MRQ_IGMA0_EPS_BIN_SEARCH )
				{
					if( in_declare_relax_infeas_if_solver_fail )
					{
						if( objCut > zl )
							zl = objCut;
					}
					else
					{
						out_return_code = MRQ_NLP_SOLVER_ERROR;
						goto termination;
					}
				}
				else
				{
					
					ret = prob.objEval(thnumber, true, global->sol, aux);
					
					if(ret != 0)
					{
						out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
						goto termination;
					}
					
					
					#if MRQ_DEBUG_MODE
						if( aux > out_best_obj )
						{
							std::cerr << MRQ_PREPRINT "IGMA1: Error! New solution is worst than previous.";
							//MRQ_getchar();
							
							out_return_code = MRQ_UNDEFINED_ERROR;
							goto termination;
						}
					#endif
					
					//printf("Funcao objetivo solucao: %f\n", aux);
					
					
					//out_best_obj = zu = aux;
					//MRQ_copyArray(n, global->sol, out_best_sol);
					
					const bool flag = tryUpdateBestSolution( thnumber, n, global->sol, aux, iter, clockStart, timeStart, in_store_history_solutions );
					
					#if MRQ_DEBUG_MODE
						assert(flag);
					#endif
				}
			}
			
		}
		
		printf("%-5ld  %+-14e  %+-14e  %+-14e  %+-14e\n", iter, global->objValue, zl, zu, zu - zl);
		
		
		if( checkTerminationCriterions(thnumber, zl, zu, iter, timeStart, clockStart, out_return_code) )
		{
			break;
		}
		
		MRQ_getchar();
		
		//end of iteration. From now, we start the next
		iter++;
		
		
		if( in_eps_updt_strategy == MRQ_IGMA0_EPS_BIN_SEARCH )
		{
			if( !improved )
			{
				/*ret = prob.objEval(thnumber, true, global->sol, aux);
				
				if(ret != 0)
				{
					out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
					goto termination;
				}
				
				objCut = MRQ_zuWithTol(aux, in_absolute_obj_cut_epsilon, in_relative_convergence_tol); */
				
				aux = MRQ_zuWithTol(zu, in_absolute_obj_cut_epsilon, in_relative_obj_cut_epsilon);
				
				if( aux > objCut )
				{
					objCut = aux;
				}
				else
				{
					out_return_code = out_best_obj < MRQ_INFINITY ? MRQ_OPTIMAL_SOLUTION : MRQ_UNDEFINED_ERROR;
					break;
				}
				
				intTolerance *= in_decrease_fator_integer_tol;
			}
			else
			{
				objCut = (1 - in_alpha)*zl +  (in_alpha)*zu;
			}
			
			//objCut = (zl + zu)/2.0;
		}
		else
		{
			objCut = MRQ_zuWithTol(zu, in_absolute_obj_cut_epsilon, in_relative_obj_cut_epsilon);
			
			/*if( objCut >= MRQ_INFINITY )
				objCut = sqrt(MRQ_INFINITY); */
		}
		
		printf("zl: %f zu: %f objCut: %f\n", zl, zu, objCut);
		
		
		//ret = global->setObjectiveCut( objCut );
		
		if( objCutSet )
		{
			ret = ocset.updateObjCut( global, objCutIndex, prob, objCut );
		}
		else
		{
			const int r = global->addConstraints(1);
			if( r != 0)
			{
				if(in_print_level > 0)
					MRQ_PRINTERRORNUMBER(r);
				
				out_return_code = MRQ_NLP_SOLVER_ERROR;
				break;
			}
			
			
			ret = ocset.setObjCut( global, objCutIndex, prob, objCut, gapMinEval );
			objCutSet = true;
		}
		
		if( ret != 0)
		{
			if(in_print_level > 0)
				MRQ_PRINTERRORNUMBER(ret);
			
			out_return_code = MRQ_NLP_SOLVER_ERROR;
			break;
		}
		
		if( in_global_solver == MRQ_IQUAD )
			ret = ( (OPT_Iquad *) global )->solveWParams(false, true, true, true, subGlobalSolverParams, NULL);
		else
			ret = global->solve(false);
		
		printf("global ret: %d obj: %f\n", ret, global->objValue);
		for(int w = 0; w < n; w++)
			printf("sol[%d]: %f\t", w, global->sol[w] );
		printf("\n");
		
		if( ret == OPT_OPTIMAL_SOLUTION )
		{
			intSol = prob.isIntegerSolution(global->sol, intTolerance);
			
			if(!intSol)
			{
				printf("Solucao nao é inteira!\n");
				if(in_eps_updt_strategy == MRQ_IGMA0_EPS_BIN_SEARCH)
				{
					if( objCut > zl )
						zl = objCut;
					
					if( !improved )
					{//we can stop
						out_return_code = MRQ_OPTIMAL_SOLUTION;
						break;
					}
				}
				else
				{
					out_return_code = MRQ_OPTIMAL_SOLUTION;
					break;
				}
			}
		}
		else if( ret == OPT_INFEASIBLE_PROBLEM )
		{
			if( objCut > zl )
				zl = objCut;
			
			if( in_eps_updt_strategy == MRQ_IGMA0_EPS_NO_UPDATE )
			{
				out_return_code = MRQ_OPTIMAL_SOLUTION; //I am not sure about it... 
				break;
			}
			
			intSol = false;
		}
		else
		{
			out_return_code = MRQ_NLP_SOLVER_ERROR;
			goto termination;
		}
		
		
		
		
		
	}
	
	
	
termination:
	
	if(in_print_level > 0)
	{
		if(out_return_code == MRQ_OPTIMAL_SOLUTION)
		{
			std::cout << MRQ_PREPRINT "Optimal solution found. Objective: " << out_best_obj << "\n";
		}
		else if(out_return_code == MRQ_INFEASIBLE_PROBLEM)
		{
			std::cout << MRQ_PREPRINT "Infeasible problem.\n";
		}
		else if(out_return_code == MRQ_MAX_ITERATIONS_STOP )
		{
			std::cout << MRQ_PREPRINT "Max iteration number reached!\n";
		}
		else if(out_return_code == MRQ_MAX_TIME_STOP)
		{
			std::cout << MRQ_PREPRINT "Max time reached!\n";
		}
		else if(out_return_code == MRQ_NLP_SOLVER_ERROR)
		{
			std::cout << MRQ_PREPRINT "Error on global solver!\n";
		}
	}
	
	
	if(intVars)	free(intVars);
	if(auxVals)	free(auxVals);
	
	if(nlp)		delete nlp;
	if(global)	delete global;
	
	if(gapMinEval)	delete gapMinEval;
	
	
	out_feasible_solution = out_best_obj < MRQ_INFINITY;
	out_number_of_iterations = iter;
	out_lower_bound = zl;
	out_upper_bound = zu;
	
	algorithmFinalization(nthreads, prob, lx, ux);
	
	out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
	out_clock_time = MRQ_getTime() - timeStart;
	
	if(in_print_level > 1)
		std::cout << "cpu time: " << out_cpu_time << "\n";
	
	return out_return_code;
}



int MRQ_IGMA0::setIntegerParameter(const char *name, const long int value)
{
	int ret = MRQ_Algorithm::setIntegerParameter(name, value);
	
	if( ret == 0 )
		return 0;
	
	ret = 0;
	
	
	if( MRQ_setAtt<bool>( MRQ_STRATT(in_declare_relax_infeas_if_solver_fail), name, value ) == 0 )
	{}
	else if( MRQ_setAtt<bool>( MRQ_STRATT(in_stop_heuristcs_on_first_solution), name, value ) == 0 )
	{}
	else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_heuristcs), name, value ) == 0 )
	{}
	else if( MRQ_setAtt<int>( MRQ_STRATT(in_number_of_threads), name, value ) == 0 )
	{}
	else
		ret = MRQ_NAME_ERROR;
	
	
	return ret;
}


int MRQ_IGMA0::setDoubleParameter(const char *name, const double value)
{
	int ret = MRQ_Algorithm::setDoubleParameter(name, value);
	
	if( ret == 0 )
		return 0;
	
	ret = 0;
	
	
	if( MRQ_setAtt( MRQ_STRATT(in_absolute_obj_cut_epsilon), name, value ) == 0 )
	{}
	else if( MRQ_setAtt( MRQ_STRATT(in_decrease_fator_integer_tol), name, value ) == 0 )
	{}
	else if( MRQ_setAtt( MRQ_STRATT(in_relative_obj_cut_epsilon), name, value ) == 0 )
	{}
	else if( MRQ_setAtt( MRQ_STRATT(in_max_time_on_heuristcs), name, value ) == 0 )
	{}
	else if( MRQ_setAtt( MRQ_STRATT(in_alpha), name, value ) == 0 )
	{}
	else
		ret = MRQ_NAME_ERROR;
	
	
	return ret;
}


int MRQ_IGMA0::setStringParameter(const char *name, const char *value)
{
	int r, ret = MRQ_Algorithm::setStringParameter(name, value);
	
	if( ret == 0 )
		return 0;
	
	ret = 0;
	
	
	if( (r = MRQ_setStrAtt( MRQ_STRATT(in_eps_updt_strategy), name, value ) ) >= 0 )
	{
		ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
	}
	else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_global_solver), name, value ) ) >= 0 )
	{
		ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
	}
	else
		ret = MRQ_NAME_ERROR;
	
	return ret;
}




























