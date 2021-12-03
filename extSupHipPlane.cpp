/*That file contains a simple implementation of
* Extnded Supporting Hyperplane algorithm.
*
* References:
*
* Kronqvist & Lundell & Westerlund, The extended supporting hyperplane algorithm for convex mixed-integer nonlinear programming. Journal of Global Optimization 64 (2016), pages 249-272.
*
*
* Author: Wendel Melo
*
* Date: 10-Sept-2016
*/


#include <cassert>
#include <climits>
#include <cmath>
#include <cstdio>
#include <ctime>


#include <iostream>
#include <new>

#include "BBL_tools.hpp"

#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_advanced.hpp"



using namespace optsolvers;
using namespace muriqui;



MRQ_ExtSupHypPlane::MRQ_ExtSupHypPlane():MRQ_LinearApproxAlgorithm()
{
    resetParameters();
    resetOutput();
    out_algorithm = MRQ_ESH_ALG;
}


MRQ_ExtSupHypPlane::~MRQ_ExtSupHypPlane()
{
}


int MRQ_ExtSupHypPlane::checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx, const double *ux)
{
    const int m = prob.m;
    
    const bool *nlConstr = prob.nlConstr;
    const double *lc = prob.lc, *uc = prob.uc;
    MRQ_SparseMatrix *QC = prob.QC;
    
    
    for( int i = 0; i < m; i++ )
    {
        //by now, we are not handling nonlinear double bounded constraint
        if( lc[i] > -MIP_INFINITY && uc[i] < MIP_INFINITY && (nlConstr[i] || QC[i].getNumberOfElements() > 0) )
        {
            std::cerr << MRQ_PREPRINT "Problem has double bound nonlinear constraints. By now, Extended Suporting Hyperplane cannot be applied\n";
            
            return MRQ_ALG_NOT_APPLICABLE;
        }
    }
    
    return 0;
}


void MRQ_ExtSupHypPlane::printParameters(std::ostream &out) const
{
    char strValue[100];
    
    strValue[0] = '\0'; //only by safe if some parameter does not find the string value...
    
    MRQ_LinearApproxAlgorithm::printParameters(out);
    out << "\n"
    
    
    //MRQ_STRFFATT(in_delete_linearizations_from_int_vars_fixing) << "\n"
    
    MRQ_STRFFATT(in_fix_int_vars_to_try_improve_cut) << "\n"
    MRQ_STRFFATT(in_linearize_on_interior_sol) << "\n"
    MRQ_STRFFATT(in_try_solve_interior_problem_if_cont_relax_fail) << "\n"
    MRQ_STRFFATT(in_max_lp_subiters) << "\n"
    MRQ_STRFFATT(in_max_lp_subiter_to_improve_obj_app) << "\n";
    
    MRQ_enumToStr( in_lp_constr_linearisation_strategy, strValue );
    out << MRQ_STR(in_lp_constr_linearisation_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr( in_interior_point_strategy, strValue );
    out << MRQ_STR(in_starting_point_strategy) " " << strValue << "\n"
    
    MRQ_STRFFATT(in_cont_relax_absolute_convergence_tol) << "\n"
    MRQ_STRFFATT(in_cont_relax_relative_convergence_tol) << "\n"
    MRQ_STRFFATT(in_absolute_tol_to_check_previous_sol) << "\n"
    MRQ_STRFFATT(in_relative_tol_to_check_previous_sol) << "\n"
    MRQ_STRFFATT(in_delta_to_inc_eps_to_active_constraints_to_linearization) << "\n"
    MRQ_STRFFATT(in_eps_to_enforce_interior_sol_on_cont_relax_sol) << "\n"
    MRQ_STRFFATT(in_eps_to_interior_sol) << "\n"
    MRQ_STRFFATT(in_eps_to_line_search) << "\n"
    ;
}


void MRQ_ExtSupHypPlane::resetOutput()
{
    MRQ_LinearApproxAlgorithm::resetOutput();
    
    out_number_of_lp_iterations = 0;
}


void MRQ_ExtSupHypPlane::resetParameters()
{
    MRQ_LinearApproxAlgorithm::resetParameters();
    
    in_printing_frequency = 10;
    in_constr_linearisation_strategy = MRQ_CLS_ONLY_INFEAS_AND_ACTIVE;
    in_eps_to_active_constr_to_linearisation = 0.2;
    
    
    in_fix_int_vars_to_try_improve_cut = false;
    in_linearize_on_interior_sol = false;
    //in_delete_linearizations_from_int_vars_fixing = false;
    in_try_solve_interior_problem_if_cont_relax_fail = true;
    in_max_lp_subiters = 50;
    in_max_lp_subiter_to_improve_obj_app = 10;
    
    in_lp_constr_linearisation_strategy = MRQ_CLS_ONLY_INFEAS_AND_ACTIVE;
    in_interior_point_strategy = muriqui::MRQ_ESHP_IPS_MOST_INTERIOR;
    
    in_cont_relax_absolute_convergence_tol = in_absolute_convergence_tol;
    in_cont_relax_relative_convergence_tol = in_relative_convergence_tol;
    
    in_absolute_tol_to_check_previous_sol = 1e-6;
    in_relative_tol_to_check_previous_sol = 1e-6;
    in_delta_to_inc_eps_to_active_constraints_to_linearization = 0.05;
    
    in_eps_to_enforce_interior_sol_on_cont_relax_sol = 0.025;
    in_eps_to_interior_sol = 0.01;
    in_eps_to_line_search = 1e-6;
}


int MRQ_ExtSupHypPlane::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_LinearApproxAlgorithm::setIntegerParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    //if( MRQ_setAtt<bool>( MRQ_STRATT(in_delete_linearizations_from_int_vars_fixing), name, value ) == 0 ); else 
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_fix_int_vars_to_try_improve_cut), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_linearize_on_interior_sol), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_try_solve_interior_problem_if_cont_relax_fail), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_max_lp_subiters), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    return ret;
}


int MRQ_ExtSupHypPlane::setDoubleParameter(const char *name, const double value)
{
    int ret = MRQ_LinearApproxAlgorithm::setDoubleParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt( MRQ_STRATT(in_absolute_tol_to_check_previous_sol), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_relative_tol_to_check_previous_sol), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_delta_to_inc_eps_to_active_constraints_to_linearization), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_eps_to_enforce_interior_sol_on_cont_relax_sol), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_eps_to_interior_sol), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_eps_to_line_search), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    return ret;
}


int MRQ_ExtSupHypPlane::setStringParameter(const char *name, const char *value)
{
    int ret = MRQ_LinearApproxAlgorithm::setStringParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    if( (ret = MRQ_setStrAtt( MRQ_STRATT(in_interior_point_strategy), name, value ) ) >= 0 )
    {
        ret = ret == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else if( (ret = MRQ_setStrAtt( MRQ_STRATT(in_lp_constr_linearisation_strategy), name, value ) ) >= 0 )
    {
        ret = ret == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}



static inline int MRQ_changeToMILPLoop( const int nI, const int *intVars, OPT_LPSolver &master )
{
    //std::cout << "******************************************************************\n";
    //MRQ_getchar();
    for(int i = 0; i < nI; i++)
    {
        const int r = master.setVariableType( intVars[i], OPT_VT_INTEGER );
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    }
    
    return 0;
}



int MRQ_ExtSupHypPlane::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    const double timeStart = MRQ_getTime();
    clock_t clockStart = clock();
    
    
    const int n = prob.n;
    const int m = prob.m;
    const int nI = prob.getNumberOfIntegerVars();
    
    const double &epsToIntSol = in_eps_to_interior_sol;
    const double &deltaEpsActLin = in_delta_to_inc_eps_to_active_constraints_to_linearization;
    const double &absTolToTestPreviousSol = in_absolute_tol_to_check_previous_sol; //1e-6;
    const double &relTolToTestPreviousSol = in_relative_tol_to_check_previous_sol;
    
    bool linearizeObj, updtConstrBounds, linOnInteriorSol = in_linearize_on_interior_sol;
    bool milpLoop = false, feasSol;
    bool solveNLPIPProb = in_interior_point_strategy == muriqui::MRQ_ESHP_IPS_MOST_INTERIOR;
    int ret;
    int constrLinearisationStrategy = in_lp_constr_linearisation_strategy;
    
    unsigned long int iter = 0, iterlpfeas = 0;
    double robj;
    double epsToActConstrToLin = in_eps_to_active_constr_to_linearisation;
    
    double zu_cont = MRQ_INFINITY; //upper bound of continuous relaxation...
    
    
    bool *auxCEval = NULL;
    int *intVars = NULL, *auxInds = NULL;
    double *plc = NULL, *puc = NULL;
    double *isol = NULL, *auxConstr = NULL; //isol is the interior solution
    double *rsol, *lsSol, *auxVals = NULL;
    
    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux; 
    
    double NLPCpuTime, NLPClockTime;
    double *pNLPCpuTime = in_measure_nlp_time ? &NLPCpuTime : NULL;
    double *pNLPClockTime = in_measure_nlp_time ? &NLPClockTime : NULL;
    
    OPT_LPSolver *master, *pnlp;
    MRQ_MasterMILPProb masterMilp;
    MRQ_NLPIPProblem *ipNLP = NULL; //interior point problem
    MRQ_NLPSolver *nlp = NULL;
    MRQ_GradientsEvaluation  gradEval;
    MRQ_Preprocessor preprocessor(&prob);
    
    
    nthreads = 1;
    nthreads_lazy = in_number_of_threads > 0 ? in_number_of_threads : branchAndBound::BBL_getNumCores() ;
    if( in_milp_solver == MRQ_GUROBI )
        nthreads_lazy = 1; //gurobi apply multithreading to solve the problem, but only thread 0 is called to add lazy constraints.
    
    
    {
        auto ret = algorithmInitialization(1, (in_preprocess_lin_constr || in_preprocess_obj_function || in_preprocess_quad_constrs), milpSolverParams, nlpSolverParams, prob, lx, ux, &preprocessor, &updtConstrBounds, &plc, &puc ); //that algorithm is monothread...
        MRQ_IFERRORGOTOLABEL(ret, out_return_code, ret, termination);
    }
    
    //we do not need preprocessor more...
    MRQ_secFree(plc);
    preprocessor.deallocateMemory();
    
    if( in_print_level > 1 )
    {
        std::cout << "\n";
        MRQ_PRINTMSG("Starting Extended Supporting Hyperplane Algorithm\n\n");
    }
    
    if( in_print_level > 3 )
        printSubSolvers(true, true, false);
    
    ret = gradEval.initialize(thnumber, &prob);
    MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MEMORY_ERROR, termination);
    
    {
        int sizeAuxInds = n + 1;
        #if 0
        if( in_delete_linearizations_from_int_vars_fixing )
            sizeAuxInds = MRQ_max(sizeAuxInds, m - prob.getNumberOfLinearConstraints() + 1);
        #endif
        MRQ_malloc(auxInds,  sizeAuxInds);
    }
    
    MRQ_malloc(intVars, nI); 	// intVars = (int *) malloc( nI * sizeof(int) );
    MRQ_malloc(isol, 3*n); 		//isol = (double *) malloc( 3*n * sizeof(double) );
    MRQ_malloc(auxVals, n+1);	//auxVals = (double *) malloc( (n+1) *sizeof(double) );
    MRQ_malloc(auxCEval, m); //auxCEval = (bool *) malloc(m * sizeof(bool) );
    MRQ_malloc(auxConstr, m); //auxConstr = (double *) malloc(m * sizeof(double) );
    MRQ_IFMEMERRORGOTOLABEL(!auxInds || !intVars || !isol || !auxVals || !auxCEval || !auxConstr, out_return_code, termination);
    
    rsol = &isol[n];
    lsSol = &isol[2*n];
    
    
    prob.getIntegerIndices(intVars);
    
    
    for(int i = 0; i < m; i++)
        auxCEval[i] = prob.nlConstr[i] || prob.QC[i].getNumberOfElements() > 0;
    
    
    if( in_interior_point_strategy == MRQ_ESHP_IPS_CLOSED_TO_CONT_RELAX_SOL )
    {
        nlp = OPT_newNLPSolver(in_nlp_solver);
        MRQ_IFMEMERRORGOTOLABEL(!nlp, out_return_code, termination);
        
        ret = MRQ_setNLPRelaxProb(prob, lx, ux, plc, puc, nlp, true, true, true, false, thnumber, in_set_special_nlp_solver_params, nlpSolverParams, 1, in_max_cpu_time, in_max_time, 0, 0);
        MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
        
        
        const bool *nlConstr = prob.nlConstr;
        const double *lc = plc ? plc : prob.lc, *uc = plc ? puc : prob.uc;
        const MRQ_SparseMatrix *QC = prob.QC;
        
        for( int i = 0; i < m; i++ )
        {
            if( nlConstr[i] || ( QC[i].getNumberOfElements() > 0 && !in_set_quadratics_in_master_problem ) )
            {
                double nlci, nuci;
                
                if(lc[i] > -MIP_INFINITY)
                {
                    const double eps = MRQ_abs(lc[i]) > 1.0 ? MRQ_abs(lc[i])*in_eps_to_enforce_interior_sol_on_cont_relax_sol : in_eps_to_enforce_interior_sol_on_cont_relax_sol ;
                    
                    nlci = lc[i] + eps;
                }
                else
                    nlci = -OPT_INFINITY;
                
                
                if(uc[i] < MIP_INFINITY)
                {
                    const double eps = MRQ_abs(uc[i]) > 1.0 ? MRQ_abs(uc[i])*in_eps_to_enforce_interior_sol_on_cont_relax_sol : in_eps_to_enforce_interior_sol_on_cont_relax_sol ;
                    
                    nuci = uc[i] - eps;
                }
                else
                    nuci = OPT_INFINITY;
                
                const int r = nlp->setConstraintBounds(i, nlci, nuci);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
            }
        }
        
        
        pnlp = nlp;
        
        if( run_by_inside )
        {
            if( !std::isinf(insideSolverMaxTime) )
            {
                int r = pnlp->setMaxTime(insideSolverMaxTime);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
            }
        }
        
        
        pnlp->solveAndGetTime(pNLPCpuTime, pNLPClockTime, false);
        
        if(in_measure_nlp_time)
        {
            out_cpu_time_of_nlp_solving += *pNLPCpuTime;
            out_clock_time_of_nlp_solving += *pNLPClockTime;
        }
        
        out_number_of_nlp_probs_solved++;
        
        if( pnlp->feasSol == false )
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORMSG("Error at solving problem to get an interior solution");
            
            if( in_try_solve_interior_problem_if_cont_relax_fail )
            {
                solveNLPIPProb = true;
            }
            else
            {
                out_return_code = MRQ_NLP_SOLVER_ERROR;
                goto termination;
            }
        }
        #if 0
        else
        {
            //if in_starting_point_strategy == MRQ_ESHP_SPS_CLOSED_TO_CONT_RELAX_SOL, if we have a feasible solution, we got success to have an interior sol...
            if( in_interior_point_strategy == MRQ_ESHP_IPS_MOST_INTERIOR && pnlp->objValue > -epsToIntSol )
            {
                //if(in_print_level > 0)
                    //MRQ_PRINTERRORMSG("Failure to obtain a valid interior solution");
                
                if( in_try_solve_interior_problem_if_cont_relax_fail )
                {
                    solveNLPIPProb = true;
                }
                else
                {
                    out_return_code = MRQ_INITIAL_SOLUTION_ERROR;
                    goto termination;
                }
            }
        }
        #endif
    }
    
    
    if( solveNLPIPProb )
    {
        ipNLP = new (std::nothrow) MRQ_NLPIPProblem;
        MRQ_IFMEMERRORGOTOLABEL(!ipNLP, out_return_code, termination);
        

        ret = ipNLP->setProblem( in_nlp_solver, prob , lx, ux, nlpSolverParams, thnumber, in_set_special_nlp_solver_params, !in_set_quadratics_in_master_problem, in_number_of_threads, in_max_cpu_time, in_max_time);
        MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
        
        
        pnlp = ipNLP->solver;
        
        if( run_by_inside )
        {
            if( !std::isinf(insideSolverMaxTime) )
            {
                ret = pnlp->setMaxTime( insideSolverMaxTime );
                MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
            }
        }
        
        pnlp->solveAndGetTime(pNLPCpuTime, pNLPClockTime, false);
        
        if(in_measure_nlp_time)
        {
            out_cpu_time_of_nlp_solving += *pNLPCpuTime;
            out_clock_time_of_nlp_solving += *pNLPClockTime;
        }
        
        out_number_of_nlp_probs_solved++;
        
        if( !pnlp->feasSol || pnlp->objValue > -epsToIntSol )
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORMSG("Failure to obtain a valid interior solution");
            
            out_return_code = MRQ_NLP_SOLVER_ERROR;
            goto termination;
        }
    }
    
    
    if( in_print_level > 4 )
    {
        double *p = pnlp->sol;
        std::cout << MRQ_PREPRINT  "Strategy: " << (solveNLPIPProb ? MRQ_ESHP_IPS_MOST_INTERIOR : in_interior_point_strategy) << " Interior solution gotten: \n"
        MRQ_PREPRINT "optsolver status: " << pnlp->retCode << " original solver status: " << pnlp->origSolverRetCode << " Objective: " << pnlp->objValue << "\n";
        
        if( in_print_level > 5 )
        {
            for(int i = 0; i < n; i++)
                std::cout << MRQ_PREPRINT "x["<<i<<"]: " << p[i] << "\n";
            
            if( solveNLPIPProb )
                std::cout << MRQ_PREPRINT "aux variable: " << p[n] << "\n";
        }
    }
    
    
    
    MRQ_copyArray(n, pnlp->sol, isol);//do not use getsolution here, because ipNLP has one variable more...
    
    if( nlp )
    {
        delete nlp;
        nlp = NULL;
    }
    if( ipNLP )
    {
        delete ipNLP;
        ipNLP = NULL;
    }
    
    
    
    
    ret = masterMilp.setProblemBase( thnumber, prob, in_milp_solver, true, in_set_quadratics_in_master_problem, in_max_lp_subiters == 0, lx, ux, 1, milpSolverParams, in_number_of_threads );
    MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    
    master = masterMilp.master;
    
    if( run_by_inside )
    {
        if( !std::isinf(insideSolverMaxTime) )
        {
            ret = master->setMaxTime( insideSolverMaxTime );
            MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
    }
    
    
    linearizeObj = prob.hasNlObj || (prob.Q.getNumberOfElements() > 0 && !in_set_quadratics_in_master_problem);
    
    
    
    
    if( !linOnInteriorSol )
    {
        //checking is some variable is unbounded. We need that because if it is true, we have to add linearizations to all constraints to guarantee lp relaxation will not be unbounded in the first iteration
        
        for(int i = 0; i < n; i++)
        {
            if( lx[i] <= -MIP_INFINITY || ux[i] >= MIP_INFINITY )
            {
                linOnInteriorSol = true;
                break;
            }
        }
    }
    
    
    if( linOnInteriorSol )
    {
        const int r = masterMilp.addLinearizedNLConstraints(true, isol, in_set_quadratics_in_master_problem, auxCEval, auxInds, auxVals, auxConstr);
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    }
    
    
    if(linearizeObj)
    {
        const int r = masterMilp.addLinearizedObjFunction( !(linOnInteriorSol && prob.hasNlConstrs), isol, in_set_quadratics_in_master_problem, auxInds, auxVals, NULL );
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    }
    
    
    if( in_max_lp_subiters <= 0 )
    {
        milpLoop = true;
    }
    
    
    while(true)
    {
        //int lastLinearizationIndex = INT_MAX;
        iter++;
        
        
        //master->generateModelFile("eshp.lp");
        //std::cout << "Gerei eshp.lp\n";
        
        //when w = 0, we are in a regular iteration. When w = 1, is because we got a feasible solution before and we are fixing integer variables.
        for(signed char w = 0; w < 2; w++)
        {
            
            if( in_max_cpu_time < INFINITY || in_max_time < INFINITY )
            {
                const double rcputime = in_max_cpu_time - ( ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC);
                const double rtime = in_max_time - (MRQ_getTime() - timeStart);
                int r = 0;
                
                if(rcputime <= 0.0 || rtime <= 0.0)
                {
                    out_return_code = MRQ_MAX_TIME_STOP;
                    goto endLoop;
                }
                
                if(in_max_cpu_time < INFINITY)
                    r = master->setMaxCPUTime(rcputime);
                
                if(in_max_time < INFINITY)
                    r += master->setMaxTime(rtime);
                
                if(r != 0)
                {
                    if(in_print_level > 0)
                        MRQ_PRINTERRORMSG("Error to set MILP maximum time!");
                    //out_return_code = MRQ_MILP_SOLVER_ERROR;
                    //break;
                }
            }
            
            
            ret = master->solve( false );
            
            { //couting number of iterations
                long unsigned int masterIters;
                int r = master->getNumberOfIterations(masterIters);
                if(r == 0)
                    out_number_of_milp_solver_iters += masterIters;
                else if(in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(r);
            }
            
            if( w == 1 ) //we have to unfix the integer variables
                MRQ_unfixIntegerVarsByList(nI, intVars, lx, ux, *master);
            
            
            if( in_print_level > 6 )
            {
                std::cout << MRQ_PREPRINT "w: " << (int) w << " master solving. opstovers code: " << master->retCode << " original solver code: " << master->origSolverRetCode << " obj value: " << master->objValue << "\n";
                
                if( in_print_level > 7 )
                {
                    const double *x = master->sol;
                    
                    for(int i = 0; i < n; i++)
                        std::cout << "x["<<i<<"]: " <<  x[i] << "\n";
                    std::cout << "auxiliary objectuve varable: " << x[n] << "\n";
                }
            }
            
            
            if( ret == OPT_OPTIMAL_SOLUTION )
            {
                if( master->objValue > zl && w == 0)
                    zl = master->objValue;
            }
            else if( (w == 1 && ret == OPT_INFEASIBLE_PROBLEM) )
            {
                //we fix integer variables and problem gets infeasible. it is ok, just let algorithm continue its work. Note, in this case, we need the last linearization (because it avoid this integer solution), and, so, we do not delete it
                break;
            }
            else 
            {
                
                #if MRQ_DEBUG_MODE
                    assert( milpLoop || ret != OPT_INFEASIBLE_PROBLEM ); //we already got an interior solution. So at least the continuous lp problem cannot be infeasible
                #endif
                
                if( ret == OPT_UNBOUNDED_PROBLEM )
                {
                    #if MRQ_DEBUG_MODE
                        assert( iter == 1 );
                    #endif
                    
                    if( in_print_level > 0 )
                        MRQ_PRINTMSG("Unbounded master problem!\n");
                    
                    out_return_code = MRQ_UNBOUNDED_PROBLEM;
                    goto endLoop;
                }
                else if( ret == OPT_MAX_TIME )
                {
                    out_return_code = MRQ_MAX_TIME_STOP;
                    goto endLoop;
                }
                
                if( in_print_level > 0 )
                    std::cerr << MRQ_PREPRINT "Warning: Error at master problem solving. optsolver code: " << ret << " original solver code: " << master->origSolverRetCode << "\n";
                
                if( !master->feasSol )
                {
                    out_return_code = MRQ_MILP_SOLVER_ERROR;
                    goto endLoop;
                }
            }
            
            #if 0
            if(in_delete_linearizations_from_int_vars_fixing && w == 1)
            { //if master is not infeasible when we fix integer vars, we can try remove the last linearization
                
                int nconstrs, nconstrToDel;
                int r = master->getNumberOfConstraints(nconstrs);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
                
                nconstrToDel = nconstrs - lastLinearizationIndex;
                
                #if MRQ_DEBUG_MODE
                {
                    int maxLinearizations = m - prob.getNumberOfLinearConstraints() + 1;
                    
                    if(in_set_quadratics_in_master_problem)
                        maxLinearizations -= prob.getNumberOfQuadConstraints();
                    
                    assert( lastLinearizationIndex > 0  && nconstrToDel > 0  &&  nconstrToDel <= maxLinearizations );
                }
                #endif
                
                #pragma ivdep
                #pragma GCC ivdep
                for(int k = 0; k < nconstrToDel; k++)
                    auxInds[k] = k + lastLinearizationIndex;
                    
                r = master->removeConstraints(nconstrToDel, auxInds);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
            }
            #endif
            
            if(iter > 1 && w == 0 && MRQ_abs(robj - master->objValue) <= MRQ_abs(relTolToTestPreviousSol*robj) +  absTolToTestPreviousSol)
            { //rsol has the solution in the previous iterations
                //we have two subsequenet solutions having the same objective value. Maybe, we are ciclyng because the solution is the same. Let's check:
                
                bool sameSol = true;
                double *csol = master->sol;
                
                for(int i = 0; i < n; i++)
                {
                    if(  MRQ_abs( rsol[i] - csol[i] ) > + MRQ_abs(relTolToTestPreviousSol*rsol[i]) + absTolToTestPreviousSol )
                    {
                        sameSol = false;
                        break;
                    }
                }
                
                if(sameSol)
                {
                    if( in_print_level > 5 )
                        MRQ_PRINTMSG("Same solution from previous iterations obtained. Increasing epsilon to active constraint to linearization\n");
                    
                    epsToActConstrToLin += deltaEpsActLin;
                }
            }
            
            
            if( !in_fix_int_vars_to_try_improve_cut )
                w = 9; //we do not fix integer variables to improve cut
            
            robj = master->objValue;
            MRQ_copyArray(n, (const double *) master->sol, rsol);
            
            ret = prob.isFeasibleToConstraints( thnumber, rsol, true, auxCEval, in_absolute_feasibility_tol, in_relative_feasibility_tol, feasSol, auxConstr );
            MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_CALLBACK_FUNCTION_ERROR, endLoop);
            
            
            if( !milpLoop || (feasSol && !linearizeObj) )
            {
                //if we are not in the milpLoop, we can break the for. If solution is feasible and obj is linear, we can break the for since feasible solution is optimal...
                w = 9;
            }
            else  if(w == 0)
            {
                //so we fix integer variables to solve problem again
                int r = MRQ_fixIntVarsOnSolByList(nI, intVars, rsol, *master);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
                
                /*r = master->getNumberOfConstraints(lastLinearizationIndex);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);*/
                
                if( !feasSol )
                {
                    //for(int i = 0; i < prob.m; i++)
                        //printf("c[%d]: %f\n", i, auxConstr[i]);
                    
                    //if the solution is already feasible, I believe we do not need linearize constarints. Anyway, linearization of objective can move the optimal solution out of feasible region. However, I think it is better avoid linearization in this case 
                    int ret = masterMilp.addLinearizedNLConstraintsByStrategy( epsToActConstrToLin, &out_number_of_constr_linears_saved, false, rsol, in_set_quadratics_in_master_problem, constrLinearisationStrategy, auxCEval, NULL, NULL, NULL, auxConstr, true );
                    MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
                    
                }
                
                if(linearizeObj)
                {
                    const bool newx = !prob.hasNlConstrs;
                    
                    ret = masterMilp.addLinearizedObjFunction(newx, rsol, in_set_quadratics_in_master_problem, auxInds, auxVals, NULL);
                    MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
                }
                
            }
            
            
            if( feasSol )
            {
                if( milpLoop )
                {
                    //so, we have a feasible integer solution at hands...
                    double objValue;
                    
                    int ret = prob.objEval( thnumber, !prob.hasNlConstrs, rsol, objValue );
                    MRQ_IFCALLBACKERRORGOTOLABEL(ret, out_return_code, endLoop);
                    
                    if( objValue < zu )
                    {
                        //updateBestSol( n, objValue, rsol, clockStart, timeStart, iter );
                        tryUpdateBestSolution(thnumber, n, rsol, objValue, iter, clockStart, timeStart, in_store_history_solutions);
                    }
                    
                    //since solution is feasible, even if if in_fix_int_vars_to_try_improve_cut is false, there is no sense in perform linear search because linear search is for infeasible solution. So, we linearize only objective function if it s nonlinear
                    
                    if(w > 0 && linearizeObj) //if w == 0, objective already linearized in the if above. Note if in_fix_int_vars_to_try_improve_cut is true, we have w > 0
                    {
                        ret = masterMilp.addLinearizedObjFunction(true, rsol, in_set_quadratics_in_master_problem, auxInds, auxVals, NULL);
                        MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
                    }
                    
                }
                else
                {
                    bool change = true;
                    //changing to milp loop 
                    //w = 9; //we break the for
                    
                    if( linearizeObj )
                    {
                        //we add some objective linearization
                        double obj;
                        
                        int ret = prob.objEval(thnumber, !prob.hasNlConstrs, rsol, obj );
                        MRQ_IFCALLBACKERRORGOTOLABEL(ret, out_return_code, endLoop);
                        
                        //std::cout << "Solucao viavel! zl: " << zl << " zu_cont: " << zu_cont << " obj:" << obj  << "\n";
                        //MRQ_getchar();
                        
                        ret = masterMilp.addLinearizedObjFunction(false, rsol, in_set_quadratics_in_master_problem, auxInds, auxVals, &obj);
                        MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
                        
                        
                        if( iterlpfeas == 0 )
                            iterlpfeas = iter;
                        
                        if( obj < zu_cont )
                            zu_cont = obj;
                        
                        if( iter - iterlpfeas < in_max_lp_subiter_to_improve_obj_app )
                        {
                            if( zu_cont - zl > in_cont_relax_absolute_convergence_tol && zu_cont - zl > MRQ_abs(zu_cont)* in_cont_relax_relative_convergence_tol )
                                change = false;
                        }
                    }
                    
                    
                    if( change )
                    {
                        ret = MRQ_changeToMILPLoop( nI, intVars, *master );
                        MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
                        
                        std::cout << MRQ_PREPRINT " iter: " << iter << " changing to milp by (nlp) feasible solution. Objective: " << master->objValue << "\n";
                        
                        if( master->retCode == OPT_OPTIMAL_SOLUTION && !linearizeObj ) //we have a linear (or quadratic) objective function. So, if solution is feasible and we have an optimal lp solution, that solutios is solution of constinuous nonlinear relaxation
                            out_obj_opt_at_continuous_relax = master->objValue;
                        
                        constrLinearisationStrategy = in_constr_linearisation_strategy;
                        out_number_of_lp_iterations = iter;
                        milpLoop = true;
                    }
                    
                    //oficial algorithm does not requires it, but maybe it is usefull linearize active constraint in this feasible solution
                    
                    ret = masterMilp.addLinearizedNLConstraintsByStrategy( epsToActConstrToLin, &out_number_of_constr_linears_saved, false, rsol, in_set_quadratics_in_master_problem, constrLinearisationStrategy, auxCEval, NULL, NULL, NULL, auxConstr, true );
                    MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
                    
                }
                
            }
            else if( w > 0)
            {
                //w = 9;//we break the for
                bool newx = true;
                double lslambda;
                
                //performing the line search
                
                ret = MRQ_lineSearch2(thnumber, prob, in_eps_to_line_search, isol, rsol, in_absolute_feasibility_tol, in_relative_feasibility_tol, auxCEval, auxConstr, lsSol, &lslambda );
                MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_UNDEFINED_ERROR, endLoop);
                
                #if 0
                
                //checking if the line search solution is feasible
                if( milpLoop && in_test_linear_search_solution_as_best_solution && lslambda <= 0.001 )
                {
                    bool feasSol;
                    double objValue;
                    double *tsol = auxVals;
                    
                    MRQ_copyArray(n, lsSol, tsol);
                    
                    for(int i = 0; i < nI; i++)
                    {
                        const int ind = intVars[i];
                        tsol[ind] = rsol[ind];
                    }
                    
                    
                    ret = prob.objEval(thnumber, true, tsol, objValue);
                    if( ret != 0 )
                    {
                        if( in_print_level > 0 )
                            MRQ_PRINTERRORNUMBER(ret);
                        
                        out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
                        goto endLoop;
                    }
                    
                    
                    //for(int i = 0; i < n; i++)
                        //std::cout << "tsol["<<i<<"]: " << tsol[i] << "\n";
                    std::cout << "obj: " << objValue << " zl: " << zl <<  " zu: " << zu << "\n";
                    //MRQ_getchar();
                    
                    
                    
                    if( in_test_linear_search_solution_as_best_solution && milpLoop && objValue < zu ) // && MRQ_isIntegerSol( nI, intVars, lsSol, in_integer_tol) )
                    {
                        const bool newSol = !prob.hasNlObj;
                        
                        ret = prob.isFeasibleToConstraints(thnumber, tsol, newSol, auxCEval, in_absolute_feasibility_tol, in_relative_feasibility_tol, feasSol, auxConstr);
                        
                        if( ret != 0 )
                        {
                            if( in_print_level > 0 )
                                MRQ_PRINTERRORNUMBER(ret);
                            
                            out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
                            goto endLoop;
                        }
                        
                        for(int i = 0; i < prob.m; i++)
                        {
                            if( !auxCEval[i] )
                                continue;
                            
                            std::cout << "i: " << i << " lc: " << prob.lc[i] <<  " c: "<< auxConstr[i] << " uc: " << prob.uc[i] << "\n";
                        }
                        
                        std::cout << "feas: " << feasSol << "\n";
                        MRQ_getchar();
                        
                        if( feasSol )
                        {
                            //std::cout << "atualizei best sol. antigo zu: " << zu << " novo zu: " << objValue << "\n";
                            
                            tryUpdateBestSolution(thnumber, n, tsol, objValue, iter, clockStart, timeStart, in_store_history_solutions);
                            
                            MRQ_getchar();
                        }
                        
                    }
                    
                }
                
                #endif
                
                
                if(linearizeObj)
                {
                    //we evaluate constraints on linear search
                    ret = masterMilp.addLinearizedObjFunction( true, lsSol, in_set_quadratics_in_master_problem, auxInds, auxVals, NULL );
                    MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
                
                    newx = !prob.hasNlObj;
                }
                
                ret = masterMilp.addLinearizedNLConstraintsByStrategy( epsToActConstrToLin, &out_number_of_constr_linears_saved, newx, lsSol, in_set_quadratics_in_master_problem, constrLinearisationStrategy, auxCEval, NULL, NULL, NULL, NULL, false );
                MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
                
            }
            
            
            if( in_print_level > 1 && (iter-1)%in_printing_frequency == 0 )
            {
                printf(MRQ_PREPRINT "%lu %s  %+-14e  %+-14e  %+-14e\n", iter, milpLoop ? "milp" : "lp",  zl, zu, zu-zl );
                //MRQ_getchar();
            }
            
            
            if( checkTerminationCriterions(thnumber, zl, zu, iter, timeStart, clockStart, out_return_code) )
            {
                goto endLoop;
            }
        }
        
        
        if( !milpLoop && iter >= in_max_lp_subiters )
        {
            milpLoop = true;
            ret = MRQ_changeToMILPLoop(nI, intVars, *master);
            MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, endLoop);
            
            out_number_of_lp_iterations = iter;
        }
        
    }
    
    
endLoop:
    
    if(in_print_level > 0)
    {
        std::cout << MRQ_PREPRINT;
        if(out_return_code == MRQ_OPTIMAL_SOLUTION)
            std::cout << "Optimal solution found. Objective: " << out_best_obj << "\n";
        else if(out_return_code == MRQ_INFEASIBLE_PROBLEM)
            std::cout << "Infeasible problem.\n";
        else if(out_return_code == MRQ_MAX_ITERATIONS_STOP )
            std::cout << "Max iteration number reached!\n";
        else if(out_return_code == MRQ_MAX_TIME_STOP)
            std::cout << "Max time reached!\n";
    }
    
    
    
    
termination:
    
    if(ipNLP)
        delete ipNLP;
    
    if(plc)			free(plc);
    if(intVars)		free(intVars);
    if(isol)		free(isol);
    if(auxInds)		free(auxInds);
    if(auxVals)		free(auxVals);
    if(auxCEval)	free(auxCEval);
    if(auxConstr)	free(auxConstr);
    
    
    if( !milpLoop )
        out_number_of_lp_iterations = iter; //we did not ch
    
    
    out_number_of_iterations = iter;
    
    
    algorithmFinalization(1, prob, lx, ux);
    
    out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
    out_clock_time = MRQ_getTime() - timeStart;
    
    if(in_print_level > 1)
        std::cout << MRQ_PREPRINT "cpu time: " << out_cpu_time << "\n";
    
    return out_return_code;
}





/*bool MRQ_ExtSupHypPlane::updateBestSol( const int n, double objValue, double *sol, const clock_t& clockStart, const int timeStart, const long unsigned int iter )
{
    bool updt = true;
    double time, cpuTime;
    
    
    if( in_store_history_solutions )
    {
        cpuTime = MRQ_calcCPUTtime(clockStart, clock());
        time = MRQ_getTime() - timeStart;
    }
    
    
    if( in_user_callbacks )
    {
        const int r = in_user_callbacks->updatingBestSolution( out_algorithm, thnumber, sol, objValue, zu, iter );
        
        if(r != 0)
        {
            if( in_print_level > 0 )
                std::cerr << MRQ_PREPRINT "Callback function error " << r << "on in_user_callbacks->updatingBestSolution" << MRQ_GETFILELINE << "\n";
            
            updt = false;
        }
    }
    
    //objValue can have been changed
    if( objValue >= zu )
        updt = false;
    
    
    if( updt )
    {
        out_best_obj = objValue;
        MRQ_copyArray(n, sol, out_best_sol);
        
        if( in_store_history_solutions )
        {
            
            const int r = out_sol_hist.addSolution(n, iter, time, cpuTime, sol, objValue);
            
            if( r != 0 )
            {
                if( in_print_level > 0 )
                    std::cerr << MRQ_PREPRINT "Warning: Error " << r << " on out_sol_hist.addSolution." << MRQ_GETFILELINE << "\n";
            }
        }
    }
    
    return updt;
} */


