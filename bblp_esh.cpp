/*That file contains a simple implementation of
* LP-BB extended cutting plane based algorithm
*
* References:
*
* 1 - Kronqvist & Lundell & Westerlund, The extended supporting hyperplane algorithm for convex mixed-integer nonlinear programming. Journal of Global Optimization 64 (2016), pages 249-272.
* 
* 2 - Lundell & Kronqvist & Westerlund, Improvements to the Supporting Hyperplane Optimization Toolkit Solver for Convex MINLP, Proceedings of the XIII Global Optimization Workshop, GOW'16, 4-8 September 2016.
* 
* 3 - My head.
*
* Author: Wendel Melo
*
* Date: 10-June-2018
* 
* 
* * Reference to set cplex callback before solve: admipex1.c
* 
*/


#include <cstdlib>
#include <cmath>
#include <new>


#include "BBL_tools.hpp"

#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_milpCallbacks.hpp"


using namespace optsolvers;
using namespace muriqui;




MRQ_LPBBExtSupHypPlane::MRQ_LPBBExtSupHypPlane():MRQ_ExtSupHypPlane()
{
    resetParameters();
    out_algorithm = MRQ_LP_BB_ESH_BASED_ALG;
}


MRQ_LPBBExtSupHypPlane::~MRQ_LPBBExtSupHypPlane()
{
}


int MRQ_LPBBExtSupHypPlane::checkAlgorithmRequirements( MRQ_MINLPProb &prob, const double *lx, const double *ux)
{
    if(!MRQ_isLazyConstraintsAvaliable(in_milp_solver) )
    {
        //char solverName[30];
        //MRQ_enumToStr( in_milp_solver, solverName );
        
        const std::string &solverName = OPT_getSolverName(in_milp_solver);
        
        std::cerr << MRQ_PREPRINT "We are so sorry. Algorithm " <<  getAlgorithmName() << " (" << out_algorithm << ") does not work with milp solver " << solverName  << " (" << in_milp_solver << ")  \n";
        
        return MRQ_BAD_PARAMETER_VALUES;
    }
    
    return MRQ_ExtSupHypPlane:: checkAlgorithmRequirements(prob, lx, ux);
}


void MRQ_LPBBExtSupHypPlane::resetOutput() 
{
    out_number_of_master_relaxation_solved = 0;
}


void MRQ_LPBBExtSupHypPlane::resetParameters() 
{
    MRQ_ExtSupHypPlane::resetParameters();
    in_fix_int_vars_to_try_improve_cut = false;
    in_max_lp_subiters = 0;
}


int MRQ_LPBBExtSupHypPlane:: solverCallbackLazyConstraints( MRQ_MILPCallbackData &data)
{
    MRQ_MINLPProb &prob = *data.prob;
    
    const bool linearizeObj = data.linearizeObj;
    const int n = prob.n;
    const int thnumber = data.thnumber;
    const bool incQuadsInMaster = data.setQuadsInMaster;
    
    bool feasSol, newx = true;
    int r, retCode;
    long int iter = -1;
    double objValue = NAN;
    
    
    
    const int *indices = data.indices;
    const bool *constrEval = data.constrEval;
    
    bool *auxConstrEval2 = data.auxConstrEval2;
    double *auxVars = data.auxVars;
    double *masterSol = data.auxVars2;
    double *lineSearchSol = data.auxVars3;
    double *masterConstrValues = data.auxConstr;
    double *constrValues = data.auxConstr2;
    double *interiorSol = data.refSol;
    double *plc = data.plc, *puc = data.puc;
    
    MRQ_GradientsEvaluation &gradEval = data.gradEval;
    MRQ_MILPSolverCallbackInterface *callbackSolver = data.callbackSolver;
    
    
    
    r = callbackSolver->getNodeSolution(n+1, masterSol);
    MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_MILP_SOLVER_ERROR, termination);
    
    r = prob.isFeasibleToConstraints(thnumber, masterSol, true, constrEval,  in_absolute_feasibility_tol, in_relative_feasibility_tol, feasSol, masterConstrValues);
    MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_CALLBACK_FUNCTION_ERROR, termination);
    
    if( prob.hasNlConstrs )
        newx = false;
    
    
    if(feasSol || linearizeObj)
    {
        r = prob.objEval(thnumber, newx, masterSol, objValue);
        MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_CALLBACK_FUNCTION_ERROR, termination);
        
        if( prob.hasNlObj )
            newx = false;
    }
    
    if(feasSol && objValue < zu)
    {
        callbackSolver->getNumberOfIterations(iter);
        
        SEMAPH_updtSol.lock(nthreads_lazy);
        {
            tryUpdateBestSolution(thnumber, n, masterSol, objValue, iter, data.clockStart, data.timeStart, false);
        }
        SEMAPH_updtSol.unlock(nthreads_lazy);
        
        if( zl >= zu )
        { //we could, but, by now, w do not considere optimal tolerances here
            if(in_print_level > 2)
                MRQ_PRINTMSG("A feasible solution better than lower bound was found. stopping");
            
            data.out_sol_lower_than_zl = true;
            
            retCode = MRQ_LAZY_MILP_BB_SOLUTION_LOWER_THAN_ZL;
            goto termination;
        }
    }
    
    if( feasSol )
    { //if solution is feasible, we only linearize objetive function since there is no sense in performing linear search (linear search is just to find a feasible solution in the edge of feasible reagion from an infeasible solution)
        
        if(linearizeObj)
        {
            r = MRQ_addLazyConstraintObjectiveLinearizationOnSolution(thnumber, callbackSolver, prob, incQuadsInMaster, newx, masterSol, &objValue, indices, auxVars);
            MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
        }
    }
    else
    {
        newx = true;
        double lslambda;
        
        r = MRQ_lineSearch2(thnumber, prob, in_eps_to_line_search, interiorSol, masterSol, in_absolute_feasibility_tol, in_relative_feasibility_tol, constrEval, constrValues, lineSearchSol, &lslambda );
        
        
        
        //we need calculate constraint values. So, we taking advantage and use isFeasibleToConstraints. Warning: probabily, thi solution is not integer. So, we do not try update best solution...
        r = prob.isFeasibleToConstraints(thnumber, lineSearchSol, true, constrEval,  in_absolute_feasibility_tol, in_relative_feasibility_tol, feasSol, constrValues);
        MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        
        r = addLazyConstraintsLinearizationOnSolution( thnumber, callbackSolver, prob, gradEval, incQuadsInMaster, linearizeObj, constrEval, lineSearchSol, NULL, constrValues, masterConstrValues, indices, plc, puc, auxVars, auxConstrEval2 );
        MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
    }
    
    
    
    
    retCode = 0;
    
termination:
    
    return retCode;
}


int MRQ_LPBBExtSupHypPlane::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    const double timeStart = MRQ_getTime();
    const clock_t clockStart = clock();
    
    const bool binProblem = prob.isBinaryProblem();
    const int n = prob.n;
    const int m = prob.m;
    //const int nI = prob.getNumberOfIntegerVars();
    const bool preproc = in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_obj_function;
    
    const bool setQuadsInMaster = in_set_quadratics_in_master_problem;
    const bool linearizeObj = prob.hasObjNLTerm() || (prob.Q.getNumberOfElements() > 0 && !setQuadsInMaster);
    
    bool updtConstrBounds;
    bool linOnInteriorSol = in_linearize_on_interior_sol;
    bool solveNLPIPProb = in_interior_point_strategy == muriqui::MRQ_ESHP_IPS_MOST_INTERIOR;
    int r;
    
    bool *constrEval = NULL, *auxConstrEval; 
    //int *intVars = NULL;
    int *indices = NULL;
    double *constrValues = NULL;
    
    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux;
    double *plc = NULL, *puc = NULL;
    double *interiorSol = NULL;
    
    MRQ_NLPSolver *nlp = NULL;
    MRQ_NLPIPProblem *ipNLP = NULL; //interior point problem
    
    
    MRQ_Preprocessor preprocessor(&prob);
    MRQ_LAAPointsStoring laps(n);
    MRQ_GradientsEvaluation gradEval;
    
    MRQ_LPSolver *master;
    MRQ_MasterMILPProb masterMilp;
    
    MRQ_MILPCallbackData *callbackData = NULL;
    MRQ_PointsStoreKeeper *pointsKeeper = NULL;
    MRQ_Mutex SEMAPH_addLazy;
    
    
    nthreads = in_number_of_threads > 0 ? in_number_of_threads : branchAndBound::BBL_getNumCores() ;
    
    nthreads_lazy = nthreads;
    if( in_milp_solver == MRQ_GUROBI )
        nthreads_lazy = 1; //gurobi apply multithreading to solve the problem, but only thread 0 is called to add lazy constraints.
    
    {
        auto ret = algorithmInitialization(nthreads_lazy, preproc, milpSolverParams, nlpSolverParams, prob, lx, ux, &preprocessor, &updtConstrBounds, &plc, &puc);
        
        if(ret != MRQ_SUCCESS)
        {
            if(in_print_level > 0)
            {
                if(ret == MRQ_INFEASIBLE_PROBLEM)
                    std::cout << MRQ_PREPRINT << "Preprocessor detected infeasible problem\n";
                else
                    MRQ_PRINTERRORMSG("Error at algorithm initialization\n");
            }
            
            out_return_code = ret;
            goto termination;
        }
    }
    
    
    if(in_print_level > 1)
    {
        std::cout << "\n";
        MRQ_PRINTMSG("Starting LP Branch-And-Bound based on Extended Supporting Hyperplane\n\n");
    }
    
    
    if(in_print_level > 1)
        printSubSolvers(true, true, false);
    
    //MRQ_malloc(intVars, nI);
    MRQ_malloc(indices, n+1);
    MRQ_malloc(constrEval, 2*m);
    MRQ_malloc(constrValues, m);
    MRQ_malloc(interiorSol, n+1);
    callbackData = new (std::nothrow) MRQ_MILPCallbackData[nthreads_lazy];
    MRQ_IFMEMERRORGOTOLABEL(!indices || !constrEval || !constrValues || !interiorSol || !callbackData, out_return_code, termination);
    
    
    auxConstrEval = &constrEval[m];
    
    //prob.getIntegerIndices(intVars);
    
    {
        const bool *nlConstr = prob.nlConstr;
        const minlpproblem::MIP_SparseMatrix *QC = prob.QC;
        
        #pragma ivdep
        #pragma GCC ivdep
        for(int i = 0; i < m; i++)
            constrEval[i] = nlConstr[i] || (QC[i].getNumberOfElements() > 0);
    }
    
    if( in_fix_int_vars_to_try_improve_cut )
    {
        pointsKeeper = new (std::nothrow) MRQ_PointsStoreKeeper[nthreads_lazy];
        MRQ_IFMEMERRORGOTOLABEL(!pointsKeeper, out_return_code, termination);
        
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < nthreads_lazy; i++)
            pointsKeeper[i].nthreads = nthreads_lazy;
    }
    
    for(unsigned int i = 0; i < nthreads_lazy; i++)
    {
        const bool useMilpRelax = false; // in_fix_int_vars_to_try_improve_cut;
        
        
        int r = callbackData[i].allocateBase(thnumber + i, nthreads_lazy, &prob, in_milp_solver, useMilpRelax, NULL, &SEMAPH_addLazy);
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MEMORY_ERROR, termination);
        
        callbackData[i].binProblem = binProblem;
        callbackData[i].linearizeObj = linearizeObj;
        callbackData[i].setQuadsInMaster = setQuadsInMaster;
        
        callbackData[i].timeStart = timeStart;
        callbackData[i].clockStart = clockStart;
        
        //callbackData[i].nI = nI;
        //callbackData[i].intVars = intVars; //actually, we do not need ontVars by now in the callbacks...
        callbackData[i].indices = indices;
        callbackData[i].constrEval = constrEval;
        
        callbackData[i].refSol = interiorSol;
        
        if(plc)
        {
            callbackData[i].plc = plc;
            callbackData[i].puc = puc;
        }
        else
        {
            callbackData[i].plc = prob.lc;
            callbackData[i].puc = prob.uc;
        }
        
        callbackData[i].alg = this;
        callbackData[i].laps = &laps;
        callbackData[i].allPointsKeepers = pointsKeeper;
        
        #if 0
        if(useMilpRelax)
        {
            /**** Setting master's LP relaxation *****/
            MRQ_MasterMILPProb*  masterMilp;
            MRQ_LPSolver *lp;
            
            
            masterMilp = callbackData[i].relaxMaster;
            
            //r = masterMilp->setProblemBase(thnumber, prob, in_milp_solver, true, setQuadsInMaster, true, lx, ux, 1, milpSolverParams); 
            
            r = setMasterProblemBase(masterMilp, thnumber, prob, in_milp_solver, true, setQuadsInMaster, true, lx, ux, 1, milpSolverParams, constrEval, NULL, 1); //note, we set a milp problems actually, but there is no problem about that because when we solve this problem, integer variable will be fixed
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
        #endif
    }
    
    
    //setar o problema de solução interior, e setar a variavel initSol abaixo (na linha com um TODO) com a solução...
    if( in_interior_point_strategy == MRQ_ESHP_IPS_CLOSED_TO_CONT_RELAX_SOL )
    {
        double NLPCpuTime, NLPClockTime;
        double *pNLPCpuTime = in_measure_nlp_time ? &NLPCpuTime : NULL;
        double *pNLPClockTime = in_measure_nlp_time ? &NLPClockTime : NULL;
        
        nlp = OPT_newNLPSolver(in_nlp_solver);
        MRQ_IFMEMERRORGOTOLABEL(!nlp, out_return_code, termination);
        
        r = MRQ_setInteriorPointProbClosedToNLPRelaxProb(prob, lx, ux, plc, puc, nlp, thnumber, in_set_special_nlp_solver_params, nlpSolverParams, 1, in_max_cpu_time, in_max_time, !in_set_quadratics_in_master_problem, in_eps_to_enforce_interior_sol_on_cont_relax_sol);
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
        
        if( run_by_inside )
        {
            if( !std::isinf(insideSolverMaxTime) )
            {
                int r = nlp->setMaxTime(insideSolverMaxTime);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
            }
        }
        
        nlp->solveAndGetTime(pNLPCpuTime, pNLPClockTime, false);
        
        if(in_measure_nlp_time)
        {
            out_cpu_time_of_nlp_solving += *pNLPCpuTime;
            out_clock_time_of_nlp_solving += *pNLPClockTime;
        }
        
        out_number_of_nlp_probs_solved++;
        
        if( nlp->feasSol )
        {
            MRQ_copyArray(n, nlp->sol, interiorSol);
        }
        else
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
        
        
        delete nlp;
        nlp = NULL;
    }
    
    
    if( solveNLPIPProb )
    {
        double NLPCpuTime, NLPClockTime;
        double *pNLPCpuTime = in_measure_nlp_time ? &NLPCpuTime : NULL;
        double *pNLPClockTime = in_measure_nlp_time ? &NLPClockTime : NULL;
        optsolvers::OPT_LPSolver *pnlp;
        
        
        ipNLP = new (std::nothrow) MRQ_NLPIPProblem;
        MRQ_IFMEMERRORGOTOLABEL(!ipNLP, out_return_code, termination);
        
        r = ipNLP->setProblem( in_nlp_solver, prob , lx, ux, nlpSolverParams, thnumber, in_set_special_nlp_solver_params, !in_set_quadratics_in_master_problem, in_number_of_threads, in_max_cpu_time, in_max_time);
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
        
        pnlp = ipNLP->solver;
        
        if( run_by_inside )
        {
            if( !std::isinf(insideSolverMaxTime) )
            {
                r = pnlp->setMaxTime( insideSolverMaxTime );
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
        
        if( !pnlp->feasSol || pnlp->objValue > -in_eps_to_interior_sol )
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORMSG("Failure to obtain a valid interior solution");
            
            out_return_code = MRQ_NLP_SOLVER_ERROR;
            goto termination;
        }
        
        
        MRQ_copyArray(n+1, pnlp->sol, interiorSol);//do not use getsolution here, because ipNLP has one variable more...
        
        delete ipNLP;
        ipNLP = NULL;
    }
    
    
    if( in_print_level > 4 )
    {
        double *p = interiorSol;
        std::cout << MRQ_PREPRINT  "Strategy: " << (solveNLPIPProb ? MRQ_ESHP_IPS_MOST_INTERIOR : in_interior_point_strategy) << " Interior solution gotten: \n";
        
        if( in_print_level > 5 )
        {
            for(int i = 0; i < n; i++)
                std::cout << MRQ_PREPRINT "x["<<i<<"]: " << p[i] << "\n";
            
            if( solveNLPIPProb )
                std::cout << MRQ_PREPRINT "aux variable: " << p[n] << "\n";
        }
    }
    
    
    
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
    
    
    
    /***** Setting master problem *****/
    if( prob.getProblemType() != minlpproblem::MIP_PT_MILP )
    {
        r = setMasterProblemBase(&masterMilp, thnumber, prob, in_milp_solver, true, setQuadsInMaster, in_max_lp_subiters == 0, lx, ux, 1, milpSolverParams, auxConstrEval, &laps, nthreads);
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    }
    
    master = masterMilp.master;
    
    if(nPoints == 0 || linOnInteriorSol )
    {
        //we have to improvise the first linearization point. We will use the interiorSol
        double objValue;
        double * const initSol = interiorSol;
        
        int r = prob.objEval(thnumber, true, initSol, objValue);
        MRQ_IFCALLBACKERRORGOTOLABEL(r, out_return_code, termination);
        
        r = prob.constraintsEval(thnumber, !prob.hasNlObj, constrEval, initSol, constrValues );
        MRQ_IFCALLBACKERRORGOTOLABEL(r, out_return_code, termination);
        
        r = masterMilp.addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, false, initSol, setQuadsInMaster, in_constr_linearisation_strategy, constrEval, NULL, NULL, NULL, constrValues, true );
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        
        
        if(linearizeObj)
        {
            int r = masterMilp.addLinearizedObjFunction( false, initSol, setQuadsInMaster, in_obj_linearisation_strategy, zu, &laps, NULL, NULL, &objValue);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
        
        //now we set this first linearization for relax master in threads
        #if 0
        if( in_fix_int_vars_to_try_improve_cut )
        {
            for(unsigned int i = 0; i < nthreads_lazy; i++)
            {
                MRQ_MasterMILPProb *masterMilp = callbackData[i].relaxMaster;
                
                int r = masterMilp->addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, false, initSol, setQuadsInMaster, in_constr_linearisation_strategy, constrEval, NULL, NULL, NULL, constrValues, true );
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                if(linearizeObj)
                {
                    int r = masterMilp->addLinearizedObjFunction( false, initSol, setQuadsInMaster, in_obj_linearisation_strategy, zu, &laps, NULL, NULL, &objValue);
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                }
            }
            
        }		
        #endif
        
    }
    
    
    if( in_max_lp_subiters > 0 )
    {
        bool feasSol;
        unsigned int iterlpfeas = UINT_MAX;
        double zulp = MRQ_INFINITY, objValue; //double lslambda;
        double *pobjValue;
        const double *psol;
        
        double *lineSearchSol = callbackData[0].auxVars3;
        double *masterConstrValues = callbackData[0].auxConstr;
        double *constrValues = callbackData[0].auxConstr2;
        double *auxVars = callbackData[0].auxVars;
        
        
        for(decltype(in_max_lp_subiters) k = 0; k < in_max_lp_subiters; k++)
        {
            feasSol = false;
            r = master->solve(false);
            
            //std::cout << MRQ_PREPRINT "k: " << (int) k << " master solving. opstovers code: " << master->retCode << " original solver code: " << master->origSolverRetCode << " obj value: " << master->objValue << "\n";
            
            const double *masterSol = master->sol;
            
            if( r == OPT_OPTIMAL_SOLUTION )
            {
                if(master->objValue > zl)
                    zl = master->objValue;
            }
            else if( r == OPT_INFEASIBLE_PROBLEM )
            {
                if( in_print_level > 0 )
                    MRQ_PRINTMSG("Infeasible master problem!\n");
                out_return_code = MRQ_INFEASIBLE_PROBLEM;
                goto termination;
            }
            else if( r == OPT_UNBOUNDED_PROBLEM )
            {
                #if MRQ_DEBUG_MODE
                    assert( k == 1 );
                #endif
                
                if( in_print_level > 0 )
                    MRQ_PRINTMSG("Unbounded master problem!\n");
                
                out_return_code = MRQ_UNBOUNDED_PROBLEM;
                goto termination;
            }
            
            if( !master->feasSol )
            {
                out_return_code = MRQ_MILP_SOLVER_ERROR;
                goto termination;
            }
            
            r = prob.isFeasibleToConstraints(thnumber, masterSol, true, constrEval, in_absolute_feasibility_tol, in_relative_feasibility_tol, feasSol, masterConstrValues);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_CALLBACK_FUNCTION_ERROR, termination);
            
            
            if( feasSol )
            {
                const bool newx = !prob.hasNlConstrs;
                psol = masterSol;
                
                r = prob.objEval(thnumber, newx, psol, objValue);
                MRQ_IFCALLBACKERRORGOTOLABEL(r, out_return_code, termination);
                
                pobjValue = &objValue;
                
                if( objValue < zulp )
                    zulp = objValue;
                
                if( iterlpfeas == UINT_MAX )
                    iterlpfeas = k;
            }
            else
            {
                r = MRQ_lineSearch2(thnumber, prob, in_eps_to_line_search, interiorSol, masterSol, in_absolute_feasibility_tol, in_relative_feasibility_tol, constrEval, constrValues, lineSearchSol, NULL);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_UNDEFINED_ERROR, termination);
                
                psol = lineSearchSol;
                pobjValue = NULL;
                
                r = masterMilp.addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, true, psol, setQuadsInMaster, in_constr_linearisation_strategy, constrEval, auxConstrEval, indices, auxVars, NULL, false );
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
            }
            
            if(linearizeObj)
            {
                r = masterMilp.addLinearizedObjFunction(true, psol, setQuadsInMaster, indices, auxVars, pobjValue);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
            }
            
            
            if( in_print_level > 2 && k%10 == 0 )
            {
                std::cout << MRQ_PREPRINT "lp iter: " << k << " zl: " << zl << " lp zu: " << zulp << " lp master obj value: " << master->objValue << "\n";
            }
            
            if(zulp - zl <= in_absolute_convergence_tol || zulp - zl <= MRQ_abs(zulp)*in_relative_convergence_tol)
            {
                std::cout << MRQ_PREPRINT "Optimal solution of lp relaxtion found. Objective: " << zulp << "\n";
                break;
            }
            else if( iterlpfeas != UINT_MAX &&  k - iterlpfeas >= in_max_lp_subiter_to_improve_obj_app)
            {
                break;
            }
        }
        
        //setting integer variables. 
        {
            auto *xtype = prob.xtype;
            for(int i = 0; i < n; i++)
            {	
                if( minlpproblem::MIP_isIntegerType(xtype[i]) )
                {
                    r = master->setVariableType( i, optsolvers::OPT_VT_INTEGER );
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                }
            }
        }
        
    }
    
    
    //do not move this to before lp loop because indices is being used...
    #pragma ivdep
    #pragma GCC ivdep
    for(int i = 0; i <= n; i++) //considering auxiliary variables also
        indices[i] = i;
    
    
    /***** setting up solver to adopt lazy constraints *****/
    r = setLazyConstraintCallbackOnMilpSolver(master, callbackData);
    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    
    if( in_user_callbacks )
    {
        if(in_call_before_solve_callback_in_milp_bb)
        {
            r = setBeforeSolveCallbackOnMilpSolver(master, callbackData);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
        
        if(in_call_branching_callback_in_milp_bb)
        {
            r = setBranchingCallbackOnMilpSolver(master, callbackData);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
    }
    
    r = master->setRelativeOptimalityTol( in_relative_convergence_tol );
    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    
    
    r = master->solve();
    
    zl = master->getDualObjValue();
    
    if( master->feasSol )
    {
        if( master->objValue < MRQ_zuWithTol(out_best_obj, in_absolute_convergence_tol,  in_relative_convergence_tol )   ||  out_best_obj >= MRQ_INFINITY  )
        {
            MRQ_PRINTMSG("Warning: milp solver found a solution slightly infeasible\n");
            
            MRQ_copyArray(n, (const double*) master->sol, out_best_sol );
            
            int r = prob.objEval( thnumber, true, out_best_sol, out_best_obj );
            MRQ_IFCALLBACKERRORGOTOLABEL(r, out_return_code, termination);
        }
    }
    
    
    if(r == OPT_OPTIMAL_SOLUTION)
    {
        //zl = master->getObjValue();
        out_return_code = MRQ_OPTIMAL_SOLUTION;
    }
    else if(r == OPT_MAX_TIME)
    {
        out_return_code = MRQ_MAX_TIME_STOP;
    }
    else if(r == OPT_INFEASIBLE_PROBLEM)
    {
        out_return_code = MRQ_INFEASIBLE_PROBLEM;
    }
    else if(r == OPT_CALLBACK_FUNCTION_ERROR || r == MRQ_CALLBACK_FUNCTION_ERROR)
    {
        out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
    }
    else if(r == OPT_MAX_ITERATIONS)
    {
        out_return_code = MRQ_MAX_ITERATIONS_STOP;
    }
    else if(r == OPT_UNBOUNDED_PROBLEM)
    {
        out_return_code = MRQ_UNBOUNDED_PROBLEM;
    }
    else if(r == MRQ_NLP_SOLVER_ERROR)
    {
        out_return_code = MRQ_NLP_SOLVER_ERROR;
    }
    else
    {
        out_return_code = MRQ_UNDEFINED_ERROR;
    }
    
    
    for(unsigned int i = 0; i < nthreads; i++)
    {
        if( callbackData[i].out_sol_lower_than_zl )
        {
            out_return_code = MRQ_OPTIMAL_SOLUTION;
            break;
        }
    }
    
    
    r = master->getNumberOfIterations(out_number_of_iterations);
    if(r != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        out_number_of_iterations = ULONG_MAX;
    }
    
    
    
termination:
    
    if(plc)			free(plc);
    
    if(constrEval)	free(constrEval);
    //if(intVars)		free(intVars);
    if(indices)		free(indices);
    if(constrValues) free(constrValues);
    if(interiorSol)	free(interiorSol);
    
    if(nlp)		delete nlp;
    if(ipNLP)	delete ipNLP;
    
    if(callbackData)	delete[] callbackData;
    if(pointsKeeper)	delete[] pointsKeeper;
    
    algorithmFinalization(nthreads, prob, lx, ux);
    
    out_number_of_threads = nthreads_lazy;
    out_number_of_milp_solver_iters = out_number_of_iterations;
    out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
    out_clock_time = MRQ_getTime() - timeStart;
    
    if(in_print_level > 1)
        std::cout << MRQ_PREPRINT "cpu time: " << out_cpu_time << "\n";
    
    
    return out_return_code;
}




