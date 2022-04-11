/*That file contains a simple implementation of
* LP-BB extended cutting plane based algorithm
*
* References:
*
* My head
*
* Author: Wendel Melo
*
* Date: 10-March-2018
* 
* 
* Reference to set cplex callback before solve: admipex1.c
* 
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
using namespace minlpproblem;
using namespace muriqui;


MRQ_LPBBExtCutPlan::MRQ_LPBBExtCutPlan(): MRQ_ExtCutPlan()
{
    resetParameters();
    //resetOutput();
    out_algorithm = MRQ_LP_BB_ECP_BASED_ALG;
}


MRQ_LPBBExtCutPlan::~MRQ_LPBBExtCutPlan()
{
}


int MRQ_LPBBExtCutPlan::checkAlgorithmRequirements( MRQ_MINLPProb &prob, const double *lx, const double *ux)
{
    if( !MRQ_isLazyConstraintsAvaliable(in_milp_solver) ) 
    {
        //char solverName[30];
        //MRQ_enumToStr( in_milp_solver, solverName );
        
        const std::string &solverName = OPT_getSolverName(in_milp_solver);
        
        std::cerr << MRQ_PREPRINT "We are so sorry. Algorithm " <<  getAlgorithmName() << " (" << out_algorithm << ") does not work with milp solver " << solverName  << " (" << in_milp_solver << ") \n";
        
        return MRQ_BAD_PARAMETER_VALUES;
    }
    
    return MRQ_ExtCutPlan:: checkAlgorithmRequirements(prob, lx, ux);
}


/*void MRQ_LPBBExtCutPlan::resetOutput() 
{
    MRQ_ExtCutPlan::resetOutput();
}*/


void MRQ_LPBBExtCutPlan::resetParameters()
{
    MRQ_ExtCutPlan::resetParameters();
    //in_fix_int_vars_to_try_improve_cut = false;
    in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 0;
}


int MRQ_LPBBExtCutPlan:: solverCallbackLazyConstraints( MRQ_MILPCallbackData &data)
{
    MRQ_MINLPProb &prob = *data.prob;
    
    const bool linearizeObj = data.linearizeObj;
    const int n = prob.n;
    const int thnumber = data.thnumber;
    const unsigned threadShift = thnumber - this->thnumber;
    
    const bool incQuadsInMaster = data.setQuadsInMaster;
    
    bool feasSol, solveLpFix = in_max_number_of_fixed_relax_master_problem_solved_per_iteration > 0;
    bool addLinearizationsOnMasterSol = true;
    int r, retCode;
    long int iter = -1;
    double objValue = NAN;
    double zlpy, zupy;
    
    double *masterSol = data.auxVars2;
    
    const int *indices = data.indices;
    const bool *constrEval = data.constrEval;
    
    bool *auxConstrEval2 = data.auxConstrEval2;
    
    double *auxVars = data.auxVars;
    double *constrValues = data.auxConstr;
    double *plc = data.plc, *puc = data.puc;
    
    MRQ_GradientsEvaluation &gradEval = data.gradEval;
    MRQ_MILPSolverCallbackInterface *callbackSolver = data.callbackSolver;
    
    MRQ_PointsStoreKeeper *allPointsKeepers = data.allPointsKeepers;
    
    
    
    
    r = callbackSolver->getNodeSolution(n+1, masterSol);
    MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_MILP_SOLVER_ERROR, termination);
    
    
    
    r = prob.isFeasibleToConstraints(thnumber, masterSol, true, constrEval,  in_absolute_feasibility_tol, in_relative_feasibility_tol, feasSol, constrValues);
    MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_CALLBACK_FUNCTION_ERROR, termination);
    
    
    if(feasSol || linearizeObj)
    {
        r = prob.objEval(thnumber, !prob.hasNlConstrs, masterSol, objValue);
        MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_CALLBACK_FUNCTION_ERROR, termination);
    }
    
    //std::cout << "entrei na callback. Obj: " << objValue << " Feas: " << feasSol << "\n";
    
    
    if(feasSol)
    {
        const double objRelax = masterSol[n];
        
        zupy = objValue; //we already have a feasible solution to Py
        
        if( MRQ_abs(objValue - objRelax) <= MRQ_max(in_absolute_convergence_tol, MRQ_abs(objValue*in_relative_convergence_tol) )  )
        {
            //so, the solution of master problem is also the optimal solution of Py. As this solution is feasible and objective value is correct, we do not need linearize in this solution!
            solveLpFix = false;
            addLinearizationsOnMasterSol = false;
        }
        
        if( objValue < zu )
        {
            callbackSolver->getNumberOfIterations(iter);
            
            SEMAPH_updtSol.lock(nthreads_lazy);
            {
                tryUpdateBestSolution(thnumber, n, masterSol, objValue, iter, data.clockStart, data.timeStart, false);
            }
            SEMAPH_updtSol.unlock(nthreads_lazy);
            
            
            if( zl >= zu )
            { //we could, but, by now, we do not considere optimal tolerances here
                if(in_print_level > 2)
                    MRQ_PRINTMSG("A feasible solution better than lower bound was found. stopping");
                
                data.out_sol_lower_than_zl = true;
                
                retCode = MRQ_LAZY_MILP_BB_SOLUTION_LOWER_THAN_ZL;
                goto termination;
            }
        }
    }
    
    
    
    
    
    
    if( solveLpFix )
    {
        double nlpBeginTime;
        clock_t nlpBeginClock;
        
        const int nI = data.nI;
        const int *intVars = data.intVars;
        const unsigned int nintloop = in_max_number_of_fixed_relax_master_problem_solved_per_iteration;
        
        bool hasLpSol = false;
        bool pyCompletelySolved = false;
        bool pyOptimalFounded = false;
        unsigned int nSubIter = 0;
        int firstSetLinearizationIndex, nMasterConstrAfterFisrtLinearization, lastSetLinearizationIndex;
        
        
        double lpSolObjValue, bestOrLastSolPyObjValue;
        double nlpTime, nlpCPUTime;
        double *bestOrLastSolPy = data.auxVars3;
        MRQ_MasterMILPProb *masterLp = data.relaxMaster;
        MRQ_LPSolver *lp = masterLp->master;
        double *lpSolConstrValues = data.auxConstr2;// = lp->constr
        
        
        if( in_measure_nlp_time )
        {
            nlpBeginTime = MRQ_getTime();
            nlpBeginClock = clock();
        }
        
        zlpy = masterSol[n]; //value of variable alpha is the objective function
        zupy = MRQ_INFINITY;
        
        if( allPointsKeepers )
        {
            MRQ_PointsStoreKeeper &myPointsKeeper = allPointsKeepers[threadShift];
            
            r = myPointsKeeper.updateLinearizationOnPoints(*masterLp, linearizeObj, in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, incQuadsInMaster, in_constr_linearisation_strategy, constrEval, MRQ_OLS_ALL_POINTS, zu, NULL);
            MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
        }
        
        r = lp->getNumberOfConstraints(firstSetLinearizationIndex);
        MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_MILP_SOLVER_ERROR, termination);
        
        
        r = masterLp->addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, false, masterSol, incQuadsInMaster, in_constr_linearisation_strategy, constrEval, NULL, NULL, NULL, lpSolConstrValues, true);
        MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
            
        if(linearizeObj)
        {
            r = masterLp->addLinearizedObjFunction(false, masterSol, incQuadsInMaster, MRQ_OLS_ALL_POINTS, zu, NULL, NULL, NULL, &objValue);
            MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
        }
        
        r = lp->getNumberOfConstraints(nMasterConstrAfterFisrtLinearization);
        MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_MILP_SOLVER_ERROR, termination);
        
        r = MRQ_fixIntVarsOnSolByList(nI, intVars, masterSol, *lp);
        MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        
        
        for(unsigned int w = 0; w < nintloop; w++)
        {
            nSubIter++; //we use this to count correclty the number of iterations in the interior loop. Note, we cannot use w beacuse we can make w = nintloop to abort this loop, and so, after this for, w cannot be the correct number of iterations.
            
            //maybe another thread update zu...
            r = lp->setVariableBounds(n, -OPT_INFINITY, zu );
            #if MRQ_DEBUG_MODE
                MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_MILP_SOLVER_ERROR, termination);
            #endif
            
            r = lp->solve(false);
            
            //printf("%d R lp. cod: %d obj: %f master obj: %f orig master obj: %f zu: %f zupy: %f\n", w, lp->retCode, lp->objValue, masterSol[n], objValue, zu, zupy);
            
            if( lp->retCode == OPT_OPTIMAL_SOLUTION )
            {
                double *sol = lp->sol;
                bool feasSol;
                int r;
                
                hasLpSol = true;
                lpSolObjValue = NAN;
                
                //MRQ_copyArray(n, sol, lpSol);
                //double *lpSol = sol;
                
                
                if(lp->objValue > zlpy)
                    zlpy = lp->objValue;
                
                
                r = prob.isFeasibleToConstraints(thnumber, sol, true, constrEval, in_absolute_feasibility_tol, in_relative_convergence_tol, feasSol, lpSolConstrValues);
                MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_CALLBACK_FUNCTION_ERROR, termination);
                
                r = prob.objEval(thnumber, !prob.hasNlConstrs, sol, lpSolObjValue);
                MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_CALLBACK_FUNCTION_ERROR, termination);
                
                
                if(feasSol)
                {
                    #if MRQ_DEBUG_MODE
                        assert( zu <= zupy ); //we must have entered in if above
                    #endif
                    
                    if( lpSolObjValue < zupy )
                    {
                        zupy = lpSolObjValue;
                        bestOrLastSolPyObjValue = lpSolObjValue;
                        MRQ_copyArray(n, sol, bestOrLastSolPy);
                    }
                    
                    if(lpSolObjValue < zu)
                    {
                        /*r = callbackSolver->setSolution(n, sol);
                        if(r != 0)
                        {
                            if(in_print_level > 0)
                                MRQ_PRINTERRORNUMBER(r);
                            
                            retCode = r;
                            goto termination;
                        } */
                        
                        if( iter == (decltype(iter)) -1 )
                            callbackSolver->getNumberOfIterations(iter);
                        
                        SEMAPH_updtSol.lock(nthreads_lazy);
                        {
                            tryUpdateBestSolution(thnumber, n, sol, lpSolObjValue, iter, data.clockStart, data.timeStart, false);
                        }
                        SEMAPH_updtSol.unlock(nthreads_lazy);
                        
                        bestOrLastSolPy[n] = lpSolObjValue;
                        #if 0
                            r = callbackSolver->setSolution(n+1, bestOrLastSolPy);
                            if(r != 0)
                            {
                                MRQ_PRINTERRORNUMBER(r); //we do not abort here if we have some error.
                                std::cout << "Encontrei solução, mas não consegui passala ao solver de milp!\n";
                                MRQ_getchar();
                            }
                        #endif
                    }
                    
                }
                
                if( zupy >= MRQ_INFINITY ) //we do not have feasible sol to Py, sol, we store the lastSol
                {
                    bestOrLastSolPyObjValue = lpSolObjValue;
                    MRQ_copyArray(n, sol, bestOrLastSolPy);
                }
                
                //printf(" obj: %f", objValue);
                
                
                if( zupy - zlpy <= MRQ_max(in_absolute_convergence_tol, MRQ_abs(zupy*in_relative_convergence_tol) ) )
                { 
                    //so, we solve the continuous nlp relaxation problem fixing y. we can stop this solving. Even if w == 0, in this case we have lower bound equal to uper bound and so, we solve the original MINLP problem
                    w = nintloop;
                    pyCompletelySolved = true;
                    pyOptimalFounded = true;
                    //std::cout << "\t\tResolvi Py!!!";
                }
                
                
                
                if( !in_delete_intermediate_linearizations_in_each_iteration )
                {
                    //adding the cuts to master problem
                    r = addLazyConstraintsLinearizationOnSolution( thnumber, callbackSolver, prob, gradEval, incQuadsInMaster, linearizeObj, constrEval, sol, &lpSolObjValue, lpSolConstrValues, lpSolConstrValues, indices, plc, puc, auxVars, auxConstrEval2);
                    MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
                }
                
                
                //adding cuts to this relax master problem
                
                r = lp->getNumberOfConstraints(lastSetLinearizationIndex);
                MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_MILP_SOLVER_ERROR, termination);
                
                r = masterLp->addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, false, sol, incQuadsInMaster, in_constr_linearisation_strategy, constrEval, NULL, NULL, NULL, lpSolConstrValues, true);
                MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
                
                if(linearizeObj)
                {
                    r = masterLp->addLinearizedObjFunction(false, sol, incQuadsInMaster, MRQ_OLS_ALL_POINTS, zu, NULL, NULL, NULL, &lpSolObjValue);
                    MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
                }
                
            }
            else if( lp->retCode == OPT_INFEASIBLE_PROBLEM )
            {
                //std::cout << "\t\tPy inviavel!!!";
                pyCompletelySolved = true;
                break;
            }
            else
            { //we have some error
                break;
            }
            
        }
        
        
        if(pyOptimalFounded)
            addLinearizationsOnMasterSol = false; //we solve Py completely. So, we do not need more the cuts on the master solution
        
        
        if(!pyCompletelySolved)
        {
            //std::cout << "\t\tPy nao foi completamente resolvido!!!";
            out_all_subproblems_py_completely_solved = false; //we do not need put it in a exclusive zone.
        }
        
        
        if( hasLpSol ) 
        {
            //deleting intermediate linearizations. Note are deleting intermediate linearizations for lp solver, even if in_delete_intermediate_linearizations_in_each_iteration is false. We choose this because the number of linearization can be very large, and in this case, they are not added like lazy constraints.
            
            
            //Here, we olny keep the original master solution and the last lp solution in the linearization point (except if we found the optimal for py. In this case, we discard master solution also)
            int startIndicesToDelete, endIndicesToDelete;
            double lpSolObjValue = bestOrLastSolPyObjValue;
            double *lpSol = bestOrLastSolPy;
            
            if( pyOptimalFounded )
            {
                //if we found the optimal solution of Py, we can discard all previous linearizations from this iteration, inclusive the original ECP linearization (linearization from integer master problem)
                startIndicesToDelete = firstSetLinearizationIndex;
                endIndicesToDelete = lastSetLinearizationIndex-1;
            }
            else
            {
                startIndicesToDelete = nMasterConstrAfterFisrtLinearization;
                endIndicesToDelete = lastSetLinearizationIndex-1;
            }
            
            if( endIndicesToDelete > startIndicesToDelete ) 
            {
                r = lp->removeConstraintsByRange(startIndicesToDelete, endIndicesToDelete);
                MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_MILP_SOLVER_ERROR, termination);
            }
            
            
            if( in_delete_intermediate_linearizations_in_each_iteration ) //so we must add lazy constraints on this solution
            {
                
                r = prob.constraintsEval(thnumber, true, constrEval, lpSol, lpSolConstrValues);
                MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_CALLBACK_FUNCTION_ERROR, termination);
                
                //adding linearizations int last lp solution in the master problems and allPointsKeepers
                r = addLazyConstraintsLinearizationOnSolution( thnumber, callbackSolver, prob, gradEval, incQuadsInMaster, linearizeObj, constrEval, lpSol, &lpSolObjValue, lpSolConstrValues, lpSolConstrValues, indices, plc, puc, auxVars, auxConstrEval2);
                MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
            }
            
            //adding cuts to other threads storing...
            for(unsigned int i = 0; i < nthreads_lazy; i++)
            {
                if(i == threadShift)
                    continue; //that is the current thread. So, we already add the cuts
                    
                r = allPointsKeepers[i].insertPoint(n, lpSol);
                MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
            }
        }
        
        
        if( in_measure_nlp_time )
        {
            nlpTime = MRQ_getTime() - nlpBeginTime;
            nlpCPUTime = MRQ_calcCPUTtime( nlpBeginClock, clock() );
        }
        
        SEMAPH_updtOut.lock(nthreads_lazy);
        {
            out_number_of_master_relaxation_solved += nSubIter;
            if( in_measure_nlp_time )
            {
                out_clock_time_of_nlp_solving += nlpTime;
                out_cpu_time_of_nlp_solving += nlpCPUTime;
            }
            if(pyCompletelySolved)
                out_number_of_subproblems_py_completely_solved++;
        }
        SEMAPH_updtOut.unlock(nthreads_lazy);
    }	
    
    
    
    
    if( addLinearizationsOnMasterSol )
    {
        r = addLazyConstraintsLinearizationOnSolution( thnumber, callbackSolver, prob, gradEval, incQuadsInMaster, linearizeObj, constrEval, masterSol, &objValue, constrValues, constrValues, indices, plc, puc, auxVars, auxConstrEval2);
        MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        
        //adding cuts to other threads storing...
        if( allPointsKeepers )
        {
            for(unsigned int i = 0; i < nthreads_lazy; i++)
            {
                if(i == threadShift)
                    continue; //that is the current thread. So, we already add the cuts
                
                r = allPointsKeepers[i].insertPoint(n, masterSol);
                MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
            }
        }
    }
    
    
    retCode = 0;
    
termination:
    
    //MRQ_getchar();
    return retCode;
}



int MRQ_LPBBExtCutPlan::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    const double timeStart = MRQ_getTime();
    const clock_t clockStart = clock();
    
    const bool binProblem = prob.isBinaryProblem();
    const int n = prob.n;
    const int m = prob.m;
    const int nI = prob.getNumberOfIntegerVars();
    const int pType = prob.getProblemType();
    const bool preproc = in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_obj_function;
    
    const bool setQuadsInMaster = in_set_quadratics_in_master_problem;
    const bool linearizeObj = prob.hasObjNLTerm() || (prob.Q.getNumberOfElements() > 0 && !setQuadsInMaster);
    
    bool updtConstrBounds;
    int r;
    
    bool *constrEval = NULL, *auxConstrEval; 
    int *intVars = NULL;
    int *indices = NULL;
    double *constrValues = NULL;
    
    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux;
    double *plc = NULL, *puc;
    
    MRQ_Preprocessor preprocessor(&prob);
    MRQ_LAAPointsStoring laps(n);
    MRQ_GradientsEvaluation gradEval;
    
    OPT_LPSolver *nlp = NULL;
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
        MRQ_PRINTMSG("Starting LP Branch-And-Bound based on Extended Cutting Plane\n\n");
    }
    
    
    if(in_print_level > 1)
        printSubSolvers(true, true, false);
    
    MRQ_malloc(intVars, nI);
    MRQ_malloc(indices, n+1);
    MRQ_malloc(constrEval, 2*m);
    MRQ_malloc(constrValues, m);
    callbackData = new (std::nothrow) MRQ_MILPCallbackData[nthreads_lazy];
    MRQ_IFMEMERRORGOTOLABEL( !intVars || !indices || !constrEval || !constrValues || !callbackData, out_return_code, termination);
    
    
    auxConstrEval = &constrEval[m];
    
    prob.getIntegerIndices(intVars);
    
    {
        const bool *nlConstr = prob.nlConstr;
        const minlpproblem::MIP_SparseMatrix *QC = prob.QC;
        
        #pragma ivdep
        #pragma GCC ivdep
        for(int i = 0; i < m; i++)
            constrEval[i] = nlConstr[i] || (QC[i].getNumberOfElements() > 0);
    }
    
    if( in_max_number_of_fixed_relax_master_problem_solved_per_iteration > 0 )
    {
        pointsKeeper = new (std::nothrow) MRQ_PointsStoreKeeper[nthreads_lazy];
        MRQ_IFMEMERRORGOTOLABEL( !pointsKeeper, out_return_code, termination);
        
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < nthreads_lazy; i++)
            pointsKeeper[i].nthreads = nthreads_lazy;
    }
    
    
    for(unsigned int i = 0; i < nthreads_lazy; i++)
    {
        const bool useMilpRelax = in_max_number_of_fixed_relax_master_problem_solved_per_iteration > 0;
        
        int r = callbackData[i].allocateBase(thnumber + i, nthreads_lazy, &prob, in_milp_solver, useMilpRelax, NULL, &SEMAPH_addLazy);
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MEMORY_ERROR, termination);
        
        callbackData[i].binProblem = binProblem;
        callbackData[i].linearizeObj = linearizeObj;
        callbackData[i].setQuadsInMaster = setQuadsInMaster;
        
        callbackData[i].timeStart = timeStart;
        callbackData[i].clockStart = clockStart;
        
        callbackData[i].nI = nI;
        callbackData[i].intVars = intVars;
        callbackData[i].indices = indices;
        callbackData[i].constrEval = constrEval;
        
        
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
        
        //printf("pointsKeeper[%u].nthreads: %u", i, pointsKeeper[i].nthreads );
        
        if(useMilpRelax)
        {
            /**** Setting master's LP relaxation *****/
            MRQ_MasterMILPProb*  masterMilp;
            
            masterMilp = callbackData[i].relaxMaster;
            
            r = setMasterProblemBase(masterMilp, thnumber, prob, in_milp_solver, true, setQuadsInMaster, false, lx, ux, 1, milpSolverParams, constrEval, NULL, 1); //note, we set a milp problems actually, but there is no problem about that because when we solve this problem, integer variable will be fixed
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
            
        }
    }
    
    
    /***** Setting master problem *****/
    
    r = setMasterProblemBase(&masterMilp, thnumber, prob, in_milp_solver, true, setQuadsInMaster, true, lx, ux, 1, milpSolverParams, auxConstrEval, &laps, nthreads);
    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    
    master = masterMilp.master;
    
    r = master->setRelativeOptimalityTol( in_relative_convergence_tol );
    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    
    
    if(nPoints == 0)
    {
        //we have to improvise the first linearization point
        
        double objValue;
        double * const initSol = master->sol; //that is not so elegant, but we use master.sol because, otherwise, we would need allocate a one more array....
        
        for(int i = 0; i < n; i++)
            initSol[i] = MRQ_min(ux[i], MRQ_max(0.0, lx[i]) ); 
        
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
        if( in_max_number_of_fixed_relax_master_problem_solved_per_iteration > 0 )
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
        
    }
    
    
    #pragma ivdep
    #pragma GCC ivdep
    for(int i = 0; i <= n; i++) //considering auxiliary variables also
        indices[i] = i;
    
    if( in_measure_nlp_time )
    {
        out_cpu_time_of_nlp_solving = 0;
        out_clock_time_of_nlp_solving = 0;
    }
    
    /***** setting up solver to adopt lazy constraints, user callbacks and pseudo prunning *****/
    
    if( prob.getProblemType() != minlpproblem::MIP_PT_MILP )
    {
        r = setLazyConstraintCallbackOnMilpSolver(master, callbackData);
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    }
    
    if( in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
    {
        r = allocateMemoryForPseudoPruning(nthreads, n, nI);
        MRQ_IFERRORGOTOLABEL(r, out_return_code, (MRQ_RETURN_CODE) r, termination);
    }
    
    
    if(in_user_callbacks || in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING)
    {
        
        if(in_call_before_solve_callback_in_milp_bb || (in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING && in_try_pseudo_pruning_before_solving_relaxations) )
        {
            
            r = setBeforeSolveCallbackOnMilpSolver(master, callbackData, in_call_before_solve_callback_in_milp_bb);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
        
        if(in_call_branching_callback_in_milp_bb || in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING)
        {
            r = setBranchingCallbackOnMilpSolver(master, callbackData, in_call_branching_callback_in_milp_bb);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
        
        
        if( in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING && in_try_pseudo_pruning_before_solving_relaxations )
        {
            r = setDeleteNodeCallbackOnMilpSolver(master, callbackData, false);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
        
        
    }
    
    
    if( in_max_number_of_fixed_relax_master_problem_solved_per_iteration > 0 )
        out_all_subproblems_py_completely_solved = true; //that is weird, but we set this flag as true before the procedure
    
    /*master->generateModelFile("bbecp.lp");
    std::cout << "Gerei bbecp.lp\n";
    MRQ_getchar(); */
    
    r = master->solve();
    
    //std::cout << "master solving ret code: " << master->retCode << " obj: " << master->objValue << " orig ret code: " << master->origSolverRetCode << "\n";
    
    {
        double myzl = master->getDualObjValue();
        if( myzl > zl )
            zl = myzl;
    }
    
    if( master->feasSol )
    {
        if( master->objValue < MRQ_zuWithTol(out_best_obj, in_absolute_convergence_tol,  in_relative_convergence_tol )   ||  out_best_obj >= MRQ_INFINITY  )
        {
            MRQ_PRINTMSG("Warning: milp solver found a solution slightly infeasible\n");
            
            MRQ_copyArray(n, (const double*) master->sol, out_best_sol );
            
            int r = prob.objEval( thnumber, true, out_best_sol, out_best_obj );
            MRQ_IFCALLBACKERRORGOTOLABEL(r, out_return_code, termination);
        }
        
        #if MRQ_DEBUG_MODE
        #endif
    }
    
    /*{
        double *sol = master->sol;
        bool feas;
        double obj;
        
        double *lc = prob.lc, *uc = prob.uc;
        
        prob.isFeasibleToConstraints(0, sol, true, constrEval, in_absolute_feasibility_tol, in_relative_feasibility_tol, feas, constrValues);
        
        prob.objEval(0, true, sol, obj);
        
        for(int i = 0; i < n; i++)
        {
            printf("x[%d]: %0.20f\n", i, sol[i]);
        }
        
        std::cout << "feas: " << feas << " obj: " << obj <<  " absTol: " << in_absolute_feasibility_tol << " relTol: " << in_relative_feasibility_tol << " \n";
        
        for(int i = 0; i < m; i++)
        {
            if( lc[i] > -MIP_INFINITY )
                printf("%0.10f <= ", lc[i]);
            
            printf("%0.10f", constrValues[i]);
            
            if( uc[i] < MIP_INFINITY )
                printf(" <= %0.10f", uc[i]);
            
            printf(" \tcEval: %d\n", (int) constrEval[i]);
        } 
    } */
    
    
    
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
        //if the problem is nonconvex, we can have found a feasible solution and, even so, solver declares infeasibility
        if( out_best_obj < MRQ_INFINITY )
            out_return_code = MRQ_UNDEFINED_ERROR;
        else
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
    
    
    if( in_refine_final_solution_using_nlp  && out_best_obj < MRQ_INFINITY )
    {
        int r;
        double NLPCpuTime, NLPClockTime;
        
        //we try to use the milp solver to refine the solution if possible
        if( pType == MIP_PT_MILP )
        {
            nlp = OPT_newLPSolver(in_milp_solver);
        }
        else
        {
            const auto milpSolverType = master->getSolverType();
            
            if( pType == MIP_PT_MIQP && OPT_isQPSolverType(milpSolverType) )
                nlp = OPT_newQPSolver(in_milp_solver);
            else if( pType == MIP_PT_MIQCP && OPT_isQCPSolverType(milpSolverType) )  
                nlp = OPT_newQCPSolver(in_milp_solver);
            else if ( pType == MIP_PT_MINLP && OPT_isNLPSolverType(milpSolverType) ) 
                nlp = OPT_newNLPSolver(in_milp_solver);
            else
                nlp = OPT_newNLPSolver(in_nlp_solver);
        }
        
        
        MRQ_IFERRORGOTOLABEL(!nlp, out_return_code, out_return_code, termination); //we do not change the return code
        
        
        //here, we do not preprocess. Since we are already fixing integer variable, we let it to nlp solver preprocesor, although some solvers do not have one.
        r = MRQ_setNLPRelaxProb( prob, lx, ux, plc, puc, nlp, true, true, true, false, thnumber, in_set_special_nlp_solver_params, nlpSolverParams, in_number_of_threads, in_max_cpu_time, in_max_time, 0, 0 );
        
        MRQ_IFERRORGOTOLABEL(r, out_return_code, out_return_code, termination); //we do not change
        
        if( run_by_inside )
        {
            if( !std::isinf(insideSolverMaxTime) )
                nlp->setMaxTime(insideSolverMaxTime );
        }
        
        MRQ_fixIntVarsOnSolByList(nI, intVars, out_best_sol, *nlp);
        
        
        if( nlp->getSolverType() == OPT_NLP )
            ( (OPT_NLPSolver*) nlp)->setInitialSolution(out_best_sol, NULL, NULL);
        
        
        
        
        if( in_max_cpu_time < INFINITY || in_max_time < INFINITY )
        {
            const double rcputime = in_max_cpu_time - ( ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC);
            const double rtime = in_max_time - (MRQ_getTime() - timeStart);
            int r = 0;
            
            if(rcputime <= 0.0 || rtime <= 0.0)
                goto termination;
            
            if(in_max_cpu_time < INFINITY)
                r = nlp->setMaxCPUTime(rcputime);
            
            if(in_max_time < INFINITY)
                r += nlp->setMaxTime(rtime);
            
            if(r != 0)
            {
                if(in_print_level > 0)
                    MRQ_PRINTERRORMSG("Error to set NLP maximum time!");
            }
        }
        
        
        
        nlp->solveAndGetTime( in_measure_nlp_time ? &NLPCpuTime : NULL, in_measure_nlp_time ? &NLPClockTime : NULL, false);
        
        if( in_measure_nlp_time )
        {
            out_cpu_time_of_nlp_solving += NLPCpuTime;
            out_clock_time_of_nlp_solving += NLPClockTime;
        }
        
        out_number_of_nlp_probs_solved++;
        
        /*if(nlp->isMyNLPClass())
        {
            ( (optsolvers::OPT_MyNLPSolver *) nlp)->generateModelFile("modelo_nlp.txt");
        }*/
        
        
        if( nlp->retCode == OPT_OPTIMAL_SOLUTION )
        {
            if(in_print_level > 4)
                std::cout << MRQ_PREPRINT "Refinement problem solved. Old Objective: " << out_best_obj << " New objective: " << nlp->objValue << "\n";
            out_best_obj = nlp->objValue;
            MRQ_copyArray(n, nlp->sol, out_best_sol);
        }
        else
        {
            MRQ_PRINTERRORMSGP("Error to solve the nlp refinement problem: ", nlp->retCode << " solver error code: " << nlp->origSolverRetCode );
        }
        
    }
    
    
    
    
termination:
    
    if(plc)			free(plc);
    
    if(intVars)			free(intVars);
    if(indices)			free(indices);
    if(constrEval) 		free(constrEval);
    if(constrValues)	free(constrValues);
    
    if(callbackData)	delete[] callbackData;
    if(pointsKeeper)	delete[] pointsKeeper;
    
    if(nlp)				delete nlp;
    
    
    algorithmFinalization(nthreads, prob, lx, ux);
    
    out_number_of_threads = nthreads_lazy;
    out_number_of_milp_solver_iters = out_number_of_iterations;
    out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
    out_clock_time = MRQ_getTime() - timeStart;
    
    if(in_print_level > 1)
        std::cout << MRQ_PREPRINT "cpu time: " << out_cpu_time << "\n";
    
    
    return out_return_code;
}

