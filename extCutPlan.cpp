#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#include <iostream>
#include <new>

#include "BBL_tools.hpp"

#include "MRQ_solvers.hpp"

#include "MRQ_dataStructures.hpp"
#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"

#if MRQ_DEBUG_MODE
    #define MRQ_CHECK_INT_SOLS_ARE_REPEATING 1
#endif

//using namespace std;

using namespace optsolvers;
using namespace muriqui;

MRQ_ExtCutPlan::MRQ_ExtCutPlan():MRQ_LinearApproxAlgorithm()
{
    resetParameters();
    resetOutput();
    
    out_algorithm = MRQ_ECP_ALG;
}


MRQ_ExtCutPlan::~MRQ_ExtCutPlan()
{
}


void MRQ_ExtCutPlan::printParameters(std::ostream& out) const
{
    MRQ_LinearApproxAlgorithm::printParameters(out);
    out << "\n"
    
    //MRQ_STRFFATT(in_fix_int_vars_to_try_improve_cut) << "\n"
    MRQ_STRFFATT(in_delete_intermediate_linearizations_in_each_iteration) << "\n"
    
    MRQ_STRFFATT(in_refine_final_solution_using_nlp) << "\n"
    MRQ_STRFFATT(in_max_number_of_fixed_relax_master_problem_solved_per_iteration) << "\n"
    MRQ_STRFFATT(in_subiter_printing_frequency) << "\n"
    ;
}


void MRQ_ExtCutPlan::resetOutput()
{
    MRQ_LinearApproxAlgorithm::resetOutput();
    
    out_all_subproblems_py_completely_solved = false;
    out_number_of_subproblems_py_completely_solved = 0;
    
    out_number_of_integer_solution_repetitions = 0;
    out_integer_solution_history.desallocate();
    
    out_number_of_master_relaxation_solved = 0;
}

void MRQ_ExtCutPlan::resetParameters()
{
    MRQ_LinearApproxAlgorithm::resetParameters();
    
    //in_fix_int_vars_to_try_improve_cut = false;
    
    in_delete_intermediate_linearizations_in_each_iteration = false;
    
    in_refine_final_solution_using_nlp = true;
    in_max_number_of_fixed_relax_master_problem_solved_per_iteration = 1;
    
    in_subiter_printing_frequency = 3;
    
    //in_max_iterations = 100000;
    //in_absolute_convergence_tol = 1.0e-3;
    //in_relative_convergence_tol = 1.0e-5;
}



int MRQ_ExtCutPlan::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_LinearApproxAlgorithm::setIntegerParameter(name, value);
    
    if(ret == 0)
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_delete_intermediate_linearizations_in_each_iteration), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_refine_final_solution_using_nlp), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_max_number_of_fixed_relax_master_problem_solved_per_iteration), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_subiter_printing_frequency), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    return ret;
}



int MRQ_ExtCutPlan::run(MRQ_MINLPProb& prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    const double timeStart = MRQ_getTime();
    const clock_t clockStart = clock();
    
    
    const bool nlObj = prob.hasNlObj;
    bool linearizeObj;
    const int n = prob.n;
    const int m = prob.m;
    
    const unsigned int nintloop = in_max_number_of_fixed_relax_master_problem_solved_per_iteration ;
    
    bool updtConstrBounds;
    bool setQuadsInMaster = in_set_quadratics_in_master_problem;
    bool feasible;
    int aux, nI;
    unsigned long int iter = 0;
    double objMasterSol;
    
    
    MRQ_NLPSolver *nlp = NULL;
    OPT_LPSolver *master, *masterLp = NULL;
    MRQ_MasterMILPProb masterMilp, masterMilpRelax;
    
    MRQ_GradientsEvaluation  gradEval;
    MRQ_LAAPointsStoring laps(n);
    MRQ_Preprocessor preprocessor(&prob);
    
    //arrays:
    bool *auxCEval = NULL;
    int *auxCols = NULL;
    int *intVars = NULL;
    double *auxVars = NULL, *auxConstr = NULL, *lpsol = NULL;
    double *masterSolConstr = NULL;
    double *plc = NULL, *puc;
    double *bestSolPy = NULL;
    
    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux; 
    
    
    #if MRQ_CHECK_INT_SOLS_ARE_REPEATING
        
        //double auxHistorySol[ prob.getNumberOfIntegerVars() ];
        
        double *auxHistorySol;
        MRQ_malloc(  auxHistorySol, prob.getNumberOfIntegerVars() );
        MRQ_IFMEMERRORGOTOLABEL( !auxHistorySol, out_return_code, termination );
    #endif
    
    nthreads = 1;
    nthreads_lazy = in_number_of_threads > 0 ? in_number_of_threads : branchAndBound::BBL_getNumCores() ;
    if( in_milp_solver == MRQ_GUROBI )
        nthreads_lazy = 1; //gurobi apply multithreading to solve the problem, but only thread 0 is called to add lazy constraints.
    
    { 
        auto r = algorithmInitialization(1, (in_preprocess_lin_constr || in_preprocess_obj_function || in_preprocess_quad_constrs), milpSolverParams, nlpSolverParams, prob, lx, ux, &preprocessor, &updtConstrBounds, &plc, &puc); //that algorithm is monothread...
        MRQ_IFERRORGOTOLABEL(r, out_return_code, r, termination);
    }
    
    
    //we do not need preprocessor more...
    MRQ_secFree(plc);
    preprocessor.deallocateMemory();
    
    if(in_print_level > 1)
        std::cout << "\nStarting Extended Cutting Plane algorithm\n\n";
    
    if(in_print_level > 3)
        printSubSolvers(true, false, false);
    
    
    aux = gradEval.initialize(thnumber, &prob);
    MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MEMORY_ERROR, termination);
    
    
    nI = prob.getNumberOfIntegerVars();
    
    MRQ_malloc(intVars, nI); //intVars = (int *) malloc( nI * sizeof(int) );
    MRQ_malloc(auxCEval, m); //auxCEval = (bool *) malloc(m * sizeof(bool) );
    MRQ_malloc(auxConstr, m); //auxConstr = (double *) malloc(m * sizeof(double) );
    MRQ_malloc(masterSolConstr, m);
    MRQ_malloc(auxCols, (n+1)); //auxCols = (int *) malloc( (n+1)*sizeof(int) );
    MRQ_malloc(auxVars, (n+1)); //auxVars = (double *) malloc( (n+1)*sizeof(double) );
    MRQ_IFMEMERRORGOTOLABEL( !intVars || !auxCEval || !auxConstr || !masterSolConstr || !auxCols || !auxVars, out_return_code, termination );
    
    if( nintloop > 0 )
    {
        MRQ_malloc(lpsol, n+1);
        MRQ_malloc(bestSolPy, n+1);
        MRQ_IFMEMERRORGOTOLABEL( !lpsol || !bestSolPy, out_return_code, termination );
    }
    
    
    prob.getIntegerIndices(intVars);
    
    
    for(int i = 0; i < m; i++)
        auxCEval[i] = prob.nlConstr[i] || prob.QC[i].getNumberOfElements() > 0;
    
    
    aux = masterMilp.setProblemBase( thnumber, prob, in_milp_solver, true, setQuadsInMaster, true, lx, ux, 1, milpSolverParams, in_number_of_threads );
    MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    
    master = masterMilp.master;
    
    if( nintloop > 0 )
    {
        aux = masterMilpRelax.setProblemBase( thnumber, prob, in_milp_solver, true, setQuadsInMaster, false, lx, ux, 1, milpSolverParams, in_number_of_threads );
        MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        
        masterLp = masterMilpRelax.master;
    }
    
    
    aux = master->getSolverType();
    if( aux == OPT_LP || (aux == OPT_QP && prob.hasQuadMatrixInSomeConstraint()  ) )
    {
        if( in_set_quadratics_in_master_problem && in_print_level > 0 )
            std::cerr << MRQ_PREPRINT "Warning: parameter in_set_quadratics_in_master_problem is set to true, but milp solver does not support quadratic constraints currently. Changing the parameter to false.\n";
        
        setQuadsInMaster = false; //integer solver does not support the quadratics
    }
    
    linearizeObj = nlObj || (prob.Q.getNumberOfElements() > 0 && !setQuadsInMaster);
    
    
    if( nPoints > 0 )
    {
        aux = masterMilp.addConstraintLinearisationPointsToMILP( in_eps_to_active_constr_to_linearisation,  &out_number_of_constr_linears_saved, nPoints, points, setQuadsInMaster, in_constr_linearisation_strategy, auxCEval);
        MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        
        aux = masterMilp.addObjLinearisationPointsToMILP( nPoints, points, zu, setQuadsInMaster, in_obj_linearisation_strategy, NULL, NULL, &laps );
        MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    }
    
    
    
    if( nPoints == 0 || aux != 0 )
    {
        //so, we need get a first linearization point
        double * const initSol = master->sol; //that is not so elegant, but we use master.sol because, otherwise, we would need allocate a one more array....
        
        for(int i = 0; i < n; i++)
            initSol[i] = MRQ_min(ux[i], MRQ_max(0.0, lx[i]) ); 
        
        
        aux = masterMilp.addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, true, initSol, setQuadsInMaster, in_constr_linearisation_strategy, auxCEval );
        MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        
        if( nlObj || (prob.Q.getNumberOfElements() > 0 && !setQuadsInMaster) )
        {
            int r, mmaster;
            
            aux = masterMilp.addLinearizedObjFunction( !prob.hasNlConstrs, initSol, setQuadsInMaster, auxCols, auxVars, NULL );
            MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
            
            if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS )
            {
                //auxVars has the coeficients of linearization to objective function and the RHS also
                r = laps.updateObjLinearizationsByNonObjCutPointsByNewPoint( *master, zu, auxVars );
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                r = laps.addPoint(n, initSol);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MEMORY_ERROR, termination);
                
                r = master->getNumberOfConstraints( mmaster );
                #if MRQ_DEBUG_MODE
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                #endif
                
                laps.indMaster[ laps.npoints-1 ] = mmaster;
            }
            
        }
    }
    
    
    if( run_by_inside )
    {
        if( !std::isinf(insideSolverMaxTime) )
        {
            aux = master->setMaxTime(insideSolverMaxTime);
            MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
    }
    
    
    if( zu < MRQ_INFINITY )
    {
        aux = master->setVariableBounds(n, -OPT_INFINITY, zu );
        MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    }
    
    
    if(in_print_level > 1)
        printf(MRQ_PREPRINT "%5s  %-14s  %-14s  %-14s\n", "iter", "    lbound", "    ubound", "     gap");
    
    if( nintloop > 0 )
        out_all_subproblems_py_completely_solved = true; //we set this variable here to be true: we consider all subproblems solved until now were complete solved. The strange thing is if we have some error, this variable can be set as true, but I do not know how to solve it. We need this information even if we do not reach optimal solution of MINLP problem, e.g., when we reach maximum number of iterations or maximum time...  :/
    
    
    //starting outer approximation loop
    while(true)
    {
        bool pyOptimalFounded = false; //true if we found a optimal soultion for Py
        bool pyCompletelySolved = false; //true if Py was completely solved: we cand find a optimal solution, or detect infeasibility
        int r, firstSetLinearizationIndex; //index to put the first set of linearization
        int nNewLpSol = 0;
        double zlpy, zupy = MRQ_INFINITY; //upper bound to solve P_y
        double *mastersol;
        
        
        iter++;
        
        /*{
            char name[100];
            
            sprintf(name, "ecp_%s_iter_%lu.lp", master->getSolverName().c_str(), iter);
            master->generateModelFile( name );
        }*/
        
        if( in_set_obj_lower_bound_on_master_problem )
            master->setVariableBounds(n, zl, zu);
        
        
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
                //we let algorithm follow even so...
            }
        }
        
        
        aux = master->solve(false);
        
        /*printf("MILP status: %d objF: %f#####################################################################\n", aux, master->objValue);
        //for(int i = 0; i < n; i++)
            //printf("s%d: %0.4f   ", i, master->sol[i]);
        //printf("\n");
        //MRQ_getchar();*/

        
        #if MRQ_CHECK_INT_SOLS_ARE_REPEATING
        if( aux == OPT_OPTIMAL_SOLUTION )
        {
            bool repeatedSol = false;
            double *msol = master->sol;
            
            for(int i = 0; i < nI; i++)
                auxHistorySol[i] = round( msol[ intVars[i] ] );
            
            //I think it is better run in reverse way because I think it is more probable repeted the last solutions
            for(unsigned int k2 = out_integer_solution_history.getnsols(); k2 > 0 ; k2-- )
            {
                unsigned int k = k2-1; //we put here because k and k2 are unsigned
                
                bool equalSols = true;
                MRQ_HistorySolution *hssol = out_integer_solution_history.getHistSolPointer(k);
                
                double *hsol = hssol->sol;
                
                for(int i = 0; i < nI; i++)
                {
                    if( hsol[i] != auxHistorySol[i] )
                    {
                        equalSols = false;
                        break;
                    }
                }
                
                if(equalSols)
                { //we use time inside history solution to count how many times solution apeared
                    hssol->time += 1;
                    repeatedSol = true;
                    out_number_of_integer_solution_repetitions++; 
                    //std::cout << "Repeti solucao da iteracao " << hssol->iter << ". Obj antigo: " << hssol->objF << " Obj novo: " << master->objValue << "\n";
                    
                    /*for(int i = 0; i < nI; i++)
                    {
                        std::cout << "i: " << i << " a: " << hsol[i] << " n: " << auxHistorySol[i] << "    ";
                        assert(hsol[i] == auxHistorySol[i]);
                    }*/
                    
                    //MRQ_getchar();
                    break;
                }
            }
            
            if( !repeatedSol )
            {
                ////we use time inside history solution to count how many times solution apeared. So, we pass 1.0 like time value
                int r = out_integer_solution_history.addSolution(nI, iter, 1.0, NAN, auxHistorySol, master->objValue);
                
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MEMORY_ERROR, termination);
            }
            
            #if 0
            {
                //here, we solve continuous relaxtions of the problem just to see if this integer solution is feasible. This is not part of algorithm. is just to evaluate the behaviour of algorithm
                MRQ_ContinuousRelax cr;
                double tlx[n], tux[n];
                
                MRQ_copyArray(n, lx, tlx);
                MRQ_copyArray(n, ux, tux);
                
                for(int i = 0; i < nI; i++)
                {
                    auto ind = intVars[i];
                    tlx[ind] = auxHistorySol[i];
                    tux[ind] = auxHistorySol[i];
                }
                
                cr.in_print_level = 0;
                
                MRQ_insideRun( &cr, prob, NULL, nlpSolverParams, thnumber, INFINITY, tlx, tux);
                
                std::cout << "iter: " << iter << " Resolvi Py. ret code: " << cr.out_return_code << " obj: " << cr.out_best_obj << "\n";
            }
            #endif
        }
        #endif
                
        {  //couting number of iterations
            long unsigned int masterIters;
            int r = master->getNumberOfIterations(masterIters);
            if(r == 0)
                out_number_of_milp_solver_iters += masterIters;
            else if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
        }
        
        
        if( aux != OPT_OPTIMAL_SOLUTION )
        {
            if( aux == OPT_INFEASIBLE_PROBLEM )
                out_return_code = out_best_obj < MRQ_INFINITY ? MRQ_OPTIMAL_SOLUTION : MRQ_INFEASIBLE_PROBLEM ;
            else if ( aux == OPT_MAX_TIME )
                out_return_code = MRQ_MAX_TIME_STOP;
            else if( aux == OPT_UNBOUNDED_PROBLEM)
                out_return_code = MRQ_UNBOUNDED_PROBLEM;
            else
                out_return_code = MRQ_MILP_SOLVER_ERROR;
            
            break;
        }
        
        
        if( master->objValue > zl )
            zl = master->objValue;
        
        
        mastersol = master->sol;
        
        aux = prob.isFeasibleToConstraints(thnumber, mastersol, true, auxCEval, in_absolute_feasibility_tol, in_relative_feasibility_tol, feasible, masterSolConstr);
        MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_CALLBACK_FUNCTION_ERROR, termination);
        
        aux = prob.objEval(thnumber, !prob.hasNlConstrs, mastersol, objMasterSol);
        MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_CALLBACK_FUNCTION_ERROR, termination);
                
        if(feasible)
        {
            if( objMasterSol < zupy )
            {
                zupy = objMasterSol;
                if(bestSolPy)
                    MRQ_copyArray(n+1, mastersol, bestSolPy);
            }
            
            const bool updt = tryUpdateBestSolution( thnumber, n, mastersol, objMasterSol, iter, clockStart, timeStart, in_store_history_solutions );
            
            if( updt )
            {
                aux = master->setVariableBounds(n, -OPT_INFINITY, zu );
                #if MRQ_DEBUG_MODE
                    MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                #endif
                
                if( masterLp )
                {
                    aux = masterLp->setVariableBounds(n, -OPT_INFINITY, zu );
                    #if MRQ_DEBUG_MODE
                        MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                    #endif
                }
            }
        
        }
        
        
        //we already have a feasible solution...
        if(zu - zl <= in_absolute_convergence_tol || zu - zl <= MRQ_abs(zu)*in_relative_convergence_tol)
        {
            if(out_best_obj < MRQ_INFINITY)
                out_return_code = MRQ_OPTIMAL_SOLUTION;
            else
                out_return_code = MRQ_INFEASIBLE_PROBLEM;
            break;
        }
        
        
        if( nintloop > 0 )
        {
            unsigned int w;
            double objLpSol;
            double nlpBeginTime;
            clock_t nlpBeginClock;
            double *pLinSol = NULL, *pObjLinSol;
            
            
            if( in_measure_nlp_time )
            {
                nlpBeginTime = MRQ_getTime();
                nlpBeginClock = clock();
            }
            
            
            zlpy = master->objValue;
            
            r = master->getNumberOfConstraints( firstSetLinearizationIndex);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
            
            { //udating the relaxation of master contraints: we will copy linearization from master problem instead of recalulated them
                int m_masterLp;
                
                r = masterLp->getNumberOfConstraints(m_masterLp);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                #if MRQ_DEBUG_MODE
                    assert( firstSetLinearizationIndex >= m_masterLp);
                #endif
                
                if( firstSetLinearizationIndex - m_masterLp > 0 )
                {
                    r = masterLp->addConstraints( firstSetLinearizationIndex - m_masterLp );
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                    
                    r = OPT_copyConstraintLinearParts(m_masterLp, firstSetLinearizationIndex-1, m_masterLp, *master, *masterLp);
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                }
                
                //linearizing in the master solution
                aux = masterMilpRelax.addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, false, mastersol,  setQuadsInMaster, in_constr_linearisation_strategy, auxCEval, NULL, NULL, NULL, masterSolConstr, true );
                MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                if( linearizeObj ) //nonlinear objective
                {
                    aux = masterMilpRelax.addLinearizedObjFunction( false, mastersol, setQuadsInMaster, auxCols, auxVars, &objMasterSol );
                    MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                }
                
            }
            
            /*master->generateModelFile("modelo1.lp");
            masterLp->generateModelFile("modelo2.lp");
            
            std::cout << "Gerei arquivos de saida!\n";
            MRQ_getchar();*/
            
            r = MRQ_fixIntVarsOnSolByList(nI, intVars, mastersol, *masterLp);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
            
            //objLpSol = objMasterSol;
            //MRQ_copyArray(n+1, mastersol, lpsol);
            //MRQ_copyArray(m, masterSolConstr, auxConstr);
            
            
            
            for(w = 0; w < nintloop; w++)
            {
                masterLp->solve( false );
                
                //std::cout << "\t" << w << " resolvi master lp. retcode: " << masterLp->retCode << " obj: " << masterLp->objValue << "\n";
                
                if( masterLp->retCode == OPT_INFEASIBLE_PROBLEM )
                {
                    //std::cout << "Py inviavel!\n";
                    pyCompletelySolved = true;
                    break;
                }
                else if( masterLp->retCode != OPT_OPTIMAL_SOLUTION )
                {
                    break;
                }
                
                nNewLpSol++;
                
                MRQ_copyArray(n+1, masterLp->sol, lpsol);
                
                if( masterLp->objValue > zlpy )
                    zlpy = masterLp->objValue;
                
                #if MRQ_DEBUG_MODE
                    assert(zlpy <= zu + MRQ_abs(zu* 1e-4) + 1e-4); //we set zu as upper bound of alpha variable. So, the lp relax problem cannot give us lower bound greater 
                #endif
                
                
                aux = prob.isFeasibleToConstraints(thnumber, lpsol, true, auxCEval, in_absolute_feasibility_tol, in_relative_feasibility_tol, feasible, auxConstr);
                MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_CALLBACK_FUNCTION_ERROR, termination);
                
                
                aux = prob.objEval(thnumber, !prob.hasNlConstrs, lpsol, objLpSol);
                MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_CALLBACK_FUNCTION_ERROR, termination);
                
                
                if(feasible)
                {
                    if( objLpSol < zupy )
                    {
                        zupy = objLpSol;
                        MRQ_copyArray(n+1, lpsol, bestSolPy);
                    }
                    
                    const bool updt = tryUpdateBestSolution( thnumber, n, lpsol, objLpSol, iter, clockStart, timeStart, in_store_history_solutions );
                    
                    if( updt )
                    {
                        aux = master->setVariableBounds(n, -OPT_INFINITY, zu );
                        #if MRQ_DEBUG_MODE
                            MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                        #endif
                        
                        aux = masterLp->setVariableBounds(n, -OPT_INFINITY, zu );
                        #if MRQ_DEBUG_MODE
                            MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                        #endif
                    }
                    
                }
                
                
                aux = masterMilpRelax.addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, false, lpsol,  setQuadsInMaster, in_constr_linearisation_strategy, auxCEval, NULL, NULL, NULL, auxConstr, true );
                MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                
                if( linearizeObj  ) //nonlinear objective
                {
                    aux = masterMilpRelax.addLinearizedObjFunction( false, lpsol, setQuadsInMaster, auxCols, auxVars, &objLpSol );
                    MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                }
                
                
                //adding this point to  master milp.
                if( !in_delete_intermediate_linearizations_in_each_iteration )
                {
                    aux = masterMilp.addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, false, lpsol,  setQuadsInMaster, in_constr_linearisation_strategy, auxCEval, NULL, NULL, NULL, auxConstr, true );
                    MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                    
                    if( linearizeObj  ) //nonlinear objective
                    {
                        //by now, we always linearize the last point in the master problem
                        
                        aux = masterMilp.addLinearizedObjFunction( false, lpsol, setQuadsInMaster, auxCols, auxVars, &objLpSol );
                        MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                        
                        if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS )
                        {
                            int r, mmaster;
                            
                            //now, we check if we can discard the linearization on previous point...
                            
                            //auxVars has the gradient of obj linearization, i.e., the coeficients of obj linearizations constraint and the rhs at the end...
                            r = laps.updateObjLinearizationsByNonObjCutPointsByNewPoint( *master, zu, auxVars );
                            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                            
                            r = laps.addPoint(n, lpsol);
                            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MEMORY_ERROR, termination);
                            
                            r = master->getNumberOfConstraints( mmaster);
                            #if MRQ_DEBUG_MODE 
                                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                            #endif
                            
                            laps.indMaster[ laps.npoints-1 ] = mmaster-1;
                        }
                        
                    }
                }
                
                
                
                if( zupy - zlpy <= MRQ_max(in_absolute_convergence_tol, MRQ_abs(zupy*in_relative_convergence_tol) ) )
                { //we solved Py. We linearize milp amster on optimal solution
                    //std::cout << "Resolvi Py\n";
                    pyCompletelySolved = true;
                    pyOptimalFounded = true;
                    break;
                }
                
                
            }
            
            
            //deleting intermediate linearizations from Py. It is easier in this way. In the next iteration, we will update masterLp with constraints added to Master
            if( pyOptimalFounded || (in_delete_intermediate_linearizations_in_each_iteration && nNewLpSol > 1 ) )
            {
                int endm;
                aux = masterLp->getNumberOfConstraints(endm);
                MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                if( endm > firstSetLinearizationIndex )
                {
                    aux = masterLp->removeConstraintsByRange(firstSetLinearizationIndex, endm-1);
                    MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                }
            }
            
            
            if( in_measure_nlp_time )
            {
                out_clock_time_of_nlp_solving += MRQ_getTime() - nlpBeginTime;
                out_cpu_time_of_nlp_solving += MRQ_calcCPUTtime(nlpBeginClock, clock());
            }
            
            
            if( pyCompletelySolved )
                out_number_of_subproblems_py_completely_solved++;
            else
                out_all_subproblems_py_completely_solved = false; //when we detect optimal solution or other termination condition, we skip this lines, but I think there is no problem about that.
            
            
            
            if( pyOptimalFounded )
            {
                if( !in_delete_intermediate_linearizations_in_each_iteration )
                { //since we found an optimal solution, we delete intermediate linearizations even so.
                    int mmaster;
                    
                    r = master->getNumberOfConstraints(mmaster);
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                    
                    if( mmaster > firstSetLinearizationIndex )
                    {
                        r = master->removeConstraintsByRange(firstSetLinearizationIndex, mmaster-1);
                        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                    }
                }
                
                #if MRQ_DEBUG_MODE
                    assert(nNewLpSol > 0);
                #endif
                
                pLinSol = bestSolPy;
                pObjLinSol = &zupy;
                
                //pLinSol = lpsol;
                //pObjLinSol = &objLpSol;
            }
            else if( in_delete_intermediate_linearizations_in_each_iteration )
            {
                //we must linearize on las t sol and master milp solution also
                if( nNewLpSol > 0 )
                {
                    pLinSol = lpsol;
                    pObjLinSol = &objLpSol;
                }
            }
            
            
            if( pLinSol )
            {
                aux = masterMilp.addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, true, pLinSol,  setQuadsInMaster, in_constr_linearisation_strategy, auxCEval, NULL, NULL, NULL, auxConstr, false );
                MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                if( linearizeObj  ) //nonlinear objective
                {
                    aux = masterMilp.addLinearizedObjFunction( !prob.hasNlConstrs, pLinSol, setQuadsInMaster, auxCols, auxVars, pObjLinSol );
                    MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                    
                    if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS )
                    {
                        int r, mmaster;
                        
                        //now, we check if we can discard the linearization on previous point...
                        
                        //auxVars has the gradient of obj linearization, i.e., the coeficients of obj linearizations constraint and the rhs at the end...
                        r = laps.updateObjLinearizationsByNonObjCutPointsByNewPoint( *master, zu, auxVars );
                        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                        
                        r = laps.addPoint(n, pLinSol);
                        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MEMORY_ERROR, termination);
                        
                        r = master->getNumberOfConstraints( mmaster);
                        #if MRQ_DEBUG_MODE 
                            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                        #endif
                        
                        laps.indMaster[ laps.npoints-1 ] = mmaster-1;
                    }
                }
            }
            
            
        }
        
        
        
        if(in_print_level > 1)// && w % in_subiter_printing_frequency == 0)
            printf(MRQ_PREPRINT "%5lu.%-2u  %+-14e  %+-14e  %+-14e\n", iter, 0, zl, zu, zu-zl); //TODO: replace the 0 by w
        
        
        if( checkTerminationCriterions(thnumber, zl, zu, iter, timeStart, clockStart, out_return_code) )
        {
            break;
        }
        
        
        if( !pyOptimalFounded )
        { //so, we must linearize also in the master problem solution
            aux = masterMilp.addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, true, mastersol,  setQuadsInMaster, in_constr_linearisation_strategy, auxCEval, NULL, NULL, NULL, masterSolConstr, true );
            MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
            
            if( linearizeObj  ) //nonlinear objective
            {
                aux = masterMilp.addLinearizedObjFunction( !prob.hasNlConstrs, mastersol, setQuadsInMaster, auxCols, auxVars, &objMasterSol );
                MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS )
                {
                    int r, mmaster;
                    
                    //now, we check if we can discard the linearization on previous point...
                    
                    //auxVars has the gradient of obj linearization, i.e., the coeficients of obj linearizations constraint and the rhs at the end...
                    r = laps.updateObjLinearizationsByNonObjCutPointsByNewPoint( *master, zu, auxVars );
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                    
                    r = laps.addPoint(n, mastersol);
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MEMORY_ERROR, termination);
                    
                    r = master->getNumberOfConstraints( mmaster);
                    #if MRQ_DEBUG_MODE 
                        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                    #endif
                    
                    laps.indMaster[ laps.npoints-1 ] = mmaster-1;
                }
            }
            
            #if MRQ_DEBUG_MODE
            if( masterLp && ( !in_delete_intermediate_linearizations_in_each_iteration || nNewLpSol <= 1) )
            {
                int mlp, mmilp;
                
                masterLp->getNumberOfConstraints(mlp);
                master->getNumberOfConstraints(mmilp);
                
                if(mlp != mmilp)
                    printf("mlp: %d mmilp: %d\n", mlp, mmilp);
                
                assert(mlp == mmilp);
            }
            #endif
        }
        
        
        #if 0
        for(unsigned int w = 0; w < nintloop; w++)
        {
            //we just need evaluate nonlinear constraint to check feasibility. Linear constraints are already enforced...
            
            objCalculated = false;
            
            
            aux = prob.isFeasibleToConstraints(thnumber, mastersol, true, auxCEval, in_absolute_feasibility_tol, in_relative_feasibility_tol, feasible, auxConstr);
            
            if(aux != 0)
            {
                if( in_print_level > 0 )
                    MRQ_PRINTCALLBACKERRORNUMBER(aux);
                
                if(w > 0)
                    break;
                
                out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
                goto termination;
            }
            
            
            if(feasible)
            {
                const int r = prob.objEval(thnumber, !prob.hasNlConstrs, mastersol, objTemp);
                
                if( r == 0 )
                {
                    objCalculated = true;
                    
                    #if 0
                    if( objTemp < zu )
                    {
                        zu = out_best_obj = objTemp;
                        
                        //for(i = 0; i < n; i++)
                            //out_best_sol[i] = master->sol[i];
                        
                        MRQ_copyArray(n, master->sol, out_best_sol);
                        
                        if(in_store_history_solutions)
                            out_sol_hist.addSolution(n, iter, MRQ_getTime() - timeStart, clock() - clockStart, out_best_sol, out_best_obj);
                        
                        
                        
                        aux = master->setVariableBounds(n, -OPT_INFINITY, zu );
                        
                        #if MRQ_DEBUG_MODE
                            if(aux != 0)
                            {
                                if( in_print_level > 0 )
                                    MRQ_PRINTERRORNUMBER(aux);
                                
                                out_return_code = MRQ_MILP_SOLVER_ERROR;
                                goto termination;
                            }
                        #endif
                    
                    }
                    #endif
                    
                    if( objTemp < zupy )
                        zupy = objTemp;
                    
                    const bool updt = tryUpdateBestSolution( thnumber, n, master->sol, objTemp, iter, clockStart, timeStart, in_store_history_solutions );
                    
                    if( updt )
                    {
                        aux = master->setVariableBounds(n, -OPT_INFINITY, zu );
                        #if MRQ_DEBUG_MODE
                            MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                        #endif
                    }
                    
                }
            }
            
            
            if(in_print_level > 1 && w % in_subiter_printing_frequency == 0)
                printf(MRQ_PREPRINT "%5lu.%-2u  %+-14e  %+-14e  %+-14e\n", iter, w, zl, zu, zu-zl);
            
            if( zupy - zlpy <= MRQ_max(in_absolute_convergence_tol, MRQ_abs(zupy*in_relative_convergence_tol) ) ) //we perform this test after print iteratinprogress because we change the value of w. If we performed this test before printing, w would be printed incorrect
            { 
                //so, we solve the continuous nlp relaxation problem fixing y. we can stop this solving. Even if w == 0, in this case we have lower bound equal to uper bound and so, we solve the original MINLP problem
                w = nintloop;
                //std::cout << "*************************************Resolvi Py!!!\n";
                pyCompletelySolved = true;
                pyOptimalFounded = true;
            }
            
            if( checkTerminationCriterions(thnumber, zl, zu, iter, timeStart, clockStart, out_return_code) )
            {
                goto endLoop;//break;
            }
            
            
            #if 0
            {
                int inds[m];
                int r, mmaster;
                
                r = master->getNumberOfConstraints(mmaster);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                for(int i = nMasterConstrBeforeLoop; i < mmaster; i++)
                    inds[i - nMasterConstrBeforeLoop] = i;
                
                if( mmaster - nMasterConstrBeforeLoop > 0 )
                {
                    r = master->removeConstraints( mmaster - nMasterConstrBeforeLoop, inds );
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                }
            }
            #endif
            
            aux = master->getNumberOfConstraints(lastSetLinearizationIndex);
            MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
            
            aux = masterMilp.addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, false, mastersol,  setQuadsInMaster, in_constr_linearisation_strategy, auxCEval, NULL, NULL, NULL, auxConstr, true );
            MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
            
            aux = masterMilpRelax.addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, false, mastersol,  setQuadsInMaster, in_constr_linearisation_strategy, auxCEval, NULL, NULL, NULL, auxConstr, true );
            MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
            
            
            if( linearizeObj  ) //nonlinear objective
            {
                //by now, we always linearize the last point in the master problem
                
                aux = masterMilp.addLinearizedObjFunction( false, master->sol, setQuadsInMaster, auxCols, auxVars, (objCalculated ? &objTemp : NULL) );
                MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                aux = masterMilpRelax.addLinearizedObjFunction( false, master->sol, setQuadsInMaster, auxCols, auxVars, (objCalculated ? &objTemp : NULL) );
                MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                
                if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS )
                {
                    int r, mmaster;
                    
                    //now, we check if we can discard the linearization on previous point...
                    
                    //auxVars has the gradient of obj linearization, i.e., the coeficients of obj linearizations constraint and the rhs at the end...
                    r = laps.updateObjLinearizationsByNonObjCutPointsByNewPoint( *master, zu, auxVars );
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                    
                    r = laps.addPoint(n, master->sol);
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MEMORY_ERROR, termination);
                    
                    r = master->getNumberOfConstraints( mmaster);
                    #if MRQ_DEBUG_MODE 
                        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                    #endif
                    
                    laps.indMaster[ laps.npoints-1 ] = mmaster-1;
                }
                
            }
            
            
            if(w == 0)
            {
                aux = master->getNumberOfConstraints(nMasterConstrAfterFisrtLinearization);
                MRQ_IFERRORGOTOLABEL(aux, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
                //if we are meausring nlp time and user choose fix int vars in master, we count the time. We considr nlp roslution start from after linearizing master cuts and finish after this internal loop
                if(in_measure_nlp_time && nintloop > 1)
                {
                    nlpBeginTime = MRQ_getTime();
                    nlpBeginClock = clock();
                }
            }
            
            
            if( w < nintloop -1 ) //in the last iteration of this for, we do not have to solve master problem
            {
                //fixing integer variables
                MRQ_fixIntVarsOnSolByList(nI, intVars, mastersol, *masterLp);
                
                aux = masterLp->solve(false);
                
                /*if(aux == OPT_INFEASIBLE_PROBLEM)
                    std::cout << "problema mestre inviavel na fixacao de int vars. iter: " << iter << "\n"; */
                
                //unfixing integer variables
                //MRQ_unfixIntegerVarsByList( nI, intVars, lx, ux, *master );
                
                if( aux == OPT_INFEASIBLE_PROBLEM)
                {
                    //std::cout << "*************************************Py Inviavel!!!\n";
                    pyCompletelySolved = true;
                    break;
                }
                else if( aux != OPT_OPTIMAL_SOLUTION )//&& aux != OPT_INFEASIBLE_PROBLEM )
                {
                    MRQ_PRINTERRORMSGP(" lp approximation for py returned code: ", aux);
                    break;
                }
                
                zlpy = masterLp->objValue;
                mastersol = masterLp->sol;
            }
            
        }
        
        if( !isnan(nlpBeginTime) ) //so, we are measuring nlp time
        {
            out_clock_time_of_nlp_solving += MRQ_getTime() - nlpBeginTime;
            out_cpu_time_of_nlp_solving += MRQ_calcCPUTtime(nlpBeginClock, clock());
            
            nlpBeginTime = NAN;
        }
        
        
        if( pyCompletelySolved )
            out_number_of_subproblems_py_completely_solved++;
        else
            out_all_subproblems_py_completely_solved = false; //when we detect optimal solution or other termination condition, we skip this lines, but I think there is no problem about that.
        
        
        if( (pyOptimalFounded || in_delete_intermediate_linearizations_in_each_iteration) && lastSetLinearizationIndex > nMasterConstrAfterFisrtLinearization )
        {
            int startIndicesToDelete, endIndicesToDelete; 
            
            if( pyOptimalFounded )
            {
                //if we found the optimal solution of Py, we can discard all previous linearization from this iteration, inclusive the original ECP linearization (linearization from integer master problem)
                startIndicesToDelete = firstSetLinearizationIndex;
                endIndicesToDelete = lastSetLinearizationIndex-1;
            }
            else if( in_delete_intermediate_linearizations_in_each_iteration )
            {
                //if we do not reach the optimal solution of Py, we must keep the original ECP linearization (linearization from integer master problem)
                startIndicesToDelete = nMasterConstrAfterFisrtLinearization;
                endIndicesToDelete = lastSetLinearizationIndex-1;
            }
            else
            {
                assert(false); //if we reach this point, you forget to set ninds for an aditional condition
            }
            
            /*const int ninds = endIndicesToDelete - startIndicesToDelete;
            int inds[ninds];
            
            MRQ_setSequentialValuesArray(inds, startIndicesToDelete, endIndicesToDelete);
            
            int r = master->removeConstraints(ninds, inds);*/
            
            
            int r = master->removeConstraintsByRange(startIndicesToDelete, endIndicesToDelete);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
            
            r = masterLp->removeConstraintsByRange(startIndicesToDelete, endIndicesToDelete);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
        #endif
    }
    
    
    
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
    
    
    if( in_refine_final_solution_using_nlp && out_best_obj < MRQ_INFINITY )
    {
        int r;
        double NLPCpuTime, NLPClockTime;
        
        //here, we do not preprocess. Since we are already fixing integer variable, we let it to nlp solver preprocesor, although some solvers do not have one.
        
        nlp = OPT_newNLPSolver(in_nlp_solver);
        MRQ_IFERRORGOTOLABEL(!nlp, out_return_code, out_return_code, termination); //we do not change the return code
        
        r = MRQ_setNLPRelaxProb( prob, lx, ux, plc, puc, nlp, true, true, true, false, thnumber, in_set_special_nlp_solver_params, nlpSolverParams, in_number_of_threads, in_max_cpu_time, in_max_time,  0, 0 );
        
        MRQ_IFERRORGOTOLABEL(r, out_return_code, out_return_code, termination); //we do not change
        
        if( run_by_inside )
        {
            if( !std::isinf(insideSolverMaxTime) )
                nlp->setMaxTime(insideSolverMaxTime );
        }
        
        MRQ_fixIntVarsOnSolByList(nI, intVars, out_best_sol, *nlp);
        
        
        nlp->setInitialSolution(out_best_sol, NULL, NULL);
        
        nlp->solveAndGetTime( in_measure_nlp_time ? &NLPCpuTime : NULL, in_measure_nlp_time ? &NLPClockTime : NULL, false);
        
        if( in_measure_nlp_time )
        {
            out_cpu_time_of_nlp_solving += NLPCpuTime;
            out_clock_time_of_nlp_solving += NLPClockTime;
        }
        
        out_number_of_nlp_probs_solved++;
        
        
        if( nlp->retCode == OPT_OPTIMAL_SOLUTION )
        {
            if(in_print_level > 4)
                std::cout << MRQ_PREPRINT "Refinement problem solved. Old Objective: " << out_best_obj << " New objective: " << nlp->objValue << "\n";
            out_best_obj = nlp->objValue;
            MRQ_copyArray(n, nlp->sol, out_best_sol);
        }
        else
        {
            MRQ_PRINTERRORMSGP("Error to solve the nlp refinement problem: ", nlp->retCode);
        }
        
    }
    
    
    
    
    
termination:
    
    
    
    #if MRQ_CHECK_INT_SOLS_ARE_REPEATING
    if( out_return_code == MRQ_OPTIMAL_SOLUTION && out_number_of_integer_solution_repetitions > 0)
    { //if the optmal solution is repetead, we desconsiderate this, because optimal solution can repeat
        
        for(int i = 0; i < nI; i++)
            auxHistorySol[i] = round( out_best_sol[ intVars[i] ]  );
        
        //I think it is better run in reverse way because I think it is more probable repeted the last solutions
        for(unsigned int k2 = out_integer_solution_history.getnsols(); k2 > 0 ; k2-- )
        {
            unsigned int k = k2-1; //we put here because k and k2 are unsigned
            
            bool equalSols = true;
            MRQ_HistorySolution *hssol = out_integer_solution_history.getHistSolPointer(k);
            
            double *hsol = hssol->sol;
            
            for(int i = 0; i < nI; i++)
            {
                if( hsol[i] != auxHistorySol[i] )
                {
                    equalSols = false;
                    break;
                }
            }
            
            if(equalSols)
            {
                //we found optimal solution in the history
                if( hssol->time > 1.0 )
                    out_number_of_integer_solution_repetitions--;
                //std::cout << "Desconsiderando uma repetição da solução otima!\n";
                break;
            }
            
        }
    }
    
    if(auxHistorySol)   free(auxHistorySol);
    
    #endif
    
    
    if(plc)				free(plc);
    if(intVars)			free(intVars);
    if(auxCols)			free(auxCols);
    if(auxVars)			free(auxVars);
    if(auxConstr)		free(auxConstr);
    if(masterSolConstr)	free(masterSolConstr);
    if(auxCEval)		free(auxCEval);
    
    if(lpsol)			free(lpsol);
    if(bestSolPy)		free(bestSolPy);
    
    if(nlp)				delete nlp;
    
    out_feasible_solution = out_best_obj < MRQ_INFINITY;
    out_number_of_iterations = iter;
    out_lower_bound = zl;
    out_upper_bound = zu;
    
    algorithmFinalization(1, prob, lx, ux);
    
    out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
    out_clock_time = MRQ_getTime() - timeStart;
    
    if(in_print_level > 1)
        std::cout << "cpu time: " << out_cpu_time << "\n";
    
    
    return out_return_code;
}



