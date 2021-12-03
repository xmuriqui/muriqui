/*That file contains a simple implementation of
* Outer Approximation Algorithm, in its basic formulation.
*
* References:
*
* Bonami et al, An algorithm framework for convex mixed integer nonlinear
* programs. Discrete Optimization 5 (2008), pages 186-204.
*
* Duran & Grossmann, An outer-approximation algorithm for a class of
* mixed-integer nonlinear programs, Mathematical Programming 36 (1986),
* pages 307-339.
*
* Author: Wendel Alexandre Melo
*
* Date: 06-Sept-2013
*/

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>
#include <iostream>
#include <new>

#include "BBL_tools.hpp"

#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_solvers.hpp"





//using namespace std;

using namespace optsolvers;
using namespace muriqui;





MRQ_OuterApp::MRQ_OuterApp():MRQ_LinearApproxAlgorithm()
{
    resetParameters();
    resetOutput();
    out_algorithm = MRQ_OA_ALG;
}



MRQ_OuterApp::~MRQ_OuterApp()
{
}


void MRQ_OuterApp::printParameters(std::ostream &out) const
{
    MRQ_LinearApproxAlgorithm::printParameters(out);
    out << "\n"
    
    MRQ_STRFFATT(in_binarie_cut_when_nlp_infeasible) << "\n"
    MRQ_STRFFATT(in_round_first_nlp_relaxation_solution) << "\n"
    MRQ_STRFFATT(in_use_first_nlp_relaxation) << "\n"
    ;
}


void MRQ_OuterApp::resetParameters()
{
    MRQ_LinearApproxAlgorithm::resetParameters();
    
    in_binarie_cut_when_nlp_infeasible = false;
    in_round_first_nlp_relaxation_solution = false;
    in_use_first_nlp_relaxation = true;
    //in_max_iterations = 1000;
    
}



int MRQ_OuterApp::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_LinearApproxAlgorithm::setIntegerParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_binarie_cut_when_nlp_infeasible), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_round_first_nlp_relaxation_solution), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_first_nlp_relaxation), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}



//note: we do not perform test here, just update the values...
void MRQ_OuterApp::updateBestSol( const int n, double *newsol, const double *dualSol, const double objValue, optsolvers::OPT_LPSolver *master, MRQ_LAAPointsStoring *laps, const bool linearizeObj, MRQ_QUAD_APP_MASTER_STRATEGY quadAppStrategy, MRQ_IntegratedHessian *intHess, double *auxVars, const clock_t &clockStart, const int timeStart, const long unsigned int iter )
{
    
    //out_best_obj = objValue;
    //zu = objValue;
    
    //MRQ_copyArray(n, newsol, out_best_sol);
    
    
    //note we do not store history solutions here. This store is made by run method...
    const bool updt = tryUpdateBestSolution( thnumber, n, newsol, objValue, iter, clockStart, timeStart, false );
    
    
    if( updt )
    {
        int aux;
        double zuaux;
        
        if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS && linearizeObj)
        {
            aux = laps->updateObjLinearizationsByNonObjCutPoints2(*master, zu, auxVars);
            
            if(aux != 0)
            {
                std::cerr << MRQ_PREPRINT << "Error " << aux << MRQ_GETFILELINE << std::endl;
            }
        }
        
        zuaux = quadAppStrategy == MRQ_QAMS_NO_QUAD_APP ? zu : MRQ_zuWithTol(zu, in_absolute_convergence_tol, in_relative_convergence_tol); //to quadratic app work, we need to set a value smaller than zu...
        
        aux = master->setVariableBounds(n, -OPT_INFINITY, zuaux );//do not pass zl here because it damages branch-and-bound performance
        
        #if MRQ_DEBUG_MODE
            if(aux != 0)
            {
                if( in_print_level > 0 )
                    MRQ_PRINTERRORNUMBER(aux);
            
                //out_return_code = MRQ_UNDEFINED_ERROR;
                //goto termination;
            }
        #endif
        
        if( quadAppStrategy == MRQ_QAMS_NO_QUAD_APP )
        {
            master->setObjCutUpperBound(zu);
        }
        else if( quadAppStrategy == MRQ_QAMS_ON_BEST_POINT )
        {
            OPT_QPSolver *qmaster = (OPT_QPSolver *) master;
            
            
            aux = intHess->evalCompleteHessian( thnumber, true, newsol, 1.0, dualSol );
            
            if( aux != 0 )
            {
                if( in_print_level > 0 )
                    MRQ_PRINTCALLBACKERRORNUMBER(aux);
            }
            else
            {
                aux = qmaster->setObjQuadMatrix( intHess->nzs, intHess->rows, intHess->cols, intHess->values );
                
                if( aux != 0 )
                {
                    if( in_print_level > 0 )
                        MRQ_PRINTERRORNUMBER(aux);
                }
                
                //std::cout << "Atualizei a matrix quadratica da funcao objetivo no master!" << std::endl;
            }
            
        }
        
    }
    
}



int MRQ_OuterApp::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    const bool nlObj = prob.hasNlObj;
    const bool binProblem = prob.isBinaryProblem();
    const int n = prob.n;
    const int m = prob.m;
    const int nI= prob.getNumberOfIntegerVars();
    const bool preproc = in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_obj_function;
    
    
    bool setQuadsInMaster = in_set_quadratics_in_master_problem;
    bool useFirstNlp = in_use_first_nlp_relaxation;
    bool linearizeObj;
    bool feasProbSol = false;
    bool objCalculated = false, constrCalculated = false;
    bool binaryCut = false;
    bool updtConstrBounds;
    
    int aux;
    MRQ_QUAD_APP_MASTER_STRATEGY quadAppStrategy = in_quad_app_master_strategy;
    unsigned long int iter = 0;
    int *intVars = NULL;
    double *flx = NULL, *fux;
    double *plc = NULL, *puc = NULL;
    double timeStart, objTemp, zuaux;
    double NLPCpuTime, NLPClockTime;
    double *pNLPCpuTime = in_measure_nlp_time ? &NLPCpuTime : NULL;
    double *pNLPClockTime = in_measure_nlp_time ? &NLPClockTime : NULL;
    
    double *auxPSol, *auxPConstr;
    double *newCalclc = NULL, *newCalcuc;
    clock_t clockStart;
    
    
    MRQ_NLPSolver *nlp = NULL;
    MRQ_NLPFeasProb nlpFeas2;
    //MRQ_NLPSolver *nlpFeas = NULL;
    MRQ_GradientsEvaluation  gradEval;
    MRQ_LAAPointsStoring laps(n);
    MRQ_Preprocessor preprocessor(&prob);
    
    MRQ_IntegratedHessian *intHess = NULL;
    
    MRQ_MasterMILPProb masterMilp;
    OPT_LPSolver *master;
    OPT_NLPSolver *constrBoxCalcNLP = NULL;
    
    
    //araays:
    bool *auxConstrEval = NULL;
    int *auxCols = NULL;
    double *auxVars = NULL, *auxConstr = NULL, *dualSol = NULL;
    double *masterSolConstr = NULL;
    
    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux; 
    
    
    timeStart = MRQ_getTime();
    clockStart = clock();
    
    
    nthreads = 1;
    nthreads_lazy = in_number_of_threads > 0 ? in_number_of_threads : branchAndBound::BBL_getNumCores() ;
    if( in_milp_solver == MRQ_GUROBI )
        nthreads_lazy = 1; //gurobi apply multithreading to solve the problem, but only thread 0 is called to add lazy constraints.
    
    
    {
        auto r = algorithmInitialization(1, (in_preprocess_lin_constr || in_preprocess_obj_function || in_preprocess_quad_constrs), milpSolverParams, nlpSolverParams, prob, lx, ux, &preprocessor, &updtConstrBounds, &plc, &puc); //that algorithm is monothread...
        if(r != MRQ_SUCCESS)
        {
            if(in_print_level > 0)
            {
                if( r == MRQ_INFEASIBLE_PROBLEM )
                    std::cout << MRQ_PREPRINT "Preprocessor detected infeasible problem\n";
                else
                    MRQ_PRINTERRORMSG("Error at algorithm initialization\n");
            }
                
            out_return_code = r;
            goto termination;
        }
        
    }
    
    
    
    if(in_print_level > 1)
    {
        std::cout << "\n";
        MRQ_PRINTMSG("Starting Outer Approximation algorithm\n\n");
    }
    
    if(in_print_level > 3)
        printSubSolvers(true, true, false);
    
    
    aux = gradEval.initialize(thnumber, &prob);
    if(aux != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(aux);
        
        out_return_code = MRQ_MEMORY_ERROR;
        goto termination;
    }
    
    
    intVars = (int *) malloc( nI * sizeof(int) );
    auxConstrEval = (bool *) malloc( m * sizeof(bool) );
    //auxConstrEval2 = (bool *) malloc( m * sizeof(bool) );
    auxConstr = (double *) malloc( m * sizeof(double) );
    auxCols = (int *) malloc( (n+1)*sizeof(int) );
    auxVars = (double *) malloc( (n+1) *sizeof(double) );
    dualSol = (double *) malloc( m*sizeof(double) );
    //nlp = MRQ_newNlpSolver(in_nlp_solver, thnumber, &prob);
    
    nlp = OPT_newNLPSolver(in_nlp_solver);
    
    if(!nlp || !intVars || !auxConstrEval || !auxConstr || !auxCols || !auxVars || !dualSol)
    {
        if(in_print_level > 0)
            MRQ_PRINTMEMERROR;
        
        out_return_code = MRQ_MEMORY_ERROR;
        goto termination;
    }
    
    
    prob.getIntegerIndices(intVars);
    
    for(int i = 0; i < m; i++)
        auxConstrEval[i] = prob.nlConstr[i] || prob.QC[i].getNumberOfElements() > 0;
    
    
    if( in_constr_linearisation_strategy == MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_MASTER_SOL_ALSO )
    {
        masterSolConstr = (double *) malloc( m * sizeof(double) );
        if( !masterSolConstr )
        {
            if(in_print_level > 0)
                MRQ_PRINTMEMERROR;
            
            out_return_code = MRQ_MEMORY_ERROR;
            goto termination;
        }
    }
    else if( in_constr_linearisation_strategy == MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_BY_BOX_FILL )
    {
        constrBoxCalcNLP = OPT_newNLPSolver(in_nlp_solver);
        
        newCalclc = (double *) malloc( 2*m * sizeof(double) );
        
        if(!constrBoxCalcNLP || !newCalclc)
        {
            if(in_print_level > 0)
                MRQ_PRINTMEMERROR;
            
            out_return_code = MRQ_MEMORY_ERROR;
            goto termination;
        }
        
        newCalcuc = &newCalclc[m];
        
        clock_t cs = clock();
        
        //calculating bounds for all nonlinear constraints
        int r = OPT_ConstrsBoundsCalculator:: calculate(prob, plc, puc, constrBoxCalcNLP, nlpSolverParams, false, auxConstrEval, newCalclc, newCalcuc);
        
        out_cpu_time_on_box_to_constr_calculation = MRQ_calcCPUTtime(cs, clock());
        
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
            
            out_return_code = MRQ_NLP_SOLVER_ERROR;
            goto termination;
        }
        
        /*std::cout << "Caixa das restricoes: \n";
        
        for(int i = 0; i < m; i++)
        {
            if( auxConstrEval[i] )
            {
                std::cout << "i: " << i << " lc: " << newCalclc[i] << " uc: " << newCalcuc[i] << "\n";
            }
        }
        
        MRQ_getchar(); */
        
        delete constrBoxCalcNLP;
        constrBoxCalcNLP = NULL;
    }
    
    
    if( preproc )
    {
        flx = (double *) malloc( 2*n * sizeof(double) );
        if( !flx )
        {
            if(in_print_level > 0)
                MRQ_PRINTMEMERROR;
            
            out_return_code = MRQ_MEMORY_ERROR;
            goto termination;
        }
        
        fux = &flx[n];
    }
    
    
    aux = masterMilp.setProblemBase( thnumber, prob, in_milp_solver, true, setQuadsInMaster, true, lx, ux, 1, milpSolverParams, in_number_of_threads );
    if( aux != 0 )
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(aux);
        
        out_return_code = MRQ_MILP_SOLVER_ERROR;
        goto termination;
    }
    
    master = masterMilp.master;
    
    
    aux = master->getSolverType();
    if( aux == OPT_LP || ( aux == OPT_QP && prob.hasQuadMatrixInSomeConstraint() )  )
    {
        if( in_set_quadratics_in_master_problem && in_print_level > 0 )
            std::cerr << MRQ_PREPRINT "Warning: parameter in_set_quadratics_in_master_problem is set to true, but milp solver does not support quadratic constraints currently. Changing the parameter to false.\n";
        
        if( in_quad_app_master_strategy != MRQ_QAMS_NO_QUAD_APP && in_print_level > 0 )
            std::cerr << MRQ_PREPRINT "Warning: parameter in_quad_app_master_strategy is set to use quadratic appoximation, but milp solver does not support quadratic constraints currently. Changing to linear master problem.\n";
        
        setQuadsInMaster = false; //integer solver does not support the quadratics
        quadAppStrategy = MRQ_QAMS_NO_QUAD_APP;
    }
    
    linearizeObj = nlObj || (prob.Q.getNumberOfElements() > 0 && !setQuadsInMaster);
    
    
    if( quadAppStrategy != MRQ_QAMS_NO_QUAD_APP  )
    {
        intHess = new (std::nothrow) MRQ_IntegratedHessian (&prob);
        
        if( !intHess )
        {
            if( in_print_level > 0 )
                MRQ_PRINTMEMERROR;
            
            out_return_code = MRQ_MEMORY_ERROR;
            goto termination;
        }
        
        aux = intHess->buildStructures(true, true);
        if( aux != 0 )
        {
            if( in_print_level > 0 )
                MRQ_PRINTERRORNUMBER(aux);
            
            out_return_code = MRQ_MEMORY_ERROR;
            goto termination;
        }
        
        useFirstNlp = true; //we need lambda multipliers
    }
    
    aux = MRQ_setNLPRelaxProb( prob, lx, ux, plc, puc, nlp, true, true, true, false, thnumber, in_set_special_nlp_solver_params, nlpSolverParams, in_number_of_threads, in_max_cpu_time, in_max_time, 0, 0 );
    
    if( aux != 0 )
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(aux);
        
        out_return_code = MRQ_NLP_SOLVER_ERROR;
        goto termination;
    }
    
    if( run_by_inside )
    {
        if( !std::isinf(insideSolverMaxTime) )
            aux += nlp->setMaxTime(insideSolverMaxTime );
    }
    
    
    /* aux = nlp->setnVariablesBounds(n, lx, ux);
    if( updtConstrBounds )
    {
        for(int i = 0; i < m; i++)
            aux += nlp->setConstraintBounds(i, plc[i], puc[i]);
    }*/
    
    
    if(aux != 0)
    {
        if( in_print_level > 0 )
            MRQ_PRINTERRORMSG("Error to set NLP relaxation!");
        
        out_return_code = MRQ_NLP_SOLVER_ERROR;
        goto termination;
    }
    
    
    if( useFirstNlp || nPoints == 0 )
    {
        if( nlp->getSolverCode() != optsolvers::OPT_MOSEK )
        {
            //aux += nlp->setInitialSolution();
            
            if(in_use_initial_solution)
            {
                aux = nlp->setInitialSolution( xInit, NULL, NULL);
            
                if(aux != 0)
                {
                    if( in_print_level > 0 )
                        std::cerr << MRQ_PREPRINT "Error " << aux << " at setting initial solution" << MRQ_GETFILELINE << "\n";
                    
                    //out_return_code = MRQ_UNDEFINED_ERROR;
                    //goto termination;
                }
            }
        }
        
        if(in_print_level > 2)
            std::cout << MRQ_PREPRINT "Solving NLP relaxation\n";
        
        
        
        nlp->solveAndGetTime(pNLPCpuTime, pNLPClockTime, false);
        
        if( in_measure_nlp_time )
        {
            out_cpu_time_of_nlp_solving += *pNLPCpuTime;
            out_clock_time_of_nlp_solving += *pNLPClockTime;
        }
        
        out_number_of_nlp_probs_solved++;
        
        
        if( nlp->retCode == OPT_OPTIMAL_SOLUTION )
        {
            out_obj_opt_at_continuous_relax = nlp->objValue;
            zl = MRQ_max( zl,  nlp->objValue);
            
            if(in_print_level > 1)
                std::cout << MRQ_PREPRINT "NLP relaxation solution: " << nlp->objValue << "\n";
            
            //for(int i = 0; i < n; i++)
                //printf("nlp[%d]->sol: %f\n", i, nlp->sol[i]);
            
            
            if(in_store_history_solutions)
                out_sol_hist.addSolution(n, iter, MRQ_getTime() - timeStart, MRQ_calcCPUTtime(clockStart), nlp->sol, nlp->objValue);
            
            if( zl > zu )
            {
                if(in_print_level > 0)
                    MRQ_PRINTMSG("Solution of NLP relaxation is greater than upper_bound ");

                out_return_code = MRQ_INFEASIBLE_PROBLEM;
                goto termination;
            }
            
            //if( prob.isIntegerSolution(nlp->sol, in_integer_tol)  )
            if( MRQ_isIntegerSol( nI, intVars, nlp->sol, in_integer_tol ) )
            {
                if( in_print_level > 1 )
                    MRQ_PRINTMSG("An integer optimal solution was gotten as NLP relaxation solution\n");
                
                //MRQ_copyArray( n, nlp->sol, out_best_sol);
                //out_best_obj = zu = nlp->objValue;
                
                tryUpdateBestSolution( thnumber, n, nlp->sol, nlp->objValue, 0, clockStart, timeStart, false ); //we already store this solution in the history if in_store_history_solutions is true
                
                
                out_return_code = MRQ_OPTIMAL_SOLUTION;
                
                goto termination;
            }
            
        }
        else if( nlp->retCode == OPT_INFEASIBLE_PROBLEM)
        {
            if( in_print_level > 2 )
                std::cout << MRQ_PREPRINT "Continuous relaxation is infeasible!\n";
            
            out_return_code = MRQ_INFEASIBLE_PROBLEM;
            goto termination;
        }
        else if( nlp->retCode == OPT_UNBOUNDED_PROBLEM )
        {
            if( in_print_level > 2 )
                std::cout << MRQ_PREPRINT "Continuous relaxation is unbounded!\n";
            
            out_return_code = MRQ_UNBOUNDED_PROBLEM;
            goto termination;
        }
        else
        {
            if(in_print_level > 0)
                std::cerr << MRQ_PREPRINT "Failure at solving continuous relaxation!\n";
        }
        
        
        if( in_round_first_nlp_relaxation_solution )
        {
            if( nlp->feasSol ) //( nlp->retCode == OPT_OPTIMAL_SOLUTION || nlp->retCode == OPT_FEASIBLE_SOLUTION )
            {
                double* const sol = nlp->sol;
                
                for(int i = 0; i < nI; i++)
                {
                    const int &ind = intVars[i];
                    sol[ind] = round(sol[ind]);
                }
            }
        }
        
    }
    
    
    
    if( nPoints > 0 )
    {
        aux += masterMilp. addConstraintLinearisationPointsToMILP( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, nPoints, points, setQuadsInMaster, in_constr_linearisation_strategy, auxConstrEval);
        
        if(linearizeObj)
            aux += masterMilp.addObjLinearisationPointsToMILP( nPoints, points, zu, setQuadsInMaster, in_obj_linearisation_strategy, NULL, NULL, &laps );
    }
    
    
    if(useFirstNlp || nPoints == 0)
    {
        bool newx = true;
        int r;
        
        for(int i = 0; i < m; i++)
            //auxConstrEval[i] = true;
            auxConstrEval[i] = prob.nlConstr[i] || (prob.QC[i].getNumberOfElements() > 0 && !setQuadsInMaster) ;
        
        if( !nlp->feasSol || in_round_first_nlp_relaxation_solution ) //if( ( nlp->retCode != OPT_OPTIMAL_SOLUTION && nlp->retCode != OPT_FEASIBLE_SOLUTION )  ||  in_round_first_nlp_relaxation_solution  )
        {
            
            if( !nlp->feasSol ) //( nlp->retCode != OPT_OPTIMAL_SOLUTION && nlp->retCode != OPT_FEASIBLE_SOLUTION )
            {
                double *sol = nlp->sol;
                
                #pragma ivdep
                #pragma GCC ivdep
                for(int i = 0; i < n; i++)
                    sol[i] = MRQ_min(ux[i], MRQ_max(0.0, lx[i]) );
            }
            
            
            r = prob.objEval( thnumber, newx, nlp->sol, nlp->objValue );
            if( r != 0 )
            {
                if( in_print_level > 0 )
                    std::cerr << MRQ_PREPRINT <<"Callback function error " << r << MRQ_GETFILELINE << std::endl;
                
                out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
                goto termination;
            }
            
            
            if( prob.hasNlObj )
                newx = false;
            
            
            r = prob.constraintsEval(thnumber, newx, auxConstrEval, nlp->sol, nlp->constr);
            if( r != 0 )
            {
                if( in_print_level > 0 )
                    std::cerr << MRQ_PREPRINT <<"Callback function error " << r << MRQ_GETFILELINE << std::endl;
                
                out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
                goto termination;
            }
            
            
            if( prob.hasNlConstrs )
                newx = false;
        }
        
        aux +=  masterMilp .addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, newx, nlp->sol, setQuadsInMaster, in_constr_linearisation_strategy, auxConstrEval, NULL, NULL, NULL, nlp->constr, true, NULL, newCalclc, newCalcuc );
        
        
        if(prob.hasNlConstrs)
            newx = false;
        
        
        if( linearizeObj  )
        {
            int mmaster;
            
            aux += masterMilp. addLinearizedObjFunction( newx, nlp->sol, setQuadsInMaster, auxCols, auxVars, &(nlp->objValue) );
            
            if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS )
            {
                //auxVars has the coeficients of linearization to objective function and the RHS also
                
                
                r = laps.updateObjLinearizationsByNonObjCutPointsByNewPoint( *master, zu, auxVars );
                if( r != 0 )
                {
                    if( in_print_level > 0 )
                        MRQ_PRINTERRORNUMBER(r);
                    
                    out_return_code = MRQ_MILP_SOLVER_ERROR;
                    goto termination;
                }
                
                
                r = laps.addPoint(n, nlp->sol );
                if( r != 0 )
                {
                    if( in_print_level > 0 )
                        MRQ_PRINTERRORNUMBER(r);
                    
                    out_return_code = MRQ_MEMORY_ERROR;
                    goto termination;
                }
                
                
                r = master->getNumberOfConstraints( mmaster );
                #if MRQ_DEBUG_MODE
                    if(r != 0)
                    {
                        if(in_print_level > 0)
                            MRQ_PRINTERRORNUMBER(r);
                        
                        out_return_code = MRQ_MILP_SOLVER_ERROR;
                        goto termination;
                    }
                #endif
                
                laps.indMaster[ laps.npoints-1 ] = mmaster;
            }
        }
    }

    
    if( run_by_inside )
    {
        if( !std::isinf(insideSolverMaxTime) )
            aux += master->setMaxTime( insideSolverMaxTime );
    }
    
    
    if( zu < MRQ_INFINITY )
    {
        zuaux = quadAppStrategy == MRQ_QAMS_NO_QUAD_APP ? zu : MRQ_zuWithTol(zu, in_absolute_convergence_tol, in_relative_convergence_tol); //to quadratic app work, we need to set a value smaller than zu...
        
        aux += master->setVariableBounds( n, -OPT_INFINITY, zuaux); //auxiliary variable. do not pass zl here because it damages branch-and-bound performance
        
        if( quadAppStrategy == MRQ_QAMS_NO_QUAD_APP )
            master->setObjCutUpperBound( zu );
    }
    
    
    if(aux != 0)
    {
        if( in_print_level > 0 )
            MRQ_PRINTERRORMSG("Error at setting MILP relaxation!");
        
        out_return_code = MRQ_MILP_SOLVER_ERROR;
        goto termination;
    }
    
    
    if( quadAppStrategy != MRQ_QAMS_NO_QUAD_APP )//( quadAppStrategy == MRQ_QAMS_ON_LAST_POINT )
    {
        if( nlp->retCode == OPT_OPTIMAL_SOLUTION )
        {
            OPT_QPSolver *qmaster = (OPT_QPSolver *) master;
            
            
            nlp->getDualSolution(auxConstr, NULL, true);
            
            aux = intHess->evalCompleteHessian( thnumber, true, nlp->sol, 1.0,  auxConstr );
            
            if( aux != 0 )
            {
                if( in_print_level > 0 )
                    std::cerr << MRQ_PREPRINT "Callback function error " << aux << MRQ_GETFILELINE << "\n";
                
                #if MRQ_DEBUG_MODE
                    out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
                    goto termination;
                #endif
            }
            
            aux = qmaster->setObjQuadMatrix( intHess->nzs,  intHess->rows, intHess->cols, intHess->values ); //it is missing the therm -x⁰ * H * x⁰, but that is a constant and it is not really necessary
            
            if( aux != 0 )
            {
                if( in_print_level > 0 )
                    std::cerr << MRQ_PREPRINT "Callback function error " << aux << MRQ_GETFILELINE << "\n";
                
                #if MRQ_DEBUG_MODE
                    out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
                    goto termination;
                #endif
            }
            
        }
        else
        {
            if( in_print_level > 0 )
                MRQ_PRINTERRORMSG("Failure to solve nlp relaxation. Adopting the first master problem as linear.");
        }
        
    }
    
    
    if(in_print_level > 1)
        printf("%-5s  %-14s  %-14s  %-14s\n", "iter", "    lbound", "    ubound", "     gap");
    
    
    
    for(int i = 0; i < m; i++)
        auxConstrEval[i] = prob.nlConstr[i] || (prob.QC[i].getNumberOfElements() > 0 && !setQuadsInMaster) ;
    
    
    //starting outer approximation loop
    while(true)
    {
        iter++;
        
        #if MRQ_DEBUG_MODE
            auxPConstr = NULL;
            auxPSol  = NULL;
        #endif
        
        //master->generateModelFile("master2.lp");
        //std::cout << "gerei master2.lp" << std::endl;
        //getchar();
        
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
                
                //out_return_code = MRQ_MILP_SOLVER_ERROR;
                //break;
            }
        }
        
        
        
        
        
        if( in_set_obj_lower_bound_on_master_problem )
            master->setVariableBounds(n, zl, zu);
        
        aux = master->solve(false);
        
        {  //couting number of iterations
            long unsigned int masterIters;
            int r = master->getNumberOfIterations(masterIters);
            if(r == 0)
                out_number_of_milp_solver_iters += masterIters;
            else if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
        }
        
        if( in_constr_linearisation_strategy == MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_MASTER_SOL_ALSO )
        {
            const int r = prob.constraintsEval( thnumber, true, auxConstrEval, master->sol, masterSolConstr );
            
            if( r != 0 )
            {
                if(in_print_level > 0)
                    MRQ_PRINTCALLBACKERRORNUMBER(r);
                
                out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
                goto termination;
            }
        }
        
        
        //printf("MILP status: %d orig status: %d objF: %f\n", aux, master->origSolverRetCode, master->objValue);
        //for(int i = 0; i < n; i++)
            //printf("isol[%d]: %f \t", i, master->sol[i]);
        //getchar();
        
        if( aux != OPT_OPTIMAL_SOLUTION )
        {
            if( aux == OPT_INFEASIBLE_PROBLEM )
                out_return_code = out_best_obj < MRQ_INFINITY ? MRQ_OPTIMAL_SOLUTION : MRQ_INFEASIBLE_PROBLEM ;
            else if( aux == OPT_UNBOUNDED_PROBLEM)
                out_return_code = MRQ_UNBOUNDED_PROBLEM;
            else if( aux == OPT_UNDEFINED_ERROR )
                out_return_code = out_best_obj < MRQ_INFINITY ? MRQ_OPTIMAL_SOLUTION : MRQ_MILP_SOLVER_ERROR ;
            else if( aux == OPT_MAX_TIME )
                out_return_code = MRQ_MAX_TIME_STOP;
            else
            {
                if( in_print_level > 0 )
                    std::cerr << MRQ_PREPRINT " Error at milp solver. Optsolvers return code: " << aux << " Original return code: " << master->origSolverRetCode << "\n";
                out_return_code = MRQ_MILP_SOLVER_ERROR;
            }
            
            break;
        }
        
        
        if( quadAppStrategy == MRQ_QAMS_NO_QUAD_APP )
        {
            if( master->objValue > zl )
            {
                zl = master->objValue; //if we are using quadratic master problem, there is no valid lower bound...
                
                //we check here if we reach the optimal solution to avoid solve a NLP problem
                if(zu - zl <= in_absolute_convergence_tol || zu - zl <= MRQ_abs(zu)*in_relative_convergence_tol)
                {
                    //zl is equal to zu. We can stop
                    if(in_print_level > 1)
                        printf("%-5ld  %+-14e  %+-14e  %+-14e\n", iter, zl, zu, zu-zl);
                    out_return_code = MRQ_OPTIMAL_SOLUTION;
                    break;
                }
                
            }
        }
        
        
        
        if( preproc )
        {
            const double *sol = master->sol;
            bool updtv, updtc;
            
            
            MRQ_copyArray(n, lx, flx);
            MRQ_copyArray(n, ux, fux);
            
            #pragma GCC ivdep
            for(int i = 0; i < nI; i++)
            {
                const int ind = intVars[i];
                flx[ind] = fux[ind] = round( sol[ind] );
            }
            
            
            aux = preprocessor.preprocess( in_preprocess_quad_constrs, in_preprocess_obj_function && zu < MRQ_INFINITY, zu, flx, fux, updtv, updtc, NULL, NULL, plc, puc );
            
            if( aux == minlpproblem::MIP_INFEASIBILITY )
            {
                aux = OPT_INFEASIBLE_PROBLEM;
                nlp->retCode = OPT_INFEASIBLE_PROBLEM;
                nlp->feasSol = false;
            }
            else
            {
                int r = nlp->setnVariablesBounds( n, flx, fux );
                
                for(int i = 0; i < m; i++)
                    r += nlp->setConstraintBounds( i, plc[i], puc[i] );
                
                #if MRQ_DEBUG_MODE
                    if( r != 0 )
                    {
                        if( in_print_level )
                            MRQ_PRINTERROR;
                        out_return_code = MRQ_NLP_SOLVER_ERROR;
                        break;
                    }
                #endif
            }
            
        }
        else
        {
            MRQ_fixIntVarsOnSolByList(nI, intVars, master->sol, *nlp);
            aux = 0;
        }
        
        
        if( aux != OPT_INFEASIBLE_PROBLEM )
        { //preprocessor can detect infeasibility
            nlp->setInitialSolution( master->sol, NULL, NULL );
            
            aux = nlp->solveAndGetTime(pNLPCpuTime, pNLPClockTime, false);
            
            if(in_measure_nlp_time)
            {
                out_cpu_time_of_nlp_solving += *pNLPCpuTime;
                out_clock_time_of_nlp_solving += *pNLPClockTime;
            }
            
            out_number_of_nlp_probs_solved++;
        }
        
        
        /*printf("Resolvi nlp. code: %d nlp->code: %d objValue: %f feas: %d\n ", aux, nlp->retCode, nlp->objValue, (int) nlp->feasSol );
        MRQ_getchar();*/
        
        
        //note, if we only have a feasible solution to nlp, we just linearize on them if they improve the lower bound. Otherwise, we will let add the ecp cut. We do it because if we linearize in a feasible solution that does not improve zu, we have the danger of cycling (ok, we could avoid this cycling, but it would be too many work), so we adotp this strategy 
        if( aux == OPT_OPTIMAL_SOLUTION || (nlp->feasSol && nlp->objValue < zu) )
        {
            #if OPT_DEBUG_MODE
                assert( aux != OPT_INFEASIBLE_PROBLEM );
            #endif
            
            nlp->getDualSolution(dualSol, NULL, true);
            
            if( nlp->objValue < zu )
            {
                updateBestSol(n, nlp->sol, dualSol, nlp->objValue, master, &laps, linearizeObj, quadAppStrategy, intHess, auxVars, clockStart, timeStart, iter);
            }
            
            objCalculated = true;
            constrCalculated = true;
            
            objTemp = nlp->objValue;
            auxPConstr = nlp->constr;
            
            auxPSol = nlp->sol;
        }
        else if( aux == OPT_INFEASIBLE_PROBLEM )
        {
            if( in_print_level > 6 )
                std::cout << MRQ_PREPRINT "integer fixed NLP problem is infeasible.\n";
            
            if( in_binarie_cut_when_nlp_infeasible && binProblem )
            {
                binaryCut = true;
            }
            else
            {
                int r;
                
                if( in_print_level > 4 )
                    std::cout << MRQ_PREPRINT "Solving feasibility problem.\n";
                
                if( !nlpFeas2.solver )
                {
                    r = nlpFeas2.setProblem( in_nlp_solver, prob, lx, ux, nlpSolverParams, thnumber, in_set_special_nlp_solver_params, in_number_of_threads, in_max_cpu_time, in_max_time);
                    
                    if( r != 0 )
                    {
                        if( in_print_level > 0 )
                            MRQ_PRINTERRORNUMBER(r);
                        
                        out_return_code = MRQ_NLP_SOLVER_ERROR;
                        break;
                    }
                    
                    //r = nlp->setnVariablesBounds(n, lx, ux);
                    
                    if( run_by_inside )
                    {
                        if(!std::isinf(insideSolverMaxTime))
                            r += nlpFeas2.solver->setMaxTime(insideSolverMaxTime);
                        
                        if(r != 0)
                        {
                            if( in_print_level > 0 )
                                MRQ_PRINTERROR;
                            
                            out_return_code = MRQ_MEMORY_ERROR;
                            break;
                        }
                    }
                }
                
                
                
                MRQ_fixIntVarsOnSolByList(nI, intVars, master->sol, *nlpFeas2.solver);
                
                
                {//setting initial solution to nlp feasible problem from master solution
                    int nfeas, r;
                    double *xInitNLPFeas = nlpFeas2.solver->sol;
                    
                    r = nlpFeas2.solver->getNumberOfVars(nfeas);
                    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
                    
                    MRQ_copyArray(n, master->sol, xInitNLPFeas);
                    MRQ_setAllArray(nfeas - n, &xInitNLPFeas[n], 0.0 );
                    
                    ( (MRQ_NLPSolver *) nlpFeas2.solver)->setInitialSolution( xInitNLPFeas, NULL, NULL );
                }
                
                r = nlpFeas2.solver->solveAndGetTime( pNLPCpuTime, pNLPClockTime, false);
                
                if( in_measure_nlp_time )
                {
                    out_cpu_time_of_nlp_solving += *pNLPCpuTime;
                    out_clock_time_of_nlp_solving += *pNLPClockTime;
                }
                
                out_number_of_nlp_probs_solved++;
                
                //nlpFeas2.solver->generateModelFile( "oanlpfeasmodel.out" );
                
                /*printf("Resolvi nlpfeas by optsolvers. code: %d orig code: %d objValue: %f\n", r, nlpFeas2.solver->origSolverRetCode,  nlpFeas2.solver->objValue);
                
                for(int i = 0; i < prob.n + prob.m; i++)
                    printf("osol[%d]: %f ", i, nlpFeas2.solver->sol[i] );
                printf("\n");
                
                MRQ_getchar(); */
                
                if( r == OPT_OPTIMAL_SOLUTION )
                {
                    auxPSol = nlpFeas2.solver->sol;
                    nlpFeas2.solver->getDualSolution( dualSol, NULL, true );
                    feasProbSol = true;
                }
                else
                {
                    //we got failure to solve feasibility problem. We try add binary cut.
                    
                    if( in_print_level > 4 )
                        std::cerr << MRQ_PREPRINT "Failure to solve NLP feasibility relaxation. \n";
                    
                    if( binProblem )
                    {
                        binaryCut = true;
                    }
                    else
                    {
                        //our unique choice is add a cut on the milp solution, like ecp
                        auxPSol = master->sol;
                    }
                }
                
                //if( in_print_level > 4 )
                //printf("\n" );
            }
        }
        else
        {
            bool feasible;
            int r;
            //printf("Adicionando corte ECP\n");
            
            if( in_print_level > 4 )
                std::cerr << MRQ_PREPRINT "Failure to solve NLP relaxation fixing integer vars.\n" ;
            
            //check if milp solution is feasible
            
            r = prob.isFeasibleToConstraints(thnumber, master->sol, true, auxConstrEval, in_absolute_feasibility_tol, in_relative_feasibility_tol, feasible, auxConstr);
            
            if(r == 0)
            {
                auxPConstr = auxConstr;
                constrCalculated = true;
            }
            
            
            //our unique choice is add a cut on the milp solution, like ecp
            auxPSol = master->sol;
            MRQ_setAllArray(m, dualSol, 1.0);
            
            
            if(feasible)
            {
                int r = prob.objEval( thnumber, !prob.hasNlConstrs, master->sol, objTemp);
                
                if(r == 0)
                {
                    objCalculated = true;
                    if( objTemp < zu )
                    {
                        updateBestSol(n, master->sol, dualSol, objTemp, master, &laps, linearizeObj, quadAppStrategy, intHess, auxVars, clockStart, timeStart, iter);
                        
                        /*zu = out_best_obj = objTemp;
                        
                        MRQ_copyArray(n, master->sol, out_best_sol);
                        
                        //if(in_store_history_solutions)
                            //out_sol_hist.addSolution(n, iter, MRQ_getTime() - timeStart, clock() - clockStart, out_best_sol, out_best_obj);
                        
                        
                        if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS && linearizeObj)
                        {
                            aux = laps.updateObjLinearizationsByNonObjCutPoints2(*master, zu, auxVars, auxVars2);
                            
                            if(aux != 0)
                            {
                                std::cerr << MRQ_PREPRINT << "Error " << aux << MRQ_GETFILELINE << std::endl;
                            }
                            
                        }
                        
                        
                        zuaux = quadAppStrategy == MRQ_QAMS_NO_QUAD_APP ? zu : MRQ_zuWithTol(zu, in_absolute_convergence_tol, in_relative_convergence_tol); //to quadratic app work, we need to set a value smaller than zu... 
                        
                        aux = master->setVariableBounds(n, -OPT_INFINITY, zuaux); //do not pass zl here because it damages branch-and-bound performance
                        
                        if( quadAppStrategy == MRQ_QAMS_NO_QUAD_APP )
                            master->setUpperObjCut(zu);
                        
                        #if MRQ_DEBUG_MODE
                            if(aux != 0)
                            {
                                if( in_print_level > 0 )
                                    fprintf(stderr, "Error at setting auxiliary variable bound!\n");
                            }
                        #endif */
                    }
                }
            }
        }
        
        
        if(in_print_level > 1)
            printf("%-5ld  %+-14e  %+-14e  %+-14e\n", iter, zl, zu, zu-zl);
        
        if( checkTerminationCriterions(thnumber, zl, zu, iter, timeStart, clockStart, out_return_code) )
        {
            break;
        }
        
        
        if(binaryCut)
        {
            //aux = milp->addBinaryCut( milp2->sol );
            
            aux = masterMilp.addBinaryCut(nI, intVars, master->sol, auxVars);
            binaryCut = false;
        }
        else
        {
            int r;
            
            if(in_store_history_solutions)
                out_sol_hist.addSolution(n, iter, MRQ_getTime() - timeStart, clock() - clockStart, auxPSol, (objCalculated ? objTemp : INFINITY) );
            
            aux = masterMilp. addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, true, auxPSol, setQuadsInMaster, in_constr_linearisation_strategy, auxConstrEval, NULL, NULL, NULL, (constrCalculated ? auxPConstr : auxConstr), constrCalculated, masterSolConstr, newCalclc, newCalcuc );
            
            
            if( linearizeObj ) //nonlinear objective
            {
                
                if( in_obj_linearisation_strategy != MRQ_OLS_NON_OBJ_CUT_POINTS || !feasProbSol )
                {
                    //we always linearize the last point in the master problem
                    aux += masterMilp.addLinearizedObjFunction( true, auxPSol, setQuadsInMaster, auxCols, auxVars, (objCalculated ? &objTemp : NULL) );
                    
                    
                    
                    if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS )
                    {
                        int r;
                        int mmaster;
                        
                        //now, we check if we can discard the linearization on previous point...
                        
                        //auxVars has the gradient of obj linearization, i.e., the coefiicents of obj linearizations constraint and the rhs at the end...
                        r = laps.updateObjLinearizationsByNonObjCutPointsByNewPoint( *master, zu, auxVars );
                        
                        if( r != 0 )
                        {
                            if( in_print_level > 0 )
                                MRQ_PRINTERRORNUMBER(r);
                            
                            out_return_code = MRQ_MILP_SOLVER_ERROR;
                            goto termination;
                        }
                        
                        
                        r = laps.addPoint(n, auxPSol);
                        if( r != 0 )
                        {
                            if( in_print_level > 0 )
                                MRQ_PRINTERRORNUMBER(r);
                            
                            out_return_code = MRQ_MEMORY_ERROR;
                            goto termination;
                        }
                        
                        r = master->getNumberOfConstraints( mmaster);
                        #if MRQ_DEBUG_MODE 
                            if(r != 0)
                            {
                                if(in_print_level > 0)
                                    MRQ_PRINTERRORNUMBER(r);
                                
                                out_return_code = MRQ_MILP_SOLVER_ERROR;
                                goto termination;
                            }
                        #endif
                        
                        laps.indMaster[ laps.npoints-1 ] = mmaster-1;
                        
                        
                        //milp2->generateModelFile("master_readd_cortes.lp");
                        //std::cout << "Readicionei os cortes!" << std::endl;
                        //getchar();
                    }
                }
                
                
                if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS && feasProbSol)
                    out_number_of_obj_linears_saved_by_infeasibility++;
                
            }
            
            
            if( quadAppStrategy == MRQ_QAMS_ON_LAST_POINT ) // && auxPSol == nlp->sol )
            {
                OPT_QPSolver *qmaster = (OPT_QPSolver*) master;
                
                //if auxPSol == nlp->sol, we have a feasible solution gotten by feasible nlp problem.
                
                //for(int i = 0; i < m; i++)
                    //auxConstr[i] = -nlp->dualSol[i];
                
                r = intHess->evalCompleteHessian( thnumber, true, nlp->sol, 1.0, dualSol );
                
                if( r != 0 )
                {
                    if( in_print_level > 0 )
                        std::cerr << MRQ_PREPRINT "Callback function error " << aux << MRQ_GETFILELINE << "\n";
                    
                    out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
                    goto termination;
                }
                
                r = qmaster->setObjQuadMatrix( intHess->nzs, intHess->rows, intHess->cols, intHess->values ); //it is missing the therm -x⁰ * H * x⁰, but that is a constant and it is not really necessary
                
                if( r != 0 )
                {
                    if( in_print_level > 0 )
                        MRQ_PRINTERRORNUMBER(r);
                    
                    out_return_code = MRQ_MILP_SOLVER_ERROR;
                    goto termination;
                }
            }
        }
        
        
        if( aux )
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORMSG("Error to add cutting to MILP\n");
            
            out_return_code = MRQ_MILP_SOLVER_ERROR;
            break;
        }
        
        objCalculated = false;
        constrCalculated = false;
        feasProbSol = false;
        //getchar();
    }
    
    
    
    if(in_print_level > 1)
    {
        if(out_return_code == MRQ_OPTIMAL_SOLUTION)
        {
            std::cout << MRQ_PREPRINT "Optimal solution found. Objective: " << out_best_obj;
        }
        else if(out_return_code == MRQ_INFEASIBLE_PROBLEM)
        {
            std::cout << MRQ_PREPRINT "Infeasible problem. ";
        }
        else if(out_return_code == MRQ_MAX_ITERATIONS_STOP)
        {
            std::cout << MRQ_PREPRINT "Max number of iterations reached! ";
        }
        else if(out_return_code == MRQ_MAX_TIME_STOP)
        {
            std::cout << MRQ_PREPRINT "Max time reached! ";
        }
    }
    
    
    
    
termination:
    
    if(plc)				free(plc);
    if(flx)				free(flx);
    if(intVars)			free(intVars);
    if(auxCols)			free(auxCols);
    if(auxVars)			free(auxVars);
    if(auxConstr)		free(auxConstr);
    if(auxConstrEval)	free(auxConstrEval);
    if(dualSol)			free(dualSol);
    //if(auxConstrEval2)	free(auxConstrEval2);
    //if(master)			delete master;
    if(nlp)				delete nlp;
    //if(nlpFeas)			delete nlpFeas;
    if(intHess)			delete intHess;
    
    if(masterSolConstr)	free(masterSolConstr);
    
    if(newCalclc)		free(newCalclc);
    if(constrBoxCalcNLP)	delete constrBoxCalcNLP;
    
    
    out_feasible_solution = out_best_obj < MRQ_INFINITY;
    out_number_of_iterations = iter;
    out_algorithm = muriqui::MRQ_OA_ALG;
    out_lower_bound = zl;
    out_upper_bound = zu;
    
    out_number_of_obj_linears_saved_by_zu = laps.nObjLinRem;
    
    algorithmFinalization(1, prob, lx, ux);
    
    out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
    out_clock_time = MRQ_getTime() - timeStart;
    
    
    if(in_print_level > 1)
        std::cout << " cpu time: " << out_cpu_time << "\n";
    
    //std::cout << "in_measue_nlp_time: " << in_measure_nlp_time <<  " out_clock_time_of_nlp_solving: " << out_clock_time_of_nlp_solving << " out_cpu_time_of_nlp_solving: " << out_cpu_time_of_nlp_solving << "\n";
    
    return out_return_code;
}











