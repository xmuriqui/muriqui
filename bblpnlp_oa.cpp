/*That file contains a simple implementation of
* LP/NLP BB outter approximation based algorithm
*
* References:
*
* Bonami et al, An algorithm framework for convex mixed integer nonlinear
* programs. Discrete Optimization 5 (2008), pages 186-204.
*
* Quesada & Grossmann, An LP/NLP based branch and bound algorithm for convex MINLP optimization problems, Comput. Optim. Appl. 18 (1992),
* pages 937-947.
*
* Author: Wendel Melo
*
* Date: 18-Sept-2016
* 
* 
* Reference to set cplex callback before solve: admipex1.c
* 
*/


#include <cstdlib>
#include <cmath>
#include <climits>

#include <new>
#include <thread>
#include <typeinfo>   // operator typeid

#include "BBL_tools.hpp"

#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_milpCallbacks.hpp"


using namespace optsolvers;
using namespace muriqui;




MRQ_LPNLPBBOuterApp::MRQ_LPNLPBBOuterApp():MRQ_LinearApproxAlgorithm()
{
    resetParameters();
    //resetOutput();
    out_algorithm = MRQ_LP_NLP_BB_OA_BASED_ALG;
}


MRQ_LPNLPBBOuterApp::~MRQ_LPNLPBBOuterApp()
{
}


void MRQ_LPNLPBBOuterApp::printParameters(std::ostream &out) const
{
    MRQ_LinearApproxAlgorithm::printParameters(out);
    out << "\n"
    
    MRQ_STRFFATT(in_use_first_nlp_relaxation) << "\n"
    MRQ_STRFFATT(in_binarie_cut_when_nlp_infeasible) << "\n"
    MRQ_STRFFATT(in_linearize_obj_in_nl_feas_solutions) << "\n"
    ;
}


void MRQ_LPNLPBBOuterApp::resetParameters()
{
    MRQ_LinearApproxAlgorithm::resetParameters();
    
    in_binarie_cut_when_nlp_infeasible = false;
    in_use_first_nlp_relaxation = true;
    in_linearize_obj_in_nl_feas_solutions = true;
    in_number_of_threads = 1;
}


int MRQ_LPNLPBBOuterApp::checkAlgorithmRequirements( MRQ_MINLPProb &prob, const double *lx, const double *ux)
{
    if(!MRQ_isLazyConstraintsAvaliable(in_milp_solver))
    {
        //char solverName[30];
        //MRQ_enumToStr( in_milp_solver, solverName );
        
        std::string solverName = optsolvers::OPT_getSolverName(in_milp_solver);
        
        std::cerr << MRQ_PREPRINT "We are so sorry. Algorithm " <<  getAlgorithmName() << " (" << out_algorithm << ") does not work with milp solver " << solverName  << " (" << in_milp_solver <<")  \n";
        
        return MRQ_BAD_PARAMETER_VALUES;
    }
    
    return MRQ_LinearApproxAlgorithm::checkAlgorithmRequirements(prob, lx, ux);
}



int MRQ_LPNLPBBOuterApp::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_LinearApproxAlgorithm::setIntegerParameter(name, value);
    
    if(ret == 0)
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_first_nlp_relaxation), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_binarie_cut_when_nlp_infeasible), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_linearize_obj_in_nl_feas_solutions), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    return ret;
}




int MRQ_LPNLPBBOuterApp:: solverCallbackLazyConstraints( MRQ_MILPCallbackData &data)
{
    MRQ_MINLPProb &prob = *data.prob;
    
    const int n = prob.n;
    const int thnumber = data.thnumber;
    const int nI = data.nI;
    const int *intVars = data.intVars;
    
    long int iter = -1;
    
    bool binaryCut = false;
    bool feasible = false;
    //bool constr_added = false; //maybe we do not need this
    const bool incQuadsInMaster = data.setQuadsInMaster;
    int r, retCode = 0;
    double objsol = NAN;
    
    const int *indices = data.indices;
    const bool *auxConstrEval = data.constrEval;
    bool *auxConstrEval2 = data.auxConstrEval2;
    double *auxVars = data.auxVars;
    double *auxConstr = data.auxConstr;
    double *masterSol = data.auxVars2;
    
    double NLPCpuTime, NLPClockTime;
    double *pNLPCpuTime = in_measure_nlp_time ? &NLPCpuTime : NULL;
    double *pNLPClockTime = in_measure_nlp_time ? &NLPClockTime : NULL;
    
    double *plc = data.plc, *puc = data.puc;
    
    
    double *psol = NULL, *pconstr = NULL;
    
    
    MRQ_NLPSolver *nlp = (MRQ_NLPSolver *) data.nlp;
    MRQ_NLPSolver *nlpFeas = (MRQ_NLPSolver*) data.nlpFeas->solver;
    MRQ_GradientsEvaluation &gradEval = data.gradEval;
    MRQ_MILPSolverCallbackInterface *callbackSolver = data.callbackSolver;
    
    double &totalNLPCpuTime = data.out_cpu_time_of_nlp_solving;
    double &totalNLPClockTime = data.out_clock_time_of_nlp_solving;
    long unsigned int &nlpProbsSolved = data.out_number_of_nlp_probs_solved;
    
    
    //*useraction_p = CPX_CALLBACK_DEFAULT;
    //callbackSolver->initializeSolverData(params);
    
    
    
    
    /*double obj;
    
    r = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &obj);
    if(r != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        
        retCode = MRQ_MILP_SOLVER_ERROR;
        goto termination;
    }*/
    
    
    
    
    
    
    
    //r = CPXgetcallbacknodex(env, cbdata, wherefrom, masterSol, 0, n); //solution is a n+1 array
    r = callbackSolver->getNodeSolution(n+1, masterSol);
    if(r != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        
        retCode = MRQ_MILP_SOLVER_ERROR;
        goto termination;
    }
    
    
    r = MRQ_fixIntVarsOnSolByList(nI, intVars, masterSol, *nlp);
    if(r != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        
        retCode = MRQ_MILP_SOLVER_ERROR;
        goto termination;
    }
    
    
    nlp->setInitialSolution(masterSol, NULL, NULL);
    
    nlp->solveAndGetTime(pNLPCpuTime, pNLPClockTime, false);
    
    //std::cout << "nlpfix - retCode: " << nlp->retCode << " obj: " << nlp->objValue << " feas: " << nlp->feasSol << "\n";
    
    if(in_measure_nlp_time)
    {
        totalNLPCpuTime += *pNLPCpuTime;
        totalNLPClockTime += *pNLPClockTime;
    }
    
    nlpProbsSolved++;
    
    if(nlp->retCode == OPT_OPTIMAL_SOLUTION || (nlp->feasSol))       //&& nlp->objValue < zu) )
    {
        psol = nlp->sol;
        pconstr = nlp->constr;
        objsol = nlp->objValue;
        feasible = true;
    }
    else if(nlp->retCode == OPT_INFEASIBLE_PROBLEM)
    {
        bool binProblem = data.binProblem;
        
        if( in_print_level > 4 )
            std::cout << MRQ_PREPRINT "NLP problem is infeasible.\n";
        
        if(in_binarie_cut_when_nlp_infeasible && binProblem)
        {
            binaryCut = true;
        }
        else
        {
            MRQ_fixIntVarsOnSolByList(nI, intVars, masterSol, *nlpFeas);
            
            {//setting initial solution to nlp feasible problem from master solution
                int nfeas;
                double *xInitNLPFeas = nlpFeas->sol;
                
                r = nlpFeas->getNumberOfVars(nfeas);
                MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_NLP_SOLVER_ERROR, termination);
                
                MRQ_copyArray(n, masterSol, xInitNLPFeas);
                MRQ_setAllArray(nfeas - n, &xInitNLPFeas[n], 0.0 );
                
                nlpFeas->setInitialSolution(xInitNLPFeas, NULL, NULL);
            }
            
            int r = nlpFeas->solveAndGetTime(pNLPCpuTime, pNLPClockTime, false);
            
            //std::cout << "nlpfeas - code: " << nlpFeas->retCode << " obj: " << nlpFeas->objValue << " feasSol: " << nlpFeas->feasSol << "\n";
            
            if(in_measure_nlp_time)
            {
                totalNLPCpuTime += *pNLPCpuTime;
                totalNLPClockTime += *pNLPClockTime;
            }
            
            nlpProbsSolved++;
            
            if(r == OPT_OPTIMAL_SOLUTION)
            {
                psol = nlpFeas->sol;
            }
            else
            {
                if(in_print_level > 4)
                    std::cerr << MRQ_PREPRINT "Failure to solve NLP feasibility relaxation. \n";
                
                if(binProblem)
                {
                    binaryCut = true;
                }
                else
                {
                    //our unique choice is add a cut on the milp solution, like ecp
                    psol = masterSol;
                }
            }
        }
        
    }
    else
    {
        //check if milp solution is feasible
        psol = masterSol;
    }
    
    
    
    if(binaryCut)
    {
        //here, we add the binary cut. Note, if we add binary cut, we do not linearise objective function.
        double rhs;
        int r;
        
        MRQ_calculateBinCut(nI, intVars, masterSol, auxVars, rhs);
        
        //SEMAPH_addLazy.lock(nthreads);
        {
            r = callbackSolver->setLazyConstraint(nI, intVars, auxVars, -INFINITY, rhs);
        }
        //SEMAPH_addLazy.unlock(nthreads);
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
            
            retCode = MRQ_MILP_SOLVER_ERROR;
            goto termination;
        }
        
        //std::cout << "Adicionei corte binario! nI: " << nI << " rhs: " << rhs << "\n";
        /*for(int i = 0; i < nz; i++)
            std::cout << "i: " << i << " col: " << auxCols[i] << " var: " << auxVars[i] << " \t";
        std::cout << "\n";*/
        
        /*std::cout << "integer sol: \n";
        for(int i = 0; i < nI; i++)
            std::cout << " x_" << intVars[i] << ": " << masterSol[intVars[i]];
        std::cout << "\n"; */
        
        
        //constr_added = true;
    }
    else
    {
        bool newx = true;
        int r;
        
        
        if(pconstr == NULL)
        {
            int r = prob.isFeasibleToConstraints(thnumber, psol, newx, auxConstrEval, in_absolute_feasibility_tol, in_relative_feasibility_tol, feasible, auxConstr);
            
            if(r != 0)
            {
                if(in_print_level > 0)
                    MRQ_PRINTCALLBACKERRORNUMBER(r);
                
                retCode = MRQ_CALLBACK_FUNCTION_ERROR;
                goto termination;
            }
            
            if(prob.hasNlConstrs)
                newx = false;
            
            //here, we only check by feasibility considering nonlinear constraints. We can have a situation where we have a solution from nlpFeas only infeasible by linear constraints. So, we put this if here..
            if(psol == nlpFeas->sol)
                feasible = false;
            
            pconstr = auxConstr;
        }
        
        
        if(std::isnan(objsol))
        {
            int r = prob.objEval(thnumber, newx, psol, objsol);
            if(r != 0)
            {
                if(in_print_level > 0)
                    MRQ_PRINTCALLBACKERRORNUMBER(r);
                
                retCode = MRQ_CALLBACK_FUNCTION_ERROR;
                goto termination;
            }
            
            if(prob.hasNlConstrs)
                newx = false;
        }
        
        if(feasible)
        {
            if(objsol < zu)
            {
                //std::cout << "Indo atualizar sol. obj: " << objsol << " zu: " << zu << "\n";
                
                callbackSolver->getNumberOfIterations(iter);
                
                SEMAPH_updtSol.lock(nthreads);
                {
                    tryUpdateBestSolution(thnumber, n, psol, objsol, iter, data.clockStart, data.timeStart, false);
                }
                SEMAPH_updtSol.unlock(nthreads);
                
                if( zu <= zl )
                {//we could, but, by now, w do not considere optimal tolerances here
                    if(in_print_level > 2)
                        MRQ_PRINTMSG("A feasible solution better than lower bound was found. stopping");
                    
                    data.out_sol_lower_than_zl = true;
                    retCode = MRQ_LAZY_MILP_BB_SOLUTION_LOWER_THAN_ZL;
                    goto termination;
                }
            }
        }
        
        
        #if 0
        
        //if this solution is from nlpfeas, we save the obj linearization since this solution is infeasible
        if(data.linearizeObj && (psol != nlpFeas->sol || in_linearize_obj_in_nl_feas_solutions) )
        {
            double rhs;
            
            int r = MRQ_calculateObjLinearizedConstraint(prob, thnumber, newx, psol, incQuadsInMaster, &objsol, NULL, auxVars, rhs);
            if(r != 0)
            {
                if(in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(r);
                
                retCode = MRQ_MILP_SOLVER_ERROR;
                goto termination;
            }
            
            if(prob.hasNlObj)
                newx = false;
            
            SEMAPH_addLazy.lock(nthreads);
            {
                //r = CPXcutcallbackadd(env, cbdata, wherefrom, n+1, rhs, 'L', indices, auxVars, CPX_USECUT_PURGE);
                
                r = callbackSolver->setLazyConstraint(n+1, indices, auxVars, -INFINITY, rhs);
            }
            SEMAPH_addLazy.unlock(nthreads);
            if(r != 0)
            {
                if(in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(r);
                
                retCode = MRQ_MILP_SOLVER_ERROR;
                goto termination;
            }
            
            //std::cout << "linearizei funcao objetivo\n";
            
            constr_added = true;
        }
        
        
        //linearising constraints
        {
            const int m = prob.m;
            const bool *nlConstr = prob.nlConstr;
            MRQ_SparseMatrix *QC = prob.QC;
            double *grad = auxVars;
            
            //TODO: remove linearizations tarategy by active and infeas in master solution...
            MRQ_calculateConstraintsToBeLinearizedByStrategy(prob, in_constr_linearisation_strategy, in_eps_to_active_constr_to_linearisation, auxConstrEval, pconstr, pconstr, plc, puc, &out_number_of_constr_linears_saved, auxConstrEval2);
            
            
            if(prob.hasNlConstrs)
            {
                int r = gradEval.evaluateJacobian(newx, auxConstrEval2, psol);
                
                if(r != 0)
                {
                    #if MRQ_DEBUG_MODE
                        MRQ_PRINTERRORNUMBER(r);
                    #endif
                    
                    retCode = MRQ_CALLBACK_FUNCTION_ERROR;
                    goto termination;
                }
                
                if( prob.hasNlConstrs )
                    newx = false;
            }
            
            
            
            for(int i = 0; i < m; i++)
            {
                if(auxConstrEval2[i] == false)
                    continue;
                
                if( nlConstr[i] == false && ( QC[i].getNumberOfElements() == 0u || incQuadsInMaster ) )
                    continue;
                
                
                gradEval.constraintCompleteGradient(i, psol, grad);
                
                const double deltaRHS = MRQ_calculateDeltaRHSToConstraintLinearisation(pconstr[i], n, grad, psol);
                
                const double lc = plc[i] > -MIP_INFINITY ? plc[i] + deltaRHS : -INFINITY;
                const double uc = puc[i] < MIP_INFINITY ? puc[i] + deltaRHS : INFINITY;
                
                //const char sense = MRQ_cplexConstraintSense(plc[i], puc[i]);
                
                
                
                SEMAPH_addLazy.lock(nthreads);
                {
                    //r = CPXcutcallbackadd(env, cbdata, wherefrom, n, rhs, sense, indices, grad, CPX_USECUT_PURGE);
                    
                    r = callbackSolver->setLazyConstraint(n, indices, grad, lc, uc);
                }
                SEMAPH_addLazy.unlock(nthreads);
                
                if(r != 0)
                {
                    if(in_print_level > 0)
                        MRQ_PRINTERRORNUMBER(r);
                    
                    retCode = MRQ_MILP_SOLVER_ERROR;
                    goto termination;
                }
            
            }
            
        }
        
        #endif
        
        
        r = addLazyConstraintsLinearizationOnSolution( thnumber, callbackSolver, prob, gradEval, incQuadsInMaster, data.linearizeObj && (psol != nlpFeas->sol || in_linearize_obj_in_nl_feas_solutions), auxConstrEval, psol, &objsol, pconstr, pconstr, indices, plc, puc, auxVars, auxConstrEval2);
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
            
            retCode = r;
            goto termination;
        }
        
        
        
        if(in_store_history_solutions)
        {
            if(iter == -1)
                callbackSolver->getNumberOfIterations(iter);
            
            SEMAPH_history.lock(nthreads);
            {
                out_sol_hist.addSolution(n, iter, MRQ_getTime() - data.timeStart, clock() - data.clockStart, psol, objsol);
            }
            SEMAPH_history.unlock(nthreads);
        }
            
    }
    
    
    
    
    
    
    
    /* Tell CPLEX that cuts have been created 
    if(constr_added)
    {
        *useraction_p = CPX_CALLBACK_SET;
    }*/
    
termination:
    
    //if(retCode != 0)
        //*useraction_p = CPX_CALLBACK_FAIL;
    
    return retCode;
}



int MRQ_LPNLPBBOuterApp::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    const double timeStart = MRQ_getTime();
    const clock_t clockStart = clock();
    
    const bool binProblem = prob.isBinaryProblem();
    const int n = prob.n;
    const int m = prob.m;
    const int nI = prob.getNumberOfIntegerVars();
    const bool preproc = in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_obj_function;
    
    const bool setQuadsInMaster = in_set_quadratics_in_master_problem;
    const bool linearizeObj = prob.hasObjNLTerm() || (prob.Q.getNumberOfElements() > 0 && !setQuadsInMaster);
    const bool useFirstNlp = in_use_first_nlp_relaxation || nPoints == 0;
    
    bool updtConstrBounds;
    int ret;
    
    bool *constrEval = NULL, *auxConstrEval = NULL; 
    int *intVars = NULL;
    int *indices = NULL;
    
    
    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux; 
    double *plc = NULL, *puc = NULL; //puc must be initialized due to MRQ_setNLPRelaxProb
    
    MRQ_NLPSolver *nlp;
    //MRQ_NLPFeasProb nlpFeas;
    MRQ_Preprocessor preprocessor(&prob);
    MRQ_LAAPointsStoring laps(n);
    MRQ_GradientsEvaluation gradEval;
    
    OPT_LPSolver *master;
    MRQ_MasterMILPProb masterMilp;
    
    MRQ_MILPCallbackData *callbackData = NULL;
    MRQ_Mutex SEMAPH_addLazy;
    
    
    
    nthreads = in_number_of_threads > 0 ? in_number_of_threads : branchAndBound::BBL_getNumCores() ;
    
    nthreads_lazy = nthreads;
    if( in_milp_solver == MRQ_GUROBI )
        nthreads_lazy = 1; //gurobi apply multithreading to solve the problem, but only thread 0 is called to add lazy constraints.
    
    
    {
        auto ret = algorithmInitialization(nthreads, preproc, milpSolverParams, nlpSolverParams, prob, lx, ux, &preprocessor, &updtConstrBounds, &plc, &puc);
        
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
        MRQ_PRINTMSG("Starting LP/NLP Branch-And-Bound based on Outer Approximation\n\n");
    }
    
    
    
    
    /*if(in_milp_solver == MRQ_GUROBI)
        nthreads = 1; //I am sorry, but due to solver limitations, we cannot run this algorithm in more than 1 thread... */
    
    
    if(in_print_level > 3)
        printSubSolvers(true, true, false);
    
    MRQ_malloc(intVars, nI); //intVars = (int *) malloc(nI *sizeof(int));
    MRQ_malloc(constrEval, 2*m); //auxConstrEval = (bool *) malloc(2*m *sizeof(bool));
    MRQ_malloc(indices, n+1); //auxCols = (int *) malloc((n+1) *sizeof(int));
    callbackData = new (std::nothrow) MRQ_MILPCallbackData[nthreads];
    
    if(!intVars || !constrEval || !indices || !callbackData)
    {
        if(in_print_level > 0)
            MRQ_PRINTMEMERROR;
        
        out_return_code = MRQ_MEMORY_ERROR;
        goto termination;
    }
    
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
    
    
    for(unsigned int i = 0; i < nthreads; i++)
    {
        MRQ_NLPSolver* nlp;
        MRQ_NLPFeasProb *nlpFeas;
        
        int r = callbackData[i].allocateBase(thnumber + i, nthreads, &prob, in_milp_solver, false, &in_nlp_solver, &SEMAPH_addLazy);
        
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
            
            out_return_code = muriqui::MRQ_MEMORY_ERROR;
            goto termination;
        }
        
        
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
        
        
        /************* Setting NLP relaxation *****************/
        
        nlp = (MRQ_NLPSolver* &) callbackData[i].nlp;
        
        
        ret = MRQ_setNLPRelaxProb(prob, lx, ux, plc, puc, nlp, true, true, true, false, thnumber + i, in_set_special_nlp_solver_params, nlpSolverParams, in_number_of_threads, in_max_cpu_time, in_max_time, 0, 0);
        MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
        
        if( run_by_inside )
        {
            if( !std::isinf(insideSolverMaxTime) )
            {
                ret = nlp->setMaxTime(insideSolverMaxTime );
                MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
            }
        }
        
        /*ret = nlp->setnVariablesBounds(n, lx, ux);
        MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
        
        if( updtConstrBounds )
        {
            for(int i = 0; i < m; i++)
            {
                ret = nlp->setConstraintBounds(i, plc[i], puc[i]);
                MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
            }
        }*/
        
        if( nlp->getSolverCode() != optsolvers::OPT_MOSEK && in_use_initial_solution )
        {
            const int r = nlp->setInitialSolution(xInit, NULL, NULL);
            
            if(r != 0)
            {
                if( in_print_level > 0 )
                    std::cerr << MRQ_PREPRINT "Error " << r << " at setting initial solution" << MRQ_GETFILELINE << "\n";
            }
        }
        
        
        
        /***** Setting NLP feasibility problem *****/
        nlpFeas = callbackData[i].nlpFeas;
        
        r = nlpFeas->setProblem(in_nlp_solver, prob, lx, ux, nlpSolverParams, thnumber + i, in_set_special_nlp_solver_params, in_number_of_threads, in_max_cpu_time, in_max_time);
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
        
        
        if(run_by_inside)
        {
            if(!std::isinf(insideSolverMaxTime))
            {
                r = nlpFeas->solver->setMaxTime(insideSolverMaxTime);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
            }
        }
        
        if( updtConstrBounds )
        {
            for(int k = 0; k < m; k++)
            {
                r = nlp->setConstraintBounds(k, plc[k], puc[k]);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
            }
        }
        
        /***** end of NLP feasibility problem seting *****/
    }
    
    
    
    nlp = (MRQ_NLPSolver*) callbackData[0].nlp; //taking advantage object allocated for first thread
    
    /**** end of NLP relaxation seting *****/
    
    
    
    /**** Solving NLP relaxation  *****/
    if(useFirstNlp)
    {
        double NLPCpuTime, NLPClockTime;
        double *pNLPCpuTime = in_measure_nlp_time ? &NLPCpuTime : NULL;
        double *pNLPClockTime = in_measure_nlp_time ? &NLPClockTime : NULL;
        
        
        if(in_print_level > 2)
            MRQ_PRINTMSG("Solving NLP relaxation\n");
        
        const int retCode = nlp->solveAndGetTime(pNLPCpuTime, pNLPClockTime, false);
        
        if( in_measure_nlp_time )
        {
            out_cpu_time_of_nlp_solving += *pNLPCpuTime;
            out_clock_time_of_nlp_solving += *pNLPClockTime;
        }
        
        out_number_of_nlp_probs_solved++;
        
        if( retCode == OPT_OPTIMAL_SOLUTION )
        {
            out_obj_opt_at_continuous_relax = nlp->objValue;
            if(nlp->objValue > zl)
                zl = nlp->objValue;
            
            if(in_print_level > 1)
                std::cout << MRQ_PREPRINT "NLP relaxation solution: " << nlp->objValue << "\n";
            
            if(in_print_level > 5)
            {
                double *sol = nlp->sol;
                for(int i = 0; i < n; i++)
                    std::cout << MRQ_PREPRINT << "sol["<<i<<"]: " << sol[i] << "\n";
            }
            
            if(in_store_history_solutions)
                out_sol_hist.addSolution(n, 0, MRQ_getTime() - timeStart, MRQ_calcCPUTtime(clockStart), nlp->sol, nlp->objValue);
            
            if(zl > zu)
            {
                if(in_print_level > 0)
                    MRQ_PRINTMSG("Solution of NLP relaxation is greater than upper_bound ");

                out_return_code = MRQ_INFEASIBLE_PROBLEM;
                goto termination;
            }
            
            if( MRQ_isIntegerSol(nI, intVars, nlp->sol, in_integer_tol) )
            {
                if( in_print_level > 1 )
                    MRQ_PRINTMSG("An integer optimal solution was gotten as NLP relaxation solution\n");
                
                tryUpdateBestSolution(thnumber, n, nlp->sol, nlp->objValue, 0, clockStart, timeStart, false); //we already store this solution in the history if in_store_history_solutions is true
                
                out_return_code = MRQ_OPTIMAL_SOLUTION;
                
                goto termination;
            }
            
        }
        else if( nlp->retCode == OPT_INFEASIBLE_PROBLEM)
        {
            if( in_print_level > 2 )
                MRQ_PRINTMSG("Continuous relaxation is infeasible!\n");
            
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
    }
    
    /**** end of NLP relaxation solving *****/
    
    
    /***** Setting master problem *****/
    
    
    //if solver does not suport quadratics, setProblemBase will return an error
    
    //ret = masterMilp.setProblemBase(thnumber, prob, in_milp_solver, true, setQuadsInMaster, true, lx, ux, 1, milpSolverParams);
    ret = setMasterProblemBase(&masterMilp, thnumber, prob, in_milp_solver, true, setQuadsInMaster, true, lx, ux, 1, milpSolverParams, auxConstrEval, &laps, nthreads);
    MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    
    master = masterMilp.master;
    
    ret = master->setRelativeOptimalityTol( in_relative_convergence_tol );
    MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    
    
    if(useFirstNlp)
    {
        bool newx = true;
        
        if(!nlp->feasSol)
        {
            //so, we have to improvise a point to have the first point to linearize
            double *sol = nlp->sol;
            
            for(int i = 0; i < n; i++)
                sol[i] = MRQ_min(ux[i], MRQ_max(0.0, lx[i]) );
            
            int r = prob.objEval(thnumber, newx, nlp->sol, nlp->objValue);
            if(r != 0)
            {
                if( in_print_level > 0 )
                    MRQ_PRINTCALLBACKERRORNUMBER(r);
                
                out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
                goto termination;
            }
            
            if(prob.hasNlObj)
                newx = false;
            
            r = prob.constraintsEval(thnumber, newx, constrEval, nlp->sol, nlp->constr);
            if(r != 0)
            {
                if( in_print_level > 0 )
                    MRQ_PRINTCALLBACKERRORNUMBER(r);
                
                out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
                goto termination;
            }
            
            if(prob.hasNlConstrs)
                newx = false;
        }
        
        int ret = masterMilp.addLinearizedNLConstraintsByStrategy( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, newx, nlp->sol, setQuadsInMaster, in_constr_linearisation_strategy, constrEval, NULL, NULL, NULL, nlp->constr, true );
        
        if(ret != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(ret);
            
            out_return_code = muriqui::MRQ_MILP_SOLVER_ERROR;
            goto termination;
        }
        
        
        if(prob.hasNlConstrs)
            newx = false;
        
        
        if(linearizeObj)
        {
            double *auxVars = callbackData[0].auxVars;
            
            int ret = masterMilp.addLinearizedObjFunction( newx, nlp->sol, setQuadsInMaster, indices, auxVars, &(nlp->objValue));
            
            if(ret != 0)
            {
                if(in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(ret);
                
                out_return_code = MRQ_MILP_SOLVER_ERROR;
                goto termination;
            }
            
            
            if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS )
            {
                int mmaster;
                
                //auxVars has the coeficients of linearization to objective function and the RHS also
                
                ret = laps.updateObjLinearizationsByNonObjCutPointsByNewPoint( *master, zu, auxVars );
                if(ret != 0)
                {
                    if( in_print_level > 0 )
                        MRQ_PRINTERRORNUMBER(ret);
                    
                    out_return_code = muriqui::MRQ_MILP_SOLVER_ERROR;
                    goto termination;
                }
                
                
                ret = laps.addPoint(n, nlp->sol);
                if(ret != 0)
                {
                    if( in_print_level > 0 )
                        MRQ_PRINTERRORNUMBER(ret);
                    
                    out_return_code = muriqui::MRQ_MEMORY_ERROR;
                    goto termination;
                }
                
                
                ret = master->getNumberOfConstraints( mmaster);
                if(ret != 0)
                {
                    if( in_print_level > 0 )
                        MRQ_PRINTERRORNUMBER(ret);
                    
                    out_return_code = MRQ_MILP_SOLVER_ERROR;
                    goto termination;
                }
                
                laps.indMaster[laps.npoints-1] = mmaster;
            }
        }
        
        
    }
    
    /***** end of master problem setting *****/
    
    
    
    /***** setting callback data to milp solver *****/
    
    #pragma ivdep
    #pragma GCC ivdep
    for(int i = 0; i <= n; i++) //considering auxiliary variables also
        indices[i] = i;
    
    
    
    
    /***** end of callback data to milp solver setting *****/
    
    
    
    
    /***** setting up solver to adopt lazy constraints *****/
    //if the problem is linear, we do not need lazy constraint
    if( prob.getProblemType() != minlpproblem::MIP_PT_MILP )
    {
        ret = setLazyConstraintCallbackOnMilpSolver(master, callbackData);
        MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
    }
    
    
    if( in_user_callbacks )
    {
        if(in_call_before_solve_callback_in_milp_bb)
        {
            ret = setBeforeSolveCallbackOnMilpSolver(master, callbackData);
            MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
        
        if(in_call_branching_callback_in_milp_bb)
        {
            ret = setBranchingCallbackOnMilpSolver(master, callbackData);
            MRQ_IFERRORGOTOLABEL(ret, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
        }
    }
    
    
    /***** setting up solver to use callback to add cuts. We use that for bonmin hybrid algorithm  *****/
    
    if( typeid(*this) == typeid(MRQ_BonminHybrid) )
    {
        
        switch(master->getSolverCode())
        {
            case optsolvers::OPT_CPLEX:
            {
            #if OPT_HAVE_CPLEX
                optsolvers::OPT_Cplex *cplex = (OPT_Cplex*) master;
                
                CPXENVptr cplex_env = cplex->env;
                
                int r = CPXsetusercutcallbackfunc(cplex_env, MRQ_bonminHybCplexCallback, callbackData);
                MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MILP_SOLVER_ERROR, termination);
                
            #endif
                
                break;
            }
            default:
            {
                MRQ_PRINTERRORMSGP("Invalid solver code: ", master->getSolverCode());
                out_return_code = MRQ_NONIMPLEMENTED_ERROR;
                goto termination;
            }
        }
        
    }
    
    //std::cout <<  "cpu time: " << MRQ_calcCPUTtime(clockStart) << "\n";
    
    /***** end of master problem setting *****/
    
    
    /*master->generateModelFile("bboa.lp");
    std::cout << "Gerei bboa.lp\n";
    MRQ_getchar(); */
    
    
    ret = master->solve();
    
    /*std::cout << "master solving ret code: " << master->retCode << " obj: " << master->objValue << " orig ret code: " << master->origSolverRetCode << "\n";
    
    {
        double *masterSol = master->sol;
        
        std::cout << "integer sol: \n";
        for(int i = 0; i < nI; i++)
            std::cout << " x_" << intVars[i] << ": " << masterSol[intVars[i]];
        std::cout << "\n";
    }*/
    
    zl = master->getDualObjValue();
    
    if(ret == OPT_OPTIMAL_SOLUTION)
    {
        //zl = master->getObjValue();
        out_return_code = MRQ_OPTIMAL_SOLUTION;
    }
    else if(ret == OPT_MAX_TIME)
    {
        out_return_code = MRQ_MAX_TIME_STOP;
    }
    else if(ret == OPT_INFEASIBLE_PROBLEM)
    {
        out_return_code = MRQ_INFEASIBLE_PROBLEM;
    }
    else if(ret == OPT_CALLBACK_FUNCTION_ERROR || ret == MRQ_CALLBACK_FUNCTION_ERROR)
    {
        out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
    }
    else if(ret == OPT_MAX_ITERATIONS)
    {
        out_return_code = MRQ_MAX_ITERATIONS_STOP;
    }
    else if(ret == OPT_UNBOUNDED_PROBLEM)
    {
        out_return_code = MRQ_UNBOUNDED_PROBLEM;
    }
    else if(ret == MRQ_NLP_SOLVER_ERROR)
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
    
    
    for(unsigned int i = 0; i < nthreads; i++)
    {
        out_cpu_time_of_nlp_solving += callbackData[i].out_cpu_time_of_nlp_solving;
        out_clock_time_of_nlp_solving += callbackData[i].out_clock_time_of_nlp_solving;
        out_number_of_nlp_probs_solved += callbackData[i].out_number_of_nlp_probs_solved;
    }
    
    
    ret = master->getNumberOfIterations(out_number_of_iterations);
    if(ret != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(ret);
        out_number_of_iterations = ULONG_MAX;
    }
    
    
    
termination:
    
    if(plc)			free(plc);
    
    if(intVars)		free(intVars);
    if(constrEval)	free(constrEval);
    if(indices)		free(indices);
    
    if(callbackData)	delete[] callbackData;
    
    
    
    algorithmFinalization(nthreads, prob, lx, ux);
    
    out_number_of_threads = nthreads;
    out_number_of_milp_solver_iters = out_number_of_iterations;
    out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
    out_clock_time = MRQ_getTime() - timeStart;
    
    if(in_print_level > 1)
        std::cout << MRQ_PREPRINT "cpu time: " << out_cpu_time << "\n";
    
    
    return out_return_code;
}








