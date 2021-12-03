/*That file contains a implementation of
* Bonmin hybrid algorithm
*
* References:
*
* Bonami et al, An algorithm framework for convex mixed integer nonlinear
* programs. Discrete Optimization 5 (2008), pages 186-204.
*
*
* Author: Wendel Melo
*
* Date: 7-Jan-2016
*/


#include <climits>
#include <cstdio>
#include <iostream>


#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_milpCallbacks.hpp"


using namespace optsolvers;
using namespace muriqui;




#if OPT_HAVE_CPLEX

int CPXPUBLIC muriqui::MRQ_bonminHybCplexCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
    MRQ_MILPCallbackData *data = (MRQ_MILPCallbackData*) cbhandle;
    MRQ_BonminHybrid *alg = NULL;
    
    unsigned int thnumber; // I am not sure about the correct type to thread number. Unfortunatelly, CPXgetcallbackinfo receives a void pointer.
    
    const unsigned int nparams = 4; //size of array params
    void* params[] = {&env, cbdata, &wherefrom, useraction_p};
    
    
    const int r = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &thnumber);
    if(r != 0)
    {
        if(data->alg->in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        
        return MRQ_MILP_SOLVER_ERROR;
    }
    
    
    alg = (MRQ_BonminHybrid *) data[thnumber].alg;
    
    //std::cout << "thnumber: " << thnumber << "\n";
    
    data[thnumber].callbackSolver->initializeSolverData(nparams, params);
    
    return alg->solverCallbackCutToBonminHyb( data[thnumber]);
}
#endif









MRQ_BonminHybrid::MRQ_BonminHybrid():MRQ_LPNLPBBOuterApp()
{
    resetParameters();
    resetOutput();
    out_algorithm = MRQ_BONMIN_HYBRID_ALG;
}


MRQ_BonminHybrid::~MRQ_BonminHybrid()
{
}


int MRQ_BonminHybrid::checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx, const double *ux)
{
    if(in_milp_solver != MRQ_CPLEX)
    {
        char solverName[30];
        
        MRQ_enumToStr( in_milp_solver, solverName  );
        
        std::cerr << MRQ_PREPRINT "We are so sorry. Algorithm " <<  getAlgorithmName() << "(" << out_algorithm << ")" << " does not work with milp solver " << solverName  << "(" << in_milp_solver <<")  \n";
        
        return MRQ_BAD_PARAMETER_VALUES;
    }
    
    return MRQ_LPNLPBBOuterApp::checkAlgorithmRequirements(prob, lx, ux);
}


void MRQ_BonminHybrid::printParameters(std::ostream &out) const
{
    MRQ_LPNLPBBOuterApp::printParameters();
    out << "\n"
    
    MRQ_STRFFATT(in_out_app_max_iterations) << "\n"
    MRQ_STRFFATT(in_nlp_relaxation_solving_frequence) << "\n"
    MRQ_STRFFATT(in_outer_app_max_time) << "\n";
    
    out << "\ninside Outer Approximation parameters:\n\n";
    
    in_outer_app.printParameters(out);
}


void MRQ_BonminHybrid::resetParameters()
{
    in_out_app_max_iterations = ULONG_MAX;
    in_nlp_relaxation_solving_frequence= 10;
    
    in_outer_app_max_time =INFINITY;
    in_outer_app_max_cpu_time = 30;
}


int MRQ_BonminHybrid::setDoubleParameter(const char *name, const double value)
{
    int ret = MRQ_LPNLPBBOuterApp::setDoubleParameter(name, value);
    
    if(ret == 0)
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt<double>( MRQ_STRATT(in_outer_app_max_time), name, value ) == 0 );
    else if( MRQ_setAtt<double>( MRQ_STRATT(in_outer_app_max_cpu_time), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    return ret;
}


int MRQ_BonminHybrid::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_LPNLPBBOuterApp::setIntegerParameter(name, value);
    
    if(ret == 0)
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt<long unsigned int>( MRQ_STRATT(in_out_app_max_iterations), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_nlp_relaxation_solving_frequence), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    return ret;
}



int MRQ_BonminHybrid::solverCallbackCutToBonminHyb( MRQ_MILPCallbackData &data)
{
    MRQ_MINLPProb &prob = *data.prob;
    
    const int n = prob.n;
    const int thnumber = data.thnumber;
    
    const bool incQuadsInMaster = data.setQuadsInMaster;
    
    const int nI = data.nI;
    const int *intVars = data.intVars;
    const int *indices = data.indices;
    
    const bool *auxConstrEval = data.constrEval;
    bool *auxConstrEval2 = data.auxConstrEval2;
    
    double *auxVars = data.auxVars;
    double *auxVars2 = data.auxVars2;
    
    double NLPCpuTime, NLPClockTime;
    double *pNLPCpuTime = in_measure_nlp_time ? &NLPCpuTime : NULL;
    double *pNLPClockTime = in_measure_nlp_time ? &NLPClockTime : NULL;
    
    double *plc = data.plc, *puc = data.puc;
    
    MRQ_NLPSolver *nlp = (MRQ_NLPSolver *) data.nlp;
    MRQ_GradientsEvaluation &gradEval = data.gradEval;
    MRQ_MILPSolverCallbackInterface *callbackSolver = data.callbackSolver;
    
    int r, retCode;
    long int myiter;
    double &totalNLPCpuTime = data.out_cpu_time_of_nlp_solving;
    double &totalNLPClockTime = data.out_clock_time_of_nlp_solving;
    long unsigned int &nlpProbsSolved = data.out_number_of_nlp_probs_solved;
    
    //callbackSolver->initializeSolverData(params);
    
    
    r = callbackSolver->getNumberOfIterations(myiter);
    if(r != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        
        retCode = MRQ_MILP_SOLVER_ERROR;
        goto termination;
    }
    
    
    
    //std::cout << "myiter: " << myiter << "\n";
    
    if(myiter%in_nlp_relaxation_solving_frequence != 0 || data.lastiter == myiter)
    {
        retCode = 0;
        goto termination;
    }
    
    data.lastiter = myiter;
    
    
    
    r = callbackSolver->getVarsLowerBoundsOnNode(n, auxVars);
    r += callbackSolver->getVarsUpperBoundsOnNode(n, auxVars2);
    if(r != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        
        retCode = MRQ_MILP_SOLVER_ERROR;
        goto termination;
    }
    
    
    //I am not sure if cplex change continuous variable bounds. So, we only fix integer variables to do not damage the lazzy constraint procces...
    for(int i = 0; i < nI; i++)
    {
        const int ind = intVars[i];
        r += nlp->setVariableBounds(ind, auxVars[ind], auxVars2[ind]);
    }
    
    
    r += callbackSolver->getNodeSolution(n, auxVars);
    
    r += nlp->setInitialSolution(auxVars, NULL, NULL);
    
    if(r != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        
        retCode = MRQ_MILP_SOLVER_ERROR;
        goto termination;
    }
    
    
    nlp->solveAndGetTime(pNLPCpuTime, pNLPClockTime, true);
    
    //std::cout << "resolvi nlp. myiter: " << myiter << " retCode: " << nlp->retCode << " obj: " << nlp->objValue << "\n";
    
    if(in_measure_nlp_time)
    {
        totalNLPCpuTime += *pNLPCpuTime;
        totalNLPClockTime += *pNLPClockTime;
    }
    
    nlpProbsSolved++;
    
    if( nlp->retCode == OPT_OPTIMAL_SOLUTION )
    {
        const double objsol = nlp->objValue;
        const double* const psol = nlp->sol;
        const double* const pconstr = nlp->constr;
        
        //cplex uses the same function to add global cuts and lazzy contsraints...
        
        r = addLazyConstraintsLinearizationOnSolution( thnumber, callbackSolver, prob, gradEval, incQuadsInMaster, data.linearizeObj, auxConstrEval, psol, &objsol, pconstr, pconstr, indices, plc, puc, auxVars, auxConstrEval2);
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
            
            retCode = r;
            goto termination;
        }
        
        
        if(in_store_history_solutions)
        {
            SEMAPH_history.lock(nthreads);
            {
                out_sol_hist.addSolution(n, myiter, MRQ_getTime() - data.timeStart, clock() - data.clockStart, psol, objsol);
            }
            SEMAPH_history.unlock(nthreads);
        }
        
    }
    
    
    
    retCode = 0;
    
termination:
    
    return retCode;
}





int MRQ_BonminHybrid::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    const double timeStart = MRQ_getTime();
    const clock_t clockStart = clock();
    
    const int n = prob.n;
    //decltype(nPoints) original_nPoints = nPoints;
    
    int ret;
    long unsigned int oaNLPProbsSolved = 0;
    long unsigned int oaMILPIters = 0;
    double oaNLPCpuTime = 0.0, oaNLPClockTime = 0.0;
    
    
    resetOutput();
    
    //iterCount = 0;
    
    if(in_measure_nlp_time) //since we do not call algorithmInitialization before outer app, i thin it is better initilizae nlp time counters by myself
    {
        out_cpu_time_of_nlp_solving = 0.0;
        out_clock_time_of_nlp_solving = 0.0;
    }
    
    
    //theorically, we should call algorithm initialization before computations because this function will initialize the user evaluation callbaks. But, since we are not performign any problem function evaluation here, we do not call this function and we let our algorithms call them...
    
    if( in_outer_app_max_time > 0.0 && in_outer_app_max_cpu_time > 0.0 && in_out_app_max_iterations > 0 )
    {
        in_outer_app.deletePointsToLinearisation();
        
        in_outer_app.nPoints = nPoints;
        in_outer_app.points = points;
        
        ////we overwrite  some options. Unfortunatelly we need do it...
        in_outer_app.copyParametersFrom(*this);
        
        in_outer_app.in_store_history_solutions = true;
        
        in_outer_app.in_max_iterations = in_out_app_max_iterations;
        
        in_outer_app.in_number_of_threads = in_number_of_threads;
        
        in_outer_app.in_max_time = MRQ_min( in_max_time, in_outer_app_max_time );
        
        in_outer_app.in_max_cpu_time = MRQ_min( in_max_cpu_time, in_outer_app_max_cpu_time );
        
        
        
        in_outer_app.run(prob, milpSolverParams, nlpSolverParams);
        
        in_outer_app.nPoints = 0; //we need do it to avoid oa free our points 
        
        out_obj_opt_at_continuous_relax = in_outer_app.out_obj_opt_at_continuous_relax;
        oaNLPProbsSolved += in_outer_app.out_number_of_nlp_probs_solved;
        oaMILPIters += in_outer_app.out_number_of_milp_solver_iters;
        
        if(in_store_history_solutions)
        {
            const unsigned int nsols = out_sol_hist.getnsols();
            
            for(unsigned int i = 0; i < nsols; i++)
            {
                MRQ_HistorySolution *hsol = out_sol_hist.getHistSolPointer(i);
                
                out_sol_hist.addSolution(n, 0, hsol->gettime(), hsol->getcputime(), hsol->sol, hsol->getobjvalue());
            }
        }
        
        if( in_measure_nlp_time)
        {
            oaNLPCpuTime += in_outer_app.out_cpu_time_of_nlp_solving;
            oaNLPClockTime += in_outer_app.out_clock_time_of_nlp_solving;
        }
        
        if( in_outer_app.out_return_code == MRQ_OPTIMAL_SOLUTION || in_outer_app.out_return_code == MRQ_INFEASIBLE_PROBLEM || in_outer_app.out_return_code == MRQ_UNBOUNDED_PROBLEM )
        {
            int r = allocateBestSol(n);
            if(r != 0)
            {
                if(in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(r);
                
                out_return_code = muriqui::MRQ_MEMORY_ERROR;
                goto termination;
            }
            
            out_return_code = in_outer_app.out_return_code;
            out_best_obj = in_outer_app.out_best_obj;
            out_lower_bound = in_outer_app.out_lower_bound;
            
            if(out_return_code != MRQ_INFEASIBLE_PROBLEM)
                MRQ_copyArray(n, in_outer_app.out_best_sol, out_best_sol);
            
            goto termination;
        }
        
    }
    
    
    ret = addPointsToLinearisation(n, in_outer_app.out_sol_hist );
    
    if(ret != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(ret);
    }
    
    in_outer_app.out_sol_hist.desallocate();
    
    
    if(in_max_cpu_time < INFINITY)
    {
        if( MRQ_calcCPUTtime(clockStart, clock()) >= in_max_cpu_time )
        {
            out_return_code = MRQ_MAX_TIME_STOP;
            goto termination;
        }
    }
    
    
    if(in_max_time < INFINITY)
    {
        if( MRQ_getTime() - timeStart >= in_max_time )
        {
            out_return_code = MRQ_MAX_TIME_STOP;
            goto termination;
        }
    }
    
    
    MRQ_LPNLPBBOuterApp::run(prob, milpSolverParams, nlpSolverParams);
    
    
    
    
termination:
    
    //since we did not call algorithm initialization, we do not cll algogorithmFinalization
    
    in_outer_app.out_sol_hist.desallocate();
    
    out_feasible_solution = out_best_obj < MRQ_INFINITY;
    
    out_cpu_time_of_nlp_solving += oaNLPCpuTime;
    out_clock_time_of_nlp_solving += oaNLPClockTime;
    out_number_of_nlp_probs_solved += oaNLPProbsSolved;
    out_number_of_milp_solver_iters += oaMILPIters;
    
    out_cpu_time = MRQ_calcCPUTtime(clockStart, clock());
    out_clock_time = MRQ_getTime() - timeStart;
    
    if(in_print_level > 1)
        printf("cpu time: %f\n", out_cpu_time);
    
    return out_return_code;
}









