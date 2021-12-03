
#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cassert>



#include <new>

#include "MRQ_algClasses.hpp"
//#include "MRQ_nlpSolvers.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_tools.hpp"




using namespace optsolvers;
using namespace muriqui;

//using namespace std;



MRQ_FeasibilityPump::MRQ_FeasibilityPump():MRQ_Heuristic()
{
    resetParameters();
    resetOutput();
    out_algorithm = MRQ_FP_HEUR_ALG;
}




MRQ_FeasibilityPump::~MRQ_FeasibilityPump()
{ 
    
}





int MRQ_FeasibilityPump::checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx, const double *ux)
{
    return 0;
    //return  MRQ_isBinarieProblemAtRegion(prob, lx, ux) ? 0 : MRQ_ALG_NOT_APPLICABLE;
}


void MRQ_FeasibilityPump::printParameters(std::ostream &out) const
{
    MRQ_Heuristic::printParameters(out);
    out << "\n"
    
    MRQ_STRFFATT(in_set_linear_obj_term_on_bin_vars_at_nlp) << "\n"
    MRQ_STRFFATT(in_set_norm1_on_nlp) << "\n"
    MRQ_STRFFATT(in_last_iters_considered_to_cycles) << "\n"
    MRQ_STRFFATT(in_max_cycle_subiters) << "\n"
    MRQ_STRFFATT(in_lower_bound_to_pi) << "\n"
    MRQ_STRFFATT(in_upper_bound_to_pi) << "\n"
    ;
}


void MRQ_FeasibilityPump::resetParameters()
{
    MRQ_Heuristic::resetParameters();
    
    in_set_norm1_on_nlp = true;
    in_set_linear_obj_term_on_bin_vars_at_nlp = true;
    in_last_iters_considered_to_cycles = 5;
    in_max_cycle_subiters = 20;
    
    in_lower_bound_to_pi = -0.3;
    in_upper_bound_to_pi = 0.7;
    
    //ok, it is not a parameter, but i think it is a good idea set it here...
    //out_algorithm = MRQ_FP_HEUR_ALG;
}


int MRQ_FeasibilityPump::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_Heuristic::setIntegerParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_linear_obj_term_on_bin_vars_at_nlp), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_norm1_on_nlp), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_last_iters_considered_to_cycles), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_max_cycle_subiters), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


int MRQ_FeasibilityPump::setDoubleParameter(const char *name, const double value)
{
    int ret = MRQ_Heuristic::setDoubleParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt( MRQ_STRATT(in_lower_bound_to_pi), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_upper_bound_to_pi), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}





/*void inline MRQ_roundSol(const int ninds, const unsigned int* inds, double* sol)
{
    int i;
    
    for(i = 0; i < ninds; i++)
        sol[inds[i]] = round( sol[inds[i]] );
} */


void inline MRQ_roundSol(const int ninds, const int* inds, const double* solin, double *solout)
{
    int i;
    
    for(i = 0; i < ninds; i++)
        solout[inds[i]] = round( solin[inds[i]] );
}






/*inline void MRQ_convertBinSolToUIntArray(const int ninds, const unsigned int* inds, const double* sol, const int sizeConvArray, unsigned int *convArray)
{
    const int bytesBin = sizeof( unsigned int );
    int aux;
    int i, j, k;
    
    
    for(i = k = 0; i < sizeConvArray; i++)
    {
        convArray[i] = 0;
        
        aux = MRQ_min( bytesBin, ninds - k);
        
        for(j = 0; j < aux; j++, k++)
        {
            
            if( sol[ inds[k] ] > 0.5 ) //we consider 1.0
                convArray[i] |= 1 << j; //we can use | (bitwise or) operator because we are sure we do not have "go 1 (vai um)"
        }
    }
    
    
} */


bool inline MRQ_checkSolLastSols(const int nI, const int *intVars, const double *sol, const int nLastSols, int ** const lastIntSols)
{
    bool equal;
    int i, j, aux;
    
    for(i = 0; i < nLastSols; i++)
    {
        equal = true;
        
        for(j = 0; j < nI; j++)
        {
            aux = (int) sol[ intVars[j] ];
            
            if( lastIntSols[i][j] != aux )
            {
                equal = false;
                break;
            }
        }
        
        if(equal)
            return true;
    }
    
        
    return false;
}


void inline MRQ_copySolToIntSol(const int nI, const int *intVars, const double *sol, int *intSol)
{
    int i;
    
    for(i = 0; i < nI; i++)
        intSol[i] = sol[ intVars[i] ];
}

int inline MRQ_flipIntSol( MRQ_Random &random, const double lbpi, const double ubpi, const double *lx, const double *ux, const int nI, const int *intVars, const double *nlpSol, double *sol )
{
    long int j;
    unsigned int ind, k;
    double aux;
    
    unsigned int w = 0;
    
    //choosing a coordinate to guarantee the fliping
    
    do
    {
        k = intVars[ random.randInt(0, nI -1) ]; //only this line should be in this loop. The block below is just to test if there is some bug in our implementation
        
        
        if( w >= 20lu*nI ) //this block is just to test if there is some bug in our code i.e., the solution is already fixed
        {
            int i = 0;
            //probably we have a bug here. We run integer variables to check if the solution is already fixed
            for(i = 0; i < nI; i++)
            {
                k = intVars[i];
                if( ceil(lx[k]) != floor(ux[k]) )
                    break;
            }
            
            if(i == nI)
            {
                //if we reach this point, it is because this solution is already fixed. It is a bug, I think
                MRQ_PRINTERRORMSG("Fixed solution at MRQ_flipIntSol. Please send this message and instance to " MRQ_EMAIL);
                
                return MRQ_UNDEFINED_ERROR;
            }
        }
        
        w++;
    }while( ceil(lx[k]) == floor(ux[k]) );
    
    
    
    for(int i = 0; i < nI; i++)
    {
        ind = intVars[i];
        
        if( ceil(lx[ind]) == floor(ux[ind]) )
            continue;
        
        
        aux = MRQ_max( random.random(lbpi, ubpi), 0.0);
        
        
        if( MRQ_abs(nlpSol[ind] - sol[ind]) + aux > 0.5 || ind == k)
        {
            do
            {
                j = random.randInt( lx[ind], ux[ind] );
            }while( j == (long int) sol[ind] );
            
            sol[ind] = j;
        }
    }
    
    return 0;
}



void inline MRQ_getHighestGapBetweenSols(const double *lx, const double *ux, const int nI, int *intVars, const double *sol1, const double *sol2, unsigned int &indhighestGap, double &highestGap )
{
    int i;
    unsigned int j;
    double gap;
    
    
    highestGap = -1.0;
    
    for(i = 0; i < nI; i++)
    {
        j = intVars[i];
        
        gap = MRQ_abs( sol1[j] - sol2[j] );
        
        if( gap > highestGap && lx[j] != ux[j] )
        {
            highestGap = gap;
            indhighestGap = j;
        }
    }
    
    
}












int MRQ_FeasibilityPump::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    const int n = prob.n;
    const int m = prob.m;
    const int nI= prob.getNumberOfIntegerVars();
    const bool preprocess = in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_obj_function;
    
    bool flag, updtConstrBounds;
    int i, j, ret, nlast;
    unsigned int ind;
    unsigned int indLastSol;
    unsigned int nSubIter;
    unsigned long int iter = 0;
    double timeStart, highestGap;
    
    clock_t clockStart;
    
    MRQ_NLPSolver *nlp = NULL;
    //MRQ_NLPSolver *nlpfp = NULL;
    MRQ_NLPFeasPumpProb nlpfp;
    MRQ_Random random;
    MRQ_Preprocessor preprocessor(&prob);
    
    //araays:
    bool *auxEval = NULL;
    double *auxVars = NULL, *auxConstr = NULL;
    
    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux;
    double *plc = NULL, *puc = NULL;
    
    int *intVars = NULL;
    
    int **lastSols = NULL; //we will code the last solutions obtained in that matrix. Unfortunatelly, i could not do shift using long int...
    
    
    
    timeStart = MRQ_getTime();
    clockStart = clock();
    
    //zl = in_lower_bound;
    //zu = in_upper_bound;
    
    
    {
        auto ret = algorithmInitialization(1, (in_preprocess_lin_constr || in_preprocess_obj_function || in_preprocess_quad_constrs), milpSolverParams, nlpSolverParams, prob, lx, ux, &preprocessor, &updtConstrBounds, &plc, &puc); //that algorithm is monothread...
        if(ret != 0)
        {
            if(ret == MRQ_ALG_NOT_APPLICABLE)
            {
                if(in_print_level > 1)
                    std::cerr << MRQ_PREPRINT "Feasibility pump is not applicable to that problem currently!\n";
            }
            else
            {
                if(in_print_level > 0 )
                    std::cerr << MRQ_PREPRINT "Error " << ret << " at algorithm initialization\n";
            }
            
            out_return_code = ret;
            goto termination;
        }
    }
    
    
    if(in_print_level > 1)
        std::cout << "\n" MRQ_PREPRINT "Starting Feasibility Pump heuristic\n";
    
    if(in_print_level > 3)
        printSubSolvers(true, true, false);
    
    
    
    if( in_last_iters_considered_to_cycles <= 0)
        in_last_iters_considered_to_cycles = 1;
    
    
    MRQ_malloc(auxEval, m); 
    MRQ_malloc(auxConstr, m); 
    MRQ_malloc(auxVars, n); 
    MRQ_malloc(intVars, nI); 
    MRQ_calloc(lastSols, in_last_iters_considered_to_cycles); 
    
    if( !auxEval || !auxConstr || !auxVars || !intVars || !lastSols )
    {
        if(in_print_level > 0)
            MRQ_PRINTMEMERROR;
        
        out_return_code = MRQ_MEMORY_ERROR;
        goto termination;
    }
    
    
    for(ind = 0; ind < in_last_iters_considered_to_cycles; ind++)
    {
        MRQ_malloc(lastSols[ind], nI); //lastSols[ind] = (int *) malloc(nI * sizeof(int));
        
        if(!lastSols[ind])
        {
            if(in_print_level > 0)
                MRQ_PRINTMEMERROR;
            
            out_return_code = MRQ_MEMORY_ERROR;
            goto termination;
        }
    }
    
    
    for(i = 0; i < m; i++)
        auxEval[i] = true;
    
    prob.getIntegerIndices(intVars);
    
    
    
    if( !in_use_initial_solution || in_solve_nlp_as_local_search_at_end )
    {
        nlp = OPT_newNLPSolver( in_nlp_solver );
        
        if( !nlp )
        {
            if(in_print_level > 0)
                MRQ_PRINTMEMERROR;
        
            out_return_code = MRQ_MEMORY_ERROR;
            goto termination;
        }
        
        ret = MRQ_setNLPRelaxProb(prob, lx, ux, plc, puc, nlp, true, true, true, false, thnumber, in_set_special_nlp_solver_params, nlpSolverParams, in_number_of_threads, in_max_cpu_time, in_max_time, 0, 0);
        
        if( ret != 0 )
        {
            if( in_print_level > 0 )
                MRQ_PRINTERRORNUMBER(ret);
            
            out_return_code = MRQ_NLP_SOLVER_ERROR;
            goto termination;
        }
        
        /*ret = nlp->setnVariablesBounds(n, lx, ux);
        if( updtConstrBounds )
        {
            for(int i = 0; i < m; i++)
                ret += nlp->setConstraintBounds(i, plc[i], puc[i]);
        } */
        
        if( run_by_inside )
        {
            if( !std::isinf(insideSolverMaxTime) )
                ret += nlp->setMaxTime(insideSolverMaxTime );
        }
        
        
        if( ret != 0 )
        {
            if( in_print_level > 0 )
                std::cerr << MRQ_PREPRINT  "Error " MRQ_GETFILELINE << "\n";
        
            out_return_code = MRQ_NLP_SOLVER_ERROR;
            goto termination;
        }
        
        
        if( !in_use_initial_solution )
        {
            ret = nlp->solve(false);
            
            if( ret == OPT_OPTIMAL_SOLUTION )
            {
                out_obj_opt_at_continuous_relax = nlp->objValue;
                zl = MRQ_max( zl,  nlp->objValue);
                
                
                if(in_print_level > 1)
                    std::cout << MRQ_PREPRINT  "NLP relaxation solution: " << nlp->objValue << "\n";
                
                
                //if( prob.isIntegerSolution(nlp->sol, in_integer_tol) )
                if( MRQ_isIntegerSol(nI, intVars, nlp->sol, in_integer_tol) )
                {
                    if( in_print_level > 1 )
                        std::cout << MRQ_PREPRINT  "An integer optimal solution was gotten as NLP relaxation solution\n";
                    
                    //MRQ_copyArray(n, nlp->sol, out_best_sol);
                    //out_best_obj = zu = nlp->objValue;
                    
                    const bool r = tryUpdateBestSolution( thnumber, n, nlp->sol, nlp->objValue, iter, clockStart, timeStart, in_store_history_solutions );
                    
                    #if MRQ_DEBUG_MODE
                        assert(r == true);
                    #endif
                
                    out_return_code = MRQ_OPTIMAL_SOLUTION;
                    
                    goto termination;
                }
                
                MRQ_copyArray(n, nlp->sol, auxVars);
            }
            else if( nlp->feasSol )
            {
                MRQ_copyArray(n, nlp->sol, auxVars);
            }
            else if( ret == OPT_INFEASIBLE_PROBLEM )
            {
                out_return_code = MRQ_INFEASIBLE_PROBLEM;
                goto termination;
            }
            else
            {
                //we need a first solution...
                
                for(i = 0; i < n; i++)
                    auxVars[i] = MRQ_min( ux[i], MRQ_max(0.0, lx[i]) );
            }
            
        }
        
    }
    
    
    if(in_use_initial_solution)
    {
        MRQ_copyArray(n, xInit, auxVars);
    }
    
    
    MRQ_roundSol(nI, intVars, auxVars, auxVars);
    
    
    random.setSeed( in_use_random_seed_to_random_numbers ? NULL : &in_seed_to_random_numbers );
    
    
    while( true )
    {
        indLastSol = iter % in_last_iters_considered_to_cycles;
        
        iter++;
        
        
        
        if( iter > 1 )
        {
            nlast = MRQ_min<unsigned int>( iter-1, in_last_iters_considered_to_cycles);
            
            nSubIter = 0;
            
            while(true)
            {
                flag = MRQ_checkSolLastSols(nI, intVars, auxVars, nlast, lastSols);
                
                
                if(flag)
                {
                    //we have flip the solution 
                    
                    int r = MRQ_flipIntSol(random, in_lower_bound_to_pi, in_upper_bound_to_pi, lx, ux, nI, intVars, nlpfp.solver->sol, auxVars);
                    
                    if(r != 0)
                    {
                        out_return_code = MRQ_UNDEFINED_ERROR;
                        goto termination;
                    }
                }
                else
                {
                    break;
                }
                
                nSubIter++;
                
                if( nSubIter >= in_max_cycle_subiters )
                {
                    out_return_code = MRQ_HEURISTIC_FAIL;
                    goto termination;
                }
                
            }
        }
        
        
        
        MRQ_copySolToIntSol(nI, intVars, auxVars, lastSols[indLastSol] );
        
        
        //check if the rounded solution is already feasible.
        
        prob.isFeasibleToConstraints(thnumber, auxVars, true, auxEval, in_absolute_feasibility_tol, in_relative_feasibility_tol, flag, auxConstr);
        
        if( flag )
        {
            MRQ_copyArray( n, auxVars, out_best_sol );
            
            
            ret = prob.objEval(thnumber, !prob.hasNlConstrs, out_best_sol, out_best_obj);
            if(ret == 0)
            {
                if( out_best_obj < zu )
                    zu = out_best_obj;
                
                out_return_code = MRQ_HEURISTIC_SUCCESS;
                goto end_loop;
            }
        }
        
        
        //solving the nlp fp problem
        if( nlpfp.solver == NULL )
        {
            
            ret = nlpfp.setProblemBase( in_nlp_solver, prob, lx, ux, plc, puc, nlpSolverParams, thnumber, in_set_special_nlp_solver_params, in_set_linear_obj_term_on_bin_vars_at_nlp, in_set_norm1_on_nlp, in_number_of_threads, in_max_cpu_time, in_max_cpu_time );
            
            if( ret != 0 )
            {
                if(in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(ret);
                
                out_return_code = MRQ_MEMORY_ERROR;
                goto termination;
            }
            
            
            
            if( run_by_inside )
            {
                if( !std::isinf(insideSolverMaxTime) )
                    ret += nlpfp.solver->setMaxTime(insideSolverMaxTime );
                
                if( ret != 0 )
                {
                    if(in_print_level > 0)
                        MRQ_PRINTERROR;
                    
                    out_return_code = MRQ_MEMORY_ERROR;
                    goto termination;
                }
            }
            
        }
        
        
        ret = nlpfp.setObjective( nI, intVars, lx, ux, auxVars );
        
        if( ret != 0 )
        {
            if(in_print_level > 0)
                std::cerr << MRQ_PREPRINT "Error " << ret << MRQ_GETFILELINE << "\n";
            
            out_return_code = MRQ_MEMORY_ERROR;
            goto termination;
        }
        
        
        /*for(i = 30; i < n; i++)
            printf("y[%2d]: %0.1f \t", i-30, auxVars[i]);
        printf("\n"); */
        
        ret = nlpfp.solver->solve(false);
        
        if( ret == OPT_OPTIMAL_SOLUTION || nlpfp.solver->feasSol )
        {
            if( prob.isIntegerSolution( nlpfp.solver->sol, in_integer_tol ) )
            {
                ret = prob.objEval(thnumber, true, nlpfp.solver->sol, out_best_obj);
                
                if(ret == 0)
                {
                    if( out_best_obj < zu )
                        zu = out_best_obj;
                    
                    MRQ_copyArray( n, nlpfp.solver->sol, out_best_sol );
                    
                    out_return_code = MRQ_HEURISTIC_SUCCESS;
                    
                    goto end_loop;
                }
            }
            
            //printf("iter: %lu objF: %f\n", iter, nlpfp->objValue);
            //getchar();
            
            /*for(i = 0; i < n; i++)
                printf("x[%2d]: %0.3f \t", i-12, nlpfp->sol[i]);
            printf("\n");
            
            getchar(); */
            
            //MRQ_copyArray(n, nlpfp->sol, auxVars);
        }
        else
        {
            if(in_print_level > 1)
                std::cout << MRQ_PREPRINT "Error at NLP solving!\n";
            
            out_return_code = MRQ_NLP_SOLVER_ERROR;
            goto termination;
        }
        
        
        
        if( checkTerminationCriterions(thnumber, zl, zu, iter, timeStart, clockStart, out_return_code) )
        {
            goto termination;
        }
        
        
        
        
        MRQ_roundSol(nI, intVars, nlpfp.solver->sol, auxVars);
        
        //check if the rounded sol is equal to the last rounded sol, i. e., the current initial solution of nlpfp...
        
        flag = MRQ_checkSolLastSols(nI, intVars, auxVars, 1, &lastSols[indLastSol] );
        
        if(flag)
        {
            //"soft" flipping
            
            //calculating highest gap...
            
            
            MRQ_getHighestGapBetweenSols( lx, ux, nI, intVars, nlpfp.solver->sol, auxVars, ind, highestGap);
            
            
            #if MRQ_DEBUG_MODE
                //highestGap should be greater than in_integer_tol. Otherwise, we would have stoped after nlpfp->solve...
                assert( highestGap > in_integer_tol );
            #endif
            
            i = auxVars[ind];
                
            do
            {
                j = random.randInt( lx[ind], ux[ind] );
            }while( j == i );
            
            
            auxVars[ind] = j;
        }
        
    }
    
    
    
    
    
end_loop:

    
    if( out_return_code == MRQ_HEURISTIC_SUCCESS )
    {
        if( in_solve_nlp_as_local_search_at_end )
        {
            if(in_print_level > 2)
                std::cout << MRQ_PREPRINT "Current obj fucntion: " << out_best_obj << ".Solving nlp problem as local search.\n";
            
            
            if( preprocess )
            {
                bool updtvb, updtcb;
                
                //lx and ux are new arrays different of prob.lx and nlx
                
                for(int i = 0; i < nI; i++)
                {
                    const int ind = intVars[i];
                    lx[ind] = ux[ind] = out_best_sol[ind];
                }
                
                ret = preprocessor.preprocess( in_preprocess_quad_constrs, in_preprocess_obj_function, zu, lx, ux, updtvb, updtcb, plc, puc, plc, puc );
                
                
                ret = nlp->setnVariablesBounds( n, lx, ux );
                
                if( updtcb )
                {
                    for(int i = 0; i < m; i++)
                        ret += nlp->setConstraintBounds(i, plc[i], puc[i]);
                }
                
            }
            else
            {
                //MRQ_fixIntVarsOnSol( prob.n, prob.xtype, out_best_sol, *nlp );
                
                MRQ_fixIntVarsOnSolByList( nI, intVars, out_best_sol, *nlp);
            }
            
            ret = nlp->solve(false);
            
            if( ret == MRQ_OPTIMAL_SOLUTION || ret == MRQ_NLP_FEASIBLE_SOLUTION )
            {
                if(in_print_level > 2)
                    std::cout << MRQ_PREPRINT  "New objective function: " << nlp->objValue << "\n";
                
                if( nlp->objValue < out_best_obj )
                {
                    out_best_obj = nlp->objValue;
                    
                    if( out_best_obj < zu )
                        zu = out_best_obj;
                    
                    MRQ_copyArray(n, nlp->sol, out_best_sol);
                }
            }
            else
            {
                if(in_print_level > 2)
                    std::cout << "Failure.\n";
            }
            
        }
    }
    
    
    
    
    
    
    
    
termination:
    
    
    if( in_print_level > 1 )
    {
        if( out_return_code == MRQ_HEURISTIC_SUCCESS )
            std::cout << MRQ_PREPRINT  "Feasibility Pump heuristic found a feasible solution! ";
        else
            std::cout << MRQ_PREPRINT  "Feasibility Pump did not find a feasible solution! ";
    }
    
    
    if(nlp)	delete nlp;
    
    if(plc)			free(plc);
    if(auxEval)		free(auxEval);
    if(auxVars)		free(auxVars);
    if(auxConstr)	free(auxConstr);
    if(intVars)		free(intVars);
    
    if(lastSols)
    {
        for(ind = 0; ind < in_last_iters_considered_to_cycles; ind++)
        {
            if( lastSols[ind] )
                free(lastSols[ind]);
            else
                break;
        }
        
        free(lastSols);
    }
    
    
    
    
    out_feasible_solution = out_best_obj < MRQ_INFINITY;
    out_number_of_iterations = iter;
    //out_algorithm = muriqui::MRQ_FP_HEUR_ALG;
    out_lower_bound = zl;
    out_upper_bound = zu;
    
    algorithmFinalization(1, prob, lx, ux);
    
    out_cpu_time = MRQ_calcCPUTtime(clockStart);
    out_clock_time = MRQ_getTime() - timeStart;
    
    if(in_print_level > 1)
        std::cout << MRQ_PREPRINT "cpu time: " << out_cpu_time << "\n";
    
    
    
    return out_return_code;
}





































