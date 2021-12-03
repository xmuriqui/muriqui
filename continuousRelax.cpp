#include <cassert>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <cmath> //we use that to avoid problems on cmath on std 2011

#include <new>

#include "MRQ_algClasses.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_tools.hpp"


//using namespace std;

using namespace muriqui;
using namespace optsolvers;


MRQ_ContinuousRelax::MRQ_ContinuousRelax():MRQ_Algorithm()
{
    out_dual_sol = NULL;
    out_constraint_values = NULL;
    
    resetParameters();
    resetOutput();
    
    //ok, it is not a parameter, but I think it is a good idea set it here...
    out_algorithm = MRQ_CONT_RELAX_ALG;
}


MRQ_ContinuousRelax::~MRQ_ContinuousRelax()
{
    desallocateDualAndConstraintsArrays();
}



int MRQ_ContinuousRelax:: allocateDualAndConstraintsArrays( const int ndual, const int m)
{
    desallocateDualAndConstraintsArrays();
    
    MRQ_malloc(out_dual_sol, ndual); //out_dual_sol = (double *) malloc( ndual * sizeof(double) );
    MRQ_malloc(out_constraint_values, m); //out_constraint_values = (double *) malloc( m * sizeof(double) );
    MRQ_IFMEMERRORRETURN(!out_dual_sol || !out_constraint_values);
    
    MRQ_setAllArray<double>(ndual, out_dual_sol, NAN);
    MRQ_setAllArray<double>(m, out_constraint_values, NAN);
    
    return 0;
}



void MRQ_ContinuousRelax:: desallocateDualAndConstraintsArrays()
{
    MRQ_secFree( out_dual_sol);
    MRQ_secFree( out_constraint_values );
}



void MRQ_ContinuousRelax:: printParameters( std::ostream &out) const
{
    MRQ_Algorithm::printParameters(out);
    
    out << "\n"
    
    MRQ_STRFFATT(in_set_integer_vars_as_integers) << "\n";
}


void MRQ_ContinuousRelax::resetParameters()
{
    MRQ_Algorithm::resetParameters();
    
    in_set_integer_vars_as_integers = false;
    in_number_of_threads = 1;
}


int MRQ_ContinuousRelax::setIntegerParameter( const char *name, const long int value)
{
    int ret = MRQ_Algorithm::setIntegerParameter( name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_integer_vars_as_integers), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


void MRQ_ContinuousRelax::resetOutput()
{
    MRQ_Algorithm::resetOutput();
    
    
    desallocateDualAndConstraintsArrays();
    
    out_original_solver_return_code = INT_MIN;
    out_nlp_feasible_solution = false;
}




int MRQ_ContinuousRelax::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    const double timeStart = MRQ_getTime();
    const clock_t clockStart = clock();
    
    //const int n = prob.n;
    //const int m = prob.m;
    
    bool updtConstrBounds, intSol = false;
    int r;
    double *plc = NULL, *puc = NULL; //we do not allocate puc. Even so, we must initialzie it with NULL beacuse MRQ_setNLPRelaxProb. Take care to do not free puc!
    //MRQ_NLPSolver *nlp = NULL;
    MRQ_NLPSolver *nlp = NULL;
    MRQ_Preprocessor preprocessor(&prob);
    
    
    //arrays:
    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux;
    
    
    {
        auto ret = algorithmInitialization(1, (in_preprocess_lin_constr || in_preprocess_obj_function || in_preprocess_quad_constrs), milpSolverParams, nlpSolverParams, prob, lx, ux, &preprocessor, &updtConstrBounds, &plc, &puc); //that algorithm is monothread...
        if(ret != 0)
        {
            if(in_print_level > 0 && ret != MRQ_INFEASIBLE_PROBLEM )
                MRQ_PRINTERRORNUMBER(ret);
            
            out_return_code = ret;
            goto termination;
        }
    }
    
    //we do not need preprocessor more...
    preprocessor.deallocateMemory();
    
    
    /*for(int i = 0; i < m; i++)
    {
        std::cout << "lc["<<i<<"]: " << plc[i] << " uc["<<i<<"]: " << puc[i] << "  " << std::endl;
        
        if( plc[i] <= -MIP_INFINITY && puc[i] >= MIP_INFINITY )
        {
            std::cout << "Detectei restricao " << i << " redundante! " << std::endl;
            MRQ_getchar();
        }
        
    }
    MRQ_getchar(); */
    
    
    
    
    if(in_print_level > 1)
    {
        std::cout << "\nStarting Continuous Relaxation solving \n\n";
    
        //if(in_print_level > 3)
        printSubSolvers(false, true, false);
    }
    
    
    nlp = OPT_newNLPSolver( in_nlp_solver );
    MRQ_IFMEMERRORGOTOLABEL(!nlp, out_return_code, termination);
    
    
    r = MRQ_setNLPRelaxProb(prob, lx, ux, plc, puc, nlp, true, true, true, in_set_integer_vars_as_integers, thnumber, in_set_special_nlp_solver_params, nlpSolverParams, in_number_of_threads, in_max_cpu_time, in_max_time, 0, 0);
    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination)
    
    
    if( in_print_level > 5 )
        nlp->setOutputLevel(4);
    
    r = allocateDualAndConstraintsArrays( prob.m + 2*prob.n, prob.m );
    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_MEMORY_ERROR, termination);
    
    
    if( in_max_cpu_time < MRQ_INFINITY )
    {
        r = nlp->setMaxCPUTime( in_max_cpu_time );
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
    }
    
    if( in_max_time < MRQ_INFINITY )
    {
        r = nlp->setMaxTime( in_max_time );
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
    }
    
    if( run_by_inside )
    {
        if( !std::isinf(insideSolverMaxTime) )
        {
            r = nlp->setMaxTime(insideSolverMaxTime );
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination);
        }
    }
    
    
    if( in_use_initial_solution && xInit)
    {
        r = nlp->setInitialSolution(xInit, NULL, NULL);
        MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_NLP_SOLVER_ERROR, termination)
    }
    
    MRQ_secFree(plc);
    
    
    r = nlp->solve(false);
    
    if( nlp->retCode == OPT_OPTIMAL_SOLUTION )
    {
        if( !std::isnan(nlp->dualObjValue) )
            zl = out_lower_bound = nlp->dualObjValue;
        
        out_obj_opt_at_continuous_relax = nlp->objValue;
        
        if( prob.isIntegerSolution(nlp->sol, in_integer_tol) )
        {
            zu = out_upper_bound = nlp->objValue;
            intSol = true;
        }
        
        if(in_fix_int_vars_from_nlp_relax_sol)
        {
        }
    }
    
    out_original_solver_return_code = nlp->origSolverRetCode;
    out_best_obj = nlp->objValue;
    nlp->getSolution(out_best_sol, out_constraint_values);
    nlp->getDualSolution(out_dual_sol, &out_dual_sol[prob.m]);
    
    out_nlp_feasible_solution = nlp->feasSol;
    
    {
        r = nlp->getNumberOfIterations(out_number_of_iterations);
        if(r != 0)
        {
            //maybe some solvers does not provide number of iterations on optoslvers
            out_number_of_iterations = ULONG_MAX;
        }
    }
    
    
    switch(nlp->retCode)
    {
        case OPT_OPTIMAL_SOLUTION:
            out_return_code = intSol ? MRQ_OPTIMAL_SOLUTION : MRQ_CONT_RELAX_OPTIMAL_SOLUTION;
            break;
            
        case OPT_INFEASIBLE_PROBLEM:
            out_return_code = MRQ_INFEASIBLE_PROBLEM;
            break;
            
        case OPT_FEASIBLE_SOLUTION:
            out_return_code = MRQ_NLP_FEASIBLE_SOLUTION;
            break;
            
        case OPT_UNBOUNDED_PROBLEM:
            out_return_code = MRQ_UNBOUNDED_PROBLEM;
            break;
            
        case OPT_MAX_ITERATIONS:
            out_return_code = MRQ_MAX_ITERATIONS_STOP;
            break;
            
        case OPT_MAX_TIME:
            out_return_code = MRQ_MAX_TIME_STOP;
            break;
            
        case OPT_CALLBACK_FUNCTION_ERROR:
            out_return_code = MRQ_CALLBACK_FUNCTION_ERROR;
            break;
            
        default:
            out_return_code = MRQ_NLP_SOLVER_ERROR;
    }
    
    
    
termination:
    
    algorithmFinalization(1, prob, lx, ux);
    
    if(plc)		free(plc); //do not free puc
    if(nlp)		delete nlp;
    
    out_feasible_solution = zu == zl;
    //out_number_of_iterations = 1;
    
    out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
    out_clock_time = MRQ_getTime() - timeStart;
    
    if(in_print_level > 1)
        printf("cpu time: %f\n", out_cpu_time);
    
    
    return out_return_code;
    
}
