
#include <math.h>
#include <ctime>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <new>


#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_milpCallbacks.hpp"
#include "MRQ_advanced.hpp"

#include "OPT_solvers.hpp"



#define MRQ_MAX_SIZE_OFSTRING_TO_READ_PARAMS 512




using namespace muriqui;
using namespace optsolvers;
using namespace minlpproblem;



MRQ_Algorithm::MRQ_Algorithm()
{
    xInit = NULL;
    
    nlx = NULL;
    nux = NULL;
    
    out_best_sol = NULL;
    out_best_obj = INFINITY;
    
    
    run_by_inside = false;
    thnumber = 0;
    //insideNumberOfThreads = 1;
    insideSolverMaxTime = INFINITY;
    
    dcs0 = NULL;
    dcs1 = NULL;
    
    ndcs0 = 0;
    ndcs1 = 0;
    
    //resetParameters();
    resetOutput();
}


MRQ_Algorithm::~MRQ_Algorithm()
{
    desallocateBestSol();
    desallocateInitialSolution();
    desallocateDCS0();
    desallocateDCS1();
}


int MRQ_Algorithm::allocateDCS0(const unsigned int size)
{
    MRQ_malloc(dcs0, size); //dcs0 = (MRQ_DynConstrSetUnity*) malloc( size* sizeof(MRQ_DynConstrSetUnity) );
    
    if( !dcs0 )
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTMEMERROR;
        #endif
        return MRQ_MEMORY_ERROR;
    }
    
    return 0;
}


int MRQ_Algorithm::allocateDCS1(const unsigned int size)
{
    MRQ_malloc(dcs1, size); //dcs1 = (MRQ_DynConstrSetUnity*) malloc( size* sizeof(MRQ_DynConstrSetUnity) );
    
    if( !dcs1 )
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTMEMERROR;
        #endif
        return MRQ_MEMORY_ERROR;
    }
    
    return 0;
}


int MRQ_Algorithm::allocateBestSol(const unsigned int n)
{
    desallocateBestSol();
    
    MRQ_malloc(out_best_sol, n); //out_best_sol = (double *) malloc( n * sizeof(double) );
    if(!out_best_sol)
    {
        if( in_print_level > 0 )
            MRQ_PRINTMEMERROR;
        return MRQ_MEMORY_ERROR;
    }
    
    MRQ_setAllArray<double>(n, out_best_sol, NAN);
    
    out_best_obj = INFINITY;
    
    return 0;
}


int MRQ_Algorithm::allocateInitialSolution(const unsigned int n)
{
    //desallocateInitialSolution();
    //double *p;
    int r;
    
    //we use realloc. So, if user set another solution to same problem (i.e, having the same sine n), we try avoid reallocation;
    r = MRQ_realloc(xInit, n); //p = (double *) realloc( xInit, n*sizeof(double) );
    if(r != 0)//(!p)
    {
        if( in_print_level > 0 )
            MRQ_PRINTMEMERROR;
        return MRQ_MEMORY_ERROR;
    }
    
    //xInit = p;
    
    
    return 0;
}



int MRQ_Algorithm::checkAlgorithmRequirements(MRQ_MINLPProb& prob, const double* lx, const double* ux)
{
    return 0;
}


void MRQ_Algorithm::desallocateBestSol()
{
    
    if(out_best_sol)
    {
        free(out_best_sol);
        out_best_sol = NULL;
    }
}


void MRQ_Algorithm::desallocateDCS0()
{
    if( dcs0 )
    {
        free(dcs0);
        dcs0 = NULL;
        ndcs0 = 0;
    }
}


void MRQ_Algorithm::desallocateDCS1()
{
    if( dcs1 )
    {
        free(dcs1);
        dcs1 = NULL;
        ndcs1 = 0;
    }
}


void MRQ_Algorithm::desallocateInitialSolution()
{
    if(xInit)
    {
        free(xInit);
        xInit = NULL;
        in_use_initial_solution = false;
    }
}


void MRQ_Algorithm::copyParametersFrom(const MRQ_Algorithm& source)
{
    in_assume_convexity = source.in_assume_convexity;
    in_call_end_of_iteration_callback = source.in_call_end_of_iteration_callback;
    in_call_new_best_solution_callback = source.in_call_new_best_solution_callback;
    in_call_update_best_sol_callback= source.in_call_update_best_sol_callback;
    in_fix_int_vars_from_nlp_relax_sol = source.in_fix_int_vars_from_nlp_relax_sol;
    
    in_preprocess_lin_constr = source.in_preprocess_lin_constr;
    in_preprocess_quad_constrs = source.in_preprocess_quad_constrs;
    in_preprocess_obj_function = source.in_preprocess_obj_function;
    in_print_parameters_values= source.in_print_parameters_values;
    in_set_special_nlp_solver_params= source.in_set_special_nlp_solver_params;
    in_store_history_solutions= source.in_store_history_solutions;
    in_use_initial_solution= source.in_use_initial_solution;
    in_use_dynamic_constraint_set= source.in_use_dynamic_constraint_set;
    
    in_max_iterations= source.in_max_iterations;
    in_number_of_threads = source.in_number_of_threads;
    in_printing_frequency= source.in_printing_frequency;
    in_print_level= source.in_print_level;

    in_milp_solver= source.in_milp_solver;
    in_nlp_solver= source.in_nlp_solver;
    
    in_absolute_convergence_tol= source.in_absolute_convergence_tol;
    in_absolute_feasibility_tol= source.in_absolute_feasibility_tol;
    in_integer_tol= source.in_integer_tol;
    in_lower_bound= source.in_lower_bound;
    in_max_time= source.in_max_time;
    in_max_cpu_time= source.in_max_cpu_time;
    in_relative_convergence_tol= source.in_relative_convergence_tol;
    in_relative_feasibility_tol= source.in_relative_feasibility_tol;
    in_upper_bound= source.in_upper_bound;
    
    in_user_callbacks= source.in_user_callbacks;
}


void MRQ_Algorithm::getDCS0Array(int* size, MRQ_DynConstrSetUnity* dcs) const
{
    if(size)
        *size = ndcs0;
    
    if( dcs )
        MRQ_copyArray(ndcs0, dcs0, dcs);
}


//size or dcs can be NULL; We just return over arguments havinga valid pointers
void MRQ_Algorithm::getDCS1Array(int *size, MRQ_DynConstrSetUnity *dcs) const
{
    if(size)
        *size = ndcs1;
    
    if( dcs )
        MRQ_copyArray(ndcs1, dcs1, dcs);
}


bool MRQ_Algorithm::isLinearApproximationAlgorithm() const
{
    return false;
}


void MRQ_Algorithm::resetParameters()
{
    in_assume_convexity = true;
    in_call_end_of_iteration_callback = false;
    in_call_new_best_solution_callback = false;
    in_call_update_best_sol_callback = false;
    in_fix_int_vars_from_nlp_relax_sol = false;
    
    in_preprocess_lin_constr = true;
    in_preprocess_quad_constrs = true;
    in_preprocess_obj_function = true;
    in_print_parameters_values = false;
    in_store_history_solutions = false;
    in_use_initial_solution = false;
    
    in_set_special_nlp_solver_params = false;
    in_use_dynamic_constraint_set = false;
    
    in_print_level = 3;
    in_number_of_threads = 0;
    in_milp_solver = MRQ_getDefaultMILPSolverCode();
    in_nlp_solver = MRQ_getDefaultNLPSolverCode();
    
    in_integer_tol = 1.0e-4;
    in_absolute_feasibility_tol = 1.0e-3;
    in_relative_feasibility_tol = 1.0e-6;
    in_lower_bound = -MRQ_INFINITY; //do not use INFINITY
    in_upper_bound = MRQ_INFINITY;  //do not use INFINITY
    
    in_max_iterations = ULONG_MAX;//1000000000000;
    in_printing_frequency = 1;
    
    in_max_time = INFINITY;
    in_max_cpu_time = INFINITY;
    
    in_absolute_convergence_tol = 1.0e-3;
    in_relative_convergence_tol = 1.0e-3;
    
    in_user_callbacks = NULL;
}



void MRQ_Algorithm::printParameters(std::ostream &out) const
{
    char strValue[100];
    
    strValue[0] = '\0'; //only by safe...
    
    
    
    //we print algorithm code
    out << 
    // MRQ_STRFFATT(out_algorithm) << "\n"
    
    MRQ_STRFFATT(in_assume_convexity) << "\n"
    MRQ_STRFFATT(in_call_end_of_iteration_callback) << "\n"
    MRQ_STRFFATT(in_call_new_best_solution_callback) << "\n"
    MRQ_STRFFATT(in_call_update_best_sol_callback) << "\n"
    MRQ_STRFFATT(in_fix_int_vars_from_nlp_relax_sol) << "\n"
    
    MRQ_STRFFATT(in_preprocess_lin_constr) << "\n"
    MRQ_STRFFATT(in_preprocess_quad_constrs) << "\n"
    MRQ_STRFFATT(in_preprocess_obj_function) << "\n"
    MRQ_STRFFATT(in_print_parameters_values) << "\n"
    MRQ_STRFFATT(in_set_special_nlp_solver_params) << "\n"
    MRQ_STRFFATT(in_store_history_solutions) << "\n"
    MRQ_STRFFATT(in_use_dynamic_constraint_set) << "\n"
    MRQ_STRFFATT(in_use_initial_solution) << "\n"
    
    MRQ_STRFFATT(in_number_of_threads) << "\n"
    MRQ_STRFFATT(in_max_iterations) << "\n"
    MRQ_STRFFATT(in_printing_frequency) << "\n"
    MRQ_STRFFATT(in_print_level) << "\n";
    
    //MRQ_STRFFATT(in_milp_solver) << "\n"
    //MRQ_STRFFATT(in_nlp_solver) << "\n"
    
    MRQ_enumToStr(in_milp_solver, strValue) ;
    out << MRQ_STR(in_milp_solver) ": " << strValue << "\n";
    
    MRQ_enumToStr(in_nlp_solver, strValue) ;
    out << MRQ_STR(in_nlp_solver) ": " << strValue << "\n"
    
    
    MRQ_STRFFATT(in_absolute_convergence_tol) << "\n"
    MRQ_STRFFATT(in_relative_feasibility_tol) << "\n"
    MRQ_STRFFATT(in_integer_tol) << "\n"
    MRQ_STRFFATT(in_absolute_feasibility_tol) << "\n"
    MRQ_STRFFATT(in_lower_bound) << "\n"
    MRQ_STRFFATT(in_max_time) << "\n"
    MRQ_STRFFATT(in_max_cpu_time) << "\n"
    MRQ_STRFFATT(in_relative_convergence_tol) << "\n"
    MRQ_STRFFATT(in_upper_bound) << "\n"
    ;
}



int MRQ_Algorithm::setInitialSolution(const int n, const double *xI)
{
    int i; 
    
    if(n > 0)
    {
        i = allocateInitialSolution(n);
        
        if(i != 0)
            return i;
        
        //for(i = 0; i < n; i++)
            //this->xInit[i] = xI[i];
        
        MRQ_copyArray(n, xI, xInit);
        
        in_use_initial_solution = true;
    }
    else
    {
        desallocateInitialSolution();
    }
    
    return 0;
}



int MRQ_Algorithm::readParametersFromFile(const char *fileName, const bool printErrorMsgs)
{
    unsigned int lineCounter = 0;
    int r;
    FILE *file;
    char param[MRQ_MAX_SIZE_OFSTRING_TO_READ_PARAMS], value[MRQ_MAX_SIZE_OFSTRING_TO_READ_PARAMS];
    
    const int sline = 2*MRQ_MAX_SIZE_OFSTRING_TO_READ_PARAMS + 3;
    char line[sline];	
    
    
    file = fopen(fileName, "r");
    if( !file )
    {
        if(printErrorMsgs && in_print_level > 0)
            std::cerr << MRQ_PREPRINT "Error at opening paremeter file " << fileName << MRQ_GETFILELINE << "\n";
            
        return MRQ_NAME_ERROR;
    }
    
    
    
    
    
    while( !feof(file) )
    {
        bool set = false;
        bool valueError = false;
        char *p;
        
        lineCounter++;
        
        p = fgets(line, sline, file);
        if(p == NULL) //so, the end of file was reached before read any character. That is so strange because was not detect at while test...
            break;
        
        
        if( line[0] == MRQ_CHARAC_COMENT_ON_PARAMS_FILE ) //we assume there is only one peer to line
            continue;
        
        if( OPT_isEmptyString(line) ) //empty line
            continue;
        
        
        param[0] = '\0';
        value[0] = '\0';
        
        
        
        r = sscanf(line, "%" MRQ_EXPSTR(MRQ_MAX_SIZE_OFSTRING_TO_READ_PARAMS) "s " "%" MRQ_EXPSTR(MRQ_MAX_SIZE_OFSTRING_TO_READ_PARAMS) "s", param, value );
        
        if( r < 3 )
        {
            if( printErrorMsgs )
                std::cerr << MRQ_PREPRINT "Error to process line " << lineCounter << " in the file " << fileName << ": " << line; //do not \n here because line already have \n
            continue;
        }
        
        
        
        const int rs = setStringParameter(param, value);
        
        if( rs == 0 )
        {
            set = true;
        }
        else if( rs == MRQ_VALUE_ERROR )
        {
            valueError = true;
            set = true;
        }
        
        
        if( !set )
        {
            long int ivalue;
            
            sscanf(value, "%ld", &ivalue);
            const int ri = setIntegerParameter(param, ivalue);
            
            if( ri == 0 )
            {
                set = true;
            }
            else if( rs == MRQ_VALUE_ERROR )
            {
                valueError = true;
                set = true;
            }
        }
        
        
        if( !set )
        {
            double dvalue;
            
            sscanf(value, "%lf", &dvalue);
            const int rd = setDoubleParameter(param, dvalue);
            
            if( rd == 0 )
            {
                set = true;
            }
            else if( rd == MRQ_VALUE_ERROR )
            {
                valueError = true;
                set = true;
            }
        }
        
        
        
        if( in_print_level > 0 )
        {
            if(valueError)
                std::cerr << MRQ_PREPRINT "Value error at setting parameter " << param << "to value: " << value << "\n";
            else if(set)
                std::cout << MRQ_PREPRINT "Parameter " << param << " set to value " << value << "\n";
            else
                std::cerr << MRQ_PREPRINT "Name error at setting parameter " << param << "\n";
        }
    }
    
    
    fclose(file);
    
    return 0;
}




int MRQ_Algorithm::readParametersWithTypeFromFile(const char *fileName, const bool printErrorMsgs, const bool printFileOpenError)
{
    MRQ_AlgorithmParameterSetter algSetter(this);
    
    return OPT_readParametersWithTypeFromFile(fileName, printErrorMsgs, printFileOpenError, algSetter);
}



void MRQ_Algorithm::resetOutput()
{
    desallocateBestSol();
    out_best_obj = INFINITY;
    
    out_feasible_solution = false;
    
    out_number_of_feas_sols = 0;
    out_number_of_threads = 1;
    out_number_of_iterations = 0;
    out_user_callback_error_code = 0;
    out_return_code = MRQ_UNDEFINED_ERROR;
    //out_algorithm = MRQ_UNDEFINED_ALG; do not do it, since it overwrite real information inside each class
    
    out_number_of_iterations_to_first_feas_sol = -1; //we set the highest value
    out_number_of_iterations_to_best_sol = -1; //we set the highest value
    
    out_clock_time = 0.0;
    out_cpu_time = 0.0;
    out_cpu_time_to_first_feas_sol = -1.0;
    out_clock_time_to_first_feas_sol = -1.0;
    out_cpu_time_to_best_sol = -1.0;
    out_clock_time_to_best_sol = -1.0;
    out_lower_bound = -MRQ_INFINITY;
    out_upper_bound = MRQ_INFINITY;
    
    out_obj_opt_at_continuous_relax = NAN;
    
    out_sol_hist.desallocate();
}



int MRQ_Algorithm::setIntegerParameter(const char *name, const long int value)
{
    int ret = 0;
    
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_assume_convexity), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_call_end_of_iteration_callback), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_call_new_best_solution_callback), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_call_update_best_sol_callback), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_fix_int_vars_from_nlp_relax_sol), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_preprocess_lin_constr), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_preprocess_quad_constrs), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_preprocess_obj_function), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_print_parameters_values), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_store_history_solutions), name, value ) == 0 );
    
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_special_nlp_solver_params), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_number_of_threads), name, value) == 0 );
    else if( MRQ_setAtt<unsigned long int>( MRQ_STRATT(in_max_iterations), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_printing_frequency), name, value ) == 0 );
    else if( MRQ_setAtt<int>( MRQ_STRATT(in_print_level), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


int MRQ_Algorithm::setDoubleParameter(const char *name, const double value)
{
    int ret = 0;
    
    
    if( MRQ_setAtt( MRQ_STRATT(in_absolute_convergence_tol), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_relative_feasibility_tol), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_integer_tol), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_absolute_feasibility_tol), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_lower_bound), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_max_time), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_max_cpu_time), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_relative_convergence_tol), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_upper_bound), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}



int MRQ_Algorithm::setParameters(const MRQ_GeneralSolverParams &params)
{
    int r, ret = 0;
    
    for( auto &pairDbl : params.dblParams )
    {
        r = setDoubleParameter( pairDbl.first.c_str(), pairDbl.second );
        if( r != 0 )
        {
            if(in_print_level > 0)
                std::cerr << MRQ_PREPRINT "Error to set double parameter " << pairDbl.first << " to value " << pairDbl.second << "\n";
            
            ret = r;
        }
    }
    
    for( auto &pairInt : params.intParams )
    {
        r = setIntegerParameter( pairInt.first.c_str(), pairInt.second );
        if( r != 0 )
        {
            if(in_print_level > 0)
                std::cerr << MRQ_PREPRINT "Error to set integer parameter " << pairInt.first << " to value " << pairInt.second << "\n";
            
            ret = r;
        }
    }
    
    for( auto &pairStr : params.strParams )
    {
        r = setStringParameter( pairStr.first.c_str(), pairStr.second.c_str() );
        if( r != 0 )
        {
            if(in_print_level > 0)
                std::cerr << MRQ_PREPRINT "Error to set string parameter " << pairStr.first << " to value " << pairStr.second << "\n";
            
            ret = r;
        }
    }
    
    return ret;
}



//that method should be used to set enumeration parameters
int MRQ_Algorithm::setStringParameter(const char *name, const char *value)
{
    int r, ret;
    
    
    if( (r = MRQ_setStrAtt( MRQ_STRATT(in_milp_solver), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_nlp_solver), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}



















MRQ_RETURN_CODE MRQ_Algorithm::algorithmInitialization(const int nthreads, const bool allocBoundsCopy, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, MRQ_MINLPProb& prob, double* &lx, double* &ux, MRQ_Preprocessor *preprocessor, bool *updtSomeConstrBounds, double** preproclc, double** preprocuc, bool callUserCallbackBeforeAll)
{
    const int n = prob.n;
    const int m = prob.m;
    int r;
    
    MRQ_helloOnceTime();
    
    resetOutput();
    
    if( in_print_parameters_values )
    {
        if( in_print_level > 1 && !run_by_inside )
        {
            MRQ_PRINTMSG("parameters values:\n");
            printParameters();
            std::cout << "\n";
            MRQ_PRINTMSG("end of parameters values\n");
            
            std::cout << "\n";
            MRQ_PRINTMSG("first solver parameters values:\n\n");
            
            if(milpSolverParams)
                milpSolverParams->print();
            
            std::cout << "\n";
            MRQ_PRINTMSG("end of first solver parameters values\n"
            "second solver parameters values:\n\n");
            
            if(nlpSolverParams)
                nlpSolverParams->print();
            
            std::cout << "\n";
            MRQ_PRINTMSG("end of second solver parameters values\n");
        }
    }
    
    
    //this->prob = &prob;
    
    zl = MRQ_max(in_lower_bound, -MRQ_INFINITY);
    zu = MRQ_min(in_upper_bound, MRQ_INFINITY);
    
    
    if(updtSomeConstrBounds)
        *updtSomeConstrBounds = false;
    
    
    if( allocBoundsCopy )
    {
        double *oldlx = lx;
        double *oldux = ux;
        
        MRQ_malloc(lx, 2*n);
        MRQ_IFMEMERRORRETURN( !lx );
        
        
        ux = &lx[n];
        
        MRQ_copyArray( n, oldlx, lx );
        MRQ_copyArray( n, oldux, ux );
    }
    
    
    if( in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_obj_function )
    {
        bool updt, updtConstr;
        int r;
        
        
        if( preproclc )
        {
            MRQ_malloc( *preproclc, 2*m );
            MRQ_IFMEMERRORRETURN( !*preproclc);
            
            *preprocuc = &( (*preproclc)[m] );
        }
        
        
        r = preprocessor->allocateMemory( n, m );
        MRQ_IFERRORRETURN(r, MRQ_intToReturnCode(r));
        
        
        
        if( in_print_level >= 3 )
            std::cout << MRQ_PREPRINT << "preprocessing...\n";
        
        r = preprocessor->preprocess( in_preprocess_quad_constrs, in_preprocess_obj_function && zu < MRQ_INFINITY, zu, lx, ux, updt, updtConstr, NULL, NULL, preproclc ? *preproclc : NULL, preprocuc ? *preprocuc : NULL );
        
        
        if( updtSomeConstrBounds )
            *updtSomeConstrBounds = updtConstr;
        
        if( r == MIP_INFEASIBILITY )
        {
            std::cout << MRQ_PREPRINT << " minlpproblem preprocessor detected infeasibility\n";
            return MRQ_INFEASIBLE_PROBLEM;
        }
        
    }
    
    
    r = checkAlgorithmRequirements(prob, lx, ux);
    MRQ_IFERRORRETURN(r, MRQ_ALG_NOT_APPLICABLE);
    
    
    if( prob.nlEval && prob.hasNLTerm() && run_by_inside == false )
    {
        //we assume if we run by inside, i.e., by a call from other algorithm, the initialize method was already called
        r = prob.nlEval->initialize( nthreads, prob.n, prob.m, prob.J.getNumberOfElements(), prob.lagH.getNumberOfElements() );
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTCALLBACKERRORNUMBER(r);
            
            return MRQ_CALLBACK_FUNCTION_ERROR;
        }
    }
    
    if( in_user_callbacks )
    {
        //we do it because we can run some algorithms inside others (e.g. oa inside bb), and we must keep alg pointing to correct
        if( run_by_inside == false || in_user_callbacks->alg == NULL )
            in_user_callbacks->alg = this;
        
        if( callUserCallbackBeforeAll )
        {
            r = in_user_callbacks->beforeAll(out_algorithm, nthreads);
            MRQ_IFERRORRETURN(r, MRQ_STOP_REQUIRED_BY_USER);
        }
    }
    
    r = _algorithmInitialization();
    MRQ_IFERRORRETURN(r, MRQ_intToReturnCode(r));
    
    r = allocateBestSol(prob.n);
    MRQ_IFERRORRETURN(r, MRQ_intToReturnCode(r));
    
    return MRQ_SUCCESS;
}


int MRQ_Algorithm::_algorithmInitialization()
{
    return 0;
}


void MRQ_Algorithm::algorithmFinalization( const int nthreads, MRQ_MINLPProb& prob, double* &lx, double* &ux )
{
    const int n = prob.n;
    const int m = prob.m;
    
    
    out_lower_bound = zl;
    out_upper_bound = zu;
    out_feasible_solution = out_best_obj < MRQ_INFINITY;
    
    
    if( lx != prob.lx && lx != nlx )
    { //so, we allocate a new array to lx and ux in algorithm initialization, probably to do preprocessing
        free(lx);
        lx = NULL;
        ux = NULL;
    }
    
    if( prob.nlEval && prob.hasNLTerm() && run_by_inside == false )
    {
        prob.nlEval->finalize(1, n, m, prob.J.getNumberOfElements(), prob.lagH.getNumberOfElements());
    }
    
    prob.objValue = out_best_obj;
    if( out_best_obj < MRQ_INFINITY )
        MRQ_copyArray(n, out_best_sol, prob.x);
    
    
    if( in_user_callbacks )
    {
        in_user_callbacks->afterAll( out_algorithm, out_number_of_iterations, out_cpu_time, out_clock_time, out_lower_bound, out_upper_bound);
        
        if( in_user_callbacks->alg == this )
            in_user_callbacks->alg = NULL;
    }
    
}



bool MRQ_Algorithm::checkTerminationCriterions(const int threadNumber, const double zl, const double zu, long unsigned int iter, const double timeStart, const clock_t &clockStart, MRQ_RETURN_CODE& retCode) const
{
    bool timesComputed = false;
    double cpuTime;
    double wallTime;
    
    
    if(in_call_end_of_iteration_callback)
    {
        timesComputed = true;
        
        //we try avoid calculating times unnecessary (I believe it call system functions and that is not so good...)
        cpuTime = MRQ_calcCPUTtime( clockStart, clock());//(double(clock()-clockStart) )/CLOCKS_PER_SEC;
        wallTime = MRQ_getTime() - timeStart;
        
        if( in_user_callbacks && in_user_callbacks->endOfIteration( out_algorithm, threadNumber, iter, cpuTime, wallTime, zl, zu) != 0)
        {
            retCode = MRQ_STOP_REQUIRED_BY_USER;
            return true;
        }
    }
    
    
    
    //we already have a feasible solution...
    if(zu - zl <= in_absolute_convergence_tol || zu - zl <= MRQ_abs(zu)*in_relative_convergence_tol)
    {
        if(out_best_obj < MRQ_INFINITY)
        {
            retCode = MRQ_OPTIMAL_SOLUTION;
            //return true;
        }
        else
        {
            retCode = MRQ_INFEASIBLE_PROBLEM;
            //return true;
        }
        
        return true;
    }
    
    
    
    if( iter >= in_max_iterations )
    {
        retCode = MRQ_MAX_ITERATIONS_STOP;
        return true;
    }
    
    if( in_max_cpu_time < INFINITY )
    {
        if(!timesComputed)
            cpuTime = MRQ_calcCPUTtime( clockStart, clock());//(double(clock()-clockStart) )/CLOCKS_PER_SEC;
        
        
        if( cpuTime >= in_max_cpu_time )
        {
            retCode = MRQ_MAX_TIME_STOP;
            return true;
        }
    }
    
    if( in_max_time < INFINITY )
    {
        if(!timesComputed)
            wallTime = MRQ_getTime() - timeStart;
        
        if( wallTime >= in_max_time )
        {
            retCode = MRQ_MAX_TIME_STOP;
            return true;
        }
    }
    
    
    return false;
}


std::string MRQ_Algorithm::getAlgorithmInitials() const 
{
    return MRQ_getAlgInitials(out_algorithm);
}


std::string MRQ_Algorithm::getAlgorithmName() const
{
    return MRQ_getAlgName(out_algorithm);
}


int MRQ_Algorithm::getBestSolutionCopy(const int n, double* solution, double& fsolution) const
{
    //const int n = prob->n;
    
    fsolution = out_best_obj;
    MRQ_copyArray(n, out_best_sol, solution);
    
    return 0;
}


/*
//method to try get a good feasible solution. Return  MRQ_HEURISTIC_SUCCESS or MRQ_OPTIMAL_SOLUTION if a feasible solution is found and false otherwise
int MRQ_Algorithm::getFeasibleSolution(MRQ_MINLPProb &prob, double *lx, double *ux, const double *xInitial, const double maxTime, const bool stopOnFirst, double *sol, double &obj)
{
    MRQ_Diving diving;
    MRQ_FeasibilityPump fp;
    MRQ_OAFeasibilityPump oafp;
    
    obj = INFINITY;
    
    
    if(xInitial)
        fp.setInitialSolution(prob.n, xInitial);
    
    
    fp.insideRun(prob, NULL, NULL, thnumber, 1, maxTime, lx, ux);
    
    if( fp.out_return_code == MRQ_HEURISTIC_SUCCESS || fp.out_return_code == MRQ_OPTIMAL_SOLUTION )
    {
        obj = fp.out_best_obj;
        MRQ_copyArray(prob.n,  fp.out_best_sol, sol);
        
        if( stopOnFirst || fp.out_return_code == MRQ_OPTIMAL_SOLUTION )
            return fp.out_return_code;
    }
    
    
    if(xInitial)
        diving.setInitialSolution(prob.n, xInitial);
    
    
    diving.insideRun(prob, NULL, NULL, thnumber, 1, maxTime, lx, ux);
    
    if( diving.out_return_code == MRQ_HEURISTIC_SUCCESS || diving.out_return_code == MRQ_OPTIMAL_SOLUTION )
    {
        
        if( diving.out_best_obj < obj )
        {
            obj = diving.out_best_obj;
            MRQ_copyArray(prob.n, diving.out_best_sol, sol);
            
            if( stopOnFirst || diving.out_return_code == MRQ_OPTIMAL_SOLUTION )
                return diving.out_return_code;
        }
    }
    
    
    if(xInitial)
        oafp.setInitialSolution(prob.n, xInitial);
    
    
    oafp.insideRun(prob, NULL, NULL, thnumber, 1, maxTime, lx, ux);
    
    if( oafp.out_return_code == MRQ_HEURISTIC_SUCCESS || oafp.out_return_code == MRQ_OPTIMAL_SOLUTION )
    {
        
        if( oafp.out_best_obj < obj )
        {
            obj = oafp.out_best_obj;
            MRQ_copyArray(prob.n, oafp.out_best_sol, sol);
            
            if( stopOnFirst || oafp.out_return_code == MRQ_OPTIMAL_SOLUTION )
                return oafp.out_return_code;
        }
    }
    
    
    return obj < INFINITY ? MRQ_HEURISTIC_SUCCESS : MRQ_HEURISTIC_FAIL;
}
*/








int MRQ_Algorithm::insideRun(MRQ_MINLPProb& prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, const int thnumber, const double insideSolverMaxTime, double* nlx, double* nux)
{
    this->thnumber = thnumber;
    run_by_inside = true;
    //this->insideNumberOfThreads = nThreads;
    this->insideSolverMaxTime = insideSolverMaxTime;
    
    this->nlx = nlx;
    this->nux = nux;
    
    run(prob, milpSolverParams, nlpSolverParams);
    
    this->thnumber = 0;
    run_by_inside = false;
    //this->insideNumberOfThreads = 1;
    this->insideSolverMaxTime = INFINITY;
    
    this->nlx = NULL;
    this->nux = NULL;
    
    return out_return_code;
}



void MRQ_Algorithm::printSubSolvers(const bool printMilp, const bool printNlp, const bool printGlobal)
{
    std::cout << MRQ_PREPRINT ;
    if(printMilp)
    {
        std::cout << "milp solver: " << optsolvers::OPT_getSolverName(in_milp_solver);
        /*switch(in_milp_solver)
        {
            case MRQ_CPLEX:
                std::cout << "cplex";
                break;
            
            case MRQ_GUROBI:
                std::cout << "gurobi";
                break;
            
            case MRQ_XPRESS:
                std::cout << "xpress";
                break;
                
            case MRQ_MILP_MOSEK:
                std::cout << "mosek";
                break;
            case MRQ_MILP_KNITRO:
                std::cout << "knitro";
                break;
            case MRQ_GLPK:
                std::cout << "glpk";
                break;
            
            default:
                std::cout << "unknow";
        }
        */
        std::cout << ", ";
    }
    
    if( printNlp )
    {
        std::cout << "nlp solver: " << optsolvers::OPT_getSolverName(in_nlp_solver);
        /*switch(in_nlp_solver)
        {
            case MRQ_NLP_MOSEK:
                std::cout << "mosek";
                break;
                
            case MRQ_IPOPT:
                std::cout << "ipopt";
                break;
            case MRQ_NLP_KNITRO:
                std::cout << "knitro";
                break;
            case MRQ_WORHP:
                std::cout << "worhp";
                break;
            case MRQ_ALGENCAN:
                std::cout << "algencan";
                break;
            default:
                std::cout << "unknow";
        }*/
    }
    
    std::cout << "\n";
}





bool MRQ_Algorithm::tryUpdateBestSolution(const int unsigned threadNumber, const int n, double* sol, double objValue, const long unsigned int iter, const clock_t &clockStart, const double timeStart, const bool storeOnHistory)
{
    return tryUpdateBestSolutionWithLock(threadNumber, 1, NULL, n, sol, objValue, iter, clockStart, timeStart, storeOnHistory);
}



bool MRQ_Algorithm::tryUpdateBestSolutionWithLock(const unsigned int threadNumber, const unsigned int nthreads, MRQ_Mutex *mutex, const int n, double* sol, double objValue, const long unsigned int iter, const clock_t &clockStart, const double timeStart, const bool storeOnHistory)
{
    bool updt = true, change = false;
    decltype(out_number_of_feas_sols) nfeas;
    double oldBestSol;
    double clockTime, cpuTime;
    //const int n = prob->n;
    
    
    if( in_call_update_best_sol_callback )
    {
        const int r = in_user_callbacks->updatingBestSolution( out_algorithm, threadNumber, sol, objValue, zu, iter);
        
        if( r != 0 )
        {
            if( in_print_level > 0 )
                std::cerr << MRQ_PREPRINT "Callback function error " << r << "on in_user_callbacks->updatingBestSolution" << MRQ_GETFILELINE << "\n";
            
            updt = false;
        }
    }
    
    if( updt && objValue < zu ) //note, the strict correct is to compare objValue with zu inside the mutual exclusion zone. However, zu never increases its value. So, if objValue is already greater than zu, so, we do not need set semaphore to check in a strict correct way.
    {
        cpuTime = MRQ_calcCPUTtime(clockStart, clock()); //(double(clock()-clockStart) )/CLOCKS_PER_SEC;
        clockTime = MRQ_getTime() - timeStart;
        
        if(mutex)
            mutex->lock(nthreads);
        { //mutual exclusion zone
            if( updt && objValue < zu )
            {
                oldBestSol = out_best_obj;
                out_best_obj = zu = objValue;
                MRQ_copyArray(n, sol, out_best_sol);
                
                out_number_of_feas_sols++;
                
                out_cpu_time_to_best_sol = cpuTime;
                out_clock_time_to_best_sol = clockTime;
                out_number_of_iterations_to_best_sol = iter;
                
                nfeas = out_number_of_feas_sols;
                change = true;
            }
            
        }
        if(mutex)
            mutex->unlock(nthreads);
    }
    
    
    if(change)
    {
        
        if( nfeas == 1 )
        {
            out_cpu_time_to_first_feas_sol = cpuTime;
            out_clock_time_to_first_feas_sol = clockTime;
            out_number_of_iterations_to_first_feas_sol = iter;
        }
        
        
        if( storeOnHistory )
        {
            const int r = out_sol_hist.addSolution( n, iter, clockTime, cpuTime, sol, objValue );
            
            if( r != 0 )
            {
                if( in_print_level > 0 )
                    std::cerr << MRQ_PREPRINT "Warning: Error " << r << " on out_sol_hist.addSolution." << MRQ_GETFILELINE << "\n";
            }
        
        }
    
        if( in_call_new_best_solution_callback && in_user_callbacks )
        {
            in_user_callbacks->newBestSolution( threadNumber, sol, oldBestSol, objValue, iter);
        }
    }
    
    
    
    
    return change;
}









MRQ_LinearApproxAlgorithm::MRQ_LinearApproxAlgorithm():MRQ_Algorithm()
{
    nPoints = 0;
    points = NULL;
    
    
    /********** for pseudo prunning *******/
    pp_nvars = 0;
    pp_nI = 0;
    pp_nthreads = 0;
    pp_intVars = NULL;
    pp_nBranchs = NULL;
    pp_counter1 = NULL;
    pp_counter2 = NULL;
    pp_lowestBoundPseudoPruned = NULL;
    pp_lupcosts = NULL;
    pp_sols = NULL;
    /********** end for pseudo prunning *******/
    
    resetOutput();
}


//destructors of base class are ever called by its descendats
MRQ_LinearApproxAlgorithm::~MRQ_LinearApproxAlgorithm()
{
    deletePointsToLinearisation();
    deallocateMemoryForPseudoPruning();
}


void MRQ_LinearApproxAlgorithm::printParameters(std::ostream& out) const
{
    char strValue[100];
    
    strValue[0] = '\0'; //only by safe if some parameter does not find the string value...
    
    
    MRQ_Algorithm::printParameters(out);
    out << "\n"
    
    MRQ_STRFFATT(in_call_before_solve_callback_in_milp_bb) << "\n"
    MRQ_STRFFATT(in_call_branching_callback_in_milp_bb) << "\n"
    MRQ_STRFFATT(in_measure_nlp_time) << "\n"
    MRQ_STRFFATT(in_set_obj_lower_bound_on_master_problem) << "\n"
    MRQ_STRFFATT(in_set_quadratics_in_master_problem) << "\n";
    
    
    MRQ_enumToStr( in_obj_linearisation_strategy, strValue );
    out << MRQ_STR(in_obj_linearisation_strategy) ": " << strValue << "\n";
    
    MRQ_enumToStr( in_constr_linearisation_strategy, strValue );
    out << MRQ_STR(in_constr_linearisation_strategy) ": " << strValue << "\n";
    
    MRQ_enumToStr( in_quad_app_master_strategy, strValue );
    out << MRQ_STR(in_quad_app_master_strategy) ": " << strValue << "\n"
    
    MRQ_STRFFATT(in_eps_to_active_constr_to_linearisation) << "\n"
    
    MRQ_STRFFATT(in_try_pseudo_pruning_before_solving_relaxations) << "\n";
    
    //MRQ_STRFFATT(in_min_number_of_brachings_per_var_before_pseudo_pruning) << "\n";
    
    
    MRQ_enumToStr( in_pseudo_pruning_strategy, strValue );
    out << MRQ_STR(in_pseudo_pruning_strategy) ": " << strValue << "\n"
    
    MRQ_STRFFATT(in_alpha_to_balance_estimative_in_pseudo_pruning) << "\n"
    MRQ_STRFFATT(in_absolute_convergence_tol_for_pseudo_pruning) << "\n"
    MRQ_STRFFATT(in_relative_convergence_tol_for_pseudo_pruning) << "\n"
    MRQ_STRFFATT(in_absolute_upper_bound_slack_for_pseudo_pruning) << "\n"
    MRQ_STRFFATT(in_relative_upper_bound_slack_factor_for_pseudo_pruning) << "\n"
    ;
}


void MRQ_LinearApproxAlgorithm::resetOutput()
{
    MRQ_Algorithm::resetOutput();
    
    out_number_of_constr_linears_saved = 0;
    out_number_of_obj_linears_saved_by_infeasibility = 0;
    out_number_of_obj_linears_saved_by_zu = 0;
    
    out_number_of_nlp_probs_solved = 0;
    out_number_of_milp_solver_iters = 0;
    
    out_cpu_time_on_box_to_constr_calculation = 0;
    
    out_clock_time_of_nlp_solving = NAN;
    out_cpu_time_of_nlp_solving = NAN;
    
    out_number_of_pseudo_prunes_on_solving = 0;
    out_number_of_pseudo_prunes_on_branching = 0;
}


void MRQ_LinearApproxAlgorithm::resetParameters()
{
    MRQ_Algorithm::resetParameters();
    
    in_call_before_solve_callback_in_milp_bb = false;
    in_call_branching_callback_in_milp_bb = false;
    in_measure_nlp_time = false;
    in_set_obj_lower_bound_on_master_problem = false;
    in_set_quadratics_in_master_problem = false;
    in_quad_app_master_strategy = MRQ_QAMS_NO_QUAD_APP;
    in_obj_linearisation_strategy = MRQ_OLS_ALL_POINTS;
    in_constr_linearisation_strategy = MRQ_CLS_ALL_CONSTRS;
    in_eps_to_active_constr_to_linearisation = 1.0e-3;
    
    /* for pseudo prunning */
    in_try_pseudo_pruning_before_solving_relaxations = false;
    
    //in_min_number_of_brachings_per_var_before_pseudo_pruning = 0;
    
    in_pseudo_pruning_strategy = MRQ_BB_PPS_NO_PSEUDO_PRUNING;
    
    in_alpha_to_balance_estimative_in_pseudo_pruning = 0.9;
    
    in_absolute_convergence_tol_for_pseudo_pruning = 1.0;
    in_relative_convergence_tol_for_pseudo_pruning = 0.05;
    
    in_absolute_upper_bound_slack_for_pseudo_pruning = 1e-3;
    in_relative_upper_bound_slack_factor_for_pseudo_pruning = 1e-3;
    /* end for pseudo prunning */
}


int MRQ_LinearApproxAlgorithm::setIntegerParameter( const char *name, const long int value)
{
    int ret = MRQ_Algorithm::setIntegerParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_obj_lower_bound_on_master_problem), name, value ) == 0 );
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_measure_nlp_time), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_quadratics_in_master_problem), name, value ) == 0 );
    /********  for pseudo pruning *********/
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_try_pseudo_pruning_before_solving_relaxations), name, value ) == 0 );
    /********  end of for pseudo pruning *********/
    else
        ret = MRQ_NAME_ERROR;
    
    return ret;
}



int MRQ_LinearApproxAlgorithm::setDoubleParameter( const char *name, const double value)
{
    int ret = MRQ_Algorithm::setDoubleParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt( MRQ_STRATT(in_eps_to_active_constr_to_linearisation), name, value ) == 0 );
    /********  for pseudo pruning *********/
    else if( MRQ_setAtt( MRQ_STRATT(in_alpha_to_balance_estimative_in_pseudo_pruning), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_absolute_convergence_tol_for_pseudo_pruning), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_relative_convergence_tol_for_pseudo_pruning), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_absolute_upper_bound_slack_for_pseudo_pruning), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_relative_upper_bound_slack_factor_for_pseudo_pruning), name, value ) == 0 );
    /******** end for pseudo pruning *********/
    else
        ret = MRQ_NAME_ERROR;
    
    return ret;
}



int MRQ_LinearApproxAlgorithm::setStringParameter( const char *name, const char *value)
{
    int r, ret = MRQ_Algorithm::setStringParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    
    if( (r = MRQ_setStrAtt( MRQ_STRATT(in_obj_linearisation_strategy), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_constr_linearisation_strategy), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_quad_app_master_strategy), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_pseudo_pruning_strategy), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}



int MRQ_LinearApproxAlgorithm:: addLazyConstraintsLinearizationOnSolution( const unsigned int thnumber, MRQ_MILPSolverCallbackInterface *callbackSolver, MRQ_MINLPProb &prob, MRQ_GradientsEvaluation &gradEval, const bool incQuadsInMaster, const bool linearizeObj, const bool *constrEval, const double *psol, const double *objSol, const double *pconstr, const double *masterConstr, const int *indices, const double *plc, const double *puc, double *auxVars, bool *auxConstrEval2)
{
    bool newx = true;
    int r, retCode = 0;
    
    
    if(linearizeObj)
    {
        r = MRQ_addLazyConstraintObjectiveLinearizationOnSolution(thnumber, callbackSolver, prob, incQuadsInMaster, newx, psol, objSol, indices, auxVars);
        MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
        
    }
    
    
    //linearising constraints
    r = MRQ_addLazyConstraintConstraintLinearizationOnSolution(thnumber, callbackSolver, prob, gradEval, incQuadsInMaster, constrEval, newx, psol, pconstr, masterConstr, indices, plc, puc, auxVars, auxConstrEval2, in_constr_linearisation_strategy, in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved);
    MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    
    
termination:
    
    return retCode;
}


void MRQ_LinearApproxAlgorithm:: deallocateMemoryForPseudoPruning()
{
    for(unsigned int i = 0; i < pp_nthreads; i++)
    {
        if( pp_lupcosts[i] )
            free( pp_lupcosts[i] );
    }
    
    for(unsigned int i = 0; i < pp_nthreads; i++)
    {
        if( pp_sols[i] )
            free(pp_sols[i]);
    }
    
    MRQ_secFree(pp_intVars);
    MRQ_secFree(pp_lupcosts);
    MRQ_secFree(pp_nBranchs);
    MRQ_secFree(pp_counter1);
    MRQ_secFree(pp_counter2);
    MRQ_secFree(pp_lowestBoundPseudoPruned);
    MRQ_secFree(pp_sols);
    
    pp_nthreads = 0;
    pp_nvars = 0;
    pp_nI = 0;
}


//nvars and nI are used to have a maximum value for the arrays since presolve can eliminate variablesi (if lazy constraints are not nin use, for example, in linear problems)
int MRQ_LinearApproxAlgorithm:: allocateMemoryForPseudoPruning( unsigned int nthreads, unsigned int nvars, unsigned int nI)
{
    deallocateMemoryForPseudoPruning();
    
    pp_nthreads = nthreads;
    
    //cplex creates a one more variable to objective (I think), since we will get values from cplex getting the number of varition from cplex lp pointer problems, we have to take this additional variable
    nvars++;
    
    MRQ_malloc( pp_intVars, nI );
    MRQ_calloc( pp_nBranchs, nvars );
    MRQ_calloc( pp_lupcosts, nthreads );
    MRQ_calloc( pp_sols, nthreads );
    MRQ_calloc( pp_counter1, nthreads );
    MRQ_calloc( pp_counter2, nthreads );
    MRQ_malloc( pp_lowestBoundPseudoPruned, nthreads );
    
    MRQ_IFMEMERRORRETURN( !pp_nBranchs || !pp_lupcosts || !pp_sols || !pp_counter1 || !pp_counter2 || !pp_lowestBoundPseudoPruned );
    
    
    //printf("nvars: %d\n", nvars);
    //MRQ_getchar();
    
    for(unsigned int i = 0; i < nthreads; i++)
    {
        MRQ_malloc( pp_lupcosts[i], 2*nvars );
        MRQ_malloc( pp_sols[i], nvars );
        MRQ_IFMEMERRORRETURN( !pp_lupcosts[i] || !pp_sols[i] );
    }
    
    MRQ_setAllArray<double>(nthreads, pp_lowestBoundPseudoPruned, INFINITY);
    
    //pp_nvars and pp_nI will be set by the callbacks procedures
    
    return 0;
}


int MRQ_LinearApproxAlgorithm::checkAlgorithmRequirements( MRQ_MINLPProb& prob, const double* lx, const double* ux)
{
    if( prob.getNumberOfNLEqualityConstraints() > 0 )
    {
        if(in_print_level > 1)
            MRQ_PRINTERRORMSG("Error: Linear Approximation algorithm does not address problems having nonlinear equality constraints!\n");
            
        return MRQ_ALG_NOT_APPLICABLE;
    }
    
    return 0;
}


int MRQ_LinearApproxAlgorithm::_algorithmInitialization()
{
    int r = MRQ_Algorithm::_algorithmInitialization();
    
    if(in_measure_nlp_time)
    {
        out_clock_time_of_nlp_solving = 0.0;
        out_cpu_time_of_nlp_solving = 0.0;
    }
    
    return r;
}



void MRQ_LinearApproxAlgorithm::algorithmFinalization( const int nthreads, MRQ_MINLPProb& prob, double* &lx, double* &ux )
{
    
    if( in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
    {
        double minlb = INFINITY;
        
        for(unsigned int i = 0; i < pp_nthreads; i++)
            out_number_of_pseudo_prunes_on_solving += pp_counter1[i];
        
        for(unsigned int i = 0; i < pp_nthreads; i++)
            out_number_of_pseudo_prunes_on_branching += pp_counter2[i];
        
        for(unsigned int i = 0; i < pp_nthreads; i++)
        {
            if( pp_lowestBoundPseudoPruned[i] < minlb )
                minlb = pp_lowestBoundPseudoPruned[i];
        }
        
        
        if( !std::isinf(minlb) )
        {
            if( minlb < zl )
                zl = minlb;
        }
        
        deallocateMemoryForPseudoPruning();
    }
    
    
    return MRQ_Algorithm::algorithmFinalization(nthreads, prob, lx, ux);
}


int MRQ_LinearApproxAlgorithm::deletePointsToLinearisation()
{
    if(points)
    {
        for(int i = 0; i < nPoints; i++)
            free(points[i]);

        free(points);
        points = NULL;
    }
    
    nPoints = 0;
    return 0;
}


bool MRQ_LinearApproxAlgorithm:: isLinearApproximationAlgorithm() const
{
    return true;
}


int MRQ_LinearApproxAlgorithm::addPointToLinearisation(const int dim,  const double* point)
{
    void *p = (void *) &point; //this is ridiculous, but C has serius problems about casting betwen bilevel pointers and const bilevel pointers...
    
    return addPointsToLinearisation(1, dim, (double **) p);
}


int MRQ_LinearApproxAlgorithm::addPointsToLinearisation(const int npoints, const int dim, double **points)
{
    int i, j, k, r;
    //double **aux;

    
    
    //aux = (double**)realloc( this->points, (nPoints + npoints)*sizeof(double*));
    
    r = MRQ_realloc(this->points, (nPoints + npoints));
    if(r != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTMEMERROR;
        
        return MRQ_MEMORY_ERROR;
    }
    //this->points = aux;
    

    for(k =0, i = nPoints; i < npoints + nPoints; k++, i++)
    {
        
        MRQ_malloc(this->points[i], dim); //this->points[i] = (double*) malloc( dim * sizeof(double) );
        if( !this->points[i] )
        {
            if(in_print_level > 0)
                MRQ_PRINTMEMERROR;
            
            #pragma GCC ivdep
            #pragma ivdep
            for(j = nPoints; j < i; j++ )
                free( this->points[j] );
                
            return MRQ_MEMORY_ERROR;
        }

        //for(j = 0; j < dim; j++)
        //this->points[i][j] = points[k][j];
        
        MRQ_copyArray(dim, points[k], this->points[i]);
    }

    nPoints += npoints;
    
    return 0;
}



int MRQ_LinearApproxAlgorithm::addPointsToLinearisation(const int dim, MRQ_SolutionHistory& points)
{
    const int npoints = points.getnsols();
    int r;
    //double **aux;
    MRQ_HistorySolution *solAux;
    
    if(npoints == 0)
        return 0;
    
    
    //aux = (double**)realloc( this->points, (nPoints + npoints)*sizeof(double*));
    r = MRQ_realloc(this->points, (nPoints + npoints) );
    if(r != 0)
    {
        if(in_print_level > 0)
            MRQ_PRINTMEMERROR;
        
        return MRQ_MEMORY_ERROR;
    }
    //this->points = aux;
    
    
    for(int k = 0, i = nPoints; i < npoints + nPoints; k++, i++)
    {
        //fprintf(stderr, "k: %d\n", k);
        
        MRQ_malloc(this->points[i], dim); //this->points[i] = (double*) malloc( dim * sizeof(double) );
        if( !this->points[i] )
        {
            if(in_print_level > 0)
                MRQ_PRINTMEMERROR;
            
            for(int j = nPoints; j < i; j++ )
                free( this->points[j] );
            
            return MRQ_MEMORY_ERROR;
        }
        
        solAux = points.getHistSolPointer(k);
        
        solAux->getsolution(this->points[i]) ;
    }
    
    
    nPoints += npoints;
    
    return 0;
}


int MRQ_LinearApproxAlgorithm::getNPoints()
{
    return nPoints;
}


void MRQ_LinearApproxAlgorithm::copyParametersFrom(const MRQ_LinearApproxAlgorithm &source)
{
    MRQ_Algorithm::copyParametersFrom(source);
    
    in_measure_nlp_time = source.in_measure_nlp_time;
    
    in_set_obj_lower_bound_on_master_problem = source.in_set_obj_lower_bound_on_master_problem;
    in_set_quadratics_in_master_problem = source.in_set_quadratics_in_master_problem;
    in_obj_linearisation_strategy = source.in_obj_linearisation_strategy;
    in_constr_linearisation_strategy = source.in_constr_linearisation_strategy;
    in_quad_app_master_strategy = source.in_quad_app_master_strategy;
    in_eps_to_active_constr_to_linearisation = source.in_eps_to_active_constr_to_linearisation;  //relative tolerance to decide if a constraint will be linearized
}


//for lazy constraints, we have to enforce  the mapping to keep original variables in callbacks
int MRQ_LinearApproxAlgorithm:: setLazyConstraintCallbackOnMilpSolver( MRQ_LPSolver *master, void *callbackData)
{
    /***** setting up solver to adopt lazy constraints *****/
    switch(master->getSolverCode())
    {
    case optsolvers::OPT_CPLEX:
        {
        #if OPT_HAVE_CPLEX
            optsolvers::OPT_Cplex *cplex = (OPT_Cplex*) master;
            
            CPXENVptr cplex_env = cplex->env;
            
            
            /* Set up to use MIP lazyconstraint callback. */
            int r = CPXsetlazyconstraintcallbackfunc(cplex_env, MRQ_labbLazyConstraintCplexCallback, callbackData);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            /* Assure linear mappings between the presolved and original
        models (presolver can eliminate variabes...) */
            r = CPXsetintparam(cplex_env, CPXPARAM_Preprocessing_Linear, 0);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            //Turn on traditional search for use with control callbacks
            r = CPXsetintparam(cplex_env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            /*r = CPXsetdblparam(cplex_env, CPXPARAM_Simplex_Tolerances_Feasibility, in_absolute_feasibility_tol);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR); */
            
            /* Let MIP callbacks work on the original model */
            r = CPXsetintparam(cplex_env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
        #endif
            break;
        }
    case optsolvers::OPT_GUROBI:
        {
        #if OPT_HAVE_GUROBI
            OPT_Gurobi *gurobi = (OPT_Gurobi*) master;
            
            //GRBenv *env = gurobi->env;
            GRBmodel *model = gurobi->prob;
            
            int r = GRBsetcallbackfunc(model, MRQ_bboaGurobiCallback, callbackData);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            
            r = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_LAZYCONSTRAINTS, 1);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            
        #endif
            break;
        }
    case optsolvers::OPT_GLPK:
        {
        #if OPT_HAVE_GLPK
            OPT_Glpk *glpk = (OPT_Glpk*) master;
            glp_iocp *parm = glpk->intParam;
            
            parm->cb_func = muriqui::MRQ_bboaGlpkCallback;
            parm->cb_info = callbackData;
            
            parm->binarize = GLP_OFF;
            parm->presolve = GLP_OFF;
            
        #endif
            break;
        }
    default:
        {
            MRQ_PRINTERRORMSG("Invalid solver code");
            return MRQ_NONIMPLEMENTED_ERROR;
        }
    }
    
    return 0;
}




int MRQ_LinearApproxAlgorithm:: setBeforeSolveCallbackOnMilpSolver( MRQ_LPSolver *master, void *callbackData, bool enforceOriginalIndices)
{
    const auto solverCode = master->getSolverCode();
    
    /***** setting up solver to before solve callback *****/
    switch(solverCode)
    {
    case optsolvers::OPT_CPLEX:
        {
        #if OPT_HAVE_CPLEX
            optsolvers::OPT_Cplex *cplex = (OPT_Cplex*) master;
            
            CPXENVptr cplex_env = cplex->env;
            
            
            /* Set up to use MIP solve callback. */
            int r = CPXsetsolvecallbackfunc(cplex_env, MRQ_labbCplexBeforeSolveCallback, callbackData);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            
            if( enforceOriginalIndices )
            {
                /* Assure linear mappings between the presolved and original
            models (presolver can eliminate variabes...) */
                r = CPXsetintparam(cplex_env, CPXPARAM_Preprocessing_Linear, 0);
                MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
                
                /* Let MIP callbacks work on the original model */
                r = CPXsetintparam(cplex_env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
                MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            }
            
            
            //Turn on traditional search for use with control callbacks
            r = CPXsetintparam(cplex_env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            
        #else
            return MRQ_LIBRARY_NOT_AVAILABLE;
        #endif
            break;
        }
    default:
        {
            MRQ_PRINTERRORMSGP("Invalid solver code ", solverCode);
            return MRQ_NONIMPLEMENTED_ERROR;
        }
    }
    
    return 0;
}



int MRQ_LinearApproxAlgorithm:: setBranchingCallbackOnMilpSolver( MRQ_LPSolver *master, void *callbackData, bool enforceOriginalIndices)
{
    const auto solverCode = master->getSolverCode();
    
    /***** setting up solver to before solve callback *****/
    switch(solverCode)
    {
    case optsolvers::OPT_CPLEX:
        {
        #if OPT_HAVE_CPLEX
            optsolvers::OPT_Cplex *cplex = (OPT_Cplex*) master;
            
            CPXENVptr cplex_env = cplex->env;
            
            
            /* Set up to use branching callback. */
            int r = CPXsetbranchcallbackfunc(cplex_env, MRQ_labbCplexBranchingCallback, callbackData);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            if( enforceOriginalIndices )
            {
                /* Assure linear mappings between the presolved and original
        models (presolver can eliminate variabes...) */
                r = CPXsetintparam(cplex_env, CPXPARAM_Preprocessing_Linear, 0);
                MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
                
                /* Let MIP callbacks work on the original model */
                r = CPXsetintparam(cplex_env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
                MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            }
            
            
            //Turn on traditional search for use with control callbacks
            r = CPXsetintparam(cplex_env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            
        #else
            return MRQ_LIBRARY_NOT_AVAILABLE;
        #endif
            break;
        }
    default:
        {
            MRQ_PRINTERRORMSGP("Invalid solver code ", solverCode);
            return MRQ_NONIMPLEMENTED_ERROR;
        }
    }
    
    return 0;
}



int MRQ_LinearApproxAlgorithm:: setDeleteNodeCallbackOnMilpSolver( MRQ_LPSolver *master, void *callbackData, bool enforceOriginalIndices)
{
    const auto solverCode = master->getSolverCode();
    
    /***** setting up solver to before solve callback *****/
    switch(solverCode)
    {
    case optsolvers::OPT_CPLEX:
        {
        #if OPT_HAVE_CPLEX
            optsolvers::OPT_Cplex *cplex = (OPT_Cplex*) master;
            
            CPXENVptr cplex_env = cplex->env;
            
            
            /* Set up to use branching callback. */
            int r = CPXsetdeletenodecallbackfunc( cplex_env, MRQ_labbCplexDeleteNodeCallback, callbackData);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            if( enforceOriginalIndices )
            {
                /* Assure linear mappings between the presolved and original
        models (presolver can eliminate variabes...) */
                r = CPXsetintparam(cplex_env, CPXPARAM_Preprocessing_Linear, 0);
                MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
                
                /* Let MIP callbacks work on the original model */
                r = CPXsetintparam(cplex_env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
                MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            }
            
            
            //Turn on traditional search for use with control callbacks
            r = CPXsetintparam(cplex_env, CPXPARAM_MIP_Strategy_Search, CPX_MIPSEARCH_TRADITIONAL);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            
        #else
            return MRQ_LIBRARY_NOT_AVAILABLE;
        #endif
            break;
        }
    default:
        {
            MRQ_PRINTERRORMSGP("Invalid solver code ", solverCode);
            return MRQ_NONIMPLEMENTED_ERROR;
        }
    }
    
    return 0;
}



int MRQ_LinearApproxAlgorithm::setMasterProblemBase( MRQ_MasterMILPProb *masterMilp, const int thnumber,  MRQ_MINLPProb& prob, const int solver, const bool setLinearObj, const bool setQuad, const bool setVarTypes, const double* lx, const double* ux, const int nauxvars, MRQ_GeneralSolverParams* params, bool *auxConstrEval, MRQ_LAAPointsStoring *laps, const int nthreads  )
{
    MRQ_LPSolver *master;
    
    int ret = masterMilp->setProblemBase(thnumber, prob, solver, setLinearObj, setQuad, setVarTypes, lx, ux, nauxvars, params, nthreads);
    
    if( ret != 0 )
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(ret);
        
        return ret;
    }
    
    master = masterMilp->master;
    
    //we prefer do not set master relative optimality here because only lazy constarints methods like LP-BB, LP/NLP-BB, ESH-BB will use that
    
    if(in_max_cpu_time < INFINITY)
    {
        int r = master->setMaxCPUTime(in_max_cpu_time);
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
            
            return MRQ_MILP_SOLVER_ERROR;
        }
    }
    
    if(in_max_time < INFINITY)
    {
        int r = master->setMaxTime(in_max_time);
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
            
            return MRQ_MILP_SOLVER_ERROR;
        }
    }
    
    if(run_by_inside && !std::isinf(insideSolverMaxTime))
    {
        int r = master->setMaxTime( insideSolverMaxTime );
    
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORMSG("Error to set master probem!");
            
            return MRQ_MILP_SOLVER_ERROR;
        }
    }
    
    if( zu < MRQ_INFINITY )
    {
        const int n = prob.n;
        
        int r = master->setVariableBounds(n, -OPT_INFINITY, zu );
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
            
            return MRQ_MILP_SOLVER_ERROR;
        }
    }
    
    
    if(nPoints > 0)
    {
        ret = masterMilp->addConstraintLinearisationPointsToMILP( in_eps_to_active_constr_to_linearisation, &out_number_of_constr_linears_saved, nPoints, points, setQuad, in_constr_linearisation_strategy, auxConstrEval);
        
        
        ret = masterMilp->addObjLinearisationPointsToMILP( nPoints, points, zu, setQuad, in_obj_linearisation_strategy, NULL, NULL, laps);
        
        if(ret != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(ret);
            
            return MRQ_MILP_SOLVER_ERROR;
        }
        
    }
    
    
    return 0;
}


int MRQ_LinearApproxAlgorithm:: solverCallbackLazyConstraints( MRQ_MILPCallbackData &data)
{
    if(in_print_level > 0)
        MRQ_PRINTERRORMSGP ("Error:  solverCallbackLazyConstraints is not implemented at ", getAlgorithmName() );
    
    return MRQ_NONIMPLEMENTED_ERROR;
}


//I think the behaviour for this method is the same for all linear approximation based branch-and-bound algorithms. So, we implement it here in MRQ_LinearApproxAlgorithm
int MRQ_LinearApproxAlgorithm:: solverCallbackBeforeSolve( MRQ_MILPCallbackData &data)
{
    MRQ_MILPSolverCallbackInterface *milpSolverCallbackInterface = data.callbackSolver;
    
    
    if( in_try_pseudo_pruning_before_solving_relaxations && in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
    {
        const unsigned int thnumber = data.thnumber;
        const int maxBranchIndices = 10;
        
        int r, nBranchIndices;
        int branchIndices[maxBranchIndices];
        MRQ_SizeIndexSol *sizeIndexSol = NULL;
        
        
        
        milpSolverCallbackInterface->getNodeBranchVarIndex(nBranchIndices, branchIndices);
        
        
        #if MRQ_DEBUG_MODE
            assert( nBranchIndices <= maxBranchIndices );
        #endif
        
        
        r = milpSolverCallbackInterface->getUserPointerAtCurrentNode( (void**) &sizeIndexSol);
        MRQ_IFERRORRETURN(r, r);
        
        //in some situatons, solvers do not perform  branching over variables, but in SOS, constraints, etc. So, we must treat these case. In some cases, cplex information diverges from my information... we just try pseudoprunne when both agree. By now, we are only considering the braching CPX_TYPE_VAR where only one branching is done. TODO: if some day we could have branching on more than one variable, fix the if below.
        if(sizeIndexSol && (unsigned int) branchIndices[0] == sizeIndexSol->indexSol[0].index )
        { 
            const unsigned int nBranchIndexSol = sizeIndexSol->size;
            const MRQ_IndexSol *branchIndexSol = sizeIndexSol->indexSol;
            
            double ub = INFINITY;
            double objNode, nodeGap;
            
            //printf("branchIndex: %d\n", branchIndex);
            
            
            #if MRQ_DEBUG_MODE
                if( (unsigned int) branchIndices[0] != branchIndexSol[0].index )
                {
                    printf("branchIndices[0]: %d  branchIndexSol[0].index: %u\n", branchIndices[0], branchIndexSol[0].index);
                }
                assert( (unsigned int) branchIndices[0] == branchIndexSol[0].index ); //this should only works because are threat only branchings under one variable
            #endif
            
            
            //maybe we should put a semaphore here, but i think that too expensive for the benefits
            for(unsigned int i = 0; i < nBranchIndexSol; i++)
                pp_nBranchs[ branchIndexSol[i].index ]++;
            
            
            /* By now, we do not have way to get solution in the parent node to estimate objective increasing. Our alternative would be to store a user handle with the node having this information, but I am not sure about this proccess.*/
              
            r = milpSolverCallbackInterface->getBestSolutionObjValue(ub);
            MRQ_IFERRORRETURN(r, r);
        
            r = milpSolverCallbackInterface->getCurrentNodeLowerBound(objNode);
            MRQ_IFERRORRETURN(r, r);
            
            nodeGap = ub - objNode;
            
            if( nodeGap < in_absolute_convergence_tol_for_pseudo_pruning || nodeGap/MRQ_abs(ub) < in_relative_convergence_tol_for_pseudo_pruning )
            {
                const double ub_slack =  ub + in_relative_upper_bound_slack_factor_for_pseudo_pruning * MRQ_abs(ub) + in_absolute_upper_bound_slack_for_pseudo_pruning;
                
                double deltaObj = 0.0;
                
                for(unsigned int i = 0; i < nBranchIndexSol; i++)
                {
                    const int ind = branchIndexSol[i].index;
                    const double solInd = branchIndexSol[i].sol;
                    
                    double lx, ux;
                    double dowpcost, uppcost, realpcost;
                    double step;
                    
                    
                    r = milpSolverCallbackInterface->getVarBoundsOnNode(ind, lx, ux);
                    MRQ_IFERRORRETURN(r, r);
                    
                    //we only are interested in variables getting fixed in this branch
                    if( lx != ux )
                        continue; 
                    
                    r = milpSolverCallbackInterface->getVarPseudoCosts(ind, dowpcost, uppcost);
                    MRQ_IFERRORRETURN(r, r);
                    
                    
                    if( solInd < lx )
                    { //upper rounding (right branching )
                        step = lx - solInd;
                        realpcost = uppcost;
                    }
                    else
                    { //down rounding (left branching)
                        step = solInd - lx;
                        realpcost = dowpcost;
                    }
                    
                    deltaObj +=  in_alpha_to_balance_estimative_in_pseudo_pruning * realpcost * step;
                }
                
                if(objNode + deltaObj >= ub_slack)
                {
                    //pseudo-prunning!!!!!
                    
                    //printf("pseudo prunning antes de solve (tipo 1)\n");
                    
                    r = milpSolverCallbackInterface->pruneCurrentNode();
                    MRQ_IFERRORRETURN(r, r);
                    
                    if( objNode < pp_lowestBoundPseudoPruned[thnumber] )
                        pp_lowestBoundPseudoPruned[thnumber] = objNode;
                                    
                    pp_counter1[thnumber]++;
                }
                
            }
            
        }
        
        
    }
    
    
    if( !milpSolverCallbackInterface->prunedNode && in_user_callbacks && in_call_before_solve_callback_in_milp_bb)
    {
        const unsigned int thnumber = data.thnumber;
        MRQ_MILPSolverCallbackInterface *milpSolverCallbackInterface = data.callbackSolver;
        MRQ_SolutionStorer solsToBuildLazyConstraints;
        
        int ret = in_user_callbacks->linearApp_beforeSolveInMILPBB( out_algorithm, thnumber, *milpSolverCallbackInterface, solsToBuildLazyConstraints );
        
        if(ret != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTCALLBACKERRORNUMBERWITHNAME(ret, "in_user_callbacks->linearApp_beforeSolveInMILPBB");
            return ret;
        }
        
        //adding lazy constraints 
        if( solsToBuildLazyConstraints.sols.size() > 0 )
        {
            MRQ_MINLPProb &prob = *data.prob;
            
            const int *indices = data.indices;
            const bool *constrEval = data.constrEval;
            
            const bool incQuadsInMaster = data.setQuadsInMaster;
            const bool linearizeObj = data.linearizeObj;
            
            bool *auxConstrEval2 = data.auxConstrEval2;
            
            double *auxVars = data.auxVars;
            double *constrValues = data.auxConstr;
            double *plc = data.plc, *puc = data.puc;
            
            MRQ_GradientsEvaluation &gradEval = data.gradEval;
            MRQ_MILPSolverCallbackInterface *callbackSolver = data.callbackSolver;
            
            MRQ_PointsStoreKeeper *allPointsKeepers = data.allPointsKeepers;
            
            
            for( auto pairSol: solsToBuildLazyConstraints.sols )
            {
                const double objValue = pairSol.first;
                const double *sol = pairSol.second;
                
                ret = prob.constraintsEval(thnumber, true, constrEval, sol, constrValues);
                MRQ_IFERRORRETURN(ret, MRQ_CALLBACK_FUNCTION_ERROR);
                
                ret = addLazyConstraintsLinearizationOnSolution(thnumber, callbackSolver, prob, gradEval, incQuadsInMaster, linearizeObj, constrEval, sol, &objValue, constrValues, constrValues, indices, plc, puc, auxVars, auxConstrEval2);
                MRQ_IFERRORRETURN(ret, ret);
                
                
                //adding solution in allPointsKeepers, for MRQ_LP_BB_ECP_BASED_ALG
                if( allPointsKeepers )
                {
                    const int n = prob.n;
                    
                    for(unsigned int i = 0; i < nthreads_lazy; i++)
                    { //note here, we must add the point in all points keepers...
                        ret = allPointsKeepers[i].insertPoint(n, sol);
                        MRQ_IFERRORRETURN(ret, ret)
                    }
                }
            }
        }
        
    }
    
    
    return 0;
}





int MRQ_LinearApproxAlgorithm:: solverCallbackBranching( MRQ_MILPCallbackData &data)
{
    MRQ_MILPSolverCallbackInterface *milpSolverCallbackInterface = data.callbackSolver;
    
    
    milpSolverCallbackInterface->prunedNode = false;
    
    if( in_user_callbacks && in_call_branching_callback_in_milp_bb )
    {
        const unsigned int thnumber = data.thnumber;
        
        int ret = in_user_callbacks->linearApp_branchingInMILPBB( out_algorithm, thnumber, *milpSolverCallbackInterface);
        
        if(ret != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTCALLBACKERRORNUMBERWITHNAME(ret, "in_user_callbacks->linearApp_beforeSolveInMILPBB");
            return ret;
        }
        
    }
    
    
    /*if the current node was not pruned by used, we proced to pseudo prunning*/
    if( !milpSolverCallbackInterface->prunedNode )
    {
        if( in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
        {
            const unsigned int thnumber = data.thnumber;
            MRQ_MILPSolverCallbackInterface *milpSolverCallbackInterface = data.callbackSolver;
            
            int r;
            double ub = INFINITY;
            double objNode, nodeGap;
            
            
            if( pp_nvars <= 0 )
            {
                //presolve can have eliminated variables. So, we have to ask
                r = milpSolverCallbackInterface->getNumberOfVariables(pp_nvars);
                MRQ_IFERRORRETURN(r, r);
            }
            
            if( pp_nI <= 0 )
            {
                r = milpSolverCallbackInterface->getIntegerVariableIndices(pp_nI, pp_intVars);
                MRQ_IFERRORRETURN(r, r);
            }
            
            
            r = milpSolverCallbackInterface->getBestSolutionObjValue(ub);
            MRQ_IFERRORRETURN(r, r);
            
            r = milpSolverCallbackInterface->getNodeObjValue(objNode);
            MRQ_IFERRORRETURN(r, r);
            
            
            
            //in general, ub should be greater than objNode. However, since we are in the multihtreading environment, ub can be updated after the bound prunning, and so,in this case, objNode will be greater than ub. Due to it, we perform the test if nodeGap is positive.
            
            if( objNode < ub )
            {
                double *sol = NULL;
                /*#if MRQ_DEBUG_MODE
                    assert( ub > objNode || std::isnan(ub) );
                #endif */
                
                nodeGap = ub - objNode;
                
                
                if( nodeGap < in_absolute_convergence_tol_for_pseudo_pruning || nodeGap/MRQ_abs(ub) < in_relative_convergence_tol_for_pseudo_pruning )
                {
                    double *lpcosts = pp_lupcosts[thnumber];
                    double *upcosts = &lpcosts[pp_nvars];
                    
                    
                    const double ub_slack =  ub + in_relative_upper_bound_slack_factor_for_pseudo_pruning * MRQ_abs(ub) + in_absolute_upper_bound_slack_for_pseudo_pruning;
                    
                    
                    
                    //getting the solution
                    sol = pp_sols[thnumber];
                    r = milpSolverCallbackInterface->getNodeSolution(pp_nvars, sol);
                    MRQ_IFERRORRETURN(r, r);
                    
                    
                    //getting the pseudo costs
                    r = milpSolverCallbackInterface->getPseudoCosts(pp_nvars, lpcosts, upcosts);
                    MRQ_IFERRORRETURN(r, r);
                    
                    /*printf("pp_nvars: %d\n", pp_nvars);
                    
                    for(int i = 0; i < pp_nvars; i++)
                    {
                        printf("i: %d\n", i);
                        fflush(stdout);
                        lpcosts[i] = 777;
                        upcosts[i] = 888;
                    } */
                    
                    
                    /*for(decltype(pp_nI) i = 0; i < pp_nI; i++)
                        printf("intVars[%d]: %d\n", i, pp_intVars[i]);
                    
                    for(decltype(pp_nvars) i = 0; i < pp_nvars; i++)
                        printf("i: %d sol: %lf lpc: %lf upc: %lf\n", i, sol[i], lpcosts[i], upcosts[i]);
                    MRQ_getchar(); */
                    
                    
                    for(int i = 0; i < pp_nI; i++)
                    {
                        const int ind = pp_intVars[i];
                        
                        
                        //if( pp_nBranchs[ind]  >= in_min_number_of_brachings_per_var_before_pseudo_pruning ) //TODO: evaluate if we realy must impose it in cplex branching...
                        {
                            double lx, ux;
                            
                            r = milpSolverCallbackInterface->getVarBoundsOnNode(ind, lx, ux);
                            MRQ_IFERRORRETURN(r, r);
                            
                            //we only are interest in non fixed variables that will be fixed in the next branching to try apply the pseudo prunning. If the difference between upper and lower bound is greater than 1.0, variable can be not fixed in the branchings 
                            if( ux - lx == 1.0 )
                            {
                                double step = ceil(sol[ind]) - sol[ind];
                                
                                const double deltaUpObj = in_alpha_to_balance_estimative_in_pseudo_pruning *  upcosts[ind] * step;
                                
                                
                                
                                if( objNode + deltaUpObj >= ub_slack )
                                {
                                    //up rounding (right branch) branching would be pruned by bound accord to pseudo-cost.
                                    step = sol[ind] - floor(sol[ind]);
                                    
                                    const double deltaDownObj = in_alpha_to_balance_estimative_in_pseudo_pruning * lpcosts[ind] * step;
                                    
                                    if(objNode + deltaDownObj >= ub_slack)
                                    {
                                        //pseudo-prunning!!!!!
                                        //printf("pseudo-prunning! var: %d sol: %lf obj: %lf lpc: %lf upc: %lf deltal: %lf deltau: %lf ub_slack: %lf\n", ind, sol[ind], objNode, lpcosts[ind], upcosts[ind], deltaDownObj, deltaUpObj, ub_slack );
                                        //MRQ_getchar();
                                        
                                        r = milpSolverCallbackInterface->pruneCurrentNode();
                                        MRQ_IFERRORRETURN(r, r);
                                        
                                        if( objNode < pp_lowestBoundPseudoPruned[thnumber] )
                                            pp_lowestBoundPseudoPruned[thnumber] = objNode;
                                        
                                        pp_counter2[thnumber]++;
                                        
                                        break;
                                    } //end of if(objNode + deltaDownObj >= ub_slack)
                                } //end of if( objNode + deltaUpObj >= ub_slack )
                            } //end of if( ux - lx == 1.0 )
                        
                        } //end of if( pp_nBranchs[ind]  >= in_min_number_of_brachings_per_var_before_pseudo_pruning )
                        
                    } //end of for(int i = 0; i < pp_nI; i++)
                    
                } //end of  if( nodeGap < in_absolute_convergence_tol_for_pseudo_pruning || nodeGap/MRQ_abs(ub) < in_relative_convergence_tol_for_pseudo_pruning )
                
                
                
                if( in_try_pseudo_pruning_before_solving_relaxations )
                {
                    //performing branching by ourselves to attach information about current relaxation solution at the branching vars.
                    
                    if( !milpSolverCallbackInterface->prunedNode )
                    {
                        if(!sol)
                        {
                            sol = pp_sols[thnumber];
                            //getting the solution
                            r = milpSolverCallbackInterface->getNodeSolution(pp_nvars, sol);
                            MRQ_IFERRORRETURN(r, r);
                        }
                        
                        
                        r = milpSolverCallbackInterface->tryGenerateChildNodesSavingTheParentSolOnBranchVars(sol); //do not return here if you got an error, except MRQ_MEMORY_ERROR
                        if( r == MRQ_MEMORY_ERROR )
                        {
                            MRQ_PRINTMEMERROR;
                            return MRQ_MEMORY_ERROR;
                        }
                    }
                }
                
            } //end of if( objNode < ub )
            
            
        } //end of if( in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
    
    } //end of if( !milpSolverCallbackInterface->prunedNode )
    
    
    
    return 0;
}





MRQ_Heuristic::MRQ_Heuristic():MRQ_Algorithm()
{
}



MRQ_Heuristic::~MRQ_Heuristic()
{
}


void MRQ_Heuristic::copyParametersFrom(const MRQ_Heuristic &source)
{
    MRQ_Algorithm::copyParametersFrom(source);
    
    in_solve_nlp_as_local_search_at_end = source.in_solve_nlp_as_local_search_at_end;
    in_use_random_seed_to_random_numbers = source.in_use_random_seed_to_random_numbers;
    in_seed_to_random_numbers = source.in_seed_to_random_numbers;
    out_seed_to_random_numbers = source.out_seed_to_random_numbers;
}


void MRQ_Heuristic::printParameters(std::ostream &out ) const
{
    MRQ_Algorithm::printParameters(out);
    
    
    out << "\n"
    MRQ_STRFFATT(in_solve_nlp_as_local_search_at_end) << "\n"
    MRQ_STRFFATT(in_use_random_seed_to_random_numbers) << "\n"
    MRQ_STRFFATT(in_seed_to_random_numbers) << "\n"
    ;
}


void MRQ_Heuristic::resetParameters()
{
    MRQ_Algorithm::resetParameters();
    
    in_max_iterations = 1000;
    
    in_solve_nlp_as_local_search_at_end = true;
    in_use_random_seed_to_random_numbers = false;
    
    in_seed_to_random_numbers = 1986;
}



void MRQ_Heuristic::resetOutput()
{
    MRQ_Algorithm::resetOutput();
    
    out_seed_to_random_numbers = in_seed_to_random_numbers;
}



int MRQ_Heuristic::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_Algorithm::setIntegerParameter(name, value);
    
    if(ret == 0)
        return 0;
    
    ret = 0;
    
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_solve_nlp_as_local_search_at_end), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_random_seed_to_random_numbers), name, value ) == 0 );
    else if( MRQ_setAtt<decltype(in_seed_to_random_numbers) >( MRQ_STRATT(in_seed_to_random_numbers), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}




/*
bool MRQ_Heuristic::checkTerminationCriterions(const int threadNumber, const double zl, const double zu, long unsigned int iter, const double timeStart, const clock_t clockStart, int& retCode)
{
    
    const double cpuTime = (double(clock()-clockStart) )/CLOCKS_PER_SEC;
    const double wallTime = MRQ_getTime() - timeStart;
    
    
    
    if(in_user_calbacks)
    {
        if( in_user_calbacks->endOfIteration(threadNumber, iter, cpuTime, wallTime, zl, zu) != 0)
        {
            retCode = MRQ_STOP_REQUIRED_BY_USER;
            return true;
        }
    }
    
    
    
    if( iter >= in_max_iterations )
    {
        retCode = MRQ_MAX_ITERATIONS_STOP;
        return true;
    }
    
    if( cpuTime > in_max_cpu_time )
    {
        retCode = MRQ_MAX_TIME_STOP;
        return true;
    }
    
    if( wallTime > in_max_time )
    {
        retCode = MRQ_MAX_TIME_STOP;
        return true;
    }
    
    
    return false;
}

*/




//to delcare the bool in_use_* references, we use the same order of __getAlgPointer
MRQ_HeuristicExecutor::MRQ_HeuristicExecutor()
{
    __resetMyParameters();
}


bool* MRQ_HeuristicExecutor::getAlgUseFlagPointer( const int number )
{
    return algs[number].first;
}


MRQ_Algorithm* MRQ_HeuristicExecutor::getAlgPointer( const int number )
{
    return algs[number].second;
}


int MRQ_HeuristicExecutor::getNumberOfAlgs() const
{
    return algs.size();
}


void MRQ_HeuristicExecutor::__resetMyParameters()
{
    ai = 0;
    
    in_stop_on_first_sol = true;
    
    in_use_diving = true;
    in_use_fp = true;
    in_use_oafp = true;
    in_use_igma1 = false;
    in_use_igma2 = true;
    in_use_rens = false;
    in_use_ssr = true;
    
    in_max_total_clock_time = INFINITY;
    in_max_total_cpu_time = INFINITY;
    
    
    igma1.in_heuristic_mode = true;
    igma1.in_use_general_feas_heuristics = false;
    igma1.in_max_iterations = 100;
    igma1.in_number_of_threads = 1;
}


bool MRQ_HeuristicExecutor::hasAlgPointer(const MRQ_Algorithm *alg)
{
    const int nAlgs = getNumberOfAlgs();
    
    for(int i = 0; i < nAlgs; i++)
    {
        if(getAlgPointer(i) == alg)
            return true;
    }
    
    return false;
}


void MRQ_HeuristicExecutor::printParameters( std::ostream &out) const
{
    out << MRQ_STRFFATT(in_stop_on_first_sol) << "\n"
    MRQ_STRFFATT(in_use_diving) << "\n"
    MRQ_STRFFATT(in_use_fp) << "\n"
    MRQ_STRFFATT(in_use_oafp) << "\n"
    MRQ_STRFFATT(in_use_igma1) << "\n"
    MRQ_STRFFATT(in_use_igma2) << "\n"
    MRQ_STRFFATT(in_use_rens) << "\n"
    MRQ_STRFFATT(in_max_total_cpu_time) << "\n"
    MRQ_STRFFATT(in_max_total_clock_time) << "\n"
    ;
}



void MRQ_HeuristicExecutor::resetParameters()
{
    const int nAlgs = getNumberOfAlgs();
    
    for(int i = 0; i < nAlgs; i++)
        getAlgPointer(i)->resetParameters();
    
    __resetMyParameters();
}



int MRQ_HeuristicExecutor::setInitialSolution( const int n, const double *sol )
{
    const int nAlgs = getNumberOfAlgs();
    int code = 0;
    
    for(int i = 0; i < nAlgs; i++)
    {
        MRQ_Algorithm *alg = getAlgPointer(i);
        
        const int r = alg->setInitialSolution(n, sol);
        
        if( r != 0 )
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERRORNUMBER(r);
            #endif
            code = r;
        }
    }
    
    return code;
}


void MRQ_HeuristicExecutor::setSolvers( const MRQ_MILP_SOLVER milpSolver,  const MRQ_NLP_SOLVER nlpSolver)
{
    const int nAlgs = getNumberOfAlgs();
    
    for(int i = 0; i < nAlgs; i++)
    {
        MRQ_Algorithm *alg = getAlgPointer(i);
        
        alg->in_milp_solver = milpSolver;
        alg->in_nlp_solver = nlpSolver;
    }
}


//set all heuristics to use the same maximum cpu time
void MRQ_HeuristicExecutor::setMaxCPUTimes( const double maxTime)
{
    const int nAlgs = getNumberOfAlgs();
    
    for(int i = 0; i < nAlgs; i++)
        getAlgPointer(i)->in_max_cpu_time = maxTime;
}



void MRQ_HeuristicExecutor::setMaxTimes( const double maxTime)
{
    const int nAlgs = getNumberOfAlgs();
    
    for(int i = 0; i < nAlgs; i++)
        getAlgPointer(i)->in_max_time = maxTime;
}


void MRQ_HeuristicExecutor::setMaxIters(const long unsigned int maxIters)
{
    const int nAlgs = getNumberOfAlgs();
    
    for(int i = 0; i < nAlgs; i++)
        getAlgPointer(i)->in_max_iterations = maxIters;
}


void MRQ_HeuristicExecutor::setNumberOfThreads(const unsigned int numberOfThreads)
{
    const int nAlgs = getNumberOfAlgs();
    
    for(int i = 0; i < nAlgs; i++)
        getAlgPointer(i)->in_number_of_threads = numberOfThreads;
}


void MRQ_HeuristicExecutor::setPrintLevels(const int level)
{
    const int nAlgs = getNumberOfAlgs();
    
    for(int i = 0; i < nAlgs; i++)
        getAlgPointer(i)->in_print_level = level;
}





int MRQ_HeuristicExecutor::insideRun( MRQ_MINLPProb& prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, double zu, double &obj, double *sol, bool cycleRunning, MRQ_Algorithm **successAlg, const int thnumber, const int nThreads, const double insideSolverMaxTime, double* nlx, double* nux)
{
    return __run(prob, milpSolverParams, nlpSolverParams, zu, obj, sol, cycleRunning, successAlg, true, thnumber, nThreads, insideSolverMaxTime, nlx, nux);
}





int MRQ_HeuristicExecutor::run( MRQ_MINLPProb &prob, MRQ_GeneralSolverParams *milpSolverParams, MRQ_GeneralSolverParams *nlpSolverParams, double zu, double &obj, double *sol, bool cycleRunning, MRQ_Algorithm **successAlg )
{
    //note, arguments after insideRun will not be used
    return __run(prob, milpSolverParams, nlpSolverParams, zu, obj, sol, cycleRunning, successAlg, false, 0, 0, 0.0, NULL, NULL);
}




int MRQ_HeuristicExecutor::__run( MRQ_MINLPProb& prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, double zu, double &obj, double *sol, bool cycleRunning, MRQ_Algorithm **successAlg, const bool insideRun, const int thnumber, const int nThreads, const double insideSolverMaxTime, double* nlx, double* nux)
{
    const double timeStart = MRQ_getTime();
    const clock_t clockStart = clock();
    
    const int nAlgs = getNumberOfAlgs();
    
    int code = MRQ_HEURISTIC_FAIL, r;
    //double time;
    MRQ_Algorithm *alg;
    
    
    
    
    obj = INFINITY;
    
    
    if( !cycleRunning || ai >= nAlgs )
        ai = 0;
    
    
    for( ; ai < nAlgs; ai++)
    {
        if( *getAlgUseFlagPointer(ai) == false )
            continue;
        
        alg = getAlgPointer(ai);
        
        
        
        alg->in_upper_bound = zu;
        alg->in_number_of_threads = nThreads;
        
        if(insideRun)
        {
            //r = alg->insideRun( prob, milpSolverParams, nlpSolverParams, thnumber, nThreads, insideSolverMaxTime, nlx, nux );
            
            r = MRQ_insideRun(alg, prob, milpSolverParams, nlpSolverParams, thnumber, insideSolverMaxTime, nlx, nux);
        }
        else
        {
            r = alg->run( prob, milpSolverParams, nlpSolverParams );
        }
        
        if( r == MRQ_HEURISTIC_SUCCESS || r == MRQ_OPTIMAL_SOLUTION )
        {
            if( alg->out_best_obj <= zu )
            {
                obj = alg->out_best_obj;
                
                if( sol )
                    MRQ_copyArray( prob.n, alg->out_best_sol, sol);
                
                
                zu = obj;
                
                if( successAlg )
                    *successAlg = alg;
                
                if( in_stop_on_first_sol || r == MRQ_OPTIMAL_SOLUTION )
                {
                    //MRQ_getchar();
                    return r;
                }
                else
                    code = MRQ_HEURISTIC_SUCCESS;
            }
        }
        
        
        if( MRQ_getTime() - timeStart >= in_max_total_clock_time  ||  (clock() - clockStart)/CLOCKS_PER_SEC >= in_max_total_cpu_time   )
        {
            if( code != MRQ_HEURISTIC_SUCCESS )
            {
                ai++; //we have to increase ai, otherwise we can always run the same heuristic if we stop by maxtime...
                code = MRQ_MAX_TIME_STOP;
                break;
            }
        }
        
        
    }
    
    
    
    return code;
}









int MRQ_HeuristicExecutor::setIntegerParameter( const char *name, const long int value)
{
    int ret = 0;
    
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_stop_on_first_sol), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_diving), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_fp), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_oafp), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_igma1), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_igma2), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}



int MRQ_HeuristicExecutor::setDoubleParameter( const char *name, const double value)
{
    int ret = 0;
    
    
    if( MRQ_setAtt( MRQ_STRATT(in_max_total_cpu_time), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_max_total_clock_time), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


int MRQ_HeuristicExecutor:: setIntegerParameterToAlgs(const char *name, const long int value)
{
    const int nAlgs = getNumberOfAlgs();
    int r = 0;
    
    for(int i = 0; i < nAlgs; i++)
    {
        r += getAlgPointer(i)->setIntegerParameter(name, value);
    }
    
    return r;
}


int MRQ_HeuristicExecutor:: setDoubleParameterToAlgs(const char *name, const double value)
{
    const int nAlgs = getNumberOfAlgs();
    int r = 0;
    
    for(int i = 0; i < nAlgs; i++)
    {
        r += getAlgPointer(i)->setDoubleParameter(name, value);
    }
    
    return r;
}


int MRQ_HeuristicExecutor::setParametersToAlgs(const MRQ_GeneralSolverParams &params)
{
    const int nAlgs = getNumberOfAlgs();
    int r = 0;
    
    for(int i = 0; i < nAlgs; i++)
    {
        r += getAlgPointer(i)->setParameters(params) ;
    }
    
    return r;
}


int MRQ_HeuristicExecutor:: setStringParameterToAlgs(const char *name, const char * value)
{
    const int nAlgs = getNumberOfAlgs();
    int r = 0;
    
    for(int i = 0; i < nAlgs; i++)
    {
        r += getAlgPointer(i)->setStringParameter(name, value);
    }
    
    return r;
}



MRQ_Algorithm* muriqui::MRQ_newAlgorithm(int algCode, int numberOfNLEqualityConstraints)
{
    MRQ_Algorithm *alg;
    
    
    if(algCode == MRQ_UNDEFINED_ALG)
    { //here, we set the default algorithm. We will set according by instalation characteristcs
        
        if( numberOfNLEqualityConstraints <= 0 )
        {
            if( OPT_isSolverAvailable(OPT_CPLEX) || OPT_isSolverAvailable(OPT_GUROBI) )
            {
                algCode = MRQ_LP_BB_ECP_BASED_ALG;
            }
            else if( MRQ_getDefaultMILPSolverCode() != MRQ_UNDEFINED_MILP ) 
            {
                algCode = MRQ_ECP_ALG;
            }
            else if( MRQ_getDefaultNLPSolverCode() != MRQ_UNDEFINED_NLP )
            {
                //we have a nlp solver, but do not have milp solver. So, B&B
                algCode = MRQ_BB_ALG;
            }
            else
            {
                MRQ_PRINTERRORMSG( "Error: No MILP NOR NLP solvers configured. It is not possible run any MINLP algorithm" );
            }
        }
        else
        {
            //nonlinear equality constraints: we only can apply branch-and-bound
            algCode = MRQ_BB_ALG;
        }
        
    }
    
    
    switch(algCode)
    {
        case MRQ_OA_ALG:
            alg = new (std::nothrow) MRQ_OuterApp;
            break;
        
        case MRQ_LP_NLP_BB_OA_BASED_ALG:
            alg = new (std::nothrow) MRQ_LPNLPBBOuterApp;
            break;
            
        case MRQ_ECP_ALG:
            alg = new (std::nothrow) MRQ_ExtCutPlan;
            break;
        
        case MRQ_BB_ALG:
            alg = new (std::nothrow) MRQ_BranchAndBound;
            break;
            
        case MRQ_IGMA0_ALG:
            alg = new (std::nothrow) MRQ_IGMA0;
            break;
            
        case MRQ_IGMA1_ALG:
            alg = new (std::nothrow) MRQ_IGMA1;
            break;
            
        case MRQ_IGMA2_ALG:
            alg = new (std::nothrow) MRQ_IGMA2;
            break;
            
        case MRQ_FP_HEUR_ALG:
            alg = new (std::nothrow) MRQ_FeasibilityPump;
            break;
            
        case MRQ_DIVE_HEUR_ALG:
            alg = new (std::nothrow) MRQ_Diving;
            break;
            
        case MRQ_OA_FP_HEUR_ALG:
            alg = new (std::nothrow) MRQ_OAFeasibilityPump;
            break;
            
        case MRQ_ESH_ALG: 
            alg = new (std::nothrow) MRQ_ExtSupHypPlane;
            break;
            
        case MRQ_CONT_RELAX_ALG:
            alg = new (std::nothrow) MRQ_ContinuousRelax;
            break;
            
        case MRQ_BONMIN_HYBRID_ALG:
            alg = new (std::nothrow) MRQ_BonminHybrid;
            break;
            
        case MRQ_RENS_HEUR_ALG:
            alg = new (std::nothrow) MRQ_RENS;
            break;
            
        case MRQ_LP_BB_ECP_BASED_ALG:
            alg = new (std::nothrow) MRQ_LPBBExtCutPlan;
            break;
        
        case MRQ_UNDEFINED_ALG: //default algorithm
        case MRQ_LP_BB_ESH_BASED_ALG:
            alg = new (std::nothrow) MRQ_LPBBExtSupHypPlane;
            break;
            
        case MRQ_SSR_HEUR_ALG:
            alg = new (std::nothrow) MRQ_StructuredStochasticRounding;
            break;
            
        default:
            MRQ_PRINTERRORMSGP("Invalid algorithm code: ", algCode);
            alg = NULL;
    }
    
    
    return alg;
}


