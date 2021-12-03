
#include <cmath>
#include <ctime>
#include <climits>

#include <iostream>
#include <new>

#include "BBL_branchAndBound.hpp"
#include "BBL_tools.hpp"

#include "MRQ_solvers.hpp"

#include "MRQ_tools.hpp"
#include "MRQ_bb.hpp"
#include "MRQ_advanced.hpp"
#include "MRQ_algClasses.hpp"



using namespace branchAndBound;
using namespace optsolvers;
using namespace minlpproblem;
using namespace muriqui;





MRQ_BranchAndBound::MRQ_BranchAndBound():MRQ_Algorithm()
{
    resetParameters();
    resetOutput();
    
    out_prune_counter_by_level = NULL;
    
    out_algorithm = MRQ_BB_ALG;
}



MRQ_BranchAndBound::~MRQ_BranchAndBound()
{
    desallocatePruneCountersByLevel();
}





void MRQ_BranchAndBound::resetParameters()
{
    MRQ_Algorithm::resetParameters();

    //in_branch_even_integer_sol = false;
    in_consider_relax_infeas_if_solver_fail = true;
    in_consider_relax_infeas_if_solver_get_unbounded = false;
    in_count_total_prunes = false;

    in_only_apply_pseudo_pruning_on_fixed_integer_vars = true;
    in_repeat_strong_branching_if_get_bounds_updating = true;

    //in_reorganize_lists = true;

    //in_update_var_bounds_using_feas_sols = false;

    in_stop_multibranch_after_first_bound_prune = false;

    in_use_feas_heuristic_diving = false;
    in_use_feas_heuristic_fp = true;
    in_use_feas_heuristic_oafp = true;
    in_use_feas_heuristic_rens = false;

    //in_use_round_heuristic = true;
    in_use_outer_app = true;
    in_use_outer_app_as_heuristic = true;
    in_use_dual_obj_to_bound_prunning = false;
    in_use_early_branching = false;

    
    in_calculate_pseudo_cost_average_above_error_estimative = false;
    in_call_after_bb_loop_callback = false;
    in_call_before_bb_loop_callback = false;
    in_call_before_solving_relax_callback = false;
    in_call_after_solving_relax_callback = false;
    in_call_generate_root_node_callback = false;


    in_set_special_nlp_solver_params = false;

    in_bounds_updt_frequency = 1;
    in_igma2_frequency = 1;
    in_printing_frequency = 1000;

    in_max_number_of_branchings_in_constraint_branching = 5;
    in_max_tree_level_to_count_prunes_by_level = 0;
    in_number_of_node_sublists = 100;
    in_number_of_branching_vars = 1;
    in_feas_heuristic_frequency = 50000;
    in_lists_reorganization_frequency = 10000;

    in_min_number_of_bound_prunes_per_var_before_pseudo_pruning = 10000;

    in_seed_to_random_numbers = 1986;

    in_exp_strategy = MRQ_BB_ES_BEST_LIMIT;
    in_branching_strategy = MRQ_BB_BS_STBRANCH_PSEUDO_COSTS;
    in_variable_choosing_strategy_to_constraint_branching = MRQ_BB_BS_STBRANCH_PSEUDO_COSTS;
    //in_update_var_bounds_strategy = MRQ_BB_UBS_INT_VARS;
    in_parent_node_bounds_storage_strategy = MRQ_BB_PNBSS_SCHAR;
    in_parent_sol_storing_strategy = MRQ_BB_PSSS_NO_STORING;
    in_rounding_heuristic_strategy = MRQ_RS_SSR;
    in_linear_bounds_updt_solver = MRQ_LP_CPLEX;
    in_bounds_linear_updt_strategy = MRQ_BB_BLUS_NO_UPDT;
    in_constr_branching_strategy = MRQ_BB_CBS_NO_CONSTRAINT_BRANCH;
    in_igma2_strategy = MRQ_BB_IHS_UNTIL_FIRST_FEAS_SOL;
    in_int_feas_heurs_strategy = MRQ_BB_IHS_UNTIL_FIRST_FEAS_SOL;

    in_pseudo_pruning_strategy = MRQ_BB_PPS_NO_PSEUDO_PRUNING;

    //in_number_of_pseudo_costs_rely = 1;

    in_outer_app_frequence = 10000;
    in_outer_app_subprob_frequence = 100;
    in_rounding_heuristic_call_iter_frequence = 1000;

    in_number_of_pos_iters_to_update_var_bounds = 10000;


    in_feas_heuristic_max_time = 30;
    in_bounds_updt_factor_to_subgroup = 0.1;
    in_outer_app_time = 10;
    in_outer_app_subprob_time = 10;
    in_pseudo_cost_mu = 0.2;
    in_pseudo_cost_mu_to_variable_choosing_on_constraint_branching = 0.2;

    in_nlp_relative_opt_tol_on_early_branching = 1.0e-3;	in_nlp_relative_opt_tol_to_int_sol_on_early_branching = 1.0e-6;


    in_nlp_relative_primal_tol_on_early_branching = 1.0e-4;
    in_nlp_relative_primal_tol_to_int_sol_on_early_branching = 1.0e-6;

    in_nlp_relative_dual_tol_on_early_branching = 1.0e-4;
    in_nlp_relative_dual_tol_to_int_sol_on_early_branching = 1.0e-6;


    in_absolute_convergence_tol_for_pseudo_pruning = 1.0;
    in_relative_convergence_tol_for_pseudo_pruning = 0.1;

    in_absolute_upper_bound_slack_for_pseudo_pruning = 1e-3;
    in_relative_upper_bound_slack_factor_for_pseudo_pruning = 0.01;

    in_alpha_to_balance_estimative_in_pseudo_pruning = 0.4;
    
    in_integer_neighborhood_factor_to_ssr = 0.3;

    /********* begin of igma2 parameters **********/

    in_igma2_set_max_dist_constr_on_bin_vars = false;
    in_igma2_set_special_gap_min_solver_params = true;
    in_igma2_solve_local_search_problem_even_on_non_int_sol = true;
    //in_igma2_gap_min_solver = MRQ_getDefaultMinGapNLPSolverCode();
    in_igma2_factor_to_max_dist_constr = 0.05;
    in_igma2_factor_to_max_dist_constr_on_bin_vars = 0.05;
    in_igma2_percentual_to_rectangular_neighborhood = 0.05;

    in_igma2_neighborhood_strategy = MRQ_IGMA2_NS_SPHERIC;
    /********* end of igma2 parameters **********/




    in_OAmilpParams = NULL;
    in_OAnlpParams = NULL;

    //in_outer_app = NULL;
    //in_outer_app_sub = NULL;
}



void MRQ_BranchAndBound::resetOutput()
{
    MRQ_Algorithm::resetOutput();
    
    out_nlp_failure_in_some_relaxation = false;
    out_seed_to_random_numbers = in_seed_to_random_numbers;
    out_number_of_open_nodes = 0;
    out_number_of_iters_having_wrong_lower_bounds = 0;
    out_number_of_constraint_branchings = 0;
    out_number_of_strong_branching_calculations_to_pseudo_costs = 0;
    out_pseudo_cost_average_above_error_estimative = 0.0;
    
    out_prune_counter.reset();
}



void MRQ_BranchAndBound::desallocatePruneCountersByLevel()
{
    if(out_prune_counter_by_level)
    {
        MRQ_secDeleteArray(out_prune_counter_by_level);
    }
}



void MRQ_BranchAndBound::printParameters(std::ostream& out ) const
{
    char strValue[100];
    
    
    MRQ_Algorithm::printParameters(out);
    out << "\n"
    MRQ_STRFFATT(in_calculate_pseudo_cost_average_above_error_estimative) << "\n"
    MRQ_STRFFATT(in_call_after_bb_loop_callback) << "\n"
    MRQ_STRFFATT(in_call_before_bb_loop_callback) << "\n"
    MRQ_STRFFATT(in_call_before_solving_relax_callback) << "\n"
    MRQ_STRFFATT(in_call_after_solving_relax_callback) << "\n"
    MRQ_STRFFATT(in_call_generate_root_node_callback) << "\n"
    MRQ_STRFFATT(in_count_total_prunes) << "\n"
    MRQ_STRFFATT(in_consider_relax_infeas_if_solver_fail) << "\n"
    MRQ_STRFFATT(in_consider_relax_infeas_if_solver_get_unbounded) << "\n"
    MRQ_STRFFATT( in_only_apply_pseudo_pruning_on_fixed_integer_vars ) << "\n"
    MRQ_STRFFATT( in_repeat_strong_branching_if_get_bounds_updating ) << "\n"
    MRQ_STRFFATT(in_stop_multibranch_after_first_bound_prune) << "\n"
    MRQ_STRFFATT(in_use_dual_obj_to_bound_prunning) << "\n"
    MRQ_STRFFATT(in_use_early_branching) << "\n"
    
    MRQ_STRFFATT(in_use_feas_heuristic_diving) << "\n"
    MRQ_STRFFATT(in_use_feas_heuristic_fp) << "\n"
    MRQ_STRFFATT(in_use_feas_heuristic_oafp) << "\n"
    MRQ_STRFFATT(in_use_feas_heuristic_rens) << "\n"
    
    MRQ_STRFFATT(in_use_outer_app) << "\n"
    MRQ_STRFFATT(in_use_outer_app_as_heuristic) << "\n"
    
    MRQ_STRFFATT(in_bounds_updt_frequency) << "\n"
    MRQ_STRFFATT(in_feas_heuristic_frequency) << "\n"
    MRQ_STRFFATT(in_igma2_frequency) << "\n"
    MRQ_STRFFATT(in_lists_reorganization_frequency) << "\n"	
    MRQ_STRFFATT(in_max_number_of_branchings_in_constraint_branching) << "\n"
    MRQ_STRFFATT(in_max_tree_level_to_count_prunes_by_level) << "\n"
    MRQ_STRFFATT(in_number_of_branching_vars) << "\n"
    MRQ_STRFFATT(in_number_of_pos_iters_to_update_var_bounds) << "\n"
    //MRQ_STRFFATT(in_number_of_pseudo_costs_rely) << "\n"
    MRQ_STRFFATT(in_number_of_node_sublists) << "\n"
    MRQ_STRFFATT(in_outer_app_frequence) << "\n"
    MRQ_STRFFATT(in_outer_app_subprob_frequence) << "\n"
    MRQ_STRFFATT(in_rounding_heuristic_call_iter_frequence) << "\n"
    MRQ_STRFFATT(in_seed_to_random_numbers) << "\n"
    MRQ_STRFFATT(in_min_number_of_bound_prunes_per_var_before_pseudo_pruning) << "\n";
    
    
    MRQ_enumToStr(in_linear_bounds_updt_solver, strValue);
    out << MRQ_STRPARINTVALUE(in_linear_bounds_updt_solver) " " << strValue << "\n";
    
    MRQ_enumToStr(in_bounds_linear_updt_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_bounds_linear_updt_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr(in_exp_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_exp_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr(in_branching_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_branching_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr(in_variable_choosing_strategy_to_constraint_branching, strValue);
    out << MRQ_STRPARINTVALUE(in_variable_choosing_strategy_to_constraint_branching) " " << strValue << "\n";
    
    MRQ_enumToStr(in_parent_node_bounds_storage_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_parent_node_bounds_storage_strategy) " " << strValue << "\n";
    MRQ_enumToStr(in_parent_sol_storing_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_parent_sol_storing_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr(in_rounding_heuristic_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_rounding_heuristic_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr(in_constr_branching_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_constr_branching_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr(in_igma2_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_igma2_strategy) " " << strValue << "\n";
    
    
    MRQ_enumToStr(in_int_feas_heurs_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_int_feas_heurs_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr(in_pseudo_pruning_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_pseudo_pruning_strategy) " " << strValue << "\n";
    
    out << 
    MRQ_STRFFATT(in_feas_heuristic_max_time) << "\n"
    MRQ_STRFFATT(in_bounds_updt_factor_to_subgroup) << "\n"
    MRQ_STRFFATT(in_outer_app_subprob_time) << "\n"
    MRQ_STRFFATT(in_outer_app_time) << "\n"
    MRQ_STRFFATT(in_pseudo_cost_mu) << "\n"
    MRQ_STRFFATT(in_pseudo_cost_mu_to_variable_choosing_on_constraint_branching) << "\n"
    MRQ_STRFFATT(in_nlp_relative_opt_tol_on_early_branching) << "\n"
    MRQ_STRFFATT(in_nlp_relative_opt_tol_to_int_sol_on_early_branching) << "\n"
    MRQ_STRFFATT(in_nlp_relative_primal_tol_on_early_branching) << "\n"
    MRQ_STRFFATT(in_nlp_relative_primal_tol_to_int_sol_on_early_branching) << "\n"
    MRQ_STRFFATT(in_nlp_relative_dual_tol_on_early_branching) << "\n"
    MRQ_STRFFATT(in_nlp_relative_dual_tol_to_int_sol_on_early_branching) << "\n"
    MRQ_STRFFATT(in_absolute_convergence_tol_for_pseudo_pruning) << "\n"
    MRQ_STRFFATT(in_relative_convergence_tol_for_pseudo_pruning) << "\n"
    MRQ_STRFFATT(in_absolute_upper_bound_slack_for_pseudo_pruning) << "\n"
    MRQ_STRFFATT(in_relative_upper_bound_slack_factor_for_pseudo_pruning) << "\n"
    MRQ_STRFFATT(in_alpha_to_balance_estimative_in_pseudo_pruning) << "\n"
    MRQ_STRFFATT(in_integer_neighborhood_factor_to_ssr) << "\n"
    
    /********* begin of igma2 parameters **********/
    
    MRQ_STRFFATT(in_igma2_set_max_dist_constr_on_bin_vars) << "\n"
    MRQ_STRFFATT(in_igma2_set_special_gap_min_solver_params) << "\n"
    MRQ_STRFFATT(in_igma2_solve_local_search_problem_even_on_non_int_sol) << "\n";
    
    /*MRQ_enumToStr(in_igma2_gap_min_solver, strValue);
    out << MRQ_STRPARINTVALUE(in_igma2_gap_min_solver) ": " << strValue << "\n";*/
    
    MRQ_enumToStr(in_igma2_neighborhood_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_igma2_neighborhood_strategy) ": " << strValue << "\n"
    
    MRQ_STRFFATT(in_igma2_factor_to_max_dist_constr) << "\n"
    MRQ_STRFFATT(in_igma2_factor_to_max_dist_constr_on_bin_vars) << "\n"
    MRQ_STRFFATT(in_igma2_percentual_to_rectangular_neighborhood) << "\n"
    
    /********* end of igma3 parameters **********/
    
    ;
}



int MRQ_BranchAndBound::setDCS0Array( const int size, const MRQ_DynConstrSetUnity *dcs)
{
    if(ndcs0 > 0)
        desallocateDCS0();
    
    const int r = allocateDCS0(size);
    if( r != 0 )
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    MRQ_copyArray(size, dcs, dcs0);
    
    ndcs0 = size;
    return 0;
}


int MRQ_BranchAndBound::setDCS1Array( const int size, const MRQ_DynConstrSetUnity *dcs)
{
    if(ndcs1 > 0)
        desallocateDCS1();
    
    const int r = allocateDCS1(size);
    if( r != 0 )
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    MRQ_copyArray(size, dcs, dcs1);
    
    ndcs1 = size;
    return 0;
}


int MRQ_BranchAndBound::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_Algorithm::setIntegerParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_calculate_pseudo_cost_average_above_error_estimative), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_call_after_bb_loop_callback), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_call_before_bb_loop_callback), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_call_before_solving_relax_callback ), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_call_after_solving_relax_callback), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_call_generate_root_node_callback), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_count_total_prunes), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_consider_relax_infeas_if_solver_fail), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_consider_relax_infeas_if_solver_get_unbounded), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT( in_only_apply_pseudo_pruning_on_fixed_integer_vars ), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT( in_repeat_strong_branching_if_get_bounds_updating ), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_stop_multibranch_after_first_bound_prune), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_dual_obj_to_bound_prunning), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_early_branching), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_feas_heuristic_diving), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_feas_heuristic_fp), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_feas_heuristic_oafp), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_feas_heuristic_rens), name, value ) == 0 );
    
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_outer_app), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_outer_app_as_heuristic), name, value ) == 0 );
    
    
    else if( MRQ_setAtt<int>( MRQ_STRATT(in_bounds_updt_frequency), name, value) == 0 );
    else if( MRQ_setAtt<int>( MRQ_STRATT(in_feas_heuristic_frequency), name, value) == 0 );
    else if( MRQ_setAtt<int>( MRQ_STRATT(in_igma2_frequency), name, value) == 0 );
    else if( MRQ_setAtt<int>( MRQ_STRATT(in_lists_reorganization_frequency), name, value) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_max_number_of_branchings_in_constraint_branching), name, value) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_max_tree_level_to_count_prunes_by_level), name, value) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_min_number_of_bound_prunes_per_var_before_pseudo_pruning), name, value) == 0 );
    else if( MRQ_setAtt<int>( MRQ_STRATT(in_number_of_branching_vars), name, value) == 0 );
    else if( MRQ_setAtt<int>( MRQ_STRATT(in_number_of_node_sublists), name, value) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_number_of_pos_iters_to_update_var_bounds), name, value) == 0 );
    else if( MRQ_setAtt<int>( MRQ_STRATT(in_outer_app_frequence), name, value) == 0 );
    else if( MRQ_setAtt<int>( MRQ_STRATT(in_outer_app_subprob_frequence), name, value) == 0 );
    else if( MRQ_setAtt<int>( MRQ_STRATT(in_rounding_heuristic_call_iter_frequence), name, value) == 0);
    else if( MRQ_setAtt< decltype(in_seed_to_random_numbers) >( MRQ_STRATT(in_seed_to_random_numbers), name, value) == 0);
    
    /********* begin of igma 2 parameters **********/
    
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_igma2_set_max_dist_constr_on_bin_vars), name, value) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_igma2_set_special_gap_min_solver_params), name, value) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_igma2_solve_local_search_problem_even_on_non_int_sol), name, value) == 0);
    
    /********* end of igma3 parameters **********/
    
    else 
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}





int MRQ_BranchAndBound::setDoubleParameter(const char *name, const double value)
{
    int ret = MRQ_Algorithm::setDoubleParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt( MRQ_STRATT(in_feas_heuristic_max_time), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_bounds_updt_factor_to_subgroup), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_nlp_relative_opt_tol_on_early_branching), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_nlp_relative_opt_tol_to_int_sol_on_early_branching), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_nlp_relative_primal_tol_on_early_branching), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_nlp_relative_primal_tol_to_int_sol_on_early_branching), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_nlp_relative_dual_tol_on_early_branching), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_nlp_relative_dual_tol_to_int_sol_on_early_branching), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_absolute_convergence_tol_for_pseudo_pruning), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_relative_convergence_tol_for_pseudo_pruning), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_absolute_upper_bound_slack_for_pseudo_pruning), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_relative_upper_bound_slack_factor_for_pseudo_pruning), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_alpha_to_balance_estimative_in_pseudo_pruning), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_outer_app_subprob_time), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_outer_app_time), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_pseudo_cost_mu), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_pseudo_cost_mu_to_variable_choosing_on_constraint_branching), name, value ) == 0 );
    
    /********* begin of igma2 parameters **********/
    
    else if( MRQ_setAtt( MRQ_STRATT(in_igma2_factor_to_max_dist_constr), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_igma2_factor_to_max_dist_constr_on_bin_vars), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_igma2_percentual_to_rectangular_neighborhood), name, value ) == 0 );
    
    /********* end of igma3 parameters **********/
    
    else
        ret = MRQ_NAME_ERROR;
    
    
    
    return ret;
}




int MRQ_BranchAndBound::setStringParameter(const char *name, const char *value)
{
    int r;
    int ret = MRQ_Algorithm::setStringParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    //(r = MRQ_setStrAtt( MRQ_STRATT(in_milp_solver), name, value ) ) >= 0
    
    
    if( (r = MRQ_setStrAtt( MRQ_STRATT(in_bounds_linear_updt_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_branching_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_constr_branching_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_exp_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_igma2_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_int_feas_heurs_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_linear_bounds_updt_solver), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_parent_node_bounds_storage_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_parent_sol_storing_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_rounding_heuristic_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_variable_choosing_strategy_to_constraint_branching), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_pseudo_pruning_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    
    
    /********* begin of igma3 parameters **********/
    
    /*else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_igma2_gap_min_solver), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR; */
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_igma2_neighborhood_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    
    /********* end of igma3 parameters **********/
    
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}



bool MRQ_BranchAndBound::tryUpdateBestSolution(const unsigned int threadNumber, const int n, double* sol, double objValue, const long unsigned int iter, const clock_t& clockStart, const double timeStart, const bool storeOnHistory)
{
    //assert(false);
    
    return tryUpdateBestSolution(threadNumber, n, sol, objValue, iter, clockStart, timeStart, storeOnHistory, true);
}



bool MRQ_BranchAndBound::tryUpdateBestSolution(const unsigned int threadNumber, const int n, double* sol, double objValue, const long unsigned int iter, const clock_t& clockStart, const double timeStart, const bool storeOnHistory, const bool addSolToOALin)
{
    return _bblCallbacks->tryUpdateBestSolution( threadNumber, sol,  objValue, iter, addSolToOALin );
    //assert(false);
}



int MRQ_BranchAndBound::getBestSolutionCopy(const int n, double* solution, double& fsolution)
{
    assert(false);
    return 0;
}



int MRQ_BranchAndBound:: allocatePruneCountersByLevel(const int nlevels)
{
    desallocatePruneCountersByLevel();
    
    out_prune_counter_by_level = new (std::nothrow) MRQ_PruneCounter[nlevels];
    
    if(!out_prune_counter_by_level)
    {
        if(in_print_level > 0)
            MRQ_PRINTMEMERROR;
        
        return MRQ_MEMORY_ERROR;
    }
    
    return 0;
}




int MRQ_BranchAndBound::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    const int n = prob.n;
    const int m = prob.m;
    const int ndual = 2*n + m;


    int ret;
    int unsigned nthreads;
    double *plc = NULL, *puc = NULL; // do not free puc, because we take advatnage the malloc of plc. We just put NULL because we should sinalize if there is new bounds or no to other procedures like MRQ_BinSumConstrs::calculateIndices when we do not preprocess...


    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux; 

    BBL_BranchAndBound bb;
    MRQ_Preprocessor *preprocessor = nullptr;
    MIP_BinSumConstrsIndsByClass *ssrBinSumConstrs = nullptr;
    MRQ_BBLCallbacks bbcallbacks(&prob, this, milpSolverParams, nlpSolverParams);

    double timeStart = MRQ_getTime();
    clock_t clockStart = clock();

    BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY parentNodeBoundsStrategy;

    _bblCallbacks = &bbcallbacks;


    nthreads = in_number_of_threads <= 0 ? MRQ_getNumCores() : in_number_of_threads;


    if( in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_obj_function )
    {
        preprocessor = new (std::nothrow) MRQ_Preprocessor(&prob);
        MRQ_IFMEMERRORGOTOLABEL(!preprocessor, out_return_code, termination);
    }


    {
        auto ret = algorithmInitialization( nthreads, true, milpSolverParams, nlpSolverParams, prob, lx, ux, preprocessor, NULL, &plc, &puc, false );
        
        if( ret != 0 )
        {
            if( in_print_level > 0 )
            {
                if( ret == MRQ_ALG_NOT_APPLICABLE && (out_algorithm == MRQ_IGMA1_ALG || out_algorithm == MRQ_IGMA2_ALG) )
                    MRQ_PRINTERRORMSG("Error: Integrality Gap Minimization Algorithm only hands binary problems");
                else
                    MRQ_PRINTERRORNUMBER(ret);
            }
            out_return_code = ret;
            goto termination;
        }
        
        #if 0
        {
            const int m = prob.m;
            const int n = prob.n;
            const double *lc = prob.lc, *uc = prob.uc;
            const double *olx = prob.lx, *oux = prob.ux; 
            for(int i = 0; i < n; i++)
                printf("lx[%d]: %lf  ux[%d]: %lf  nlx[%d]: %lf nux[%d]: %lf\n", i, olx[i], i, oux[i], i, lx[i], i, ux[i]);
            
            for(int i = 0; i < m; i++)
                printf("lc[%d]: %lf  uc[%d]: %lf   plc[%d]: %lf  puc[%d]: %lf\n", i, lc[i], i, lc[i], i, plc[i], i, puc[i]);
        }
        #endif
    }

    MRQ_secDelete(preprocessor);


    //chekcing dynamic constraint set indices
    if( in_use_dynamic_constraint_set )
    {
        int r;
        MRQ_DynConstrSetSetter dcssetter;
        
        r = dcssetter.checkIndices(n, m, ndcs0, dcs0);
        if( r != 0 )
        {
            if( in_print_level > 0 )
            {
                MRQ_PRINTERRORMSG("Error on dcs0 indices\n");
                MRQ_PRINTERRORNUMBER(r);
            }
            out_return_code = MRQ_VALUE_ERROR;
            goto termination;
        }
        
        r = dcssetter.checkIndices(n, m, ndcs1, dcs1);
        if( r != 0 )
        {
            if( in_print_level > 0 )
            {
                MRQ_PRINTERRORMSG("Error on dcs1 indices\n");
                MRQ_PRINTERRORNUMBER(r);
            }
            out_return_code = MRQ_VALUE_ERROR;
            goto termination;
        }
    }


    if( in_print_level > 1 )
        std::cout << "\nStarting Branch-and-Bound algorithm\n\n";

    if(in_print_level > 2)
        printSubSolvers(true, true, false);



    if(in_max_tree_level_to_count_prunes_by_level > 0)
    {
        ret = allocatePruneCountersByLevel( in_max_tree_level_to_count_prunes_by_level );
        
        if( ret != 0 )
        {
            MRQ_PRINTERRORNUMBER(ret);
            goto termination;
        }
    }

    if( in_user_callbacks)
    {
        if( in_user_callbacks->alg == this )
            in_user_callbacks->bblCallBacks = &bbcallbacks;
    }


    bbcallbacks.oplc = plc;

    {
        int nI, *intVars;
        
        MRQ_malloc(intVars, n);
        MRQ_IFMEMERRORGOTOLABEL(!intVars, out_return_code, termination);
        
        nI = prob.getIntegerIndices(intVars);
        
        parentNodeBoundsStrategy = (BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY) in_parent_node_bounds_storage_strategy;
        
        MRQ_checkParentNodeBoundsStrategy(nI, intVars, lx, ux, parentNodeBoundsStrategy);
        
        if( in_rounding_heuristic_strategy == MRQ_RS_SSR )
        {
            ssrBinSumConstrs = new (std::nothrow) MIP_BinSumConstrsIndsByClass;
            
            int r = ssrBinSumConstrs->calculateIndices(prob, lx, ux, plc, puc, nullptr, true, true, true);
            MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_UNDEFINED_ERROR, termination);
            
            bbcallbacks.ssrBinSumConstrs = ssrBinSumConstrs;
        }
        
        free(intVars);
    }
    
    
    

    bbcallbacks.parentNodeBoundsStrategy = parentNodeBoundsStrategy;

    bb.in_call_after_bb_loop_callback = in_call_after_bb_loop_callback;
    bb.in_call_before_bb_loop_callback = in_call_before_bb_loop_callback;
    bb.in_call_before_solving_relax_callback = false;
    bb.in_call_updating_best_solution_callback = in_call_update_best_sol_callback;
    bb.in_call_new_best_solution_callback = in_call_new_best_solution_callback;

    bb.in_call_end_of_iteration_callback =  in_call_end_of_iteration_callback;

    bb.in_count_total_prunes = in_count_total_prunes;

    bb.in_branching_strategy = in_branching_strategy == MRQ_BB_BS_USER_NODE_GENERATION ? BBL_BS_USER_NODE_GENERATION : BBL_BS_USER_INDEX_CHOICE;

    bbcallbacks.origBBLBranchStrategy = bb.in_branching_strategy;

    bb.in_consider_relax_infeas_if_solver_fail = in_consider_relax_infeas_if_solver_fail;
    bb.in_use_dual_obj_to_bound_prunning = in_use_dual_obj_to_bound_prunning;
    bb.in_reorganize_lists = in_lists_reorganization_frequency > 0;
    bb.in_max_tree_level_to_count_prunes_by_level = in_max_tree_level_to_count_prunes_by_level;

    bb.in_number_of_node_sublists = in_number_of_node_sublists;
    bb.in_number_of_threads = in_number_of_threads;
    bb.in_infinity_value = MRQ_INFINITY;

    bb.in_exp_strategy = MRQ_MRQ_exp_strategy2BBL_exp_strategy(in_exp_strategy);


    bb.in_store_parent_primal_solution_on_nodes = in_parent_sol_storing_strategy == MRQ_BB_PSSS_ONLY_PRIMAL || in_parent_sol_storing_strategy == MRQ_BB_PSSS_PRIMAL_AND_DUAL;
    bb.in_store_parent_dual_solution_on_nodes = in_parent_sol_storing_strategy == MRQ_BB_PSSS_PRIMAL_AND_DUAL;

    bb.in_parent_node_bounds_strategy = parentNodeBoundsStrategy;

    bb.in_store_history_solutions = in_store_history_solutions;

    bb.in_max_iterations = in_max_iterations;
    bb.in_printing_frequency = in_printing_frequency;
    bb.in_print_level = in_print_level;

    bb.in_absolute_convergence_tol = in_absolute_convergence_tol;
    bb.in_lower_bound = in_lower_bound;
    bb.in_max_time = in_max_time;
    bb.in_max_cpu_time = in_max_cpu_time;
    bb.in_relative_convergence_tol = in_relative_convergence_tol;
    bb.in_upper_bound = in_upper_bound;


    bb.run( n, m, ndual, lx, ux, bbcallbacks );


    out_number_of_open_nodes = bb.out_number_of_open_nodes;

    out_number_of_iterations = bb.out_number_of_iterations;

    out_number_of_threads = bb.out_number_of_threads;

    if( bb.out_best_obj < zu )
    {
        out_best_obj = zu = bb.out_best_obj;
        MRQ_copyArray(n, bb.out_best_sol, out_best_sol);
    }

    out_number_of_iterations_to_first_feas_sol = bb.out_first_sol_iter;
    out_number_of_iterations_to_best_sol = bb.out_best_sol_iter;
    out_cpu_time_to_first_feas_sol = bb.out_cpu_time_to_fisrt_sol;
    out_clock_time_to_first_feas_sol = bb.out_clock_time_to_fisrt_sol;
    out_cpu_time_to_best_sol = bb.out_cpu_time_to_best_sol;
    out_clock_time_to_best_sol = bb.out_clock_time_to_best_sol;


    out_number_of_feas_sols = bb.out_number_of_feas_sols;

    out_lower_bound = bb.out_lower_bound;

    if( in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
    {
        bool updt = false;
        
        //we have to correct lower bound output
        
        auto ppruningLowestlb = bbcallbacks.pseudoPruningLowestLowerBound;
        
        for(unsigned int i = 0; i < nthreads; i++)
        {
            if( ppruningLowestlb[i] < out_lower_bound )
            {
                out_lower_bound = ppruningLowestlb[i];
                updt = true;
            }
        }
        
        if(in_print_level > 4 && updt)
            MRQ_PRINTERRORMSGP("Correcting lower bound due to pseudo pruning to ", out_lower_bound);
    }
    
    if( in_calculate_pseudo_cost_average_above_error_estimative && in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
    {
        auto pCostStdDeviationEstimative = bbcallbacks.pCostAvgAboveErrorEstimative;
        auto nPCostStdDeviationEstimative = bbcallbacks.nPCostAvgAboveErrorEstimative;
        
        long unsigned int nEstimatives = 0;
        
        out_pseudo_cost_average_above_error_estimative = 0;
        
        for(decltype(nthreads) i = 0; i < nthreads; i++)
            nEstimatives += nPCostStdDeviationEstimative[i];
        
        if( nEstimatives > 0 )
        {
            for(decltype(nthreads) i = 0; i < nthreads; i++ )
                out_pseudo_cost_average_above_error_estimative += pCostStdDeviationEstimative[i];
            
            out_pseudo_cost_average_above_error_estimative /= nEstimatives;
        }
    }
    

    zl = out_lower_bound; //let zl be updated here since algorithmFinalization will use this value
    out_upper_bound = bb.out_upper_bound;
    out_obj_opt_at_continuous_relax = bb.out_obj_opt_at_root_relax;

    {
        auto numberOfWrongLowerBounds = bbcallbacks.numberOfWrongLowerBounds;
        for(unsigned int i = 0; i < nthreads; i++)
        {
            out_number_of_iters_having_wrong_lower_bounds += numberOfWrongLowerBounds[i];
        }
        
        if( bbcallbacks.numberOfConstrBranchings )
        {
            auto numberOfConstrBranchings = bbcallbacks.numberOfConstrBranchings;
            
            for(unsigned int i = 0; i < nthreads; i++)
                out_number_of_constraint_branchings += numberOfConstrBranchings[i];
        }
    }


    switch( bb.out_return_code )
    {
        case BBL_OPTIMAL_SOLUTION:
            out_return_code = MRQ_OPTIMAL_SOLUTION;
            break;
            
        case BBL_MAX_TIME_STOP:
            out_return_code = MRQ_MAX_TIME_STOP;
            break;
            
        case BBL_MAX_ITERATIONS_STOP:
            out_return_code = MRQ_MAX_ITERATIONS_STOP;
            break;
            
        case BBL_STOP_REQUIRED_BY_USER:
            
            if( bbcallbacks.userErrorCode == 0 )
                out_return_code = MRQ_intToReturnCode(bb.out_return_subcode);
            else
            {
                out_return_code = MRQ_STOP_REQUIRED_BY_USER;
                out_user_callback_error_code = bbcallbacks.userErrorCode;
            }
            
            break;
            
        case BBL_INFEASIBLE_PROBLEM:
        case BBL_NO_FEASIBLE_SOLUTION_FOUND:
            out_return_code = MRQ_INFEASIBLE_PROBLEM;
            break;
            
        case BBL_MEMORY_ERROR:
            out_return_code = MRQ_MEMORY_ERROR;
            break;
            
        case BBL_SOLVER_ERROR:
            out_return_code = MRQ_NLP_SOLVER_ERROR;
            break;
            
        //case BBL_FEASIBLE_SOLUTION:
        default:
            out_return_code = MRQ_UNDEFINED_ERROR;
            break;
    }


    if( in_store_history_solutions )
    {
        BBL_SolutionHistory &shist = bb.out_sol_hist;
        const int nsols = shist.getnsols();
        
        const double pretime = 0;// timePrebb -timeStart;
        const double precputime = 0;//(clockPrebb - clockStart)/(double) CLOCKS_PER_SEC;
        
        for(int i = 0; i < nsols; i++)
        {
            BBL_HistorySolution *hsol = shist.getHistSolPointer(i);
            
            out_sol_hist.addSolution( n, hsol->iter, hsol->time + pretime, hsol->cputime + precputime, hsol->sol, hsol->objF );
        }
        
    }


    if( in_count_total_prunes )
    {
        out_prune_counter.accumulate( bb.out_prune_counter );
        
        auto pseudoPruningCounter1 = bbcallbacks.beforeExpoPseudoPruningCounter;
        
        auto pseudoPruningCounter2 = bbcallbacks.afterExpoPseudoPruningCounter;
        
        if(in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
        {
            for(unsigned int i = 0; i < nthreads; i++)
            {
                out_prune_counter.pseudopruning1 += pseudoPruningCounter1[i];
                out_prune_counter.pseudopruning2 += pseudoPruningCounter2[i];
            }
        }
    }


    if( in_max_tree_level_to_count_prunes_by_level > 0 )
    {
        //MRQ_copyArray(in_max_tree_level_to_count_prunes_by_level, bb.out_prune_counter_by_level, out_prune_counter_by_level );
        
        for(unsigned int k = 0; k <  in_max_tree_level_to_count_prunes_by_level; k++)
        {
            out_prune_counter_by_level[k].accumulate( bb.out_prune_counter_by_level[k]  );
        }
    }




termination:

    if(plc)		free(plc);
    if(ssrBinSumConstrs)    delete ssrBinSumConstrs;

    if( in_user_callbacks )
    {
        in_user_callbacks->bblCallBacks = NULL;
    }


    out_feasible_solution = out_best_obj < MRQ_INFINITY;
    out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
    out_clock_time = MRQ_getTime() - timeStart;

    if(in_print_level > 1)
        std::cout << "cpu time: " << out_cpu_time << "\n"; 

    if( in_user_callbacks )
    {
        if( in_user_callbacks->alg == this )
            in_user_callbacks->bblCallBacks = NULL;
    }

    algorithmFinalization(nthreads, prob, lx, ux);

    return out_return_code;
}









MRQ_BBLCallbacks::MRQ_BBLCallbacks(MRQ_MINLPProb* prob, MRQ_BranchAndBound* bb, MRQ_GeneralSolverParams* milpParams, MRQ_GeneralSolverParams* nlpParams) : BBL_UserCallbacks()
{
    initialize(prob, bb, milpParams, nlpParams);
}



MRQ_BBLCallbacks::~MRQ_BBLCallbacks()
{
    deallocate();
}



void MRQ_BBLCallbacks::deallocate()
{
    deallocateThreadStructures();
    
    MRQ_secFree(intVars); //do not desallocate reverseIntVars
    MRQ_secFree(binVars);
    MRQ_secFree(contVars);
    
    MRQ_secFree(nPCostAvgAboveErrorEstimative);
    MRQ_secFree(pCostAvgAboveErrorEstimative);
    MRQ_secDelete(oa);
    MRQ_secDelete(oaPoints);
    
    MRQ_secDelete(pcosts);
    MRQ_secDelete( pruningPcosts );
    
    MRQ_secDelete(ccstorager);
    
    binSumConstrs.deallocate();
}



void MRQ_BBLCallbacks::deallocateThreadStructures()
{
    
    MRQ_secDeleteArray(preprocessors);
    MRQ_secDeleteArray(lpBounds);
    MRQ_secDeleteArray(randomBoundUpdts);
    MRQ_secDeleteArray(chooseIndices);
    MRQ_secDeleteArray(constrChoosers);
    MRQ_secDeleteArray(gapmins);
    
    if( tplc )
    {
        for(int i = 0; i < nthreads; i++)
        {
            if( tplc[i] )
                free( tplc[i] );
        }
        
        free(tplc);
        tplc = nullptr;
    }
    
    
    if( nlps )
    {
        for(int i = 0; i < nthreads; i++)
        {
            if(nlps[i])
                delete nlps[i];
        }
        
        free(nlps);
        nlps = NULL;
    }
    
    if( oaSubs )
    {
        for(int i = 0; i < nthreads; i++)
        {
            if( oaSubs[i].points )
            {
                if( oaSubs[i].points[1] )
                    free( oaSubs[i].points[1] );
                
                free(oaSubs[i].points);
                oaSubs[i].points = nullptr;
                oaSubs[i].nPoints = 0;
            }
        }
        
        delete[] oaSubs;
        oaSubs = nullptr;
    }
    
    if(roundings)
    {
        for(int i = 0; i < nthreads; i++)
        {
            if(roundings[i])
                delete roundings[i];
        }
        
        free(roundings);
        roundings = nullptr;
    }
    
    if(ssRoundings)
    {
        for(int i = 0; i < nthreads; i++)
        {
            if(ssRoundings[i])
                delete ssRoundings[i];
        }
        
        free(ssRoundings);
        ssRoundings = nullptr;
    }
    
    if(ssrRandoms)
    {
        for(int i = 0; i < nthreads; i++)
        {
            if(ssrRandoms[i])
                delete ssrRandoms[i];
        }
        
        free(ssrRandoms);
        ssrRandoms = nullptr;
    }
    
    if( heurExecs )
    {
        for(int i = 0; i < nthreads; i++)
        {
            MRQ_HeuristicExecutor &heurExec = heurExecs[i];
            
            const int nAlgs = heurExec.getNumberOfAlgs();
            
            for(int j = 0; j < nAlgs; j++)
            {
                MRQ_Algorithm *alg = heurExec.getAlgPointer(j);
                
                alg->xInit = NULL; //we change alg->xInit to run heuristics. So, we need put NULL here to avoid a invalid free.
            }
        }
        
        delete[] heurExecs;
        heurExecs = NULL;
    }
    
    
    MRQ_secDeleteArray(userNodeGens);
    
    
    MRQ_secFree(useRoundings);
    MRQ_secFree(useHeuristics);
    MRQ_secFree(useOASubs);
    MRQ_secFree(useIGMA2s);
    
    MRQ_secFree(beforeExpoPseudoPruningCounter);
    MRQ_secFree(afterExpoPseudoPruningCounter);
    MRQ_secFree(pseudoPruningLowestLowerBound);
    
    MRQ_secFree(numberOfConstrBranchings);
    MRQ_secFree(numberOfWrongLowerBounds);
    
    if( auxInds )
    {
        for(int i = 0; i < nthreads; i++)
        {
            if(auxInds[i])
                free(auxInds[i]);
        }
        
        free(auxInds);
        auxInds = NULL;
    }
    
    
    if( auxVars )
    {
        for(int i = 0; i < nthreads; i++)
        {
            if(auxVars[i])
                free(auxVars[i]);
        }
        
        free(auxVars);
        auxVars = NULL;
    }
    
    
    if( auxConstrs )
    {
        for(int i = 0; i < nthreads; i++)
        {
            if( auxConstrs[i] )
                free( auxConstrs[i] );
        }
        
        free(auxConstrs);
        auxConstrs = NULL;
    }
    
    
}



void MRQ_BBLCallbacks::initialize(MRQ_MINLPProb* prob, MRQ_BranchAndBound* bb, MRQ_GeneralSolverParams* milpParams, MRQ_GeneralSolverParams* nlpParams)
{
    userErrorCode = 0;
    nthreads = 0;
    nI = 0;
    nbin = 0;
    nC = 0;
    
    lastOANPoints = -1;
    
    binVars = nullptr;
    nonBinVars = nullptr;
    intVars = nullptr;
    contVars = nullptr;
    reverseIntVars = nullptr;
    
    auxInds = nullptr;
    auxVars = nullptr;
    auxConstrs = nullptr;
    
    mybb = bb;
    this->prob = prob;
    milpSolverParams = milpParams;
    nlpSolverParams = nlpParams;
    
    oplc = nullptr;
    
    nPCostAvgAboveErrorEstimative = nullptr;
    pCostAvgAboveErrorEstimative = nullptr;
    
    oa = nullptr;
    oaPoints = nullptr;
    cutGen = nullptr;
    pcosts = nullptr;
    pruningPcosts = nullptr;
    ccstorager = nullptr;
    
    useRoundings = nullptr;
    useHeuristics = nullptr;
    useOASubs = nullptr;
    useIGMA2s = nullptr;
    
    beforeExpoPseudoPruningCounter = nullptr;
    afterExpoPseudoPruningCounter = nullptr;
    pseudoPruningLowestLowerBound = nullptr;
    
    numberOfConstrBranchings = nullptr;
    numberOfWrongLowerBounds = nullptr;
    
    tplc = nullptr;
    nlps = nullptr;
    lpBounds = nullptr;
    preprocessors = nullptr;
    
    thCutLists = nullptr;
    heurExecs = nullptr;
    
    
    oaSubs = nullptr;
    userNodeGens = nullptr;
    chooseIndices = nullptr;
    constrChoosers = nullptr;
    roundings = nullptr;
    ssRoundings = nullptr;
    ssrRandoms = nullptr;
    
    gapmins = nullptr;
    
    randomBoundUpdts = nullptr;
}



int MRQ_BBLCallbacks::allocateThreadStructures( const int nthreads)
{
    const bool preprocess = mybb->in_preprocess_lin_constr || mybb->in_preprocess_quad_constrs || mybb->in_preprocess_obj_function;
    const int n = prob->n, m = prob->m;
    int ret;
    
    
    this->nthreads = nthreads;
    
    
    
    if(preprocess)
    {
        preprocessors = new (std::nothrow) MRQ_Preprocessor[ nthreads];
        MRQ_calloc(tplc, nthreads); //tplc = (double **) calloc( nthreads, sizeof(double*) );
        MRQ_IFMEMERRORRETURN( !preprocessors || !tplc);
        
        for(int i = 0; i < nthreads; i++)
            preprocessors[i].initialize(prob);
    }
    
    
    MRQ_calloc(nlps, nthreads); 
    MRQ_IFMEMERRORRETURN(!nlps);
    
    
    if( mybb->in_rounding_heuristic_strategy == MRQ_RS_SSR && nbin == nI )
    {
        MRQ_calloc(ssrRandoms, nthreads);
        MRQ_calloc(ssRoundings, nthreads);
        MRQ_IFMEMERRORRETURN(!ssRoundings || !ssrRandoms);
    }
    else if( mybb->in_rounding_heuristic_strategy != MRQ_RS_NO_ROUNDING )
    {
        MRQ_calloc(roundings, nthreads);
        MRQ_IFMEMERRORRETURN(!roundings);
    }
    
    
    if( mybb->in_int_feas_heurs_strategy != MRQ_BB_IHS_NO_HEURISTICS ) //if( mybb->in_use_general_int_feas_heuristics )
    {
        heurExecs = new (std::nothrow) MRQ_HeuristicExecutor[nthreads];
        MRQ_IFMEMERRORRETURN(!heurExecs);
    }
    
    
    if( mybb->in_use_outer_app_as_heuristic )
    {
        oaSubs = new (std::nothrow) MRQ_OuterApp[nthreads];
        MRQ_IFMEMERRORRETURN(!oaSubs);
    }
    
    
    if( mybb->in_branching_strategy == MRQ_BB_BS_USER_NODE_GENERATION )
    {
        userNodeGens = new (std::nothrow) MRQ_NewUserNodeGenerator2[nthreads];
        MRQ_IFMEMERRORRETURN( !userNodeGens );
        
        //const bool use_parent_primal_sol = mybb->in_parent_sol_storing_strategy == MRQ_BB_PSSS_ONLY_PRIMAL || mybb->in_parent_sol_storing_strategy == MRQ_BB_PSSS_PRIMAL_AND_DUAL;
        //const bool use_parent_dual_sol = mybb->in_parent_sol_storing_strategy == MRQ_BB_PSSS_PRIMAL_AND_DUAL; 
        
        for(int i = 0; i < nthreads; i++)
        {
            userNodeGens[i].n = n;
            userNodeGens[i].parentBoundsStrategy = parentNodeBoundsStrategy;//userNodeGens[i].init( );
        }
    }
    else if( mybb->in_branching_strategy != MRQ_BB_BS_USER_INDEX_CHOICE )
    {
        chooseIndices = new (std::nothrow) MRQ_NewChooseIndexToBranch[nthreads];
        MRQ_IFMEMERRORRETURN(!chooseIndices);
    }
    
    if(binSumConstrs.nbinSumConstrs > 0)
    {
        constrChoosers = new (std::nothrow) MRQ_BinSumConstrsChooser[nthreads];
        MRQ_calloc( numberOfConstrBranchings, nthreads ); 
        MRQ_IFMEMERRORRETURN( !constrChoosers || !numberOfConstrBranchings );
    }
    
    
    if( mybb->in_bounds_linear_updt_strategy != MRQ_BB_BLUS_NO_UPDT )
    {
        int ml;
        
        prob->getConstraintStatistcs(&ml, NULL, NULL);
        
        
        if( ml > 0 )
        {
            lpBounds = new (std::nothrow) MRQ_LPboundsUpdater[nthreads];
            randomBoundUpdts = new (std::nothrow) MRQ_Random[nthreads];
            MRQ_IFMEMERRORRETURN(!lpBounds || !randomBoundUpdts);
            
            for(int i = 0; i < nthreads; i++)
                randomBoundUpdts[i].setSeed( &mybb->in_seed_to_random_numbers );
        }
    }
    
    
    if( mybb->in_igma2_strategy != MRQ_BB_IHS_NO_HEURISTICS  && nbin > 0 && OPT_isSolverAvailable( mybb->in_nlp_solver ) )
    {
        gapmins = new (std::nothrow) MRQ_GapMinProb[nthreads];
        MRQ_IFMEMERRORRETURN(!gapmins);
        
        
        igma2Iter.prob = prob;
        igma2Iter.nI = nI;
        igma2Iter.nC = nC;
        
        igma2Iter.intVars = intVars;
        igma2Iter.contVars = contVars;
        
        igma2Iter.in_set_max_dist_constr_on_bin_vars = mybb->in_igma2_set_max_dist_constr_on_bin_vars;
        igma2Iter.in_solve_local_search_problem_even_on_non_int_sol = mybb->in_igma2_solve_local_search_problem_even_on_non_int_sol ;
        igma2Iter.in_print_level = mybb->in_print_level;
        igma2Iter.in_neighborhood_strategy = mybb->in_igma2_neighborhood_strategy;
        igma2Iter.in_factor_to_max_dist_constr_on_bin_vars = mybb->in_igma2_factor_to_max_dist_constr_on_bin_vars;
        igma2Iter.in_percentual_to_rectangular_neighborhood = mybb->in_igma2_percentual_to_rectangular_neighborhood;
        igma2Iter.in_integer_tol = mybb->in_integer_tol;
        
    }
    
    if(  mybb->in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
    {
        MRQ_malloc(pseudoPruningLowestLowerBound, nthreads);
        MRQ_IFMEMERRORRETURN( !pseudoPruningLowestLowerBound);
        
        MRQ_setAllArray(nthreads, pseudoPruningLowestLowerBound, (double) INFINITY);
        
        if( mybb->in_count_total_prunes )
        {
            MRQ_calloc(beforeExpoPseudoPruningCounter, nthreads);
            MRQ_calloc(afterExpoPseudoPruningCounter, nthreads);
            MRQ_IFMEMERRORRETURN( !beforeExpoPseudoPruningCounter  || !afterExpoPseudoPruningCounter);
        }
        
    }
    
    if( mybb->in_calculate_pseudo_cost_average_above_error_estimative && mybb->in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
    {
        MRQ_calloc(nPCostAvgAboveErrorEstimative, nthreads);
        MRQ_calloc(pCostAvgAboveErrorEstimative, nthreads);
        MRQ_IFMEMERRORRETURN( !pCostAvgAboveErrorEstimative );
    }
    
    
    
    MRQ_calloc(useRoundings, nthreads); //useRoundings = (bool*) calloc( nthreads, sizeof(bool) );
    MRQ_calloc(useHeuristics, nthreads); //useHeuristics = (bool*) calloc( nthreads, sizeof(bool) );
    MRQ_calloc(useOASubs, nthreads); // useOASubs = (bool*) calloc( nthreads, sizeof(bool) );
    MRQ_calloc(useIGMA2s, nthreads);
    
    MRQ_calloc(auxInds, nthreads); //auxInds = (int**) calloc( nthreads, sizeof(int*) );
    MRQ_calloc(auxVars, nthreads); //auxVars = (double**) calloc( nthreads, sizeof(double*) );
    MRQ_calloc(auxConstrs, nthreads); //auxConstrs = (double**) calloc( nthreads, sizeof(double*) );
    
    MRQ_calloc(numberOfWrongLowerBounds, nthreads);
    
    MRQ_IFMEMERRORRETURN( !useRoundings || !useHeuristics || !useOASubs || !useIGMA2s || !auxInds || !auxVars || !auxConstrs || !numberOfWrongLowerBounds);
    
    
    if( mybb->in_user_callbacks )
    {//in this case, we can have cuts from user...
        thCutLists = new (std::nothrow) MRQ_NewGlobalCutList[nthreads];
        
        MRQ_IFMEMERRORRETURN(!thCutLists);
        
        for(int i = 0; i < nthreads; i++)
            thCutLists[i].initialize(nthreads);
    }
    
    
    
    //we try optimize cache usage setting NLP problem here. So we allocate memory for each thread together
    
    for(int i = 0; i < nthreads; i++)
    {
        if( preprocess )
        {
            MRQ_malloc(tplc[i], 2*m); //tplc[i] = (double *) malloc( 2*m* sizeof(double) );
            MRQ_IFMEMERRORRETURN( !tplc[i] );
        }
        
        MRQ_malloc(auxInds[i], MRQ_max(n, m)); //auxInds[i] = (int *) malloc( MRQ_max(n, m) * sizeof(int) );
        MRQ_malloc(auxVars[i], 2*n); //auxVars[i] = (double *) malloc( 2*n * sizeof(double) );
        MRQ_malloc(auxConstrs[i], m); //auxConstrs[i] = (double *) malloc( m * sizeof(double) );
        MRQ_IFMEMERRORRETURN( !auxInds[i] || !auxVars[i] || !auxConstrs[i] );
        
        
        if( preprocess )
        {
            ret = preprocessors[i].allocateMemory(n, m);
            MRQ_IFERRORRETURN( ret, MRQ_MEMORY_ERROR );
        }
        
        
        nlps[i] = OPT_newNLPSolver( mybb->in_nlp_solver );
        MRQ_IFMEMERRORRETURN(!nlps[i]);
        
        
        ret = MRQ_setNLPRelaxProb( *prob, NULL, NULL, NULL, NULL, nlps[i], true, true, true, false, mybb->thnumber + i, mybb->in_set_special_nlp_solver_params, nlpSolverParams, 1, mybb->in_max_cpu_time, mybb->in_max_time,  0, 0 );
        MRQ_IFERRORRETURN( ret, MRQ_NLP_SOLVER_ERROR );
        
        
        if( mybb->run_by_inside )
        {
            if( !std::isinf( mybb->insideSolverMaxTime) )
                nlps[i]->setMaxTime( mybb->insideSolverMaxTime );
        }
        
        if( mybb->in_use_early_branching )
        {
            nlps[i]->setRelativeOptimalityTol( mybb->in_nlp_relative_opt_tol_on_early_branching );
            
            nlps[i]->setRelativePrimalTol( mybb->in_nlp_relative_primal_tol_on_early_branching);
            
            nlps[i]->setRelativeDualTol( mybb->in_nlp_relative_dual_tol_on_early_branching);
        }
        
        if( oaSubs )
        {
            oaSubs[i].in_print_level = mybb->in_print_level - 2;
            oaSubs[i].in_use_first_nlp_relaxation = false;
            oaSubs[i].in_milp_solver = mybb->in_milp_solver;
            oaSubs[i].in_nlp_solver = mybb->in_nlp_solver;
            oaSubs[i].in_max_time = mybb->in_outer_app_subprob_time;
            oaSubs[i].in_number_of_threads = 1;
            
            //we already preprocess in each BB node...
            oaSubs[i].in_preprocess_lin_constr = false;
            oaSubs[i].in_preprocess_obj_function = false;
            oaSubs[i].in_preprocess_quad_constrs = false;
            
            oaSubs[i].in_user_callbacks = mybb->in_user_callbacks;
            oaSubs[i].in_call_update_best_sol_callback = mybb->in_call_update_best_sol_callback;
            oaSubs[i].in_call_end_of_iteration_callback = mybb->in_call_end_of_iteration_callback;
            
            
            //we preallocate space 2 linearization points in oaSub. The first is the current out_best_sol and the other is nlp->sol (solution of relaxation in the node).
            
            MRQ_malloc(oaSubs[i].points, 2); //oaSubs[i].points = (double **) malloc( 2 * sizeof(double *) );
            MRQ_IFMEMERRORRETURN( !oaSubs[i].points);
            
            
            //we allocate space for 2 array pointers, but we take advantage the solution stored at nlp->sol as the first array.
            oaSubs[i].points[0] = nlps[i]->sol;
            
            MRQ_malloc(oaSubs[i].points[1], n); //oaSubs[i].points[1] = (double *) malloc( n * sizeof(double) );
            MRQ_IFMEMERRORRETURN( !oaSubs[i].points[1] );
            
        }
        
        
        if( ssrRandoms )
        {
            ssrRandoms[i] = new (std::nothrow) MRQ_Random( mybb->in_seed_to_random_numbers );
            MRQ_IFMEMERRORRETURN( !ssrRandoms[i] );
        }
        
        if( ssRoundings )
        {
            MRQ_SSRoundingExecutor *ssr;
            
            ssRoundings[i] = new (std::nothrow) MRQ_SSRoundingExecutor;
            
            MRQ_IFMEMERRORRETURN( !ssRoundings[i] );
            
            ssr = ssRoundings[i];
            
            //ssr->in_preprocess_after_handling_constraints = true;
            ssr->in_random_order_to_threat_classes = true;
            ssr->in_random_order_to_threat_constraints_in_each_class = true;
            ssr->in_solve_minlp_as_local_search = false;
            ssr->in_max_number_of_main_iterations = 1;
            
            ssr->in_max_number_of_improvments = 0;
            ssr->in_print_level = mybb->in_print_level-2;
            ssr->in_nlp_solver = mybb->in_nlp_solver;
            
            ssr->in_absolute_feasibility_tol = mybb->in_absolute_feasibility_tol;
            ssr->in_relative_feasibility_tol = mybb->in_relative_feasibility_tol;
            
            ssr->in_integer_tol = mybb->in_integer_tol;
            ssr->in_integer_neighborhood_factor = mybb->in_integer_neighborhood_factor_to_ssr;
            
            ssr->in_cont_relax_strategy_to_stoch_rounding = MRQ_SSR_CRSSR_NO_CONTINUOUS_RELAXATION;
            
        }
        
        if( roundings )
        {
            const int strategy = mybb->in_rounding_heuristic_strategy == MRQ_RS_SSR ? MRQ_RS_PROBABILISTIC : mybb->in_rounding_heuristic_strategy ;
            
            roundings[i] = MRQ_newRounding(strategy, mybb->in_seed_to_random_numbers );
            
            MRQ_IFMEMERRORRETURN( !roundings[i] );
        }
        
        if( heurExecs )
        {
            const int MAXITERSONHEURISTICS = 20;
            
            heurExecs[i].in_max_total_clock_time = mybb->in_feas_heuristic_max_time;
            
            heurExecs[i].setMaxTimes( mybb->in_feas_heuristic_max_time );
            
            
            heurExecs[i].setPrintLevels( mybb->in_print_level - 2 );
            
            heurExecs[i].setMaxIters( MAXITERSONHEURISTICS);
            
            heurExecs[i].setNumberOfThreads(1);
            
            heurExecs[i].setSolvers( mybb->in_milp_solver, mybb->in_nlp_solver );
            
            heurExecs[i].in_use_diving = mybb->in_use_feas_heuristic_diving;
            heurExecs[i].in_use_fp = mybb->in_use_feas_heuristic_fp;
            heurExecs[i].in_use_oafp = mybb->in_use_feas_heuristic_oafp;
            heurExecs[i].in_use_rens = mybb->in_use_feas_heuristic_rens;
            
            heurExecs[i].in_use_igma1 = false;
            heurExecs[i].in_use_igma2 = false;
            heurExecs[i].in_use_ssr = false;
            
            const int nAlgs = heurExecs[i].getNumberOfAlgs();
            
            for( int j = 0; j < nAlgs; j++ )
            {
                MRQ_Algorithm *alg = heurExecs[i].getAlgPointer(j);
                
                alg->in_preprocess_lin_constr = false;
                alg->in_preprocess_quad_constrs = false;
                alg->in_preprocess_obj_function = false;
            }
            
        }
        
        
        /*if( userNodeGens )
        {
            const int r = userNodeGens[i].allocate(n);
            if( r != 0 )
            {
                if(mybb->in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(r);
                return r;
            }
        } */
        
        
        
        if( mybb->in_branching_strategy != MRQ_BB_BS_USER_NODE_GENERATION && mybb->in_branching_strategy != MRQ_BB_BS_USER_INDEX_CHOICE )
        {
            const int maxBranchVars = mybb->in_number_of_branching_vars;
            
            const int r = chooseIndices[i].allocate( maxBranchVars, n);
            MRQ_IFERRORRETURN(r, r);
        }
        
        if(binSumConstrs.nbinSumConstrs > 0)
        {
            const int r = constrChoosers[i].reallocate( binSumConstrs.nbinSumConstrs);
            MRQ_IFERRORRETURN(r, r);
        }
        
        
        if( lpBounds )
        {
            double zu = getUpperBound();
            
            const int r = lpBounds[i].buildProblem( mybb->in_nlp_solver, *prob, zu );
            MRQ_IFERRORRETURN(r, r);
            
            //lpBounds[i].lp->generateModelFile("lpmodel.lp");
            //MRQ_getchar();
        }
        
        
        if(gapmins)
        {
            const double *lx = prob->lx;
            const double *ux = prob->ux;
            
            int r = gapmins[i].setProblem(mybb->in_nlp_solver, *prob, lx, ux, nlpSolverParams, i, mybb->in_igma2_set_special_gap_min_solver_params, false, false, 1, mybb->in_max_cpu_time, mybb->in_max_time, 0,  1 + (int) mybb->in_igma2_set_max_dist_constr_on_bin_vars   );
            MRQ_IFERRORRETURN(r, r);
            
            ((OPT_NLPSolver*) gapmins[i].solver)->in_absolute_feas_tol = mybb->in_absolute_feasibility_tol;
            ((OPT_NLPSolver*) gapmins[i].solver)->in_relative_feas_tol = mybb->in_relative_feasibility_tol;
        }
        
    }
    
    
    
    
    
    
    return 0;
}







int MRQ_BBLCallbacks::beforeAll(const unsigned int numberOfThreads, double *lx, double *ux)
{
    const int n = prob->n;
    const int m = prob->m;
    const int NUMBEROFINITIALOAPOINTSALLOC = 20;
    int ret;


    #if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Entering at MRQ_BBLCallbacks::beforeAll" << std::endl;
    }
    #endif




    userErrorCode = 0;

    lastOANPoints = -1;
    constrBranchStrat = mybb->in_constr_branching_strategy;

    nI = prob->getNumberOfIntegerVars();
    nC = n - nI;

    MRQ_malloc(intVars, nI+n);
    MRQ_malloc(contVars, nC);
    MRQ_IFMEMERRORRETURN( !intVars || !contVars );

    prob->getContinuousIndices(contVars);
    prob->getIntegerIndices( intVars );

    reverseIntVars = &intVars[nI];

    prob->getReverseIntegerIndices( reverseIntVars );
    
    
    nbin = 0;
    for(int i = 0; i < nI; i++)
    {
        const int ind = intVars[i];
        if( lx[ind] > -1.0 && lx[ind] <= 1.0 && ux[ind] >= 0.0 &&  ux[ind] < 2.0  )
            nbin++;
    }
    
    

    if( mybb->in_branching_strategy == MRQ_BB_BS_BIN_FIRST_HIGHEST_INT_GAP || mybb->in_igma2_strategy != MRQ_BB_IHS_NO_HEURISTICS )
    {
        MRQ_malloc(binVars, nI);
        MRQ_IFMEMERRORRETURN( !binVars );
        
        
        nonBinVars = &binVars[nbin];
        
        
        int kb = 0, knb = 0;
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            if( lx[ind] > -1.0 && lx[ind] <= 1.0 && ux[ind] >= 0.0 && ux[ind] < 2.0 )
            {
                binVars[kb] = ind;
                kb++;
            }
            else
            {
                nonBinVars[knb] = ind;
                knb++;
            }
        }
        
        #if MRQ_DEBUG_MODE
            assert( kb == nbin );
            assert( kb + knb == nI );
        #endif
        
    }
    /*else
    {
        nbin = -1;//we do not need calculate this value
    }*/
    
    
    ccstorager = new (std::nothrow) minlpproblem:: MIP_ConstraintsByColumnsStorager;
    MRQ_IFMEMERRORRETURN(!ccstorager);
    
    {
        int r = ccstorager->storageConstraintsByColumns(*prob, mybb->in_preprocess_quad_constrs);
        MRQ_IFERRORRETURN(r, MRQ_UNDEFINED_ERROR);
    }
    
    
    if( (constrBranchStrat != MRQ_BB_CBS_NO_CONSTRAINT_BRANCH && origBBLBranchStrategy != BBL_BS_USER_NODE_GENERATION) || mybb->in_rounding_heuristic_strategy != MRQ_RS_NO_ROUNDING )
    {
        //if you change plc and puc definiton, please consider change also definiton in generateNodes method
        const double *plc = oplc;
        const double *puc = oplc ? &(plc[m]) : NULL; 
        
        int r = binSumConstrs.calculateIndices(*prob, lx, ux, plc, puc);
        MRQ_IFERRORRETURN(r, r);
        
        if( mybb->in_print_level > 5 )
        {
            std::cout << MRQ_PREPRINT "Equalities constraints of binary sum:\n";
            for(int j = 0; j < binSumConstrs.nbinSumConstrs; j++)
                std::cout << "c["<<j<<"]: " << binSumConstrs.binSumConstrs[j] << " \t";
            std::cout << "\n";
        }
        //MRQ_getchar();
    }




    if( mybb->in_branching_strategy == MRQ_BB_BS_STBRANCH_PSEUDO_COSTS || ( binSumConstrs.nbinSumConstrs > 0 && constrBranchStrat == MRQ_BB_CBS_STBRANCH_PSEUDO_COSTS ) || mybb->in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
    {
        pcosts = new (std::nothrow) MRQ_NewPseudoCostCalc;
        MRQ_IFMEMERRORRETURN( !pcosts );
        
        ret = pcosts->allocate( nI );
        MRQ_IFERRORRETURN(ret, ret);
        
        if(numberOfThreads > 1)
        {
            //MRQ_getchar();
            //we need lock other threads before bblop because we will open several threads to calculate psudeo cuusts in the first iteration of B&B
            lockAuxiliaryThreadsBeforeBBLoop();
        }
    }
    

    /*if( mybb->in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING )
    {
        pruningPcosts = new (std::nothrow) MRQ_BasePseudoCostCalc;
        MRQ_IFMEMERRORRETURN( !pruningPcosts );
        
        ret = pruningPcosts->allocate(nI);
        MRQ_IFERRORRETURN(ret, ret);
    }*/


    if( mybb->in_use_outer_app )
    {
        oa = new (std::nothrow) MRQ_OuterApp(); //if we get a memory error, we let algorithm follow without outer approximation.
        MRQ_IFMEMERRORRETURN(!oa);
        
        oa->in_print_level = mybb->in_print_level - 2;
        oa->in_use_first_nlp_relaxation = false;
        oa->in_milp_solver = mybb->in_milp_solver;
        oa->in_nlp_solver = mybb->in_nlp_solver;
        oa->in_preprocess_lin_constr = false;
        oa->in_preprocess_obj_function = false;
        oa->in_preprocess_quad_constrs = false;
        
        oa->in_user_callbacks = mybb->in_user_callbacks;
        oa->in_call_update_best_sol_callback = mybb->in_call_update_best_sol_callback;
        oa->in_call_end_of_iteration_callback = mybb->in_call_end_of_iteration_callback;
        
        oa->in_max_time = mybb->in_outer_app_time;
        
        oa->in_number_of_threads = 1;
        
        oa->in_relative_convergence_tol = mybb->in_relative_convergence_tol;
        oa->in_absolute_convergence_tol = mybb->in_absolute_convergence_tol;
        
        oa->in_store_history_solutions = true;
        
        oaPoints = new (std::nothrow) MRQ_NewPoints( NUMBEROFINITIALOAPOINTSALLOC );
        MRQ_IFMEMERRORRETURN( !oaPoints );
        
        iterNextOAApplic = 1;
    }
    else
    {
        iterNextOAApplic = ULONG_MAX;
    }



    ret = allocateThreadStructures(numberOfThreads);
    MRQ_IFERRORRETURN(ret, ret);


    if( mybb->in_user_callbacks )
    {//in this case, we can have cuts from user...
        
        cutGen = new (std::nothrow) MRQ_NewGlobalCutGenerator( numberOfThreads, thCutLists );
        MRQ_IFMEMERRORRETURN(!cutGen);
        
        
        for(unsigned int i = 0; i < numberOfThreads; i++)
            thCutLists[i].cutGen = cutGen;
        
        mybb->in_user_callbacks->cutGenerator = cutGen;
        
        
        const int r = mybb->in_user_callbacks->beforeAll( MRQ_BB_ALG, numberOfThreads );
        if( r != 0 )
        {
            if( mybb->in_print_level > 0 )
                MRQ_PRINTCALLBACKERRORNUMBER(r);
            
            userErrorCode = r;
            return MRQ_STOP_REQUIRED_BY_USER;
        }
    }



    #if 0
    if( mybb->in_int_heurs_strategy != muriqui::MRQ_BB_IHS_NO_HEURISTICS ) //if( mybb->in_use_general_int_feas_heuristics )
    {
        double zu = getUpperBound();
        double obj;
        
        double *intSol = auxVars[0];
        MRQ_HeuristicExecutor &heurExec = heurExecs[0];
        MRQ_Algorithm *alg;
        
        
        const int r = heurExec.insideRun(*prob, milpSolverParams, nlpSolverParams, zu, obj, intSol, true, &alg, 0, 1, mybb->in_max_int_heuristic_time, lx, ux);
        
        if( r == MRQ_HEURISTIC_SUCCESS )
        {
            tryUpdateBestSolution(0, intSol, obj, 0, true);
            
            if( r == MRQ_OPTIMAL_SOLUTION )
            {
                tryUpdateLowerBound(alg->in_lower_bound);
            }
        }
        
        for(int i = 0; i < nthreads; i++)
        {
            for(int j = 0; j < heurExecs->nAlgs; j++ )
            {
                heurExecs[i].__getAlgPointer(j)->in_use_initial_solution = true;
            }
        }
    } 

    #endif


    #if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Leaving MRQ_BBLCallbacks::beforeAll" << std::endl;
        
        MRQ_closeFileToThread(tid);
    }
    #endif


    return 0;
}




int MRQ_BBLCallbacks::generateRootNode( BBL_Node*& rootNode)
{
    int ret;
    MRQ_NewBBNode *myRoot = NULL;
    
    
    if( mybb->in_branching_strategy == MRQ_BB_BS_USER_NODE_GENERATION || mybb->in_call_generate_root_node_callback )
    {
        ret = mybb->in_user_callbacks->BB_generateRootNode(myRoot);
        
        rootNode = myRoot;
        
        if( ret != 0 )
        {
            if( mybb->in_print_level > 0)
                MRQ_PRINTCALLBACKERRORNUMBER(ret);
            
            userErrorCode = ret;
            return MRQ_STOP_REQUIRED_BY_USER;
        }
        
        if(myRoot == NULL)
            MRQ_PRINTMSG("Warning: user callback BB_generateRootNode did not generate root node to branch-and-bound. Generating a standard root node\n");
    }
    
    if( myRoot == NULL )
    {
        myRoot = new (std::nothrow) MRQ_NewBBNode(parentNodeBoundsStrategy, prob->n);
        rootNode = myRoot;
        MRQ_IFMEMERRORRETURN( !myRoot );
    }
    
    
    return 0;
}



int MRQ_BBLCallbacks::beforeBBLoop(const unsigned int thnumber, const double lb, const double ub)
{
    if( mybb->in_call_before_bb_loop_callback )
    {
        int r = mybb->in_user_callbacks->beforeBBLoop(thnumber, lb, ub, *nlps[thnumber] );
        
        if( r != 0 )
        {
            if( mybb->in_print_level > 0 )
                MRQ_PRINTCALLBACKERRORNUMBER(r);
            
            userErrorCode = r;
            return MRQ_STOP_REQUIRED_BY_USER;
        }
    }
    
    return 0;
}



void MRQ_BBLCallbacks::afterBBLoop(const unsigned int thnumber, const double lb, const double ub, const int threadReturnCode)
{
    if( mybb->in_call_after_bb_loop_callback )
    {
        mybb->in_user_callbacks->afterBBLoop(thnumber, lb, ub, *nlps[thnumber], threadReturnCode);
    }
}



//return false if ind is not in the array a
inline bool MRQ_tryUpdateBndsInNodeBndsSol(const unsigned int size, branchAndBound::BBL_ClassUnionNodeBoundsSolPointer &a, const unsigned int ind, const double lind, const double uind)
{
    bool found = false;
    
    unsigned int bind;
    double bl, bu;
    
    for(unsigned int i = 0; i < size; i++)
    {
        a.getArrayElement(i, &bind, &bl, &bu, NULL);
        
        if( bind == ind ) //if( a[i].ind == ind )
        {
            //a[i].l = lind;
            //a[i].u = uind;
            
            a.setArrayElement(i, NULL, &lind, &uind, NULL); //do not update sol. we need the original value...
            
            found = true;
            break;
        }
    }
    
    return found;
}


//myBoundsReallocated should be intializaed by user as false... It will be true if we need reallocate node.myBounds
int MRQ_tryUpdtBondsByLP(const int index, const double tol, const int maxAuxBndsSize, bool &myBoundsReallocated, MRQ_LPboundsUpdater &lpBound, double *nlx, double *nux, MRQ_NewBBNode &node )
{
    bool updt = false;
    int r;
    double value;
    
    
    r = lpBound.calcMinValueOnVariable( index, value);
    
    if( r == MRQ_OPTIMAL_SOLUTION && value - tol > nlx[index] )
    {
        //std::cout << "lpbounds: atualizamos limite inferior na coordenada " << index << " nlx: " << nlx[index] << " value: " << value << " \n";
        //MRQ_getchar();
            
        nlx[index] = ceil( value ); //since we already updated the bound, we prefer do not solve the maximization problem on this same variable...
        updt = true;
    }
    else if( r == MRQ_INFEASIBLE_PROBLEM )
    {
        std::cout << "lpbounds: detectei inviablidade \n";
        MRQ_getchar();
        
        //we detect infeasibility in this node
        return MRQ_INFEASIBLE_PROBLEM;
    }
    else
    {
        //since we cannot update lower bound, we try update upper bound
        
        r = lpBound.calcMaxValueOnVariable(index, value);
        
        if( r == MRQ_OPTIMAL_SOLUTION && value + tol < nux[index] )
        {
            std::cout << "lpbounds: atualizamos limite superior na coordenada " << index << " nux: " << nux[index] << " value: " << value << " \n";
            MRQ_getchar();
        
            nux[index] = floor(value);
            updt = true;
        }
        else if( r == MRQ_INFEASIBLE_PROBLEM )
        {
            //std::cout << "lpbounds: detectei inviablidade \n";
            //MRQ_getchar();
        
            //we detect infeasibility in this node
            return MRQ_INFEASIBLE_PROBLEM;
        }
        
    }
    
    
    
    if( updt )
    {
        //checkinf if parent bounds already consider this variable...
        
        const unsigned int nMybounds = node.nMyBounds;
        
        const bool found = MRQ_tryUpdateBndsInNodeBndsSol( nMybounds, node.myBounds, index, nlx[index], nux[index] );
        
        
        if( !found )
        {
            
            if( !myBoundsReallocated )
            {
                //MRQ_NodeBoundsSol *auxBnds;
                
                //auxBnds = (MRQ_NodeBoundsSol *) realloc( node.myBounds, maxAuxBndsSize* sizeof(MRQ_NodeBoundsSol) );
                
                //we cannot use MRQ_realloc because we cannot create a reference of packed structured member...
                //int r = MRQ_realloc(node.myBounds, maxAuxBndsSize);
                //MRQ_IFERRORRETURN(r, r);
                
                {
                    /*void *p = realloc(node.myBounds, maxAuxBndsSize * sizeof(*node.myBounds) );
                    MRQ_IFMEMERRORRETURN(!p);
                    node.myBounds = (decltype(node.myBounds)) p;*/
                    
                    unsigned int oNMyBounds = node.nMyBounds;
                    
                    r = node.allocateNodeBounds(maxAuxBndsSize); //allocateNodeBounds performs a realloc
                    MRQ_IFERRORRETURN(r, r);
                    
                    node.nMyBounds = oNMyBounds;
                }
                
                //node.myBounds = auxBnds;
                myBoundsReallocated = true;
            }
            
            /*MRQ_NodeBoundsSol el;
                
            el.ind = index;
            el.l = nlx[index];
            el.u = nux[index];
            el.sol = NAN; //we use NAN like a flag to do not update pseudo cost over this index...
            
            branchAndBound::BBL_addElementOnSortedNodeBoundsArray( nMybounds, node.myBounds, el );*/
            
            node.myBounds.addElementOnSortedArray(node.nMyBounds, index, nlx[index], nux[index], NAN) ;
            
            node.nMyBounds++;
        }
        
        //node.print();
        //MRQ_getchar();
    }
    
    
    
    
    
    return 0;
}



int MRQ_BBLCallbacks::linearBoundsUpdating( MRQ_BB_BOUND_LIN_UPDT_STRATEGY strategy, const double zu, MRQ_NewBBNode &node, double *nlx, double *nux, bool *auxVars, MRQ_LPboundsUpdater &lpBound, MRQ_Random &random)
{
    const int n = prob->n;
    const double tol = 1.0e-5;
    
    
    {
        //checking if there is some unfixed integer variable...
        bool flag = true;
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            if( nlx[ind] != nux[ind] )
            {
                flag = false;
                break;
            }
        }
        
        
        if( flag )
            return 0;
    }
    
    
    lpBound.setVariablesBounds(n, nlx, nux);
    lpBound.updateObjectiveCut(*prob, zu);
    
    
    if( strategy == MRQ_BB_BLUS_ONE_VAR )
    {
        bool reallocBounds = false;
        unsigned int irand;
        int r;
        
        
        do
        {
            irand = intVars[ random.randInt(0, nI-1) ];
            
        }while( nlx[irand] == nux[irand] ); //if all integer variables are fixed, we have a problem here...
        
        
        r = MRQ_tryUpdtBondsByLP( irand, tol, node.nMyBounds + 1, reallocBounds, lpBound, nlx, nux, node );
        
        return r;
    }
    else if( strategy == MRQ_BB_BLUS_ALL_VARS )
    {
        bool bndsRealloc = false;
        int r;
        
        
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            if( nlx[ind] == nux[ind] )
                continue;
            
            
            r = MRQ_tryUpdtBondsByLP( ind, tol, nI, bndsRealloc, lpBound, nlx, nux, node );
            
            if( r != 0 )
            {
                #if MRQ_DEBUG_MODE
                    if( r != MRQ_INFEASIBLE_PROBLEM )
                        MRQ_PRINTERRORNUMBER(r);
                #endif
                return r;
            }
        }
        
    }
    else if( strategy == MRQ_BB_BLUS_SUBGROUP )
    {
        const double factor = mybb->in_bounds_updt_factor_to_subgroup;
        bool bndsRealloc = false;
        int nunfix = 0;
        int irand;
        int r;
        
        
        MRQ_setAllArray(nI, auxVars, false);
        
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            if( nlx[ind] != nux[ind] )
            {
                auxVars[i] = true;
                nunfix++;
            }
        }
        
        
        const int ninds = MRQ_min<int>( ceil(factor*nI) , nunfix ); //number of index to be tried
        
        
        for(int i = 0; i < ninds; i++)
        {
            int ii;
            
            do
            {
                ii = random.randInt(0, nI-1);
            }while( !auxVars[ii] );
            
            
            irand = intVars[ii];
            
            r = MRQ_tryUpdtBondsByLP( irand, tol, nI, bndsRealloc, lpBound, nlx, nux, node );
            
            if( r != 0 )
            {
                #if MRQ_DEBUG_MODE
                    if( r != MRQ_INFEASIBLE_PROBLEM )
                        MRQ_PRINTERRORNUMBER(r);
                #endif
                return r;
            }
            
            auxVars[ii] = false; //this avriable canoot be choosen again...
        }
        
    }
    else
    {
        MRQ_PRINTERRORMSG("Bounds updatingt strategy not implemented");
    }
    
    return 0;
}








//int MRQ_BBLCallbacks::beforeSolvingRelaxation( const unsigned int thnumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, bool &pruneNode)


int MRQ_BBLCallbacks::solveSubProblem(const unsigned int thnumber, BBL_Node &node, const long unsigned int iter, const double lb, const double myub, double *nlx, double *nux, BBL_RETURN_CODES &optRetCode, double &objValue, double &dualObjValue, double *sol, double *dualSol, bool &generalFeasibleSol, bool &pruneNode, double &nodeLowerBound, BBL_BRANCH_STRATEGY &branchStrategy)
{
    #if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Entering at MRQ_BBLCallbacks::solveSubProblem" << std::endl;
    }
    #endif


    const bool useDualSol = mybb->in_parent_sol_storing_strategy == MRQ_BB_PSSS_PRIMAL_AND_DUAL;
    const int n = prob->n;
    const int m = prob->m;
    double ub = myub;
    double zu_tol = ub < MRQ_INFINITY ? MRQ_zuWithTol(ub, mybb->in_absolute_convergence_tol, mybb->in_relative_convergence_tol) : MRQ_INFINITY;

    bool branchEvenInteger = false, updtBounds;
    bool relaxIntSol;
    int r;
    int functionReturnCode;
    //double lbnode;

    MRQ_NewBBNode &mynode = (MRQ_NewBBNode&) node;


    //thread structures
    MRQ_Preprocessor *prepocess = preprocessors ? &preprocessors[thnumber] : nullptr;
    MRQ_NLPSolver *nlp = nlps[thnumber];
    int *auxInd = auxInds[thnumber];
    double *auxConstr = auxConstrs[thnumber];
    double *auxVar = this->auxVars[thnumber];
    double *auxVar2 = &auxVar[n];

    double *myplc;
    double *mypuc;


    const bool useRounding = mybb->in_rounding_heuristic_strategy != MRQ_RS_NO_ROUNDING && ( ub >= MRQ_INFINITY || (iter-1) % mybb->in_rounding_heuristic_call_iter_frequence == 0 ) ;


    const bool useHeuristic = (mybb->in_int_feas_heurs_strategy == MRQ_BB_IHS_ALWAYS || (mybb->in_int_feas_heurs_strategy == MRQ_BB_IHS_UNTIL_FIRST_FEAS_SOL && ub >= BBL_INFINITY) ) && ( (iter-1) % mybb->in_feas_heuristic_frequency == 0) ;

    const bool useOASub = mybb->in_use_outer_app_as_heuristic && (iter % mybb->in_outer_app_subprob_frequence == 0); //here, we test if iter % is equal to zero because we do not want run oaSub in the first iteration. SO, DO NOT PUT iter-1

    const bool useIGMA2 = (mybb->in_igma2_strategy == MRQ_BB_IHS_ALWAYS || (mybb->in_igma2_strategy == MRQ_BB_IHS_UNTIL_FIRST_FEAS_SOL && ub >= BBL_INFINITY ) ) && (gapmins) && ( (iter-1) % mybb->in_igma2_frequency == 0);

    const bool useBoundUpdt = mybb->in_bounds_linear_updt_strategy != MRQ_BB_BLUS_NO_UPDT && ( (iter-1) % mybb->in_bounds_updt_frequency == 0 );


    generalFeasibleSol = false;
    pruneNode = false;
    
    
    #if MRQ_DEBUG_MODE
        assert( (int) mynode.getDepth() <= nbin || nbin != nI ); // if the problem is binary, we cannot have more levels in the BB tree than the number of variables (but it is not necessraily true if the problem is no binary)
    #endif
    
    

    if( useHeuristic )
        useHeuristics[thnumber] = true;

    if( useRounding )
        useRoundings[thnumber] = true;

    if( useOASub )
        useOASubs[thnumber] = true;

    if( useIGMA2 )
        useIGMA2s[thnumber] = true;

    #if MRQ_SET_LBHEUR_ON_BB_NODE
        if( mynode.heurlb >= zu_tol )
        {
            //simulating a bound prunning
            
            objValue = dualObjValue = INFINITY;
            optRetCode = BBL_OPTIMAL_SOLUTION;
            
            //pruneNode = true;
            functionReturnCode = 0;
            goto termination;
        }
    #endif

    
    if( mybb->in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING && ub < MRQ_INFINITY)
    {
        const double nodeLowerBound = node.getLowerBound();
        const double nodeGap = ub - nodeLowerBound;
        
        if( nodeGap < mybb->in_absolute_convergence_tol_for_pseudo_pruning  ||  nodeGap/MRQ_abs(ub)  < mybb->in_relative_convergence_tol_for_pseudo_pruning  )
        {
            double objIncreasiestimativeToPseudoPruning = 0.0;
            const double ub_slack =  ub + mybb->in_relative_upper_bound_slack_factor_for_pseudo_pruning * MRQ_abs(ub) + mybb->in_absolute_upper_bound_slack_for_pseudo_pruning;
            MRQ_BasePseudoCostCalc *myPcosts = this->pcosts;
            bool allVarIncreasingEstimated;
            
            myPcosts->getObjIncreaseEstimative(mybb->in_only_apply_pseudo_pruning_on_fixed_integer_vars, reverseIntVars, mybb->in_min_number_of_bound_prunes_per_var_before_pseudo_pruning, node.nMyBounds, node.myBounds, mybb->in_alpha_to_balance_estimative_in_pseudo_pruning, objIncreasiestimativeToPseudoPruning, allVarIncreasingEstimated);
            
            if( nodeLowerBound + objIncreasiestimativeToPseudoPruning >= ub_slack ) //maybe we could use zu_tol
            {
                /*printf("\nPodando por pseudo poda! nodeLowerBound: %f total estimative: %f ub: %f ub_slack: %f\n", nodeLowerBound, nodeLowerBound + estimative, ub, ub_slack);
                node.print();
                printf("pruning pcosts: \n");
                myPcosts->print(nI, intVars, std::cout);
                MRQ_getchar();*/
                
                /*{
                    unsigned int k;
                    double l, u, s;
                    node.myBounds.getArrayElement(0, &k, &l, &u, &s);
                    
                    if( u == 0.0 )
                        MRQ_getchar();
                }*/
                
                
                if( mybb->in_count_total_prunes )
                    beforeExpoPseudoPruningCounter[thnumber]++;
                
                if( nodeLowerBound < pseudoPruningLowestLowerBound[thnumber] )
                        pseudoPruningLowestLowerBound[thnumber] = nodeLowerBound;
                
                //TODO: count pseudo prunning by level also!
                
                //pseudo pruning
                pruneNode = true;
                functionReturnCode = 0;
                goto termination;
            }
            
        }
        
    }





    /*The ideal here is use a semaphore to update mybb->zl. However, we do not use mybb->zl to anything in this algorithm. We only put it here yo in_user_callbacks have a reference to get lower bound. So, it would be too expensive set a sempahore in each B&B iteration only for that purpose. The global lower bound used to stop the algorithm is inside BBL, in this lower bound is updated in a accurated way. */
    if( lb > mybb->zl )
        mybb->zl = lb;



    if( zu_tol < MRQ_INFINITY && pcosts ) // && mybb->in_branching_strategy == MRQ_BB_BS_STBRANCH_PSEUDO_COSTS )
    {
        bool prune;
        
        prune = pcosts->updateVarBoundsByUpperBound(nI, intVars, zu_tol, nlx, nux);
        
        if(prune)
        {
            #if MRQ_DEBUG_MODE
                //std::cout << "Pruning node by bound second bound from first iteration pseudocust." << MRQ_GETFILELINE << "\n"  ;
                //MRQ_getchar();
            #endif
            
            //simulating a prune by bound
            
            objValue = dualObjValue = INFINITY;
            optRetCode = BBL_OPTIMAL_SOLUTION;
            
            functionReturnCode = 0;
            goto termination;
        }
        //MRQ_getchar();
        //we should increase user bound prune counter, but BBL will count like a user prune... :(
    }



    if(prepocess)
    {
        bool updtvb, updtcb;
        //const double *opuc = oplc ? &oplc[m] : NULL;
        const double *opuc = &oplc[m];
        
        myplc = tplc[thnumber];
        mypuc = &myplc[m];
        
        
        /*double nlx2[n], nux2[n], myplc2[m], mypuc2[m];
                
        MRQ_copyArray(n, nlx, nlx2);
        MRQ_copyArray(n, nux, nux2);
        MRQ_copyArray(m, myplc, myplc2);
        MRQ_copyArray(m, mypuc, mypuc2); */
        
        
        
        
        
        //r = prepocess->preprocess( mybb->in_preprocess_quad_constrs, mybb->in_preprocess_obj_function, zu_tol, nlx, nux, updtvb, updtcb, oplc, opuc, myplc, mypuc );
        
        r = prepocess->preprocess(nI, intVars, *ccstorager, mybb->in_preprocess_quad_constrs, mybb->in_preprocess_obj_function, zu_tol, nlx, nux, updtvb, updtcb, oplc, opuc, myplc, mypuc );
        
        #if 0
        {
            r = prepocess->preprocess(0, NULL, *ccstorager, mybb->in_preprocess_quad_constrs, mybb->in_preprocess_obj_function, zu_tol, nlx2, nux2, updtvb, updtcb, oplc, opuc, myplc2, mypuc2 );
            
            {
                for(int i = 0; i < n; i++)
                {
                    if( MRQ_abs(nlx[i] - nlx2[i] ) > 1e-6 )
                        printf("nlx[%d]: %0.16lf    nlx2[%d]: %0.16lf\n", i, nlx[i], i, nlx2[i]);
                    
                    if( MRQ_abs( nux[i] - nux2[i] ) > 1e-6 )
                        printf("nux[%d]: %0.16lf    nux2[%d]: %0.16lf\n", i, nux[i], i, nux2[i]);
                    
                    assert( MRQ_abs(nlx[i] - nlx2[i]) <= 1e-6 );
                    assert( MRQ_abs(nux[i] - nux2[i]) <= 1e-6 );
                }
            }
        }
        #endif
        
        if( r == minlpproblem::MIP_INFEASIBILITY )
        {
            //MRQ_PRINTERRORMSG("Podei no por preprocessamento!");
            optRetCode = BBL_INFEASIBLE_PROBLEM;
            //pruneNode = true;
            functionReturnCode = 0;
            goto termination;
        }
        
        for( int i = 0; i < m; i++ )
            r += nlp->setConstraintBounds(i, myplc[i], mypuc[i]);
        
        #if MRQ_DEBUG_MODE
            if( r!= 0 )
            {
                if( mybb->in_print_level > 0 )
                    MRQ_PRINTERRORNUMBER(r);
                functionReturnCode = MRQ_NLP_SOLVER_ERROR;
                goto termination;
            }
        #endif
    }
    else
    {
        myplc = oplc;
        mypuc = &myplc[m];
    }
    

    if( useBoundUpdt )
    {
        MRQ_LPboundsUpdater &lpBound = lpBounds[thnumber];
        MRQ_Random &random = randomBoundUpdts[thnumber];
        
        bool *auxFlags = (bool *) auxInd;
        
        ub = getBestSolutionObj();
        zu_tol = ub < MRQ_INFINITY ? MRQ_zuWithTol(ub, mybb->in_absolute_convergence_tol, mybb->in_relative_convergence_tol) : MRQ_INFINITY;
        
        
        int r = linearBoundsUpdating( mybb->in_bounds_linear_updt_strategy, zu_tol, mynode, nlx, nux, auxFlags, lpBound, random );
        
        if( r == MRQ_INFEASIBLE_PROBLEM )
        {
            optRetCode = BBL_INFEASIBLE_PROBLEM;
            functionReturnCode = 0;
            goto termination;
        }
        
    }


    
    /*we believe the best moment to aply heuristic is after solving the currente relaxation. However, it can be advantageous play before solving current relaxation in the first iteration since having a feasible soltuion can bring better results aplying strong branching pseudo costs calculation because it allow fixing some variable by bound*/
    if( iter == 1u && pcosts && useHeuristics[thnumber] )
    {
        double obj;
        double *intSol = auxVar;
        MRQ_HeuristicExecutor &heurExec = heurExecs[thnumber];
        
        
        const int r = heurExec.insideRun( *prob, milpSolverParams, nlpSolverParams, getBestSolutionObj(), obj, intSol, true, NULL, mybb->thnumber + thnumber, 1, mybb->in_feas_heuristic_max_time, nlx, nux );
        
        if(r == MRQ_HEURISTIC_SUCCESS || r == MRQ_OPTIMAL_SOLUTION)
        {
            tryUpdateBestSolution(thnumber, intSol, obj, iter, true);
        }
        
        useHeuristics[thnumber] = false;
    }
            
            
    


    if( mybb->in_call_before_solving_relax_callback )
    {
        r = mybb->in_user_callbacks->BB_beforeSolvingRelaxation(thnumber, mynode, iter, lb, ub, nlx, nux, *nlp, pruneNode );
        
        if(r != 0)
        {
            if( mybb->in_print_level > 0 )
                MRQ_PRINTCALLBACKERRORNUMBER(r);
            
            userErrorCode = r;
            functionReturnCode = MRQ_STOP_REQUIRED_BY_USER;
            goto termination;
        }
        
        if( pruneNode )
        {
            functionReturnCode = 0;
            goto termination;
        }
    }


    //Now, we start the procedures to solve current relaxation


    //checking by global cuts to add...
    if( thCutLists )
    {
        if( !thCutLists[thnumber].isempty )
            thCutLists[thnumber].addCutsOnSolver(nlp);
    }



    if( node.getParentInfo()->a.xParent )
    {
        const double *initSol = node.getParentInfo()->a.xParent;
        nlp->setInitialSolution( initSol, useDualSol ? &initSol[n] : NULL, useDualSol ? &initSol[n+m] : NULL );
    }



    if( mybb->in_use_dynamic_constraint_set )
    {
        MRQ_DynConstrSetSetter dcsSetter;
        dcsSetter.setDyncConstrSet( mybb->ndcs0, mybb->dcs0, mybb->ndcs1, mybb->dcs1, nlx, nux, prob->nlConstr, nlp );
    }
    
    


    do
    {
        bool solveAgain;
        updtBounds = false;
        
        if( zu_tol < MRQ_INFINITY )
            nlp->setObjCutUpperBound(zu_tol);
        
        r = nlp->setnVariablesBounds(n, nlx, nux);
        
        
        do
        {
            solveAgain = false;
            
            nlp->solve(false);
            
            
            #if MRQ_DEBUG_MODE
            {
                const double nodelb = mynode.getLowerBound();
                if( nlp->retCode == optsolvers::OPT_OPTIMAL_SOLUTION  &&  nlp->objValue < nodelb - 0.05 - MRQ_abs(0.001*nodelb) )
                {
                    numberOfWrongLowerBounds[thnumber]++;
                    
                    //std::cout << MRQ_PREPRINT "warning: node relaxation lower than lower bound. iter: " << iter << " nlp.objValue: " << nlp->objValue << " lb: " << lb << " nlp.retCode: " << nlp->retCode << " nlp->origSolverRetCode: " << nlp->origSolverRetCode << " node lb: " << nodelb << "\n";
                    #if MRQ_DEBUG_IGMA2_BUG
                        std::cout << "igma2OnAncestral: " << mynode.igma2OnAncestral << "\n" ;
                    #endif
                    /*for(int i = 0; i < n; i++)
                        printf("i: %d nlx: %f nux: %f\n", i, nlx[i], nux[i]);
                    for(int i = 0; i < m; i++)
                        printf("i: %d myplcx: %f mypuc: %f\n", i, myplc[i], mypuc[i] );*/
                    //MRQ_getchar();
                }
            }
            #endif
            
            //printf("iaia 2.6 Thread: %d\n", thnumber);
            //fflush(stdout);
            
            //std::cout << "nlp.retCode: " << nlp->retCode << " objValue: " << nlp->objValue << " dualObjValue: " << nlp->dualObjValue << " nlp.origSolverRetCode: " << nlp->origSolverRetCode << "\n";
            
            
            if( mybb->in_call_after_solving_relax_callback )
            {
                const int r = mybb->in_user_callbacks->BB_afterSolvingRelaxation( thnumber, mynode, iter, lb, ub, nlx, nux, *nlp, nlp->retCode, nlp->objValue, nlp->dualObjValue, nlp->sol, nlp->constr, nlp->dualSolC, nlp->dualSolV, pruneNode, solveAgain, branchEvenInteger );
                
                if( r != 0 )
                {
                    if( r < 0 && mybb->in_print_level > 0 ) //we assume if r is positive, user wants finish the proccess with success
                        MRQ_PRINTCALLBACKERRORNUMBER(r);
                    
                    userErrorCode = r;
                    functionReturnCode = MRQ_STOP_REQUIRED_BY_USER;
                    goto termination;
                }
                
                
                if( pruneNode )
                {
                    functionReturnCode = 0;
                    goto termination;
                }
            }
            
            
            
        }while( solveAgain );
        
        
        relaxIntSol = nlp->feasSol && MRQ_isIntegerSol(nI, intVars, nlp->sol, mybb->in_integer_tol);
        
        /*std::cout << "nlp.retCode: " << nlp->retCode << " objValue: " << nlp->objValue << " dualObjValue: " << nlp->dualObjValue << " nlp.origSolverRetCode: " << nlp->origSolverRetCode << " integer sol: " << relaxIntSol << "\n";

        for(int i = 0; i < n; i++)
            std::cout << "sol["<<i<<"]: " << nlp->sol[i] << " \t";
        //std::cout << "\n"; */
        
        
        if( relaxIntSol ) 
        {
            if( nlp->objValue < ub )
            {
                tryUpdateBestSolution(thnumber, nlp->sol, nlp->objValue, iter, true);
                ub = getUpperBound();
                
                zu_tol = MRQ_zuWithTol( ub, mybb->in_absolute_convergence_tol, mybb->in_relative_convergence_tol );
            }
            
        }
        
        
        if( pcosts )
        {
            
            if( nlp->retCode == OPT_OPTIMAL_SOLUTION && !relaxIntSol )
            {
                
                if( !pcosts->allPCostInit )
                {
                    const bool calculateEvenIfMaxComptSBranchReached = mybb->in_repeat_strong_branching_if_get_bounds_updating && iter == 1;
                    const bool multiThreads = iter == 1;
                    bool intSolFound;
                    const int nlpRetCode = nlp->retCode;
                    const double objValue = nlp->objValue;
                    const bool feasSol = nlp->feasSol;
                    double zuaux = getBestSolutionObj();
                    double *olx, *oux;
                    double *sol = auxVar;
                    double *intSol = auxVar2;
                    getVariableBoundsArrayPointers(olx, oux);
                    
                    MRQ_copyArray(n, nlp->sol, sol);
                    
                    //std::cout << "\n\n\n";
                    
                    if(mybb->in_print_level > 4)
                        MRQ_PRINTMSG("Calculating strong branching for pseudocosts\n");
                    
                    pcosts->calculateStrongBranching( nthreads, multiThreads ? nthreads : 1, *prob, nI, intVars, multiThreads, mybb->in_consider_relax_infeas_if_solver_fail, mybb->in_print_level, mybb->in_integer_tol, olx, oux, multiThreads ? nlps : &nlp, &mynode, nlx, nux, sol, nlp->objValue, calculateEvenIfMaxComptSBranchReached, updtBounds, pruneNode, zuaux, intSolFound, intSol ) ;
                    
                    
                    if( mybb->in_print_level > 2 && updtBounds )
                        MRQ_PRINTMSG("Variable bounds updated by strong branching for pseudocosts\n");
                    
                    if( iter == 1 )
                        mybb->out_number_of_strong_branching_calculations_to_pseudo_costs ++ ;
                    
                    if( intSolFound ) //if an integer solution was found, it is already better than upper bound
                    {
                        tryUpdateBestSolution( thnumber, intSol, zuaux, iter, true );
                        
                        ub = getUpperBound();
                        zu_tol = MRQ_zuWithTol( ub, mybb->in_absolute_convergence_tol, mybb->in_relative_convergence_tol );
                    }
                    
                    pcosts->checkIfAllPCostsAreCalculated(nI, intVars, olx, oux);
                    
                    //restoring values in nlp
                    nlp->retCode = nlpRetCode;
                    nlp->feasSol = feasSol;
                    MRQ_copyArray(n, sol, nlp->sol);
                    nlp->objValue = objValue;
                    
                    
                    if( pruneNode )
                    {
                        functionReturnCode = 0;
                        goto termination;
                    }
                    
                }
                
            }
            
        }
        
        
        //std::cout << "updtBounds: " << updtBounds << "\n";
        //MRQ_getchar();
        
    }while(updtBounds);


    if( iter == 1 && pcosts )
    {//unlocking other threads, so they could enter in the branch and bound loop
        unlockAllAuxiliaryThreads(); //note: this fucntion can be called more than one time because we have a loop here. It only works because BBL prevents and no try unlock mutex if it is not locked.
    }
    
    
    
    
    
    if( mybb->in_use_early_branching )
    {
        
        if( relaxIntSol && nlp->dualObjValue < zu_tol )
        {
            //we must refine the solution. 
            
            MRQ_fixIntVarsOnSolByList(nI, intVars, nlp->sol, *nlp);
            
            nlp->solve(false);
            
            relaxIntSol = nlp->feasSol;
            
            //We do not need do it actually, but we resume the bounds here...
            nlp->setnVariablesBounds(n, nlx, nux);
        }
        
    }


    if( nlp->retCode == OPT_OPTIMAL_SOLUTION )
    {
        
        if(iter > 1) //we just update pseudocosts if this is not the first iteration
        {
            if( pcosts )
            {
                
                if( mybb->in_calculate_pseudo_cost_average_above_error_estimative && mybb->in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING  )
                {
                    //here, we eval the quality of our estimative
                    const unsigned int minNumberOfPCostsCalculations = 1;
                    double objIncreasiestimativeToPseudoPruning = 0.0;
                    const double nodeLowerBound = node.getLowerBound();
                    
                    bool allVarIncreasingEstimated;
                    
                    MRQ_BasePseudoCostCalc *myPcosts = this->pcosts;
                    
                    myPcosts->getObjIncreaseEstimative(mybb->in_only_apply_pseudo_pruning_on_fixed_integer_vars, reverseIntVars, minNumberOfPCostsCalculations, node.nMyBounds, node.myBounds, mybb->in_alpha_to_balance_estimative_in_pseudo_pruning, objIncreasiestimativeToPseudoPruning, allVarIncreasingEstimated);
                    
                    if( allVarIncreasingEstimated )
                    {
                        const double diff = MRQ_max( (nodeLowerBound + objIncreasiestimativeToPseudoPruning) - nlp->objValue, 0.0  );  //ok, instead off calculate the error of the estimative, we only consider the error if the estimativa is above the real value.So we can get a better idea how pseudo costs are estimating the objective...
                        
                        if( diff > 0 )
                        {
                            pCostAvgAboveErrorEstimative[thnumber] += diff;
                            nPCostAvgAboveErrorEstimative[thnumber]++;
                        }
                    }
                    
                }
                
                pcosts->updatePCosts( nthreads, reverseIntVars, mynode.nMyBounds, mynode.myBounds, mynode.getLowerBound(), nlp->objValue, nlp->sol );
            }
            
            #if 0
            //checking if this node will be pruned by BBL. If we are using psedo prunning, we have to update the pseudocosts
            if( mybb->in_pseudo_pruning_strategy != MRQ_BB_PPS_NO_PSEUDO_PRUNING  )//&&  nodeLowerBound >= zu_tol )
            {
                pruningPcosts->updatePCosts(nthreads, reverseIntVars, mynode.nMyBounds, mynode.myBounds, mynode.getLowerBound(), nlp->objValue, nlp->sol);
                
                /*node.print();
                pruningPcosts->print(nI, intVars);
                std::cout << "nlp.retCode: " << nlp->retCode << " objValue: " << nlp->objValue << " dualObjValue: " << nlp->dualObjValue << " nlp.origSolverRetCode: " << nlp->origSolverRetCode  << " parent lower bound: " << mynode.getBestLowerBound() << " ub: " << ub << " zu_tol: " << zu_tol << "\n";
                
                MRQ_getchar();*/
                
            }
            #endif
        }
        
        nodeLowerBound = mybb->in_use_dual_obj_to_bound_prunning ? nlp->dualObjValue : nlp->objValue;
        
        optRetCode = BBL_OPTIMAL_SOLUTION;
        
        if( nlp->objValue < zu_tol && ub < MRQ_INFINITY )
        {
            if( mybb->in_fix_int_vars_from_nlp_relax_sol && nbin > 0 )
            {
                int nFixed;
                MRQ_BinVarsOptNlpRelaxSolFixer fixer;
                
                
                /*{
                    double *sol = nlp->sol;
                    double *duallx = nlp->dualSolV, *dualux = &nlp->dualSolV[n];
                    
                    for(int i = 0; i < nI; i++)
                    {
                        int ind = intVars[i];
                        
                        if( nlx[ind] != nux[ind] )
                        {
                            if( MRQ_gap( sol[ind] ) < 0.001 )
                            {
                                std::cout << ind << " - sol: " << sol[ind] << " dlx: " << duallx[ind] << " dux: " << dualux[ind] << "\n";
                                MRQ_getchar();
                            }
                        }
                    }
                }*/
                
                
                r = fixer.fixBinVarsFromNlpRelaxSol(nI, intVars, ub, nlp->objValue, nlp->dualSolV, &(nlp->dualSolV[n]), nFixed, nlx, nux, &mynode );
                MRQ_IFERRORGOTOLABEL(r, functionReturnCode, r, termination);
                
                /*if(nFixed > 0)
                {
                    std::cout << "Fixei " << nFixed << " variaveis!\n";
                    MRQ_getchar();
                }*/
                
            }
            
            
            if( mybb->in_pseudo_pruning_strategy == MRQ_BB_PPS_ON_NODE_EXPLORATION_AND_BRANCHING )
            {
                const double nodeGap = ub - nlp->objValue;
                
                if( nodeGap < mybb->in_absolute_convergence_tol_for_pseudo_pruning  ||  nodeGap/MRQ_abs(ub)  < mybb->in_relative_convergence_tol_for_pseudo_pruning  )
                {
                    bool prune;
                    const double ub_slack =  ub + mybb->in_relative_upper_bound_slack_factor_for_pseudo_pruning * MRQ_abs(ub) + mybb->in_absolute_upper_bound_slack_for_pseudo_pruning;
                    MRQ_BasePseudoCostCalc *myPcosts = this->pcosts;
                    
                    prune = myPcosts->checkIfNodeCanBePrunedByEstimative(mybb->in_only_apply_pseudo_pruning_on_fixed_integer_vars, nI, intVars, mybb->in_min_number_of_bound_prunes_per_var_before_pseudo_pruning, nlx, nux, nlp->objValue, nlp->sol, mybb->in_integer_tol, mybb->in_alpha_to_balance_estimative_in_pseudo_pruning, ub_slack);
                    
                    if(prune)
                    {
                        //pruningPcosts->print(nI, intVars);
                        //for(int k = 0; k < nI; k++)
                            //printf("x_%d: %f\t", intVars[k], nlp->sol[intVars[k]] );
                        //std::cout << "Podando no por pseudo poda no branching! nlp->objValue: " << nlp->objValue << "  ub_slack: " <<  ub_slack << "\n";
                        //MRQ_getchar();
                        
                        if( mybb->in_count_total_prunes )
                            afterExpoPseudoPruningCounter[thnumber]++;
                        
                        if( nlp->objValue < pseudoPruningLowestLowerBound[thnumber] )
                            pseudoPruningLowestLowerBound[thnumber] = nlp->objValue;
                        
                        pruneNode = true;
                        functionReturnCode = 0;
                        goto termination;
                    }
                }
            }
            
        }
    }
    else if( nlp->retCode == OPT_INFEASIBLE_PROBLEM )
    {
        nodeLowerBound = NAN;
        optRetCode = BBL_INFEASIBLE_PROBLEM;
    }
    else if( nlp->feasSol )
    {
        nodeLowerBound = nlp->dualObjValue;
        if( std::isnan(nodeLowerBound) )
            nodeLowerBound = -MRQ_INFINITY;
        
        optRetCode = BBL_FEASIBLE_SOLUTION;
    }
    else if( nlp->retCode == OPT_UNBOUNDED_PROBLEM )
    {
        MRQ_PRINTERRORMSG("Subproblem is unbouded!");
        
        pruneNode = true;
        
        if( nlp->getSolverCode() == optsolvers::OPT_MOSEK || mybb->in_consider_relax_infeas_if_solver_get_unbounded )
        {
            nodeLowerBound = NAN;
            optRetCode = BBL_INFEASIBLE_PROBLEM;
        }
        else
        {
            #if MRQ_DEBUG_MODE
                assert(iter == 1);
            #endif
            functionReturnCode = MRQ_UNBOUNDED_PROBLEM; //it will stop branch and bound...
            goto termination;
        }
    }
    else
    {
        mybb->out_nlp_failure_in_some_relaxation = true;
        nodeLowerBound = mynode.getLowerBound();
        optRetCode = BBL_UNDEFINED;
    }

    

    dualObjValue = nlp->dualObjValue;
    if( std::isnan(dualObjValue) ) //some solvers has no dualObjValue. In this case, the value is set as NaN
        dualObjValue = nlp->objValue;

    #if MRQ_SET_LBHEUR_ON_BB_NODE
        if( nodeLowerBound > mynode.heurlb )
            mynode.heurlb = nodeLowerBound;
    #endif


    if( nlp->feasSol )
    {
        objValue = nlp->objValue;
        MRQ_copyArray(n, nlp->sol, sol);
        
        if( useDualSol )
        {
            MRQ_copyArray( m, nlp->dualSolC, dualSol );
            MRQ_copyArray( 2*n, nlp->dualSolV, &dualSol[m] );
        }
        
        
        
        if( nodeLowerBound < zu_tol && !relaxIntSol )
        {
            
            if( useRoundings[thnumber] )
            {
                if( ssRoundings )
                {
                    const auto nlp_feasSol = nlp->feasSol;
                    const auto nlp_retCode = nlp->retCode;
                    const auto nlp_origSolverRetCode = nlp->origSolverRetCode;
                    const auto nlp_objValue = nlp->objValue;
                    const auto nlp_dualObjValue = nlp->dualObjValue;
                    
                    int algRetCod;
                    double outObj = INFINITY, *outSol = auxVar2;
                    
                    MRQ_Random &random = *ssrRandoms[thnumber];
                    MRQ_SSRoundingExecutor &ssRounding = *ssRoundings[thnumber];
                    
                    
                    ssRounding.run( *prob, *ssrBinSumConstrs, ccstorager, random, *nlp, thnumber, nthreads, INFINITY, myplc, mypuc, nlx, nux, nI, intVars, reverseIntVars, nC, contVars, nlp->sol, algRetCod, outObj, outSol);
                    
                    //std::cout << "algRetCod: " << algRetCod << " outObj: " << outObj << "\n";
                    //MRQ_getchar();
                    
                    //restoring the nlp atributes
                    nlp->feasSol = nlp_feasSol;
                    nlp->retCode = nlp_retCode;
                    nlp->origSolverRetCode = nlp_origSolverRetCode;
                    nlp->objValue = nlp_objValue;
                    nlp->dualObjValue = nlp_dualObjValue;
                    
                    MRQ_copyArray(n, (const double*) sol, nlp->sol ); //restoring the solution. We take advantage we have already copy the solution on so array.
                    
                    if( outObj < zu_tol  )
                    {
                        if( mybb->in_print_level > 5 )
                            std::cout << MRQ_PREPRINT "structured rounding heuristic found a feasible solution. obj: " << outObj << "\n";
                        
                        tryUpdateBestSolution(thnumber, outSol, outObj, iter, true);
                    }
                }
                else
                {
                    double objrsol;
                    //bool *auxEval = (bool *) auxInd;
                    double *rsol = auxVar2;
                    MRQ_Rounding* rounding = roundings[thnumber];
                    
                    
                    const bool r = rounding->roundSolution( thnumber, *prob, nI, intVars, mybb->in_absolute_feasibility_tol, mybb->in_relative_feasibility_tol, nlp->sol, ub, true, false, nlx, nux, auxConstr, auxVar, *nlp, rsol, objrsol );
                    
                    if( r )
                    {
                        if( mybb->in_print_level > 5 )
                            std::cout << MRQ_PREPRINT "rounding heuristic found a feasible solution. obj: " << objrsol << "\n";
                        
                        tryUpdateBestSolution(thnumber, rsol, objrsol, iter, true);
                        
                        ub = getUpperBound();
                        zu_tol = MRQ_zuWithTol( ub, mybb->in_absolute_convergence_tol, mybb->in_relative_convergence_tol );
                    }
                    
                }
                useRoundings[thnumber] = false;
            }
            
            
            
            if( useIGMA2s[thnumber] )
            {
                const int nnonbin = nI - nbin;
                bool apply = true;
                
                
                for(int i = 0; i < nnonbin; i++)
                {
                    const int ind = nonBinVars[i];
                    
                    if( nlx[ind] != nux[ind] )
                    {
                        //we have a integer nonbinary variable non fixed. So, we cannot apply igma
                        apply= false;
                        break;
                    }
                }
                
                if(apply)
                {
                    int r;//, r2;
                    
                    
                    MRQ_GapMinProb &gapminsolver = gapmins[thnumber];
                    auto *gapmin = gapminsolver.solver;
                    
                    const double maxDistance = 1.0 + ceil(mybb->in_igma2_factor_to_max_dist_constr *n);
                    
                    const int distConstIndex = m; //that only works because we are not adopting the objective cut. If we adopt the objective cut, this index is not more this value...
                    
                    if( prepocess )
                    {
                        for( int i = 0; i < m; i++ )
                        {
                            int r = gapmin->setConstraintBounds(i, myplc[i], mypuc[i]);
                            
                            #if MRQ_DEBUG_MODE
                            if(r != 0)
                            {
                                if(mybb->in_print_level > 0)
                                    MRQ_PRINTERRORNUMBER(r);
                            }
                            #endif
                        }
                    }
                    
                    
                    const bool nlp_feasSol = nlp->feasSol;
                    const int nlp_retCode = nlp->retCode;
                    const int nlp_origSolverRetCode = nlp->origSolverRetCode;
                    const double nlp_objValue = nlp->objValue;
                    const double nlp_dualObjValue = nlp->dualObjValue;
                    const double *nlp_sol = sol; //we take advantage we have already copied the solution to sol array //double *nlp_sol = auxVar;
                    //double nlp_constr[m], nlp_dualSolC[m], nlp_dualSolV[2*n];
                    
                    bool gapMinIntSol;
                    double objOutSol;
                    double *outSol = auxVar2;
                    
                    //MRQ_copyArray(n, (const double*) nlp->sol, nlp_sol);
                    /*MRQ_copyArray(m, (const double*) nlp->constr, nlp_constr);
                    MRQ_copyArray(2*n, (const double*) nlp->dualSolV, nlp_dualSolV);
                    MRQ_copyArray(m, (const double*) nlp->dualSolC, nlp_dualSolC); */
                    
                    
                    //std::cout << "indo executar igma2 iter:" << iter << "\n";
                    
                    /*for(int i = 0; i < n; i++)
                    {
                        double slx, sux;
                        
                        //std::cout << "i: " << i << " nlx: " << nlx[i] << " nux: " << nux[i] << " \t";
                        
                        nlp->getVariableBounds(i, slx, sux);
                        
                        if( (slx != nlx[i] &&  (slx > -MRQ_INFINITY && nlx[i] > -MRQ_INFINITY) )  ||  (sux != nux[i] && (sux < MRQ_INFINITY && nux[i] < MRQ_INFINITY) ) )
                        {
                            std::cout << "\nDeu diferença antes! i: " << i << " slx: " << slx << " nlx: " << nlx[i] << " sux: " << sux << " nux: " << nux[i] << "\n" ;
                            MRQ_getchar();
                        }
                    }*/
                    //std::cout << "\n";
                    
                    r = igma2Iter.run(thnumber, gapminsolver, nlp, nlx, nux, distConstIndex, maxDistance, nlp_sol, nlp->dualSolC, nlp->dualSolV, outSol, objOutSol, prepocess, ccstorager, &gapMinIntSol);
                    
                    nlp->feasSol = nlp_feasSol;
                    nlp->retCode = nlp_retCode;
                    nlp->origSolverRetCode = nlp_origSolverRetCode;
                    nlp->objValue = nlp_objValue;
                    nlp->dualObjValue = nlp_dualObjValue;
                    
                    //restoring solution
                    MRQ_copyArray(n, (const double*) nlp_sol, nlp->sol);
                    /*MRQ_copyArray(m, (const double*) nlp_constr, nlp->constr);
                    MRQ_copyArray(m, (const double*) nlp_dualSolC, nlp->dualSolC);
                    MRQ_copyArray(2*n, (const double*) nlp_dualSolV, nlp->dualSolV); */
                    
                    #if MRQ_DEBUG_IGMA2_BUG
                    {
                        mynode.igma2OnAncestral = 1;
                    }
                    #endif
                    
                    //std::cout << "apos executar igma3 iter:" << iter << "\n";
                    /*for(int i = 0; i < n; i++)
                    {
                        double slx, sux;
                        
                        //std::cout << "i: " << i << " nlx: " << nlx[i] << " nux: " << nux[i] << " \t";
                        
                        nlp->getVariableBounds(i, slx, sux);
                        
                        if( (slx != nlx[i] &&  (slx > -MRQ_INFINITY && nlx[i] > -MRQ_INFINITY) )  ||  (sux != nux[i] && (sux < MRQ_INFINITY && nux[i] < MRQ_INFINITY) ) )
                        {
                            std::cout << "\nDeu diferença depois! i: " << i << " slx: " << slx << " nlx: " << nlx[i] << " sux: " << sux << " nux: " << nux[i] << "\n" ;
                            MRQ_getchar();
                        }
                    } */
                    //std::cout << "\n";
                    
                    //maybe it is not necessary, but we restore the node bounds
                    /*r = nlp->setnVariablesBounds(n, nlx, nux);
                    if( r != 0)
                    {
                        if(mybb->in_print_level > 0)
                            MRQ_PRINTERRORNUMBER(r);
                        //we report error, but we do not stop algorithm
                    }
                    */
                    
                    if(r == MRQ_HEURISTIC_SUCCESS)
                    {
                        /*std::cout << "igma 3 encontrou solucao viavel. obj: " << objOutSol << " iter: " << iter << "\n";
                        printf("%0.16f\n", objOutSol);
                        
                        for(int i = 0; i < n; i++)
                            std::cout << "sol["<<i<<"]: " << outSol[i] << " \t"; */
                        
                        
                        //MRQ_getchar();
                        
                        tryUpdateBestSolution(thnumber, outSol, objOutSol, iter, true);
                    }
                    else if(r != MRQ_HEURISTIC_FAIL)
                    { //if we got an heuristic fail, we just continue the computation. For any other code, we stop the algorithm
                        functionReturnCode = r;
                        goto termination;
                    }
                    
                }
                
                useIGMA2s[thnumber] = false;
            }
            
            
            
            if( useHeuristics[thnumber] )
            {
                double obj;
                double *intSol = auxVar;
                MRQ_HeuristicExecutor &heurExec = heurExecs[thnumber];
                
                
                const int nAlgs = heurExec.getNumberOfAlgs();
                for(int i = 0; i < nAlgs; i++)
                {
                    MRQ_Algorithm *alg = heurExec.getAlgPointer(i);
                    
                    alg->xInit = nlp->sol; //we do it here every time when we apply heuristic because maybe some day optsolvers can realloc nlp.sol, and, so, pointer can change...
                }
                
                
                
                const int r = heurExec.insideRun( *prob, milpSolverParams, nlpSolverParams, ub, obj, intSol, true, NULL, mybb->thnumber + thnumber, 1, mybb->in_feas_heuristic_max_time, nlx, nux );
                
                if(r == MRQ_HEURISTIC_SUCCESS || r == MRQ_OPTIMAL_SOLUTION)
                {
                    tryUpdateBestSolution(thnumber, intSol, obj, iter, true);
                }
                
                useHeuristics[thnumber] = false;
            }
            
            
            
            if( useOASubs[thnumber] )
            {
                int r;
                double obj;
                double heurlb;
                MRQ_OuterApp &oaSub = oaSubs[thnumber];
                
                useOASubs[thnumber] = false;
                
                #if MRQ_SET_LBHEUR_ON_BB_NODE
                    heurlb = mynode.heurlb;
                #else
                    heurlb = mynode.getLowerBound();
                #endif
                
                oaSub.in_lower_bound = MRQ_max( lb, MRQ_max(mynode.getLowerBound(), heurlb) );
                oaSub.in_upper_bound = ub;
                
                oaSub.points[0] = nlp->sol; // in the past, we only do this setting in beforeAll and after, never more. However, user can add new variables in nlp and so, nlp.sol can change. So, by the safe, we reset it here...
                
                obj = getBestSolutionObj();
                
                if( obj < MRQ_INFINITY )
                {
                    getBestSolutionCopy( oaSub.points[1], obj );
                    
                    oaSub.nPoints = 2;
                }
                else
                {
                    oaSub.nPoints = 1;
                }
                
                
                if(mybb->in_print_level > 4)
                    MRQ_PRINTMSG("Applying Outer Approximation on a subproblem\n");
                
                r = oaSub.insideRun(*prob, mybb->in_OAmilpParams, mybb->in_OAnlpParams, mybb->thnumber + thnumber, mybb->in_outer_app_subprob_time, nlx, nux);
                
                
                if( r == MRQ_OPTIMAL_SOLUTION )
                {
                    //if everything is fine, this node will be pruned by optimallity
                    optRetCode = BBL_OPTIMAL_SOLUTION;
                    relaxIntSol = true;
                    MRQ_copyArray(n, oaSub.out_best_sol, sol);
                    objValue = dualObjValue = oaSub.out_best_obj;
                }
                else if( r == MRQ_INFEASIBLE_PROBLEM )
                {
                    objValue = dualObjValue = MRQ_INFINITY;
                }
                else
                {
                    #if MRQ_DEBUG_MODE
                        assert( oaSub.out_lower_bound < ub );
                    #endif
                    
                    #if MRQ_SET_LBHEUR_ON_BB_NODE
                        if( oaSub.out_lower_bound > ( (MRQ_NewBBNode&) node).heurlb )
                            ( (MRQ_NewBBNode&) node ).heurlb  = oaSub.out_lower_bound;
                    #else
                        if( oaSub.out_lower_bound > nodeLowerBound )
                            nodeLowerBound = oaSub.out_lower_bound;
                    #endif
                    
                    if( oaSub.out_best_obj < ub )
                    {
                        tryUpdateBestSolution(thnumber, oaSub.out_best_sol, oaSub.out_best_obj, iter, true);
                        
                        ub = getUpperBound();
                        zu_tol = MRQ_zuWithTol( ub, mybb->in_absolute_convergence_tol, mybb->in_relative_convergence_tol );
                    }
                    
                }
                
            }
            
        }
        
    }
    
    

    generalFeasibleSol = relaxIntSol;


    if( branchEvenInteger )
    {
        //if we let generalFeasibleSol as true, BBL will prune the node by otlimallity. So, we update best soltuion here
        generalFeasibleSol = false;
    }


    //applying outer approximation
    if( (int) thnumber == (nthreads - 1)/2 )
    { //only midle thread can apply outer approximation. In this way, thread 0 can go ahead and generate nodes to other thread explore them while midle thread apply outer approximation. We perform in this way to let last thread free to muriqui dcserver...
        
        
        //we just apply OA if we have a new linearization point or the last OA application could terminate an iteration with success. Note if OA get MRQ_MAX_TIME_STOP before finish the first iteration, we will apply OA several times having the same result (MRQ_MAX_TIME_STOP before finish the first iteration). So, we only apy OA again if BB found a new linearization point or last OA application could find some solution to put in the solution history, i.e, the next oa applciation will be different from previous. Note lastOANPoints is initialized with -1. So, we will enter in this if in the first time to apply OA. An alternative stregy could be increase the time to apply OA automatically here but that is little dangerous...
        
        if( iter >= iterNextOAApplic && (oaPoints->nPoints > 0 || lastOANPoints != 0)  )
        {
            double *olx, *oux;
            
            iterNextOAApplic += mybb->in_outer_app_frequence;
            
            
            getVariableBoundsArrayPointers(olx, oux);
            
            
            SEMAPH_OAPoints.lock(nthreads);
                oa->addPointsToLinearisation( oaPoints->nPoints, n, oaPoints->points );
                
                oaPoints->desallocate();
            SEMAPH_OAPoints.unlock(nthreads);
            
            
            if( mybb->in_print_level > 2 )
                MRQ_PRINTMSG("Applying Outer Approximation on the root\n");
            
            ub = getUpperBound(); //we can have updated zu
            zu_tol = MRQ_zuWithTol( ub, mybb->in_absolute_convergence_tol, mybb->in_relative_convergence_tol );
            
            oa->in_lower_bound = lb;
            oa->in_upper_bound = zu_tol;
            
            const int r = oa->insideRun(*prob, mybb->in_OAmilpParams, mybb->in_OAnlpParams, thnumber + mybb->thnumber,  mybb->in_outer_app_time, olx, oux );
            
            //std::cout << "oa ret code: " << oa->out_return_code << " obj: " << oa->out_best_obj << " Thread: " << thnumber << "\n";
            
            if( r == MRQ_OPTIMAL_SOLUTION )
            {
                //we have to find a way to avoid exploitation continue
                tryUpdateLowerBound( oa->out_upper_bound ); //we use upper bound to bbl declare optimal solution
                pruneNode = true;
            }
            else
            {
                double bestobj = getBestSolutionObj();
                
                if( r == MRQ_INFEASIBLE_PROBLEM)
                {
                    if( bestobj < MRQ_INFINITY )
                    {
                        //we have to find a way to avoid exploitation continue
                        tryUpdateLowerBound( bestobj ); //we use upper bound to bbl declare optimal solution
                        pruneNode = true;
                    }
                    else
                    {
                        functionReturnCode = MRQ_INFEASIBLE_PROBLEM;
                        goto termination;
                    }
                }
                else
                {
                    lastOANPoints = oa->out_sol_hist.getnsols();
                    
                    const int r = oa->addPointsToLinearisation( n, oa->out_sol_hist );
                    
                    if(r != 0 )
                    {
                        if(mybb->in_print_level > 0)
                            MRQ_PRINTERRORNUMBER(r);
                        //here, we continue if we get an error
                    }
                    
                    if( oa->out_lower_bound > lb )
                        tryUpdateLowerBound( oa->out_lower_bound );
                    
                    //std::cout << "oa.iters: " << oa->out_number_of_iterations << " lastOANPoints: " << lastOANPoints << "\n";
                }
                
            }
            
            
            if( oa->out_feasible_solution )
            {
                tryUpdateBestSolution(thnumber, oa->out_best_sol, oa->out_best_obj, iter, false);
                
                ub = getUpperBound();
                zu_tol = MRQ_zuWithTol( ub, mybb->in_absolute_convergence_tol, mybb->in_relative_convergence_tol );
            }
            
        }
        
    }



    #if MRQ_CHECK_VAR_BOUNDS_ON_SOLVER_SOL_IN_BB
    {
        //that is not ok, but we have to check if solver respect variable bounds... ridiculous...
        if( nlp->feasSol && !generalFeasibleSol )
        {
            const double intTol = mybb->in_integer_tol;
            bool change = false;
            
            
            for(int i = 0; i < nI; i++)
            {
                const int ind = intVars[i];
                
                if( sol[ind] > nux[ind] + intTol )
                {
                    sol[ind] = nux[ind];
                    change = true;
                    //MRQ_getchar();
                }
                else if( sol[ind] < nlx[ind] - intTol )
                {
                    sol[ind] = nlx[ind];
                    change = true;
                    //MRQ_getchar();
                }
            }
            
            
            //considere a possibilidade fazer aqui um código que fixe variaveis inteiras e resolva o problema novamente...
            
            if(change)
            {
                /*std::cout << "Precisei mudar a solucao!\n";
                MRQ_getchar();
                
                for(int i = 0; i < nI; i++)
                {
                    const int ind = intVars[i];
                    
                    std::cout << "nlpx["<<ind<<"]: " << nlp->sol[ ind ] << " myx["<<ind<<"]: " << sol[ind] << "\n";
                } */
                
                
                retCode = BBL_FEASIBLE_SOLUTION;// we put that to allow branching even on integer solution...
                nlp->retCode = optsolvers::OPT_UNDEFINED;
                
                bool relaxIntSol = MRQ_isIntegerSol(nI, intVars, sol, mybb->in_integer_tol);
                
                
                if( relaxIntSol )
                {
                    
                    bool feas;
                    bool *auxEval = (bool *) auxInd;
                    
                    
                    MRQ_setAllArray(m, auxEval, true);
                    
                    prob->isFeasibleToConstraints(thnumber, sol, true, auxEval, mybb->in_absolute_feasibility_tol, mybb->in_relative_feasibility_tol, feas, nlp->constr );
                    
                    if( feas )
                    {
                        double obj;
                        
                        const int r = prob->objEval(thnumber, false, sol, obj);
                        
                        if( r == 0 && obj < ub )
                        {
                            objValue = obj;
                            tryUpdateBestSolution(thnumber, sol, obj, iter, true);
                            
                            //that is the last ub updation. we do not need ub more, so we comment it. If some day you put more precedures after that in this method, uncomment this...
                            //ub = getUpperBound();
                            //zu_tol = MRQ_zuWithTol( ub, mybb->in_absolute_convergence_tol, mybb->in_relative_convergence_tol );
                        }
                        
                        generalFeasibleSol = true;
                    }
                    
                }
                
            }
        }
    }
    #endif

    #if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        std::ostream &thOut = *MRQ_thsOut[std::this_thread::get_id()];
        thOut << "\tMRQ_BBLCallbacks::solveSubProblem 10" << std::endl;
    }
    #endif

    branchStrategy = origBBLBranchStrategy; //we can have already change the branchstarategy in the previous iteration... we have to check if we must put BBL_BS_USER_NODE_GENERATION to constraint branching

    if( constrBranchStrat != MRQ_BB_CBS_NO_CONSTRAINT_BRANCH && binSumConstrs.nbinSumConstrs > 0 && origBBLBranchStrategy != BBL_BS_USER_NODE_GENERATION )
    {
        if( mybb->in_stop_multibranch_after_first_bound_prune &&  getSomeBoundPrune() )
        {
            constrBranchStrat = MRQ_BB_CBS_NO_CONSTRAINT_BRANCH;
            
            if( mybb->in_print_level > 4 )
                MRQ_PRINTMSG("Fisrt prune by bound detected. Avoiding multibranching\n");
            //MRQ_getchar();
        }
        else
        {
            if( !pruneNode && objValue < ub && optRetCode != BBL_INFEASIBLE_PROBLEM )
            {
                MRQ_BinSumConstrsChooser &constrChooser = constrChoosers[thnumber];
                const double *plc = prepocess ? tplc[thnumber] : prob->lc;
                const double *puc = prepocess ? &plc[m] : prob->uc;
                
                
                if( constrChooser.calculateCandidateConstrsToBranch( binSumConstrs.nbinSumConstrs, binSumConstrs.binSumConstrs, nlx, nux, prob->xtype, prob->A, plc, puc, mybb->in_integer_tol, sol ) > 0 )
                    branchStrategy = BBL_BS_USER_NODE_GENERATION;
            }
            
        }
        
    }

    //std::cout << "retCode: " << retCode << " objValue: " << objValue << " dualObjValue: " << dualObjValue << " generalFeasibleSol: " << generalFeasibleSol << "\n";

    //for(int i = 0; i < n; i++)
        //std::cout << "sol["<<i<<"]: " << sol[i] << " \t";

    //std::cout << "\n"; 

    functionReturnCode = 0;

    #if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        std::ostream &thOut = *MRQ_thsOut[std::this_thread::get_id()];
        thOut << "\tMRQ_BBLCallbacks::solveSubProblem 11" << std::endl;
    }
    #endif


    termination:

    if( iter == 1 && pcosts )
    {
        //TODO: check if we can erase the previous call and only call unlockAllAuxiliaryThreads here
        
        //unlocking other threads, so they could enter in the branch and bound loop
        unlockAllAuxiliaryThreads(); //note: this fucntion can be called more than one time because we could call it before. It only works because BBL prevents and no try unlock mutex if it is not locked.
    }


    #if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Leaving MRQ_BBLCallbacks::solveSubProblem" << std::endl;
        
        //we already know the current bug which we are looking for is in this function. So, when we reach this point, we know this executions was succesfull and we can close the debug file to open it as empty file in the next execution:
        MRQ_closeFileToThread(tid);
    }
    #endif


    return functionReturnCode;
}



int MRQ_BBLCallbacks::chooseIndexToBranch( const int thnumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, BBL_RETURN_CODES retCode, double objValue, double dualObjValue, double *sol, double *dualSol, unsigned int &sizeIndices, unsigned int *indices, double *breakValues1, double *breakValues2, BBL_Node* &nodes)
{
    /*#if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Entering at MRQ_BBLCallbacks::chooseIndexToBranch" << std::endl;
    }
    #endif */



    MRQ_NewBBNode &mynode = (MRQ_NewBBNode&) node;

    bool branchOnInt = false;
    int nBranchVars;
    MRQ_NLPSolver *nlp = nlps[thnumber];
    const double *nlpsol = sol;
    double *auxVar = auxVars[thnumber];


    if( mybb->in_branching_strategy == MRQ_BB_BS_USER_INDEX_CHOICE )
    {
        const int r = mybb->in_user_callbacks->BB_chooseIndexToBranch( thnumber, mynode, iter, lb, ub, nlx, nux, nlp->retCode, nlp->objValue, nlpsol, nlp->constr, nlp->dualSolC, nlp->dualSolV, nBranchVars, indices );
        
        if(r != 0)
        {
            if(mybb->in_print_level > 0)
                MRQ_PRINTCALLBACKERRORNUMBER(r);
            
            userErrorCode = r;
            return r;
        }
        
    }
    else
    {
        const int maxBranchVars = mybb->in_number_of_branching_vars;
        MRQ_NewChooseIndexToBranch &chooseInd = chooseIndices[thnumber];
        
        chooseInd.chooseIndices( pcosts ? pcosts->pcost : NULL, mybb->in_branching_strategy, mybb->in_integer_tol, nI, intVars, reverseIntVars, nbin, binVars, nonBinVars, mybb->in_pseudo_cost_mu, maxBranchVars, nlpsol, nlp->retCode == OPT_OPTIMAL_SOLUTION, nlx, nux, prob->xprior, auxVar, nBranchVars, indices, branchOnInt );
        
    }

    sizeIndices = nBranchVars;

    if(branchOnInt == false)
    {
        for(unsigned int i = 0; i < sizeIndices; i++)
        {
            const int k = indices[i];
            
            //that is so strange, but we can have cases when solver give us solution out of box constraints, and the gap is greater than integer tolerance. So, we let branch on this values
            if( nlx[k] <= nlpsol[k] && nlpsol[k] <= nux[k] )
            {
                breakValues1[i] = floor( nlpsol[k] );
                breakValues2[i] = breakValues1[i] + 1;
            }
            else if( nlpsol[k] < nlx[k])
            {
                #if MRQ_DEBUG_MODE
                    assert( nlx[k] != nux[k] );
                #endif
                
                breakValues1[i] = nlx[k];
                breakValues2[i] = nlx[k] + 1;
            }
            else if( nlpsol[k] > nux[k] )
            {
                #if MRQ_DEBUG_MODE
                    assert( nlx[k] != nux[k] );
                #endif
                
                breakValues1[i] = nux[k] - 1;
                breakValues2[i] = nux[k];
            }
            else
            {
                assert(false);
            }
            
            //std::cout << "branch - index: " << indices[i] << " sol: " << nlpsol[indices[i]] << " lx: " << nlx[indices[i]] << " bv1: " << breakValues1[i] << " bv2: " << breakValues2[i]  << " ux: " << nux[indices[i]] << "\n";
        }
    }
    else
    {
        for(unsigned int i = 0; i < sizeIndices; i++)
        {
            const int ind = indices[i];
            
            const double solind = round( nlpsol[ind] ); //nlpsol[ind] is already an integer value, or almost an integer value
            
            if( solind == nux[ind] )
            {
                breakValues1[i] = solind - 1;
                breakValues2[i] = solind;
            }
            else
            {
                breakValues1[i] = solind;
                breakValues2[i] = solind + 1;
            }
            
            #if MRQ_DEBUG_MODE
                assert( nlx[ind] < nux[ind] ); //variable could not be fixed...
            #endif
        }
    }


    {
        const unsigned int nnodes = 1u << sizeIndices;
        int code;
        
        
        #if MRQ_SET_LBHEUR_ON_BB_NODE
            double heurlb = MRQ_max(objValue, node.getLowerBound() ); //we already updated node.lb in solveSubProblem with the relaxation value...
            if( ((MRQ_NewBBNode&) node).heurlb  > heurlb )
                heurlb = ((MRQ_NewBBNode&) node).heurlb ;
        #endif
        
        #if MRQ_DEBUG_IGMA2_BUG
            unsigned int igma2OnAncestral = ((MRQ_NewBBNode&) node).igma2OnAncestral;
            if( igma2OnAncestral > 0 )
                igma2OnAncestral++;
        #endif
        
        
        MRQ_NewBBNode *auxNodes, *first = new (std::nothrow) MRQ_NewBBNode(parentNodeBoundsStrategy, prob->n);
        
        auxNodes = first;
        
        if( !first )
        {
            if(mybb->in_print_level > 0)
                MRQ_PRINTMEMERROR;
            
            code = MRQ_MEMORY_ERROR;
            goto termination;
        }
        
        
        
        for( unsigned int i = 0;  ;  )
        {
            #if MRQ_SET_LBHEUR_ON_BB_NODE
                auxNodes->heurlb = heurlb;
            #else
                
            #endif
            //auxNodes->nBranchVars = nBranchVars;
            
            #if MRQ_DEBUG_IGMA2_BUG
                auxNodes->igma2OnAncestral = igma2OnAncestral;
            #endif
            
            i++;
            if( i >= nnodes )
                break;
            
            auxNodes->next = new (std::nothrow) MRQ_NewBBNode(parentNodeBoundsStrategy, prob->n);
            
            if( !auxNodes->next )
            {
                if(mybb->in_print_level > 0)
                    MRQ_PRINTMEMERROR;
                
                code = MRQ_MEMORY_ERROR;
                goto termination;
            }
            
            auxNodes->next->previous = auxNodes;
            
            auxNodes = (MRQ_NewBBNode *) auxNodes->next;
            
            //auxNodes->heurlb = heurlb;
        }
        
        nodes = first;
        
        code = 0;
        
        
    termination:
        
        if(code != 0)
        {
            if( auxNodes )
            {
                while( true )
                {
                    if( auxNodes->previous )
                    {
                        auxNodes = (MRQ_NewBBNode *) auxNodes->previous;
                        
                        delete auxNodes->next;
                    }
                    else
                    {
                        delete auxNodes;
                        break;
                    }
                }
            }
            
            return code;
        }
    }


    //std::cout << "Numero de variaveis de branching: " << sizeIndices << "\n";
    //for(unsigned int i = 0; i < sizeIndices; i++)
        //std::cout << "indices["<<i<<"]: " << indices[i] << " bvalue1[" << i <<"]: " << breakValues1[i] << " bvalue2[" << i <<"]: " << breakValues2[i] << "\n";


    /*#if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        //just main thread pass by here. Anyway, we create a file 
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Leaving MRQ_BBLCallbacks::chooseIndexToBranch" << std::endl;
    }
    #endif */


    return 0;
}



int MRQ_BBLCallbacks::generateNodes(const int thnumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, BBL_RETURN_CODES retCode, const double objValue, double *sol, double *dualSol, BBL_UserNodeGenerator &userNodeGenerator )
{
    /*#if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Entering at MRQ_BBLCallbacks::generateNodes" << std::endl;
        
    }
    #endif */

    MRQ_NewBBNode &mynode = (MRQ_NewBBNode&) node;
    unsigned int *varsSelected = NULL;
    BBL_NodeBoundsSol *zeroNodeBounds = nullptr;
    int rcode;
    
    
    if( mybb->in_branching_strategy == MRQ_BB_BS_USER_NODE_GENERATION )
    {
        double *olx, *oux;
        MRQ_NewUserNodeGenerator2 &userNodeGen =  userNodeGens[thnumber];
        MRQ_NLPSolver *nlp = nlps[thnumber];
        
        
        getVariableBoundsArrayPointers(olx,oux);
        
        userNodeGen.initialize(&mynode, &userNodeGenerator);
        
        const int r = mybb->in_user_callbacks->BB_generateNodes( thnumber, mynode, iter, lb, ub, nlx, nux, nlp->retCode, nlp->objValue, nlp->sol, nlp->constr, nlp->dualSolC, nlp->dualSolV, userNodeGen );
        
        if(r != 0)
        {
            if(mybb->in_print_level > 0)
                MRQ_PRINTCALLBACKERRORNUMBER(r);
            
            userErrorCode = r;
            rcode = r;
            goto termination;
        }
        
    }
    else //we are branching on constraints...
    {
        #if MRQ_DEBUG_MODE
            assert( constrBranchStrat != MRQ_BB_CBS_NO_CONSTRAINT_BRANCH || nthreads > 1 ); //another thread can already have changed constrBranchStrat
        #endif
        
        int *varsToFixOne = auxInds[thnumber];
        double *auxVar = this->auxVars[thnumber];
        
        const int m = prob->m;
        //const double ONE = 1.0;
        const unsigned int maxVarsToBranching = MRQ_max(2u, mybb->in_max_number_of_branchings_in_constraint_branching);
        
        
        int index;
        int r;
        
        #if MRQ_SET_LBHEUR_ON_BB_NODE
            const double heurlb = MRQ_max(objValue, MRQ_max( node.getLowerBound(), ( (MRQ_NewBBNode&) node ).heurlb ) ); //we already updated node.lb in solveSubProblem with the relaxation value...
        #endif
        
        //plc and puc pointers should point to same array was passed in binSumConstrs.calculateIndices
        const double *plc = oplc ? oplc : prob->lc;
        const double *puc = oplc ? &oplc[m] : prob->uc; 
        
        const minlpproblem::MIP_SparseMatrix &A = prob->A;
        BBL_NodeBoundsSol nodeBounds;
        MRQ_BinSumConstrsChooser &constrChooser = constrChoosers[thnumber];
        
        //here, it is better use the original mybb->_in_constr_branching_strategy because other thread can have change constrBranchStrat to new value before this thread take it new value in account
        constrChooser.chooseIndexToBranch( A, mybb->in_integer_tol, sol, mybb->in_constr_branching_strategy, reverseIntVars, mybb->in_pseudo_cost_mu, pcosts ? pcosts->pcost : NULL, index ); //we already call calculateCandidateConstrsToBranch in solveProblem method...
        
        /*{
            MRQ_NewPseudoCost *pcost = pcosts->pcost;
            
            for(int i = 0; i < nI; i++)
            {
                printf("i: %d intVars: %d pl: %f npl: %d nrealPl: %d pr: %f npr: %d nrealPr: %d \n", i, intVars[i], pcost[i].pl, (int) pcost[i].nPl, (int) pcost[i].nrealPl, pcost[i].pr, (int) pcost[i].nPr, (int) pcost[i].nrealPr  );
            }
        }
        
        std::cout << "index of constraint to branch: " << index << "\n";
        
        MRQ_getchar();*/
        
        #if MRQ_DEBUG_MODE
            assert( index >= 0 );
        #endif
        
        
        
        
        //minlpproblem::MIP_SparseRow &row = A[index];
        //unsigned int maxnodes = row.getNumberOfElements();
        const int *acols = A[index];
        const double *avalues = A(index);
        const unsigned int maxvars = A.getNumberOfElementsAtRow(index);
        
        const auto xtype = prob->xtype;
        
        unsigned int nNegCoef = 0, nVarsToFixOne = 0;
        unsigned int negCoefIndex = -1;
        double fixed = 0.0;
        
        for(unsigned int k = 0; k < maxvars; k++)
        {
            const int col = acols[k];
            const double val = avalues[k];
            
            if( nlx[col] == nux[col] )
            { //variable is already fixed
                fixed += nux[col] * val;
                continue;
            }
            #if MRQ_DEBUG_MODE
            else
            {
                assert( val == -1.0 || val == 1.0 );
            }
            #endif
            
            
            #if MRQ_DEBUG_MODE
                assert( minlpproblem::MIP_isIntegerType(xtype[col]) );
                assert( nlx[col] > -1.0  &&  nux[col] < 2.0 );
            #endif
            
            
            
            if(val == -1.0)
            {
                #if MRQ_DEBUG_MODE
                    assert(nNegCoef == 0); //we can have at most one varibale having negative coefficient in thisconstraint
                #endif
                nNegCoef++;
                negCoefIndex = col;
            }
            else
            {
                varsToFixOne[nVarsToFixOne] = col;
                nVarsToFixOne++;
            }
            
        }
        
        const double rhs = puc[index] - fixed;
        const double lhs = plc[index] - fixed;
        
        
        #if MRQ_DEBUG_MODE
            assert( rhs == 0.0 || rhs == 1.0 );
            assert( lhs <= 1.0 );
        #endif
        
        
        
        if( nVarsToFixOne <= maxVarsToBranching )
        {
            if( plc[index] != puc[index] ) //we have lower equal constraint
            { /*we have a constraint in one of these formats:
                
                y_{k_1} + y{k_2} + ... + y_{k_n-1} <= 1         (1)
                
                y_{k_1} + y{k_2} + ... + y_{k_n-1} <=  y{k_n}   (2) */
                
                #if MRQ_DEBUG_MODE
                    assert( (rhs == 1.0 && lhs <= 0.0 && nNegCoef == 0) || (rhs == 0.0 && lhs < 0.0  && nNegCoef == 1) );
                #endif
                
                
                //in this case, we need generate a node fixing all nonfixed variables having coefficient 1.0 at zero in this constraints. We will left a possible variable having coefficient -1.0 free to take any value (maybe it can be fixed in a future branchif if it get a fractional value in continuous relaxation).
                
                //BBL_NodeBoundsSol zeroNodeBounds[nVarsToFixOne];
                MRQ_malloc(zeroNodeBounds, nVarsToFixOne);
                MRQ_IFMEMERRORGOTOLABEL(!zeroNodeBounds, rcode, termination);
                
                for(unsigned int i = 0; i < nVarsToFixOne; i++)
                {
                    const unsigned int col = varsToFixOne[i];
                    
                    #if MRQ_DEBUG_MODE
                        assert( nlx[col] != nux[col] );
                    #endif
                    
                    //note we fix all (binarie) variables having coefficiente 1.0 to zero. We will left a possible variable having coefficient -1.0 free to take any value (maybe it can be fixed in a future branchif it get a fractional value in continuous relaxation).
                    zeroNodeBounds[i].ind = col;
                    zeroNodeBounds[i].l = 0.0;
                    zeroNodeBounds[i].u = 0.0;
                    zeroNodeBounds[i].sol = sol[col];
                }
                
                
                MRQ_NewBBNode *auxNode = new (std::nothrow) MRQ_NewBBNode(parentNodeBoundsStrategy, prob->n);
                
                MRQ_IFMEMERRORGOTOLABEL(!auxNode, rcode, termination);
                
                #if MRQ_SET_LBHEUR_ON_BB_NODE
                    auxNode->heurlb = heurlb;
                #endif
                
                r = userNodeGenerator.generateNode(nVarsToFixOne, zeroNodeBounds, true, false, auxNode );
                
                MRQ_IFERRORGOTOLABEL(r, rcode, MRQ_MEMORY_ERROR, termination);
            }
            else if(rhs == 0 && lhs == 0) //( (plc[index] == 0.0) && (puc[index] == 0.0) )
            {
                /*we (probably, :) ) have a constraint in the form:
                    y_{k_1} + y{k_2} + ... + y_{k_n-1} = y{k_n}  */
                
                #if MRQ_DEBUG_MODE
                    assert( nNegCoef == 1 );
                #endif
                
                
                //we generate a node fixing all variables to zero
                //it is enough fix the variable having the negative signal to zero
            
                BBL_NodeBoundsSol zeroNodeBounds;
            
                zeroNodeBounds.l = 0.0;
                zeroNodeBounds.u = 0.0;
                zeroNodeBounds.ind = negCoefIndex;
                zeroNodeBounds.sol = sol[negCoefIndex];
                
                
                MRQ_NewBBNode *auxNode = new (std::nothrow) MRQ_NewBBNode(parentNodeBoundsStrategy, prob->n);
                
                MRQ_IFMEMERRORGOTOLABEL(!auxNode, rcode, termination);
                
                
                #if MRQ_SET_LBHEUR_ON_BB_NODE
                    auxNode->heurlb = heurlb;
                #endif
                
                r = userNodeGenerator.generateNode(1, &zeroNodeBounds, true, true, auxNode );
                
                MRQ_IFERRORGOTOLABEL(r, rcode, MRQ_MEMORY_ERROR, termination);
            }
            
        }
        else
        { //so, we have to choose just maxVarsToBranching variables to fix on 1. We will generate a extra node fixing these variables to zero
            MRQ_NewChooseIndexToBranch &chooseInd = chooseIndices[thnumber];
            MRQ_NLPSolver *nlp = nlps[thnumber];
            const double *nlpsol = sol;
            
            bool branchOnInt;
            int nVarsSelected;
            
            
            MRQ_malloc(varsSelected, maxVarsToBranching);
            MRQ_IFMEMERRORGOTOLABEL(!varsSelected, rcode, termination);
            
            /*we taking advantage chooseInd to choose the subsets of indices in the constrant choosen to branching. To do it, we pass nVarsToFix and varsToFix replacing nI and intVars and -INFINITY as integer tolerance because we must choose maxVarsToBranching even if some of these varoables gets integer values in the current relaxation */
            chooseInd.chooseIndices( pcosts->pcost, mybb->in_variable_choosing_strategy_to_constraint_branching, -INFINITY, nVarsToFixOne, varsToFixOne, reverseIntVars, nbin, binVars, nonBinVars, mybb->in_pseudo_cost_mu_to_variable_choosing_on_constraint_branching, maxVarsToBranching, nlpsol, nlp->retCode == OPT_OPTIMAL_SOLUTION, nlx, nux, prob->xprior, auxVar, nVarsSelected, varsSelected, branchOnInt ); //we pass infinity as inttol because we have to allow variables having integer values being selected to branching
            
            #if MRQ_DEBUG_MODE
                assert( nVarsSelected == (int) maxVarsToBranching );
            #endif
            
            //BBL_NodeBoundsSol zeroNodeBounds[nVarsSelected];
            
            MRQ_malloc(zeroNodeBounds, nVarsSelected);
            MRQ_IFMEMERRORGOTOLABEL(!zeroNodeBounds, rcode, termination);
                
            //Since we will not generate a new node fixing each noxifed variables to one, we must generate node fixing selected variables to zero
            for(int i = 0; i < nVarsSelected; i++)
            {
                const unsigned int col = varsSelected[i];
                
                zeroNodeBounds[i].ind = col;
                zeroNodeBounds[i].l = 0.0;
                zeroNodeBounds[i].u = 0.0;
                zeroNodeBounds[i].sol = sol[col];
            }
            
            MRQ_NewBBNode *auxNode = new (std::nothrow) MRQ_NewBBNode(parentNodeBoundsStrategy, prob->n);
            MRQ_IFMEMERRORGOTOLABEL(!auxNode, rcode, termination);
                
            #if MRQ_SET_LBHEUR_ON_BB_NODE
                auxNode->heurlb = heurlb;
            #endif
                
            r = userNodeGenerator.generateNode(nVarsSelected, zeroNodeBounds, true, false, auxNode );
            MRQ_IFERRORGOTOLABEL(r, rcode, MRQ_MEMORY_ERROR, termination);
            
            
            //copying varsSelected to varsToFix to free varsSelected
            MRQ_copyArray(nVarsSelected, varsSelected, varsToFixOne);
            
            nVarsToFixOne = nVarsSelected;
            
            free(varsSelected); varsSelected = NULL;
        }
        
        
        
        //now, we generate one new node for each nonfixed binary variable having coeficient 1.0 in varsToFixOne. Note, we update this array when are are branching just a subset of free variables in this constraint
        {
            
            nodeBounds.l = 1.0;
            nodeBounds.u = 1.0;
            
            for(unsigned int i = 0; i < nVarsToFixOne; i++)
            {
                const unsigned int col = varsToFixOne[i];
                
                //note here, we do not have to encadeate nodes, the encadeating will be done by userNodeGenerator
                MRQ_NewBBNode *auxNode = new (std::nothrow) MRQ_NewBBNode(parentNodeBoundsStrategy, prob->n);
                MRQ_IFMEMERRORGOTOLABEL(!auxNode, rcode, termination);
                
                
                #if MRQ_SET_LBHEUR_ON_BB_NODE
                    auxNode->heurlb = heurlb;
                #endif
                
                nodeBounds.ind = col;
                nodeBounds.sol = sol[col];
                
                
                r = userNodeGenerator.generateNode(1, &nodeBounds, true, true, auxNode );
                MRQ_IFERRORGOTOLABEL(r, rcode, MRQ_MEMORY_ERROR, termination);
            }
            
        }
        
        numberOfConstrBranchings[thnumber]++;
        
    }


    rcode = 0;
    
termination:

    if(varsSelected)    free(varsSelected);
    if(zeroNodeBounds)  free(zeroNodeBounds);
    
    /*#if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Leaving MRQ_BBLCallbacks::generateNodes" << std::endl;
    }
    #endif */

    return rcode;
}



int MRQ_BBLCallbacks::endOfIteration(const int thnumber, const long unsigned int iter, const double cpuTime, const double wallTime, const double lb, const double ub, BBL_Node &node, const double *nlx, const double *nux)
{
    /*#if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Entering at MRQ_BBLCallbacks::endOfIteration" << std::endl;
    }
    #endif */
    
    
    //we assume if we enter here, if because we have to call endOfIteration
    
    
    if( mybb->in_user_callbacks )
    {
        const int r = mybb->in_user_callbacks->endOfIteration( MRQ_BB_ALG, thnumber, iter, cpuTime, wallTime, lb, ub);
        
        if( r != 0 )
        {
            if(mybb->in_print_level > 0)
                MRQ_PRINTCALLBACKERRORNUMBER(r);
            
            userErrorCode = r;
            return MRQ_STOP_REQUIRED_BY_USER;
        }
    }
    else
    {
        #if MRQ_DEBUG_MODE
        if( mybb->in_print_level > 5 )
        {
            printOpenNodesList();
            MRQ_getchar();
        }
        #endif
    }
    
    /*#if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Leaving MRQ_BBLCallbacks::endOfIteration" << std::endl;
    }
    #endif */
    
    return 0;
}


int MRQ_BBLCallbacks::updatingBestSolution(const int thnumber, double* sol, double &objValue, const double ub, const long unsigned int iter)
{
    return mybb->in_user_callbacks->updatingBestSolution( MRQ_BB_ALG, thnumber, sol, objValue, ub, iter );
}



void MRQ_BBLCallbacks::newBestSolution( const int thnumber, const double *newSol, const double oldBestObj, const double newBestObj, const long unsigned int iter )
{
    return mybb->in_user_callbacks->newBestSolution( thnumber, newSol, oldBestObj, newBestObj, iter);
}




void MRQ_BBLCallbacks::afterAll(const long unsigned int iters, const double cpuTime, const double wallTime, const double lb, const double ub)
{
    #if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Entering at MRQ_BBLCallbacks::afterAll" << std::endl;
    }
    #endif
    //desallocate(); //do dot deallocate memory here
    
    #if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "Leaving MRQ_BBLCallbacks::afterAll" << std::endl;
    }
    #endif
}






bool MRQ_BBLCallbacks::tryUpdateBestSolution(const int threadNumber, double* solution, const double fsolution, const long unsigned int iter, const bool addSolToOALin)
{
    #if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "\tEntering at MRQ_BBLCallbacks::tryUpdateBestSolution" << std::endl;
    }
    #endif
    
    const bool r = BBL_UserCallbacks::tryUpdateBestSolution( threadNumber, solution, fsolution, iter);
    
    if( r && addSolToOALin && oa )
    {
        const int n = prob->n;
        
        SEMAPH_OAPoints.lock(nthreads);
            oaPoints->addPoints(1, n, &solution);
        SEMAPH_OAPoints.unlock(nthreads);
    }
    
    
    #if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        MRQ_createFileToThread( tid );
        std::ostream &thOut = *MRQ_thsOut[tid];
        
        thOut << "\tLeaving at MRQ_BBLCallbacks::tryUpdateBestSolution" << std::endl;
    }
    #endif
    
    return r;
}



MRQ_Preprocessor* MRQ_BBLCallbacks::getPreprocessorPointer( unsigned int threadNumber)
{
    return preprocessors ? &preprocessors[threadNumber] : nullptr;
}




