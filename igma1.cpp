/*
* That file implements the Integrality Gap Minimization Algorithm, version 2, by Wendel Melo, Marcia Fampa and Fernanda Raupp
* 
* Author: Wendel Melo
* 
* Date: 01-June-2015
* 
* */

#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <climits>

#include <iostream>
#include <new>


#include "MRQ_igma1.hpp"


#if OPT_HAVE_IPOPT
    #include "IpIpoptData.hpp"
    #include "IpIpoptCalculatedQuantities.hpp"
    #include "IpTNLPAdapter.hpp"
    #include "IpOrigIpoptNLP.hpp"
    
    //#include "IpRestoIpoptNLP.hpp" 
    
#endif


using namespace branchAndBound;
using namespace optsolvers;
using namespace muriqui;






#if OPT_HAVE_IPOPT



MRQ_IpoptIntermediateCallback::MRQ_IpoptIntermediateCallback()
{
    initialize();
}

MRQ_IpoptIntermediateCallback::~MRQ_IpoptIntermediateCallback()
{
    desallocate();
}

void MRQ_IpoptIntermediateCallback::initialize(const int thnumber, const MRQ_MINLPProb *prob, const double integerTol, const double absFeasTol, const double relFeasTol, const int nI, const int* intVars)
{
    this->thnumber = thnumber;
    this->prob = prob;
    this->integer_tol = integerTol;
    this->nI = nI;
    this->intVars = intVars;
    this->absFeasTol = absFeasTol;
    this->relFeasTol = relFeasTol;
    
    this->minimumItersToabortIfNoprogress = 200;
    this->minImprovToObj =integerTol;
    this->minImprovToInfeas = 0.1;
    
    this->ipopt = NULL;
    this->sol = NULL;
    this->constrEval = NULL;
    this->constrs = NULL;
}

int MRQ_IpoptIntermediateCallback::allocate(const int n, const int m)
{
    constrEval = (bool*) malloc(m * sizeof(bool));
    sol = (double*) malloc(n * sizeof(double));
    constrs = (double *) malloc(m * sizeof(double));
    if(!constrEval || !sol || !constrs)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTMEMERROR;
        #endif
        return MRQ_MEMORY_ERROR;
    }
    
    MRQ_setAllArray(m, constrEval, true);
    
    return 0;
}

void MRQ_IpoptIntermediateCallback::desallocate()
{
    MRQ_secFree(sol);
    MRQ_secFree(constrEval);
    MRQ_secFree(constrs);
}


/* That callback is to be used in the Integrality gap minimization process with Ipopt. We just check if we can stop since we have a integer feasible solution
    *
    */ 
bool MRQ_IpoptIntermediateCallback::intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value, Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    //const double MAX_CONSTR_VIOL = 1e-8;
    //double constr_viol = ip_cq->curr_nlp_constraint_violation( Ipopt::NORM_MAX);
    
    
    //std::cout << "ipopt int call iter: " << iter << " mode: " << mode << " obj_value: " << obj_value << " inf_pr: " << inf_pr << "\n";
    
    
    
    //getting the solution. from http://www.coin-or.org/Ipopt/documentation/node23.html;
    if( mode == Ipopt::RegularMode )
    {
        Ipopt::TNLPAdapter* tnlp_adapter = NULL;
        if( ip_cq != NULL )
        { 
            Ipopt::OrigIpoptNLP* orignlp = NULL;
            orignlp = dynamic_cast<Ipopt::OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
            
            if( orignlp != NULL )
                tnlp_adapter = dynamic_cast<Ipopt::TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
        }
        
        //If ipopt is in the restoration phase, tnlp_adapter will be NULL
        
        if( tnlp_adapter && obj_value <= nI*integer_tol )// && constr_viol <= MAX_CONSTR_VIOL )
        {
            tnlp_adapter->ResortX(*ip_data->curr()->x(), sol);
            
            //now, sol has the solution
            
            //for(int i = 0; i < n; i++)
                //std::cout << "x" << i << ": " << sol[i] << "  ";
            //std::cout << "\n";
            
            
            bool intSol = MRQ_isIntegerSol(nI, intVars, sol, integer_tol);
            
            if(intSol)
            {
                bool feasible = false;
                
                int r = prob->isFeasibleToConstraints(thnumber, sol, true, constrEval, absFeasTol, relFeasTol, feasible, constrs );
                
                if(r == 0 && feasible)
                {
                    
                    double realObj, zucut, lb;
                    
                    if(ipopt->prob.m > prob->m) //so, we have the objective cut...
                    {
                        int r = ipopt->getConstraintBounds(prob->m, lb, zucut); //constraint having index m is the objective cut
                        
                        if(r != 0)
                        {
                            #if MRQ_DEBUG_MODE
                                MRQ_PRINTERRORNUMBER(r);
                            #endif
                            assert(false);
                        }
                        
                        r = prob->objEval(thnumber, !prob->hasNlConstrs, sol, realObj);
                        
                        feasible = r == 0 && realObj <= zucut;
                        
                        //std::cout << "realObj: " << realObj << "\n";
                        //MRQ_getchar();
                    }
                    
                    if(feasible)
                    {
                        //std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Interrompendo ipopt integer_tol: " << integer_tol << " absFeasTol: " << absFeasTol << " relFeasTol: " << relFeasTol << "\n";
                        
                        //for(int i = 0; i < prob->n; i++)
                            //printf("x[%d]: %0.12f\n", i, sol[i]);
                        
                        //for(int i = 0; i < prob->m; i++)
                            //printf("g[%d]: %0.12f\n", i, constrs[i]);
                        
                        
                        //MRQ_getchar();
                        return false; //requiring ipopt to stop
                    }
                }
            }
        }
    }
    
    
    if( iter >= minimumItersToabortIfNoprogress )
    {
        if( MRQ_abs(obj_value - lastObj) < minImprovToInfeas && MRQ_abs(inf_pr - lastInf_pr) < minImprovToInfeas )
        {
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$Requisitando ipopt a parar por falta de progresso\n";
            
            
            //MRQ_getchar();
            return false; //requiring ipopt to stop
        }
        
    }
    
    
    lastObj = obj_value;
    lastInf_pr = inf_pr;
    
    
    return true;
}




bool MRQ_IpoptIntermediateCallback2:: intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value, Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    Ipopt::TNLPAdapter* tnlp_adapter = NULL;
    
    
    if( obj_value <= nI*integer_tol )
    {
        if( ip_cq != NULL )
        {
            Ipopt::OrigIpoptNLP* orignlp = NULL;
            orignlp = dynamic_cast<Ipopt::OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
            
            
            #if 0
            //just put this portion of code if you want stop in the solutions from restoration phase. 
            if(orignlp == NULL)
            {
                //If ipopt is in the restoration phase, orignlp will be null ( http://list.coin-or.org/pipermail/ipopt/2009-June/001581.html ). So, we must get a RestoIpoptNLP object
                
                #if MRQ_DEBUG_MODE
                    assert(mode == Ipopt::RestorationPhaseMode);
                #endif
                
                std::cerr << "orignlp e nulo! ipopt esta na fase de restauracao\n";
                
                Ipopt::RestoIpoptNLP* restonlp = dynamic_cast<Ipopt::RestoIpoptNLP*> (GetRawPtr(ip_cq->GetIpoptNLP()));
                
                
                orignlp = dynamic_cast<Ipopt::OrigIpoptNLP*> (&restonlp->OrigIpNLP());
            }
            #endif
            
            
            
            if( orignlp != NULL )
            {
                tnlp_adapter = dynamic_cast<Ipopt::TNLPAdapter*> (GetRawPtr(orignlp->nlp()));
                
                //If ipopt is in the restoration phase, tnlp_adapter will be NULL
            }
            
            
            
            if(tnlp_adapter)
            {
                tnlp_adapter->ResortX(*ip_data->curr()->x(), sol);
                
                const bool intSol = MRQ_isIntegerSol(nI, intVars, sol, integer_tol);
                
                if(intSol)
                {
                    //we stop ipopt because we have an integer solution, even if the solution is infeasible. We hope local search can find a good solution ;)
                    
                    //std::cout << "Interrompendo Ipopt por solucao inteira. integer_tol: " << integer_tol << "\n";
                    return false;
                }
            }
            
            
        }
        
        
        
    }
    
    return true;
}



MRQ_IpoptIntermediateCallback3:: MRQ_IpoptIntermediateCallback3() :MRQ_IpoptIntermediateCallback()
{
    nlpSolved = false;
    zu = INFINITY;
}


bool MRQ_IpoptIntermediateCallback3:: intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value, Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size, Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq)
{
    Ipopt::TNLPAdapter* tnlp_adapter = NULL;
    
    
    if( obj_value <= nI*integer_tol )
    {
        if( ip_cq != NULL )
        {
            Ipopt::OrigIpoptNLP* orignlp = NULL;
            orignlp = dynamic_cast<Ipopt::OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
            
            
            #if 0
            //just put this portion of code if you want stop in the solutions from restoration phase. 
            if(orignlp == NULL)
            {
                //If ipopt is in the restoration phase, orignlp will be null ( http://list.coin-or.org/pipermail/ipopt/2009-June/001581.html ). So, we must get a RestoIpoptNLP object
                
                #if MRQ_DEBUG_MODE
                    assert(mode == Ipopt::RestorationPhaseMode);
                #endif
                
                std::cerr << "orignlp e nulo! ipopt esta na fase de restauracao\n";
                
                Ipopt::RestoIpoptNLP* restonlp = dynamic_cast<Ipopt::RestoIpoptNLP*> (GetRawPtr(ip_cq->GetIpoptNLP()));
                
                
                orignlp = dynamic_cast<Ipopt::OrigIpoptNLP*> (&restonlp->OrigIpNLP());
            }
            #endif
            
            
            
            if( orignlp != NULL )
            {
                tnlp_adapter = dynamic_cast<Ipopt::TNLPAdapter*> (GetRawPtr(orignlp->nlp()));
                
                //If ipopt is in the restoration phase, tnlp_adapter will be NULL
            }
            
            
            
            if(tnlp_adapter)
            {
                
                
                tnlp_adapter->ResortX(*ip_data->curr()->x(), sol);
                
                const bool intSol = MRQ_isIntegerSol(nI, intVars, sol, integer_tol);
                
                if(intSol)
                {
                    double origObj = MRQ_INFINITY;
                    
                    if( !nlpSolved || zu < MRQ_INFINITY )
                    {
                        int r;
                        r = prob->objEval(thnumber, !prob->hasNlConstrs, sol, origObj);
                        
                        if(r != 0)
                        {
                            #if MRQ_DEBUG_MODE
                                MRQ_PRINTCALLBACKERRORNUMBER(r);
                            #endif
                            return true;
                        }
                        
                    }
                    
                    
                    /*if( origObj <= zu )
                    {
                        //if we alredy solved the local search problem, we only solve again if solution is feasible and objective is better than zu. Note, if zu is MRQ_INFINITY, objetive was not calculated, but there is ni problem about that
                        
                        bool feas = false;
                        double *cValues = nlp->constr;
                        
                        prob->isFeasibleToConstraints(thnumber, sol, true, NULL, absFeasTol, relFeasTol, feas, cValues );
                        
                        if( feas )
                        {
                            //std::cout << "vou resolver busca local pq encontrei solucao melhor! zu: " << zu << " obj: " << origObj << "\n";
                            
                            //for(int i = 0; i < nI; i++)
                                //std::cout << "sol["<<intVars[i]<<"]: " << sol[intVars[i]] << "\n";
                            
                            return false;
                        }
                    }*/
                    
                    
                    if(!nlpSolved || origObj <= zu)
                    {
                        
                        //std::cout << "Interrompendo Ipopt por solucao inteira. integer_tol: " << integer_tol << "\n";
                        
                        
                        MRQ_fixIntVarsOnSolByList( nI, intVars, sol, *nlp );
                        
                        nlp->setInitialSolution( sol, NULL, NULL);
                        
                        nlp->solve();
                        
                        if( nlp->feasSol && nlp->objValue <= zu )
                        {
                            //std::cout << "Resolvi busca local e consegui melhorar solucao. zu: " << zu << " novo obj: " << nlp->objValue << "\n";
                            return false;
                        }
                        else
                        {
                            //std::cout << "Resolvi busca local, mas nao adiantou. feas: " << nlp->feasSol << " zu: " << zu << " novo obj: " << nlp->objValue << "\n";
                            nlpSolved = true;
                        }
                        
                    }
                }
            }
            
            
        }
        
        
        
    }
    
    return true;
}




#endif




MRQ_IGMA1::MRQ_IGMA1():MRQ_Algorithm()
{
    resetParameters();
    resetOutput();
    out_algorithm = MRQ_IGMA1_ALG;
}



MRQ_IGMA1::~MRQ_IGMA1()
{
}



int MRQ_IGMA1::checkAlgorithmRequirements( MRQ_MINLPProb &prob, const double *lx, const double *ux)
{
    return MRQ_isBinarieProblemAtRegion(prob, lx, ux) ? 0 : MRQ_ALG_NOT_APPLICABLE;
}


void MRQ_IGMA1::printParameters(std::ostream &out) const
{
    char strValue[100];
    
    MRQ_Algorithm::printParameters(out);
    out << "\n"
    
    MRQ_STRFFATT(in_adopt_obj_cut) << "\n"
    MRQ_STRFFATT(in_consider_relax_infeas_if_solver_fail) << "\n"
    MRQ_STRFFATT(in_enable_gap_min_solver_premature_stoping) << "\n"
    MRQ_STRFFATT(in_heuristic_mode) << "\n"
    MRQ_STRFFATT(in_reorganize_lists) << "\n"
    MRQ_STRFFATT(in_set_special_gap_min_solver_params) << "\n"
    MRQ_STRFFATT(in_set_gap_exp_on_constr) << "\n"
    MRQ_STRFFATT(in_set_gap_ubound_constr) << "\n"
    MRQ_STRFFATT(in_use_random_initial_sols) << "\n"
    MRQ_STRFFATT(in_use_general_feas_heuristics) << "\n"
    MRQ_STRFFATT(in_lists_reorganization_frequency) << "\n"
    MRQ_STRFFATT(in_max_improvments_of_best_sol) << "\n"
    MRQ_STRFFATT(in_max_nonimprovment_integer_sols) << "\n"
    MRQ_STRFFATT(in_min_number_of_iters_on_gap_min_premature_stop) << "\n"
    MRQ_STRFFATT(in_seed_to_random_numbers) << "\n";
    
    MRQ_enumToStr(in_constr_branching_strategy, strValue);
    out << MRQ_STR(in_constr_branching_strategy) ": " << strValue << '\n';
    
    MRQ_enumToStr(in_exp_strategy, strValue);
    out << MRQ_STR(in_exp_strategy) ": " << strValue << '\n';
    
    MRQ_enumToStr(in_gap_min_solver, strValue);
    out << MRQ_STR(in_gap_min_solver) ": " << strValue << '\n';
    
    MRQ_enumToStr(in_gap_min_obj_strategy, strValue);
    out << MRQ_STR(in_gap_min_obj_strategy) ": " << strValue << '\n';
    
    
    out << MRQ_STRFFATT(in_absolute_obj_cut_eps) << "\n"
    MRQ_STRFFATT(in_relative_obj_cut_eps) << "\n"
    
    MRQ_STRFFATT(in_factor_to_increase_abs_obj_cut_eps) << "\n"
    MRQ_STRFFATT(in_factor_to_increase_rel_obj_cut_eps) << "\n"
    MRQ_STRFFATT(in_factor_to_increase_abs_obj_cut_eps_on_zero) << "\n"
    
    MRQ_STRFFATT(in_absolute_obj_tol_to_active_obj_cut) << "\n"
    MRQ_STRFFATT(in_relative_obj_tol_to_active_obj_cut) << "\n"
    MRQ_STRFFATT(in_factor_to_min_gap_obj_by_avg_gap) << "\n"
    
    //MRQ_STRFFATT(in_integer_tol_on_gap_min_premature_stop) << "\n"
    MRQ_STRFFATT(in_lower_bound_to_random_sol) << "\n"
    MRQ_STRFFATT(in_upper_bound_to_random_sol) << "\n"
    MRQ_STRFFATT(in_min_improv_to_obj_on_gap_min_premature_stop) << "\n"
    MRQ_STRFFATT(in_min_improv_to_infeas_on_gap_min_premature_stop) << "\n";
    
}


void MRQ_IGMA1::resetOutput()
{
    MRQ_Algorithm::resetOutput();
    
    out_number_of_open_nodes = 0;
    out_number_of_iterations_to_best_sol = ULONG_MAX;
    out_number_of_inner_iterations = 0;
    out_number_of_prunes_by_obj_cut_active = 0;
}


void MRQ_IGMA1::resetParameters()
{
    MRQ_Algorithm::resetParameters();
    
    in_gap_min_solver = MRQ_IPOPT; //we cannot use mosek to solve the nonconvex problems...
    
    in_adopt_obj_cut = true;
    in_consider_relax_infeas_if_solver_fail = true;
    in_heuristic_mode = false;
    in_reorganize_lists = true;
    in_enable_gap_min_solver_premature_stoping = false;
    in_use_random_initial_sols = true;
    in_use_general_feas_heuristics = true;
    
    in_set_gap_exp_on_constr = false;
    in_set_gap_ubound_constr = false;
    in_set_special_nlp_solver_params = false; //to avoid numerical problems. Note this is only to problems where integer variables are fixed...
    in_set_special_gap_min_solver_params = true;
    //in_stop_on_first_sol = false;
    
    in_lists_reorganization_frequency = 10000;
    in_max_improvments_of_best_sol = UINT_MAX;
    in_max_nonimprovment_integer_sols = 10;
    in_min_number_of_iters_on_gap_min_premature_stop = 200;
    in_seed_to_random_numbers = 1986;
    
    in_constr_branching_strategy = MRQ_BB_CBS_NO_CONSTRAINT_BRANCH;
    in_exp_strategy = MRQ_BB_ES_BEST_LIMIT;
    in_gap_min_obj_strategy = MRQ_IGMA_GMOS_BY_GAP_AVERAGE; //MRQ_IGMA2_GMOS_SAME_WEIGHT;
    
    in_absolute_obj_cut_eps = 1.0e-4;
    in_relative_obj_cut_eps = 1.0e-3;
    
    in_factor_to_increase_abs_obj_cut_eps = 1.5;
    in_factor_to_increase_rel_obj_cut_eps = 1.5;
    in_factor_to_increase_abs_obj_cut_eps_on_zero = 100;
    
    in_absolute_obj_tol_to_active_obj_cut = 1.0e-3;
    in_relative_obj_tol_to_active_obj_cut = 1.0e-3;
    
    in_factor_to_min_gap_obj_by_avg_gap = 10.0;
    
    in_min_improv_to_infeas_on_gap_min_premature_stop = 1e-3;
    in_min_improv_to_obj_on_gap_min_premature_stop = 1e-1;
    //in_integer_tol_on_gap_min_premature_stop = 1e-2;
    
    in_lower_bound_to_random_sol = -100;
    in_upper_bound_to_random_sol = 100;
    
    in_printing_frequency = 1000;
    //in_print_level = 6;
}



int MRQ_IGMA1::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_Algorithm::setIntegerParameter( name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_adopt_obj_cut), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_consider_relax_infeas_if_solver_fail), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_enable_gap_min_solver_premature_stoping), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_heuristic_mode), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_reorganize_lists), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_special_gap_min_solver_params), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_gap_exp_on_constr), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_gap_ubound_constr), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_random_initial_sols), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_use_general_feas_heuristics), name, value ) == 0 );
    else if( MRQ_setAtt<int>( MRQ_STRATT(in_lists_reorganization_frequency), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_max_improvments_of_best_sol), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_max_nonimprovment_integer_sols), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_min_number_of_iters_on_gap_min_premature_stop), name, value ) == 0 );
    else if( MRQ_setAtt<decltype(in_seed_to_random_numbers)>( MRQ_STRATT(in_seed_to_random_numbers), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}



int MRQ_IGMA1::setDoubleParameter(const char *name, const double value)
{
    int ret = MRQ_Algorithm::setDoubleParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    
    if( MRQ_setAtt( MRQ_STRATT(in_absolute_obj_cut_eps), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_relative_obj_cut_eps), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_factor_to_increase_abs_obj_cut_eps), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_factor_to_increase_rel_obj_cut_eps), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_factor_to_increase_abs_obj_cut_eps_on_zero), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_absolute_obj_tol_to_active_obj_cut), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_relative_obj_tol_to_active_obj_cut), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_factor_to_min_gap_obj_by_avg_gap), name, value ) == 0 );
    //else if( MRQ_setAtt( MRQ_STRATT(in_integer_tol_on_gap_min_premature_stop), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_lower_bound_to_random_sol), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_upper_bound_to_random_sol), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_min_improv_to_obj_on_gap_min_premature_stop), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_min_improv_to_infeas_on_gap_min_premature_stop), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}



int MRQ_IGMA1::setStringParameter(const char *name, const char *value)
{
    int r, ret = MRQ_Algorithm::setStringParameter( name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    
    if( (r = MRQ_setStrAtt( MRQ_STRATT(in_constr_branching_strategy), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_exp_strategy), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_gap_min_solver), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_gap_min_obj_strategy), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else
        ret = MRQ_NAME_ERROR;
    
    
    
    return ret;
}





int MRQ_IGMA1::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* gapMinSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    const int ndual = prob.m + prob.n + prob.n;
    
    bool updtSomeConstrBond;
    int nthreads;
    
    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux; 
    double *plc = NULL, *puc;
    
    
    MRQ_IGMA1BBCallbacks bbcallbacks(&prob, this, gapMinSolverParams, nlpSolverParams);
    BBL_BranchAndBound bb;
    
    
    MRQ_Preprocessor *preprocessor = NULL;
    
    clock_t clockStart, clockPrebb;
    double timeStart, timePrebb;
    
    
    
    timeStart = MRQ_getTime();
    clockStart = clock();
    
    
    nthreads = in_number_of_threads <= 0 ? BBL_getNumCores() : in_number_of_threads;
    
    
    
    
    if( in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_obj_function )
    {
        preprocessor = new (std::nothrow) MRQ_Preprocessor(&prob);
        
        if( !preprocessor )
        {
            if( in_print_level > 0 )
                MRQ_PRINTMEMERROR;
            out_return_code = MRQ_MEMORY_ERROR;
            goto termination;
        }
    }
    
    
    {
        auto ret = algorithmInitialization( nthreads, (in_preprocess_lin_constr || in_preprocess_obj_function || in_preprocess_quad_constrs), gapMinSolverParams, nlpSolverParams, prob, lx, ux, preprocessor, &updtSomeConstrBond, &plc, &puc );
        
        if( ret != 0 )
        {
            if( in_print_level > 0 )
            {
                if( ret == MRQ_ALG_NOT_APPLICABLE )
                    MRQ_PRINTERRORMSG("Error: Integrality Gap Minimization Algorithm only hands binary problems");
                else
                    MRQ_PRINTERRORNUMBER(ret);
            }
            out_return_code = ret;
            goto termination;
        }
    }
    
    MRQ_secDelete(preprocessor);
    
    
    if( updtSomeConstrBond )
        bbcallbacks.oplc = plc;
    else
    {
        free(plc);
        plc = NULL;
    }
    
    
    if( in_print_level > 1 )
        std::cout << "\n" MRQ_PREPRINT "Starting Integrality Gap Minimization Algorithm - version 1 (IGMA1)\n\n";
    
    if( in_print_level > 3 )
        printSubSolvers(false, true, false);
    
    
    if( !in_heuristic_mode && in_use_general_feas_heuristics )
    {
        int r;
        double obj;
        const double insideSolverMaxTime = 30;
        MRQ_Algorithm *alg;
        MRQ_HeuristicExecutor heurExec;
        
        
        heurExec.in_stop_on_first_sol = true;
        heurExec.in_use_igma1 = false;
        
        heurExec.setNumberOfThreads( in_number_of_threads );
        
        heurExec.setSolvers( in_milp_solver, in_nlp_solver );
        
        r = heurExec.insideRun( prob, NULL, nlpSolverParams, zu, obj, NULL, false, &alg, thnumber, 1, insideSolverMaxTime, lx, ux );
        
        if( r == MRQ_HEURISTIC_SUCCESS || r == MRQ_OPTIMAL_SOLUTION )
        {
            const bool flag = tryUpdateBestSolution(thnumber, prob.n, alg->out_best_sol, obj, 0, clockStart, timeStart, in_store_history_solutions );
            
            #if MRQ_DEBUG_MODE
                assert(flag);
            #endif
        }
    }
    
    bbcallbacks.clockStart = clockStart;
    bbcallbacks.timeStart = timeStart;
    
    
    bb.in_call_before_solving_relax_callback = in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_obj_function;
    bb.in_call_new_best_solution_callback = false;
    bb.in_call_updating_best_solution_callback = false;
    bb.in_call_end_of_iteration_callback = false;
    
    //bb.in_consider_relax_infeas_if_solver_fail = in_consider_relax_infeas_if_solver_fail;
    
    bb.in_consider_relax_infeas_if_solver_fail = false;
    
    bb.in_number_of_threads = nthreads;
    bb.in_max_iterations = in_max_iterations;
    bb.in_absolute_convergence_tol = -INFINITY; //we cannot stop by limits because we deturp lower bound concept
    bb.in_relative_convergence_tol = -INFINITY;
    bb.in_infinity_value = MRQ_INFINITY;
    bb.in_use_dual_obj_to_bound_prunning = true;
    bb.in_prune_nodes_by_bound = false;
    bb.in_lower_bound = zl;
    bb.in_upper_bound = zu;
    
    bb.in_store_parent_dual_solution_on_nodes = false;
    bb.in_store_parent_primal_solution_on_nodes = false;
    
    bb.in_print_level = in_print_level;
    bb.in_store_history_solutions = in_store_history_solutions;
    
    bb.in_reorganize_lists = in_reorganize_lists;
    bb.in_lists_reorganization_frequency = in_lists_reorganization_frequency;
    
    bb.in_printing_frequency = in_printing_frequency;
    
    
    bb.in_exp_strategy = MRQ_MRQ_exp_strategy2BBL_exp_strategy(in_exp_strategy);
    
    
    
    
    
    bb.in_max_cpu_time = in_max_cpu_time - (clock() - clockStart)/(double) CLOCKS_PER_SEC;
    bb.in_max_time = in_max_time - (MRQ_getTime() - timeStart);
    
    
    bbcallbacks.origBBLBranchStrategy = bb.in_branching_strategy;
    
    
    clockPrebb = clock();
    timePrebb = MRQ_getTime();
    
    
    bb.run(prob.n, prob.m, ndual, lx, ux, bbcallbacks);
    
    
    switch( bb.out_return_code )
    {
        case BBL_OPTIMAL_SOLUTION:
            out_return_code = in_heuristic_mode ? MRQ_HEURISTIC_SUCCESS : MRQ_OPTIMAL_SOLUTION;
            break;
            
        case BBL_FEASIBLE_SOLUTION:
            out_return_code = in_heuristic_mode ? MRQ_HEURISTIC_SUCCESS : MRQ_UNDEFINED_ERROR;
            break;
        
        case BBL_INFEASIBLE_PROBLEM:
            out_return_code = MRQ_INFEASIBLE_PROBLEM;
            break;
        
        case BBL_MAX_TIME_STOP:
            out_return_code = MRQ_MAX_TIME_STOP;
            break;
        
        case BBL_MAX_ITERATIONS_STOP:
            out_return_code = MRQ_MAX_ITERATIONS_STOP;
            break;
        
        case BBL_STOP_REQUIRED_BY_USER:
            out_return_code = MRQ_intToReturnCode( bb.out_return_subcode );
            break;
        
        default:
            out_return_code = MRQ_UNDEFINED_ERROR;
            break;
    }
    
    out_number_of_open_nodes = bb.out_number_of_open_nodes;
    
    out_number_of_iterations = bb.out_number_of_iterations;
    if( bb.out_feasible_sol )
    {
        if( bb.out_best_obj < zu )
        {
            out_best_obj = zu = bb.out_best_obj;
            MRQ_copyArray(prob.n, bb.out_best_sol, out_best_sol);
            
            if( zu <= -MRQ_INFINITY )
                out_return_code = MRQ_UNBOUNDED_PROBLEM;
        }
    }
    
    out_number_of_iterations_to_first_feas_sol = bb.out_first_sol_iter;
    out_number_of_iterations_to_best_sol = bb.out_best_sol_iter;
    out_cpu_time_to_first_feas_sol = bb.out_cpu_time_to_fisrt_sol;
    out_clock_time_to_first_feas_sol = bb.out_clock_time_to_fisrt_sol;
    out_cpu_time_to_best_sol = bb.out_cpu_time_to_best_sol;
    out_clock_time_to_best_sol = bb.out_clock_time_to_best_sol;
    
    out_number_of_feas_sols = bb.out_number_of_feas_sols;
    
    
    for(int i = 0; i <nthreads; i++)
        out_number_of_prunes_by_obj_cut_active += bbcallbacks.prunesByObjCutActive[i];
    
    for(int i = 0; i < nthreads; i++)
        out_number_of_inner_iterations = bbcallbacks.nsubiters[i];
    
    if( in_store_history_solutions )
    {
        BBL_SolutionHistory &shist = bb.out_sol_hist;
        const int nsols = shist.getnsols();
        
        const double pretime = timePrebb -timeStart;
        const double precputime = (clockPrebb - clockStart)/(double) CLOCKS_PER_SEC;
        
        for(int i = 0; i < nsols; i++)
        {
            BBL_HistorySolution *hsol = shist.getHistSolPointer(i);
            
            out_sol_hist.addSolution( prob.n, hsol->iter, hsol->time + pretime, hsol->cputime + precputime, hsol->sol, hsol->objF );
        }
        
    }
    
    
    
termination:
    
    if(preprocessor)	delete preprocessor;
    
    
    if(plc)		free(plc);
    
    
    out_feasible_solution = out_best_obj < MRQ_INFINITY;
    out_lower_bound = zl;
    out_upper_bound = zu;
    out_number_of_threads = nthreads;
    out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
    out_clock_time = MRQ_getTime() - timeStart;
    
    
    algorithmFinalization(nthreads, prob, lx, ux);
    
    
    
    return out_return_code;
}





