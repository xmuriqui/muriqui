
#include <ctime>
#include <cmath>
#include <climits>

#include <iostream>
#include <new>


#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_advanced.hpp"


using namespace branchAndBound;
using namespace optsolvers;
using namespace muriqui;




static inline double MRQ_quadEuclideanDistance(const unsigned int ninds, const int *inds, const double *sol1, const double *sol2)
{
    double d = 0.0;
    
    for(int unsigned i = 0; i < ninds; i++)
    {
        const int ind = inds[i];
        const double f = sol1[ind] - sol2[ind];
        
        d += f*f;
    }
    
    return d;
}




class MRQ_IGMA2UserCallbacks : public MRQ_UserCallbacks
{
    MRQ_MINLPProb *prob;
    MRQ_IGMA2 *igma2;
    MRQ_GeneralSolverParams *gapMinSolverParams;
    MRQ_GeneralSolverParams *nlpSolverParams;
    
    int nI, nthreads;
    int nC; //number of continuous vars
    int *intVars, *contVars;
    
    MRQ_GapMinProb *gapmins;
    //MRQ_NLPSolver2 **nlps;
    double **auxVars;
    
    MRQ_IGMA2Iteration igma2Iter;
    
    
    int allocateThreadStructures(const unsigned int nthreads);
    
    void desallocateThreadStructures();
    
    void desallocate();
    
public:
    
    double obj_opt_at_continuous_relax;
    
    MRQ_IGMA2UserCallbacks(MRQ_MINLPProb *prob, MRQ_IGMA2 *alg, MRQ_GeneralSolverParams *gapMinSolverParams, MRQ_GeneralSolverParams *nlpSolverParams) : MRQ_UserCallbacks()
    {
        this->prob = prob;
        this->igma2 = alg;
        this->gapMinSolverParams = gapMinSolverParams;
        this->nlpSolverParams = nlpSolverParams;
        
        nthreads = 0;
        nI = 0;
        nC = 0;
        intVars = NULL;
        contVars = NULL;
        
        gapmins = NULL;
        //nlps = NULL;
        auxVars = NULL;
        
        obj_opt_at_continuous_relax = -INFINITY;
    }
    
    
    
    virtual int beforeAll(const MRQ_ALG_CODE algCode, const unsigned int numberOfThreads) override;
    
    
    virtual int BB_afterSolvingRelaxation(const unsigned int threadNumber, MRQ_NewBBNode &node, const long unsigned int iter, const double lb, const double ub, const double *nlx, const double *nux, MRQ_NLPSolver &nlpSolver, const int status, double &objFSolution, double &dualObjFSolution, double *sol, double *constrs, double *dualSolC, double *dualSolV, bool &pruneNode, bool &solveAgain, bool &branchEvenIntegerSol) override;
    
};






MRQ_IGMA2::MRQ_IGMA2() : MRQ_BranchAndBound()
{
    resetParameters();
    resetOutput();
    
    out_algorithm = MRQ_IGMA2_ALG;
}


MRQ_IGMA2::~MRQ_IGMA2()
{
}


void MRQ_IGMA2::printParameters(std::ostream &out) const
{
    char strValue[100];
    
    MRQ_BranchAndBound::printParameters(out);
    out << "\n"
    
    MRQ_STRFFATT(in_set_max_dist_constr_on_bin_vars) << "\n"
    MRQ_STRFFATT(in_set_special_gap_min_solver_params) << "\n"
    MRQ_STRFFATT(in_solve_local_search_problem_even_on_non_int_sol) << "\n";
    
    MRQ_enumToStr(in_neighborhood_strategy, strValue);
    out << MRQ_STR(in_neighborhood_strategy) ": " << strValue << "\n"
    
    MRQ_STRFFATT(in_factor_to_max_dist_constr) << "\n"
    MRQ_STRFFATT(in_factor_to_max_dist_constr_on_bin_vars) << "\n"
    MRQ_STRFFATT(in_percentual_to_rectangular_neighborhood) << "\n";
    
}


void MRQ_IGMA2::resetOutput()
{
    MRQ_BranchAndBound::resetOutput();
    
    out_feas_sol_on_gap_min_problem = false;
    out_feas_sol_found_by_igma2_procedure = false;
}


void MRQ_IGMA2::resetParameters()
{
    MRQ_BranchAndBound::resetParameters();
    
    /*TODO: 
    in_number_of_node_sublists = 1;
    in_lists_reorganization_frequency = INT_MAX; */
    
    in_int_feas_heurs_strategy = MRQ_BB_IHS_NO_HEURISTICS; //in_use_general_int_feas_heuristics =false;
    in_use_outer_app = false;
    in_use_outer_app_as_heuristic = false;
    in_rounding_heuristic_strategy = MRQ_RS_NO_ROUNDING; //in_use_round_heuristic = false;
    
    in_branching_strategy = MRQ_BB_BS_BIN_FIRST_HIGHEST_INT_GAP; //TODO: MRQ_BB_BS_HIGHEST_INT_GAP
    in_constr_branching_strategy = MRQ_BB_CBS_NO_CONSTRAINT_BRANCH;
    in_exp_strategy = MRQ_BB_ES_DEPTH;
    
    
    
    in_set_max_dist_constr_on_bin_vars = false;
    in_set_special_gap_min_solver_params = true;
    in_solve_local_search_problem_even_on_non_int_sol = true;
    in_nlp_solver = MRQ_getDefaultMinGapNLPSolverCode();
    in_neighborhood_strategy = MRQ_IGMA2_NS_SPHERIC;
    in_factor_to_max_dist_constr = 0.05;
    in_factor_to_max_dist_constr_on_bin_vars = 0.05;
    in_percentual_to_rectangular_neighborhood = 0.05;
}


int MRQ_IGMA2:: checkAlgorithmRequirements( MRQ_MINLPProb &prob, const double *lx, const double *ux)
{
    return MRQ_isBinarieProblemAtRegion(prob, lx, ux) ? 0 : MRQ_ALG_NOT_APPLICABLE;
}


int MRQ_IGMA2::setDoubleParameter(const char *name, const double value)
{
    int ret = MRQ_BranchAndBound::setDoubleParameter( name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    
    if( MRQ_setAtt<double>( MRQ_STRATT(in_factor_to_max_dist_constr), name, value ) == 0 );
    else if( MRQ_setAtt<double>( MRQ_STRATT(in_factor_to_max_dist_constr_on_bin_vars), name, value ) == 0 );
    else if( MRQ_setAtt<double>( MRQ_STRATT(in_percentual_to_rectangular_neighborhood), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


int MRQ_IGMA2::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_BranchAndBound::setIntegerParameter( name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_max_dist_constr_on_bin_vars), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_set_special_gap_min_solver_params), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_solve_local_search_problem_even_on_non_int_sol), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


int MRQ_IGMA2::setStringParameter(const char *name, const char *value)
{
    int r, ret = MRQ_BranchAndBound::setStringParameter( name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    
    
    if( (r = MRQ_setStrAtt( MRQ_STRATT(in_neighborhood_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else
        ret = MRQ_NAME_ERROR;
    
    
    
    return ret;
}


int MRQ_IGMA2::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* gapMinSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    MRQ_IGMA2UserCallbacks bbcallbacks(&prob, this, gapMinSolverParams, nlpSolverParams); //remember: those callbacks are muriqui callbacks, not bbl callbakcs;
    
    
    if( in_print_level > 1 )
        std::cout << "\n" MRQ_PREPRINT "Starting Integrality Gap Minimization Algorithm - version 2 (IGMA2)\n\n";
    
    
    in_call_after_solving_relax_callback = true;
    MRQ_BranchAndBound::in_user_callbacks = &bbcallbacks;
    
    in_igma2_strategy = muriqui::MRQ_BB_IHS_NO_HEURISTICS; //strategy for igma2 inside B&B
    
    
    //we just pass gapMinSolverParams to be printed if user wants...
    /*if(run_by_inside)
    {
        MRQ_BranchAndBound::insideRun(prob, gapMinSolverParams, nlpSolverParams, thnumber, insideNumberOfThreads, insideSolverMaxTime, nlx, nux);
    }
    else */
    {
        MRQ_BranchAndBound::run(prob, gapMinSolverParams, nlpSolverParams);
    }
    
    //std::cout << "out_return_code: " << out_return_code << "\n";
    
    if( out_return_code == MRQ_STOP_REQUIRED_BY_USER && out_best_obj < MRQ_INFINITY )
    {
        out_return_code = MRQ_HEURISTIC_SUCCESS;
        
        if(in_print_level > 1)
            std::cout << "IGMA2 found a feasible solution! Obj Function: " << out_best_obj << " \n";
        
    }
    
    MRQ_BranchAndBound::in_user_callbacks = NULL;
    
    
    if( !(bbcallbacks.obj_opt_at_continuous_relax >= out_obj_opt_at_continuous_relax) ) // we have to take care of nan
        out_obj_opt_at_continuous_relax = bbcallbacks.obj_opt_at_continuous_relax;
    
    if( out_lower_bound <= -MRQ_INFINITY && out_obj_opt_at_continuous_relax > -MRQ_INFINITY )
        out_lower_bound = out_obj_opt_at_continuous_relax;
    
    
    return out_return_code;
}





int  MRQ_IGMA2UserCallbacks::allocateThreadStructures (const unsigned int nthreads)
{
    const int n = prob->n;
    const double *lx = prob->lx;
    const double *ux = prob->ux;
    
    this->nthreads = nthreads;
    
    
    gapmins = new (std::nothrow) MRQ_GapMinProb[nthreads];
    
    //MRQ_calloc(nlps, nthreads);
    MRQ_calloc(auxVars, nthreads);
    
    if(!gapmins || !auxVars)
    {
        if(igma2->in_print_level > 0)
            MRQ_PRINTMEMERROR;
        return MRQ_MEMORY_ERROR;
    }
    
    
    for(unsigned int i = 0; i < nthreads; i++)
    {
        //const int naddconstrs = (igma2->in_neighborhood_strategy == MRQ_IGMA2_NS_RECTANGULAR ? 0 : 1) + (int) igma2->in_set_max_dist_constr_on_bin_vars ;
        
        int r = gapmins[i].setProblem(igma2->in_nlp_solver, *prob, lx, ux, gapMinSolverParams, i, igma2->in_set_special_gap_min_solver_params, false, false, 1, igma2->in_max_cpu_time, igma2->in_max_time, 0,  1 + (int) igma2->in_set_max_dist_constr_on_bin_vars);
        
        if( r != 0 )
        {
            if(igma2->in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
            return r;
        }
        
        ((OPT_NLPSolver*) gapmins[i].solver)->in_absolute_feas_tol = igma2->in_absolute_feasibility_tol;
        ((OPT_NLPSolver*) gapmins[i].solver)->in_relative_feas_tol = igma2->in_relative_feasibility_tol;
        
        
        MRQ_malloc(auxVars[i], 2*n );
        if(!auxVars[i])
        {
            if(igma2->in_print_level > 0)
                MRQ_PRINTMEMERROR;
            return MRQ_MEMORY_ERROR;
        }
        
    }
    
    return 0;
}


int MRQ_IGMA2UserCallbacks::beforeAll(const MRQ_ALG_CODE algCode, const unsigned int numberOfThreads)
{
    const int n = prob->n;
    int r;
    
    nI = prob->getNumberOfIntegerVars();
    nC = n - nI;
    
    MRQ_malloc(intVars, nI);
    MRQ_malloc(contVars, nC);
    
    if(!intVars || !contVars)
    {
        if(igma2->in_print_level > 0)
            MRQ_PRINTMEMERROR;
        return MRQ_MEMORY_ERROR;
    }
    
    prob->getIntegerIndices(intVars);
    prob->getContinuousIndices(contVars);
    
    r = allocateThreadStructures(numberOfThreads);
    if(r != 0)
    {
        if(igma2->in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        return r;
    }
    
    
    igma2Iter.prob = prob;
    igma2Iter.nI = nI;
    igma2Iter.nC = nC;
    
    igma2Iter.intVars = intVars;
    igma2Iter.contVars = contVars;
    
    igma2Iter.in_set_max_dist_constr_on_bin_vars = igma2->in_set_max_dist_constr_on_bin_vars;
    igma2Iter.in_solve_local_search_problem_even_on_non_int_sol = igma2->in_solve_local_search_problem_even_on_non_int_sol ;
    igma2Iter.in_print_level = igma2->in_print_level;
    igma2Iter.in_neighborhood_strategy = igma2->in_neighborhood_strategy;
    igma2Iter.in_factor_to_max_dist_constr_on_bin_vars = igma2->in_factor_to_max_dist_constr_on_bin_vars;
    igma2Iter.in_percentual_to_rectangular_neighborhood = igma2->in_percentual_to_rectangular_neighborhood;
    igma2Iter.in_integer_tol = igma2->in_integer_tol;
    
    return 0;
}


void MRQ_IGMA2UserCallbacks::desallocate()
{
    MRQ_secFree(intVars);
    nI = 0;
    MRQ_secFree(contVars);
    nC = 0;
}


void MRQ_IGMA2UserCallbacks:: desallocateThreadStructures()
{
    MRQ_secDeleteArray(gapmins);
    
    if(auxVars)
    {
        for(int i = 0; i < nthreads; i++)
        {
            if(auxVars[i])
                free(auxVars[i]);
        }
        
        free(auxVars);
        auxVars = NULL;
    }
}


#if 1

int MRQ_IGMA2UserCallbacks::BB_afterSolvingRelaxation(const unsigned int threadNumber, MRQ_NewBBNode &node, const long unsigned int iter, const double lb, const double ub, const double *nlx, const double *nux, MRQ_NLPSolver &nlpSolver, const int status, double &objFSolution, double &dualObjFSolution, double *sol, double *constrs, double *dualSolC, double *dualSolV, bool &pruneNode, bool &solveAgain, bool &branchEvenIntegerSol)
{
    pruneNode = false;
    solveAgain = false;
    branchEvenIntegerSol = false;
    
    
    if( ub < igma2->in_upper_bound && ub < MRQ_INFINITY )
    { //some other procedure in B&B found a feasible solution. We can stop
        return MRQ_HEURISTIC_SUCCESS;
    }
    
    
    
    if(nlpSolver.feasSol)
    {
        //checking if current solution from continuous relaxation is integer
        if( MRQ_isIntegerSol(nI, intVars, sol, igma2->in_integer_tol) )
        {
            tryUpdateBestSolution(threadNumber, prob->n, sol, objFSolution, iter);
            
            
            //std::cout << "current continuous relaxation solution is integer!\n";
            
            return MRQ_HEURISTIC_SUCCESS;
        }
        
        
        const int n = prob->n;
        const int m = prob->m;
        
        MRQ_GapMinProb &gapminsolver = gapmins[threadNumber];
        MRQ_NLPSolver *gapmin = (MRQ_NLPSolver *) gapminsolver.solver;
        
        const double maxDistance = 1.0 + ceil(igma2->in_factor_to_max_dist_constr *n);
        
        const int distConstIndex = m; //that only works because we are not adopting the objective cut. If we adopt the objective cut, this index is not more this value...
        
        
        
        {
            double *auxVars = this->auxVars[threadNumber];
            double *auxVars2 = &auxVars[n];
            int r;
            double clb, cub;
            
            
            //preprocessor can have changed constraint bounds in this node. nlpSolver already have the updated bounds
            for(int i = 0; i < m; i++)
            {
                r = nlpSolver.getConstraintBounds(i, clb, cub);
                if(r == 0)
                {
                    r = gapmin->setConstraintBounds(i, clb, cub);
                    if(r != 0)
                    {
                        if(igma2->in_print_level > 0)
                            MRQ_PRINTERRORNUMBER(r);
                    }
                }
                else
                {
                    if(igma2->in_print_level > 0)
                        MRQ_PRINTERRORNUMBER(r);
                }
            }
            
            const bool nlpSolver_feasSol = nlpSolver.feasSol;
            const int nlpSolver_retCode = nlpSolver.retCode;
            const int nlpSolver_origSolverRetCode = nlpSolver.origSolverRetCode;
            const double nlpSolver_objValue = nlpSolver.objValue;
            const double nlpSolver_dualObjValue = nlpSolver.dualObjValue;
            double *nlpSolver_sol = auxVars;
            
            bool gapMinIntSol;
            double objOutSol;
            double *outSol = auxVars2;
            
            MRQ_Preprocessor *preprocessor = BB_getPreprocessorPointer(threadNumber);
            minlpproblem:: MIP_ConstraintsByColumnsStorager *ccstorager = BB_getConstraintsByColumnsStoragerPointer();
            
            
            MRQ_copyArray(n, nlpSolver.sol, nlpSolver_sol);
            
            
            r = igma2Iter.run(threadNumber, gapminsolver, &nlpSolver, nlx, nux, distConstIndex, maxDistance, nlpSolver_sol, dualSolC, dualSolV, outSol, objOutSol, preprocessor, ccstorager, &gapMinIntSol);
            
            nlpSolver.feasSol = nlpSolver_feasSol;
            nlpSolver.retCode = nlpSolver_retCode;
            nlpSolver.origSolverRetCode = nlpSolver_origSolverRetCode;
            nlpSolver.objValue = nlpSolver_objValue;
            nlpSolver.dualObjValue = nlpSolver_dualObjValue;
            
            MRQ_copyArray(n, nlpSolver_sol, nlpSolver.sol);
            
            //maybe it is not necessary, but we restore the node bounds
            /*r2 = nlpSolver.setnVariablesBounds(n, nlx, nux);
            if( r2 != 0)
            {
                if(igma2->in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(r2);
                //we report error, but we do not stop algorithm
            }
            */
            
            if(r == MRQ_HEURISTIC_SUCCESS)
            {
                bool updt = tryUpdateBestSolution(threadNumber, n, outSol, objOutSol, iter);
                    
                if(updt)
                {
                    //std::cout << "Encontramos solucao viavel\n";
                    
                    igma2->out_feas_sol_on_gap_min_problem = gapMinIntSol;
                    
                    igma2->out_feas_sol_found_by_igma2_procedure = true;
                    
                    
                    if(iter == 1 && nlpSolver.retCode == OPT_OPTIMAL_SOLUTION)
                        obj_opt_at_continuous_relax = nlpSolver.objValue;
                    
                    return MRQ_HEURISTIC_SUCCESS;
                }
                
            }
            else if(r != MRQ_HEURISTIC_FAIL)
            { //if we got an heuristic fail, we just continue the computation. For any other code, we stop the algorithm
                return r;
            }
            
        }
        
    }
    
    
    return 0;
}

#endif


#if 0

int MRQ_IGMA2UserCallbacks::BB_afterSolvingRelaxation(const int threadNumber, MRQ_NewBBNode &node, const long unsigned int iter, const double lb, const double ub, const double *nlx, const double *nux, MRQ_NLPSolver2 &nlpSolver, const int status, double &objFSolution, double &dualObjFSolution, double *sol, double *constrs, double *dualSolC, double *dualSolV, double &nodeLowerBound, bool &pruneNode, bool &solveAgain, bool &branchEvenIntegerSol)
{
    pruneNode = false;
    solveAgain = false;
    branchEvenIntegerSol = false;
    
    
    if( ub < igma3->in_upper_bound )
    { //some other procedure in B&B found a feasible solution. We can stop
        return MRQ_HEURISTIC_SUCCESS;
    }
    
    
    //checking if current solution from continuous relaxation is integer
    if( MRQ_isIntegerSol(nI, intVars, sol, igma3->in_integer_tol) )
    {
        tryUpadteBestSolution(threadNumber, prob->n, sol, objFSolution, iter);
        
        
        return MRQ_HEURISTIC_SUCCESS;
    }
    
    
    if(nlpSolver.feasSol)
    {
        const int n = prob->n;
        const int m = prob->m;
        
        MRQ_GapMinProb &gapminsolver = gapmins[threadNumber];
        MRQ_NLPSolver2 *nlp = nlps[threadNumber];
        MRQ_NLPSolver2 *gapmin = (MRQ_NLPSolver2 *) gapminsolver.solver;
        
        const double maxDistance = 1.0 + ceil(igma3->in_factor_to_max_dist_constr *n);
        
        const int distConstIndex = m; //that only works because we are not adopting the objective cut. If we adopt the objective cut, this index is not more value...
        
        
        
        //updating variables and constraint bounds...
        {
            int r;
            double clb, cub;
            
            r = gapmin->setnVariablesBounds(n, nlx, nux);
            if(r != 0)
            {
                if(igma3->in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(r);
            }
            
            
            //preprocessor can have changed constraint bounds in this node
            for(int i = 0; i < m; i++)
            {
                r = nlpSolver.getConstraintBounds(i, clb, cub);
                if(r == 0)
                {
                    r = gapmin->setConstraintBounds(i, clb, cub);
                    if(r != 0)
                    {
                        if(igma3->in_print_level > 0)
                            MRQ_PRINTERRORNUMBER(r);
                    }
                }
                else
                {
                    if(igma3->in_print_level > 0)
                        MRQ_PRINTERRORNUMBER(r);
                }
            }
        }
        
        
        
        #if 1
        {
            OPT_SolutionDistanceConstraintSetter sdcs;
            double *auxVars = nlp->sol; //using nlp.sol as aux array..
            
            
            int r = sdcs.setDistConstraint(gapmin, distConstIndex, nC, contVars, sol, maxDistance, auxVars);
            
            if(r != 0)
            {
                if(igma3->in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(r);
                
                return MRQ_NLP_SOLVER_ERROR;
            }
            
            if( igma3->in_set_max_dist_constr_on_bin_vars )
            {
                const int distIntConstIndex = distConstIndex + 1;
                const double tolMaxIntDistance = 0.001;
                double maxIntDistance = 0.0;
                
                //maxDistance for integer vars is the total integer gap plus a tolerance plus a slack to change some variables values
                
                for(int i = 0; i < nI; i++)
                    maxIntDistance += MRQ_gap( sol[intVars[i]] );
                
                maxIntDistance += 1.0 + ceil(igma3->in_factor_to_max_dist_constr_on_bin_vars *nI) + tolMaxIntDistance; //the tolerance is to due to possible numerical errors when we sum the values to calculate maxIntDistance
                
                
                int r = sdcs.setDistConstraint(gapmin, distIntConstIndex, nI, intVars, sol, maxIntDistance, auxVars);
                
                if(r != 0)
                {
                    if(igma3->in_print_level > 0)
                        MRQ_PRINTERRORNUMBER(r);
                    
                    return MRQ_NLP_SOLVER_ERROR;
                }
            }
            
            
            //if(gapmin->isMyNLPClass())
                //((OPT_MyNLPSolver *) gapmin)->prob.print(), MRQ_getchar();
        }
        #endif
        
        
        
        
        
        
        gapmin->setInitialSolution(sol, dualSolC, dualSolV);
        
        
        gapmin->solve(false);
        
        std::cout << "obj: " << gapmin->objValue << "\n";
        
        for(int i = 0; i < n; i++)
            std::cout << "sol["<<i<<"]: " << gapmin->sol[i] << " \t"  << std::endl;
        #if 0
        {
            double dlb, dub;
            unsigned int iters = 0;
            
            gapmin->getConstraintBounds(distConstIndex, dlb, dub);
            
            double qdist = MRQ_quadEuclideanDistance(nC, contVars, sol, gapmin->sol);
            
            gapmin->getNumberOfIterations(iters);
            
            std::cout << "\n\t\tigma2 iter: " << iter << " Min Gap solving - ret: " << gapmin->retCode << " feas sol: " << gapmin->feasSol << " obj: " << gapmin->objValue << " dist: " << qdist <<  " dist ub: " << maxDistance*maxDistance << " iters: " << iters << "\n";
        }
        #endif
        
        
        if(gapmin->feasSol)
        {
            const bool gapMinIntSol = MRQ_isIntegerSol(nI, intVars, gapmin->sol, igma3->in_integer_tol);
            
            if( gapMinIntSol || igma3->in_solve_local_search_problem_even_on_non_int_sol )
            {
                double obj;
                double *psol = NULL;
                
                
                nlp->setnVariablesBounds(n, nlx, nux);
                
                for(int i = 0; i < m; i++)
                {
                    double clb, cub;
                    
                    int r = nlpSolver.getConstraintBounds(i, clb, cub);
                    if(r == 0)
                    {
                        nlp->setConstraintBounds(i, clb, cub);
                        if(r != 0)
                        {
                            if(igma3->in_print_level > 0)
                                MRQ_PRINTERRORNUMBER(r);
                        }
                    }
                    else
                    {
                        if(igma3->in_print_level > 0)
                            MRQ_PRINTERRORNUMBER(r);
                    }
                    
                }
                
                
                MRQ_fixIntVarsOnSolByList(nI, intVars, gapmin->sol, *nlp);
                
                nlp->setInitialSolution( gapmin->sol, gapmin->dualSolC, gapmin->dualSolV );
                
                nlp->solve();
                
                std::cout << "obj: " << gapmin->objValue << "\n";
                for(int i = 0; i < n; i++)
                    std::cout << "local search sol["<<i<<"]: " << nlp->sol[i] << " \t"  << std::endl;
                
                //std::cout << "\t\t\tBusca local. retCode: " << nlp->retCode << " feasSol: " << nlp->feasSol << " obj: " << nlp->objValue << "\n";
                
                if( nlp->feasSol )
                {
                    obj = nlp->objValue;
                    psol = nlp->sol;
                }
                else if(gapMinIntSol)
                {
                    psol = gapmin->sol;
                    
                    int r = prob->objEval(threadNumber, true, psol, obj);
                    
                    if(r != 0)
                    {
                        psol = NULL;
                        if(igma3->in_print_level > 0)
                            MRQ_PRINTCALLBACKERRORNUMBER(r);
                    }
                    
                }
                
                if(psol)
                {
                    bool updt = tryUpadteBestSolution(threadNumber, n, psol, obj, iter);
                    
                    if(updt)
                    {
                        //std::cout << "Encontramos solucao viavel\n";
                        
                        igma3->out_feas_sol_on_gap_min_problem = gapMinIntSol;
                        
                        igma3->out_feas_sol_found_by_igma3_procedure = true;
                        
                        
                        if(iter == 1 && nlpSolver.retCode == OPT_OPTIMAL_SOLUTION)
                            obj_opt_at_continuous_relax = nlpSolver.objValue;
                        
                        return MRQ_HEURISTIC_SUCCESS;
                    }
                }
                
            }
            
        }
        
    }
    
    return 0;
}

#endif










