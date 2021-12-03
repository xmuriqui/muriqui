/*That file contains a implementation of
* Relaxation Enforced Neighborhood Search (RENS) heuristic. 
*
* References:
*
* Berthold, RENS The Optimal rounding. Math. prog. Comp 6 (2014), pages 33-54.
*
*
* Author: Wendel  Melo
*
* Date: 14-June-2017
*/


#include <cmath>
#include <ctime>
#include <cassert>

#include <iostream>
#include <new>


#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_advanced.hpp"
#include "MRQ_rens.hpp"


using namespace optsolvers;
using namespace muriqui;






MRQ_RensExecutor::MRQ_RensExecutor()
{	
    neighConstrIndex = -1;
    rensNeighStrategy = MRQ_SNS_ORIGINAL;
    mysubprob = NULL;
    subprob = NULL;
    auxInds = NULL;
    nI = -1;
    intVars = NULL;
    nC = -1;
    contVars = NULL;
    auxValues = NULL;
    nlx = NULL;
    newEval = NULL;
    
    out_alg = NULL;
    heurs = NULL;
    
    resetParameters();
    resetOutput();
}

MRQ_RensExecutor::~MRQ_RensExecutor()
{
    deallocate();
}

void MRQ_RensExecutor::deallocate()
{
    MRQ_secDelete(mysubprob);
    subprob = NULL;
    MRQ_secDelete(newEval);
    MRQ_secFree(auxInds);
    MRQ_secFree(intVars);
    MRQ_secFree(auxValues);
    MRQ_secFree(nlx);
    contVars = NULL; //constVars is allocated togheter intVars
    nI = -1;
    nC = -1; //contVars is allocated togheter intVars
    
    deallocateOutAlg();
}


void MRQ_RensExecutor::deallocateOutAlg()
{
    if(out_alg != in_alg)
    {
        if(heurs == NULL || !heurs->hasAlgPointer(out_alg) ) //if ouf_alg is an algorithm from heuristic executor, it wiil be desallocated by  heurs destructor
            MRQ_secDelete(out_alg);
    }
    
    if( heurs != in_heurs )
        MRQ_secDelete(heurs);
    
    out_alg = NULL;
    heurs = NULL;
}


void MRQ_RensExecutor::resetOutput()
{
    deallocateOutAlg();
}


void MRQ_RensExecutor::resetParameters()
{
    in_apply_only_heuristics_on_subproblems = false;
    in_solve_continuous_relaxation_on_subproblem = true;
    in_stop_on_first_sol = false;
    in_print_level = 5;
    in_nlp_solver = MRQ_getDefaultNLPSolverCode();
    in_milp_solver = MRQ_getDefaultMILPSolverCode();
    in_integer_tol = 1e-3;
    in_integer_neighborhood_factor = 0.1;
    in_continuous_neighborhood_factor = 0.1;
    in_max_cpu_time = INFINITY;
    in_max_time = INFINITY;
    in_algorithm_to_solve_subproblem = MRQ_UNDEFINED_ALG;
    
    in_alg = NULL;
    in_heurs = NULL;
}


int MRQ_RensExecutor::setSubProb(MRQ_MINLPProb &prob, int rensNeighStrategy)
{
    const int n = prob.n;
    this->rensNeighStrategy = rensNeighStrategy;
    
    deallocate();
    
    nI = prob.getNumberOfIntegerVars();
    
    
    MRQ_malloc(nlx, 2*n);
    MRQ_IFMEMERRORRETURN(!nlx);
    
    nux = &nlx[n];
    
    
    if(rensNeighStrategy == MRQ_SNS_ORIGINAL)
    {
        subprob = &prob;
    }
    else if( rensNeighStrategy == MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD || rensNeighStrategy == MRQ_SNS_EUCLIDEAN_INTEGER_NEIGHBORHOOD || rensNeighStrategy == MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD )
    { //we will set aditional constraints to local branch neighbor
        int r;
        const int nnewConstraints = (rensNeighStrategy == MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD && nI < n) ? 2 : 1; //to set the continuous euclidean neighborhood, we must have at least one continuous var
        
        
        mysubprob = new (std::nothrow) MRQ_MINLPProb;
        MRQ_IFMEMERRORRETURN(!mysubprob);
        
        
        subprob = mysubprob;
        
        r = subprob->copyProblemFrom(prob);
        MRQ_IFERRORRETURN(r, MRQ_UNDEFINED_ERROR);
        
        //we need set another object to perform nonlinear evals becuase we are adding one more constraint in the subproblem. So, the user calbacks would receive the incorrect number of constraints if we did not do it.
        newEval = new (std::nothrow) minlpproblem::MIP_EncapsulatedNonLinearEval(&prob, prob.getNonLinearEvaluationObject());
        
        MRQ_IFMEMERRORRETURN(!newEval);
        
        subprob->setNonLinearEvaluationObject(newEval);
        
        //adding the neighborhood constraint
        neighConstrIndex = subprob->getNumberOfConstraints();
        
        r = subprob->addConstraints(nnewConstraints);
        MRQ_IFERRORRETURN(r, MRQ_UNDEFINED_ERROR);
        
    }
    else
    {
        assert(false);
    }
    
    
    
    
    MRQ_malloc(auxInds, subprob->n);
    MRQ_malloc(auxValues, subprob->n);
    MRQ_malloc(intVars, subprob->n); //we use n instead of nI because we take advatnage this array to store continous indicies also
    MRQ_IFMEMERRORRETURN(!auxInds || !auxValues || !intVars);
    
    subprob->getIntegerIndices(intVars);
    
    if( rensNeighStrategy == MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD )
    {
        contVars = &intVars[nI];
        nC = subprob->getContinuousIndices(contVars);
        
        #if MRQ_DEBUG_MODE
            assert( nI + nC == subprob->n );
        #endif
    }
    
    
    return 0;
}


/*
 * if neighStrategy ==  MRQ_SNS_ORIGINAL, we set the variable bounds to perform the local search. Original reans calculate the optimal rouding seting newlx[i] = floor( sol[i] ) and newx[i] = ceil( sol[i] ) for integer variable i.
 */ 
int muriqui::MRQ_setSubproblemNeighborhood(const double *inputlx, const double *inputux, const double *sol, const int nI, const int *intVars, const int nC, const int *contVars, const int neighStrategy, const double in_integer_neighborhood_factor, const double in_continuous_neighborhood_factor, const double in_integer_tol, const int neighConstrIndex, int *auxInds, double *auxValues, MRQ_MINLPProb *prob, double *newlx, double *newux)
{
    const int n = prob->n;
    int nINonFixed = 0;
    
    for(int i = 0; i < nI; i++)
    {
        const int ind = intVars[i];
        if( inputlx[ind] != inputux[ind] )
            nINonFixed++;
    }
    
    #if MRQ_DEBUG_MODE
        assert(nINonFixed <= nI);
    #endif
    
    
    
    MRQ_copyArray(n, inputlx, newlx);
    MRQ_copyArray(n, inputux, newux);
    
    if( neighStrategy == MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD )
    {
        //double *auxValues = auxValues;
        
        int r, rhs = 0;
        int nInds = 0;
        
        //*** here, we are addressing the binary case
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            if(inputlx[ind] == inputux[ind])
                continue;
            
            #if MRQ_DEBUG_MODE
                assert( sol[ind] >= -0.1 && sol[ind] <= 1.1 ); //we are suposing binary values
            #endif
            
            auxInds[nInds] = ind;
            
            if( sol[ind] < 0.5 )
                auxValues[nInds] = 1.0; //we would round it to 0
            else
            {	
                auxValues[nInds] = -1.0;
                rhs++;
            }
            
            nInds++;
        }
        
        #if MRQ_DEBUG_MODE
            assert(nINonFixed == nInds);
        #endif
        
        r = prob->setConstraintLinearPart(neighConstrIndex, nInds, auxInds, auxValues);
        MRQ_IFERRORRETURN(r, MRQ_UNDEFINED_ERROR);
        
        r = prob->setConstraintUpperBound(neighConstrIndex, -rhs + ceil(in_integer_neighborhood_factor*nINonFixed) );
        MRQ_IFERRORRETURN(r, MRQ_UNDEFINED_ERROR);
        
    }
    else if( neighStrategy == MRQ_SNS_EUCLIDEAN_INTEGER_NEIGHBORHOOD || neighStrategy == MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD  )
    {
        //setting now the integer solution
        double maxIntDist = 0.0;
        double maxContDist = INFINITY;
        
        int r;
        double *auxVars = auxValues; 
        optsolvers::OPT_SolutionDistanceConstraintSetter esdcs;
        
        
        
        //to set the size of neighborhood, we calculate the total integer gap to guarantee there is at least one integer solution (even if it is infeasible)
        //MRQ_copyArray(n, sol, values);
        
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            if(inputlx[ind] != inputux[ind])
            {
                const double gap = MRQ_gap( sol[ind] );
                maxIntDist += gap * gap; //we are claculating the euclidian distance. So, we have to power gap to square. After this sum, we calculate the square root...
            }
        }
        
        maxIntDist = sqrt(maxIntDist); //we get square root to get the correct distance  ro closest integer solution
        
        //std::cout << "maxIntDist so com gap: " << maxIntDist << " in_integer_neighborhood_factor: " << in_integer_neighborhood_factor << "\n";
        
        
        maxIntDist += ceil( in_integer_neighborhood_factor*nINonFixed ); //here, it is like we let, at most, ceil( in_integer_neighborhood_factor*nINonFixed ) integer variables change its values from rounding solution.
        
        
        
        //std::cout << "maxIntDist final: " << maxIntDist << "\n";
        
        r = esdcs.setDistConstraint(*prob, neighConstrIndex, nI, intVars, sol, maxIntDist, auxInds, auxVars);
        MRQ_IFERRORRETURN(r, MRQ_UNDEFINED_ERROR);
        
        
        if( neighStrategy == MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD && nC > 0 )
        {
            int r, indMaxAbsVar = -1;
            double maxAbsVar = -1.0;
            
            for(int i = 0; i < nC; i++)
            {
                const int ind = contVars[i];
                if( MRQ_abs(sol[ind]) > maxAbsVar )
                {
                    indMaxAbsVar = ind;
                    maxAbsVar = MRQ_abs(sol[ind]);
                }
            }
            
            #if MRQ_DEBUG_MODE
                assert( indMaxAbsVar >= 0 && indMaxAbsVar < n );
                assert( maxAbsVar >= 0 );
            #endif
            
            maxContDist = 1.0 + maxAbsVar*in_continuous_neighborhood_factor*nC; //we sum 1.0 because all variables can de set as 0.0
            
            maxContDist = sqrt(maxContDist); //we get square root because we are speaking about euclidian distance. 
            
            //std::cout << "nC: " << nC << " factor: " << in_continuous_neighborhood_factor << " maxAbsVar: " << max << " maxContDist: " << maxContDist << "\n";
            
            r = esdcs.setDistConstraint(*prob, neighConstrIndex+1, nC, contVars, sol, maxContDist, auxInds, auxVars);
            MRQ_IFERRORRETURN(r, MRQ_UNDEFINED_ERROR);
            
        }
        
        
        #if MRQ_DEBUG_MODE
            assert(maxIntDist >= 0);//if maxIntDist is zero, we have all integer vars fixed. So, the initial soultion should be integer 
            assert(maxContDist > 0);
        #endif
        
        
        /*for(int i = 0; i < n; i++)
            std::cout << "sol["<<i<<"]: " << sol[i] << "\t";
        std::cout << "\n";*/
        
        
        //we take advatnage to reduce box constraints also...
        for(int i = 0; i < nI; i++)
        {
            int ind = intVars[i];
            
            if( newlx[ind] != newux[ind] ) //we perform this test just to avoid numerical problems if variable is fixed and solver give a solution having small decimal places
            {
                newlx[ind] = MRQ_max( newlx[ind], ceil(sol[ind] - maxIntDist) );
                
                newux[ind] = MRQ_min( newux[ind], floor(sol[ind] + maxIntDist)  );
            }
        }
        
        if( maxContDist < INFINITY )
        {
            for(int i = 0; i < nC; i++)
            {
                int ind = contVars[i];
                
                newlx[ind] = MRQ_max( newlx[ind], sol[ind] - maxContDist );
                newux[ind] = MRQ_min( newux[ind], sol[ind] + maxContDist );
            }
        }
        
    }
    else if( neighStrategy ==  MRQ_SNS_ORIGINAL)
    {
        for(int k = 0; k < nI; k++)
        {
            const int ind = intVars[k];
            
            if( MRQ_gap(sol[ind]) < in_integer_tol )
                newlx[ind] = newux[ind] = round(sol[ind]); //we fix integer var
            else
            {
                newlx[ind] = floor( sol[ind] );
                newux[ind] = ceil( sol[ind] );
                
                #if MRQ_DEBUG_MODE
                    assert( newlx[ind] >= inputlx[ind] );
                    assert( newux[ind] <= inputux[ind] );
                #endif
            }
        }
    }
    else
    {
        assert(false);
    }
    
    
    return 0;
}




int MRQ_RensExecutor::run(unsigned int thnumber, unsigned int nThreads, double insideSolverMaxTime, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, MRQ_GeneralSolverParams* algParams, const double *lx, const double *ux, const double *sol, double upper_bound, double &outObj, double *outSol)
{
    const int n = subprob->n;
    int r, retCode = MRQ_HEURISTIC_FAIL;
    
    
    resetOutput();
    
    r = MRQ_setSubproblemNeighborhood(lx, ux, sol, nI, intVars, nC, contVars, rensNeighStrategy, in_integer_neighborhood_factor, in_continuous_neighborhood_factor, in_integer_tol, neighConstrIndex, auxInds, auxValues, subprob,                            nlx, nux);
    MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
        
    
    if( in_apply_only_heuristics_on_subproblems )
    {
        MRQ_Algorithm *pHeurAlg;
        
        if(in_heurs)
        {
            heurs = in_heurs;
        }
        else
        {
            heurs = new (std::nothrow) MRQ_HeuristicExecutor;
            MRQ_IFMEMERRORGOTOLABEL(!heurs, retCode, termination);
        }
        
        if(algParams)
            heurs->setParametersToAlgs(*algParams);
        
        //trying disable a possible gens usage inside other heuristics
        heurs->setIntegerParameterToAlgs( "in_use_int_heuristic_rens", 0);
        heurs->setIntegerParameterToAlgs( "in_use_rens", 0);
        
        heurs->setNumberOfThreads(nThreads);
        heurs->setSolvers(in_milp_solver, in_nlp_solver);
        heurs->setMaxCPUTimes(in_max_cpu_time);
        heurs->setMaxTimes(in_max_time);
        
        heurs->in_stop_on_first_sol = in_stop_on_first_sol;
        
        heurs->in_use_rens = false;
        
    
        
        r = heurs->insideRun(*subprob, milpSolverParams, nlpSolverParams, upper_bound, outObj, outSol, true, &pHeurAlg, thnumber, nThreads, insideSolverMaxTime, nlx, nux);
        
        if(r == MRQ_OPTIMAL_SOLUTION || r == MRQ_HEURISTIC_SUCCESS)
        {
            retCode = MRQ_HEURISTIC_SUCCESS;
            out_alg = pHeurAlg;
        }
        
    }
    else
    {
        if(in_alg)
        {
            out_alg = in_alg;
        }
        else
        {
            out_alg = MRQ_newAlgorithm( in_algorithm_to_solve_subproblem, subprob->getNumberOfNLEqualityConstraints()  );
            MRQ_IFMEMERRORGOTOLABEL(!out_alg, retCode, termination);
            
            
             if( out_alg->isLinearApproximationAlgorithm() && subprob->getNumberOfNLEqualityConstraints() )
             {
                 MRQ_PRINTERRORMSG(" at RENS heuristic: A linear approximation algorithm was selected to solve a problem having nonlinear constraints! Please change the rens algorithm option!");
                 retCode = MRQ_BAD_DEFINITIONS;
                 goto termination;
             }
             
            
            
            
            if(algParams)
                out_alg->setParameters(*algParams);
            
            
            //we takeadvantage the nlp relaxation solution that we already have
            out_alg->in_use_initial_solution = true;
            out_alg->setInitialSolution(n, sol);
            
            if( !in_solve_continuous_relaxation_on_subproblem )
            {
                //if algorithm can skip continuous relaxation solution, we do it...
                out_alg->setIntegerParameter("in_use_first_nlp_relaxation", 0);
                
                if( out_alg->isLinearApproximationAlgorithm() )
                {
                    MRQ_LinearApproxAlgorithm *lalg = (MRQ_LinearApproxAlgorithm*) out_alg;
                    
                    int r = lalg->addPointToLinearisation(n, sol);
                    
                    if(r != 0)
                    {
                        MRQ_PRINTERRORNUMBER(r); //we do not care if we get some error
                    }
                }
                
            }
            
            //we just set parameters if user do not send is an algorithm class pointer. Otherwise, we run user algorithm in exact way that it was passed...
            
            //if algorithm choosen can make use of rens, we have to disable it
            out_alg->setIntegerParameter("in_use_int_heuristic_rens", 0);
            out_alg->setIntegerParameter("in_use_rens", 0);
            
            out_alg->in_nlp_solver = in_nlp_solver;
            out_alg->in_milp_solver = in_milp_solver;
            out_alg->in_number_of_threads = nThreads;
            out_alg->in_max_cpu_time = in_max_cpu_time;
            out_alg->in_max_time = in_max_time;
            out_alg->in_print_level = in_print_level;
            
            if( in_stop_on_first_sol )
            {
                out_alg->in_lower_bound = MRQ_INFINITY*0.1; //we set a lerge value as lower bound to force algoruthm stop in the first feasible solution found. Do not set lower bound as MRQ_INFINITY because this value will be equal to upper bound and algorithm will declare infeasibility since lower and upper bounds are equal
            }
            
        }
        
        
        
        MRQ_insideRun(out_alg, *subprob, milpSolverParams, nlpSolverParams, thnumber, insideSolverMaxTime, nlx, nux);
        
        
        if(out_alg->out_feasible_solution)
        {
            out_alg->getBestSolutionCopy(n, outSol, outObj);
        }
        
        
        
        if( in_print_level > 5 )
        {
            std::cout << MRQ_PREPRINT << "subalg ret code: " <<    out_alg->out_return_code << " subalg iters: " << out_alg->out_number_of_iterations << " subalg cpu time: " << out_alg->out_cpu_time << "\n";
        }
        
    }
    
    
    
    
termination:
    
    return retCode;
}







MRQ_RENS::MRQ_RENS():MRQ_Heuristic()
{
    resetParameters();
    resetOutput();
    out_algorithm = MRQ_RENS_HEUR_ALG;
}


MRQ_RENS::~MRQ_RENS()
{
}


int MRQ_RENS::checkAlgorithmRequirements( MRQ_MINLPProb &prob, const double *lx, const double *ux)
{
    if(in_neighborhood_strategy == MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD ) 
    {
        if( MRQ_isBinarieProblemAtRegion(prob, lx, ux) == false )
        {
            MRQ_PRINTMSG("We are sorry, but RENS with " MRQ_STR(MRQ_RENS_NS_LOCAL_BRANCHING_STRATEGY) "is only available to be applyed at binary problems. Try original RENS neighborhood strategy.");
            
            return MRQ_ALG_NOT_APPLICABLE;
        }
    }
    
    return MRQ_Algorithm::checkAlgorithmRequirements(prob, lx, ux);
}


void MRQ_RENS::printParameters(std::ostream &out) const
{
    char strValue[100];
    MRQ_Heuristic::printParameters(out);
    
    out << "\n"
    
    MRQ_STRFFATT(in_apply_only_heuristics_on_subproblem) << "\n"
    MRQ_STRFFATT(in_solve_continuous_relaxation_on_subproblem) << "\n"
    MRQ_STRFFATT(in_stop_on_first_sol) << "\n"
    MRQ_STRFFATT(in_continuous_neighborhood_factor) << "\n"
    MRQ_STRFFATT(in_integer_neighborhood_factor) << "\n";
    
    MRQ_enumToStr(in_algorithm_to_solve_subproblem, strValue);
    out << MRQ_STR(in_algorithm_to_solve_subproblem) ": " << strValue << "\n";
    
    MRQ_enumToStr(in_neighborhood_strategy, strValue);
    out << MRQ_STR(in_neighborhood_strategy) ": " << strValue << "\n"
    MRQ_STRFFATT(in_algorithm_object_to_solve_subproblem) << "\n"
    MRQ_STRFFATT(in_heuristic_exec_object_to_solve_subproblem) << "\n";
}


void MRQ_RENS::resetParameters()
{
    MRQ_Heuristic::resetParameters();
    
    in_apply_only_heuristics_on_subproblem = false;
    in_solve_continuous_relaxation_on_subproblem = true;
    in_stop_on_first_sol = true;
    
    in_integer_neighborhood_factor = 0.25;
    in_continuous_neighborhood_factor = 0.25;
    
    in_algorithm_to_solve_subproblem = MRQ_LP_NLP_BB_OA_BASED_ALG;
    in_neighborhood_strategy = MRQ_SNS_ORIGINAL;
    
    in_algorithm_object_to_solve_subproblem = NULL;
    in_heuristic_exec_object_to_solve_subproblem = NULL;
}


void MRQ_RENS::resetOutput()
{
    out_algorithm_to_solve_subproblem = MRQ_UNDEFINED_ALG;
}


int MRQ_RENS::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_Heuristic::setIntegerParameter( name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_apply_only_heuristics_on_subproblem), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_solve_continuous_relaxation_on_subproblem), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_stop_on_first_sol), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


int MRQ_RENS::setDoubleParameter(const char *name, const double value)
{
    int ret = MRQ_Heuristic::setDoubleParameter( name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    
    if( MRQ_setAtt<double>( MRQ_STRATT(in_integer_neighborhood_factor), name, value ) == 0 );
    else if( MRQ_setAtt<double>( MRQ_STRATT(in_continuous_neighborhood_factor), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


//that method should be used to set enumeration parameters
int MRQ_RENS::setStringParameter(const char *name, const char *value)
{
    int r, ret = MRQ_Heuristic::setStringParameter( name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    
    if( (r = MRQ_setStrAtt( MRQ_STRATT(in_algorithm_to_solve_subproblem), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_neighborhood_strategy), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}




int MRQ_RENS::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    return run(prob, milpSolverParams, nlpSolverParams, NULL);
}


int MRQ_RENS::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, MRQ_GeneralSolverParams* algParams)
{
    const int n = prob.n;
    //const int m = prob.m;
    
    bool updtConstrBounds;
    int ret;
    MRQ_SUBPROBLEM_NEIGHBORHOOD_STRATEGY neighStrategy = this->in_neighborhood_strategy;
    double timeStart;
    clock_t clockStart;
    
    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux;
    double *plc = NULL, *puc = NULL;
    double *initSol;
    
    MRQ_NLPSolver *nlp = NULL;
    MRQ_Preprocessor preprocessor(&prob);
    MRQ_RensExecutor rens;
    
    
    
    
    
    
    timeStart = MRQ_getTime();
    clockStart = clock();
    
    
    {
        if( in_print_parameters_values)
        {
            std::cout << MRQ_PREPRINT "Muriqui subalgorithm parameters values:\n";
            if(algParams)
                algParams->print();
            std::cout << MRQ_PREPRINT "end of Muriqui subalgorithm parameters values:\n";
        }
        
        auto ret = algorithmInitialization(1, (in_preprocess_lin_constr || in_preprocess_obj_function || in_preprocess_quad_constrs), milpSolverParams, nlpSolverParams, prob, lx, ux, &preprocessor, &updtConstrBounds, &plc, &puc); //that algorithm is monothread...
        
        if(ret != 0)
        {
            if(in_print_level > 0 )
                std::cerr << MRQ_PREPRINT "Error " << ret << " at algorithm initialization\n";
            
            out_return_code = ret;
            goto termination;
        }
    }
    
    
    if(in_print_level > 1)
        std::cout << "\n" MRQ_PREPRINT "Starting Relaxation Enforced Neighborhood Search (RENS) heuristic\n";
    
    if(in_print_level > 3)
        printSubSolvers(true, true, false);
    
    
    if( neighStrategy == MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD && MRQ_isBinarieProblemAtRegion(prob, lx, ux) == false )
    {
        MRQ_PRINTMSG("RENS with " MRQ_STR(MRQ_RENS_NS_LOCAL_BRANCHING_STRATEGY) " is only available to be applyed at binary problems. Changing neighborhood strategy to " MRQ_STR(MRQ_RENS_NS_EUCLIDEAN_INTEGER_NEIGHBORHOOD) ".\n" );
        
        neighStrategy = MRQ_SNS_EUCLIDEAN_INTEGER_NEIGHBORHOOD;
    }
    
    
    if( !in_use_initial_solution || xInit == NULL)
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
        if( updtConstrBounds )
        {
            for(int i = 0; i < m; i++)
                ret += nlp->setConstraintBounds(i, plc[i], puc[i]);
        } */
        
       
        ret = nlp->solve(false);
        if(nlp->feasSol)
        {
            if(ret == OPT_OPTIMAL_SOLUTION)
            {
                out_obj_opt_at_continuous_relax = nlp->objValue;
                if( nlp->objValue > zl )
                    zl = nlp->objValue;
                    
                if(in_print_level > 1)
                    std::cout << MRQ_PREPRINT  "NLP relaxation solution: " << nlp->objValue << "\n";
            }
            
            initSol = nlp->sol;
        }
        else if(ret == OPT_INFEASIBLE_PROBLEM)
        {
            out_return_code = MRQ_INFEASIBLE_PROBLEM;
            goto termination;
        }
        else if(ret == OPT_UNBOUNDED_PROBLEM)
        {
            out_return_code = MRQ_UNBOUNDED_PROBLEM;
            goto termination;
        }
        else
        {
            MRQ_PRINTERRORMSGP("Error to get a feasible solution to NLP relaxation. OPT return code:", ret);
            out_return_code = MRQ_NLP_SOLVER_ERROR;
            goto termination;
        }
        
    }
    else 
    {
        initSol = this->xInit;
    }
    
    if(in_print_level > 5)
    {
        std::cout << MRQ_PREPRINT  "Rens initial solution: \n";
        for(int i = 0; i < n; i++)
            std::cout << MRQ_PREPRINT  "sol["<<i<<"]: " << initSol[i] << "\n";
    }
    
    rens.in_print_level = in_print_level;
    rens.in_apply_only_heuristics_on_subproblems = in_apply_only_heuristics_on_subproblem;
    rens.in_solve_continuous_relaxation_on_subproblem = in_solve_continuous_relaxation_on_subproblem;
    rens.in_stop_on_first_sol = in_stop_on_first_sol;
    rens.in_nlp_solver = in_nlp_solver;
    rens.in_milp_solver = in_milp_solver;
    rens.in_integer_tol = in_integer_tol;
    rens.in_continuous_neighborhood_factor = in_continuous_neighborhood_factor;
    rens.in_integer_neighborhood_factor = in_integer_neighborhood_factor;
    rens.in_max_cpu_time = in_max_cpu_time;
    rens.in_max_time = in_max_time;
    rens.in_algorithm_to_solve_subproblem = in_algorithm_to_solve_subproblem;
    rens.in_alg = in_algorithm_object_to_solve_subproblem;
    rens.in_heurs = in_heuristic_exec_object_to_solve_subproblem;
    
    
    
    ret = rens.setSubProb(prob, neighStrategy);
    if( ret != 0 )
    {
        if( in_print_level > 0 )
            MRQ_PRINTERRORNUMBER(ret);
        
        out_return_code = MRQ_NLP_SOLVER_ERROR;
        goto termination;
    }
    
    
    ret = rens.run(thnumber, in_number_of_threads, MRQ_min(insideSolverMaxTime, in_max_cpu_time), milpSolverParams, nlpSolverParams, algParams, lx, ux, initSol, zu, out_best_obj, out_best_sol );
    
    out_return_code = out_best_obj < MRQ_INFINITY ? MRQ_HEURISTIC_SUCCESS : MRQ_HEURISTIC_FAIL;
    
    if( rens.out_alg )
    {
        out_algorithm_to_solve_subproblem =  rens.out_alg->out_algorithm;
        out_number_of_iterations = rens.out_alg->out_number_of_iterations;
    }
    
    
termination:
    
    if(plc)		free(plc);
    if(nlp)		delete nlp;
    
    
    out_feasible_solution = out_best_obj < MRQ_INFINITY;
    out_lower_bound = zl;
    out_upper_bound = zu;
    
    algorithmFinalization(in_number_of_threads, prob, lx, ux);
    
    
    out_cpu_time = MRQ_calcCPUTtime(clockStart);
    out_clock_time = MRQ_getTime() - timeStart;
    
    if(in_print_level > 1)
        std::cout << MRQ_PREPRINT "cpu time: " << out_cpu_time << "\n";
    
    
    return out_return_code;
}


















