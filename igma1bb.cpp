
#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <new>


#include "MRQ_igma1.hpp"



using namespace branchAndBound;
using namespace optsolvers;
using namespace muriqui;




MRQ_IGMA1BBCallbacks::MRQ_IGMA1BBCallbacks( MRQ_MINLPProb* prob, MRQ_IGMA1* igma2, MRQ_GeneralSolverParams* gapMinParams, MRQ_GeneralSolverParams* nlpParams )
{
    initialize(prob, igma2, gapMinParams, nlpParams);
}


MRQ_IGMA1BBCallbacks::~MRQ_IGMA1BBCallbacks()
{
    desallocate();
}



int MRQ_IGMA1BBCallbacks::allocateThreadStructures( const unsigned int nthreads, double *lx, double *ux)
{
    const bool preprocess = igma1->in_preprocess_lin_constr || igma1->in_preprocess_quad_constrs || igma1->in_preprocess_obj_function;
    const int n = prob->n;
    const int m = prob->m;
    const double *opuc = oplc ? &oplc[prob->m] : NULL;
    const double zu = igma1->zu;
    
    
    this->nthreads = nthreads;
    
    gapmins = new (std::nothrow) MRQ_GapMinProb[nthreads];
    if(!gapmins)
    {
        if( igma1->in_print_level > 0 )
            MRQ_PRINTMEMERROR;
        
        return MRQ_MEMORY_ERROR;
    }
    
    #if OPT_HAVE_IPOPT
        if(igma1->in_enable_gap_min_solver_premature_stoping)
        {
            if(igma1->in_gap_min_solver == MRQ_IPOPT)
            {
                ipoptInterCall = new (std::nothrow) MRQ_IpoptIntermediateCallback3[nthreads];
                
                if(!ipoptInterCall)
                {
                    if( igma1->in_print_level > 0 )
                        MRQ_PRINTMEMERROR;
                    
                    return MRQ_MEMORY_ERROR;
                }
            }
        }
    #endif
    
    prunesByObjCutActive = (unsigned int *) calloc( nthreads, sizeof(*prunesByObjCutActive) );
    nimprovments = (unsigned int*) calloc( nthreads, sizeof(*nimprovments) );
    nsubiters = (long unsigned int *) calloc( nthreads, sizeof(*nsubiters) );
    nlps = (MRQ_NLPSolver **) calloc( nthreads, sizeof(MRQ_NLPSolver *) );
    if( !prunesByObjCutActive || !nimprovments || !nsubiters || !nlps )
    {
        if( igma1->in_print_level > 0 )
            MRQ_PRINTMEMERROR;
        
        return MRQ_MEMORY_ERROR;
    }
    
    
    
    
    if( igma1->in_use_random_initial_sols )
    {
        randoms = new (std::nothrow) MRQ_Random[nthreads];
        if( !randoms )
        {
            if( igma1->in_print_level > 0 )
                MRQ_PRINTMEMERROR;
            
            return MRQ_MEMORY_ERROR;
        }
    }
    
    
    if( preprocess )
    {
        const int offs = 2*m;
        
        tplc = (double **) malloc( nthreads * sizeof(double *) );
        if(!tplc)
        {
            if( igma1->in_print_level > 0 )
                MRQ_PRINTMEMERROR;
            
            return MRQ_MEMORY_ERROR;
        }
        
        tplc[0] = (double *) malloc( offs*nthreads * sizeof(double) );
        if( !tplc[0] )
        {
            if( igma1->in_print_level > 0 )
                MRQ_PRINTMEMERROR;
            
            return MRQ_MEMORY_ERROR;
        }
        
        for(unsigned int i = 1; i < nthreads; i++)
            tplc[i] = &(tplc[0][i*offs]);
        
        
        preprocessors = new (std::nothrow) MRQ_Preprocessor[nthreads];
        
        if( !preprocessors )
        {
            if( igma1->in_print_level > 0 )
                MRQ_PRINTMEMERROR;
            
            return MRQ_MEMORY_ERROR;
        }
    }
    
    if(binSumConstrs.nbinSumConstrs > 0)
    {
        constrChoosers = new (std::nothrow) MRQ_BinSumConstrsChooser[nthreads];
        
        if( !constrChoosers )
        {
            if(igma1->in_print_level > 0)
                MRQ_PRINTMEMERROR;
            return MRQ_MEMORY_ERROR;
        }
    }
    
    
    //we try optimize cache usage setting NLP problem here. So we allocate memory for each thread together
    
    for(unsigned int i = 0; i < nthreads; i++)
    {
        int r;
        
        if(preprocessors)
        {
            preprocessors[i].initialize(prob);
            
            r = preprocessors[i].allocateMemory(n, m);
            
            if( r != 0 )
            {
                if( igma1->in_print_level > 0 )
                    MRQ_PRINTMEMERROR;
                
                return MRQ_MEMORY_ERROR;
            }
        }
        
        
        r = gapmins[i].setProblem(igma1->in_gap_min_solver, *prob, lx, ux, gapMinSolverParams, i, igma1->in_set_special_gap_min_solver_params, igma1->in_set_gap_exp_on_constr, igma1->in_set_gap_ubound_constr, 1, igma1->in_max_cpu_time, igma1->in_max_time, 0, 0);
        
        if( r != 0 )
        {
            if(igma1->in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
            return r;
        }
        
        ((OPT_NLPSolver*) gapmins[i].solver)->in_absolute_feas_tol = igma1->in_absolute_feasibility_tol;
        ((OPT_NLPSolver*) gapmins[i].solver)->in_relative_feas_tol = igma1->in_relative_feasibility_tol;
        
        
        if( oplc )
        {
            for(int j = 0; j < m; j++)
                r += gapmins[i].solver->setConstraintBounds(j, oplc[j], opuc[j]);
            
            #if OPT_DEBUG_MODE
                if( r != 0 )
                {
                    if(igma1->in_print_level > 0)
                        MRQ_PRINTERRORNUMBER(r);
                    return MRQ_NLP_SOLVER_ERROR;
                }
            #endif
        }
        
        
        
        //we just test it here because we want try allocate data structure for the same nlps[i] in contiguous positions in the memory;
        if( zu < MRQ_INFINITY )
        {
            r = gapmins[i].updateObjCutConstr( *prob, zu );
            if( r != 0 )
            {
                if(igma1->in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(r);
                return r;
            }
        }
        
        #if OPT_HAVE_IPOPT
            if( gapmins[i].solver->getSolverCode() == OPT_IPOPT && ipoptInterCall)
            {
                ( (OPT_Ipopt*) gapmins[i].solver)->setIntermediateCallbackPointer(&ipoptInterCall[i]);
            }
        #endif
        
        
        //if( igma2->in_set_special_gap_min_solver_params )
            //MRQ_setSpecialParameters( gapmins[i].solver );
        
        
        
        
        nlps[i] = OPT_newNLPSolver( igma1->in_nlp_solver );
        
        if( !nlps[i] )
        {
            if( igma1->in_print_level > 0 )
                MRQ_PRINTMEMERROR;
            return MRQ_MEMORY_ERROR;
        }
        
        
        r = MRQ_setNLPRelaxProb(*prob, lx, ux, NULL, NULL, nlps[i], true, true, true, false, i, igma1->in_set_special_nlp_solver_params, nlpSolverParams, 1, igma1->in_max_cpu_time, igma1->in_max_time, 0, 0);
        
        if( r != 0 )
        {
            if( igma1->in_print_level > 0 )
                MRQ_PRINTERRORNUMBER(r);
            return MRQ_NLP_SOLVER_ERROR;
        }
        
        if( oplc )
        {
            for(int j = 0; j < m; j++ )
                r += nlps[i]->setConstraintBounds(j, oplc[j], opuc[j] );
            
            #if OPT_DEBUG_MODE
                if( r != 0 )
                {
                    if(igma1->in_print_level > 0)
                        MRQ_PRINTERRORNUMBER(r);
                    return MRQ_NLP_SOLVER_ERROR;
                }
            #endif
        }
        
        
        #if OPT_HAVE_IPOPT
            if(ipoptInterCall)
            {
                ipoptInterCall[i].initialize(i, prob, igma1->in_integer_tol, igma1->in_absolute_feasibility_tol, igma1->in_relative_feasibility_tol, nI, intVars);
                
                const int r = ipoptInterCall[i].allocate(n, m);
                if(r != 0)
                {
                    if(igma1->in_print_level > 0)
                        MRQ_PRINTERRORNUMBER(r);
                    
                    return r;
                }
                
                ipoptInterCall[i].nlp = nlps[i];
                
                ipoptInterCall[i].ipopt = (OPT_Ipopt*) gapmins[i].solver;
                
                
                ipoptInterCall[i].minimumItersToabortIfNoprogress = igma1->in_min_number_of_iters_on_gap_min_premature_stop;
                ipoptInterCall[i].minImprovToObj = igma1->in_min_improv_to_obj_on_gap_min_premature_stop;
                ipoptInterCall[i].minImprovToInfeas = igma1->in_min_improv_to_infeas_on_gap_min_premature_stop;
            }
        #endif
        
        if(binSumConstrs.nbinSumConstrs > 0)
        {
            const int r = constrChoosers[i].reallocate( binSumConstrs.nbinSumConstrs);
            
            if(r != 0)
            {
                if(igma1->in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(r);
                return r;
            }
        }
        
    }
    
    if( igma1->in_use_random_initial_sols )
    {
        for(unsigned int i = 0; i < nthreads; i++)
            randoms[i].setSeed( &igma1->in_seed_to_random_numbers );
    }
    
    
    
    
    return 0;
}


void MRQ_IGMA1BBCallbacks::desallocate()
{
    MRQ_secFree(intVars);
    MRQ_secDeleteArray(gapmins);
    MRQ_secDeleteArray(preprocessors);
    MRQ_secDeleteArray(randoms);
    MRQ_secDeleteArray(constrChoosers);
    #if OPT_HAVE_IPOPT
        MRQ_secDeleteArray(ipoptInterCall);
    #endif
    
    MRQ_secFree(prunesByObjCutActive);
    MRQ_secFree(nimprovments);
    MRQ_secFree(nsubiters);
    MRQ_secFree(sumGapObj);
    MRQ_secFree(nGapObj);
    
    if( tplc )
    {
        if( tplc[0] )
            free(tplc[0]);
        
        free(tplc);
        tplc = NULL;
    }
    
    if( nlps )
    {
        
        for(unsigned int i = 0; i < nthreads; i++)
        {
            if( nlps[i] )
                delete nlps[i];
        }
        
        free(nlps);
        nlps = NULL;
    }
    
    binSumConstrs.deallocate();
}


void MRQ_IGMA1BBCallbacks::initialize( MRQ_MINLPProb *prob, MRQ_IGMA1 *igma2, MRQ_GeneralSolverParams *gapMinParams, MRQ_GeneralSolverParams *nlpParams)
{
    this->prob = prob;
    this->igma1 = igma2;
    
    gapMinSolverParams = gapMinParams;
    nlpSolverParams = nlpParams;
    
    oplc = NULL;
    sumGapObj = NULL;
    nGapObj = NULL;
    
    constrBranchStrat = MRQ_BB_CBS_NO_CONSTRAINT_BRANCH;
    nI = 0;
    nnonimprovs = 0;
    intVars = NULL;
    
    
    
    #if OPT_HAVE_IPOPT
        ipoptInterCall = NULL;
    #endif
    
    prunesByObjCutActive = NULL;
    nimprovments = NULL;
    nsubiters = NULL;
    
    tplc = NULL;
    gapmins = NULL;
    nlps = NULL;
    preprocessors = NULL;
    randoms = NULL;
    constrChoosers = NULL;
}


int MRQ_IGMA1BBCallbacks::beforeAll(const unsigned int numberOfThreads, double *lx, double *ux)
{
    //const int n = prob->n;
    const int m = prob->m;
    
    
    int r;
    
    abseps = igma1->in_absolute_obj_cut_eps;
    releps = igma1->in_relative_obj_cut_eps;
    
    nnonimprovs = 0;
    
    
    nI = prob->getNumberOfIntegerVars();
    
    intVars = (int *) malloc( nI * sizeof(int) );
    if( !intVars )
    {
        if(igma1->in_print_level > 0)
            MRQ_PRINTMEMERROR;
        return BBL_MEMORY_ERROR;
    }
    
    prob->getIntegerIndices(intVars);
    
    
    
    if( igma1->in_gap_min_obj_strategy == MRQ_IGMA_GMOS_BY_GAP_AVERAGE )
    {
        sumGapObj = (double *) calloc( nI, sizeof(double) );
        nGapObj = (unsigned int *) calloc( nI, sizeof(unsigned int) );
        
        if( !sumGapObj || !nGapObj )
        {
            if(igma1->in_print_level > 0)
                MRQ_PRINTMEMERROR;
            return BBL_MEMORY_ERROR;
        }
    }
    
    if( igma1->in_constr_branching_strategy != MRQ_BB_CBS_NO_CONSTRAINT_BRANCH )
    {
        const double *plc = oplc;
        const double *puc = oplc ? &(plc[m]) : NULL;
        
        int r;
        
        
        r = binSumConstrs.calculateIndices(*prob, lx, ux, plc, puc);
        
        if( r != 0 )
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
        
        constrBranchStrat = igma1->in_constr_branching_strategy;
        
        if( constrBranchStrat == MRQ_BB_CBS_STBRANCH_PSEUDO_COSTS )
        {
            if(igma1->in_print_level)
                MRQ_PRINTERRORMSG("Warning: IGMA does not compute pseudo custs. So " MRQ_STR(MRQ_BB_CBS_STBRANCH_PSEUDO_COSTS) "is not able as constraint strategy. Adopting " MRQ_STR(MRQ_BB_CBS_LOWEST_NUMBER_OF_VARS) ".");
            
            constrBranchStrat = MRQ_BB_CBS_LOWEST_NUMBER_OF_VARS;//we do not have pseudo-custs in igma2. so, we cannot adopt this strategy
        }
        
        
        if( igma1->in_print_level > 5 )
        {
            std::cout << MRQ_PREPRINT "Constraints of binary sum:\n";
            for(int j = 0; j < binSumConstrs.nbinSumConstrs; j++)
                std::cout << "c["<<j<<"]: " << binSumConstrs.binSumConstrs[j] << " \t";
            std::cout << "\n";
        }
    }
    
    
    r = allocateThreadStructures( numberOfThreads, lx, ux );
    if( r != 0)
    {
        if(igma1->in_print_level > 0)
            MRQ_PRINTERRORNUMBER(r);
        return r;
    }
    
    
    return 0;
}


int MRQ_IGMA1BBCallbacks::beforeSolvingRelaxation( const unsigned int threadNumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, bool &pruneNode)
{
    const int m = prob->m;
    bool updtvb, updtcb;
    int r;
    
    //we assume if this method is called, is because we want perform a preprocessing.
    const double *opuc = oplc ? &oplc[m] : NULL;
    //const double *opuc = &oplc[m];
    
    
    double *myplc = tplc[threadNumber];
    double *mypuc = &myplc[m];
    
    MRQ_Preprocessor &preprocessor = preprocessors[threadNumber];
    
    MRQ_GapMinProb &gapminsolver = gapmins[threadNumber];
    MRQ_NLPSolver *gapmin = (MRQ_NLPSolver *) gapminsolver.solver;
    MRQ_NLPSolver *nlp = nlps[threadNumber];
    
    pruneNode = false;
    
    
    //std::cerr << MRQ_PREPRINT << " preprocessando." << std::endl;
    
    
    r = preprocessor.preprocess( igma1->in_preprocess_quad_constrs, igma1->in_preprocess_obj_function, ub, nlx, nux, updtvb, updtcb, oplc, opuc, myplc, mypuc );
    
    if( r == minlpproblem::MIP_INFEASIBILITY )
    {
        //if( igma2->in_print_level > 4 )
            //MRQ_PRINTMSG("Podei no por preprocessamento\n");
        
        pruneNode = true;
        return 0;
    }
    
    if( igma1->in_enable_gap_min_solver_premature_stoping )
    { //we have to update variable bounds because preprocessor can chnge some continuous variable bound and we can get a wrong solution. Note, if in_enable_gap_min_solver_premature_stoping, we only update variable bounds if we really need solve local search problem
        const int n = prob->n;
        const int r = nlp->setnVariablesBounds(n, nlx, nux);
        
        #if MRQ_DEBUG_MODE
            if( r != 0 )
            {
                MRQ_PRINTERRORNUMBER(r);
                return MRQ_NLP_SOLVER_ERROR;
            }
        #endif
    }
    
    //changing constraints bounds in knitro forces problem be reinput from the beginning. So, it is better do not update...
    if( igma1->in_gap_min_solver != MRQ_NLP_KNITRO )
    {
        for( int i = 0; i < m; i++ )
        {
            r += gapmin->setConstraintBounds(i, myplc[i], mypuc[i] );
            r += nlp->setConstraintBounds(i, myplc[i], mypuc[i] );
        }
        
        #if MRQ_DEBUG_MODE
            if( r != 0 )
            {
                MRQ_PRINTERRORNUMBER(r);
                return MRQ_NLP_SOLVER_ERROR;
            }
        #endif
    }
    
    
    //std::cerr << MRQ_PREPRINT << " preprocessei." << std::endl;
    
    return 0;
}



int MRQ_IGMA1BBCallbacks::solveSubProblem( const unsigned int threadNumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, BBL_RETURN_CODES &retCode, double &objValue, double &dualObjValue, double *sol, double *dualSol, bool &generalFeasibleSol, bool &pruneNode, double &nodeLowerBound, BBL_BRANCH_STRATEGY &branchStrategy)
{
    const int n = prob->n;
    const int m = prob->m;
    const int MAXCALLERRORS = 10;
    //bool intSol;
    int r, ncallerror = 0;
    double zu, zucut;
    
    
    MRQ_GapMinProb &gapminsolver = gapmins[threadNumber];
    MRQ_Random &random = randoms[threadNumber];
    MRQ_NLPSolver *nlp = nlps[threadNumber];
    MRQ_NLPSolver *gapmin = (MRQ_NLPSolver *) gapminsolver.solver;
    
    
    retCode = BBL_SOLVER_ERROR;
    generalFeasibleSol = false;
    pruneNode = false;
    
    
    
    //abseps = igma2->in_abs_obj_cut_eps;
    //releps = igma2->in_rel_obj_cut_eps;
    
    
    r = gapmin->setnVariablesBounds(n, nlx, nux);
    #if MRQ_DEBUG_MODE
        if( r != 0 )
        {
            if( igma1->in_print_level > 0 )
                MRQ_PRINTERRORNUMBER(r);
            return MRQ_NLP_SOLVER_ERROR;
        }
    #endif
    
    
    
    if( igma1->in_gap_min_obj_strategy == MRQ_IGMA_GMOS_BY_GAP_AVERAGE )
    {
        const double factor = igma1->in_factor_to_min_gap_obj_by_avg_gap;
        
        //updating objective function
        for( int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            if(  nlx[ind] == nux[ind] )
            {
                sol[i] = 0.0;
                continue;
            }
            
            sol[i] = nGapObj[i] > 0 ? 1.0 + factor*(sumGapObj[i]/nGapObj[i]) : 1.0;
        }
        
        
        
        r = gapmin->setObjLinearCoefs(nI, intVars, sol);
        #if MRQ_DEBUG_MODE
            if( r != 0 )
            {
                if( igma1->in_print_level > 0 )
                    MRQ_PRINTERRORNUMBER(r);
                return BBL_SOLVER_ERROR;
            }
        #endif
        
        
        #pragma ivdep
        #pragma GCC ivdep
        for( int i = 0; i < nI; i++ )
            sol[i] = -2.0*sol[i];
    }
    else 
    {
        #if MRQ_DEBUG_MODE
            assert( igma1->in_gap_min_obj_strategy == MRQ_IGMA_GMOS_SAME_WEIGHT );
        #endif
        
        for( int i = 0; i < nI; i++ )
        {
            const int ind = intVars[i];
            sol[i] = nlx[ind] != nux[ind] ? 1.0 : 0.0;
        }
        
        r = gapmin->setObjLinearCoefs(nI, intVars, sol);
        #if MRQ_DEBUG_MODE
            if( r != 0 )
            {
                if( igma1->in_print_level > 0 )
                    MRQ_PRINTERRORNUMBER(r);
                return BBL_SOLVER_ERROR;
            }
        #endif
        
        
        #pragma ivdep
        #pragma GCC ivdep
        for(int i = 0; i < nI; i++)
        {
            if( sol[i] != 0.0 )
                sol[i] = -2.0;
        }
        
    }
    
    
    //const int solverCode = gapmin->getSolverCode();
    if( gapmin->isMyNLPClass() )
    {
        //ok, it is not elegant, but we change coefficients directly. Note since we dot not change sparse matrix structure, I think there is no problem about it. It is so much faster than reallocate memory and recalculate hessian indices
        
        const int cindex = gapminsolver.intGapConstrIndex;
        OPT_MyNLPSolver* mynlp = (OPT_MyNLPSolver*) gapmin;
        
        
        minlpproblem::MIP_SparseMatrix &objQ = cindex < 0 ? mynlp->prob.Q : mynlp->prob.QC[cindex];
        
        //Q.printSparseMatrix();
        
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            //if binary vars are fixed from the beginning (e.g preprocessing or original definition), we do not set objective function
            if( objQ.getNumberOfElementsAtRow(ind) > 0 )
            {
                #if MRQ_DEBUG_MODE
                    assert( objQ[ind][0] == ind ); //assert( (int) Q[ind][0].getColumn() == ind );
                #endif
                
                objQ(ind)[0] = sol[i]; //Q[ind][0].setValue( sol[i] );
            }
            #if MRQ_DEBUG_MODE
            else
            {
                assert( sol[i] == 0.0 ); //integer variable has no quadratic term. So, it must be fixed and coefficient must be zero
            }
            #endif
        }  
        
    }
    else
    {
        const int cindex = gapminsolver.intGapConstrIndex;
        int r;
        
        if( cindex < 0 )
            r = gapmin->setObjQuadMatrix( nI, intVars, intVars, sol);
        else
            r = gapmin->setConstraintQuadMatrix(cindex, nI, intVars, intVars, sol);
        
        if( r != 0 )
        {
            if( igma1->in_print_level > 0 )
                MRQ_PRINTERRORNUMBER(r);
            return BBL_SOLVER_ERROR;
        }
    }
    
    //node.print();
    //gapmin->generateModelFile( "modelo.lp" );
    
    /*for(int i = 0; i < nI; i++)
        std::cout << "sol["<<i<<"]: " << sol[i] << " sumGapObj["<<i<<"]: " << sumGapObj[i] << " nGapObj["<<i<<"]: " << nGapObj[i] << std::endl; */
    
    //MRQ_getchar();
    
    
    //( (OPT_MyNLPSolver*) gapmin)->prob.print();
    //MRQ_getchar();
    
    
    zu = ub;
    
    
    while(true)
    {
        nsubiters[threadNumber]++;
        
        if( nthreads > 1 )
            zu = getUpperBound(); //another thread could update zu
        
        if(zu < BBL_INFINITY && nthreads != 1 && (igma1->in_adopt_obj_cut || !igma1->in_heuristic_mode) )
        {
            zucut = MRQ_zuWithTol(zu, abseps, releps); //zu -abseps -MRQ_abs( zu *releps );
            
            r = gapminsolver.updateObjCutConstr( *prob, zucut );
            if( r != 0 )
            {
                if( igma1->in_print_level > 0 )
                    MRQ_PRINTERRORNUMBER(r);
                return BBL_SOLVER_ERROR;
            }
            
            
            //gapmin->generateModelFile( "gapmin.mod" );
            //MRQ_getchar();
        } 
        
        
        if( igma1->in_use_random_initial_sols )
        {
            const double rlb = igma1->in_lower_bound_to_random_sol;
            const double rub = igma1->in_upper_bound_to_random_sol;
            const double boxSize = rub - rlb;
            
            //std::cout << "random.seed: " << random.getSeed() << "\n";
            //MRQ_getchar();
            
            
            for(int i = 0; i < n; i++)
            {
                const double lxi = nlx[i];
                const double uxi = nux[i];
                double l, u;
                
                
                if(lxi <= 0.0 && 0.0 <= uxi)
                {
                    l = MRQ_max(lxi, rlb);
                    u = MRQ_min(uxi, rub);
                }
                else if(lxi > -MRQ_INFINITY)
                {
                    l = lxi;
                    u = MRQ_min(uxi, lxi + boxSize);
                }
                else
                {
                    #if MRQ_DEBUG_MODE
                        assert(uxi < MRQ_INFINITY); //if we are here, lxi <= -MRQ_INFINITY. If uxi >= MRQ_INFINITY, we should be in the if
                    #endif
                    u = uxi;
                    l = uxi - boxSize; //here, we know lxi is -MRQ_INFINITY
                }
                
                
                
                //l = nlx[i] < rlb ? MRQ_min(0.0, nux[i]) - rlb : nlx[i];
                //u = nux[i] > rub ? MRQ_max(0.0, nlx[i]) + rub : nux[i];
                
                sol[i] = l == u ? l : random.random(l, u);
                
                //std::cout << "isol[" << i << "]: " << sol[i] << "\t";
                
                //printf("i: %d isol: %0.12f nlx: %0.3f nux: %0.3f\n", i, sol[i], nlx[i], nux[i]);
            }
            
            //std::cout << std::endl;
            
            gapmin->setInitialSolution( sol, NULL, NULL);
        }
        
        
        #if OPT_HAVE_IPOPT
        if(ipoptInterCall)
        {
            ipoptInterCall[threadNumber].zu = zu;
            ipoptInterCall[threadNumber].nlpSolved = false;
        }
        #endif
        
        
        //gapmin->generateModelFile("igma2.lp");
        //MRQ_getchar();
        
        //std::cout << "Vou resolver gapmin iter: " << iter << "\n\n";
        
        
        r = gapmin->solve(true);
        
        //for(int i = 0; i < n; i++)
            //printf("sol[%d]: %0.12f\n", i, gapmin->sol[i]);
            //std::cout << "sol["<<i<<"]: " << gapmin->sol[i] << " \t"  << std::endl;
        
        //std::cout << std::endl;
        
        {
            int nonfixedvars = 0;
            
            for(int i = 0; i < nI; i++)
            {
                const int ind = intVars[i];
                
                if( nlx[ind] != nux[ind] )
                    nonfixedvars++;
            }
            
            
            //std::cout << "\n   igma2 iter: " << iter << " Min Gap solving - ret: " << r << " feas sol: " << gapmin->feasSol << " obj: " << gapmin->objValue <<  " origRetCode: " <<  gapmin->origSolverRetCode  << " nI: " << nI << " non fixed: " << nonfixedvars << " maxgap: " << 0.25*nI	<< "\n";
            /*if( gapminsolver.objCutIndex >= 0 )
            {
                const int objCutInd = gapminsolver.objCutIndex;
                double clb, cub;
                
                gapmin->getConstraintBounds(objCutInd, clb, cub);
                
                std::cout << " real obj: " << gapmin->constr[objCutInd] << " cub: " << cub;
                std::cout << "\n";
            }*/
            
        }
        
        
        //MRQ_getchar();
        
        if( gapmin->feasSol || gapmin->retCode == OPT_STOP_REQUIRED_BY_USER )
        {
            
            if( MRQ_isIntegerSol(nI, intVars, gapmin->sol, igma1->in_integer_tol) )
            {
                int r;
                double obj;
                double *psol;
                bool increaseObjCut = false;
                
                /*{
                    double *sol = gapmin->sol;
                    
                    std::cout << "Solucao inteira encontrada:\n";
                    for(int i = 0; i < nI; i++)
                        printf("x[%d]: %0.1f\t", intVars[i], sol[intVars[i]]);
                    std::cout << "\n";
                    
                    double obj = NAN;
                    prob->objEval(threadNumber, true, sol, obj);
                    
                    std::cout << "Objetivo: " << obj << "\n";
                }*/
                
                
                nlp->setnVariablesBounds(n, nlx, nux);
                
                MRQ_fixIntVarsOnSolByList( nI, intVars, gapmin->sol, *nlp );
                
                nlp->setInitialSolution( gapmin->sol, gapmin->dualSolC, gapmin->dualSolV );
                
                
                //std::cout << "Vou resolver a busca local\n";
                
                r = nlp->solve();
                
                
                //std::cout << "\t\t\tBusca local. retCode: " << nlp->retCode << " feasSol: " << nlp->feasSol << " origCode: " << nlp->origSolverRetCode << " obj: " << nlp->objValue << "\n";
                
                /*for(int i = 0; i < n; i++)
                    std::cout << "sol["<<i<<"]: " << nlp->sol[i] << " \t"; //std::endl;
                std:: cout << std::endl; */
                
                //MRQ_getchar();
                
                if( nlp->feasSol )
                {
                    obj = nlp->objValue;
                    psol = nlp->sol;
                }
                else if( nlp->retCode == OPT_INFEASIBLE_PROBLEM)
                {
                    //we fix integer variables and problem gets infeasible. Gap min problem found a integer feasible solution due to integer tolerance. We break this loop to move on and branch this bb node. Although gap minimization found an "integer" solution in the practice, it is like it was a infeasible solution. Note this does not mean gap minimization problem is infeasible in this node.
                    if( igma1->in_print_level > 5 )
                        std::cout << MRQ_PREPRINT "Infeasible problem in the local search solving.\n";
                    break;
                }
                else if( gapmin->feasSol ) //if( !igma2->in_enable_gap_min_solver_premature_stoping || (igma2->in_integer_tol_on_gap_min_premature_stop <= igma2->in_integer_tol ) )
                {
                    //we use the gap minimization solution to update best solution. But we just do it, if we do not use premature stop on gapmin solving, or the integer tol to premature stop is lower then general integer tol.
                    
                    //std::cout << "\t\t\tTentando atualizar solucao pela solucao do problema de minimizacao de gap\n";
                    
                    psol = gapmin->sol;
                    
                    if( gapminsolver.objCutIndex >= 0 )
                    {
                        //we have the objective cut. So, objective function is already evaluated 
                        obj = gapmin->constr[gapminsolver.objCutIndex];
                    }
                    else
                    {
                        r = prob->objEval(threadNumber, true, psol, obj);
                        if( r != 0 )
                        {
                            if( igma1->in_print_level > 0 )
                                MRQ_PRINTCALLBACKERRORNUMBER(r);
                            break;
                        }
                    }
                    
                }
                else
                {
                    break;
                }
                
                
                //std::cout << MRQ_PREPRINT "Encontrei solucao viavel. Obj: " << obj << " current zu: " << zu << std::endl; 
                //MRQ_getchar();
                
                //printf("iter: %d obj: %0.16f zu: %0.16f \n", (int) iter, obj, zu);
                
                if(obj < zu)
                {
                    //printf("\t\t\tcai aqui! obj: %0.12lf zu: %0.12f\n", obj, zu);
                    
                    if( zu - obj < 1e-6 )
                    {
                        increaseObjCut = true;
                    }
                    
                    zu = obj; //update our zu
                    const bool updt = tryUpdateBestSolution( threadNumber, psol, obj, iter ); //update zu to BBL
                    
                    if( updt )
                    {
                        unsigned int nimprovs = 0;
                        nimprovments[threadNumber]++;
                        
                        for(unsigned int i = 0; i < nthreads; i++)
                            nimprovs += nimprovments[i];
                        
                        nnonimprovs = 0;
                        
                        if( nimprovs >= igma1->in_max_improvments_of_best_sol )
                        {
                            return MRQ_HEURISTIC_SUCCESS; //to finish algorithm, we return a value different of zero
                        }
                        
                    }
                    else
                    {
                        increaseObjCut = true;
                    }
                    
                }
                else if( igma1->in_heuristic_mode && !igma1->in_adopt_obj_cut )
                {
                    nnonimprovs++;
                    
                    //std::cout << "number of nonimprovmts: " << nnonimprovs << "\n";
                    
                    if( nnonimprovs > igma1->in_max_nonimprovment_integer_sols )
                    {
                        if(igma1->in_print_level > 4)
                            MRQ_PRINTMSG("IGMA 2: Stoping by maximum number of nonpimrovment intger sols");
                        
                        return MRQ_HEURISTIC_SUCCESS; //to finish algorithm, we return a value different of zero
                    }
                    
                }
                else if(gapmin->retCode == OPT_STOP_REQUIRED_BY_USER)
                { //we must avoid algorithm to repeat the same integer solution without any improvment
                    break;
                }
                else //if( obj == zu ) //note zu is a local variable
                {
                    increaseObjCut = true;
                }
                
                
                if(increaseObjCut)
                {
                    //objective cut was not effective;
                    if( igma1->in_print_level > 4 )
                        MRQ_PRINTMSG("Objective cut was not effective to cut current best solution. Increasing absolute and relative epsilon to objective cut\n");
                    
                    //MRQ_getchar();
                    
                    if( MRQ_abs(zu) > 1e-4 ) //we consider 1e-4 as zero for zu
                    {
                        abseps *= igma1->in_factor_to_increase_abs_obj_cut_eps;
                        releps *= igma1->in_factor_to_increase_rel_obj_cut_eps;
                    }
                    else
                    { //we have a zero as zu... so we have to be more aggressive with absolute eps tolerance. Note changing relative eps tolarence does not bring any results...
                        abseps *= igma1->in_factor_to_increase_abs_obj_cut_eps_on_zero;
                    }
                }
                
                
                if( nthreads == 1 && (igma1->in_adopt_obj_cut || !igma1->in_heuristic_mode) )
                { //if we have more than 1 thread, we update in the beggining of loop
                    zucut = MRQ_zuWithTol(zu, abseps, releps); //zu -abseps -MRQ_abs( zu *releps );
                    
                    //std::cout << "zucut: " << zucut << "\n";
                    
                    const int r = gapminsolver.updateObjCutConstr( *prob, zucut);
                    if( r != 0 )
                    {
                        if( igma1->in_print_level > 0 )
                            MRQ_PRINTERRORNUMBER(r);
                        return BBL_SOLVER_ERROR;
                    }
                    
                }
                
                
                /*if( igma2->in_stop_on_first_sol )
                {
                    return MRQ_HEURISTIC_SUCCESS; //to finish algorithm, we return a value different of zero
                } */
                
                
                //std::cout << " Novo zu: " << getUpperBound() << std::endl;
            }
            else
            {
                break;
            }
            
            
            if(igma1->in_heuristic_mode && !igma1->in_adopt_obj_cut)
            {
                break;
            }
            
            
        }
        else if( r == OPT_CALLBACK_FUNCTION_ERROR )
        {
            //we try again from another initial solution...
            //continue
            
            ncallerror++;
            
            if( ncallerror >= MAXCALLERRORS )
                break;
        }
        else
        {
            if( r == OPT_INFEASIBLE_PROBLEM )
                pruneNode = true;
            
            break;
        }
        
        
        if(igma1->in_max_cpu_time < INFINITY)
        {
            double cpuTime = MRQ_calcCPUTtime( clockStart, clock());
            
            if( cpuTime >= igma1->in_max_cpu_time )
            {
                return MRQ_MAX_TIME_STOP;
            }
        }
        
        if(igma1->in_max_time < INFINITY)
        {
            double wallTime = MRQ_getTime() - timeStart;
            
            if( wallTime >= igma1->in_max_time )
            {
                return MRQ_MAX_TIME_STOP;
            }
        }
        
    }
    
    
    r = gapmin->retCode;
    if( r == OPT_OPTIMAL_SOLUTION )
    {
        retCode = BBL_OPTIMAL_SOLUTION;
    }
    else if( gapmin->feasSol )
    {
        retCode = BBL_FEASIBLE_SOLUTION;
    }
    else if( r == OPT_INFEASIBLE_PROBLEM ) 
    {
        retCode = BBL_INFEASIBLE_PROBLEM;
    }
    else
    {
        //otherwise, retCode will be BBL_UNDEFINED_ERROR. Note, the node will be consider infeasible dependig of parameter in_consider_relax_infeas_if_solver_fail
        dualObjValue = nI; //we add nI to the gap value. In this way, nodes where we reached objective cut will be always after others in B&B open nodes list
    }
    
    
    if( gapmin->feasSol )
    {
        dualObjValue = gapmin->objValue;
        
        if( gapminsolver.objCutIndex >= 0 )
        {
            double zucut, zlb;
            const double absTolToObjCutActive = igma1->in_absolute_obj_tol_to_active_obj_cut; //1.0e-2;
            const double relTolToObjCutActive = igma1->in_relative_obj_tol_to_active_obj_cut;
            const double tolToObjCutActive = absTolToObjCutActive + MRQ_abs(zu*relTolToObjCutActive);
            
            //we check if objective cut is active. If yes, we add nI to the gap value. In this way, nodes where we reached objective cut will be always after others in B&B open nodes list. We belive if objective cut is active here, that can be a signal there is no solution in this partition better than zu
            
            const double obj = gapmin->constr[ gapminsolver.objCutIndex ];
            
            //const double lambdaoc = gapmin->dualSolC[gapminsolver.objCutIndex];
            
            
            gapmin->getConstraintBounds( gapminsolver.objCutIndex, zlb, zucut);
            
            
            //std::cout << "\n*******checando poda por corte de nivel obj. obj: " << obj << " zu: " << zu << " zucut: " << zucut << " tolToObjCutActive: " << tolToObjCutActive << " dif: " << zucut - obj << " poda: " << (zucut - obj <= tolToObjCutActive) << " lambda obj cut: " << lambdaoc <<  "\n";
            
            //if we do not solve until optimallity, we can have lambda nonzero even if constraint is not active...
            //if( MRQ_abs( gapmin->dualSolC[gapminsolver.objCutIndex] ) >= 1e-6) 
            if( zucut - obj <= tolToObjCutActive )
            {//just checking if dual is diferent of zero
                
                //std::cout << "Detectei corte objetivo ativo. zu: " << zu << " zucut: " << zucut << " obj: " << obj << " lambda: " << gapmin->dualSolC[gapminsolver.objCutIndex] << MRQ_GETFILELINE << std::endl;
                
                //MRQ_getchar();
                
                
                
                //assert( MRQ_abs(obj - zucut) <= tolToObjCutActive ); //we use MRQ_abs because definitin of lagrange can change depending on solver. For example, lagrange multipliers on ipopt has oposite signal compared to mosek ones...
                
                //assert( MRQ_abs( gapmin->dualSolC[gapminsolver.objCutIndex] ) >= 1e-6 ); //just checking if dual is diferent of zero
                
                
                
                if( igma1->in_heuristic_mode )
                {
                    prunesByObjCutActive[threadNumber] += 1;
                    
                    pruneNode = true; //in heuristic mode we are not interested to prove optimality. If we are here, we have already a feasible solution and the gapmin found a solution in the objective cut. So, we discard this node assuming there is no a better solution in this partition (actually it is not necessarilly true, but we are in the heuristic mode...)
                }
                
                dualObjValue += nI; //we add nI to the gap value. In this way, nodes where we reached objective cut will be always after others in B&B open nodes list
            }
        }
        
        
        
        MRQ_copyArray( n, gapmin->sol, sol );
        //MRQ_copyArray( m, gapmin->dualSolC, dualSol );
        //MRQ_copyArray( 2*n, gapmin->dualSolV, &dualSol[m]);
        
        
        if( igma1->in_gap_min_obj_strategy == MRQ_IGMA_GMOS_BY_GAP_AVERAGE )
        {
            updateSumGapObj( sol, nlx, nux );
        }
        
    }
    else
    {
        //if gapmin->retCode is OPT_STOP_REQUIRED_BY_USER is because ipopt prematurelly due to lack of progress (or we already found an integer solution). So, it is not correct prune here.
        
        if( igma1->in_consider_relax_infeas_if_solver_fail && gapmin->retCode != OPT_STOP_REQUIRED_BY_USER )
            pruneNode = true;
    }
    
    
    if( binSumConstrs.nbinSumConstrs > 0 && constrBranchStrat != MRQ_BB_CBS_NO_CONSTRAINT_BRANCH )
    {
        if( !pruneNode && retCode != BBL_INFEASIBLE_PROBLEM )
        {
            MRQ_BinSumConstrsChooser &constrChooser = constrChoosers[threadNumber];
            const double *plc = preprocessors ? tplc[threadNumber] : prob->lc;
            const double *puc = preprocessors ? &plc[m] : prob->uc;
            
            if( constrChooser.calculateCandidateConstrsToBranch( binSumConstrs.nbinSumConstrs, binSumConstrs.binSumConstrs, nlx, nux, prob->xtype, prob->A, plc, puc, igma1->in_integer_tol, sol ) > 0 )
                branchStrategy = BBL_BS_USER_NODE_GENERATION;
            else
                branchStrategy = origBBLBranchStrategy;
        }
    }
    
    
    
    
    #if 0
        std::cout << "retCode: " << retCode << " pruneNode: " << pruneNode << " objValue: " << objValue << " dualObjValue: " << dualObjValue << "\n";
        for(int i = 0; i < n; i++)
            printf("sol[%d]: %0.20f   ", i, sol[i]);//std::cout << "sol["<<i<<"]: " << sol[i] << " \t";
        std::cout << "\n";
    #endif
    
    //std::cout << "dualObjValue: " << dualObjValue << std::endl;
    
    
    //MRQ_getchar();
    
    
    return 0;
}



int MRQ_IGMA1BBCallbacks::chooseIndexToBranch( const int threadNumber, BBL_Node& node, const long unsigned int iter, const double lb, const double ub, double* nlx, double* nux, BBL_RETURN_CODES retCode, double objValue, double dualObjValue, double* sol, double* dualSol, unsigned int& sizeIndices, unsigned int* indices, double* breakValues1, double* breakValues2, BBL_Node*& nodes)
{
    //const int nI = this->nI;
    int indmax = -1;
    double maxgap = -1.0;
    
    const double gapTreeshould = 0.0;//1e-8; /* that is really strange, but something strange is happening on ipopt when igma2 finds a feasible solution minimizing gap. Follower min gap solving can give a solution having differences after 12 decimal places in different runnings, even I providing the same initial solution in the same partition... so, we try prevent algorithm get random enforcing this thressolh on maxgap calculation */
    
    
    //we assume all variables here are binary. If it changes some day, you have to fix this part of code! Do not forget breakValues1 and breakValues2
    
    
    for(int i = 0; i < nI; i++)
    {
        const int ind = intVars[i];
        
        //we can have general integer variables fixed...
        if( nlx[ind] == nux[ind] )
            continue;
        
        const double gap = MRQ_min( sol[ind], 1.0 - sol[ind] );
        
        
        if( gap - maxgap > gapTreeshould ) //it should be gap - maxgap > 0.0, but ipopt is giving crazy results...
        {
            maxgap = gap;
            indmax = ind;
        }
        
        //printf("i: %d ind: %d sol: %f gap: %0.20f maxgap: %f indmax: %d\n", i, ind, sol[ind], gap, maxgap, indmax );
        
        
        #if MRQ_DEBUG_MODE
            assert( gap >= 0.0 );
        #endif
        
    }
    
    #if MRQ_DEBUG_MODE
        assert( maxgap <= 0.5 ); // if variables are binary, this test should pass
        
        if( indmax == -1)
        {// we have a integer soluiton. It can occurs if we abort ipopt solving because lack of progrres and we have a increadible coincidence ipopt was in integer (infeasible) solution. we chhose the first integer variable nonfixed
            
            for(int i = 0; i < nI; i++)
            {
                const int ind = intVars[i];
                
                if( nlx[ind] != nux[ind] )
                {
                    indmax = ind;
                    break;
                }
            }
            
            #if MRQ_DEBUG_MODE
                
                #if 0
                {
                    const int n = prob->n;
                    MRQ_GapMinProb &gapminsolver = gapmins[threadNumber];
                    MRQ_NLPSolver2 *nlp = nlps[threadNumber];
                    MRQ_NLPSolver2 *gapmin = (MRQ_NLPSolver2 *) gapminsolver.solver;
                    
                    
                    
                    std::cout << "gapmin retCode: " << gapmin->retCode << " original retCode: " << gapmin->origSolverRetCode << " feasSol:" << gapmin->feasSol << " obj: " << gapmin->objValue << "\n";
                    
                    for(int i = 0; i < n; i++)
                        std::cout << "sol["<<i<<"]: " << gapmin->sol[i] << "\n";
                    
                    const double absTol = igma2->in_absolute_feasibility_tol;
                    const double relTol = igma2->in_relative_convergence_tol;
                    bool feas = false;
                    
                    int r = prob->isFeasibleToConstraints( threadNumber, gapmin->sol, true, NULL, absTol, relTol, feas, NULL);
                    
                    if(r != 0)
                    {
                        std::cout << "Erro na avaliacao da restricao\n";
                        MRQ_getchar();
                    }
                    
                    std::cout << "Feas para o problema original: " << feas << "\n";
                    
                    MRQ_fixIntVarsOnSolByList( nI, intVars, gapmin->sol, *nlp );
                    
                    nlp->setInitialSolution( gapmin->sol, gapmin->dualSolC, gapmin->dualSolV );
                    
                    r = nlp->solve();
                    
                    
                    std::cout << "\t\t\tBusca local. retCode: " << nlp->retCode << " feasSol: " << nlp->feasSol << " obj: " << nlp->objValue << "\n";
                    
                    /*for(int i = 0; i < n; i++)
                        std::cout << "sol["<<i<<"]: " << nlp->sol[i] << " \t"; //std::endl;
                    std:: cout << std::endl; */
                }
                #endif
                
                #if 0
                if(indmax < 0)
                {
                    for(int i = 0; i < nI; i++)
                    {
                        const int ind = intVars[i];
                        
                        std::cout << "nlx["<<ind<<"]: " << nlx[ind] << " nux["<<ind<<"]: " << nux[ind] << "\n";
                    }
                }
                #endif
                
                assert(gapmins[threadNumber].solver->retCode == OPT_STOP_REQUIRED_BY_USER || !igma1->in_consider_relax_infeas_if_solver_fail || (!igma1->in_adopt_obj_cut && igma1->in_heuristic_mode) );
                
                assert(indmax > -1 || gapmins[threadNumber].solver->retCode == OPT_STOP_REQUIRED_BY_USER || (!igma1->in_adopt_obj_cut && igma1->in_heuristic_mode));
            #endif
        }
        
    #endif
    
    
    if(indmax >= 0)
    {
        sizeIndices = 1;
        indices[0] = indmax;
        //binary problem
        breakValues1[0] = 0.0; //round( sol[indmax] );
        breakValues2[0] = 1.0;
    }
    
    //std::cout << "Indice escolhido para o branching: " << indmax << std::endl;
    
    return 0;
}



int MRQ_IGMA1BBCallbacks::generateNodes(const int thnumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, BBL_RETURN_CODES retCode, const double objValue, double *sol, double *dualSol, BBL_UserNodeGenerator &userNodeGenerator )
{
    const int m = prob->m;
    const minlpproblem::MIP_SparseMatrix &A = prob->A;
    MRQ_BinSumConstrsChooser &constrChooser = constrChoosers[thnumber];
    
    int index;
    int returnCode;
    BBL_NodeBoundsSol nodeBounds;
    BBL_NodeBoundsSol *zeroNodeBounds = nullptr;
    
    //plc and puc pointers should point to same array was passed in binSumConstrs.calculateIndices
    double *plc = oplc ? oplc : prob->lc;
    double *puc = oplc ? &oplc[m] : prob->uc; 
    
    
    
    
    //branching over a constraint
    
    constrChooser.chooseIndexToBranch( A, igma1->in_integer_tol, sol, constrBranchStrat, NULL, 0.0, NULL, index );
    
    
    #if MRQ_DEBUG_MODE
        assert( index >= 0 );
    #endif
    
    nodeBounds.l = 1.0;
    nodeBounds.u = 1.0;
    
    const int *acols = A[index];
    const unsigned int maxvars = A.getNumberOfElementsAtRow(index);
    
    
    if( plc[index] != puc[index] )
    {
        //in this case, we need generate a node fixing all nofixed variables at zero in this constraints
        int nbounds = 0;
        
        
        MRQ_malloc(zeroNodeBounds, maxvars); 
        MRQ_IFMEMERRORGOTOLABEL( !zeroNodeBounds, returnCode, termination );
        
        for(unsigned int i = 0; i < maxvars; i++)
        {
            const int col = acols[i];
            
            if( nlx[col] == nux[col] )
                continue; //variable is already fixed
            
            zeroNodeBounds[nbounds].ind = col;
            zeroNodeBounds[nbounds].l = 0.0;
            zeroNodeBounds[nbounds].u = 0.0;
            zeroNodeBounds[nbounds].sol = sol[col];
            
            nbounds++;
        }
        
        //we let userNodeGenerator allocate the new node...
        
        int r = userNodeGenerator.generateNode(nbounds, zeroNodeBounds, true, true);
        
        MRQ_IFERRORGOTOLABEL(r, returnCode, MRQ_MEMORY_ERROR, termination);
    }
    
    
    
    
    for(unsigned int i = 0; i < maxvars ; i++)
    {
        const int col = acols[i];
        
        if( nlx[col] == nux[col] )
        {
            //variable is already fixed
            continue;
        }
        
        
        nodeBounds.ind = col;
        nodeBounds.sol = sol[col];
        
        //node will be allocated by generateNode
        int r = userNodeGenerator.generateNode(1, &nodeBounds, true, true);
        MRQ_IFERRORGOTOLABEL(r, returnCode, MRQ_MEMORY_ERROR, termination);
        
    }
    
    returnCode = 0;
    
termination:

    if(zeroNodeBounds)  free(zeroNodeBounds);
    
    return returnCode;
}




int MRQ_IGMA1BBCallbacks::endOfIteration(const int threadNumber, const long unsigned int iter, const double cpuTime, const double wallTime, const double lb, const double ub, BBL_Node &node, const double *nlx, const double *nux)
{
    printOpenNodesList();
    
    MRQ_getchar();
    //assert(false);
    
    return 0;
}


int MRQ_IGMA1BBCallbacks::updatingBestSolution( const int threadNumber, double* sol, double &objValue, const double ub, const long unsigned int iter)
{
    assert(false);
    return 0;
}




void MRQ_IGMA1BBCallbacks::updateSumGapObj( const double *sol, const double *nlx, const double *nux  )
{
    
    SEMAPH_sumGapObj.lock( nthreads );
    
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            if( nlx[ind] != nux[ind] )
            {
                //we assume integer variables are binary
                sumGapObj[i] += MRQ_min( sol[ind], 1.0 -sol[ind] );
                nGapObj[i]++;
            }
        }
    
    SEMAPH_sumGapObj.unlock( nthreads );
}

