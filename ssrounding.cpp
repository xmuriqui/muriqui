/*
 * This file implements the structured stochastic rounding heurist. An experimental feasibility heuristic from our heads
 * 
 * 02, September, 2020
 * 
 * Author: Wendel Melo
 */ 

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <iostream>
#include <new>

#include "MRQ_algClasses.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_tools.hpp"

#include "MRQ_ssrounding.hpp"
#include "MRQ_rens.hpp"


using namespace muriqui;
using namespace optsolvers;
using namespace minlpproblem;


#define MRQ_TRY_FIX_FLOW_CONSTRAINTS_AT_1 0



int MRQ_SSRCore::fixFirstVarsAt0( const int nVarsToFix, const int sizea, const int *acols, double *lx, double *ux )
{
    int nfixed = 0;
    
    if( nVarsToFix > 0 )
    {
        for(int j = 0; j < sizea; j++)
        {
            const int ind = acols[j];
            
            if( lx[ind] != ux[ind] )
            {
                #if MRQ_DEBUG_MODE
                    assert( lx[ind] == 0.0 );
                    assert( ux[ind] == 1.0 );
                #endif
                ux[ind] = 0.0;
                
                nfixed++;
                
                if( nfixed >= nVarsToFix )
                    break;
            }
        }
        
    }
    
    return nfixed;
}




int MRQ_SSRCore::fixFirstVarsAt1( const int nVarsToFix, const int sizea, const int *acols, double *lx, double *ux, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *uc, const MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars )
{
    int nfixed = 0;
    
    if( nVarsToFix > 0 )
    {
        for(int j = 0; j < sizea; j++)
        {
            const int ind = acols[j];
            
            if( lx[ind] != ux[ind] )
            {
                #if MRQ_DEBUG_MODE
                    assert( lx[ind] == 0.0 );
                    assert( ux[ind] == 1.0 );
                #endif
                //lx[ind] = 1.0;
                    
                fixVarAt1(ind, lx, ux, updateVarBoundsByKnapsackConstrs, prob, uc, binSumConstrInds, reverseIntVars);
                
                nfixed++;
                
                if( nfixed >= nVarsToFix )
                    break;
            }
        }
        
    }
    
    return nfixed;
}


//Returns number of variables fixed
int MRQ_SSRCore::fixRandomVarsAt0(MRQ_Random &random, const int nVarsToFix, const int sizea, const int *acols, double *lx, double *ux, const int *reverseIntVars)
{
    return fixRandomVarsAt01(random, nVarsToFix, sizea, acols, lx, ux, true, false, nullptr, nullptr, nullptr, reverseIntVars);
}


//Returns number of variables fixed
int MRQ_SSRCore::fixRandomVarsAt1(MRQ_Random &random, const int nVarsToFix, const int sizea, const int *acols, double *lx, double *ux, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *uc, const MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars )
{
    return fixRandomVarsAt01(random, nVarsToFix, sizea, acols, lx, ux, false, updateVarBoundsByKnapsackConstrs, prob, uc, binSumConstrInds, reverseIntVars );
}


/*if fixAt0 == true,variables will be fixed on 0. Otherwise, variables will be fixed on 1. 
 * Returns number of variables fixed
*/
int MRQ_SSRCore::fixRandomVarsAt01(MRQ_Random &random, const int nVarsToFix, const int sizea, const int *acols, double *lx, double *ux, const bool fixAt0, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *uc, const MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars )
{
    const int maxRandonsFails = 1000 + nVarsToFix;
    int nRandomFails = 0;
    int nFixed = 0;
    
    while( nFixed < nVarsToFix )
    {
        const int randValue = random.randInt(0, sizea-1);
        const int col = acols[ randValue ];
        
        if( lx[col] != ux[col] )
        {
            //we assume binary vars in this constraints
            #if MRQ_DEBUG_MODE
                //assert( avalues[randValue] == 1.0 );
                assert( lx[col] == 0.0 && ux[col] == 1.0 );
            #endif
                
            if( fixAt0 )
                ux[col] = 0.0;
            else
                fixVarAt1(col, lx, ux, updateVarBoundsByKnapsackConstrs, prob, uc, binSumConstrInds, reverseIntVars);
            
            nFixed++;
        }
        else 
        {
            nRandomFails++;
            if( nRandomFails >= maxRandonsFails )
            {
                //so, we give up the random indices, and set by order
                if( fixAt0 )
                    nFixed += fixFirstVarsAt0(nVarsToFix-nFixed, sizea, acols, lx, ux);
                else
                    nFixed += fixFirstVarsAt1(nVarsToFix-nFixed, sizea, acols, lx, ux, updateVarBoundsByKnapsackConstrs, prob, uc, binSumConstrInds, reverseIntVars);
                
                #if MRQ_DEBUG_MODE
                    assert( nFixed <= nVarsToFix );
                #endif
                
                break;
            }
        }
    }
    
    return nFixed;
}


int MRQ_SSRCore:: fixAllNonFixedVarsByStochRounding(const int n, const int nI, const int *intVars, const double *startSol_, const double minimumProbToRound, const int print_level, MRQ_Random &random, int *auxInds, int &nVarsFixed, double *lx, double *ux, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager, MRQ_BoundsUpdaterSolver *boundsUpdaterSolver, double *currlc, double *curruc, const int contRelaxStrategy, optsolvers::OPT_LPSolver *nlpSolver)
{
    nVarsFixed = 0;
    const int *myIntVars;
    const bool preproc_quad = true;
    
    bool solveContRelax = contRelaxStrategy != MRQ_SSR_CRSSR_NO_CONTINUOUS_RELAXATION;
    bool solverBoundsInitialized = false;
    int retCode;
    
    
    if( preprocessor )
    {
        //if we preprocess, maybe is a good idea to run integer variables in a random order
        
        MRQ_copyArray(nI, intVars, auxInds);
        MRQ_shuffleArray(nI, auxInds, random);
        
        myIntVars = auxInds;
    }
    else
    {
        myIntVars = intVars;
    }
    
    
    
    for(int i = 0; i < nI; i++)
    {
        int ind = myIntVars[i];
        
        //std::cout << "i: " << i << " ind: " << ind << "\n";
        
        if(lx[ind] != ux[ind])
        {
            double startSolInd = startSol_[ind];
            
            if( boundsUpdaterSolver )
            {
                bool infeas;
                
                int r = boundsUpdaterSolver->calculateNewVarBounds(n, 1, &ind, lx, ux, true, infeas );
                MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
                
                /*printf("lx[%d]: %0.2f  ux[%d]: %0.2f\n", ind, lx[ind], ind, ux[ind]);
                MRQ_getchar(); */
                
                if(infeas)
                {
                    if(print_level > 5)
                        std::cout << MRQ_PREPRINT "Bounds updater sover detected infeasbility in stochastic rounding \n";
                    return MRQ_HEURISTIC_FAIL;
                }
                
                if( lx[ind] == ux[ind] )
                    continue;
                
            }
            
            
            if( solveContRelax )
            {
                
                //if( !solverBoundsInitialized )
                {
                    //int r = nlpSolver->setnVariablesBounds(n, lx, ux);
                    for(int j = 0; j < nI; j++)
                    {
                        const int jind = intVars[j];
                        int r = nlpSolver->setVariableBounds( jind, lx[jind], ux[jind] );
                        #if MRQ_DEBUG_MODE
                            MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_NLP_SOLVER_ERROR, termination);
                        #endif
                    }
                    
                    solverBoundsInitialized = true;
                }
                
                nlpSolver->solve(false, true, false, false);
                
                
                
                /* std::cout << "i: " << i << " nI: " << nI << " Resolvi relaxação continua para arredondamento. feas: " << nlpSolver->feasSol << " obj: " << nlpSolver->objValue << " minimumProbToRound: " << minimumProbToRound << "\n" ;
                
                for(int j = 0; j < nI; j++)
                    std::cout << "x_" << intVars[j] << ": " << nlpSolver->sol[ intVars[j] ] << "\t"; */
                    
                    
                if( nlpSolver->retCode == OPT_INFEASIBLE_PROBLEM )
                {
                    //MRQ_getchar();
                    retCode = MRQ_HEURISTIC_FAIL;
                    goto termination;
                }
                
                
                if( nlpSolver->feasSol )
                {
                    const double *sol = nlpSolver->sol;
                    const double integerTol = 1e-4;
                    const bool intSol = MRQ_isIntegerSol( nI, intVars, sol, integerTol);
                    
                    if( intSol )
                    { //we found a feasible solution
                        
                        
                        for(int j = 0; j < nI; j++)
                        {
                            const int jind = intVars[j];
                            ux[jind] = lx[jind] = round( sol[jind] );
                        }
                        
                        retCode = 0; 
                        goto termination;
                    }
                    
                    startSolInd = sol[ind];
                }
                else
                {
                    // threating the general integer case: we guarantee some value having 0.5 of integer gap.
                    startSolInd = (lx[ind] + ux[ind])/2.0 ;
                    startSolInd = ((int) startSolInd) + 0.5;
                }
                
                if( contRelaxStrategy == MRQ_SSR_CRSSR_ONLY_BEFORE_STARTING )
                {
                    solveContRelax = false;
                    if( nlpSolver->feasSol  )
                        startSol_ = nlpSolver->sol ;
                }
            }
            
            
            
            double prob =  startSolInd - floor(startSolInd);
            
            //we must let a minimum probability to value be rounded to up and to down
            if( prob < minimumProbToRound )
                prob = minimumProbToRound;
            else if( prob > 1 - minimumProbToRound ) 
                prob = 1 - minimumProbToRound;
            
            if( random.random() < prob )
                lx[ind] = ux[ind] = ceil( startSolInd );
            else
                lx[ind] = ux[ind] = floor( startSolInd );
            
            
            //std::cout<< "\ni: " << i <<" var: " << ind  << " prob: " << prob << " fixed: " << lx[ind] << "\n";
            
            
            if( solveContRelax )
            {
                int r = nlpSolver->setVariableBounds( ind, lx[ind], ux[ind] );
                MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_NLP_SOLVER_ERROR, termination);
            }
            
            
            nVarsFixed++; //this count can be wrong because preprocessor can fix vars also
            
            
            if( preprocessor )
            {
                bool updtVarBounds, updtConstrBounds;
                int r;
                
                //r = preprocessor->preprocess(false, false, INFINITY, lx, ux, updtVarBounds, updtConstrBounds);
                r = preprocessor->preprocess(1, &ind, *ccstorager, preproc_quad, false, INFINITY, lx, ux, updtVarBounds, updtConstrBounds, currlc, curruc, currlc, curruc );

                //MRQ_getchar();
                
                if(r == MIP_INFEASIBILITY)
                {
                    //we got a failure

                    //we try to round to the other side:

                    if( ux[ind] == ceil( startSolInd ) )
                        lx[ind] = ux[ind] = floor( startSolInd );
                    else
                        lx[ind] = ux[ind] = ceil( startSolInd );
                    
                    r = preprocessor->preprocess(1, &ind, *ccstorager, preproc_quad, false, INFINITY, lx, ux, updtVarBounds, updtConstrBounds, currlc, curruc, currlc, curruc );

                    if( r == MIP_INFEASIBILITY )
                    {
                        if(print_level > 5)
                            std::cout << MRQ_PREPRINT "Preprocessor detected infeasbility in stochastic rounding \n";
                        retCode = MRQ_HEURISTIC_FAIL;
                        goto termination;
                    }
                }
            }
            
        }
    }
    
    
    retCode = 0;
    
termination:
    
    return retCode;
}


void MRQ_SSRCore::fixAllNonFixedAt0(const int sizea, const int *acols, double *lx, double *ux)
{
    for(int j = 0; j < sizea; j++)
    {
        const int ind = acols[j];
        
        if( lx[ind] != ux[ind] )
        {
            #if MRQ_DEBUG_MODE
                assert( lx[ind] == 0.0 );
                assert( ux[ind] == 1.0 );
            #endif
            ux[ind] = 0.0;
        }
    }
}



//new version where we let some variables free in <= constraints
int MRQ_SSRCore::class0_3_4Handle(const int sizea, const int *acols, const double *avalues, const double lc, const double uc, MRQ_Random &random, int *auxInds, double *lx, double *ux,  const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *ubc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager )
{
    int *nonFixed = auxInds;
    int nToFixAt1, nToFixAt0 = 0, nNonFixed = 0;
    double rhs = uc < MIP_INFINITY ? uc : lc;
    
    
    #if MRQ_DEBUG_MODE
        assert(sizea > 0);
    #endif
    
    for(int j = 0; j < sizea; j++)
    {
        const int col = acols[j];
        const double coef = avalues[j];
        
        if( lx[col] == ux[col] )
        {
            if( ux[col] != 0.0 )
                rhs -= ux[col] * coef;
        }
        else
        {
            nonFixed[nNonFixed] = col;
            nNonFixed++;
        }
    }
    
    #if MRQ_DEBUG_MODE
        assert( MRQ_isIntegerDouble(rhs) );
    #endif
    
    if( rhs < 0.0 && lc <= -MIP_INFINITY )
        return MRQ_BAD_PARAMETER_VALUES; //we cannot satisfy this constraint
    
    if(lc == uc)
    {
        nToFixAt1 = rhs;
        if(nToFixAt1 == 0)
            nToFixAt0 = nNonFixed; //if we have nToFixAt1 > 0, we will fix the remanider variables to 0 after fix variables to 1.
    }
    else if(lc <= -MIP_INFINITY)
    { // we have a less than contsraint ... <= b.
        //nToFix = random.randInt(0, nToFix);
        
        #if MRQ_DEBUG_MODE
            assert(rhs >= 0.0);
        #endif
        
        //now, we let rhs variables free to be set by other constraints. We must fix all remainder variables to 0. So, we do not fix variables to 1, but we allow until rhs variables can be fix in the future letting those variable free.
        
        nToFixAt1 = 0;
        nToFixAt0 = MRQ_max<int>( nNonFixed - rhs, 0 );        
    }
    else
    {
        #if MRQ_DEBUG_MODE
            assert( uc >= MIP_INFINITY );
        #endif
        // we have a greater than constraint ... >= b.
        // here, we do not fix variables to 0
        nToFixAt1 = MRQ_max<int>(rhs, 0);
    }
    
    
    if(nToFixAt1 > 0)
    {
        if( nToFixAt1 > nNonFixed )
        {
            //we cannot satisfy this constraint
            return MRQ_BAD_PARAMETER_VALUES;
        }
        else if( nToFixAt1 == nNonFixed )
        {
            int r = fixFirstVarsAt1(nToFixAt1, sizea, acols, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
            if(r < nToFixAt1)
                return MRQ_BAD_PARAMETER_VALUES;
        }
        else
        {
            int r = fixRandomVarsAt1(random, nToFixAt1, nNonFixed, nonFixed, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
            if(r < nToFixAt1)
                return MRQ_BAD_PARAMETER_VALUES;
            
            //fixing the remainder variables to 0, but we do noc fix class 4, because if constraint is >= b, it is already satisfied.
            if( uc < MIP_INFINITY )
                fixAllNonFixedAt0(sizea, acols, lx, ux);
        }
    }
    else
    {
        
        if( nToFixAt0 > 0 )
        {
            if( nToFixAt0 == nNonFixed )
            {
                fixFirstVarsAt0( nToFixAt0, sizea, acols, lx, ux); //I think we do not have to test returned value here. 
            }
            else
            { //we assume nToFixAt0 < nNonFixed
                fixRandomVarsAt0(random, nToFixAt0, sizea, acols, lx, ux, reverseIntVars);
            }
        }
    }
    
    return 0;
}


//new version where we let some variables free in <= constraints
int MRQ_SSRCore::class1_2_5_6Handle( const int sizea, const int *acols, const double *avalues, const double lc, const double uc, MRQ_Random &random, int *auxInds1, int *auxInds2, double *lx, double *ux, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *ubc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager )
{
    int *pCoefInds = auxInds1, *negCoefInds = auxInds2;
    int nPCoefInds = 0, nNegCoefInds = 0;
    int nFixPCoefs, nFixNegCoefs;
    int nVarsFixed1 = 0, nVarsNegCoef = 0, nVarsPosCoef = 0;
    double rhs = uc;
    
    
    for(int j = 0; j < sizea; j++)
    {
        const int ind = acols[j];
        const double val = avalues[j];
        
        if( lx[ind] == ux[ind] )
        {
            if( lx[ind] != 0.0 )
            {
                rhs -= lx[ind] * val ;
                nVarsFixed1++;
            }
        }
        else
        {
            if( val == 1.0 )
            {
                pCoefInds[ nPCoefInds ] = ind;
                nPCoefInds++;
            }
            else if( val == -1.0 )
            {
                negCoefInds[ nNegCoefInds ] = ind;
                nNegCoefInds++;
            }
            #if MRQ_DEBUG_MODE
            else
                assert(false); //coeficient should be 1 or -1
            #endif
        }
        
        #if MRQ_TRY_FIX_FLOW_CONSTRAINTS_AT_1
        if( val > 0.0 )
            nVarsPosCoef++;
        else if ( val < 0.0 )
            nVarsNegCoef++;
        #endif
            
    }
    
    /*std::cout << "nPCoefInds: " << nPCoefInds << " nNegCoefInds: " << nNegCoefInds << " rhs: " << rhs << "\n";
    
    for(unsigned int w = 0; w < MRQ_max(nPCoefInds, nNegCoefInds) ; w++)
    {
        if(w < nPCoefInds)
            std::cout << "\tpCoefInds["<<w<<"]:" << pCoefInds[w];
        if(w < nNegCoefInds)
            std::cout << "\t\tnegCoefInds["<<w<<"]:" << negCoefInds[w];
        std::cout << "\n";
    }*/
    
    #if MRQ_DEBUG_MODE
        assert(MRQ_isIntegerDouble(rhs));
    #endif
    
    
    if( rhs == 0.0)
    {
        int maxFix = MRQ_min( nNegCoefInds, nPCoefInds );
        
        nFixNegCoefs = random.randInt(0, maxFix);
        
        #if MRQ_TRY_FIX_FLOW_CONSTRAINTS_AT_1
        if( nVarsPosCoef > 1 && nVarsNegCoef > 1 )
        { /*we assume we have a flow constraint: 
            x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} - x_{n_2} - x_{n_3} - ... - x_{n_q} {=, <=} b */
            
            if( nVarsFixed1 > 0 )
            {
                /*since rhs is zero and nVarsFixed1 > 0, we assume constraints is already satisfied          */
                nFixNegCoefs = 0;
            }
            else
            {
                if(maxFix > 0)
                    nFixNegCoefs = 1; //since we do not have variables fixed at 1, and that is a flow constraint, we try enforce the flow as 1, since the most pat of flux constraints is just to pass 1 by the flow.
            }
        }
        #endif
        
        
        if( lc == uc )
        {//classes 1 and 5
            //if( maxFix > 0) nFixNegCoefs = 1; //TODO: THIS IS A TEST REMOVE THAT
            nFixPCoefs = nFixNegCoefs;
        }
        else // we have a less than contsraint ... <= b. So, we can have more variables with coefficient -1 fixed at 1. We define a number possible lower to fix in nFixPCoefs
        {
            #if MRQ_DEBUG_MODE
                assert(lc <= -MIP_INFINITY);
            #endif
            //now, we let nFixNegCoefs free to be fixed by other contsraints...
            nFixPCoefs = nFixNegCoefs;    
            //nFixPCoefs = random.randInt(0, nFixNegCoefs);
        }
    }
    else if( rhs > 0.0 )
    { //we have more variables having coeficient -1 fixed at 1. So, we need to fix rhs more variables having coefficient 1  at 1 in equalities constraints.
        
        int maxFix = MRQ_min( nNegCoefInds, nPCoefInds - (int) rhs );
        
        if( maxFix < 0 )
        {
            //rhs is positive and we cannot fix any variable here. If constraint is <=, it is ok, constraint is  guaranteed be satisfied. 
            if( lc <= -MIP_INFINITY )
                return 0;
            else
            {
                #if MRQ_DEBUG_MODE
                    assert( lc == uc );
                #endif
                return MRQ_BAD_PARAMETER_VALUES; //equality constraint. In this case, we cannot satisfy this constraint with current set of variables fiexd
            }
        }
        
        nFixNegCoefs = random.randInt(0, maxFix);
        
        
        #if MRQ_TRY_FIX_FLOW_CONSTRAINTS_AT_1
        if( rhs == 1.0 )
        {
            if(nVarsPosCoef > 1 && nVarsNegCoef > 1)
            { /*we assume we have a flow constraint: 
            x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} - x_{n_2} - x_{n_3} - ... - x_{n_q} {=, <=} b */
                
                //we do not let more variables having negative coefs being fixed at 1, since the most part of flow constraints has 1 as flow
                nFixNegCoefs = 0;
            }
        }
        #endif
        
        
        
        if( lc == uc )
        {
            nFixPCoefs = nFixNegCoefs + rhs;
        }
        else
        {
            // we have a less than contsraint ... <= b. So, we can have more variables with coefficient -1 fixed at 1. We define a number possible lower to fix in nFixPCoefs
            #if MRQ_DEBUG_MODE
                assert(lc <= -MIP_INFINITY);
            #endif
            //nFixPCoefs = random.randInt(0, nFixNegCoefs + rhs);
            //now, we let nFixNegCoefs free to be fixed by other contsraints...
            nFixPCoefs = nFixNegCoefs + rhs;
        }
    }
    else //rhs < 0
    { //we have more variables having coeficient 1 fixed at 1. So, we need to fix rhs more variables having coefficient -1  at 1.
        int maxFix = MRQ_min( nNegCoefInds + (int) rhs, nPCoefInds ); //remember: here, rhs is negative, so we have to sum instead of subtrate, since -rhs has the number of variables having positive coefs above the numver of variables having negative coefs.
        
        if( maxFix < 0 )
            return MRQ_BAD_PARAMETER_VALUES; //we cannot satisfy this constraint with currentset of variables fiexd
        
        nFixPCoefs = random.randInt(0, maxFix);
        
        #if MRQ_TRY_FIX_FLOW_CONSTRAINTS_AT_1
        if( rhs == -1.0 )
        {
            if(nVarsPosCoef > 1 && nVarsNegCoef > 1)
            { /*we assume we have a flow constraint: 
            x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} - x_{n_2} - x_{n_3} - ... - x_{n_q} {=, <=} b */
                 //we do not let more variables having positive coefs being fixed at 1, since the most part of flow constraints has 1 as flow
                 nFixPCoefs = 0;
            }
        }
        #endif
        
        
        if( lc == uc )
        {
            nFixNegCoefs = nFixPCoefs - rhs; //here, rhs is negative, so we have to sum instead of subtrate
        }
        else
        {
            // we have a less than contsraint ... <= b. So, we can have more variables with coefficient 1 fixed at 1. We define a number possible lower to fix in nFixPCoefs
            #if MRQ_DEBUG_MODE
                assert(lc <= -MIP_INFINITY);
                assert(nFixPCoefs + rhs <= nNegCoefInds);
            #endif
                
            //nFixNegCoefs = random.randInt(nFixPCoefs - rhs, nNegCoefInds); //remember: here, rhs is negative
            
            //now, we let nFixNegCoefs free to be fixed by other constraints...
            nFixNegCoefs = nFixPCoefs - rhs; //remember: here, rhs is negative
        }
        
    }
    
    
    
    
    #if MRQ_DEBUG_MODE
        assert( nFixPCoefs >= 0 );
        assert( nFixNegCoefs >= 0 );
    #endif    
    
    if( nFixPCoefs == nPCoefInds )
    {
        if( lc == uc ) //equality constraints
        {
            int r = fixFirstVarsAt1(nFixPCoefs, nPCoefInds, pCoefInds, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
            if( r < nFixPCoefs)
                return MRQ_BAD_PARAMETER_VALUES;
        }
        else
        { //at less than constraints ( <= ), we dot not fix variables having positive  coefficients at 1
            #if MRQ_DEBUG_MODE
                assert( lc <= -MIP_INFINITY );
            #endif
        }
    }
    else
    {
        #if MRQ_DEBUG_MODE
            assert( nFixPCoefs < nPCoefInds );
        #endif
        
        if( lc == uc )
        {
            int r = fixRandomVarsAt1(random, nFixPCoefs, nPCoefInds, pCoefInds, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
            if( r < nFixPCoefs )
                return MRQ_BAD_PARAMETER_VALUES;
            
            fixAllNonFixedAt0(nPCoefInds, pCoefInds, lx, ux);
        }
        else
        { //at less than constraints ( <= ), we dot not fix variables having positive  coefficients at 1. We only fix remainder variable to 0 and let nFixPCoefs variables free to be fixed by other constraints.
            #if MRQ_DEBUG_MODE
                assert( lc <= -MIP_INFINITY );
            #endif
            
            int nFixTo0 = nPCoefInds - nFixPCoefs; //remainder coefficients to be fixed on zero. So, at most, nFixPCoefs variables could be fixed at 1 in the future
            
            int r = fixRandomVarsAt0( random, nFixTo0, nPCoefInds, pCoefInds, lx, ux, reverseIntVars );
            if( r != nFixTo0 )
                return MRQ_UNDEFINED_ERROR;
        }
    }
    
    
    if( nFixNegCoefs == nNegCoefInds )
    {
        int r = fixFirstVarsAt1(nFixNegCoefs, nNegCoefInds, negCoefInds, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
        if( r < nFixNegCoefs )
            return MRQ_BAD_PARAMETER_VALUES;
    }
    else
    {
        #if MRQ_DEBUG_MODE
            assert( nFixNegCoefs < nNegCoefInds );
        #endif
        int r = fixRandomVarsAt1(random, nFixNegCoefs, nNegCoefInds, negCoefInds, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
        if( r < nFixNegCoefs )
            return MRQ_BAD_PARAMETER_VALUES;
        
        //now, for <= constraints, we do not fix variables having negative coefficient to 0. We let them free to possibly be fixed by other constraints' handlres
        if(lc == uc)
            fixAllNonFixedAt0(nNegCoefInds, negCoefInds, lx, ux); 
    }
    
    
    return 0;
}







int MRQ_SSRCore::class8Handle(const int sizea, const int *acols, const double *avalues, const double lc, const double uc, MRQ_Random &random, int *auxInds, double *lx, double *ux,  const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *ubc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars)
{
    /*int *nonFixed = auxInds;
    int nNonFixed = 0; */
    int nRandomFails = 0;
    const int maxRandonsFails = 1000;
    double rhs = uc < MIP_INFINITY ? uc : lc;
    
    #if MRQ_DEBUG_MODE
        assert(sizea > 0);
    #endif
    
    for(int j = 0; j < sizea; j++)
    {
        const int col = acols[j];
        
        if( lx[col] == ux[col] )
        {
            const double coef = avalues[j];
            if( ux[col] != 0.0 )
                rhs -= coef * ux[col];
        }
        /*else
        {
            nonFixed[nNonFixed] = col;
            nNonFixed++;
        }*/
    }
    
    
    #if MRQ_DEBUG_MODE
        assert( MRQ_isIntegerDouble(rhs) );
    #endif
    
    
    /*std::cout << "nNonFixed: " << nNonFixed << " rhs: " << rhs << "\n";
    for( int w = 0; w < nNonFixed; w++ )
        std::cout << "nonFixed["<<w<<"]: " << nonFixed[w] << "\n"; */
    
    
    while(rhs > 0.0)
    {
        const int rind = random.randInt(0, sizea-1); //we must pick indices in the costraint array because we wiill need get the coeficient also
        
        const int rcol = acols[rind];
        
        if( lx[rcol] != ux[rcol] )
        {
            #if MRQ_DEBUG_MODE
                assert(lx[rcol] == 0.0 && ux[rcol] == 1.0);
            #endif
            
            fixVarAt1(rcol, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
            
            rhs -= avalues[rind]; //we are fixing at 1
        }
        else 
        {
            nRandomFails++;
            if( nRandomFails >= maxRandonsFails )
            { //so, we give up the random indices, and set by order
                
                for(int j = 0; j < sizea; j++)
                {
                    const int col = acols[j];
                    
                    if( lx[col] != ux[col] )
                    {
                        const double coef = avalues[j];
                        rhs -= coef;
                        
                        if( rhs <= 0.0)
                            break;
                    }
                }
                
                break;
            }
        }
        
    }
    
    if( rhs <= 0.0 )
        return 0;
    else
        return MRQ_BAD_PARAMETER_VALUES; //we could not satisfy this constraint
}



int MRQ_SSRCore::strucStochRounding(const MRQ_MINLPProb &prob, const double *startSol, const double *lx, const double *ux, const double *lc, const double *uc, const int nI, const int *intVars, const int *reverseIntVars, const bool randomOrderToThreatClasses, const bool randomOrderInsideClasses, const double minimumProbToRound, const int print_level, const minlpproblem::MIP_BinSumConstrsIndsByClass &binSumConstrInds, MRQ_Random &random, int *auxInds1, int *auxInds2, unsigned int *auxConstrInds, double *auxConstrs1, double *auxConstrs2, double *outlx, double *outux, int &nVarsFixedByStochRounding, MRQ_SSR_VARIABLES_ADDITIONAL_FIXING_STRATEGY in_vars_additional_fixing_strategy, MRQ_SSR_PREPROCESSING_STRATEGY preprocStrategy, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager, const bool preprocessAfterVarRounding, const int contRelaxStrategy, optsolvers::OPT_LPSolver &nlpSolver, MRQ_BoundsUpdaterSolver *boundsUpdaterSolver)
{
    const int m = prob.m;
    const int n = prob.n;
    int r, retCode;
    
    
    int *myAuxInds1 = nullptr, *myAuxInds2 = nullptr;
    double *myAuxConstrs1 = nullptr, *myAuxConstrs2;
    double *auxlc, *auxuc; //we use auxlc and auxuc to get results from preprocessor and store the redundant constraints...
    
    const int nConstrClassOrder = 8;
    //Test that: THE BALANCE FLOW CONSTRAINTS SHOULD BE THE LAST!
    int constrClassOrder[nConstrClassOrder] = {0, 3, 4, 8, 1, 2, 5, 6}; //order to threat classes constraints 
    
    const MIP_SparseMatrix &A = prob.A;
    const bool updtByKnapacks = true;
    const bool preproc_quad = true;
    
    nVarsFixedByStochRounding = 0;
    
    if(lc == NULL)
        lc = prob.lc;
    
    if(uc == NULL)
        uc = prob.uc;
    
    
    if(!auxInds1)
    {
        MRQ_malloc(myAuxInds1, n);
        MRQ_IFMEMERRORGOTOLABEL(!myAuxInds1, retCode, termination);
        
        auxInds1 = myAuxInds1;
    }
    
    if(!auxInds2)
    {
        MRQ_malloc(myAuxInds2, n);
        MRQ_IFMEMERRORGOTOLABEL(!myAuxInds2, retCode, termination);
        
        auxInds2 = myAuxInds2;
    }
    
    
    if(!auxConstrs1)
    {
        MRQ_malloc(myAuxConstrs1, 2*m);
        MRQ_IFMEMERRORGOTOLABEL(!myAuxConstrs1, retCode, termination);
        
        myAuxConstrs2 = &myAuxConstrs1[m];
        
        auxConstrs1 = myAuxConstrs1;
        auxConstrs2 = myAuxConstrs2;
    }
    
    auxlc = auxConstrs1;
    auxuc = auxConstrs2;
    
    MRQ_copyArray(m, lc, auxlc);
    MRQ_copyArray(m, uc, auxuc);
    
    
    
    MRQ_copyArray(n, lx, outlx);
    MRQ_copyArray(n, ux, outux);
    
    
    
    if(randomOrderToThreatClasses)
    {
        MRQ_shuffleArray(nConstrClassOrder, constrClassOrder, random );
    }
    
    
    for(int j = 0; j < nConstrClassOrder; j++ )
    {
        const int classNumber = constrClassOrder[j];
        
        
        //treating constraints class 0
        const unsigned int nClassc = binSumConstrInds.nClasses[ classNumber ];
        const unsigned *classcInds = binSumConstrInds.classes[ classNumber ];
        
        if(randomOrderInsideClasses && nClassc > 1)
        {
            //we put indices in a random order
            MRQ_copyArray( nClassc, classcInds, auxConstrInds );
            MRQ_shuffleArray(nClassc, auxConstrInds, random);
            classcInds = auxConstrInds;
        }
        
        
        for(unsigned int k = 0; k < nClassc; k++)
        {
            auto cind = classcInds[k];
            
            if( lc[cind] <= -MIP_INFINITY && uc[cind] >= MIP_INFINITY )
                continue; //this constraint is disabled
            
            const int* acols = A[cind];
            const double* avalues = A(cind);
            const unsigned int nel = A.getNumberOfElementsAtRow(cind); 
            
            
            r = classXHandle(classNumber, nel, acols, avalues, lc[cind], uc[cind], random, auxInds1, auxInds2, outlx, outux, updtByKnapacks, &prob, uc, &binSumConstrInds, reverseIntVars, preprocessor, ccstorager);
            
            /*std::cout << "j: " << j << " k: " << k << " class: " << classNumber << " constraint: " << cind;
            std::cout << " return: " << r << "\n";
            for(unsigned int i = 0; i < n; i++)
            {
                std::cout << "outlx["<<i<<"]: " << outlx[i] << " outux["<<i<<"]: " << outux[i] << " \t";
                
                if( i%2 == 1 )
                    std::cout << "\n";
            }
            //MRQ_getchar(); */
            
            if(r != 0)
            {
                //we got a failure
                if(print_level > 10)
                {
                    std::cout << MRQ_PREPRINT "Failure to satisfy constraint " << cind << "\n";
                }

                retCode = MRQ_HEURISTIC_FAIL;
                goto termination;
            }
            
            
            auxlc[cind] = -MIP_INFINITY; //we do not need preprocess more this constraint
            auxuc[cind] = MIP_INFINITY; //we do not need preprocess more this constraint
            
            if(preprocStrategy == MRQ_SSRPS_AFTER_EACH_CONSTRAINT)
            {
                /*std::cout << "########################################\nAntes de preprocessar\n";
                for(unsigned int i = 0; i < n; i++)
                {
                    std::cout << "outlx["<<i<<"]: " << outlx[i] << " outux["<<i<<"]: " << outux[i] << " \t";
                    
                    if( i%2 == 1 )
                        std::cout << "\n";
                } */
                
                #if 1
                if(in_vars_additional_fixing_strategy == MRQ_SSR_VAFS_PREPROCESSING )
                {
                    bool updtVarBounds, updtConstrBounds;
                    int r;
                    
                    r = preprocessor->preprocess(nel, acols, *ccstorager, preproc_quad, false, INFINITY, outlx, outux, updtVarBounds, updtConstrBounds, auxlc, auxuc, auxlc, auxuc);
                    
                    if(r == MIP_INFEASIBILITY)
                    {
                        //we got a failure
                        if(print_level > 10)
                            std::cout << MRQ_PREPRINT "Preprocessor detect infeasbility after constraint " << cind << " in SSR heuristic\n";
                        
                        //MRQ_getchar();
                        retCode = MRQ_HEURISTIC_FAIL;
                        goto termination;
                    }
                    
                    
                    /*if(updtVarBounds)
                    {
                        std::cout << MRQ_PREPRINT "Bounds updated by ssr preprocssing after handling constraint " << cind << "\n";
                        
                        for(unsigned int i = 0; i < n; i++)
                        {
                            std::cout << "outlx["<<i<<"]: " << outlx[i] << " outux["<<i<<"]: " << outux[i] << " \t";
                            
                            if( i%2 == 1 )
                                std::cout << "\n";
                        }
                        MRQ_getchar();
                    }*/
                }
                #endif


                
                #if 0
                else if( in_vars_additional_fixing_strategy == MRQ_SSRVAFS_LINEAR_AUXILIAR_PROBLEM || in_vars_additional_fixing_strategy == MRQ_SSRVAFS_AUXILIAR_PROBLEM )
                {
                    bool infeas;
                    
                    r = boundsUpdaterSolver->calculateNewVarBounds(n, nI, intVars, outlx, outux, true, infeas );
                    MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_NLP_SOLVER_ERROR, termination);
                    
                    if( infeas )
                    {
                        //we got a failure
                        if(print_level > 5)
                            std::cout << MRQ_PREPRINT "Bounds updater detected infeasbility after constraint class " << classNumber << "\n";
                        
                        retCode = MRQ_HEURISTIC_FAIL;
                        goto termination;
                    }
                }
                
                #endif
            
                /*std::cout << "########################################\nApós preprocessar\n";
                for(unsigned int i = 0; i < n; i++)
                {
                    std::cout << "outlx["<<i<<"]: " << outlx[i] << " outux["<<i<<"]: " << outux[i] << " \t";
                    
                    if( i%2 == 1 )
                        std::cout << "\n";
                } 
                
                MRQ_getchar();*/
            }
            
        }
        
        
        if(nClassc > 0 && preprocStrategy == MRQ_SSRPS_AFTER_EACH_CONSTRAINT_CLASS )
        {
            if( in_vars_additional_fixing_strategy == MRQ_SSR_VAFS_PREPROCESSING )
            {
                bool updtVarBounds, updtConstrBounds;
                int r;
                
                //r = preprocessor->preprocess(false, false, INFINITY, outlx, outux, updtVarBounds, updtConstrBounds);  
                r = preprocessor->preprocess(0, nullptr, *ccstorager, preproc_quad, false, INFINITY, outlx, outux, updtVarBounds, updtConstrBounds, auxlc, auxuc, auxlc, auxuc); //here, we preprocess for all variables
                if(r == MIP_INFEASIBILITY)
                {
                    //we got a failure
                    if(print_level > 5)
                        std::cout << MRQ_PREPRINT "Preprocessor detected infeasbility after constraint class " << classNumber << "\n";
                    retCode = MRQ_HEURISTIC_FAIL;
                    goto termination;
                }
            }
            
            #if 0
            
            else if( in_vars_additional_fixing_strategy == MRQ_SSRVAFS_LINEAR_AUXILIAR_PROBLEM || in_vars_additional_fixing_strategy == MRQ_SSRVAFS_AUXILIAR_PROBLEM )
            {
                bool infeas;
                
                r = boundsUpdaterSolver->calculateNewVarBounds(n, nI, intVars, outlx, outux, true, infeas );
                MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_NLP_SOLVER_ERROR, termination);
                
                if( infeas )
                {
                    //we got a failure
                    if(print_level > 5)
                        std::cout << MRQ_PREPRINT "Bounds updater detected infeasbility after constraint class " << classNumber << "\n";
                    
                    retCode = MRQ_HEURISTIC_FAIL;
                    goto termination;
                }
            }
            
            #endif
            
            
        }
    }
    
    
    {
        
        r = fixAllNonFixedVarsByStochRounding(n, nI, intVars, startSol, minimumProbToRound, print_level, random, auxInds1, nVarsFixedByStochRounding, outlx, outux, preprocessAfterVarRounding ? preprocessor : nullptr, ccstorager, boundsUpdaterSolver, auxlc, auxuc, contRelaxStrategy, &nlpSolver);
        
        /*std::cout << "after stochastic rounding: \n";
        for(unsigned int i = 0; i < n; i++)
        {
            std::cout << "outlx["<<i<<"]: " << outlx[i] << " outux["<<i<<"]: " << outux[i] << "\n";
        }
        MRQ_getchar(); */
        
        if( r != 0 )
        {
            if( r != MRQ_HEURISTIC_FAIL )
            {
                MRQ_PRINTERRORNUMBER(r);
            }
            
            retCode = r;
            goto termination;
        }
    }
    
    
    retCode = 0;
    
termination:
    
    
    if(myAuxConstrs1)       free(myAuxConstrs1);
    if(myAuxInds1)    free(myAuxInds1);
    if(myAuxInds2)    free(myAuxInds2);
    
    return retCode;
}



/* this function returns number of variable fixed in the variable nVarsFixed. This function tries fix variables in the array candidateVars aplying preprocessing after each fixing. If after fixing some variable, preprocessing
fix any additional variables, this function abort the fixing returning the number of variabless fixed in the array candidateVars
*/


static inline int MRQ_tryFixRandomVariablesWithPreprocessWithouAditionalBoundsFixing(int nCandidateVars, int *candidateVars, int nVarsToFix, double valueToFix, double *lx, double *ux, MRQ_Random &random, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager, int &nVarsFixed, bool &additionalVarFixed)
{
    bool updtVarBounds = false;
    bool updtConstrBounds;
    int r, randInd, col;
    
    nVarsFixed = 0;
    additionalVarFixed = false;

    
    while( nVarsFixed < nVarsToFix  && !updtVarBounds )
    {
        long unsigned int ntries = 0;

        do 
        {
            if( ntries > 100lu*nCandidateVars )
            {
                int w;
                //we tried a lot get a random variable to fix... so, we run in the order
                for(w = 0; w < nCandidateVars; w++)
                {
                    col = candidateVars[w];
                    if( lx[col] != ux[col] )
                        break;
                }

                if( w == nCandidateVars )
                {
                    std::cout << "nCandidateVars: " << nCandidateVars << "\n";
                    for(w = 0; w < nCandidateVars; w++)
                        std::cout << "var: " << candidateVars[w] << " lx: " << lx[candidateVars[w]] << " ux: " << ux[candidateVars[w]] << "\n" ;

                    //if we reah this pointsomethig is wrong. we do not have variables to fix. Some counting is wrong, or maybe preprocessing is not flaging a variable bound changing
                    assert(false);
                }
            }

            randInd = random.randInt(0, nCandidateVars-1);
            col = candidateVars[randInd];
            ntries++;
        }while( lx[col] == ux[col] ); 

        #if MRQ_DEBUG_MODE
            assert( lx[col] != ux[col] );
        #endif

        lx[col] = ux[col] = valueToFix;
        nVarsFixed++;

        r = preprocessor->preprocess(1, &col, *ccstorager, true, false, INFINITY, lx, ux, updtVarBounds, updtConstrBounds);
        if(r != 0)
        {
            if(r !=MIP_INFEASIBILITY )
            {
                MRQ_PRINTERRORNUMBER(r);
            }
            else
            {
                //we try fix variable in the other side
                if( lx[col] == 0.0 )
                    lx[col] = ux[col] = 1.0;
                else if( lx[col] == 1.0 )
                    lx[col] = ux[col] = 0.0;
                
                additionalVarFixed = true;
                return 0;
            }
            
            return MRQ_HEURISTIC_FAIL;
        }

        

        additionalVarFixed = additionalVarFixed || updtVarBounds;
        
    }

    return 0;
}




/* CLASS 0:  x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} = b
*  CLASS 3:  x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} <= b
*  CLASS 4:  x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} >= b
* 
* Here, all variables have coefficient +1 in the constraints 
*/
int MRQ_SSRCore2::class0_3_4Handle(const int sizea, const int *acols, const double *avalues, const double lc, const double uc, MRQ_Random &random, int *auxInds, double *lx, double *ux,  const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *ubc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager )
{
    int *nonFixed = auxInds;
    
    #if MRQ_DEBUG_MODE
        assert(sizea > 0);
    #endif
    
    int nToFixAt1, nToFixAt0, nNonFixed;
    int nfixed0, nfixed1;


    do
    {
        double rhs = uc < MIP_INFINITY ? uc : lc;
        nToFixAt0 = 0;
        nToFixAt1 = 0;
        nNonFixed = 0;
        nfixed0 = 0;
        nfixed1 = 0;

        for(int j = 0; j < sizea; j++)
        {
            const int col = acols[j];
            const double coef = avalues[j];
            
            if( lx[col] == ux[col] )
            {
                if( ux[col] != 0.0 )
                    rhs -= ux[col] * coef;
            }
            else
            {
                nonFixed[nNonFixed] = col;
                nNonFixed++;
            }
        }
        
        #if MRQ_DEBUG_MODE
            assert( MRQ_isIntegerDouble(rhs) );
        #endif
        
        if( rhs < 0.0 && lc <= -MIP_INFINITY )
            return MRQ_BAD_PARAMETER_VALUES; //we cannot satisfy this constraint
        
        if(lc == uc)
        {
            nToFixAt1 = rhs;
            if(nToFixAt1 == 0)
                nToFixAt0 = nNonFixed; //if we have nToFixAt1 > 0, we will fix the remanider variables to 0 after fix variables to 1.
        }
        else if(lc <= -MIP_INFINITY)
        { // we have a less than contsraint ... <= b.
            //nToFix = random.randInt(0, nToFix);
            
            #if MRQ_DEBUG_MODE
                assert(rhs >= 0.0);
            #endif
            
            //now, we let rhs variables free to be set by other constraints. We must fix all remainder variables to 0. So, we do not fix variables to 1, but we allow until rhs variables can be fix in the future letting those variable free.
            
            nToFixAt1 = 0;
            nToFixAt0 = MRQ_max<int>( nNonFixed - rhs, 0 );        
        }
        else
        {
            #if MRQ_DEBUG_MODE
                assert( uc >= MIP_INFINITY );
            #endif
            // we have a greater than constraint ... >= b.
            // here, we do not fix variables to 0
            nToFixAt1 = MRQ_max<int>(rhs, 0);
        }
        
        
        if(nToFixAt1 > 0)
        {
            if( nToFixAt1 > nNonFixed )
            {
                //we cannot satisfy this constraint
                return MRQ_BAD_PARAMETER_VALUES;
            }
            /*else if( nToFixAt1 == nNonFixed )
            {
                int r = fixFirstVarsAt1(nToFixAt1, sizea, acols, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
                if(r < nToFixAt1)
                    return MRQ_BAD_PARAMETER_VALUES;
            }
            else */
            {
                bool additionalVarFixed;

                //fixing variables to 1
                int r = MRQ_tryFixRandomVariablesWithPreprocessWithouAditionalBoundsFixing(nNonFixed, nonFixed, nToFixAt1, 1.0, lx, ux, random, preprocessor, ccstorager, nfixed1, additionalVarFixed );
                if( r != 0 )
                {
                    if(r != MRQ_HEURISTIC_FAIL)
                        MRQ_PRINTERRORNUMBER(r);
                    return r;
                }


                if( nfixed1 < nToFixAt1 || additionalVarFixed )
                {
                    continue;
                }
                else if( nfixed1 == nToFixAt1 && uc < MIP_INFINITY )
                {
                    //we have to fix all remainder  nonfixed vars to zero. 
                    bool updtVarBounds, updtConstrBounds;
                    
                    fixAllNonFixedAt0(nNonFixed, nonFixed, lx, ux);  //fixing the remainder variables to 0, but we do noc fix class 4, because if constraint is >= b, it is already satisfied.

                    int r = preprocessor->preprocess(nNonFixed, nonFixed, *ccstorager, true, false, INFINITY, lx, ux, updtVarBounds, updtConstrBounds);

                    if( r != 0 )
                    {
                        if( r != MIP_INFEASIBILITY )
                        {
                            MRQ_PRINTERRORNUMBER(r);
                            return MRQ_UNDEFINED_ERROR;
                        }

                        return MRQ_HEURISTIC_FAIL;
                    }
                }

                /*int r = fixRandomVarsAt1(random, nToFixAt1, nNonFixed, nonFixed, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
                if(r < nToFixAt1)
                    return MRQ_BAD_PARAMETER_VALUES;
                
                //fixing the remainder variables to 0, but we do noc fix class 4, because if constraint is >= b, it is already satisfied.
                if( uc < MIP_INFINITY )
                    fixAllNonFixedAt0(sizea, acols, lx, ux);*/
            }
        }
        
        //else
        {

            if( nToFixAt0 > 0 )
            {
                /*if( nToFixAt0 == nNonFixed )
                {
                    fixFirstVarsAt0( nToFixAt0, sizea, acols, lx, ux); //I think we do not have to test returned value here. 
                }
                else */
                {
                    bool additionalVarFixed;

                    //fixing variables to 0
                    int r = MRQ_tryFixRandomVariablesWithPreprocessWithouAditionalBoundsFixing( nNonFixed, nonFixed, nToFixAt0, 0.0, lx, ux, random, preprocessor, ccstorager, nfixed0, additionalVarFixed );
                    if( r != 0 )
                    {
                        if(r != MRQ_HEURISTIC_FAIL)
                            MRQ_PRINTERRORNUMBER(r);
                        return r;
                    }

                    if( nfixed0 < nToFixAt0 || additionalVarFixed)
                        continue;

                    //fixRandomVarsAt0(random, nToFixAt0, sizea, acols, lx, ux, reverseIntVars);
                }
            }
        }
    
    } while (nfixed0 < nToFixAt0 && nfixed1 < nToFixAt1 );
    
    return 0;
}



/*
* CLASS 1: x_{p_1} + x_{p_2} + ... + x_{p_k} - x_{n_1} = b
* CLASS 2: x_{p_1} + x_{p_2} + ... + x_{p_k} - x_{n_1} <= b
* CLASS 5: x_{p_1} + x_{p_2} + ... + x_{p_k} - x_{n_1} - x_{n_2} - ... - x_{n_q} = b  
* CLASS 6: x_{p_1} + x_{p_2} + ... + x_{p_k} - x_{n_1} - x_{n_2} - ... - x_{n_q} <= b
*/
int MRQ_SSRCore2::class1_2_5_6Handle( const int sizea, const int *acols, const double *avalues, const double lc, const double uc, MRQ_Random &random, int *auxInds1, int *auxInds2, double *lx, double *ux, const bool updateVarBoundsByKnapsackConstrs, const MRQ_MINLPProb *prob, const double *ubc, const minlpproblem::MIP_BinSumConstrsIndsByClass *binSumConstrInds, const int *reverseIntVars, MRQ_Preprocessor *preprocessor, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager )
{
    int *pCoefInds = auxInds1, *negCoefInds = auxInds2;
    int nPCoefInds, nNegCoefInds;
    int nFixPCoefs, nFixNegCoefs;
    int nVarsFixed1, nVarsNegCoef, nVarsPosCoef;
    int nPCoefFixed, nNegCoefFixed;
    double rhs;
    
    
    do
    {
        nPCoefInds = 0;
        nNegCoefInds = 0;
        nVarsFixed1 = 0;
        nVarsNegCoef = 0;
        nVarsPosCoef = 0;

        nPCoefFixed = 0;
        nNegCoefFixed = 0;

        rhs = uc;


        for(int j = 0; j < sizea; j++)
        {
            const int ind = acols[j];
            const double val = avalues[j];
            
            if( lx[ind] == ux[ind] )
            {
                if( lx[ind] != 0.0 )
                {
                    rhs -= lx[ind] * val ;
                    nVarsFixed1++;
                }
            }
            else
            {
                if( val == 1.0 )
                {
                    pCoefInds[ nPCoefInds ] = ind;
                    nPCoefInds++;
                }
                else if( val == -1.0 )
                {
                    negCoefInds[ nNegCoefInds ] = ind;
                    nNegCoefInds++;
                }
                #if MRQ_DEBUG_MODE
                else
                    assert(false); //coeficient should be 1 or -1
                #endif
            }
            
            #if MRQ_TRY_FIX_FLOW_CONSTRAINTS_AT_1
            if( val > 0.0 )
                nVarsPosCoef++;
            else if ( val < 0.0 )
                nVarsNegCoef++;
            #endif
                
        }
        
        /*std::cout << "nPCoefInds: " << nPCoefInds << " nNegCoefInds: " << nNegCoefInds << " rhs: " << rhs << "\n";
        
        for(unsigned int w = 0; w < MRQ_max(nPCoefInds, nNegCoefInds) ; w++)
        {
            if(w < nPCoefInds)
                std::cout << "\tpCoefInds["<<w<<"]:" << pCoefInds[w];
            if(w < nNegCoefInds)
                std::cout << "\t\tnegCoefInds["<<w<<"]:" << negCoefInds[w];
            std::cout << "\n";
        }*/
        
        #if MRQ_DEBUG_MODE
            assert(MRQ_isIntegerDouble(rhs));
        #endif
        
        
        if( rhs == 0.0)
        {
            int maxFix = MRQ_min( nNegCoefInds, nPCoefInds );
            
            nFixNegCoefs = random.randInt(0, maxFix);
            
            #if MRQ_TRY_FIX_FLOW_CONSTRAINTS_AT_1
            if( nVarsPosCoef > 1 && nVarsNegCoef > 1 )
            { /*we assume we have a flow constraint: 
                x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} - x_{n_2} - x_{n_3} - ... - x_{n_q} {=, <=} b */
                
                if( nVarsFixed1 > 0 )
                {
                    /*since rhs is zero and nVarsFixed1 > 0, we assume constraints is already satisfied          */
                    nFixNegCoefs = 0;
                }
                else
                {
                    if(maxFix > 0)
                        nFixNegCoefs = 1; //since we do not have variables fixed at 1, and that is a flow constraint, we try enforce the flow as 1, since the most pat of flux constraints is just to pass 1 by the flow.
                }
            }
            #endif
            
            
            if( lc == uc )
            {//classes 1 and 5
                //if( maxFix > 0) nFixNegCoefs = 1; //TODO: THIS IS A TEST REMOVE THAT
                nFixPCoefs = nFixNegCoefs;
            }
            else // we have a less than contsraint ... <= b. So, we can have more variables with coefficient -1 fixed at 1. We define a number possible lower to fix in nFixPCoefs
            {
                #if MRQ_DEBUG_MODE
                    assert(lc <= -MIP_INFINITY);
                #endif
                //now, we let nFixNegCoefs free to be fixed by other contsraints...
                nFixPCoefs = nFixNegCoefs;    
                //nFixPCoefs = random.randInt(0, nFixNegCoefs);
            }
        }
        else if( rhs > 0.0 )
        { //we have more variables having coeficient -1 fixed at 1. So, we need to fix rhs more variables having coefficient 1  at 1 in equalities constraints.
            
            int maxFix = MRQ_min( nNegCoefInds, nPCoefInds - (int) rhs );
            
            if( maxFix < 0 )
            {
                //rhs is positive and we cannot fix any variable here. If constraint is <=, it is ok, constraint is  guaranteed be satisfied. 
                if( lc <= -MIP_INFINITY )
                    return 0;
                else
                {
                    #if MRQ_DEBUG_MODE
                        assert( lc == uc );
                    #endif
                    return MRQ_BAD_PARAMETER_VALUES; //equality constraint. In this case, we cannot satisfy this constraint with current set of variables fiexd
                }
            }
            
            nFixNegCoefs = random.randInt(0, maxFix);
            
            
            #if MRQ_TRY_FIX_FLOW_CONSTRAINTS_AT_1
            if( rhs == 1.0 )
            {
                if(nVarsPosCoef > 1 && nVarsNegCoef > 1)
                { /*we assume we have a flow constraint: 
                x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} - x_{n_2} - x_{n_3} - ... - x_{n_q} {=, <=} b */
                    
                    //we do not let more variables having negative coefs being fixed at 1, since the most part of flow constraints has 1 as flow
                    nFixNegCoefs = 0;
                }
            }
            #endif
            
            
            
            if( lc == uc )
            {
                nFixPCoefs = nFixNegCoefs + rhs;
            }
            else
            {
                // we have a less than contsraint ... <= b. So, we can have more variables with coefficient -1 fixed at 1. We define a number possible lower to fix in nFixPCoefs
                #if MRQ_DEBUG_MODE
                    assert(lc <= -MIP_INFINITY);
                #endif
                //nFixPCoefs = random.randInt(0, nFixNegCoefs + rhs);
                //now, we let nFixNegCoefs free to be fixed by other contsraints...
                nFixPCoefs = nFixNegCoefs + rhs;
            }
        }
        else //rhs < 0
        { //we have more variables having coeficient 1 fixed at 1. So, we need to fix rhs more variables having coefficient -1  at 1.
            int maxFix = MRQ_min( nNegCoefInds + (int) rhs, nPCoefInds ); //remember: here, rhs is negative, so we have to sum instead of subtrate, since -rhs has the number of variables having positive coefs above the numver of variables having negative coefs.
            
            if( maxFix < 0 )
                return MRQ_BAD_PARAMETER_VALUES; //we cannot satisfy this constraint with currentset of variables fiexd
            
            nFixPCoefs = random.randInt(0, maxFix);
            
            #if MRQ_TRY_FIX_FLOW_CONSTRAINTS_AT_1
            if( rhs == -1.0 )
            {
                if(nVarsPosCoef > 1 && nVarsNegCoef > 1)
                { /*we assume we have a flow constraint: 
                x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} - x_{n_2} - x_{n_3} - ... - x_{n_q} {=, <=} b */
                    //we do not let more variables having positive coefs being fixed at 1, since the most part of flow constraints has 1 as flow
                    nFixPCoefs = 0;
                }
            }
            #endif
            
            
            if( lc == uc )
            {
                nFixNegCoefs = nFixPCoefs - rhs; //here, rhs is negative, so we have to sum instead of subtrate
            }
            else
            {
                // we have a less than constraint ... <= b. So, we can have more variables with coefficient 1 fixed at 1. We define a number possible lower to fix in nFixPCoefs
                #if MRQ_DEBUG_MODE
                    assert(lc <= -MIP_INFINITY);
                    assert(nFixPCoefs + rhs <= nNegCoefInds);
                #endif
                    
                //nFixNegCoefs = random.randInt(nFixPCoefs - rhs, nNegCoefInds); //remember: here, rhs is negative
                
                //now, we let nFixNegCoefs free to be fixed by other constraints...
                nFixNegCoefs = nFixPCoefs - rhs; //remember: here, rhs is negative
            }
            
        }
        
        
        
        
        #if MRQ_DEBUG_MODE
            assert( nFixPCoefs >= 0 );
            assert( nFixNegCoefs >= 0 );
        #endif    
        
        /*if( nFixPCoefs == nPCoefInds )
        {
            if( lc == uc ) //equality constraints
            {
                int r = fixFirstVarsAt1(nFixPCoefs, nPCoefInds, pCoefInds, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
                if( r < nFixPCoefs)
                    return MRQ_BAD_PARAMETER_VALUES;
            }
            else
            { //at less than constraints ( <= ), we dot not fix variables having positive  coefficients at 1
                #if MRQ_DEBUG_MODE
                    assert( lc <= -MIP_INFINITY );
                #endif
            }
        } 
        else */
        {
            
            if( lc == uc )
            {
                bool additionalVarFixed;
                bool updtVarBounds, updtConstrBounds;

                //fixing variables at 1
                int r = MRQ_tryFixRandomVariablesWithPreprocessWithouAditionalBoundsFixing(nPCoefInds, pCoefInds, nFixPCoefs, 1.0, lx, ux, random, preprocessor, ccstorager, nPCoefFixed, additionalVarFixed );
                if( r != 0 )
                {
                    if(r != MRQ_HEURISTIC_FAIL)
                        MRQ_PRINTERRORNUMBER(r);
                    return r;
                }

                if( nPCoefFixed < nFixPCoefs || additionalVarFixed )
                    continue;

                /*int r = fixRandomVarsAt1(random, nFixPCoefs, nPCoefInds, pCoefInds, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
                if( r < nFixPCoefs )
                    return MRQ_BAD_PARAMETER_VALUES; */
                
                fixAllNonFixedAt0(nPCoefInds, pCoefInds, lx, ux);
                r = preprocessor->preprocess( nPCoefInds, pCoefInds, *ccstorager, true, false, INFINITY, lx, ux, updtVarBounds, updtConstrBounds );

                if( r != 0 )
                {
                    if( r != MIP_INFEASIBILITY )
                    {
                         MRQ_PRINTERRORNUMBER(r);
                        return MRQ_UNDEFINED_ERROR;
                    }

                    return MRQ_HEURISTIC_FAIL;
                }

                if(updtVarBounds)
                    continue;

            }
            else
            { //at less than constraints ( <= ), we dot not fix variables having positive  coefficients at 1. We only fix remainder variable to 0 and let nFixPCoefs variables free to be fixed by other constraints.
                #if MRQ_DEBUG_MODE
                    assert( lc <= -MIP_INFINITY );
                #endif
                
                bool additionalVarFixed;
                int nfixed;
                int nFixTo0 = nPCoefInds - nFixPCoefs; //remainder coefficients to be fixed on zero. So, at most, nFixPCoefs variables could be fixed at 1 in the future

                //fixing variables at 0
                int r = MRQ_tryFixRandomVariablesWithPreprocessWithouAditionalBoundsFixing( nPCoefInds, pCoefInds, nFixTo0, 0.0, lx, ux, random, preprocessor, ccstorager, nfixed, additionalVarFixed);
                if( r != 0 )
                {
                    if(r != MRQ_HEURISTIC_FAIL)
                        MRQ_PRINTERRORNUMBER(r);
                    return r;
                }

                if( nfixed < nFixTo0 || additionalVarFixed)
                    continue;
                
                //if we reach here, we consider we fix variables to 1 to let the do while loop ends...
                nPCoefFixed = nFixPCoefs;


                /*int r = fixRandomVarsAt0( random, nFixTo0, nPCoefInds, pCoefInds, lx, ux, reverseIntVars );
                if( r != nFixTo0 )
                    return MRQ_UNDEFINED_ERROR;*/
            }
        }
        
        
        /*if( nFixNegCoefs == nNegCoefInds )
        {
            int r = fixFirstVarsAt1(nFixNegCoefs, nNegCoefInds, negCoefInds, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
            if( r < nFixNegCoefs )
                return MRQ_BAD_PARAMETER_VALUES;
        }
        else */
        {
            //fixing vars to 1
            bool additionalVarFixed;

            int r = MRQ_tryFixRandomVariablesWithPreprocessWithouAditionalBoundsFixing( nNegCoefInds, negCoefInds, nFixNegCoefs, 1.0, lx, ux, random, preprocessor, ccstorager, nNegCoefFixed, additionalVarFixed );
            if( r != 0 )
            {
                if(r != MRQ_HEURISTIC_FAIL)
                    MRQ_PRINTERRORNUMBER(r);
                return r;
            }
            

            if( nNegCoefFixed < nFixNegCoefs || additionalVarFixed )
                continue;

            /*int r = fixRandomVarsAt1(random, nFixNegCoefs, nNegCoefInds, negCoefInds, lx, ux, updateVarBoundsByKnapsackConstrs, prob, ubc, binSumConstrInds, reverseIntVars);
            if( r < nFixNegCoefs )
                return MRQ_BAD_PARAMETER_VALUES;*/
            
            //now, for <= constraints, we do not fix variables having negative coefficient to 0. We let them free to possibly be fixed by other constraints' handlres
            if(lc == uc)
            {
                bool updtVarBounds, updtConstrBounds;

                fixAllNonFixedAt0(nNegCoefInds, negCoefInds, lx, ux);
                r = preprocessor->preprocess( nNegCoefInds, negCoefInds, *ccstorager, true, false, INFINITY, lx, ux, updtVarBounds, updtConstrBounds );

                if( r != 0 )
                {
                    if( r != MIP_INFEASIBILITY )
                    {
                         MRQ_PRINTERRORNUMBER(r);
                        return MRQ_UNDEFINED_ERROR;
                    }

                    return MRQ_HEURISTIC_FAIL;
                }
            }
        }

    }while( nPCoefFixed < nFixPCoefs || nNegCoefFixed < nFixNegCoefs );
    
    
    return 0;
}












MRQ_SSRoundingExecutor::MRQ_SSRoundingExecutor()
{
    resetParameters();
    resetOutput();
    
    
    auxConstrInds = nullptr;
    
    auxVarInds1 = nullptr;
    auxVarInds2 = nullptr;
    
    auxVarValues1 = nullptr;
    roundedLx = nullptr;
    roundedUx = nullptr;
    
    auxConstrs1 = nullptr;
    auxConstrs2 = nullptr;
    
    boundsUpdaterSolver = nullptr;
    subProb = nullptr;
}


MRQ_SSRoundingExecutor::~MRQ_SSRoundingExecutor()
{
    deallocate();
}


void MRQ_SSRoundingExecutor::resetOutput()
{
    //out_alg = nullptr;
    out_vars_fixed_by_stoch_rounding = false;
    out_number_of_main_iterations = 0;
    out_number_of_improvments = 0;
    out_local_search_alg_code = MRQ_UNDEFINED_ALG;
    out_cpu_time_at_nlp_integer_fixed_sol = NAN;
    out_obj_at_nlp_integer_fixed_sol = NAN;
}


void MRQ_SSRoundingExecutor::resetParameters()
{
    //in_preprocess_after_handling_constraints = false;
    in_preprocess_after_variable_rounding = true;
    in_random_order_to_threat_classes = false;
    in_random_order_to_threat_constraints_in_each_class = true;
    in_solve_minlp_as_local_search = true;
    in_stop_local_search_solving_on_first_improvment_solution = false;
    
    in_max_number_of_main_iterations = 1000;
    in_max_number_of_improvments = 10;
    
    in_print_level = 4;
    
    in_milp_solver = MRQ_getDefaultMILPSolverCode();
    in_nlp_solver = MRQ_getDefaultNLPSolverCode();
    
    in_neighborhood_strategy = MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD;
    in_preprocessing_point_strategy = MRQ_SSRPS_AFTER_EACH_CONSTRAINT;
    in_additional_vars_fixing_strategy = MRQ_SSR_VAFS_PREPROCESSING;
    
    in_rounding_var_bounds_updt_strategy = MRQ_SSR_VBUS_NO_UPDATING;
    in_cont_relax_strategy_to_stoch_rounding = MRQ_SSR_CRSSR_ONLY_BEFORE_STARTING;
    
    in_absolute_feasibility_tol = 1e-3;
    in_relative_feasibility_tol = 1e-6;
    
    in_integer_tol = 1e-3;
    in_integer_neighborhood_factor = 0.3;
    in_continuous_neighborhood_factor = 0.1;
    in_relative_convergence_tol_to_local_search = 0.2;
    in_max_cpu_time = INFINITY;
    in_max_time = INFINITY;
    in_min_probability_to_round = 0.5;
    
    in_local_search_algorithm = MRQ_UNDEFINED_ALG;
    
    in_milp_solver_params = nullptr;
    in_nlp_solver_params = nullptr;
    in_alg_params = nullptr;
    
    in_alg = nullptr;
}


int MRQ_SSRoundingExecutor::allocateAuxMemory(const unsigned int n, const unsigned int m, const unsigned int sizeAuxConstrIndex)
{
    
    MRQ_malloc(auxConstrInds, sizeAuxConstrIndex);
    MRQ_malloc(auxVarInds1, n);
    MRQ_malloc(auxVarInds2, n);
    MRQ_malloc(auxVarValues1, n);
    MRQ_malloc(roundedLx, n);
    MRQ_malloc(roundedUx, n);
    MRQ_malloc(auxConstrs1, m);
    MRQ_malloc(auxConstrs2, m);
    
    MRQ_IFMEMERRORRETURN( !auxConstrInds || !auxVarInds1 || !auxVarInds2 || !auxVarValues1 || !roundedLx || !roundedUx || !auxConstrs1 || !auxConstrs2 );
    
    return 0;
}


void MRQ_SSRoundingExecutor::deallocate()
{
    MRQ_secFree(auxConstrInds);
    MRQ_secFree(auxVarInds1);
    MRQ_secFree(auxVarInds2);
    MRQ_secFree(auxVarValues1);
    MRQ_secFree(roundedLx);
    MRQ_secFree(roundedUx);
    MRQ_secFree(auxConstrs1);
    MRQ_secFree(auxConstrs2);
    
    MRQ_secDelete(subProb);
    MRQ_secDelete(boundsUpdaterSolver);
}


int MRQ_SSRoundingExecutor::run(const MRQ_MINLPProb &prob, const minlpproblem::MIP_BinSumConstrsIndsByClass &binSumConstrInds, const minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager, MRQ_Random &random, optsolvers::OPT_LPSolver &nlpSolver, unsigned int thnumber, unsigned int nThreads, double insideSolverMaxTime, const double *lc, const double *uc, double *lx, double *ux, const int nI, const int *intVars, const int *reverseIntVars, const int nC, const int *contVars, const double *relaxSol, int &algReturnCode, double &outObj, double *outSol)
{
    const double timeStart = MRQ_getTime();
    const clock_t clockStart = clock();
    
    const int n = prob.n;
    const int m = prob.m;
    int r, myRetCode = MRQ_UNDEFINED_ERROR;
    int nVarsFixedByStochRounding;
    unsigned int j;
    MRQ_SSRCore2 core;
    MRQ_Algorithm *pAlg, *myAlg = nullptr;
    MRQ_Preprocessor *preprocessor = nullptr;
    
    
    algReturnCode = MRQ_UNDEFINED_ERROR;
    outObj = INFINITY;
    
    resetOutput();
    
    
    if( in_additional_vars_fixing_strategy == MRQ_SSR_VAFS_PREPROCESSING )
    {
        if(in_preprocessing_point_strategy != MRQ_SSRPS_NO_PREPROCESSING)
        {
            preprocessor = new (std::nothrow) MRQ_Preprocessor(&prob);
            MRQ_IFMEMERRORGOTOLABEL(!preprocessor, myRetCode, termination);
            
            r = preprocessor->allocateMemory(n, m);
            MRQ_IFMEMERRORGOTOLABEL(r, myRetCode, termination);
        }
    }
    
    
    if( in_rounding_var_bounds_updt_strategy == MRQ_SSR_VBUS_LINEAR_AUXILIAR_PROBLEM || in_rounding_var_bounds_updt_strategy == MRQ_SSR_VBUS_AUXILIAR_PROBLEM )
    {
        if( !boundsUpdaterSolver )
        {
            const bool setLinearProblem = in_rounding_var_bounds_updt_strategy == MRQ_SSR_VBUS_LINEAR_AUXILIAR_PROBLEM;
            const int solverCode = in_rounding_var_bounds_updt_strategy == MRQ_SSR_VBUS_LINEAR_AUXILIAR_PROBLEM || prob.getProblemType() == MIP_PT_MILP ? (int) in_milp_solver : (int) in_nlp_solver ;
            
            
            
            
            boundsUpdaterSolver = new (std::nothrow) MRQ_BoundsUpdaterSolver();
            MRQ_IFMEMERRORGOTOLABEL( !boundsUpdaterSolver, myRetCode, termination );
            
            
            //NOTE: by now, we are not recieving solver parameter to set
            r = boundsUpdaterSolver->buildProblem( solverCode, prob, setLinearProblem, false, NAN, thnumber, 1, NULL, NULL, NULL, false );
            MRQ_IFMEMERRORGOTOLABEL(r, myRetCode, termination);
        }
    }
    
    
    if(!auxVarInds1)
    {
        unsigned int maxSizeMClasses = 0;
        
        for(unsigned int k = 0; k < binSumConstrInds.nBinSumConstrClasses; k++)
        {
            if( binSumConstrInds.nClasses[k] > maxSizeMClasses )
                maxSizeMClasses = binSumConstrInds.nClasses[k];
        }
        
        r = allocateAuxMemory(n, m, maxSizeMClasses);
        MRQ_IFERRORGOTOLABEL(r, myRetCode, r, termination);
    }
    
    
    for(out_number_of_main_iterations = 1; out_number_of_main_iterations <= in_max_number_of_main_iterations; out_number_of_main_iterations++)
    {
        
        //printf("iter: %u\n", out_number_of_main_iterations);
        
        int ret = core.strucStochRounding(prob, relaxSol, lx, ux, lc, uc, nI, intVars, reverseIntVars, in_random_order_to_threat_classes, in_random_order_to_threat_constraints_in_each_class, in_min_probability_to_round, in_print_level, binSumConstrInds, random, auxVarInds1, auxVarInds2, auxConstrInds, auxConstrs1, auxConstrs2, roundedLx, roundedUx, nVarsFixedByStochRounding, in_additional_vars_fixing_strategy, in_preprocessing_point_strategy, preprocessor, ccstorager, in_preprocess_after_variable_rounding, in_cont_relax_strategy_to_stoch_rounding, nlpSolver, boundsUpdaterSolver );
        
        if( in_max_cpu_time < INFINITY )
        {
            const double cpuTime = MRQ_calcCPUTtime(clockStart, clock()); 
            
            if( cpuTime >= in_max_cpu_time)
            {
                myRetCode = MRQ_HEURISTIC_FAIL;
                break;
            }
        }
        
        if( in_max_time < INFINITY )
        {
            const double wallTime = MRQ_getTime() - timeStart;
            
            if( wallTime >= in_max_time )
            {
                myRetCode = MRQ_HEURISTIC_FAIL; //we return heuristic fail because here we do not have feasible solution
                break;
            }
        }
        
        
        if( ret == MRQ_HEURISTIC_FAIL ) //if we get heuristic failure, we just try again in another iteration. 
            continue;
        MRQ_IFERRORGOTOLABEL(ret, myRetCode, MRQ_HEURISTIC_FAIL, termination);
        
        
        out_vars_fixed_by_stoch_rounding = nVarsFixedByStochRounding > 0;
            
        if( nC == 0)
        {
            bool feasSol = false;
            prob.isFeasibleToConstraints(thnumber, roundedLx, true, nullptr, in_absolute_feasibility_tol, in_relative_feasibility_tol, feasSol);
            
            if( feasSol )
            {
                r = prob.objEval(thnumber, prob.hasNLConstraints(), roundedLx, outObj);
                MRQ_IFERRORGOTOLABEL(r, myRetCode, MRQ_CALLBACK_FUNCTION_ERROR, termination);
                
                out_obj_at_nlp_integer_fixed_sol = outObj;
                MRQ_copyArray(n, (const double *) roundedLx, outSol);
                out_cpu_time_at_nlp_integer_fixed_sol = MRQ_calcCPUTtime(clockStart);
                
                break;
            }
        }
        else
        {
            //MRQ_getchar();
            
            r = MRQ_fixIntVarsOnSolByList(nI, intVars, roundedLx, nlpSolver); 
            MRQ_IFERRORGOTOLABEL(r, myRetCode, r, termination);
            
            //do not store contsraints value and dual solution because we can use ssr inside a branch-and-bound procedure, and so, we would have to backup the values from continuous relaxation
            nlpSolver.solve(false, true, false, false);
            
            r = MRQ_unfixIntegerVarsByList(nI, intVars, lx, ux, nlpSolver);
            
            //std::cout << "nlp local search - ret: " << nlpSolver.retCode << " obj: " << nlpSolver.objValue << "\n";
            
            /*for(unsigned int w = 0; w < n; w++)
                std::cout << "\tsol["<<w<<"]: " << nlpSolver.sol[w] << "\n"; */
            
            //MRQ_getchar();
            
            if( nlpSolver.feasSol )
            {
                MRQ_copyArray( n, nlpSolver.sol, outSol );
                outObj = nlpSolver.objValue;
                out_obj_at_nlp_integer_fixed_sol = nlpSolver.objValue;
                out_cpu_time_at_nlp_integer_fixed_sol = MRQ_calcCPUTtime(clockStart);
                break;
            }
        }
        
    }
    
    
    if( !std::isinf(outObj) && in_solve_minlp_as_local_search )
    {
        const int neighConstrIndex = prob.m;
        double *auxValues = auxVarValues1;
        double *nlx = roundedLx, *nux = roundedUx;
        
        
        if(in_neighborhood_strategy == MRQ_SNS_ORIGINAL)
        {
            MRQ_PRINTERRORMSG( MRQ_STRPARINTVALUE(MRQ_SNS_ORIGINAL) " cannot be used like neighborhood strategy in SSR. Changing to " MRQ_STRPARINTVALUE(MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD) "!");
            
            in_neighborhood_strategy = MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD;
        }
        
        
        if(!subProb)
        {
            const int nnewConstraints = (in_neighborhood_strategy == MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD && nI < n) ? 2 : 1; //to set the continuous euclidean neighborhood, we must have at least one continuous var
            
            subProb = new (std::nothrow) MRQ_MINLPProb;
            MRQ_IFMEMERRORGOTOLABEL(!subProb, myRetCode, termination);
            
            r = subProb->copyProblemFrom(prob);
            MRQ_IFERRORGOTOLABEL(r, myRetCode, MRQ_MEMORY_ERROR, termination);
            
            r = subProb->addConstraints(nnewConstraints);
            MRQ_IFERRORGOTOLABEL(r, myRetCode, MRQ_MEMORY_ERROR, termination);
        }
        
        
        
        if(in_alg)
        {
            pAlg = in_alg;
        }
        else
        {
            
            if(in_local_search_algorithm == MRQ_SSR_HEUR_ALG)
                in_local_search_algorithm = MRQ_UNDEFINED_ALG; //we would have a infinite recursion...
            
            myAlg = MRQ_newAlgorithm(in_local_search_algorithm, subProb->getNumberOfNLEqualityConstraints() );
            MRQ_IFMEMERRORGOTOLABEL(!myAlg, myRetCode, termination);
            
            if(in_alg_params)
                myAlg->setParameters(*in_alg_params);
            pAlg = myAlg;
        }
        
        out_local_search_alg_code = pAlg->out_algorithm;
        
        pAlg->in_number_of_threads = nThreads;
        pAlg->in_print_level = in_print_level - 2;
        pAlg->in_milp_solver = in_milp_solver;
        pAlg->in_nlp_solver = in_nlp_solver;
        
        if(in_print_level > 2)
            printf( MRQ_PREPRINT "best objective: %0.10lf\n", outObj);
        
        
        pAlg->in_relative_convergence_tol = in_relative_convergence_tol_to_local_search;
        
        
        for(j = 0; j < in_max_number_of_improvments; j++)
        {
            r = MRQ_setSubproblemNeighborhood(lx, ux, outSol, nI, intVars, nC, contVars, in_neighborhood_strategy, in_integer_neighborhood_factor, in_continuous_neighborhood_factor, in_integer_tol, neighConstrIndex, auxVarInds1, auxValues, subProb, nlx, nux); 
            MRQ_IFERRORGOTOLABEL(r, myRetCode, r, termination);
            
            if( !std::isinf(in_max_cpu_time) )
                pAlg->in_max_cpu_time = in_max_cpu_time - MRQ_calcCPUTtime(clockStart);
            if( !std::isinf(in_max_time) )
                pAlg->in_max_time = in_max_time - (MRQ_getTime() - timeStart);
            
            pAlg->in_upper_bound = outObj;
            
            if( pAlg->isLinearApproximationAlgorithm() )
            {
                MRQ_LinearApproxAlgorithm *pLA = (MRQ_LinearApproxAlgorithm*) pAlg;
                
                pLA->deletePointsToLinearisation();
                
                r = pLA->addPointsToLinearisation(1, n, &outSol);
                MRQ_IFERRORGOTOLABEL(r, myRetCode, r, termination);
            }
            else
            {
                r = pAlg->setInitialSolution(n, outSol);
                MRQ_IFERRORGOTOLABEL(r, myRetCode, r, termination);
            }
            
            
            if( in_stop_local_search_solving_on_first_improvment_solution )
            {
                //here, we imposing a solution at least 10% better to stop on first improvment
                pAlg->in_lower_bound = outObj - 0.1*MRQ_abs(outObj) - 0.1;
            }
            
            
            /*std::cout << "in_integer_neighborhood_factor: " << in_integer_neighborhood_factor << " n diff: " << (int) ceil(in_integer_neighborhood_factor*n) << "\n";
            for(int w = 0; w < n; w++)
                std:: cout << "sol["<<w<<"]: " << outSol[w] << "\n";
            
            subProb->print();*/
            
            
            MRQ_insideRun(pAlg, *subProb, in_milp_solver_params, in_nlp_solver_params, thnumber, insideSolverMaxTime, nlx, nux);
            
            if(in_print_level > 2)
                printf( MRQ_PREPRINT "Local search subproblem %d - return code: %d  best objective: %0.10lf\n", j+1, pAlg->out_return_code, pAlg->out_best_obj);
            
            if(pAlg->out_feasible_solution && pAlg->out_best_obj < outObj)
            {
                outObj = pAlg->out_best_obj;
                MRQ_copyArray(n, pAlg->out_best_sol, outSol );
                out_number_of_improvments++;
            }
            else
            {
                break;
            }
            
            
            if( !std::isinf(in_max_cpu_time) )
            {
                if( MRQ_calcCPUTtime(clockStart) >= in_max_cpu_time )
                {
                    //algReturnCode = MRQ_MAX_TIME_STOP; //this will be overwritten, but ok...
                    break;
                }
            }
            
            if( !std::isinf(in_max_time) )
            {
                if( MRQ_getTime() - timeStart >= in_max_time )
                    break;
            }
        }
    }
    
    
termination:

    if( std::isinf(outObj) )
        algReturnCode = MRQ_HEURISTIC_FAIL;
    else
        algReturnCode = MRQ_HEURISTIC_SUCCESS;


    if(myAlg)   delete myAlg;
    if(preprocessor)    delete preprocessor;

    return myRetCode;
}





MRQ_StructuredStochasticRounding:: MRQ_StructuredStochasticRounding()
{
    resetParameters();
    resetOutput();
    
    out_algorithm = MRQ_SSR_HEUR_ALG;
}


MRQ_StructuredStochasticRounding:: ~MRQ_StructuredStochasticRounding()
{
    
}


int MRQ_StructuredStochasticRounding:: checkAlgorithmRequirements( MRQ_MINLPProb &prob, const double *lx, const double *ux)
{
    if( !MRQ_isBinarieProblemAtRegion(prob, lx, ux) )
    {
        if( in_print_level > 4 )
            MRQ_PRINTERRORMSG("We are sorry, but Sturctured Stochastic Rounding only can be applyed to Binary Problems.");
        return MRQ_ALG_NOT_APPLICABLE;
    }
    
    return 0;
}


void MRQ_StructuredStochasticRounding:: printParameters( std::ostream &out) const
{
    char strValue[100];
    
    MRQ_Heuristic::printParameters(out);
    out << "\n"
    //MRQ_STRFFATT(in_preprocess_after_handling_constraints) << "\n"
    MRQ_STRFFATT(in_preprocess_after_variable_rounding) << "\n"
    MRQ_STRFFATT(in_random_order_to_threat_classes) << "\n"
    MRQ_STRFFATT(in_random_order_to_threat_constraints_in_each_class) << "\n"
    MRQ_STRFFATT(in_solve_continuous_relaxation) << "\n"
    MRQ_STRFFATT(in_solve_minlp_as_local_search) << "\n"
    MRQ_STRFFATT(in_stop_local_search_solving_on_first_improvment_solution) << "\n"
    
    MRQ_STRFFATT(in_max_number_of_improvments) << "\n"
    
    MRQ_STRFFATT(in_integer_neighborhood_factor) << "\n"
    MRQ_STRFFATT(in_continuous_neighborhood_factor) << "\n"
    MRQ_STRFFATT(in_relative_convergence_tol_to_local_search) << "\n"
    MRQ_STRFFATT(in_min_probability_to_round) << "\n"
    ;
    
    MRQ_enumToStr(in_neighborhood_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_neighborhood_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr(in_preprocessing_point_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_preprocessing_point_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr(in_additional_vars_fixing_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_additional_vars_fixing_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr(in_rounding_var_bounds_updt_strategy, strValue);
    out << MRQ_STRPARINTVALUE(in_rounding_var_bounds_updt_strategy) " " << strValue << "\n";
    
    MRQ_enumToStr(in_cont_relax_strategy_to_stoch_rounding, strValue);
    out << MRQ_STRPARINTVALUE(in_cont_relax_strategy_to_stoch_rounding) " " << strValue << "\n";
    
    MRQ_enumToStr(in_local_search_algorithm, strValue);
    out << MRQ_STRPARINTVALUE(in_local_search_algorithm) " " << strValue << "\n";
}


void MRQ_StructuredStochasticRounding:: resetParameters()
{
    MRQ_Heuristic::resetParameters();
    
    //in_preprocess_after_handling_constraints = false;
    in_preprocess_after_variable_rounding = true;
    in_random_order_to_threat_classes = false;
    in_random_order_to_threat_constraints_in_each_class = true;
    in_solve_continuous_relaxation = true;
    in_solve_minlp_as_local_search = true;
    in_stop_local_search_solving_on_first_improvment_solution = false;
    
    in_max_iterations = 1000;
    in_max_number_of_improvments = 10;
    
    in_print_level = 4;
    
    in_milp_solver = MRQ_getDefaultMILPSolverCode();
    in_nlp_solver = MRQ_getDefaultNLPSolverCode();
    
    in_neighborhood_strategy = MRQ_SNS_LOCAL_BRANCHING_NEIGHBORHOOD;
    in_preprocessing_point_strategy = MRQ_SSRPS_AFTER_EACH_CONSTRAINT;
    in_additional_vars_fixing_strategy = MRQ_SSR_VAFS_PREPROCESSING;
    in_rounding_var_bounds_updt_strategy = MRQ_SSR_VBUS_NO_UPDATING;
    in_cont_relax_strategy_to_stoch_rounding = MRQ_SSR_CRSSR_ONLY_BEFORE_STARTING;
    
    in_integer_tol = 1e-3;
    in_integer_neighborhood_factor = 0.3;
    in_continuous_neighborhood_factor = 0.1;
    in_relative_convergence_tol_to_local_search = 0.2;
    in_min_probability_to_round = 0.1;
    in_max_cpu_time = INFINITY;
    in_max_time = INFINITY;
    in_local_search_algorithm = MRQ_UNDEFINED_ALG;
    
    
    in_alg_params = nullptr;
    
    in_alg = nullptr;
}


void MRQ_StructuredStochasticRounding:: resetOutput()
{
    out_vars_fixed_by_stoch_rounding = false;
    out_cpu_time_at_nlp_integer_fixed_sol = NAN;
    out_obj_at_nlp_integer_fixed_sol = NAN;
    out_obj_at_first_sol = NAN;
    out_local_search_alg_code = MRQ_UNDEFINED_ALG;
    MRQ_Heuristic::resetOutput();
}


int MRQ_StructuredStochasticRounding:: setIntegerParameter( const char *name, const long int value)
{
    int ret = MRQ_Heuristic::setIntegerParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_preprocess_after_variable_rounding), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_random_order_to_threat_classes), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_random_order_to_threat_constraints_in_each_class), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_solve_continuous_relaxation), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_solve_minlp_as_local_search), name, value ) == 0 );
    else if( MRQ_setAtt<bool>( MRQ_STRATT(in_stop_local_search_solving_on_first_improvment_solution), name, value ) == 0 );
    else if( MRQ_setAtt<unsigned int>( MRQ_STRATT(in_max_number_of_improvments), name, value) == 0 );
    else 
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


int MRQ_StructuredStochasticRounding:: setDoubleParameter( const char *name, const double value)
{
    int ret = MRQ_Algorithm::setDoubleParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt( MRQ_STRATT(in_integer_neighborhood_factor), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_continuous_neighborhood_factor), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_relative_convergence_tol_to_local_search), name, value ) == 0 );
    else if( MRQ_setAtt( MRQ_STRATT(in_min_probability_to_round), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


int MRQ_StructuredStochasticRounding:: setStringParameter( const char *name, const char *value)
{
    int r;
    int ret = MRQ_Algorithm::setStringParameter(name, value);
    
    if( ret == 0 )
        return 0;
    
    ret = 0;
    
    if( (r = MRQ_setStrAtt( MRQ_STRATT(in_neighborhood_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_preprocessing_point_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_additional_vars_fixing_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_rounding_var_bounds_updt_strategy), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_cont_relax_strategy_to_stoch_rounding), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    else if( (r = MRQ_setStrAtt( MRQ_STRATT(in_local_search_algorithm), name, value ) ) >= 0 )
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


int MRQ_StructuredStochasticRounding:: run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams )
{
    return run(prob, milpSolverParams, nlpSolverParams, nullptr);
}


int MRQ_StructuredStochasticRounding:: run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, MRQ_GeneralSolverParams* algParams)
{
    const double timeStart = MRQ_getTime();
    const clock_t clockStart = clock();
    
    const int n = prob.n;
    //const int m = prob.m;
    int nI, nC;
    double outObj;
    
    double *plc = NULL, *puc = NULL; // do not free puc, because we take advantage the malloc of plc. We just put NULL because we should sinalize if there is new bounds or no to other procedures like MRQ_BinSumConstrs::calculateIndices when we do not preprocess...
    
    
    int *intVars = NULL, *reverseIntVars = NULL;
    int *contVars;
    double *constrValues = NULL;
    double *initSol = NULL;
    
    
    int r, algRetCod;
    unsigned int nthreads;
    //double NLPCpuTime, NLPClockTime;
    MIP_BinSumConstrsIndsByClass binSumConstrs;
    MRQ_SSRoundingExecutor ssrExecutor;
    
    MRQ_Random random;

    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux;
    
    OPT_LPSolver *nlp = NULL;
    
    MRQ_Preprocessor *preprocessor = NULL;
    minlpproblem::MIP_ConstraintsByColumnsStorager ccstorager;
    
    
    nthreads = in_number_of_threads <= 0 ? MRQ_getNumCores() : in_number_of_threads;
    
    
    if( in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_obj_function )
    {
        preprocessor = new (std::nothrow) MRQ_Preprocessor(&prob);
        MRQ_IFMEMERRORGOTOLABEL(!preprocessor, out_return_code, termination);
    }
    
    
    {
        auto ret = algorithmInitialization( nthreads, true, milpSolverParams, nlpSolverParams, prob, lx, ux, preprocessor, NULL, &plc, &puc );
        
        if( ret != MRQ_SUCCESS)
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
    }
    
    
    if(in_print_level > 1)
    {
        std::cout << "\n";
        MRQ_PRINTMSG("Starting Structured Stochastic Rounding Heuristic\n\n");
        
        printSubSolvers(true, true, false);
    }
    
    
    
    MRQ_malloc(intVars, n);
    MRQ_malloc(reverseIntVars, n);
    MRQ_malloc(initSol, n);
    
    if( prob.getProblemType() == minlpproblem::MIP_PT_MILP )
    {
        //we try use a linear solver here
        nlp = OPT_newLPSolver(in_milp_solver);
    }
    
     if(!nlp)
        nlp = OPT_newNLPSolver(in_nlp_solver); //we got a failure to instantiate a linear solver. So, we try instantiate a NLP solver.
    
    MRQ_IFMEMERRORGOTOLABEL(!intVars || !reverseIntVars || !initSol || !nlp, out_return_code, termination );
    
    nI = prob.getIntegerIndices(intVars);
    
    contVars = &intVars[nI];
    
    nC = prob.getContinuousIndices(contVars);
    
    #if MRQ_DEBUG_MODE
        assert(n == nI + nC);
    #endif
    
    
    prob.getReverseIntegerIndices(reverseIntVars);
    
    
    r = MRQ_setNLPRelaxProb( prob, lx, ux, plc, puc, nlp, true, true, true, false, thnumber, in_set_special_nlp_solver_params, nlpSolverParams, in_number_of_threads, in_max_cpu_time, in_max_time, 0, 0 );
    MRQ_IFERRORGOTOLABEL(r, out_return_code, (MRQ_RETURN_CODE) r, termination);
    
    
    if( in_solve_continuous_relaxation )
    {
        if(in_print_level > 2)
            std::cout << MRQ_PREPRINT "Solving NLP relaxation\n";
        nlp->solve( false );
        
        if(in_print_level > 4)
        {
            std::cout << MRQ_PREPRINT "Continuous relaxation. ret: " << nlp->retCode << " feasSol: " << nlp->feasSol << "\n";
            if( nlp->feasSol )
            {
                for(int i = 0; i < n; i++)
                    std::cout << "x[" << i << "]: " << nlp->sol[i] << "\n";
            }
        }
        
        
        if( nlp->retCode == OPT_OPTIMAL_SOLUTION || nlp->feasSol )
        {
            const double minGap = 0.1;
            
            if( nlp->retCode == OPT_OPTIMAL_SOLUTION )
            {
                out_obj_opt_at_continuous_relax = nlp->objValue;
                zl = MRQ_max(zl, nlp->objValue);
                
                if( zl > zu )
                {
                    if(in_print_level > 0)
                        MRQ_PRINTMSG("Solution of NLP relaxation is greater than upper_bound ");

                    out_return_code = MRQ_INFEASIBLE_PROBLEM;
                    goto termination;
                }
            }
            
            if( MRQ_isIntegerSol(nI,intVars, nlp->sol, in_integer_tol) )
            {
                out_cpu_time_to_first_feas_sol = MRQ_calcCPUTtime( clockStart );
                out_clock_time_to_first_feas_sol = MRQ_getTime() - timeStart;
                out_obj_at_first_sol = nlp->objValue;
                
                if( in_print_level > 1 )
                    MRQ_PRINTMSG("An integer optimal solution was gotten as NLP relaxation solution\n");
                
                tryUpdateBestSolution( thnumber, n, nlp->sol, nlp->objValue, 0, clockStart, timeStart, in_store_history_solutions );
                
                
                if( nlp->retCode == OPT_OPTIMAL_SOLUTION )
                {
                    out_return_code = MRQ_OPTIMAL_SOLUTION;
                }
                else
                {
                    out_return_code = MRQ_HEURISTIC_SUCCESS;
                }
                goto termination;
                /*std::cout << "UNCOMENT THE GOTO ABOVE!\n";
                MRQ_getchar();*/
            }
            
            MRQ_copyArray(n, nlp->sol, initSol);
            
            
            
            for(int i = 0; i < nI; i++)
            {
                const auto ind = intVars[i];
                
                if( MRQ_gap(initSol[ind]) < minGap )
                {
                    if(initSol[ind] - lx[ind] < minGap)
                    {
                        initSol[ind] += minGap;
                    }
                    else
                    {
                        #if MRQ_DEBUG_MODE
                            assert( ux[ind] - initSol[ind] < minGap );
                        #endif
                        initSol[ind] -= minGap;
                    }
                }
            }
            
        }
        else if( nlp->retCode == OPT_INFEASIBLE_PROBLEM )
        {
            if( in_print_level > 1 )
                MRQ_PRINTMSG("Infeasible NLP relaxation\n");
            
            out_return_code = MRQ_INFEASIBLE_PROBLEM;
            goto termination;
        }
    
    }
    else
    {
        if( in_use_initial_solution && xInit )
        {
            MRQ_copyArray(n, xInit, initSol);
        }
        else
        { //we just put 0.5 as a gap to get the rounding probability
            
            //note, we do not set initial values for continue variables (we do not need this)
            
            for(int k = 0; k < nI; k++)
            {
                const int i = intVars[k];
                initSol[i] = lx[i] == ux[i] ? lx[i] : lx[i] + 0.5;
            }
        }
    }
    
    
    random.setSeed( &in_seed_to_random_numbers );
    
    r = binSumConstrs.calculateIndices(prob, lx, ux, plc, puc, reverseIntVars, true, true, true);
    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_UNDEFINED_ERROR, termination );
    
    r = ccstorager.storageConstraintsByColumns(prob, in_preprocess_quad_constrs);
    MRQ_IFERRORGOTOLABEL(r, out_return_code, MRQ_UNDEFINED_ERROR, termination );
    
    
    //ssrExecutor.in_preprocess_after_handling_constraints =    in_preprocess_after_handling_constraints;
    ssrExecutor.in_preprocess_after_variable_rounding = in_preprocess_after_variable_rounding;
    ssrExecutor.in_random_order_to_threat_classes = in_random_order_to_threat_classes;
    ssrExecutor.in_random_order_to_threat_constraints_in_each_class = in_random_order_to_threat_constraints_in_each_class;
    ssrExecutor.in_solve_minlp_as_local_search = in_solve_minlp_as_local_search;
    ssrExecutor.in_stop_local_search_solving_on_first_improvment_solution = in_stop_local_search_solving_on_first_improvment_solution;
    
    ssrExecutor.in_max_number_of_main_iterations = in_max_iterations;
    ssrExecutor.in_max_number_of_improvments = in_max_number_of_improvments;
    ssrExecutor.in_print_level = in_print_level;
    ssrExecutor.in_nlp_solver = in_nlp_solver;
    ssrExecutor.in_milp_solver = in_milp_solver;
    ssrExecutor.in_neighborhood_strategy = in_neighborhood_strategy;
    ssrExecutor.in_preprocessing_point_strategy = in_preprocessing_point_strategy;
    ssrExecutor.in_additional_vars_fixing_strategy = in_additional_vars_fixing_strategy;
    ssrExecutor.in_rounding_var_bounds_updt_strategy = in_rounding_var_bounds_updt_strategy;
    ssrExecutor.in_cont_relax_strategy_to_stoch_rounding = in_cont_relax_strategy_to_stoch_rounding;
    
    ssrExecutor.in_absolute_feasibility_tol = in_absolute_feasibility_tol;
    ssrExecutor.in_relative_feasibility_tol = in_relative_feasibility_tol;
    
    ssrExecutor.in_integer_tol = in_integer_tol;
    
    ssrExecutor.in_integer_neighborhood_factor = in_integer_neighborhood_factor;
    ssrExecutor.in_continuous_neighborhood_factor = in_continuous_neighborhood_factor;
    ssrExecutor.in_relative_convergence_tol_to_local_search = in_relative_convergence_tol_to_local_search;
    ssrExecutor.in_max_cpu_time = in_max_cpu_time;
    ssrExecutor.in_max_time = in_max_time;
    ssrExecutor.in_min_probability_to_round = in_min_probability_to_round;
    
    ssrExecutor.in_local_search_algorithm = in_local_search_algorithm;
    ssrExecutor.in_milp_solver_params = milpSolverParams;
    ssrExecutor.in_nlp_solver_params = nlpSolverParams;
    ssrExecutor.in_alg_params = algParams;
    ssrExecutor.in_alg = in_alg;
    
    
    r = ssrExecutor.run(prob, binSumConstrs, &ccstorager, random, *nlp, thnumber, nthreads, INFINITY, plc, puc, lx, ux, nI, intVars, reverseIntVars, nC, contVars, initSol , algRetCod, outObj, out_best_sol );
    
    if( !std::isinf(outObj) )
    {
        #if MRQ_DEBUG_MODE
            assert(algRetCod == MRQ_HEURISTIC_SUCCESS);
        #endif
        
        out_best_obj = outObj;
        out_return_code = MRQ_HEURISTIC_SUCCESS;
        out_obj_at_nlp_integer_fixed_sol = ssrExecutor.out_obj_at_nlp_integer_fixed_sol;
        out_obj_at_first_sol = out_obj_at_nlp_integer_fixed_sol;
        out_cpu_time_at_nlp_integer_fixed_sol = ssrExecutor.out_cpu_time_at_nlp_integer_fixed_sol;
        out_cpu_time_to_first_feas_sol = out_cpu_time_at_nlp_integer_fixed_sol;
    }
    else
    {
        #if MRQ_DEBUG_MODE
            assert(algRetCod == MRQ_HEURISTIC_FAIL);
        #endif
        out_return_code = MRQ_HEURISTIC_FAIL;
    }
    
    out_vars_fixed_by_stoch_rounding = ssrExecutor.out_vars_fixed_by_stoch_rounding;
    out_local_search_alg_code = ssrExecutor.out_local_search_alg_code;
    out_number_of_iterations = ssrExecutor.out_number_of_main_iterations;
    
    
termination:
    
    if(plc)			free(plc);
    
    if(intVars)			free(intVars);
    if(reverseIntVars)  free(reverseIntVars);
    if(initSol)     free(initSol);
    
    
    if(constrValues)	free(constrValues);
    
    if(nlp)				delete nlp;
    if(preprocessor)    delete preprocessor;
    
    algorithmFinalization(nthreads, prob, lx, ux);
    
    out_number_of_threads = nthreads;
    out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
    out_clock_time = MRQ_getTime() - timeStart;
    
    if(in_print_level > 1)
        std::cout << MRQ_PREPRINT "cpu time: " << out_cpu_time << "\n";
    
    return out_return_code;
}























