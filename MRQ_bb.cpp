


#include <cmath>
#include <cstdlib>
#include <climits>

#include <new>
#include <iostream>


#include "BBL_branchAndBound.hpp"
#include "BBL_tools.hpp"
#include "MRQ_bb.hpp"
#include "MRQ_tools.hpp"





using namespace muriqui;
using namespace branchAndBound;





void muriqui::MRQ_checkParentNodeBoundsStrategy(const int nI, const int* intVars, const double *lx, const double *ux, branchAndBound::BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY &parentNodeBoundsStrategy)
{
    int i, end;
    long int largestValue, lowestValue;
    
    if( branchAndBound::BBL_PNBSS_SINGLE_MIN <= parentNodeBoundsStrategy && parentNodeBoundsStrategy <= branchAndBound::BBL_PNBSS_SINGLE_MAX )
        end = branchAndBound::BBL_PNBSS_SINGLE_MAX;
    else
    {
        #if MRQ_DEBUG_MODE
            assert( branchAndBound::BBL_PNBSS_ESP_MIN <= parentNodeBoundsStrategy && parentNodeBoundsStrategy <= branchAndBound::BBL_PNBSS_ESP_MAX );
        #endif
        end = branchAndBound::BBL_PNBSS_ESP_MAX;
    }
    
    
    for(i = parentNodeBoundsStrategy; i <= end; i++)
    {
        bool match = true;
        lowestValue = BBL_lowestSequence( (branchAndBound::BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY) i);
        largestValue = BBL_greatestSequence( (branchAndBound::BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY) i );
        
        for(int j = 0; j < nI; j++)
        {
            const int ind = intVars[j];
            
            if( lx[ind] < lowestValue || ux[ind] > largestValue )
            {
                match = false;
                break; //this tipe is appropriate
            }
        }
        
        if(match)
            break;
    }
    
    
    parentNodeBoundsStrategy = (branchAndBound::BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY) i;
    
}







MRQ_NewBBNode::MRQ_NewBBNode( BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY parentBoundsStrategy, const unsigned int maxVars):BBL_Node( BBL_pnbs2usfdnbp(parentBoundsStrategy, maxVars) )
{
    #if MRQ_SET_LBHEUR_ON_BB_NODE
        heurlb = -MRQ_INFINITY;
    #endif
    #if MRQ_DEBUG_IGMA2_BUG
        igma2OnAncestral = 0;
    #endif
    //nBranchVars = 0;
    //nMyBounds = 0;
    //parentBounds = NULL;
    //myBounds = NULL;
}




MRQ_NewBBNode::~MRQ_NewBBNode()
{
    //subDesallocate();
}



double MRQ_NewBBNode::getBestLowerBound() const
{
    double lb = getLowerBound();
    
    #if MRQ_SET_LBHEUR_ON_BB_NODE
        if( heurlb > lb )
            lb = heurlb;
    #endif
    return lb;
}



void MRQ_NewBBNode::print(std::ostream &out) const
{
    //unsigned int i;
    //const unsigned int nDadBounds = parentBounds ? parentBounds->size : 0;
    //BBL_NodeBounds *bounds = parentBounds ? parentBounds->a : 0;
    
    out << "depth: " << getDepth() << " lb: " << getLowerBound() 
    
    #if MRQ_SET_LBHEUR_ON_BB_NODE
        << " heurlb: " << heurlb 
    #endif
    << "\n";
    
    BBL_Node::print(out);
    
    #if 0
    getParentBounds()->print(out);
    
    /*for(i = 0; i < nDadBounds; i++)
        out << "\tvar: " << bounds[i].ind << " l: " << bounds[i].l << " u: " << bounds[i].u << " ";   */
    
    for(i = 0; i < nMyBounds; i++)
        out << "\tvar: " << myBounds[i].ind << " l: " << myBounds[i].l << " u: " << myBounds[i].u << " sol: " << myBounds[i].sol << " ";
    
    out << "\n";
    #endif
}


/*double MRQ_NewBBNode::getLowerBound() const
{
    return MRQ_max( BBL_Node::getLowerBound(), heurlb );
}*/


void MRQ_NewUserNodeGenerator2::init()
{
    parent = NULL;
    n = -1;
}


MRQ_NewUserNodeGenerator2::MRQ_NewUserNodeGenerator2()
{
    init();
}


void MRQ_NewUserNodeGenerator2::initialize( MRQ_NewBBNode *parentNode, BBL_UserNodeGenerator *userNodeGen)
{
    parent = parentNode;
    bblUserNodeGen = userNodeGen;
}


int MRQ_NewUserNodeGenerator2::generateNode( const double *nodelx, const double *nodeux, MRQ_NewBBNode *newNode, const double nodeLowerBound, const double *initSol, const double *initDual )
{
    
    if( newNode == NULL )
    {
        newNode = new (std::nothrow) MRQ_NewBBNode(parentBoundsStrategy, n );
        
        if(!newNode)
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTMEMERROR;
            #endif
            return MRQ_MEMORY_ERROR;
        }
    }
    
    #if MRQ_SET_LBHEUR_ON_BB_NODE
        newNode->heurlb = MRQ_max( nodeLowerBound, parent->heurlb );
    #endif
    
    const int r = bblUserNodeGen->generateNode( nodelx, nodeux, newNode, nodeLowerBound, initSol, initDual );
    
    if( r != 0 )
    {
        if( r == BBL_MEMORY_ERROR )
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTMEMERROR;
            #endif
            return MRQ_MEMORY_ERROR;
        }
        else
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERRORNUMBER(r);
            #endif
            return MRQ_UNDEFINED_ERROR;
        }
    }
    
    return 0;
}




int MRQ_NewUserNodeGenerator2::generateNode( const bool inheritParentBounds, const unsigned int nNewBounds, const MRQ_NodeBoundsSol *newBounds, const bool isNewBoundsAscOrdered, MRQ_NewBBNode *newNode, const double nodeLowerBound, const double *initSol, const double *initDual)
{
    
    if( newNode == NULL )
    {
        newNode = new (std::nothrow) MRQ_NewBBNode(parentBoundsStrategy, n);
        
        if(!newNode)
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTMEMERROR;
            #endif
            return MRQ_MEMORY_ERROR;
        }
    }
    
    
    const int r = bblUserNodeGen->generateNode( nNewBounds, newBounds, inheritParentBounds, isNewBoundsAscOrdered, newNode, nodeLowerBound, initSol, initDual );
    
    if( r != 0 )
    {
        if( r == BBL_MEMORY_ERROR )
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTMEMERROR;
            #endif
            return MRQ_MEMORY_ERROR;
        }
        else
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERRORNUMBER(r);
            #endif
            return MRQ_UNDEFINED_ERROR;
        }
    }
    
    return 0;
}







MRQ_NewFirstBranch::MRQ_NewFirstBranch()
{
    initialize();
}


void MRQ_NewFirstBranch::initialize()
{
    bPoint = nan(""); //Not a number
    leftlb = -INFINITY;
    rightlb= -INFINITY;
}





MRQ_BasePseudoCostCalc::MRQ_BasePseudoCostCalc()
{
    initialize();
}


MRQ_BasePseudoCostCalc::~MRQ_BasePseudoCostCalc()
{
    deallocate();
}

void MRQ_BasePseudoCostCalc::initialize()
{
    pcost = NULL;
}


int MRQ_BasePseudoCostCalc::allocate(const int nI)
{
    pcost = new (std::nothrow) MRQ_NewPseudoCost[nI];
    
    MRQ_IFMEMERRORRETURN(!pcost);
    
    return 0;
}


void MRQ_BasePseudoCostCalc::deallocate()
{
        MRQ_secDeleteArray(pcost);
}


void  MRQ_BasePseudoCostCalc::updatePCosts(const unsigned int nThreads, const int *reverseIntVars, const unsigned int sizeBounds, const branchAndBound::BBL_ClassUnionNodeBoundsSolPointer &bounds, const double fDad, const double fNode, const double* x)
{
    unsigned int nBranchVars = 0;
    const double denTol = 1.0e-8; //tolerance for denominator
    
    if(sizeBounds == 0)
        return; 
    
    for(unsigned int i = 0; i < sizeBounds; i++)
    {
        double isol;
        
        bounds.getArrayElement(i, NULL, NULL, NULL, &isol);
        
        if( !std::isnan(isol) )//if( !std::isnan( bounds[i].sol ) )
            nBranchVars++;
    }
    
    #if MRQ_DEBUG_MODE
        assert(nBranchVars > 0);
    #endif
    
    const double mean = MRQ_max(fNode - fDad, 0.0)/nBranchVars;
    
    
    SEMAPH_sem.lock(nThreads);
    {
        for(unsigned int i = 0; i < sizeBounds; i++)
        {
            unsigned int j; //const unsigned int j = bounds[i].ind;
            double solj; 	//const double solj = bounds[i].sol;
            
            bounds.getArrayElement(i, &j, NULL, NULL, &solj);
            
            
            if( std::isnan(solj) )
                continue;
            
            
            if( MRQ_abs(solj -x[j]) < denTol ) //this second condition is so strange, but, in some exceptional cases, we branch on variables having integer values in the relaxation, for example, when solver fails to find optimal solution in the relation, but finds a feasible one having integer values. In this case, we cannot update the pseudocost because the numerator and denominator are zero (-nan appears)
            {
                #if MRQ_DEBUG_MODE
                    
                    /*if( !(MRQ_gap(solj) == 0.0 || nBranchVars > 1) )
                    {
                        std::cout << "j: " << j << " solj: " << solj << " x_j: " << x[j] <<  " nBranchVars: " << nBranchVars << "\n";
                    }*/
                    
                    assert(MRQ_gap(solj) == 0.0 || nBranchVars > 1); // in this case, we branched on integer value...
                    // we can have branched over a constraint fixing several variables to zero. In this case, one of those variables could have had an integer value in the parent relaxation.
                #endif
                continue;
            }
            
            
            const int ind = reverseIntVars[j];
            
            #if MRQ_DEBUG_MODE
                assert( !std::isnan(pcost[ind].pl) );
                assert( !std::isnan(pcost[ind].pr) );
            #endif
            
            
            if( x[j] < solj )
            { //left branching
                const double pcvalue = mean/(solj - x[j]);
                
                pcost[ind].nPl++; //we first increase the number of pseudocost. do not put this increment after pl updation because we need this in this away to perform pseudo prune without mutex and in this way, we can be sure estimative will not be greather than it should be
                pcost[ind].pl += pcvalue;
                
                if( pcvalue < pcost[ind].minPl  )
                    pcost[ind].minPl = pcvalue;
            }
            else
            { //right branching
                const double pcvalue = mean/(x[j] - solj);
                
                pcost[ind].nPr++; //we first increase the number of pseudocost. do not put this increment after pr updation because we need this in this away to perform pseudo prune without mutex and in this way, we can be sure estimative will not be greather than it should be
                pcost[ind].pr += pcvalue;
                
                if( pcvalue < pcost[ind].minPr )
                    pcost[ind].minPr = pcvalue;
            }
            
            
            #if MRQ_DEBUG_MODE
                
                /*if( isnan( pcost[ind].pr ) )
                {
                    std::cout << "nan em ind: " << ind << " index: " << j << " x: " << x[j] << " solj: " << solj << " mean: "<< mean << "\n";
                    
                    MRQ_getchar();
                } */
                
                
                assert(ind >= 0); // if ind is negative, is because lastBranch[i].ind is not a integer variable...
            #endif
        }
    }
    SEMAPH_sem.unlock(nThreads);
}


//returns a estimative of objective increasing in a single variable. This method does not check the number of psedudo costs calculations, just provide the estimative based on the current value of pseudo costs.
//indexInIntVar is the order (index) that the variable appears in intVars, ie.,  if indexInIntVar has a value i, the variable is the i-th integer variable
int MRQ_BasePseudoCostCalc::getObjIncreaseEstimativeByVar ( const unsigned int indexInIntVar, const unsigned int minNumberOfPCostsCalculations, const double oldVarValue, const double newVarValue,  const double alpha, double &estimative)
{
    const int iind = indexInIntVar;
    double mean, minpc, deltax, step;
    
    estimative = 0.0;
    
        
    if( newVarValue < oldVarValue  )
    { //left branching
        
        if( pcost[iind].nPl <  minNumberOfPCostsCalculations)
            return MRQ_BAD_PARAMETER_VALUES; //this is not an anormal condition. So, we MUST NOT print an error message here
        
        mean = pcost[iind].pl/pcost[iind].nPl;
        minpc = pcost[iind].minPl;
        
        deltax = oldVarValue - newVarValue;
    }
    else
    { //right branching
        
        if( pcost[iind].nPr < minNumberOfPCostsCalculations)
            return MRQ_BAD_PARAMETER_VALUES; //this is not an anormal condition. So, we MUST NOT print an error message here
        
        mean = pcost[iind].pr/pcost[iind].nPr;
        minpc = pcost[iind].minPr;
        
        deltax = newVarValue - oldVarValue;
    }
    
    step = (1.0 - alpha)*minpc  +  alpha*mean;
    estimative =  deltax * step;
    
    
    return 0;
}




int MRQ_BasePseudoCostCalc::getObjIncreaseEstimative( const bool  onlyApplyOnFixedIntVars, const int *reverseIntVars, const unsigned int minNumberOfPCostsCalculations, const unsigned int sizeBounds, const branchAndBound::BBL_ClassUnionNodeBoundsSolPointer &bounds, const double alpha, double &estimative,bool &allVarIncreasingEstimated)
{
    unsigned int nEstimatives = 0;
    estimative = 0.0;
    
    allVarIncreasingEstimated = false;
    
    for(unsigned int i = 0; i < sizeBounds; i++)
    {
        unsigned int vind;
        double vlb, vub, vsol;
        
        bounds.getArrayElement(i, &vind, &vlb, &vub, &vsol );
        
        if( !onlyApplyOnFixedIntVars || vlb == vub  ) //if variable is not fixed, we cannot provide an estimative of its value (actually, we can, maybe we provide an estimative like teh variable were in the lower or upper bound in the future. Do not forget to change checkIfNodeCanBePrunedByEstimative)
        {
            double varEstimative;
            /*const int iind = reverseIntVars[vind];
            double mean, minpc, deltax, step;
            
            if( vub < vsol )
            { //left branching
                
                if( pcost[iind].nPl <  minNumberOfPCostsCalculations)
                    continue;
                
                mean = pcost[iind].pl/pcost[iind].nPl;
                minpc = pcost[iind].minPl;
                
                deltax = (vsol - vub);
            }
            else
            { //right branching
                #if MRQ_DEBUG_MODE
                    assert( vlb > vsol );
                #endif
                
                if( pcost[iind].nPr < minNumberOfPCostsCalculations)
                    continue;
                
                mean = pcost[iind].pr/pcost[iind].nPr;
                minpc = pcost[iind].minPr;
                
                deltax = (vlb - vsol);
            }
            
            step = alpha*mean + (1.0 - alpha)*minpc;
            estimative +=  deltax * step; */
            
            int r = getObjIncreaseEstimativeByVar (reverseIntVars[vind], minNumberOfPCostsCalculations, vsol,   vsol < vlb ? vlb : vub , alpha, varEstimative);
            if(r == 0)
                nEstimatives++;
            
            estimative += varEstimative;
        }
        
    }
    
    allVarIncreasingEstimated = (nEstimatives == sizeBounds);
    
    return 0;
}


//in this method, we check integer variablesnon fixed to see if both left anr dirgth branching will be pseudo pruned. If we find a variable in this conditions, we can assume we can prune this node.
bool MRQ_BasePseudoCostCalc:: checkIfNodeCanBePrunedByEstimative( const bool  onlyApplyOnFixedIntVars,  const int nI, const int *intVars, const unsigned int minNumberOfPCostsCalculations, const double *nlx, const double *nux, const double relaxObj, const double *relaxSol, const double intTol, const double alpha, const double zu )
{
    
    //printf("minNumberOfPCostsCalculations: %d\n", minNumberOfPCostsCalculations);
    
    for(int k = 0; k < nI; k++)
    {
        const int ind = intVars[k];
        const double nlind = nlx[ind], nuind = nux[ind];
        
        #if MRQ_DEBUG_MODE
            /*if( nlind >  relaxSol[ind] )
            {
                printf("nlind: %f nuind: %f  sol: %f\n", nlind, nuind, relaxSol[ind] );
            }*/
            assert( nlind <=  relaxSol[ind] + 1e-4 );
            assert( relaxSol[ind] - 1e-4 <= nuind );
        #endif
        
        
        if( (onlyApplyOnFixedIntVars && nuind - nlind == 1.0) || ( !onlyApplyOnFixedIntVars && nuind - nlind >= 1.0)   ) //if diference between bounds is greather than 1.0, variable can be not fixed in the branhcing. So, we do not try estimative objective increasing (by now, because we could...).
        {
            
            if( MRQ_gap(relaxSol[ind]) > intTol )
            {
                double rsol = relaxSol[ind];
                double leftEstimative = 0.0, rightEstimative = 0.0;
                
                getObjIncreaseEstimativeByVar (k, minNumberOfPCostsCalculations, rsol, floor(rsol), alpha, leftEstimative);
                
                if( relaxObj + leftEstimative >= zu)
                { //this node will be pruned by left branching. So, we check if it will be pruned by right
                    getObjIncreaseEstimativeByVar (k, minNumberOfPCostsCalculations, rsol, ceil(rsol), alpha, rightEstimative);
                    
                    if( relaxObj + rightEstimative >= zu )
                    {
                        //std::cout << "pseudo podando on branching pela variavel " << ind << " node lower bound: " << relaxObj << " lelft estimative: " << relaxObj + leftEstimative << " right estimative: " << relaxObj + rightEstimative << "\n";
                        return true; //this node will be pruned by right branching also. We, we can declare the currente node will be pruned 
                    }
                }
                
            }
        }
        
    }
    
    return false;
    
}


void MRQ_BasePseudoCostCalc::print( const int nI, const int *intVars, std::ostream &out )
{
    for(int i = 0; i < nI; i++ )
    {
        out << i << " var: " << intVars[i] << " -  pl: " << pcost[i].pl <<  "  nrealPl: " << pcost[i].nrealPl << "  nPl: "<< pcost[i].nPl << "  pr: " << pcost[i].pr <<  "  nrealPr: " << pcost[i].nrealPr << "  nPr: " << pcost[i].nPr << "  minPl: " << pcost[i].minPl << "  minPr: " << pcost[i].minPr << "\n";
    }
    
}



MRQ_NewPseudoCostCalc::MRQ_NewPseudoCostCalc():MRQ_BasePseudoCostCalc()
{
    initialize();
}



MRQ_NewPseudoCostCalc::~MRQ_NewPseudoCostCalc()
{
    deallocate();
}


void MRQ_NewPseudoCostCalc::initialize()
{
    MRQ_BasePseudoCostCalc::initialize();
    
    allPCostInit = false;
    maxComptSBranch = 1;
    highestfbranch = -INFINITY;
    fbranch = NULL;
}



int MRQ_NewPseudoCostCalc::allocate(const int nI)
{
    int r = MRQ_BasePseudoCostCalc::allocate(nI);
    MRQ_IFERRORRETURN(r, r);
    
    fbranch = new (std::nothrow) MRQ_NewFirstBranch[nI];
    
    MRQ_IFMEMERRORRETURN(!pcost || !fbranch);
    
    return 0;
}




bool MRQ_NewPseudoCostCalc::checkIfAllPCostsAreCalculated( const int nI, const int *intVars, const double* olx, const double* oux)
{
    bool answer = true;
    int i;
    unsigned int ind;
    
    
    for(i = 0; i < nI; i++)
    {
        ind = intVars[i];
        if( olx[ind] != oux[ind] )
        {
            if( pcost[i].nrealPl < maxComptSBranch || pcost[i].nrealPr < maxComptSBranch )
            {
                answer = false;
                break;
            }
        }
    }
    
    allPCostInit = answer;
    return answer;
}




void MRQ_NewPseudoCostCalc::deallocate()
{
    MRQ_BasePseudoCostCalc::deallocate();
    MRQ_secDeleteArray(fbranch);
}



int muriqui::MRQ_strongBranching( MRQ_NewPseudoCostCalc *pc, const unsigned int totalThreads, const unsigned int nThreads, const unsigned int thnumber, MRQ_Mutex *SEMAPH_indSem, MRQ_Mutex *SEMAPH_solSem, MRQ_Mutex *SEMAPH_nodeSem, int *nextInd, MRQ_MINLPProb *prob,
const int nI, const int *intVars, const bool storeFirstBranch, const bool considInfeasIfNlpFail, const int printLevel, const double intTol, const double *olx, const double *oux, MRQ_NLPSolver *nlp, MRQ_NewBBNode *node, double *nlx, double *nux, const double *nodeSol, const double objNodeSol, const bool calculateEvenIfMaxComptSBranchReached, bool *newBoundsAllocated, bool *updtBounds, bool *canFathom, double *zu, bool *intSolFound, double *intSol, int *retCode )
{
    return pc->threadStrongBranching(totalThreads, nThreads, thnumber, *SEMAPH_indSem, *SEMAPH_solSem, *SEMAPH_nodeSem, *nextInd, *prob, nI, intVars, storeFirstBranch, considInfeasIfNlpFail, printLevel, intTol, olx, oux, nlp, node, nlx, nux, nodeSol, objNodeSol, calculateEvenIfMaxComptSBranchReached, *newBoundsAllocated, *updtBounds, *canFathom, *zu, *intSolFound, intSol, *retCode);
}




int MRQ_NewPseudoCostCalc::threadStrongBranching( const unsigned int totalThreads, const unsigned int nThreads, const unsigned int thnumber, MRQ_Mutex &SEMAPH_indSem, MRQ_Mutex &SEMAPH_solSem, MRQ_Mutex &SEMAPH_nodeSem, int &nextInd, MRQ_MINLPProb &prob, const int nI, const int *intVars, const bool storeFirstBranch, const bool considInfeasIfNlpFail, const int printLevel, const double intTol, const double *olx, const double *oux, MRQ_NLPSolver *nlp, MRQ_NewBBNode *node, double *nlx, double *nux, const double *nodeSol, const double objNodeSol, const bool calculateEvenIfMaxComptSBranchReached, bool &newBoundsAllocated,  bool &updtBounds, bool &canFathom, double &zu, bool &intSolFound, double *intSol, int &retCode )
{
    const int n = prob.n;
    const int plevel = 8;
    const auto maxNMyBounds = node->getMaxNMyBounds();

    bool tryLeft, tryRight;
    unsigned int ind;
    int code = 0, r, myNextInd;


    double bleft, bright;
    double lbind, ubind;

    //BBL_ArraySize<BBL_NodeBounds> *abounds = NULL;



    while( nextInd < nI )
    {
        
        SEMAPH_indSem.lock(nThreads);
        {
            myNextInd = nextInd;
            nextInd++;
        }
        SEMAPH_indSem.unlock(nThreads);
        
        if( myNextInd >= nI)
        {
            break;
        }
        
        
        //std::cerr << "myNextInd: " << myNextInd << " nI: " << nI << " thnumber: " << thnumber << "\n";
        
        ind = intVars[myNextInd];
        
        
        if( nlx[ind] == nux[ind] )
        {
            if( printLevel >= 2*plevel )
                std::cout << MRQ_PREPRINT "Fixed variable " << ind << " at value " << nlx[ind] << "\n";
            
            if( olx[ind] == olx[ind] )
            {
                //in this case the variable is fixed in all nodes. So, we cannot calculate pseudo costs for it. we mark as if we already calculate to avoid that function be called again useless
                
                //I think we do not need worry about semaphores gere, but...
                
                SEMAPH_sem.lock(totalThreads);
                {
                    pcost[myNextInd].nPl = pcost[myNextInd].nrealPl = maxComptSBranch;
                    
                    pcost[myNextInd].nPr = pcost[myNextInd].nrealPr = maxComptSBranch;
                }
                SEMAPH_sem.unlock(totalThreads);
            }
            
            continue;
        }
        
        #if MRQ_DYN_UPDT_BNDS_ON_STRONG_BRANCH_CALCS
            if( updtBounds )
            {
                ////Maybe we should not set the variables, since it affect the pseudo curst calculation (original relaxations did not fix that variable). However, it can help us fathom by limit
                
                for( ; nextIndBoundSet < myNextInd; nextIndBoundSet++)
                {
                    const int index = intVars[nextIndBoundSet];
                    
                    r = nlp->setVariableBounds( index, nlx[index], nux[index] );
                    
                    #if MRQ_DEBUG_MODE
                        if( r != 0 )
                        {
                            MRQ_PRINTERRORNUMBER(r);
                        }
                    #endif
                }
            }
        #endif
        
        
        const double gap = MRQ_gap( nodeSol[ind] );
        
        
        if( gap <= intTol )
        {
            const double rSolInd = round( nodeSol[ind] );
            
            //we have a integer sol. So, we try strong branching if we are not in the bound of the variable... (remember variable not necessarily binary...). If the variable is binary, we only try strong branching in the opposite side of current value of solution...
            
            tryLeft = rSolInd > nlx[ind] && (pcost[myNextInd].nrealPl < maxComptSBranch || calculateEvenIfMaxComptSBranchReached);
            bleft = rSolInd - 1.0;
            
            tryRight = rSolInd < nux[ind] && (pcost[myNextInd].nrealPr < maxComptSBranch || calculateEvenIfMaxComptSBranchReached);
            bright = rSolInd + 1.0;
        }
        else
        {
            tryLeft = pcost[myNextInd].nrealPl < maxComptSBranch || calculateEvenIfMaxComptSBranchReached;
            bleft = floor( nodeSol[ind] );
            
            tryRight = pcost[myNextInd].nrealPr < maxComptSBranch || calculateEvenIfMaxComptSBranchReached;
            bright = bleft + 1.0;
        }
        
        
        //left branching
        
        lbind = nlx[ind];
        ubind = nux[ind];
        
        if( tryLeft )
        {
            nlp->setVariableBounds(ind, nlx[ind], bleft);
            
            
            r = nlp->solve( false, true, false, false);
            
            if( printLevel >= plevel )
            {
                std::cout << MRQ_PREPRINT  "Left st branching. ind: " << ind << " nodeSol: " << nodeSol[ind] << " retCode: " << r << " nlp.objF: " << nlp->objValue << " nlx: " << nlx[ind] << " nux: " << nux[ind] << " bleft: " << bleft << " " ;
                std::cout << " origretCode: " << nlp->origSolverRetCode << " ";
            }
            
            if( r == optsolvers::OPT_OPTIMAL_SOLUTION )
            {
                const double objValue = nlp->objValue;
                const double *sol = nlp->sol;
                
                if( objValue >= zu )
                {
                    //fathom by bound...
                    lbind = bleft + 1.0;
                    
                    if( printLevel >= plevel )
                        std::cout << "Fhatom by bound!";
                }
                else
                {
                    
                    if( MRQ_isIntegerSol(nI, intVars, sol, intTol) )
                    {
                        SEMAPH_solSem.lock(nThreads);
                        {
                            //we have to test again becuase the value of zu can have been changed
                            if( objValue < zu )
                            {
                                zu = objValue;
                                MRQ_copyArray(n, sol, intSol);
                                
                                intSolFound = true;
                            }
                        }
                        SEMAPH_solSem.unlock(nThreads);
                        
                        //fathom by feasibility
                        lbind = bleft + 1.0;
                        
                        if( printLevel >= plevel )
                            std::cout << "Fhatom by feasibility!";
                    }
                }
                
                if( storeFirstBranch )
                { // in general, we only calculate it on first BB iteration
                    fbranch[myNextInd].bPoint = bleft;
                    fbranch[myNextInd].leftlb = objValue;
                    
                    if( fbranch[myNextInd].leftlb > highestfbranch )
                        highestfbranch = fbranch[myNextInd].leftlb;
                }
                
                
                const double aux = MRQ_max(objValue - objNodeSol, 0.0)/( nodeSol[ind] - sol[ind] );  //we are in down branching. nodeSol[myNextInd]  is not greatter than sol[myNextInd]
                
                #if MRQ_DEBUG_MODE
                    assert(aux >= 0.0);
                #endif
                
                SEMAPH_sem.lock(totalThreads);
                {
                    pcost[myNextInd].nrealPl++;
                    pcost[myNextInd].nPl++;   //we first increase the number of pseudocost. do not put this increment after pl updation because we need this in this away to perform pseudo prune without mutex and in this way, we can be sure estimative will not be greather than it should be
                    
                    pcost[myNextInd].pl += aux;
                    
                    if( aux < pcost[myNextInd].minPl )
                        pcost[myNextInd].minPl = aux;
                }
                SEMAPH_sem.unlock(totalThreads);
                
                //std::cout << "\nobjNodeSol: " << objNodeSol << " " ;
                
            }
            else if( r == optsolvers::OPT_INFEASIBLE_PROBLEM || (!nlp->feasSol && considInfeasIfNlpFail) )
            {
                //fathom by infeasibility
                lbind = bleft + 1.0;
                
                if( printLevel >= plevel )
                    std::cout << "Fhatom by infeasibility! ";
            }
            
            if( printLevel >= plevel )
                std::cout << "\n";
            
            nlp->setVariableBounds( ind, nlx[ind], nux[ind] );
        }
        
        
        if( tryRight )
        {
            nlp->setVariableBounds( ind, bright, nux[ind] );
            
            r = nlp->solve( false, true, false, false );
            
            if( printLevel >= plevel )
            {
                std::cout << MRQ_PREPRINT << "Right st branching. ind: " << ind << " nodeSol: " << nodeSol[ind] << " retCode: " << r << " nlp.objF: " << nlp->objValue << " nlx: " << nlx[ind] << " nux: " << nux[ind] << " bright: " << bright << " ";
                std::cout << " origretCode: " << nlp->origSolverRetCode << " ";
            }
            
            
            if( r == optsolvers::OPT_OPTIMAL_SOLUTION )
            {
                const double objValue = nlp->objValue;
                const double *sol = nlp->sol;
                
                if( objValue >= zu )
                {
                    //fathom by bound...
                    ubind = bright - 1.0;
                    
                    if( printLevel >= plevel )
                        std::cout << "Fathom by bound! ";
                }
                else
                {
                    if( MRQ_isIntegerSol(nI, intVars, sol, intTol) )
                    {
                        SEMAPH_solSem.lock(nThreads);
                        {
                            //we have to test again becuase the value of zu can have been changed
                            if( objValue < zu )
                            {
                                zu = objValue;
                                MRQ_copyArray(n, sol, intSol);
                                intSolFound = true;
                            }
                        }
                        SEMAPH_solSem.unlock(nThreads);
                        
                        //fathom by feasibility
                        ubind = bright - 1.0;
                        
                        if( printLevel >= plevel )
                            std::cout << "Fathom by feasibility! ";
                    }
                    
                    if( storeFirstBranch )
                    {
                        // in general we only calculate it on first BB iteration
                        fbranch[myNextInd].bPoint = bright - 1.0; //may we are overwritting it, but ok...
                        fbranch[myNextInd].rightlb = objValue;
                        
                        if( fbranch[myNextInd].rightlb > highestfbranch )
                            highestfbranch = fbranch[myNextInd].rightlb;
                    }
                    
                    const double aux = MRQ_max(objValue - objNodeSol, 0.0)/( sol[ind] - nodeSol[ind] );
                    
                    #if MRQ_DEBUG_MODE
                        assert(aux >= 0.0);
                    #endif
                    
                    SEMAPH_sem.lock(totalThreads);
                    {
                        pcost[myNextInd].nrealPr++;
                        pcost[myNextInd].nPr++;  //we first increase the number of pseudocost. do not put this increment after pr updation because we need this in this away to perform pseudo prune without mutex and in this way, we can be sure estimative will not be greather than it should be
                        pcost[myNextInd].pr += aux;
                        
                        
                        if( aux < pcost[myNextInd].minPr )
                            pcost[myNextInd].minPr = aux;
                    }
                    SEMAPH_sem.unlock(totalThreads);
                    
                }
                
            }
            else if( r == optsolvers::OPT_INFEASIBLE_PROBLEM || ( considInfeasIfNlpFail && !nlp->feasSol) )
            {
                //fathom by infeasibility
                ubind = bright - 1.0;
                
                if( printLevel >= plevel )
                    std::cout << "Fhatom by infeasibility! ";
            }
            
            nlp->setVariableBounds( ind, nlx[ind], nux[ind] );
            
            if( printLevel >= plevel )
                std::cout << "\n";
            
        }
        
        
        if( ubind < lbind )
        {
            //we can consider the current node as infeasible...
            canFathom = true;
            //we could abort the calculations here...
        }
        else if( lbind != nlx[ind] || ubind != nux[ind] )
        {
            //we have to update the bounds vector of current node...
            
            
            //check if myBounds consider that index...
            
            SEMAPH_nodeSem.lock(nThreads);
            {
                const unsigned int nMyBounds = node->nMyBounds;
                bool found = false;
                auto &myBounds = node->myBounds;
                unsigned int nind;
                
                //here, we do not need worry about mutual exclusion because two threads do not see the same index ind
                for( unsigned int j = 0; j < nMyBounds; j++ )
                {
                    myBounds.getArrayElement(j, &nind, NULL, NULL, NULL);
                    
                    if( nind == ind )
                    {
                        //myBounds[j].l = lbind;
                        //myBounds[j].u = ubind;
                        //do not set sol here because we ned the origial value
                        
                        myBounds.setArrayElement(j, NULL, &lbind, &ubind, NULL);
                        
                        found = true;
                        break;
                    }
                }
                
                
                if( !found )
                {
                    if( !newBoundsAllocated )
                    {
                        r = myBounds.allocateElements(nI); //we allocate for nI positions to avoid have to reallocate... Note: allocate elements being done by myBounds directly does not change node_>nMyBounds
                        if(r != 0)
                        {
                            MRQ_PRINTERRORNUMBER(r);
                            code = MRQ_MEMORY_ERROR;
                            break;
                        }
                        newBoundsAllocated = true;
                    }
                    
                    
                    myBounds.addElementOnSortedArray( node->nMyBounds, ind, lbind, ubind, NAN );
                    
                    
                    if( nMyBounds < maxNMyBounds )
                        node->nMyBounds = nMyBounds + 1;
                    else
                    {
                        //so, we just let node->nMyBounds with its value. Note, we must lost some bounds fixing, but the algorithm will still be correct (I think... :/)
                        MRQ_PRINTERRORMSG("Warning: maximum number of bounds reached at node in strong branching calculation. Maybe you must set BBL_USE_SHORT_INTEGERS_TO_NUMBER_OF_MY_BOUNDS_IN_BB_NODE to 0");
                    }
                }
                
            }
            SEMAPH_nodeSem.unlock(nThreads);
            
            
            //node->print();
            //MRQ_getchar();
            
            
            nlx[ind] = lbind;
            nux[ind] = ubind;
            
            updtBounds = true;
            
        }
        
    }


    //std::cout << "pseudocust. updtBounds: " << updtBounds << "\n";


    /*if( abounds )
    {
        //just by safe, we set mutexes... abound->decPointerCounter can be called by other thread replacing the array in the current node...
        SEMAPH_nodeSem.lock(nThreads);
            abounds->decPointerCounter();
        SEMAPH_nodeSem.unlock(nThreads);
    } */


    return code;
}



int MRQ_NewPseudoCostCalc::calculateStrongBranching( const unsigned int totalThreads, const unsigned int nThreads, MRQ_MINLPProb &prob, const int nI, const int *intVars, const bool storeFirstBranch, const bool considInfeasIfNlpFail, const int printLevel, const double intTol, const double *olx, const double *oux, MRQ_NLPSolver **nlps, MRQ_NewBBNode *node, double *nlx, double *nux, const double *nodeSol, const double objNodeSol, const bool calculateEvenIfMaxComptSBranchReached, bool &updtBounds, bool &canFathom, double &zu, bool &intSolFound, double *intSol )
{
    bool newBoundsAllocated = false;
    int nextInd = 0;
    int code;
    int *thRetCodes = NULL;
    
    MRQ_Mutex indSem; //index mutex
    MRQ_Mutex solSem; //mutex fox zu and intSol updating
    MRQ_Mutex nodeSem; //for alteraltion in the node structure...
    
    
    #if MRQ_CPP_MULTITHREADING
        std::thread *myThreads = NULL;
    #endif
    
    
    updtBounds = false;
    canFathom = false;
    intSolFound = false;
    
    
    //thRetCodes = (int *) malloc( nThreads * sizeof(int) );
    MRQ_malloc(thRetCodes, nThreads);
    MRQ_IFMEMERRORGOTOLABEL(!thRetCodes, code, termination);
    
    
    
    #if MRQ_CPP_MULTITHREADING
    {
        if( nThreads > 1 )
        {
            myThreads = new (std::nothrow) std::thread[nThreads-1];
            MRQ_IFMEMERRORGOTOLABEL(!myThreads, code, termination);
            
            
            for( unsigned int i = 1; i < nThreads; i++ )
            {
                //you cannot pass arguments for thread as reference. Only copy or pointer...
                
                myThreads[i-1] = std::thread( MRQ_strongBranching, this, totalThreads, nThreads, i, &indSem, &solSem, &nodeSem, &nextInd, &prob, nI, intVars, storeFirstBranch, considInfeasIfNlpFail, printLevel, intTol, olx, oux, nlps[i], node, nlx, nux, nodeSol, objNodeSol, calculateEvenIfMaxComptSBranchReached, &newBoundsAllocated, &updtBounds, &canFathom, &zu, &intSolFound, intSol, &thRetCodes[i] );
            }
            
        }
        
        
        //using main thread also...
        threadStrongBranching( totalThreads, nThreads, 0, indSem, solSem, nodeSem, nextInd, prob, nI, intVars, storeFirstBranch, considInfeasIfNlpFail, printLevel, intTol, olx, oux, nlps[0], node, nlx, nux, nodeSol, objNodeSol, calculateEvenIfMaxComptSBranchReached, newBoundsAllocated, updtBounds, canFathom, zu, intSolFound, intSol, thRetCodes[0]);
        
        
        for( unsigned int i = 1; i < nThreads; i++)
            myThreads[i-1].join();
        
    }
    #elif MRQ_OMP_MULTITHREADING
    {
        omp_set_num_threads( nThreads );
        
        #pragma omp parallel
        {
            threadStrongBranching( totalThreads, nThreads, omp_get_thread_num(), indSem, solSem, nodeSem, nextInd, prob, nI, intVars, storeFirstBranch, considInfeasIfNlpFail, printLevel, intTol, olx, oux, nlps[omp_get_thread_num()], node, nlx, nux, nodeSol, objNodeSol, newBoundsAllocated, updtBounds, canFathom, zu, intSolFound, intSol, thRetCodes[omp_get_thread_num()] );
        }
    }
    #else
    {
        #if MRQ_DEBUG_MODE
            assert( nThreads == 1 );
        #endif
        
        threadStrongBranching(1, 1, 0, indSem, solSem, nodeSem, nextInd, prob, nI, intVars, storeFirstBranch, considInfeasIfNlpFail, printLevel, intTol, olx, oux, nlps[0], node, nlx, nux, nodeSol, objNodeSol, newBoundsAllocated, updtBounds, canFathom, zu, intSolFound, intSol, thRetCodes[0]);
    }
    #endif
    
    
    
termination:
    
    
    if(thRetCodes)	free(thRetCodes);
    
    #if MRQ_CPP_MULTITHREADING
        if(myThreads)
            delete[] myThreads;
            
    #endif
    
    
    return code;
}



//return true if the node can be prunned
bool MRQ_NewPseudoCostCalc:: updateVarBoundsByUpperBound( const int nI, const int *intVars, const double zu, double *nlx, double *nux )
{
    
    if( highestfbranch >= zu )
    {
        double b;
        
        for(int  i = 0; i < nI; i++)
        {
            int ind = intVars[i];
            
            if( fbranch[i].leftlb >= zu )
            { //we can prune the left ramification
                
                b = fbranch[i].bPoint + 1.0;
                
                if( b > nlx[ind] )
                {
                    #if MRQ_DEBUG_MODE
                        //std::cout << MRQ_PREPRINT "updating l bound on MRQ_NewPseudoCostCalc:: updateVarBoundsByUpperBound index: " << ind << " oldl: " << nlx[ind] << " newl: " << b << MRQ_GETFILELINE << "\n";
                    #endif
                        
                    nlx[ind] = b;
                    //MRQ_getchar();
                }
            }
            
            if( fbranch[i].rightlb >= zu )
            {
                b = fbranch[i].bPoint;
                
                if( b < nux[ind] )
                {
                    #if MRQ_DEBUG_MODE
                        //std::cout << MRQ_PREPRINT "updating u bound on MRQ_NewPseudoCostCalc:: updateVarBoundsByUpperBound index: " << ind << " oldu: " << nux[ind] << " newu: " << b << MRQ_GETFILELINE << "\n";
                    #endif
                    
                    nux[ind] = b;
                    //MRQ_getchar();
                }
            }
            
            if( nlx[ind] > nux[ind] )
                return true;
        }
        
    }
    
    
    return false;
}








MRQ_NewChooseIndexToBranch::MRQ_NewChooseIndexToBranch()
{
    //this->prob = prob;
    //this->intTol = intTol;
    //this->strategy = strategy;
    //this->pCosts = pcost;
    
    flagVars = NULL;
    //intVars = NULL;
    //binVars = NULL;
    //gap = NULL;
    //pcostEst = NULL;
    
}




MRQ_NewChooseIndexToBranch::~MRQ_NewChooseIndexToBranch()
{
    desallocate();
}



void MRQ_NewChooseIndexToBranch::desallocate()
{
    MRQ_secFree(flagVars);
}



int MRQ_NewChooseIndexToBranch::allocate(const int maxBranchVars, const int n)
{
    //if( maxBranchVars > 1 )
    {
        flagVars = (bool *) calloc( n, sizeof(bool) );
        
        if( !flagVars )
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTMEMERROR;
            #endif
            return MRQ_MEMORY_ERROR;
        }
    }
    
    return 0;
}



void MRQ_NewChooseIndexToBranch::getHighestPriorities( const int maxInds, const double intTol, const int* inds, const int sizeInds, const double *sol, const int *priorities, double *auxVars, int& ninds, unsigned int* choosenInds )
{
    double *gaps = auxVars;


    for(int i = 0; i < sizeInds; i++)
    {
        const int ind = inds[i];
        gaps[i] = MRQ_gap( sol[ind] );
    }


    for(ninds = 0; ninds < maxInds; )
    {
        int indMaxPrior = -1;
        int maxPrior = INT_MIN;
        
        for(int i = 0; i < sizeInds; i++)
        {
            const int ind = inds[i];
            
            if( !flagVars[ ind ] && gaps[i] > intTol )
            {
                if( priorities[ind] >= maxPrior)
                {
                    maxPrior = priorities[ind];
                    indMaxPrior = ind;
                }
            }
        }
        
        
        if( indMaxPrior >= 0 )
        {
            choosenInds[ninds] = indMaxPrior;
            ninds++;
            if( maxInds > ninds )
                flagVars[indMaxPrior] = true;
        }
        else
        { //no more variables having fractional vaules
            break;
        }
    }

}


void MRQ_NewChooseIndexToBranch::getHighestGaps( const int maxInds, const double intTol, const int* inds, const int sizeInds, const double *sol, double *auxVars, int& ninds, unsigned int* choosenInds )
{
    int indMaxGap;
    double maxGap;
    double *gaps = auxVars;


    for(int i = 0; i < sizeInds; i++)
    {
        const int ind = inds[i];
        gaps[i] = MRQ_gap( sol[ind] );
    }



    for(ninds = 0; ninds < maxInds ;  )
    {
        maxGap = 0.0;
        
        for( int i = 0; i < sizeInds; i++ )
        {
            const int ind = inds[i];
            
            if( !flagVars[ ind ] )
            {
                const double gap = gaps[i];
                
                if( gap > maxGap )
                {
                    maxGap = gap;
                    indMaxGap = ind;
                }
            }
        }
        
        if( maxGap > intTol )
        {
            choosenInds[ninds] = indMaxGap;
            ninds++;
            if( maxInds > ninds )
                flagVars[indMaxGap] = true;
        }
        else
        { //no more variables having fractional vaules
            break;
        }
    }


    for(int j = 0; j < ninds; j++)
        flagVars[ choosenInds[j] ] = false; //flagsVars should be set on false...
}



void MRQ_NewChooseIndexToBranch::getHighestPcosts( const int maxInds, const double intTol, const int nI, const int* intVars, const int *reverseIntVars, const double *sol, const double* nlx, const double* nux, const double pcost_mu, const MRQ_NewPseudoCost *pcost, double *auxVars, int& ninds, unsigned int* choosenInds )
{
    int indMaxPC;
    double maxPC;
    double *cpcost = auxVars;


    for(int k = 0; k < nI; k++)
    {
        const int ind = intVars[k];
        
        const int i = reverseIntVars[ind]; //we use reverseIntVars here instead of k directly because we can call this method from constraint branching passing a subset of intVars 
        const double gap = MRQ_gap( sol[ind] );
        
        if( nlx[ind] != nux[ind] && gap > intTol )
        {
            const double fsolind = floor(sol[ind]);
            
            const double fl = ( pcost[i].pl *(sol[ind] - fsolind) )/pcost[i].nPl ;
            
            const double fr = ( pcost[i].pr *(fsolind+1 - sol[ind]) )/pcost[i].nPr ;
            
            const double pcost = (1.0 - pcost_mu)*MRQ_min(fl, fr) + pcost_mu*MRQ_max(fl, fr) ;
            
            cpcost[k] = pcost;
        }
        else
        {
            cpcost[k] = -MRQ_INFINITY;
        }
        //std::cout << "cpcost["<<i<<"]: " << cpcost[i] << "\n";
    }

    //std::cout << "maxInds: " << maxInds << "\n";

    for(ninds = 0; ninds < maxInds; )
    {
        maxPC = -1.0;
        
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            if( !flagVars[ind] )
            {
                //const double gap = MRQ_gap( sol[ind] );
                //if( gap > intTol )
                
                if( cpcost[i] >= 0.0 )
                {
                    const double pcost =cpcost[i] ;
                    
                    if( pcost > maxPC )
                    {
                        maxPC = pcost;
                        indMaxPC = ind;
                    }
                }
            }
        }
        
        if( maxPC > -1.0 )
        {
            choosenInds[ninds] = indMaxPC;
            ninds++;
            if( maxInds > ninds )
                flagVars[indMaxPC] = true;
        }
        else
        { //no more variables having fractional vaules
            break;
        }
        
    }


    for(int j = 0; j < ninds; j++)
        flagVars[ choosenInds[j] ] = false; //flagsVars should be set on false..

}


void MRQ_NewChooseIndexToBranch::getHighestPcostsEvenIntegerValues( const int maxInds, const int nI, const int* intVars, const int *reverseIntVars, const double *sol, const double* nlx, const double* nux, const double pcost_mu, const MRQ_NewPseudoCost *pcost, double *auxVars, int& ninds, unsigned int* choosenInds )
{
    
    return getHighestPcosts(maxInds, -INFINITY, nI, intVars, reverseIntVars, sol, nlx, nux, pcost_mu, pcost, auxVars, ninds, choosenInds); //we just call getHighestGaps passing -Infinity as integer tol
    
    
    #if 0
    int indMaxPC;
    double maxPC;
    double *cpcost = auxVars;


    for(int k = 0; k < nI; k++)
    {
        const int ind = intVars[k];
        
        const int i = reverseIntVars[ind];
        if( nlx[ind] != nux[ind] )
        {
            const double fsolind = floor(sol[ind]);
            
            const double fl = (pcost[i].pl *( sol[ind] - fsolind ) )/pcost[i].nPl ;
            
            const double fr = (pcost[i].pr *( fsolind + 1.0 - sol[ind] ))/pcost[i].nPr ;
            
            const double pcost = (1.0 - pcost_mu)*MRQ_min(fl, fr) + pcost_mu*MRQ_max(fl, fr) ;
            
            cpcost[k] = pcost;
        }
    }



    for(ninds = 0; ninds < maxInds; )
    {
        maxPC = -1.0;
        
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            if( !flagVars[ind] )
            {
                //const double gap = MRQ_gap( sol[ind] );
                
                if( nlx[ind] != nux[ind] )
                {
                    const double pcost = cpcost[i];
                    
                    if( pcost > maxPC )
                    {
                        maxPC = pcost;
                        indMaxPC = ind;
                    }
                }
            }
        }
        
        if( maxPC > -1.0 )
        {
            choosenInds[ninds] = indMaxPC;
            ninds++;
            if( maxInds > ninds )
                flagVars[indMaxPC] = true;
        }
        else
        { //no more variables having fractional vaules
            break;
        }
        
    }


    for(int j = 0; j < ninds; j++)
        flagVars[ choosenInds[j] ] = false; //flagsVars should be set on false..
    #endif
}



void MRQ_NewChooseIndexToBranch::chooseIndices( const MRQ_NewPseudoCost *pcost, const int strategy, const double intTol, const int nI, const int *intVars, const int *reverseIntVars, const int nbin, const int *binVars, const int *nonBinVars, const double pcost_mu, const int maxInds, const double* relaxSol, const bool nlpSucess, const double* nlx, const double* nux, const int *priorities, double *auxVars, int& ninds, unsigned int* choosenInds, bool &branchOnInt )
{
    int i, j, sizeIndices;
    const int *indices;
    
    branchOnInt = false;
    
    /*printf("nI: %d nbin: %d\n", nI, nbin);
    for(i = 0; i < nbin; i++)
        printf("binVars[%2d]: %d\t", i, (int) binVars[i]);
    printf("\n"); */
    
    
    if( strategy == MRQ_BB_BS_STBRANCH_PSEUDO_COSTS )
    {
        //pseudocost strategy
        
        indices = intVars; //do not put binVars here, because in pseudo cust calculations, we assume psudo custs are in the same order that intVars
        sizeIndices = nI;
        
        
        getHighestPcosts(maxInds, intTol, nI, intVars, reverseIntVars, relaxSol, nlx, nux, pcost_mu, pcost, auxVars, ninds, choosenInds );
        
        //now, we take advatange of getHighestGaps
        
        
        if( ninds == 0 )
        {
            if( nlpSucess )
            {
                MRQ_PRINTERRORMSG("warning: non fractional variable to perform branching."); //we still have an strange case where we have a integer solution, but objValue is greater than zu and dualObjValue is lower than zu (in this case, if we are consider dualObjValue to prune by bound, node will not be pruno). This strange case lead us here and can be normal, for example if difference between objValue and dualObjValue is high (in early_branching mode, for example)
                #if MRQ_DEBUG_MODE
                {
                    double *cpcost = auxVars; //auxVars has the calculated psudo cost for variables...
                    
                    for(int i = 0; i < nI; i++)
                    {
                        int ind = intVars[i];
                        
                        std::cout << "sol["<<ind<<"]: " << relaxSol[ind] << " pcost["<<i<<"]:: pl: " << pcost[i].pl << " pr: " << pcost[i].pr << " nPl: " << pcost[i].nPl << " nPr: " << pcost[i].nPr << " cpcost["<<i<<"]: " << cpcost[i] <<  "  \t";
                    }
                    std::cout << "\n";
                    
                    assert(false); // we just put this assert here in the debug mode to me see this error and so evaluate the situation.
                }
                #endif
            }
            
            //we choose a non fixed integer variable. It occours because nlp solver got only a feasible solution without certified optimality. So, we choose a variable taking an integer value
            
            getHighestPcostsEvenIntegerValues(maxInds, nI, intVars, reverseIntVars, relaxSol, nlx, nux, pcost_mu, pcost, auxVars, ninds, choosenInds);
            
            branchOnInt = true;
        }
        
    }
    else if( strategy == MRQ_BB_BS_HIGHEST_INT_GAP || strategy == MRQ_BB_BS_BIN_FIRST_HIGHEST_INT_GAP )
    {
        
        if( strategy == MRQ_BB_BS_HIGHEST_INT_GAP )
        {
            indices = intVars;
            sizeIndices = nI;
        }
        else //strategy == MRQ_BB_BS_BIN_FIRST_HIGHEST_INTGAP
        {
            indices = binVars;
            sizeIndices = nbin;
        }
        
        
        getHighestGaps(maxInds, intTol, indices, sizeIndices, relaxSol, auxVars, ninds, choosenInds);
        
        
        
        //printf("ninds: %d maxInds: %d\n", ninds, maxInds);
        
        if( ninds < maxInds && strategy == MRQ_BB_BS_BIN_FIRST_HIGHEST_INT_GAP )
        {
            //we need search on general integer variables
            
            indices = nonBinVars;
            sizeIndices = nI - nbin;
            
            getHighestGaps(maxInds-ninds, intTol, indices, sizeIndices, relaxSol, auxVars, i, &choosenInds[ninds]);
            
            
            ninds += i;
        }
        
        
        
        if( ninds == 0 )
        {
            if( nlpSucess )
            {
                MRQ_PRINTERRORMSG("warning: non fractional variable to perform branching\n"); //we still have an strange case where we have a integer solution, but objValue is greater than zu and dualObjValue is lower than zu (in this case, if we are consider dualObjValue to prune by bound, node will not be pruno). This strange case lead us here and can be normal, for example if difference between objValue and dualObjValue is high (in early_branching mode, for example)
                #if MRQ_DEBUG_MODE
                    assert(false); // we just put this assert here in the debug mode to me see this error and so evaluate the situation.
                #endif
            }
            
            //we choose a non fixed integer variable. It occours because nlp solver got only a feasible solution without certified optimality. So, we choose a variable taking an integer value 
            
            for(i = 0; i < sizeIndices; i++)
            {
                j = indices[i];
                if( nlx[j] != nux[j] )
                {
                    choosenInds[ninds] = j;
                    ninds = 1;
                    break; //we  take only 1 indice
                }
            }
            
            if( ninds == 0 && strategy == MRQ_BB_BS_BIN_FIRST_HIGHEST_INT_GAP)
            {
                sizeIndices = nI - nbin;
                
                for(i = 0; i < sizeIndices; i++)
                {
                    j = nonBinVars[i];
                    if( nlx[j] != nux[j] )
                    {
                        choosenInds[ninds] = j;
                        ninds = 1;
                        break; //we  take only 1 indice
                    }
                }
            }
            
            branchOnInt = true;
        }
        
        
        
    }
    else if( strategy == MRQ_BB_BS_VAR_PRIORITIES )
    {
        getHighestPriorities( maxInds, intTol, intVars, nI, relaxSol, priorities, auxVars, ninds, choosenInds );
        
        if( ninds == 0 )
        {
            if( nlpSucess )
            {
                MRQ_PRINTERRORMSG("warning: non fractional variable to perform branching\n"); //we still have an strange case where we have a integer solution, but objValue is greater than zu and dualObjValue is lower than zu (in this case, if we are consider dualObjValue to prune by bound, node will not be pruno). This strange case lead us here and can be normal, for example if difference between objValue and dualObjValue is high (in early_branching mode, for example)
                #if MRQ_DEBUG_MODE
                    assert(false); // we just put this assert here in the debug mode to me see this error and so evaluate the situation.
                #endif
            }
            
            int indMaxPrior = -1, maxPrior = INT_MIN;
            
            
            for(int i = 0; i < nI; i++)
            {
                const int ind = intVars[i];
                
                if( nlx[ind] != nux[ind] )
                {
                    if( priorities[ind] >= maxPrior )
                    {
                        maxPrior = priorities[ind];
                        indMaxPrior = ind;
                    }
                }
            }
            
            #if MRQ_DEBUG_MODE
                assert( indMaxPrior >= 0 );
            #endif
            
            ninds = 1;
            choosenInds[0] = indMaxPrior;
            
            branchOnInt = true;
            
            
        }
    }
    else
    {
        MRQ_PRINTERRORMSG("branching strategy not implemented");
        assert(false);
    }
    
    
    /*#if MRQ_DEBUG_MODE
    if( ninds <= 0 )
    {//error
        
        std::cout << "nlpSucess: " << nlpSucess << std::endl;
        
        for(int i = 0; i < nI; i++)
        {
            std::cout << "i: " << i << " intVar: " << intVars[i] << " sol["<<intVars[i]<<"]: " << relaxSol[intVars[i]] << " pcost.nPl: " << pcost[i].nPl << " pcost.pl: " << pcost[i].pl << " pcost.nPr: " << pcost[i].nPr << " pcost.pr: " << pcost[i].pr << " cpcost: " << auxVars[i] << std::endl;
        }
        
        
        assert(false); //ninds should be > 0
    }
    #endif*/
    
    
}






void MRQ_NewPoints::initialize(const int nPreAlloc)
{
    nPoints = 0;
    nPointsAllocated = 0;
    this->nPreAlloc = nPreAlloc;
    points = NULL;
}



MRQ_NewPoints::MRQ_NewPoints(const int nPreAlloc)
{
    initialize(nPreAlloc);
}

MRQ_NewPoints::~MRQ_NewPoints()
{
    desallocate();
}


int MRQ_NewPoints::preAllocate(const int *numberOfPreAlloc)
{
    double **auxP;
    
    int nn = numberOfPreAlloc ? *numberOfPreAlloc : nPreAlloc;
    
    
    //if points is NULL, realloc behaves like malloc...
    auxP = (double **) realloc( points, (nPointsAllocated + nn) * sizeof(double *) );
    
    if(!auxP)
    {
        #if MRQ_MEMORY_ERROR
            MRQ_PRINTMEMERROR;
        #endif
        return MRQ_MEMORY_ERROR;
    }
    
    points = auxP;
    nPointsAllocated += nn;
    
    return 0;
}


//warning: that function returns the number of points allocated...
int MRQ_NewPoints::addPoints(const int nPoints, const int dim, double** x)
{
    int i, row;
    int *nn;
    
    
    if( this->nPoints + nPoints > nPointsAllocated )
    {
        if( nPoints > nPreAlloc )
        {
            i = nPoints;
            nn = &i;
        }
        else
            nn = NULL;
        
        i = preAllocate(nn);
        if(i != 0)
            return 0;
    }
    
    
    for(i = 0; i < nPoints ; i++)
    {
        row = this->nPoints + i;
        points[row] = (double *) malloc( dim * sizeof(double) );
        
        if(!points[row])
            break;
        
        MRQ_copyArray( dim, x[i], points[row] );
        
        this->nPoints++;
    }
    
    
    return i;
}


void MRQ_NewPoints::desallocate()
{
    int i;
    
    if( points )
    {
        for(i = 0; i < nPoints; i++)
            free( points[i] );
        
        free(points);
        //points = NULL;
    }
    
    initialize(nPreAlloc);
}












MRQ_NewQuadraticCut::MRQ_NewQuadraticCut()
{
    anz = 0;
    Qnz = 0;
    
    acols = NULL;
    avals = NULL;
    
    Qrows = NULL;
    Qcols = NULL;
    Qvals = NULL;
}



MRQ_NewQuadraticCut::~MRQ_NewQuadraticCut()
{
    desallocateLinearPart();
    desallocateQuadraticPart();
}



int MRQ_NewQuadraticCut::addCutOnSolver( MRQ_NLPSolver* nlp)
{
    bool added = false;
    int ind, r, code;
    
    
    r = nlp->getNumberOfConstraints( ind );
    #if MRQ_DEBUG_MODE
        if( r != 0 )
        {
            MRQ_PRINTERRORNUMBER(r);
            code = MRQ_NLP_SOLVER_ERROR;
            goto termination;
        }
    #endif
    
    
    r = nlp->addConstraints(1);
    if( r != 0 )
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
        code = MRQ_NLP_SOLVER_ERROR;
        goto termination;
    }
    
    added = true;
    
    r = nlp->setConstraintLinearCoefs(ind, anz, acols, avals);
    if( r != 0 )
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
        code = MRQ_NLP_SOLVER_ERROR;
        goto termination;
    }
    
    
    r = nlp->setConstraintQuadMatrix(ind, Qnz, Qrows, Qcols, Qvals);
    if( r != 0 )
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
        code = MRQ_NLP_SOLVER_ERROR;
        goto termination;
    }
    
    
    r = nlp->setConstraintBounds(ind, lb, ub);
    #if MRQ_DEBUG_MODE
        if( r != 0 )
        {
            MRQ_PRINTERRORNUMBER(r);
            code = MRQ_NLP_SOLVER_ERROR;
            goto termination;
        }
    #endif
    
    
    code = 0;
    
termination:
    
    if(code != 0 && added)
    {
        nlp->removeConstraints(1, &ind);
    }
    
    return 0;
}



int MRQ_NewQuadraticCut::allocateLinearPart( const int nz)
{
    desallocateLinearPart();
    
    MRQ_malloc(acols, nz); //acols = (int *) malloc(nz * sizeof(int));
    MRQ_malloc(avals, nz); //avals = (double *) malloc(nz * sizeof(double));
    
    if(!acols || !avals)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTMEMERROR;
        #endif
        return MRQ_MEMORY_ERROR;
    }
    
    return 0;
}



int MRQ_NewQuadraticCut::allocateQuadraticPart( const int nz)
{
    desallocateQuadraticPart();
    
    MRQ_malloc(Qrows, nz); //Qrows = (int *) malloc( nz * sizeof(int) );
    MRQ_malloc(Qcols, nz); //Qcols = (int *) malloc( nz * sizeof(int) );
    MRQ_malloc(Qvals, nz); //Qvals = (double *) malloc(nz * sizeof(double) );
    
    if(!Qrows || !Qcols || !Qvals)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTMEMERROR;
        #endif
        return MRQ_MEMORY_ERROR;
    }
    
    return 0;
}



void MRQ_NewQuadraticCut::desallocateLinearPart()
{
    MRQ_secFree(acols);
    MRQ_secFree(avals);
}


void MRQ_NewQuadraticCut::desallocateQuadraticPart()
{
    MRQ_secFree(Qrows);
    MRQ_secFree(Qcols);
    MRQ_secFree(Qvals);
}



void MRQ_NewQuadraticCut::setBounds(const double lb, const double ub)
{
    this->lb = lb;
    this->ub = ub;
}



int MRQ_NewQuadraticCut::setLinearPart(const int nz, const int* cols, const double* vals)
{
    int code;
    
    code = allocateLinearPart(nz);
    
    if(code != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(code);
        #endif
        return code;
    }
    
    anz = nz;
    MRQ_copyArray(nz, cols, acols);
    MRQ_copyArray(nz, vals, avals);
    
    
    return 0;
}


int MRQ_NewQuadraticCut::setQuadraticPart(const int nz, const int* rows, const int* cols, const double* vals)
{
    int code;
    
    code = allocateQuadraticPart(nz);
    
    if(code != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(code);
        #endif
        return code;
    }
    
    anz = nz;
    
    MRQ_copyArray(nz, rows, Qrows);
    MRQ_copyArray(nz, cols, Qcols);
    MRQ_copyArray(nz, vals, Qvals);
    
    return 0;
}





MRQ_NewQuadCutStorer::MRQ_NewQuadCutStorer(int nThreads, unsigned int cutId)
{
    remainThreads = nThreads;
    this->cutId = cutId;
}



int MRQ_NewQuadCutStorer::decThreadCounter()
{
    int aux;
    
    remainThreads--;
    
    aux = remainThreads;
    
    //suicide. No threads have to set that cut more...
    if(aux <= 0)
        delete this;
    
    return aux; //there is no remainThreads more...
}




int MRQ_NewQuadCutStorer::setQuadraticCut(const int nza, const int* acols, const double* avalues, const int nzQ, const int* Qrows, const int* Qcols, const double* Qvalues, const double lb, const double ub)
{
    int r;
    
    
    cut.setBounds(lb, ub);
    
    r = 0;
    if(nza > 0)
    {
        r = cut.setLinearPart(nza, acols, avalues);
        if(r != 0)
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERRORNUMBER(r);
            #endif
        }
    }
    
    if(nzQ > 0)
    {
        r = cut.setQuadraticPart(nzQ, Qrows, Qcols, Qvalues);
        
        if(r != 0)
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERRORNUMBER(r);
            #endif
            return MRQ_MEMORY_ERROR;
        }
    }
    
    return 0;
}









MRQ_NewGlobalCutList::MRQ_NewGlobalCutList()
{
}


MRQ_NewGlobalCutList::MRQ_NewGlobalCutList(const int numberOfThreads)
{
    initialize(numberOfThreads);
}


int MRQ_NewGlobalCutList::initialize(const int numberOfThreads)
{
    isempty= true;
    nThreads = numberOfThreads;
    //sem = NULL;
    //mutexInit = false;
    
    return 0;//initializeMutex();
}




MRQ_NewGlobalCutList::~MRQ_NewGlobalCutList()
{
    eraseList();
    //deleteMutex();
}


int MRQ_NewGlobalCutList::addCutsOnSolver( MRQ_NLPSolver *nlp )
{
    std::list< MRQ_NewQuadCutStorer* >:: iterator it;
    
    //MRQ_lockMutex(nThreads, sem);
    SEMAPH.lock(nThreads);
        for( it = cutList.begin(); it != cutList.end(); it++ )
        {
            (*it)->cut.addCutOnSolver(nlp);
            cutGen->decThreadCounter( *it );
            
            //printf("Adicionei corte global no solver nlp\n");
        }
        
        cutList.clear();
        
        isempty = true;
    SEMAPH.unlock(nThreads);
    //MRQ_unlockMutex(nThreads, sem);
    
    return 0;
}




int MRQ_NewGlobalCutList::insertCutStorer( MRQ_NewQuadCutStorer* cutStorer)
{
    //insert the cut at the end of the encadeate std::list
    //MRQ_lockMutex(nThreads, SEMAPH);
    SEMAPH.lock( nThreads );
        cutList.insert( cutList.end() , cutStorer );
        isempty = false;
    SEMAPH.unlock( nThreads );
    //MRQ_unlockMutex(nThreads, SEMAPH);
    
    return 0;
}




void MRQ_NewGlobalCutList::eraseList()
{
    std::list< MRQ_NewQuadCutStorer* >::iterator it;
    
    
    //MRQ_lockMutex(nThreads, SEMAPH);
    SEMAPH.lock(nThreads);
        for( it = cutList.begin(); it != cutList.end(); it++ )
        {
            cutGen->decThreadCounter( *it );
            //delete (*it);
        }
        
        cutList.clear();
        isempty = true;
    SEMAPH.unlock(nThreads);
    //MRQ_unlockMutex(nThreads, SEMAPH);
}















MRQ_NewGlobalCutGenerator::MRQ_NewGlobalCutGenerator(const int nThreads, MRQ_NewGlobalCutList *cutLists)
{
    nextId = 0;
    //mutexInit = false;
    this->nThreads = nThreads;
    
    
    /*if(nThreads > 1)
        initializeMutex(nThreads);
    else
        SEMAPH = NULL; */
    
    this->cutLists = cutLists;
}


MRQ_NewGlobalCutGenerator::~MRQ_NewGlobalCutGenerator()
{
    //deleteMutex();
}



int MRQ_NewGlobalCutGenerator::decThreadCounter( MRQ_NewQuadCutStorer* cutStorer)
{
    int aux;
    
    //MRQ_lockMutex(nThreads, SEMAPH);
    SEMAPH.lock(nThreads);
        aux = cutStorer->decThreadCounter(); //when the counter get 0, the cutStorer is deleted.
    SEMAPH.unlock(nThreads);
    //MRQ_unlockMutex(nThreads, SEMAPH);
    
    return aux;
}



int MRQ_NewGlobalCutGenerator::setQuadraticCut( const int nza, const int *acols, const double *avalues, const int nzQ, const int *Qrows, const int* Qcols, const double *Qvalues, const double lb, const double ub)
{
    MRQ_NewQuadCutStorer *cutStorer;
    unsigned int id;
    int r;
    
    
    //MRQ_lockMutex(nThreads, SEMAPH);
    SEMAPH.lock(nThreads);
        id = nextId;
        nextId++;
    SEMAPH.unlock(nThreads);
    //MRQ_unlockMutex(nThreads, SEMAPH);
    
    
    cutStorer = new (std::nothrow) MRQ_NewQuadCutStorer(nThreads, id);
    
    if(!cutStorer)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTMEMERROR;
        #endif
        return MRQ_MEMORY_ERROR;
    }
    
    
    r = cutStorer->setQuadraticCut(nza, acols, avalues, nzQ, Qrows, Qcols, Qvalues, lb, ub);
    if(r != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
        delete cutStorer;
        return MRQ_MEMORY_ERROR;
    }
    
    
    for(int i = 0; i < nThreads; i++)
        cutLists[i].insertCutStorer(cutStorer);
    
    
    return 0;
}




MRQ_BinSumConstrsInds::MRQ_BinSumConstrsInds()
{
    initialize();
}


MRQ_BinSumConstrsInds::~MRQ_BinSumConstrsInds()
{
    deallocate();
}


void MRQ_BinSumConstrsInds::deallocate()
{
    MRQ_secFree(binSumConstrs);
    nbinSumConstrs = 0;
}


void MRQ_BinSumConstrsInds::initialize()
{
    nbinSumConstrs = 0;
    binSumConstrs = NULL;
}


int MRQ_BinSumConstrsInds::_calculateIndices(const MRQ_MINLPProb &prob, const double *lx, const double *ux, const double *lc, const double *uc )
{
    //int r;
    const int m = prob.m;
    const bool *nlConstr = prob.nlConstr;
    const int *xtype = prob.xtype;

    const minlpproblem::MIP_SparseMatrix *QC = prob.QC;
    const minlpproblem::MIP_SparseMatrix &A = prob.A;

    const double ONE = 1.0;
    


    //detecting constraints: sum_{i \in I} x_i = 1
    //for(int i = 0; i < m; i++)
    for( minlpproblem::MIP_SparseMatrixRowIndexIterator rowait = A.beginRowIndex() ;  *rowait < (unsigned int) m ; ++rowait )
    {
        const int i = *rowait;
        
        if( lc[i] <= -MIP_INFINITY && uc[i] >= MIP_INFINITY )  //free constraint. This constraint was discarded by user or by the preprocessor.
            continue;
        
        if(nlConstr[i] || QC[i].getNumberOfElements() > 0)
            continue;
        
        //minlpproblem::MIP_SparseRow &row = A[i];
        const int* acols = A[i];
        const double* avalues = A(i);
        const unsigned int nel = A.getNumberOfElementsAtRow(i); // row.getNumberOfElements();
        
        unsigned int nPossibleBinVarWithCoefOne = 0;
        
        
        //if( row.getNumberOfElements() == 0 )
            //continue;
        
        #if MRQ_DEBUG_MODE
            assert( nel > 0 ); //there must be at least some element in this row of A (we are running rowIndexIterator)
        #endif
        
        
        bool binSumConstr = true;
        
        for(unsigned int j = 0; j < nel; j++)
        {
            const int var = acols[j]; //row[j].getColumn();
            
            
            //if variable is fixed, we can consider this constraint even if this variable is integer or continuous
            if( lx[var] != ux[var] )
            {
                const bool intVar = minlpproblem::MIP_isIntegerType(xtype[var]);
                //by now, we only consider nonfixed variables having coeficient 1.0 in the constraint
                
                const double coef = avalues[j];
                
                //we just enable constraint branching consideration for this constraint if this variable is integer and has integer coeficient in this constraint
                if( intVar &&  MRQ_isIntegerDouble(coef) )
                {
                    //we are trying detect something general here. If a integer variavel has lower bound lower than 2.0 and upper bound greater tnah -1.0, we have an integer variable that could be binary in some(s) (sub)partitions of the space. Our objective is making possible use constraint branhching over this integer variable in partitions where this integer variable will be binary
                    if( coef == ONE  &&  lx[var] < 2.0 && ux[var] > -1.0 )
                        nPossibleBinVarWithCoefOne++;
                    
                    //even if a general integer variable cannot be treat like binary in any partition, we still can perform constraint branching over this constraint if this varibale gets fixed and updated rhs be one. So, we do not discard this contsraint 
                }
                else
                {
                    binSumConstr = false;
                    break;
                }
                
            }
            
        }
        
        if( binSumConstr && nPossibleBinVarWithCoefOne > 1 ) //we must have at least 2 possible binary variable in constraint to perform constraint branching
        {
            binSumConstrs[nbinSumConstrs] = i;
            nbinSumConstrs++;
        }
        
    }

    return 0;
}


int MRQ_BinSumConstrsInds::calculateIndices( const MRQ_MINLPProb& prob, const double* lx, const double* ux, const double* lc, const double* uc )
{
    int r;
    const int m = prob.m;
    
    if( !lx )
        lx = prob.lx;
    if( !ux )
        ux = prob.ux;
    if( !lc )
        lc = prob.lc;
    if( !uc )
        uc = prob.uc;
    
    nbinSumConstrs = 0;
    
    
    r = reallocate(m, binSumConstrs);
    MRQ_IFERRORRETURN(r, r);
    
    
    r = _calculateIndices(prob, lx, ux, lc, uc);
    MRQ_IFERRORRETURN(r, r);
    
    
    if( nbinSumConstrs != m )
    {
        reallocate(nbinSumConstrs, binSumConstrs);
    }
    
    
    return 0;
}



int MRQ_BinSumConstrsInds::reallocate( unsigned int size, int* &inds )
{
    return MRQ_reallocateArray(size, inds);
}


#if 0
MRQ_BinSumConstrsInds2::MRQ_BinSumConstrsInds2() : MRQ_BinSumConstrsInds()
{
}


MRQ_BinSumConstrsInds2::~MRQ_BinSumConstrsInds2()
{
}


int MRQ_BinSumConstrsInds2::_calculateIndices(const MRQ_MINLPProb &prob, const double *lx, const double *ux, const double *lc, const double *uc )
{
    int r;
    const int m = prob.m;
    //const bool *nlConstr = prob.nlConstr;
    const int *xtype = prob.xtype;
    
    //const minlpproblem::MIP_SparseMatrix *QC = prob.QC;
    const minlpproblem::MIP_SparseMatrix &A = prob.A;
    
    const double ONE = 1.0;
    
    
    
    r = MRQ_BinSumConstrsInds::_calculateIndices(prob, lx, ux, lc, uc);
    if(r != 0)
    {
        MRQ_PRINTERRORNUMBER(r);
        return r;
    }
    
    
    //now, we detect the constraints in the form: 
    //y_{k_1} + y{k_2} + ... + y_{k_n-1} - y{k_n} =  0
    
    
    for( minlpproblem::MIP_SparseMatrixRowIndexIterator rowait = A.beginRowIndex() ;  *rowait < (unsigned int) m ; ++rowait )
    {
        const int i = *rowait;
        
        if(prob.isLinearConstraint(i) == false)
            continue;
        
        
        if( lc[i] != 0.0 ||  uc[i] != 0.0 )
            continue; //both lc[i] and uc[i] should be zero
        
        int negVars = 0;
        
        const int* acols = A[i];
        const double* avalues = A(i);
        const unsigned int nel = A.getNumberOfElementsAtRow(i); 
        
        bool binSumConstr = true;
        
        for(unsigned int j = 0; j < nel; j++)
        {
            const int var = acols[j];
            const double value = avalues[j];
            
            
            //it is ok if we have a non binary variable fixed at zero here
            if( lx[var] == 0.0 && ux[var] == 0.0 )
                continue;
            
            
            if( !minlpproblem::MIP_isBinaryVar(xtype[var], lx[var], ux[var]) )
            {
                binSumConstr = false;
                break;
            }
            
            if( value == -ONE )
            {
                /* now, we detect constraint having any number of unitary negative coeeficites
                if( negVars >= 1 )
                {//we can have only one variable having negative index
                    binSumConstr = false;
                    break;
                }*/
                negVars++; //we do not need more this var, but ok until we be sure about this strategy
            }
            else if( value != ONE )
            {
                binSumConstr = false;
                break;
            }
            
        }
        
        if(binSumConstr)
        {
            binSumConstrs[nbinSumConstrs] = i;
            nbinSumConstrs++;
        }
    }
    
    
    
    
    return 0;
}

#endif




MRQ_BinSumConstrsChooser::MRQ_BinSumConstrsChooser()
{
    initialize();
}


MRQ_BinSumConstrsChooser::~MRQ_BinSumConstrsChooser()
{
    desallocate();
}


void MRQ_BinSumConstrsChooser::desallocate()
{
    MRQ_secFree(candInds);
    nCandInds = 0;
}


void MRQ_BinSumConstrsChooser::initialize()
{
    nCandInds = 0;
    candInds = NULL;
}


int MRQ_BinSumConstrsChooser::reallocate(unsigned int size)
{
    return MRQ_reallocateArray(size, candInds);
}


int MRQ_BinSumConstrsChooser:: calculateCandidateConstrsToBranch( const int nbinSumConstrs, const int *binSumConstrs, const double *lx, const double *ux, const int *varType, const minlpproblem::MIP_SparseMatrix& A, const double *lc, const double *uc, const double intTol, const double* sol)
{
    nCandInds = 0;

    for(int i = 0; i < nbinSumConstrs; i++)
    {
        bool fracBinaryVar = false;
        bool candidateConstr = true;
        
        const int c = binSumConstrs[i];
        
        int nonFixedNegCoefVars = 0;
        //double rhs = uc[c]; //we only consider equality constraints or lower equal constraints
        double fixed = 0.0;
        
        const int *acols = A[c];
        const double *avalues = A(c);
        const unsigned int nel = A.getNumberOfElementsAtRow(c);
        
        for(unsigned int j = 0; j < nel; j++)
        {
            const int col = acols[j];
            
            
            if( lx[col] == ux[col] )
            { //fixed variable
                fixed += avalues[j] *lx[col];
            }
            else
            {
                const double s = sol[col];
                
                
                #if MRQ_DEBUG_MODE
                    assert( minlpproblem::MIP_isIntegerType(varType[col])  ); //this variable should be integer
                #endif
                
                if( !(lx[col] > -1.0 && ux[col] < 2.0) )
                { //we have a non-fixed non-binary variable (general integer variable). So, we cannot perform branching over this constraint
                    candidateConstr = false;
                    break;
                }
                
                
                if( avalues[j] == -1.0 )
                {
                    nonFixedNegCoefVars++;
                }
                else if( avalues[j] != 1.0 )
                {
                    //we have a nonfixed variable having coeficient different of -1.0 and 1.0. We cannot perform constraint branching
                    candidateConstr = false;
                    break;
                }
                
                
                if( MRQ_binGap(s) > intTol )
                {
                    #if MRQ_DEBUG_MODE
                        assert( lx[col] > -1.0 && ux[col] < 2.0 );
                    #endif
                    
                    fracBinaryVar = true; //we cannot break here because we need calculate the fixed value
                }
            }
        }
        
        
        if( candidateConstr && fracBinaryVar )
        {
            double rhs = uc[c] - fixed;
            double lhs = lc[c] - fixed;
            
            
            //here, we are trusting we have no numerical errors, but maybe it is not true...
                
            /* We can perform contsraint branchinhg over the following type of constraints:
            * 
            * y_{k_1} + y{k_2} + ... + y_{k_n-1} = 1        (1)
            * y_{k_1} + y{k_2} + ... + y_{k_n-1} <= 1       (2)
            * y_{k_1} + y{k_2} + ... + y_{k_n-1} = y{k_n}   (3)
            * y_{k_1} + y{k_2} + ... + y_{k_n-1} <= y{k_n}  (4)
            */
            if( (rhs == 1.0 && lhs == 1.0 && nonFixedNegCoefVars == 0)  ||
                (rhs == 1.0 && lhs < 1.0  && nonFixedNegCoefVars == 0)  ||
                (rhs == 0.0 && lhs == 0.0 && nonFixedNegCoefVars == 1)  ||  
                (rhs == 0.0 && lhs < 0.0  && nonFixedNegCoefVars == 1)  )
            { //so we can perform branch over this constraint since we have a sum of binary vars, at least one binary var is fractional and final rhs is 1.0
                candInds[nCandInds] = c;
                nCandInds++;
            }
        }
        
    }


    return nCandInds;
}


//to choose index in strategies MRQ_BB_CBS_LOWEST_NUMBER_OF_VARS and MRQ_BB_CBS_HIGHEST_NUMBER_OF_VARS
inline bool MRQ_BinSumConstrsChooser:: chooseIndexToBranchLowHigh( const minlpproblem::MIP_SparseMatrix& A, const int nbinSumEqConstrs, const int* binSumEqConstrs, const double intTol, const double* sol, const int branchStrategy, int& index ) const
{
    const bool compLow = branchStrategy == MRQ_BB_CBS_LOWEST_NUMBER_OF_VARS;
    unsigned int highNFracs = 0, lowNFracs = INT_MAX;
    
    
    index = -1;
    
    //checking if some constraint can be branched
    
    for(int i = 0; i < nbinSumEqConstrs; i++)
    {
        const int c = binSumEqConstrs[i];
        
        unsigned int nFracs = 0;
        
        
        const int *acols = A[c];
        const unsigned int nel = A.getNumberOfElementsAtRow(c);
        
        for(unsigned int j = 0; j < nel; j++)
        {
            const int col = acols[j];//row[j].getColumn();
            
            if( MRQ_binGap(sol[col]) > intTol )
                nFracs++;
            
        }
        
        if( nFracs > 0 )
        {
            if( compLow )
            {
                if( nFracs < lowNFracs )
                {
                    lowNFracs = nFracs;
                    index = c;
                }
            }
            else //if( branchStrategy == MRQ_BB_CBS_HIGHEST_NUMBER_OF_VARS )
            {
                if( nFracs > highNFracs )
                {
                    highNFracs = nFracs;
                    index = c;
                }
            }
        }
    }
    
    
    /*#if MRQ_DEBUG_MODE
        assert( compLow ? lowNFracs < INT_MAX : highNFracs > 0 );
    #endif*/
    
    
    #if MRQ_DEBUG_MODE
        assert(index >= 0);
    #endif
    
    return index >= 0;
}




inline bool MRQ_BinSumConstrsChooser ::chooseIndexToBranchPseudoCost( const minlpproblem::MIP_SparseMatrix& A, const int nbinSumEqConstrs, const int* binSumEqConstrs, const double intTol, const double* sol, const int* reverseIntVars, const double pcost_mu, const MRQ_NewPseudoCost* pcost, int& index ) const
{
    double maxpc = -1.0;

    index = -1;

    for(int k = 0; k < nbinSumEqConstrs; k++)
    {
        const int c = binSumEqConstrs[k];
        
        double lowfr = INFINITY, higfr = -INFINITY, sumfr = 0.0;
        
        unsigned int ninds = 0;
        //minlpproblem::MIP_SparseRow &row = A[c];
        //const unsigned int nel = row.getNumberOfElements();
        
        const int *acols = A[c];
        const unsigned int nel = A.getNumberOfElementsAtRow(c);
        
        
        for(unsigned int j = 0; j < nel; j++)
        { 
            const int ind = acols[j]; //row[j].getColumn();
            
            if( MRQ_binGap(sol[ind]) > intTol )
            {
                const int rind = reverseIntVars[ind];
                
                //we have binary variable. So, the diference is 1.0 - sol on right
                const double fr = (pcost[rind].pr *( 1.0 - sol[ind] ))/pcost[rind].nPr ;
                
                if( fr < lowfr )
                    lowfr = fr;
                if( fr > higfr )
                    higfr = fr;
                
                sumfr += fr;
                ninds++;
                
            }
            
        }
        
        if( ninds > 0 )
        {
            const double pcost = (1.0 - pcost_mu)*lowfr + pcost_mu *higfr;
            
            //printf("constraint: %d pcost: %f\n", c, pcost);
            
            if(pcost > maxpc)
            {
                maxpc = pcost;
                index = c;
            }
        }
        
    }

    #if MRQ_DEBUG_MODE
        assert(index >= 0);
    #endif

    return index >= 0;
}



bool MRQ_BinSumConstrsChooser ::chooseIndexToBranch( const minlpproblem::MIP_SparseMatrix& A, const double intTol, const double* sol, const int branchStrategy, const int* reverseIntVars, const double pcost_mu, const MRQ_NewPseudoCost* pcost, int& index )
{
    bool code;

    #if MRQ_DEBUG_MODE
        //if( nCandInds <= 0  )
            //std::cout << "nCandInds: " << nCandInds << "\n";
        assert(nCandInds > 0);
    #endif

    if( branchStrategy == MRQ_BB_CBS_STBRANCH_PSEUDO_COSTS )
    {
        code = chooseIndexToBranchPseudoCost( A, nCandInds, candInds, intTol, sol, reverseIntVars, pcost_mu, pcost, index );
    }
    else
    {
        #if MRQ_DEBUG_MODE
            assert( branchStrategy == MRQ_BB_CBS_HIGHEST_NUMBER_OF_VARS || branchStrategy == MRQ_BB_CBS_LOWEST_NUMBER_OF_VARS );
        #endif
        
        code = chooseIndexToBranchLowHigh( A, nCandInds, candInds, intTol, sol, branchStrategy, index );
    }

    #if MRQ_DEBUG_MODE
        nCandInds = -1; //We must call calculateCandidateConstrsToBranch before call this method. So, we put -1 here to be sure calculateCandidateConstrsToBranch was called by the first assert
    #endif

    return code;
}




MRQ_IGMA2Iteration::MRQ_IGMA2Iteration()
{
    prob = nullptr;
    nI = UINT_MAX;
    nC = UINT_MAX;
    intVars = nullptr;
    contVars = nullptr;
    
    fixnlx = nullptr;
    fixnux = nullptr;
    
    resetParameters();
}


MRQ_IGMA2Iteration::~MRQ_IGMA2Iteration()
{
    MRQ_secFree(fixnlx);
}


void MRQ_IGMA2Iteration::resetParameters()
{
    in_set_max_dist_constr_on_bin_vars = false;
    in_solve_local_search_problem_even_on_non_int_sol = false;
    in_print_level = 4;
    in_neighborhood_strategy = MRQ_IGMA2_NS_RECTANGULAR;
    in_factor_to_max_dist_constr_on_bin_vars = 0.05;
    in_percentual_to_rectangular_neighborhood = 0.05;
    in_integer_tol = 0.001;
}


int MRQ_IGMA2Iteration::run(unsigned int thnumber, MRQ_GapMinProb &gapminsolver, MRQ_NLPSolver *nlp, const double *nlx, const double *nux, const int distConstIndex, const double maxDistance, const double *sol, const double *dualSolC, const double *dualSolV,  double *outputSol, double &objOutputSol, MRQ_Preprocessor *prepocessor, minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager, bool *out_feas_sol_on_gap_min_problem)
{
    const int n = prob->n;
    const int m = prob->m;
    
    MRQ_NLPSolver *gapmin = (MRQ_NLPSolver *) gapminsolver.solver;
    
    if(out_feas_sol_on_gap_min_problem)
        *out_feas_sol_on_gap_min_problem = false;
    
    
    //updating variables and constraint bounds...
    {
        int r;
        double clb, cub;
        
        r = gapmin->setnVariablesBounds(n, nlx, nux);
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
        }
        
        #if 1
        //seting bounds of continuous variables to consider maxDist (I hope that helps nlp solver)
        for(unsigned int i = 0; i < nC; i++)
        {
            const int ind = contVars[i];
            const double lb = MRQ_max(sol[ind] - maxDistance, nlx[ind]);
            const double ub = MRQ_min(sol[ind] + maxDistance, nux[ind]);
            
            int r = gapmin->setVariableBounds( ind, lb, ub);
            if(r != 0)
            {
                if(in_print_level > 0)
                    MRQ_PRINTERRORNUMBER(r);
            }
        }
        #endif
        
        
        
        #if 1
        if(nlp)
        {
            //preprocessor can have changed constraint bounds in this node
            for(int i = 0; i < m; i++)
            {
                r = nlp->getConstraintBounds(i, clb, cub);
                if(r == 0)
                {
                    r = gapmin->setConstraintBounds(i, clb, cub);
                    if(r != 0)
                    {
                        if(in_print_level > 0)
                            MRQ_PRINTERRORNUMBER(r);
                    }
                }
                else
                {
                    if(in_print_level > 0)
                        MRQ_PRINTERRORNUMBER(r);
                }
            }
        }
        #endif
        
    }
    
    
    if( in_neighborhood_strategy == MRQ_IGMA2_NS_SPHERIC )
    {
        optsolvers::OPT_SolutionDistanceConstraintSetter sdcs;
        //int auxInds[n];
        //double auxVars[n];
        //double *auxVars = nlp->sol; //using nlp.sol as aux array... note, we will overwrite nlp->sol below anyway if gapmin->solve gets success
        
        
        int r = sdcs.setDistConstraint(*gapmin, distConstIndex, nC, contVars, sol, maxDistance);
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        
        if( in_set_max_dist_constr_on_bin_vars )
        {
            const int distIntConstIndex = distConstIndex + 1;
            const double tolMaxIntDistance = 0.001;
            double maxIntDistance = 0.0;
            
            //maxDistance for integer vars is the total integer gap plus a tolerance plus a slack to change some variables values
            
            for(unsigned int i = 0; i < nI; i++)
                maxIntDistance += MRQ_gap( sol[intVars[i]] );
            
            maxIntDistance += 1.0 + ceil(in_factor_to_max_dist_constr_on_bin_vars *nI) + tolMaxIntDistance; //the tolerance is to due to possible numerical errors when we sum the values to calculate maxIntDistance
            
            
            int r = sdcs.setDistConstraint(*gapmin, distIntConstIndex, nI, intVars, sol, maxIntDistance);
            MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        }
        
        
        //if(gapmin->isMyNLPClass())
            //((OPT_MyNLPSolver *) gapmin)->prob.print(), MRQ_getchar();
    }
    else if( in_neighborhood_strategy == MRQ_IGMA2_NS_RECTANGULAR )
    {
        optsolvers::OPT_SolutionDistanceBoxConstraintsSetter sdbcs;
        
        int r = sdbcs.setDistConstraint(*gapmin, nlx, nux, nC, contVars, sol, in_percentual_to_rectangular_neighborhood);
        
        if(r != 0)
        {
            if(in_print_level > 0)
                MRQ_PRINTERRORNUMBER(r);
            
            return MRQ_NLP_SOLVER_ERROR;
        }
        
        /*for(int i = 0; i < nC; i++)
        {
            const int ind = contVars[i];
            std::cout << "sol["<<ind<<"]: " << sol[ind] << " nlx["<<ind<<"]: " << nlx[ind] << " nux["<<ind<<"]: " << nux[ind] << "\n";
        }
        
        gapmin->generateModelFile("model.mod");
        std::cout << "Gerei model.mod\n";
        MRQ_getchar(); */
    }
    else
    {
        assert(false);
    }
    
    
    
    gapmin->setInitialSolution(sol, NULL, dualSolV); //we cannot use constraints dual solution because gapmin has one additional constraint
    
    gapmin->solve(false);
    
    //std::cout << "igma - gap min solver retrun code: " << gapmin->retCode <<  " obj: " << gapmin->objValue << " orig ret code: " << gapmin->origSolverRetCode << "\n";
    //for(int i = 0; i < n; i++)
        //std::cout << "gap min sol["<<i<<"]: " << gapmin->sol[i] << " \t"  << std::endl;
    
    //MRQ_getchar();
    
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
    
    
    
    if(gapmin->feasSol && nlp)
    {
        const bool gapMinIntSol = MRQ_isIntegerSol(nI, intVars, gapmin->sol, in_integer_tol);
        
        
        if( gapMinIntSol || in_solve_local_search_problem_even_on_non_int_sol )
        {
            double obj;
            double *psol = NULL;
            //const double *plx, *pux;
            int r;
            
            
            //we prepocess to try detect infeasibility. We only preprocessif the solution from Gap Min Problem is not feasible
            if(prepocessor && !gapMinIntSol )
            {
                bool varBoundsUpdt, constrBoundsUpdt;
                double *gsol = gapmin->sol;
                
                if( !fixnlx )
                {
                    MRQ_malloc(fixnlx, 2*n);
                    MRQ_IFMEMERRORRETURN(!fixnlx); //we only deallocate in the destructor
                    
                    fixnux = &fixnlx[n];
                }
                
                
                MRQ_copyArray(n, nlx, fixnlx);
                MRQ_copyArray(n, nux, fixnux);
                
                
                for(decltype(nI) i = 0; i < nI; i++)
                {
                    auto ind = intVars[i];
                    const double fix = round( gsol[ind] );
                    
                    fixnlx[ind] = fix;
                    fixnux[ind] = fix;
                }
                
                
                
                //we only preprocess linear constraints to try detect infeasibility. 
                
                int r = prepocessor->preprocess(nI, intVars, *ccstorager, false, false, INFINITY, fixnlx, fixnux, varBoundsUpdt, constrBoundsUpdt);
                
                if( r == minlpproblem::MIP_INFEASIBILITY )
                {
                    return MRQ_HEURISTIC_FAIL;
                }
                
                MRQ_IFERRORRETURN(r, MRQ_UNDEFINED_ERROR);
                
                
                /*plx = fixnlx;
                pux = fixnux;*/
                
                nlp->setnVariablesBounds(n, fixnlx, fixnux); //integer vars were already rounded
            }
            else
            {
                /*plx = nlx;
                pux = nlx;*/
                
                nlp->setnVariablesBounds(n, nlx, nux);
                
                r = MRQ_fixIntVarsOnSolByList(nI, intVars, gapmin->sol, *nlp);
                if(r != 0)
                {
                    MRQ_PRINTERRORNUMBER(r);
                    return n;
                }
                
            }
                
            
            /*for(int i = 0; i < m; i++)
            {
                double clb, cub;
                
                int r = nlpSolver.getConstraintBounds(i, clb, cub);
                if(r == 0)
                {
                    nlp->setConstraintBounds(i, clb, cub);
                    if(r != 0)
                    {
                        if(in_print_level > 0)
                            MRQ_PRINTERRORNUMBER(r);
                    }
                }
                else
                {
                    if(in_print_level > 0)
                        MRQ_PRINTERRORNUMBER(r);
                }
                
            } */
            
            
            nlp->setInitialSolution( gapmin->sol, gapmin->dualSolC, gapmin->dualSolV );
            
            nlp->solve(false, true, false, false);
            
            r = MRQ_unfixIntegerVarsByList(nI, intVars, nlx, nux, *nlp);
            if(r != 0)
            {
                MRQ_PRINTERRORNUMBER(r);
                return r;
            }
            
            //std::cout << "\t\tnlp - feasSol: " << nlp->feasSol << " obj: " << nlp->objValue << "\n";
            //for(int i = 0; i < n; i++)
                //std::cout << "\t\tlocal search sol["<<i<<"]: " << nlp->sol[i] << " \t"  << std::endl;
            
            
            
            /*{
                bool feas;
                double objValue;
                
                prob->objEval(thnumber, true, nlp->sol, objValue);
                
                prob->isFeasibleToConstraints(thnumber, nlp->sol, true, nullptr, 1e-5, 1e-5, feas);
            }
            
            MRQ_getchar(); */
            
            
            //std::cout << "\t\t\tBusca local. retCode: " << nlp->retCode << " feasSol: " << nlp->feasSol << " obj: " << nlp->objValue << "\n";
            
            if( nlp->feasSol )
            {
                obj = nlp->objValue;
                psol = nlp->sol;
            }
            else if(gapMinIntSol)
            {
                psol = gapmin->sol;
                
                int r = prob->objEval(thnumber, true, psol, obj);
                
                if(r != 0)
                {
                    psol = NULL;
                    if(in_print_level > 0)
                        MRQ_PRINTCALLBACKERRORNUMBER(r);
                }
                
            }
            
            
            if(psol)
            {
                MRQ_copyArray(n, psol, outputSol);
                
                objOutputSol = obj;
                
                if(out_feas_sol_on_gap_min_problem)
                    *out_feas_sol_on_gap_min_problem = gapMinIntSol;
                
                //MRQ_getchar();
                return MRQ_HEURISTIC_SUCCESS;
            }
            
            
            #if 0 
                coloque isso na igma3
                if(psol)
                {
                    bool updt = tryUpadteBestSolution(thnumber, n, psol, obj, iter);
                    
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
            #endif
            
        }
        
    }
    
    
    return MRQ_HEURISTIC_FAIL;

}







