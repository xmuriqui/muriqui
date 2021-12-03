/*
* Functions to implement new branch-and-bound procedures at Muriqui
* 
* 
* 
*/



#ifndef NEW_MRQ_BB_HPP_
#define NEW_MRQ_BB_HPP_

#include <iostream>
#include <list>

#include "BBL_branchAndBound.hpp"
//#include "OPT_solvers.hpp"
#include "MRQ_config.hpp"
#include "MRQ_algClasses.hpp"
#include "MRQ_solvers.hpp"

#include "MRQ_constants.hpp"
#include "MRQ_dataStructures.hpp"
#include "MRQ_ssrounding.hpp"


//#define MRQ_CPP_MULTITHREADING BBL_CPP_MULTITHREADING

//#define MRQ_OMP_MULTITHREADING BBL_OMP_MULTITHREADING



namespace muriqui
{


//flag to use bounds from Outer Approximation as official lower bounds in bb node
#define MRQ_USE_OA_BOUND_AS_OFICIAL_LB_IN_BBNODE 0


#define MRQ_DEBUG_IGMA2_BUG 0



typedef branchAndBound::BBL_NodeBounds 	MRQ_NodeBounds;

typedef branchAndBound::BBL_Array<double>	MRQ_DoubleArray;

typedef branchAndBound::BBL_Mutex MRQ_Mutex;

typedef branchAndBound::BBL_NodeBoundsSol MRQ_NodeBoundsSol;



//check if the value parentNodeBoundsStrategy is appropriate to problem being addressed. If parentNodeBoundsStrategy indicates some float type, we have to make sure it all vaiable bounds can be stored in a 32 bits float. Otherwise, we have to change to a double type. 
void MRQ_checkParentNodeBoundsStrategy(const int nI, const int* intVars, const double *lx, const double *ux, branchAndBound::BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY &parentNodeBoundsStrategy);








/*class MRQ_NewBBNode : public branchAndBound::BBL_NodeBase
{
    
public:
    
    
    double heurlb; //heuristic value to node. We set here the outer approximation lower bound, for example... That value is not take in account to sort the nodes in best limit strategy...
    
    branchAndBound::BBL_ArraySize <branchAndBound::BBL_NodeBounds> *parentBounds;
    
    unsigned int nMyBounds;
    MRQ_IndexBranch *myBounds;
    
    
    MRQ_NewBBNode();
    
    virtual ~MRQ_NewBBNode();
    
    virtual void print( std::ostream &out = std::cout );
    
    void setDadNodeBoundsPointer( branchAndBound::BBL_ArraySize< branchAndBound::BBL_NodeBounds >* p );
    
    void subDesallocate();
    
    //that functions set olny positions correspondets to new bounds
    void getVarBoundsOnNode(double *nlx, double *nux);
    
    //void restoreOriginalVarBounds(const double *lx, const double *ux, double *nlx, double *nux);
};*/





BBL_PRE_PACK class MRQ_NewBBNode : public branchAndBound::BBL_Node
{
    
public:
    
    #if MRQ_SET_LBHEUR_ON_BB_NODE
        double heurlb; //heuristic value to node. We set here the outer approximation lower bound, for example... That value is not take in account to sort the nodes in best limit strategy...
    #endif
    
    #if MRQ_DEBUG_IGMA2_BUG
        unsigned int igma2OnAncestral; //if this variable is greather than zero, igma2 was aplied in some ancestral of the current node. We add 1 in this variable in each level
    #endif
    
    
    //unsigned int nBranchVars; //number of variables used to original branching. This value can be different from nMyBounds because we can run some heuristics to add new bounds to other variables...
    
    MRQ_NewBBNode( branchAndBound:: BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY parentBoundsStrategy, const unsigned int maxVars);
    
    #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
    #endif 
    ~MRQ_NewBBNode();
    
    
    #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
    #endif 
    double getBestLowerBound() const;
    
    
    #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
    #endif 
    void print(std::ostream &out = std::cout) const 
    #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        override
    #endif 
    ;
    
    //virtual double getLowerBound() const override;
}BBL_POS_PACK ;






class MRQ_NewPseudoCost
{
        
public:
    
    unsigned int nrealPl, nrealPr; //number of times that real strong branching was computed for left and right.
    unsigned long int nPl, nPr; //total number of terms computed in pl and pr
    long double pl, pr; //pl -> left branch     pr ->right branch (here we have the absolute sum, not the mean...)
    double minPl, minPr; //lowest left and right pseudocost in the mean...
    
    MRQ_NewPseudoCost()
    {
        nrealPl = nrealPr = 0;
        nPl = nPr = 0;
        pl = pr = 0.0;
        minPl = minPr = INFINITY;
    }
    
    
};



/*Class to alow user generate branch-and-bound nodes and so, performs a customized branch using callback object. Users can creade nodes passing bounds using one of the methods in the class */

class MRQ_NewUserNodeGenerator2
{
    MRQ_NewBBNode *parent;
    branchAndBound::BBL_UserNodeGenerator *bblUserNodeGen;
    branchAndBound::BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY parentBoundsStrategy;
    
    void initialize(MRQ_NewBBNode *parentNode, branchAndBound::BBL_UserNodeGenerator *userNodeGen);
    
    int n;
    
public:
    
    
    MRQ_NewUserNodeGenerator2();
    
    
    void init();
    
    
    /*
    * create a BB node. That node will be added to open nodes list. nodelx and nodeux describes variable bounds for the new node
    * 
    * newNode pointer can points to an object where Muriqui will represent the node. It can be useful if you are using your own class derived from MRQ_BBNode to represent nodes. If you pass a NULL pointer, Muriqui will generate that object using MRQ_BBNode class. Even if you generate newNode by yourself, Muriqui will be responsable by delete the newNode in the appropriate moment.
    * 
    * You can use nodeLowerBound if you want calculate a new lower bound to node. BB procedure will use as node lower bound the maximum between that value and the bound gotten from node's parent (curretn node beinx exploited). So, if you do not want to calculate, you can pass a huge negative value as -MRQ_INFINITY.  
    * 
    * If either initSol or initDual is NULL, the new  node will use the optimal solution from its parent to build its initial solution. Note if you wish change it, you need provide both primal and dual solution.
    * 
    * data will be copy inside that procedure to generate the BB node. So, you can free memory pointed by 
    * MRQ_UserNodeData *data after call this fucntion.
    * 
    */
    
    
    int generateNode( const double *nodelx, const double *nodeux, MRQ_NewBBNode *newNode = NULL, const double nodeLowerBound = -MRQ_INFINITY, const double *initSol = NULL, const double *initDual = NULL );
    
    
    
    /*
    * create a BB node. That node will be added to open nodes list.
    * 
    * If flag inheritParentBounds is true, the new node inherits bounds from its parent (the current node being exploited) automatically beyond its own bounds. Note bounds described in newBounds can overwrites bounds from parents when both reference a same subset of variables.
    * 
    * newNode pointer can points to an object where Muriqui will represent the node. It can be useful if you are using your own class derived from MRQ_BBNode to represent nodes. If you pass a NULL pointer, Muriqui will generate that object using MRQ_BBNode class. Even if you generate newNode by yourself, Muriqui will be responsoble by delete the newNode in the appropriate moment.
    * 
    * You can use nodeLowerBound if you want calculate a new lower bound to node. BB procedure will use as node lower bound the maximum between that value and the bound gotten from node's parent (curretn node beinx exploited). So, if you do not want to calculate, you can pass a huge negative value as -MRQ_INFINITY.  
    * 
    * If either initSol or initDual is NULL, the new  node will use the optimal solution from its parent to build its initial solution. Note if you wish change it, you need provide both primal and dual solution.
    * 
    * data will be copy inside that procedure to generate the BB node. So, you can free memory pointed by 
    * MRQ_UserNodeData *data after call this fucntion.
    * 
    */
    
    int generateNode( const bool inheritParentBounds, const unsigned int nNewBounds, const muriqui::MRQ_NodeBoundsSol* newBounds, const bool isNewBoundsAscOrdered = false, muriqui::MRQ_NewBBNode* newNode = 0, const double nodeLowerBound = -1.0e20, const double* initSol = 0, const double* initDual = 0);
    
    
friend class MRQ_BBLCallbacks;
};




/* During Pseudo-cost calculations, we initialize using strong branching. In that strong branching procedure, we simulate branching for all fractional vars. So we use that class to store the lower bounds on this ramifications. Our goal is use that information to perform future prunes when upper bounds are gotten.
    */
class MRQ_NewFirstBranch
{
public:
    
    double bPoint; //the point choosem to perform branhc
    double leftlb; //lower bound to branch [lx bPoint]
    double rightlb; //upper bound to branch [bPoint+1 ux] on case of integer vars. Maybe in the future we can use that strategy o continuous variables also
    
    MRQ_NewFirstBranch();
    
    void initialize();
};


/*a class to storage pseudo costs*/
class MRQ_BasePseudoCostCalc
{
    
private:
    
    int getObjIncreaseEstimativeByVar (const unsigned int indexInIntVar, const unsigned int minNumberOfPCostsCalculations, const double oldVarValue, const double newVarValue,  const double alpha, double &estimative);
    
public:
    
    MRQ_NewPseudoCost *pcost;
    MRQ_Mutex SEMAPH_sem;
    
    
    MRQ_BasePseudoCostCalc();
    
    virtual ~MRQ_BasePseudoCostCalc();
    
    virtual int allocate(const int nI);
    
    virtual void deallocate();
    
    virtual void initialize();
    
    void updatePCosts( const unsigned int nThreads, const int *reverseIntVars, const unsigned int sizeBounds, const branchAndBound::BBL_ClassUnionNodeBoundsSolPointer &bounds, const double fDad, const double fNode, const double* x );
    
    int getObjIncreaseEstimative(const bool  onlyApplyOnFixedIntVars, const int *reverseIntVars, const unsigned int minNumberOfPCostsCalculations, const unsigned int sizeBounds, const branchAndBound::BBL_ClassUnionNodeBoundsSolPointer &bounds, const double alpha, double &estimative, bool &allVarIncreasingEstimated);
    
    bool checkIfNodeCanBePrunedByEstimative( const bool  onlyApplyOnFixedIntVars, const int nI, const int *intVars, const unsigned int minNumberOfPCostsCalculations, const double *nlx, const double *nux, const double relaxObj, const double *relaxSol, const double intTol, const double alpha, const double zu );
    
    void print( const int nI, const int *intVars, std::ostream &out = std::cout);
};



class MRQ_NewPseudoCostCalc : public MRQ_BasePseudoCostCalc
{
    //totalThreads is the total number of threads in our B&B. nThreads is the number of threads being used to calculate strong branching...
    int threadStrongBranching( const unsigned int totalThreads, const unsigned int nThreads, const unsigned int thnumber, MRQ_Mutex& SEMAPH_indSem, MRQ_Mutex& SEMAPH_solSem, MRQ_Mutex &SEMAPH_nodeSem, int& nextInd, MRQ_MINLPProb &prob, const int nI, const int* intVars, const bool storeFirstBranch, const bool considInfeasIfNlpFail, const int printLevel, const double intTol, const double* olx, const double* oux, MRQ_NLPSolver* nlp, MRQ_NewBBNode* node, double* nlx, double* nux, const double* nodeSol, const double objNodeSol, const bool calculateEvenIfMaxComptSBranchReached, bool &newBoundsAllocated, bool& updtBounds, bool& canFathom, double& zu, bool& intSolFound, double* intSol, int& retCode );
    
    
    
public:
    
//vars to pseudocost
    bool allPCostInit;
    unsigned int maxComptSBranch;
    
    //int *intVars;
    //int *reverseIntVars; //reverse indices in intVars. That vector has the position where a variable appear in intVars
    
    MRQ_NewFirstBranch *fbranch; //bounds on first branching
    double highestfbranch;
    //MRQ_MINLPProb *prob;
    
    
    
    MRQ_NewPseudoCostCalc();
    
    virtual ~MRQ_NewPseudoCostCalc();
    
    virtual void initialize() override;
    
    virtual int allocate(const int nI) override;
    
    
    //totalThreads is the total number of threads in our B&B. nThreads is the number of threads being used to calculate strong branching...
    int calculateStrongBranching( const unsigned int totalThreads, const unsigned int nThreads, MRQ_MINLPProb &prob, const int nI, const int* intVars, const bool storeFirstBranch, const bool considInfeasIfNlpFail, const int printLevel, const double intTol, const double* olx, const double* oux, MRQ_NLPSolver** nlps, MRQ_NewBBNode* node, double* nlx, double* nux, const double* nodeSol, const double objNodeSol, const bool calculateEvenIfMaxComptSBranchReached, bool& updtBounds, bool& canFathom, double& zu, bool& intSolFound, double* intSol );
    
    
    bool checkIfAllPCostsAreCalculated( const int nI, const int *intVars, const double* olx, const double* oux);
    
    
    virtual void deallocate() override;
    
    
    
    //return true if the node can be prunned
    bool updateVarBoundsByUpperBound( const int nI, const int *intVars, const double zu, double *nlx, double *nux );
    
    
    
    friend int MRQ_strongBranching(MRQ_NewPseudoCostCalc *pc, const unsigned int totalThreads, const unsigned int nThreads, const unsigned int thnumber, MRQ_Mutex *SEMAPH_indSem, MRQ_Mutex *SEMAPH_solSem, MRQ_Mutex *SEMAPH_nodeSem, int *nextInd, MRQ_MINLPProb *prob, const int nI, const int *intVars, const bool storeFirstBranch, const bool considInfeasIfNlpFail, const int printLevel, const double intTol, const double *olx, const double *oux, MRQ_NLPSolver *nlp, MRQ_NewBBNode *node, double *nlx, double *nux, const double *nodeSol, const double objNodeSol, const bool calculateEvenIfMaxComptSBranchReached, bool *newBoundsAllocated, bool *updtBounds, bool *canFathom, double *zu, bool *intSolFound, double *intSol, int *retCode ); 
    
};



int MRQ_strongBranching( MRQ_NewPseudoCostCalc *pc, const unsigned int totalThreads, const unsigned int nThreads, const unsigned int thnumber, MRQ_Mutex *SEMAPH_indSem, MRQ_Mutex *SEMAPH_solSem, MRQ_Mutex *SEMAPH_nodeSem, int *nextInd, MRQ_MINLPProb *prob, const int nI, const int *intVars, const bool storeFirstBranch, const bool considInfeasIfNlpFail, const int printLevel, const double intTol, const double *olx, const double *oux, MRQ_NLPSolver *nlp, MRQ_NewBBNode *node, double *nlx, double *nux, const double *nodeSol, const double objNodeSol, const bool calculateEvenIfMaxComptSBranchReached, bool *newBoundsAllocated, bool *updtBounds, bool *canFathom, double *zu, bool *intSolFound, double *intSol, int *retCode ); 





class MRQ_NewChooseIndexToBranch
{
    
protected:
    
    //void inline getHighestGaps(const int maxGaps, const int* inds, const int sizeInds, const double* gap, const double *gap2, int& ninds, unsigned int* choosenInds);
    
    void getHighestGaps( const int maxInds, const double intTol, const int* inds, const int sizeInds, const double *sol, double *auxVars, int& ninds, unsigned int* choosenInds );
    
    void getHighestPcosts( const int maxInds, const double intTol, const int nI, const int* intVars, const int *reverseIntVars, const double *sol, const double* nlx, const double* nux, const double pcost_mu, const MRQ_NewPseudoCost *pcost, double *auxVars, int& ninds, unsigned int* choosenInds );
    
    void getHighestPcostsEvenIntegerValues( const int maxInds, const int nI, const int* intVars, const int *reverseIntVars, const double *sol, const double* nlx, const double* nux, const double pcost_mu, const MRQ_NewPseudoCost *pcost, double *auxVars, int& ninds, unsigned int* choosenInds );
    
    void getHighestPriorities( const int maxInds, const double intTol, const int* inds, const int sizeInds, const double *sol, const int *priorities, double *auxVars, int& ninds, unsigned int* choosenInds );
    
    
public:
    
    bool *flagVars;
    //MRQ_MINLPProb *prob;
    //MRQ_NewPseudoCostCalc *pCosts;
    
    //int nI; //number of integer variables
    //int nbin; //number of binary varaibles...
    //int strategy;
    //int *intVars, *binVars, *nonBinVars;
    //double intTol;
    //double *gap, *pcostEst;
    
    
    
    
    
    MRQ_NewChooseIndexToBranch();
    
    
    ~MRQ_NewChooseIndexToBranch();
    
    
    int allocate(const int maxBranchVars, const int n);
    
    
    void desallocate();
    
    
    void chooseIndices( const MRQ_NewPseudoCost *pcost, const int strategy, const double intTol, const int nI, const int *intVars, const int *reverseIntVars, const int nbin, const int *binVars, const int *nonBinVars, const double pcost_mu, const int maxInds, const double* relaxSol, const bool nlpSucess, const double* nlx, const double* nux, const int *priorities, double *auxVars, int& ninds, unsigned int* choosenInds, bool &branchOnInt );
    
};



//class to save points to oa linearization...
class MRQ_NewPoints
{
    
    int nPointsAllocated; //number of positions allocated to first array (*points) 
    
    
    void initialize(const int nPreAlloc);
    
    
    
public:
    
    int nPreAlloc; //number of points that will be preallocated in a preAllocate call...
    int nPoints;
    double **points;
    
    
    
    MRQ_NewPoints(const int nPreAlloc);
    
    ~MRQ_NewPoints();
    
    
    int preAllocate(const int *numberOfPreAlloc = NULL);
    
    //warning: that function returns the number of points allocated...
    int addPoints(const int nPoints, const int dim, double **x);
    
    void desallocate();
    
};




//class to store a quadratic cut in a temporary way:
//  lb <=  a'x + 0.5x'Qx <= ub
// only put the lower triangle of Q
class MRQ_NewQuadraticCut
{
    int allocateLinearPart(const int nz);
    
    int allocateQuadraticPart(const int nz);
    
public:
    double lb, ub;
    
    int anz; //number of nonzeros in a
    int *acols; //indices to nonzeros in a
    double *avals; // nonzero values in a
    
    int Qnz; //number of nonzeros in Q
    int *Qrows; //rows indices of nonzeros in Q (triple sparse format)
    int *Qcols; //cols indices of nonzeros in Q (triple sparse format)
    double *Qvals; //vals indices of nonzeros in Q (triple sparse format)
    
    
    MRQ_NewQuadraticCut();
    
    ~MRQ_NewQuadraticCut();
    
    int addCutOnSolver(MRQ_NLPSolver *nlp);
    
    void desallocateLinearPart();
    
    void desallocateQuadraticPart();
    
    void setBounds(const double lb, const double ub);
    
    //warning: that method erase all line part coeficients already stored
    int setLinearPart(const int nz, const int* cols, const double* vals);
    
    int setQuadraticPart(const int nz, const int* rows, const int* cols, const double* vals);
};




class MRQ_NewQuadCutStorer
{
    int remainThreads; //A counter to Threads that do not still set this cut
    unsigned int cutId; //An unique (sequential) number to identify the cutId. (By now, it will be not used)
    
public:
    
    MRQ_NewQuadraticCut cut;
    
    //the first argument is the number of threads that will set the cut
    MRQ_NewQuadCutStorer(int nThreads, unsigned int cutId);
    
    int decThreadCounter();
    
    // that function set the cut:
    //  lb <=  a'x + 0.5x'Qx <= ub
    // only put the lower triangle of Q
    int setQuadraticCut(const int nza, const int *acols, const double *avalues, const int nzQ, const int *Qrows, const int* Qcols, const double *Qvalues, const double lb, const double ub);
    
    
    friend class MRQ_NewGlobalCutList;
};




class MRQ_NewGlobalCutGenerator;
    
//each thread will have an object like that to store the cuts...
class MRQ_NewGlobalCutList
{
    MRQ_Mutex SEMAPH;
    
    
public:
    
    bool isempty;
    int nThreads;
    MRQ_NewGlobalCutGenerator *cutGen;
    std::list< MRQ_NewQuadCutStorer* > cutList;
    
    MRQ_NewGlobalCutList( );
    
    MRQ_NewGlobalCutList(const int numberOfThreads);
    
    ~MRQ_NewGlobalCutList();
    
    int addCutsOnSolver( MRQ_NLPSolver *nlp );
    
    //void deleteMutex();
    
    void eraseList();
    
    int initialize(const int numberOfThreads);
    
    int insertCutStorer( MRQ_NewQuadCutStorer *cutStorer );
};




class MRQ_NewGlobalCutGenerator
{
    unsigned int nextId; //number used to set ID's to cuts
    
    MRQ_Mutex SEMAPH; 
    
public:
    
    int nThreads;
    MRQ_NewGlobalCutList *cutLists;
    
    
    
    MRQ_NewGlobalCutGenerator(const int nThreads, MRQ_NewGlobalCutList *cutLists);
    
    ~MRQ_NewGlobalCutGenerator();
    
    int decThreadCounter(MRQ_NewQuadCutStorer *cutStorer);
    
    //void deleteMutex();
    
    //int initializeMutex(const int nThreads);
    
    // that function set the cut:
    //  lb <=  a'x + 0.5x'Qx <= ub
    //only put the lower triangle of Q
    int setQuadraticCut(const int nza, const int *acols, const double *avalues, const int nzQ, const int *Qrows, const int* Qcols, const double *Qvalues, const double lb, const double ub);
    
    
};






/* Class to store indices of constraints. Here, we are interested in constraints in the form:
* y_{k_1} + y{k_2} + ... + y{k_n} = b   (1)
* or
* y_{k_1} + y{k_2} + ... + y{k_n} <= b  (2)
*
* where all y variables are binaries and b >= 1.
*
* this class will generate a special branching based on a constraints in one of these classes.
*
* For a constraint in the case (1),  if b== 1, we generate at most n new nodes. In each one of them, we have one diferent variable fixed at 1 (we assume there is no variable fixed at 1 in this constraint, off course). 
*
* For a constraint in the case (2), if b == 1, we generate the set of nodes for the case (1) plus an aditional node fixing all nonfixed variables at zero.
* 
* For constraints having b > 1, we subtrate b from the number of binary variables fixed at 1. If this value is one, we proceed to branch using the rules above.
* 
* UPDATED: we have updated to detect the following constarints also
* 
* y_{k_1} + y{k_2} + ... + y{k_n} + sum_i{ a_{j_i}*w{j_i} }  = b  (3)
* or
* y_{k_1} + y{k_2} + ... + y{k_n} + sum_i{ a_{j_i}*w{j_i} } <= b  (4)
* 
* 
* Where all y variable are binaries and all w variable are integer. 
* 
* For a constraint in the case (3), we can perform a constraint branching if all w variable are fixed and the result rhs: b - sum_i{ a_{j_i}*w{j_i} } = 1. In this case, like case (1) above, we generate at most n new nodes. In each one of them, we have a diferent y variable fixed at 1.
* 
* For a constraint in the case (4), we we can perform a constraint branching if all w variable are fixed and the result rhs: b - sum_i{ a_{j_i}*w{j_i} } = 1. In this case, like case (2), we ganerate the set of nodes for the case (3) plus an aditional node fixing all nonfixed y variables at zero.
* 
*/


class MRQ_BinSumConstrsInds
{
protected:
    
    int reallocate(unsigned int size, int*& inds );
    
    int _calculateIndices(const MRQ_MINLPProb &prob, const double *lx, const double *ux, const double *lc, const double *uc );
    
public:
    
    int nbinSumConstrs; //number of binary sum equality or inequality (<=) constraints having rhs > 1. 
    int *binSumConstrs;
    
    
    
    MRQ_BinSumConstrsInds();
    
    ~MRQ_BinSumConstrsInds();
    
    
    int calculateIndices(const MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL, const double *lc = NULL, const double *uc = NULL );
    
    void deallocate();
    
    void initialize();
    
};



#if 0 //deprecated
/* Class to store indices of constraints. Here, we are interested in constraints in the form:
* y_{k_1} + y{k_2} + ... + y{k_n} = b   (1)
* or
* y_{k_1} + y{k_2} + ... + y{k_n} <= b  (2)
* or
* y_{k_1} + y{k_2} + ... + y_{k_n-1} - y{k_n} =  0  (3) 
*
* where all variables are binaries and b >= 1.
*
* this class will generate a special branching based on a constraints in one of these classes.
*
* For a constraint in the case (1), if b== 1, we generate at most n new nodes. In each one of them, we have a diferent variable fixed at 1. 
*
* For a constraint in the case (2), if b == 1, we ganerate the set of nodes for the case (1) plus an aditional node fixing all variables at zero.
* 
* For a constraint in the case (3), we generate a node fixing all variables in zero (note it enough fix only y{k_n} at zero) and n-1 nodes fixing y_{k_n} and each variable y_{k_i}, i = 1,...,n-1 at 1 (note we do not need fix y_{k_n} at 1 in those nodes).
* 
* For constraints having b > 1, we subtrate b from the number of binary variables fixed at 1. If this value is one, we proceed to branch using the rules above.
* 
*/


class MRQ_BinSumConstrsInds2 : public MRQ_BinSumConstrsInds
{
    virtual int _calculateIndices(const MRQ_MINLPProb &prob, const double *lx, const double *ux, const double *lc, const double *uc ) override;
    
public:
    
    MRQ_BinSumConstrsInds2();
    
    virtual ~MRQ_BinSumConstrsInds2();
    
    
};

#endif







/* Class to choose a index of constraint to perform branching. Here, we are interested in constraints in the form:
* y_{k_1} + y{k_2} + ... + y{k_n} = b   (1)
* or
* y_{k_1} + y{k_2} + ... + y{k_n} <= b  (2)
* or
* y_{k_1} + y{k_2} + ... + y_{k_n-1} - y{k_n} =  0  (3) 
*
* where all variables are binaries and b >= 1.
*
* this class will generate a special branching based on a constraints in one of these classes.
*
* For a constraint in the case (1), if b== 1, we generate at most n new nodes. In each one of them, we have a diferent variable fixed at 1. 
*
* For a constraint in the case (2), if b == 1, we ganerate the set of nodes for the case (1) plus an aditional node fixing all variables at zero.
* 
* For a constraint in the case (3), we generate a node fixing all variables in zero (note it enough fix only y{k_n} at zero) and n-1 nodes fixing y_{k_n} and each variable y_{k_i}, i = 1,...,n-1 at 1 (note we do not need fix y_{k_n} at 1 in those nodes).
* 
* For constraints having b > 1, we subtrate b from the number of binary variables fixed at 1. If this value is one, we proceed to branch using the rules above.
* 
*/
class MRQ_BinSumConstrsChooser
{
    
    inline bool chooseIndexToBranchLowHigh(const minlpproblem::MIP_SparseMatrix &A, const int nbinSumEqConstrs, const int *binSumEqConstrs, const double intTol, const double *sol, const int branchStrategy, int &index ) const;
    
    inline bool chooseIndexToBranchPseudoCost(const minlpproblem::MIP_SparseMatrix &A, const int nbinSumEqConstrs, const int *binSumEqConstrs, const double intTol, const double *sol, const int *reverseIntVars, const double pcost_mu, const MRQ_NewPseudoCost *pcost, int &index ) const;
    
    
    
public:
    
    
    
    int nCandInds; //number of binary sum equality or inequality (<=) canidates to be branched
    int *candInds; //indices of candidates constraints to be branches
    
    
    
    
    MRQ_BinSumConstrsChooser();
    
    ~MRQ_BinSumConstrsChooser();
    
    
    void desallocate();
    
    void initialize();
    
    
    int reallocate(unsigned int size);
    
    
    //calculate candidate constraint indices to be branched in the current solution. This method returns the number of candidate indices
    int calculateCandidateConstrsToBranch( const int nbinSumConstrs, const int* binSumConstrs, const double* lx, const double* ux, const int *varType, const minlpproblem::MIP_SparseMatrix& A, const double* lc, const double* uc, const double intTol, const double* sol);
    
    
    //chose an index of constraint to perform branch... before call this method, it is necessary call calculateCandidateConstrsToBranch
    bool chooseIndexToBranch(const minlpproblem::MIP_SparseMatrix &A, const double intTol, const double *sol, const int branchStrategy, const int *reverseIntVars, const double pcost_mu, const MRQ_NewPseudoCost *pcost, int &index);
    
    
};






/* 
 * Class to perform random rounding based on constraints over
 * binary variables.
 * 
 */ 
class MRQ_ConstraintRouding
{
};




static inline branchAndBound::BBL_EXT_EXP_STRATEGY MRQ_MRQ_exp_strategy2BBL_exp_strategy( const int value )
{
    switch(value)
    {
        case MRQ_BB_ES_DEPTH_BEST_LIMIT:
            return branchAndBound::BBL_EES_DEPTH_BEST_LIMIT;
            
        case MRQ_BB_ES_BEST_LIMIT:
            return branchAndBound::BBL_EES_BEST_LIMIT;
            
        case MRQ_BB_ES_DEPTH:
            return branchAndBound::BBL_EES_DEPTH;
            
        case MRQ_BB_ES_WIDTH:
            return branchAndBound::BBL_EES_WIDTH;
            
        default:
            return branchAndBound::BBL_EES_DEPTH_BEST_LIMIT;
    }
    
}











//class to run a iteration of IGMA3. We separate it to use as a separated procedure IGMA3 and in other algorithms like BB
class MRQ_IGMA2Iteration
{
    
    double *fixnlx, *fixnux; //auxiliary arrays to 
    
public:
    
    MRQ_MINLPProb *prob;
    unsigned int nI, nC;
    const int *intVars;
    const int *contVars;
    
    
    bool in_set_max_dist_constr_on_bin_vars;
    bool in_solve_local_search_problem_even_on_non_int_sol;
    int in_print_level;
    MRQ_IGMA2_NEIGHBORHOOD_STRATEGY in_neighborhood_strategy;
    double in_factor_to_max_dist_constr_on_bin_vars;
    double in_percentual_to_rectangular_neighborhood;
    double in_integer_tol;
    
    
    MRQ_IGMA2Iteration();
    
    ~MRQ_IGMA2Iteration();
    
    void resetParameters();
    
    
    //if prepocessor is not null preprocessor will be used before solve nlp local search problem. You only need pass ccstorager if you pass prepocessor
    int run(unsigned int thnumber, MRQ_GapMinProb &gapminsolver, MRQ_NLPSolver *nlp, const double *nlx, const double *nux, const int distConstIndex, const double maxDistance, const double *sol, const double *dualSolC, const double *dualSolV,  double *outputSol, double &objOutputSol, MRQ_Preprocessor *prepocessor, minlpproblem::MIP_ConstraintsByColumnsStorager *ccstorager, bool *out_feas_sol_on_gap_min_problem = NULL);


};















class MRQ_BBLCallbacks : public branchAndBound::BBL_UserCallbacks
{
    long unsigned int iterNextOAApplic;


    int allocateThreadStructures(const int nthreads);

    int linearBoundsUpdating( MRQ_BB_BOUND_LIN_UPDT_STRATEGY strategy, const double zu, MRQ_NewBBNode &node, double *nlx, double *nux, bool *auxVars, MRQ_LPboundsUpdater &lpBound, MRQ_Random &random);

    //return true if solution updates the best solution
    bool tryUpdateBestSolution(const int threadNumber, double* solution, const double fsolution, const long unsigned int iter, const bool addSolToOALin);
    
    
    MRQ_Preprocessor* getPreprocessorPointer( unsigned int threadNumber);


    public:

    void deallocate();

    void deallocateThreadStructures();

    void initialize(MRQ_MINLPProb *prob = NULL, MRQ_BranchAndBound *bb = NULL, MRQ_GeneralSolverParams *milpParams = NULL, MRQ_GeneralSolverParams *nlpParams = NULL);


    protected:

    int userErrorCode; //error code returned by user...

    int nthreads;
    branchAndBound::BBL_BRANCH_STRATEGY origBBLBranchStrategy;
    int constrBranchStrat;
    branchAndBound::BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY parentNodeBoundsStrategy;


    int nI, nbin, nC;
    int lastOANPoints;
    int *binVars, *nonBinVars, *intVars, *contVars;
    int *reverseIntVars;


    int **auxInds;
    double **auxVars;
    double **auxConstrs;


    MRQ_BranchAndBound *mybb;
    MRQ_MINLPProb *prob;
    MRQ_GeneralSolverParams *milpSolverParams;
    MRQ_GeneralSolverParams *nlpSolverParams;
    
    
    double *oplc; //we do not save a pointer to puc. We assume puc is allocated m positions after plc
    
    long unsigned int *nPCostAvgAboveErrorEstimative;
    double *pCostAvgAboveErrorEstimative;

    MRQ_OuterApp *oa;
    MRQ_Mutex SEMAPH_OAPoints;
    MRQ_NewPoints *oaPoints;
    MRQ_NewGlobalCutGenerator *cutGen;
    MRQ_NewPseudoCostCalc *pcosts; //pseudo costs to drive variable choosing in the branching time
    MRQ_BasePseudoCostCalc *pruningPcosts; //pseudo costs to perform pseudo prune

    MRQ_BinSumConstrsInds binSumConstrs;

    MRQ_IGMA2Iteration igma2Iter;
    
    minlpproblem::MIP_BinSumConstrsIndsByClass *ssrBinSumConstrs;
    
    minlpproblem:: MIP_ConstraintsByColumnsStorager *ccstorager;
    


    //thread structures

    bool *useHeuristics;
    bool *useRoundings;
    bool *useOASubs;
    bool *useIGMA2s;

    long unsigned int *beforeExpoPseudoPruningCounter;
    long unsigned int *afterExpoPseudoPruningCounter;

    long unsigned int *numberOfConstrBranchings;
    long unsigned int *numberOfWrongLowerBounds;

    double *pseudoPruningLowestLowerBound;

    double **tplc;
    MRQ_NLPSolver **nlps;
    MRQ_LPboundsUpdater *lpBounds;
    MRQ_Preprocessor *preprocessors;


    MRQ_NewGlobalCutList *thCutLists;
    MRQ_HeuristicExecutor *heurExecs;


    MRQ_OuterApp *oaSubs;
    MRQ_NewUserNodeGenerator2 *userNodeGens;
    MRQ_NewChooseIndexToBranch *chooseIndices;
    MRQ_BinSumConstrsChooser *constrChoosers;
    MRQ_Rounding **roundings;
    MRQ_SSRoundingExecutor **ssRoundings; 
    MRQ_Random **ssrRandoms; //we have a separated random objects to ssr to keep ssr usage deterministic. If we share random objects with other procedures, it can change the behaviour of solution found between two executions having differen paraeters

    MRQ_GapMinProb *gapmins;

    MRQ_Random *randomBoundUpdts;




    MRQ_BBLCallbacks(MRQ_MINLPProb *prob = NULL, MRQ_BranchAndBound *bb = NULL, MRQ_GeneralSolverParams *milpParams = NULL, MRQ_GeneralSolverParams *nlpParams = NULL);

    virtual ~MRQ_BBLCallbacks();


    virtual int beforeAll(const unsigned int numberOfThreads, double *lx, double *ux) override;


    virtual int generateRootNode( branchAndBound::BBL_Node* &rootNode ) override;


    virtual int beforeBBLoop(const unsigned int thnumber, const double lb, const double ub) override;

    virtual void afterBBLoop(const unsigned int thnumber, const double lb, const double ub, const int threadReturnCode) override;


    //virtual int beforeSolvingRelaxation( const unsigned int thnumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, bool &pruneNode);


    virtual int solveSubProblem(const unsigned int thnumber, branchAndBound::BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, branchAndBound::BBL_RETURN_CODES &retCode, double &objValue, double &dualObjValue, double *sol, double *dualSol, bool &generalFeasibleSol, bool &pruneNode, double &nodeLowerBound, branchAndBound::BBL_BRANCH_STRATEGY &branchStrategy) override;


    virtual int chooseIndexToBranch(const int thnumber, branchAndBound::BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, branchAndBound::BBL_RETURN_CODES retCode, double objValue, double dualObjValue, double *sol, double *dualSol, unsigned int &sizeIndices, unsigned int *indices, double *breakValues1, double *breakValues2, branchAndBound::BBL_Node* &nodes) override;


    virtual int generateNodes(const int threadNumber, branchAndBound::BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, branchAndBound::BBL_RETURN_CODES retCode, const double objValue, double *sol, double *dualSol, branchAndBound::BBL_UserNodeGenerator &userNodeGenerator ) override;


    virtual int endOfIteration(const int threadNumber, const long unsigned int iter, const double cpuTime, const double wallTime, const double lb, const double ub, branchAndBound::BBL_Node &node, const double *nlx, const double *nux) override;


    virtual int updatingBestSolution(const int threadNumber, double* sol, double &objValue, const double ub, const long unsigned int iter) override;


    virtual void newBestSolution( const int threadNumber, const double *newSol, const double oldBestObj, const double newBestObj, const long unsigned int iter ) override;


    virtual void afterAll(const long unsigned int iters, const double cpuTime, const double wallTime, const double lb, const double ub) override;


    friend MRQ_BranchAndBound;
    friend MRQ_UserCallbacks;
};



















}

#endif
