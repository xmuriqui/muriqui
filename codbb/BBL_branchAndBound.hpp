
#ifndef _BBL_BRANCHANDBOUND_HPP
#define _BBL_BRANCHANDBOUND_HPP

#include <ctime>

#include <new>
#include <iostream>


#include "BBL_constants.hpp"
#include "BBL_node.hpp"






#if BBL_CPP_MULTITHREADING
    #include <thread>
    #include <mutex>
    #include <condition_variable>
#endif


#if BBL_OMP_MULTITHREADING
    #include <omp.h>
#endif





namespace branchAndBound{
    
    
    
    class BBL_PruneCounter
    {
    public:
        long unsigned int bound;	 //to count prunes by bound
        long unsigned int infeas; //to count prunes by infeasibility
        long unsigned int opt;	 //to count prunes by optimality, i.e., because a feasible was found
        long unsigned int user;	 //to count prunes requested by user
        
        BBL_PruneCounter();
        
        void reset();
        
        void accumulate( BBL_PruneCounter& other );
    };
    
    
    class BBL_MUTEX_BASIS_TYPE;
    
    //mutex to implement mutual exclusion on multithread proceduring. In this way, the rest of code do not need worry about haw library to multithreading is being used...
    class BBL_Mutex
    {
        
    public:
        
        #if BBL_CPP_MULTITHREADING
            std::mutex mymutex;
        #else
            #if BBL_OMP_MULTITHREADING
                omp_lock_t mutex;
            #endif
        #endif
        
        
        BBL_Mutex();
        
        void initialize();
        
        int lock( const unsigned int nthreads );
        
        int tryLock( const unsigned int nthreads );
        
        int unlock( const unsigned int nthreads );
        
        void destroy();
        
        ~BBL_Mutex();
    };
    
    
    //class to implement multi_thread proceduring
    
    
    class BBL_MTNodeListManager;
    
    class BBL_NodeList
    {   
        int expStrat;
        unsigned long int nNodes;
        BBL_Node *head, *tail, *iter; //iter is used only to insert nodes...
        double firstlb; //we use that variable to allow getFirstNodeLowerBound gets first lower bound whitout acces head pointer and so, avoid semaphore use to do it on superior levels
        
    public:   
        
        BBL_NodeList(const int strategy = BBL_ES_BEST_LIMIT);
        
        void initialize(const int strategy = BBL_ES_BEST_LIMIT);
        
        //Return the number of nodes pruned. If maxLevelPruneCounter > 0, that method accumulates in pruneLevelCounter the number of prunes by level
        long unsigned int pruneNodesByBound(const double zu, double& zl, const unsigned int maxLevelPruneCounter, unsigned int* pruneLevelCounter = 0);
        
        //return the number of nodes pruned. If maxLevelPruneCounter > 0, that method accumulates in pruneLevelCounter the number of prunes by level
        long unsigned int getNodePointer(const double zu, branchAndBound::BBL_Node* &node, const unsigned int maxLevelPruneCounter = 0, unsigned int* pruneLevelCounter = NULL);
        
        //we assume that the nodes are encadeate and they have the same parent
        void insertNodes(const unsigned int numberOfNodes, BBL_Node* nodes);
        
        unsigned long int deallocateAllNodes(void);
        
        unsigned long int countNodes(void) const;
        
        //that method return a nonzero value if the list is empty
        int getFirstNodeDepth(const double zu, unsigned int& depth, long unsigned int& nNodesPruned);
        
        //that method return a nonzero value if the list is empty
        int getFirstNodeLowerBound(double& lb) const;
        
        //int getFirstNodefDad(double &fdad);
        
        //that function returns the number of nodes encadeated in variable nodes;
        long unsigned int getLastNodes(const long unsigned int numberOfNodes, BBL_Node* &nodes);
        
        //that function returns the number of nodes encadeated in variable nodes;
        long unsigned int getFirstNodes(const long unsigned int numberOfNodes, BBL_Node* &nodes);
        
        inline unsigned long int getNumberOfNodes() const
        {
            return nNodes;
        }
        
        void print(void) const;
        
        void reorganizeToBestLimit(void);
        
        //void reorganizeTofDad(void);
        
        ~BBL_NodeList();
        
        friend class BBL_MTNodeListManager;
    };




    //that class implements an array of BBL_NodeList's. Each BBL_NodeList is responsable to save nodes whose lower bound is in a specific interval. The main objective is remediating the needed time to insert nodes in the list at best limit strategies.  If either deep or width startegy is used, the object implement only one BBL_NodeList since the access in those cases is done in a constant time.
    class BBL_NodeListManager
    {
        int expStrat;
        unsigned int nNodeLists;
        BBL_NodeList *nodeLists; //vector with nodes list. Width and dept strategy only use one nodelist. Best limit strategy use several nodes sublists sorted by lower bound.
        double *llb; //"lowest lower bound" for each node list
        double reorgFactorToLastEmptyList;
        
        unsigned long int nNodes;
        
        unsigned int *pruneLevelCounter;
        
        /*
        * calculate the intervals of lower bounds to list. Each list i will have equal intervals unless the last, that will get all nodes with lb greater than endLBInterval. The list i will get nodes in the interval begLBInterval + i*(endLBInterval - begLBInterval)/(nNodeLists - 1). Moreover, the first lower bound is set to -INFINITY.
        */
        void calculateLBintervalsForLists(const double begLBInterval, const double endLBInterval);
        
        
        
        
    public:
        
        
        
        BBL_NodeListManager();
        
        void initialize();
        
        
        
        void acumulatePruneLevelCounters(const unsigned int maxLevelPruneCounter, branchAndBound::BBL_PruneCounter* counter) const;
        
        
        /*
        * begLBInterval and endLBInterval are estimatives about the interval where lower bounds of nodes will be. There is no problem if a node will have lower bound out of this interval. Even so, it will work.
        */
        int allocateLists(const int strategy, const unsigned int nNodeLists, const double begLBInterval, const double endLBInterval, const unsigned int maxLevelPruneCounter = 0);
        
        unsigned long int countNodes(void) const;
        
        unsigned long int deallocateAllNodes(void);
        
        void deallocateLists();
        
        
        //that function returns the number of nodes encadeated in variable nodes;
        long unsigned int getLastNodes(long unsigned int numberOfNodes, BBL_Node* &nodes);
        
        //that function returns the number of nodes encadeated in variable nodes;
        long unsigned int getFirstNodes(long unsigned int numberOfNodes, BBL_Node* &nodes);
        
        
        //that method give the depth of first node in the list whose lower bound is above zu and it returns a nonzero value if the list is empty
        int getFirstNodeDepth(const double zu, unsigned int& depth, long unsigned int& nNodesPruned);
        
        //that method returns a nonzero value if the list is empty
        int getFirstNodeLowerBound(double& lb) const;
        
        long unsigned int getNodePointer(const double zu, branchAndBound::BBL_Node* &p, const unsigned int maxLevelPruneCounter = 0);
        
        inline unsigned long int getNumberOfNodes() const
        {
            return nNodes;
        }
        
        //we assume that the nodes are encadeate and they have the same parent
        void insertNodes(const unsigned int numberOfNodes, BBL_Node* nodes);
        
        void print(void) const;
        
        long unsigned int pruneNodesByBound(const double zu, double& zl, const unsigned int maxLevelPruneCounter = 0);
        
        void reorganizeToBestLimit(const double begLBInterval, const double endLBInterval);
        
        void reorganizeNodeLists();
        
        ~BBL_NodeListManager();
    };
    
    
    /*
    * That class manager a branch-and-bound list node to paralel (multithread) environments. It has its own mutex (semaphores to manage the access to nodes) 
    */
    class BBL_MTNodeListManager
    {
        int expStrat;
        unsigned int nThreads;
        unsigned int nNodeLists; //on the first versions, nNodeLists = nThreads
        BBL_NodeListManager *nodeLists;
        
        
        BBL_Mutex *SEMAPH_lists; //to protect each list about removing and insertion
        
        BBL_Mutex *SEMAPH_listQueries; //to protect each list about the verification on getNode...
        
        
        
    public:
        
        BBL_MTNodeListManager();
        
        void initialize();
        
        void acumulatePruneLevelCounters(const unsigned int maxLevelPruneCounter, BBL_PruneCounter *counter) const;
        
        /*
        * begLBInterval and endLBInterval are estimatives about the interval where lower bounds of nodes will be. There is no problem if a node will have lower bound out of this interval. Even so, it will work. nNodeSubList is the number of paralel lists for each node (e.g. 100 or 1000).
        */
        int allocateLists(const int strategy, const unsigned int numberOfThreads, const unsigned int nNodeSubLists, const double begLBInterval, const double endLBInterval, const unsigned int maxLevelPruneCounter = 0);
        
        
        unsigned long int countNodes(void) const;
        
        unsigned long int desallocateAllNodes(void);
        
        void deallocate();
        
        
        //get the first node whose lower bound is above zu and return the number of pruned nodes.
        long unsigned int getNodePointer(const double zu, BBL_Node* &p, const unsigned int maxLevelPruneCounter = 0);
        
        /*that method gets a number of nodes in the list. There is no any statement about the nodes. Just return n first nodes founded. Note, not necessarilly are  the nodes having lowest lower bound or something like this. This function just return a number of arbitrary nodes. This method does not prune nodes by bound...
        that method returns the number of nodes encadeated in variable nodes;
        */
        long unsigned int getNodes(long unsigned int numberOfNodes, BBL_Node* &nodes, bool sortNodesByBound = false);
        
        
        long unsigned int getNumberOfNodes(void) const;
        
        
        /* we assume that the nodes are encadeate and they have the same parent. Here, nodes can have different lower bounds. To insert nodes, you need specify the threadNumber (starting from zero). Each threads must use unique a integer number in [0  (nThreas-1)]. Two different threads must not use the same number!
        */
        void insertNodesDifsLowerBounds(const unsigned int listNumber, const unsigned int numberOfNodes, BBL_Node* nodes);
        
        /* we assume that the nodes are encadeate and they have the same parent and the same lower bound. To insert nodes, you need specify the threadNumber (starting from zero). Each threads must use unique a integer number in [0  (nThreas-1)]. Two different threads must not use the same number!
        */
        void insertNodes(const unsigned int listNumber, const unsigned int numberOfNodes, BBL_Node* nodes);
        
        //if listNumber < 0, we print all lists...
        void print(const int listNumber = -1) const;
        
        unsigned int pruneNodesByBound(const double zu, double& zl, const unsigned int maxLevelPruneCounter = 0);
        
        void reorganizeToBestLimit(const double begLBInterval, const double endLBInterval);
        
        void reorganizeNodeLists(const unsigned int listNumber);
        
        ~BBL_MTNodeListManager();
    };



    class BBL_HistorySolution
    {
        
        int allocateSol(const int n);
        
    public:
        int nvars; //total number of variables
        long unsigned int iter; //iteration when that solution was gotten
        double time; //time (wall) when that solution was gotten
        double cputime; //cpu time when that solution was gotten
        double objF;
        double *sol; //solution
        
        BBL_HistorySolution();
        
        BBL_HistorySolution(const int n, const long unsigned int iter, const double time, const double cputime, const double* sol, const double objF);
        
        ~BBL_HistorySolution();
        
        void freeSolution();
        
        inline int getnvars() const
        {
            return nvars;
        }
        
        inline int getiter() const
        {
            return iter;
        }
        
        inline double gettime() const
        {
            return time;
        }
        
        int getsolution(double *solution) const;
    };



    class BBL_SolutionHistory
    {
        int nsols;
        BBL_HistorySolution **hsols; //we use a double pointer to turn easier and faster the reallocation process...
        
    public:
        
        BBL_SolutionHistory();
        
        ~BBL_SolutionHistory();
        
        void desallocate();
        
        inline int getnsols() const
        {
            return nsols;
        }
        
        int addSolution(const int n, const long unsigned int iter, const double time, const double cputime, const double *sol, const double objF);
        
        //it is only a pointer, not a copy. Do not free. If you want a copy of solution use the method getsol of MRQ_HistorySolution pointer
        BBL_HistorySolution * getHistSolPointer(const int index) const;
        
    };
    
    
    
    class BBL_BranchAndBound;
    
    

    class BBL_UserNodeGenerator
    {
    private:
        unsigned int n, ndual;
        unsigned int nnodes;
        
        BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY nodesType;
        
        BBL_Node *parent;
        BBL_Array<double> *x;
        BBL_BasePointer<BBL_ParentNodeInfo> *parentBounds; 
        BBL_BasePointer<BBL_ParentNodeInfo> *parentBoundsNoInherit;
        //BBL_FloatOrDoubleNodeBoundsPointer parentBounds;//BBL_ArraySize <BBL_NodeBounds> *parentBounds;
        
        double lb; //lower bound for the nodes
        double *nlx, *nux;
        double *olx, *oux;
        double *auxlx, *auxux;
        double *sol, *dualSol;
        
        int allocate(const unsigned int n);
        
        void deallocate();
        
        void initialize(BBL_Node *parentNode, double *olx, double *oux, double *nlx, double *nux, double lb, double *sol, double *dualSol);
        
        inline int insertNode( branchAndBound::BBL_Node* node );
        
        inline int generateBoundsArray(const double* nodelx, const double* nodeux, BBL_Node &node) const;
        
        int setx(branchAndBound::BBL_Array< double >*& x, const double* sol, const double* dualSol) const;
        
        
        
    public:
        
        bool usePrimalSol;
        bool useDualSol;
        BBL_Node *nodes;
        BBL_Node *tail;
        
        
        
        BBL_UserNodeGenerator(const unsigned int n, const unsigned int ndual, const bool use_parent_primal_sol, const bool use_parent_dual_sol, const BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY nodesType);
        
        ~BBL_UserNodeGenerator();
        
        
        
        /*
        * create a BB node. That node will be added to open nodes lits. nodelx and nodeux describes variable bounds for the new node
        * 
        * Note: the generated node will not inherit parent bounds directly. All bounds are described by nodelx, nodeux. Due to it, it spends more memory to save the node. Other method generateNode is more eficient to save nodes.   
        * 
        * newNode pointer can points to an object where BranchAndBound will represent the node. It can be useful if you are using your own class derived from BBL_Node to represent nodes. If you pass a NULL pointer, BranchAndBound will generate that object using BBL_Node class. Even if you generate newNode by yourself, BranchAndBound will be responsable by delete the newNode in the appropriate moment.
        * 
        * You can use nodeLowerBound if you want calculate a new lower bound to node. BB procedure will use as node lower bound the maximum between that value and the bound gotten from node's parent (curretn node being exploited). So, if you do not want to calculate, you can pass a huge negative value as -BBL_INFINITY.  
        * 
        * 
        * 
        * data will be copy inside that procedure to generate the BB node. So, you can free memory pointed by * nodelx and nodeux.
        * 
        */
        
        
        int generateNode( const double *nodelx, const double *nodeux, BBL_Node *newNode = NULL, double parentLowerBound = -BBL_INFINITY, const double *parentSol = NULL, const double *parentDual = NULL );
        
        
        
        /*
        * create a BB node. That node will be added to open nodes lits.
        * 
        * If flag inheritParentBounds is true, the new node inherits bounds from its parent (the current node being exploited) automatically beyond its own bounds. Note bounds described in newBounds can overwrites bounds from parents when both reference a same subset of variables.
        * 
        * newNode pointer can points to an object where BBL will represent the node. It can be useful if you are using your own class derived from BBL_Node to represent nodes. If you pass a NULL pointer, BBL will generate that object using BBL_Node class. Even if you generate newNode by yourself, BBL will be responsable by delete the newNode in the appropriate moment.
        * 
        * You can use nodeLowerBound if you want calculate a new lower bound to node. BB procedure will use as node lower bound the maximum between that value and the bound gotten from node's parent (curretn node beinx exploited). So, if you do not want to calculate, you can pass a huge negative value as -MRQ_INFINITY.  
        * 
        * If either initSol or initDual is NULL, the new  node will use the optimal solution from its parent to build its initial solution. Note if you wish change it, you need provide both primal and dual solution.
        * 
        * data will be copy inside that procedure to generate the BB node. So, you can free memory pointed by * BBL_NodeBounds newBounds after call this fucntion.
        * 
        */
        
        int generateNode( const unsigned int nNewBounds, const BBL_NodeBoundsSol *newBounds, const bool inheritParentBounds = true, const bool isNewBoundsAscOrdered=false, BBL_Node *newNode = NULL, double nodeLowerBound = -BBL_INFINITY, const double *parentSol = NULL, const double *parentDual = NULL);
        
        
        unsigned int getNumberOfNodes() const;
        
        BBL_Node * getNodesPointer() const;
        
        
        
        friend class BBL_BranchAndBound;
    };
    
    
    
    class BBL_UserCallbacks
    {
        
    private:
        
        BBL_BranchAndBound *bb;
        
        
        friend class BBL_BranchAndBound;
        
        
    protected:
        
        //return true if solution updates the best solution
        bool tryUpdateBestSolution(const int threadNumber, double* solution, const double fsolution, const long unsigned int iter);
        
        bool tryUpdateLowerBound(const double lb);
        
        int getBestSolutionCopy(double* solution, double& fsolution) const;
        
        double getBestSolutionObj() const;
        
        double getLowerBound() const;
        
        long unsigned int getNodes(long unsigned int numberOfNodes, BBL_Node* &nodes, bool sortNodesByBound);
        
        long unsigned int getNumberOfOpenNodes() const;
        
        bool getSomeBoundPrune() const;
        
        bool getSomeInfeasPrune() const;
        
        bool getSomeOptPrune() const;
        
        bool getSomeUserPrune() const;
        
        double getUpperBound() const;
        
        void getVariableBoundsArrayPointers( double* &lx, double* &ux ) const;
        
        //that function lock (put to sleep) all threads unless the main thread (Thread 0) before they enter in the branch-and-bound loop. Note, that function has no effect if you call after those threads alredy had entered in the loop. So, you should use this function before those threads be created. Best places are beforeAll and generateRootNode... To unlock the threads and let them run the BB loop, call the function to unlock those threads
        void lockAuxiliaryThreadsBeforeBBLoop();
        
        void unlockAllAuxiliaryThreads();
        
        bool hasFeasibleSolution() const;
        
        
    public:
        
        
        void printOpenNodesList() const;
        
        
        
        /*
        * That method is called in the beggining of the algorithm
        */
            
        virtual int beforeAll(const unsigned int numberOfThreads, double *lx, double *ux)
        {
            return 0;
        }
        
        
        /* that method is to generate the root node in the branch-and-bound. It can be useful if you are generating nodes from your own class derived from BBL_Node. If that is not the case, you do not have to do anything and BBL code will generate the rootNode from BBL_Node class. All you have to do is generate a object from a class derived from BBL_Node on variable rootNode. You do not need set any attribute or call any method of BBL_Node. If you let rootNode == NULL, rootNode will be generated automatically. You can generate multiple root nodes encadeating them in rootNode.
        * 
        */
        virtual int generateRootNode( BBL_Node* &rootNode )
        {
            return 0;
        }
        
        
        /*
        * that method is called, for each thread, just once time by run, before respective thread enter in the branch-and-bound loop. It can be useful to perform some initialization for each thread before they explore nodes.
        */ 
        virtual int beforeBBLoop(const unsigned int threadNumber, const double lb, const double ub)
        {
            return 0;
        }
        
        
        /*
        * that method is called, for each thread, just once time by run, after respective thread exits from the branch-and-bound loop. It can be useful to perform some initialization for each thread before they explore nodes.
        */ 
        virtual void afterBBLoop(const unsigned int threadNumber, const double lb, const double ub, const int threadReturnCode)
        {
            //return 0;
        }
        
        
        
        /*method called before BB solve the relaxation in the current node of Branch-and-Bound algorithm. That method is called if the parameter in_call_before_solving_relax_callback is set to true.
        *
        * nlx and nux points to current bounds to variables in the node and will be used in solving process. So, you can change those bounds, for example, to fix aditional variables. Take care to do it to do not damaging the BB correctness. The changes of bounds in a node will not be propagated to its node descendants (If you want do it, you need set again the respective bounds in node descendants).
        * 
        * 
        * 
        * pruneNode is a flag that allow the user prune the current node before solving the relaxation. Make sure set that flag appropriately.
        * 
        * 
        */
        
        virtual int beforeSolvingRelaxation( const unsigned int threadNumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, bool &pruneNode)
        {
            //pruneNode = false;
            //return 0;
            return BBL_USER_CALLBACK_NOT_IMPLEMENTED; //we put it here to bbl return a error if user do not overload this method correctly
        }
        
        
        /* Method called to solve subproblems in the Branch-&-Bound (B&B) tree.
        * 
        * threadNumber: number of thread that is calling the callback function
        * 
        * node: curren node being exploited
        * 
        * iter: current iteration (in this thread)
        * 
        * lb: global lower bound in the B&B procedure
        * 
        * ub: global upper bound in the B&B procedure
        * 
        * nlx: lower bounds for variables in the current node
        * 
        * nux: upper bounds for variables in the current node
        * 
        * retCode: (output) final status of subproblem solving
        * 
        * objValue: (output) objective function value in the solution sol
        * 
        * dualObjValue: (output) dual objective value for the current node. Note that is considered as a strict lower bound to the current node. If you set parameter in_use_dual_obj_to_bound_prunning, this value will be consider to bound prune instead of objValue
        * 
        * sol: (output) solution being returnded. Off course, if the current node is infeasible, you will nos set this array
        * 
        * dualSol (output): dual solution being returned. You only need that if you are storing parent dual solution in the nodes by means of parameter in_store_parent_dual_solution_on_nodes
        * 
        * generalFeasibleSol: (output) set this flag if the current solution is feasible to the problem being addreesed by the B&B procedure. Note if this flag is true and retCode is BBL_OPTIMAL_SOLUTION, the current node will be pruned by optimality.
        * 
        * pruneNode: (output) set this flag if you want prune the current node by some reason. Note you do not need be worry about optimality, bound and infeasibility prunes because they are already done by the BBL_BranchAndBound 
        * 
        * nodeLowerBound: (output) value to node lower bound. If you do not set this, parent lower bound, objValue or dualObjValue will be used. 
        * 
        * branchStrategy: (output) this variable is initialized with the value of bb parameter in_branching_strategy. It is usefull if user wishes change the branching strategy during B&B...
        * 
        */ 
        
        
        virtual int solveSubProblem(const unsigned int threadNumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, BBL_RETURN_CODES &retCode, double &objValue, double &dualObjValue, double *sol, double *dualSol, bool &generalFeasibleSol, bool &pruneNode, double &nodeLowerBound, BBL_BRANCH_STRATEGY &branchStrategy) = 0;
        
        
        
        
        
        
        /* That method can be used to choose a index to perform branching in the appropriate time. That method only will be called if the branching strategy is set to BBL_BS_USER_INDEX_CHOICE. Here, we generate 2^g new nodes, where g is the number of indices taken to branch. Note that you can perform branch on several indices at same time. Branch-and-Bound will generate all nodes covering all 2^g partitions of space. Note for each index, the procedure will generate two partitions whose bounds are determined by breakValue1 and breakValue2.
        * 
        * 
        * threadNumber: number of thread that is calling the callback function
        * iter: number of the current iteration
        * lb: current lower bound to problem
        * nlx: array with lower bounds in the current node
        * nux: array with upper bounds in the current node
        * retCode: BBL_RETRUN_CODE of the subproblem solving process
        * objValue: objective function on solution
        * dualObjValue: dual objective function
        * sol: solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * dualSol: dual solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * 
        * sizeIndices: (output) number of indices choosen to branch (> 0). Note that the number of new nodes is 2^sizeIndices
        * 
        * indices: (output) array containing the indices to perform branch.
        * 
        * breakValue1: (output) break values used to perform branch. Let j be the index in indices[i], first partition will stay in (nlx[j]  breakValue1[i] ). Off course, nlx[j] <= breakValue1[i] <= nux[j]. Dimension of breakValue1 is the same of indices.
        * 
        * breakValue2: (output) break values used to perform branch. Let j be the index in indices[i], second partition will stay in (breakValue1[i]  nux[j]). Off course, nlx[j] <= breakValue1[i] <= nux[j]. Dimension of breakValue2 is the same of indices.
        * 
        * 
        * nodes: (output) optional argument to pass a encadeate list of nodes generated by user. it is useful in the cases where user are working with their own class (dervedo from BBL_Node) to represent branch-and-bound nodes. BBL initialize this argument as NULL. If user does not generate BBL_Node in this pointer, BBL will perform it.
        * 
        */
        
        virtual int chooseIndexToBranch(const int threadNumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, BBL_RETURN_CODES retCode, double objValue, double dualObjValue, double *sol, double *dualSol, unsigned int &sizeIndices, unsigned int *indices, double *breakValues1, double *breakValues2, BBL_Node* &nodes)
        {
            return BBL_USER_CALLBACK_NOT_IMPLEMENTED; //we put it here to bbl return a error if user do not overload this method correctly
        }
        
        
        
        /*
        * That method can be used to users generate new nodes in branching procedure by themselves. In this way, an user can apply his own strategy to branching the space. That method only will be called if the branching strategy is set to BBL_BB_BS_USER_NODE_GENERATION
        *
        *
        * threadNumber: number of thread that is calling the callback function
        * iter: number of the current iteration
        * lb: current lower bound to problem
        * ub: current upper bound to problem
        * nlx: array with lower bounds in the current node
        * nux: array with upper bounds in the current node
        * status: MRQ_RETRUN_CODE of the relaxation solving process
        * * fsolution: objective function on solution
        * sol: solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * constrs: constraint values in the solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * dualSol: dual solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * userNodeGenerator: object that allow user generate new nodes. User must call the method generateNode to do it. Atention: nodeBounds of new node should be ordered by index of variables in a ascendet way.
        * 
        */
        virtual int generateNodes(const int threadNumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, BBL_RETURN_CODES retCode, const double objValue, double *sol, double *dualSol, BBL_UserNodeGenerator &userNodeGenerator )
        {
            return BBL_USER_CALLBACK_NOT_IMPLEMENTED; //we put it here to bbl return a error if user do not overload this method correctly
        }
        
        
        
        /*
        * That method is called at the end of each iteration.
        */
        
        virtual int endOfIteration(const int threadNumber, const long unsigned int iter, const double cpuTime, const double wallTime, const double lb, const double ub, BBL_Node &node, const double *nlx, const double *nux)
        {
            return BBL_USER_CALLBACK_NOT_IMPLEMENTED; //we put it here to bbl return a error if user do not overload this method correctly
        }
        
        
        
        /*
        * That method is called when imediatelly before algorithms try update the best solution found. You can change the objective value or some variable values under your own responsibility. 
        * 
        * After that callback, algorihms will perform the test on objValue to update or no the best solution found. Note, if objValue is not better than best solution, the best solution will not be updated!
        * 
        * sol: solution being used to try update best solution
        * objValue: objective value in sol. Note you can change that value. If you do not change, the original objective value is used
        *
        * ub: current upper bound. Note when that callback is called, the respective algorithm still did not try update the solution.
        * 
        * 
        */
        virtual int updatingBestSolution(const int threadNumber, double* sol, double &objValue, const double ub, const long unsigned int iter)
        {
            return 0;
        }
        
        
        /*
        * That method is called when the best solution is effectvelly updated. Note, it is called after the updating.
        * 
        */ 
        
        virtual void newBestSolution( const int threadNumber, const double *newSol, const double oldBestObj, const double newBestObj, const long unsigned int iter )
        {
        }
        
        
        /*
        *That method is called at the end of the Algorithm
        */
        virtual void afterAll(const long unsigned int iters, const double cpuTime, const double wallTime, const double lb, const double ub)
        {
            
        }
        
        
    };
    
    
    /*
    * That class implements a framework for paralell branch-and-bound.
    */
    class BBL_BranchAndBound
    {
        
    public:
        
        bool in_call_after_bb_loop_callback;
        bool in_call_before_bb_loop_callback;
        bool in_call_before_solving_relax_callback;
        bool in_call_end_of_iteration_callback;
        bool in_call_new_best_solution_callback;
        bool in_call_updating_best_solution_callback;
        
        bool in_count_total_prunes;
        bool in_consider_relax_infeas_if_solver_fail;
        bool in_prune_nodes_by_bound;
        bool in_reorganize_lists;
        bool in_store_history_solutions;
        bool in_store_parent_dual_solution_on_nodes;
        bool in_store_parent_primal_solution_on_nodes;
        bool in_use_dual_obj_to_bound_prunning;
        
        
        unsigned int in_lists_reorganization_frequency;
        unsigned int in_max_tree_level_to_count_prunes_by_level;
        unsigned int in_number_of_node_sublists;
        unsigned int in_number_of_threads;
        unsigned int in_print_level;
        unsigned int in_printing_frequency;
        
        BBL_BRANCH_STRATEGY in_branching_strategy;
        BBL_EXT_EXP_STRATEGY in_exp_strategy;
        BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY in_parent_node_bounds_strategy;
        
        long unsigned int in_max_iterations;
        
        double in_absolute_convergence_tol;
        double in_infinity_value;
        double in_lower_bound;
        double in_max_time;
        double in_max_cpu_time;
        double in_relative_convergence_tol;
        double in_upper_bound;
        
        
        
        
        
        bool out_feasible_sol;
        unsigned long int out_number_of_feas_sols;
        unsigned long int out_number_of_iterations;
        unsigned long int out_number_of_open_nodes;
        unsigned int out_number_of_threads;
        int out_return_code;
        int out_return_subcode;
        
        long unsigned int out_first_sol_iter;
        long unsigned int out_best_sol_iter;
        
        double out_clock_time;
        double out_cpu_time;
        double out_clock_time_to_fisrt_sol;
        double out_cpu_time_to_fisrt_sol;
        double out_clock_time_to_best_sol;
        double out_cpu_time_to_best_sol;
        double out_lower_bound;
        double out_upper_bound;
        double out_obj_opt_at_root_relax;
        
        double *out_best_sol; //best solution
        double out_best_obj;
        
        BBL_SolutionHistory out_sol_hist;
        
        BBL_PruneCounter out_prune_counter;
        BBL_PruneCounter *out_prune_counter_by_level;
        
        
        
        BBL_BranchAndBound();
        
        virtual ~BBL_BranchAndBound();
        
        virtual void resetOutput();
        
        virtual void resetParameters();
        
        virtual int run(const unsigned int nPrimalVars, const unsigned int ncons, const unsigned int nDualVars, const double *lx, const double *ux, branchAndBound::BBL_UserCallbacks& userCalbacks);
        
        
        
        
    protected:
        
        bool __lockAuxThreadsBeforeLoop; //this flag is used to lock threads (except thread 0) before they enter in the main loop
        
        BBL_Mutex SEMAPH_lockAuxThs;
        
        
        bool someBoundPrune, someInfeasPrune, someOptPrune, someUserPrune;
        
        unsigned int n, m, ndual;
        
        int expStrategy;
        double zl, zu;
        double *lx, *ux;
        
        //thread structures:
        unsigned int nthreads;
        bool endThreads;
        signed char *thRunning;
        int *thReturnCodes;
        int *thReturnSubCodes;
        double *thlbCurrentNode;
        //end of thread structures
        
        long unsigned int iter;
        clock_t clockStart;
        double timeStart;
        
        BBL_MTNodeListManager *nodes;
        
        
        BBL_Mutex SEMAPH_SolBounds; //to protect best sol, lower and upper bound
        BBL_Mutex SEMAPH_Attrib; //to protect class' attributes in general (unless bounds)
        BBL_Mutex SEMAPH_Print; //to protect printing
        BBL_Mutex SEMAPH_DelNodes; //to protect nodes deletion
        
        
        
        void lockAuxThreadsBeforeLoop();
        
        void unlockAuxThreads();
        
        
        virtual int algorithmInitialization(const unsigned int n, const unsigned int m, const unsigned int ndual, const unsigned int nthreads, const double* lx, const double* ux);
            
        virtual void algorithmFinalization();
        
        virtual int allocateBoundArrays(const unsigned int n);
        
        virtual int allocateBestSol(const unsigned int n, const unsigned int m);
        
        virtual int allocatePruneCountersByLevel(const int nlevels);
        
        virtual int allocateThreadData(const unsigned int nthreads);
        
        virtual bool checkTerminationCriterions(const int threadNumber, const double zl, const double zu, long unsigned int iter, int& retCode);
        
        virtual int bbLoop( branchAndBound::BBL_UserCallbacks& userCalbacks, const unsigned int thnumber);
        
        virtual void desallocateBestSol();
        
        virtual void desallocateBoundArrays();
        
        virtual void desallocatePruneCountersByLevel();
        
        virtual void desallocateThreadData();
        
        virtual int generateNodes(const unsigned int nBranchVars, unsigned int* branchVars, const double* breakValues1, const double* breakValues2, branchAndBound::BBL_Node* parent, const double* lxp, const double* uxp, const double nodelb, const double* solp, const double* dualSolp, const bool use_parent_primal_sol, const bool use_parent_dual_sol, branchAndBound::BBL_Node*& nodes );
        
        virtual double getBestSolutionObj();
        
        virtual int getBestSolutionCopy(double *solution, double &fsolution);
        
        long unsigned int getNodes(long unsigned int numberOfNodes, BBL_Node* &nodes, bool sortNodesByBound);
        
        inline void incBoundPruneCounters( BBL_PruneCounter &pcounter, BBL_PruneCounter *pcounterlevel, const unsigned int depth );
        
        inline void incInfeasPruneCounters( BBL_PruneCounter &pcounter, BBL_PruneCounter *pcounterlevel, const unsigned int depth );
        
        inline void incOptPruneCounters( BBL_PruneCounter &pcounter, BBL_PruneCounter *pcounterlevel, const unsigned int depth );
        
        inline void incUserPruneCounters( BBL_PruneCounter &pcounter, BBL_PruneCounter *pcounterlevel, const unsigned int depth );
        
        virtual void initializeProtectedData();
        
        virtual bool tryUpdateBestSolution(const unsigned int threadNumber, const long unsigned int iter, double* sol, double objValue, BBL_UserCallbacks& user_calbacks);
        
        virtual bool tryUpdateLowerBound(const double lb);
        
        
        friend class BBL_UserCallbacks;
        
        friend int BBL_bbLoop( BBL_UserCallbacks* userCalbacks, BBL_BranchAndBound *bb, const unsigned int number );
        
    };
    
    
    int BBL_bbLoop(  branchAndBound::BBL_UserCallbacks* userCalbacks, branchAndBound::BBL_BranchAndBound* bb, const unsigned int number ); //that is ridiculus but you cannot pass arguments for thread as reference. Only copy or pointer...


}



#endif




