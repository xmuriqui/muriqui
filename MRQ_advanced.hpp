

#ifndef MRQ_ADVANCED_HPP_
#define MRQ_ADVANCED_HPP_

#include <ctime>

#include <vector>

#include "BBL_branchAndBound.hpp"

#include "MRQ_constants.hpp"
#include "MRQ_dataStructures.hpp"
#include "MRQ_algClasses.hpp"
#include "MRQ_bb.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_milpCallbacks.hpp"



namespace muriqui
{
    /*class to allow users intervene directly on branch-and-bound and other algorithms. Objects from that class implements callback functions called in all iterations (or almost all). In this way, user can give directives to change the standard behavior of branch-and-bound procedure, for example, fixing aditional variables, adding and removing cuts, pruning nodes by itself, etc.
    * 
    * If you need a copy of the current best solution, call the method get a best solution copy. If you want update the current upper bound, call the method tryUpadteUpperBound. Make sure pass correct values to soution and its objective function value.
    * 
    * Note: If some of those functions returns a number different of zero, the algorithm execution is aborted.
    * 
    * 
    */ 
    
    typedef MRQ_NewBBNode MRQ_BBNode;
    typedef MRQ_NewUserNodeGenerator2 MRQ_UserNodeGenerator;
    typedef MRQ_NLPSolver MRQ_NLPSolver;
    
    
    
    class MRQ_SolutionStorer
    {
    public:
        
        std::vector< std::pair<double, double*>  > sols;  //we store the objective value and an array to solution
        
        ~MRQ_SolutionStorer();
        
        int addSolution(const unsigned int n, const double objValue, const double *solution);
        
        void deallocate();
    };
    
    
    
    class MRQ_UserCallbacks
    {
        
    private:
        
        MRQ_Algorithm *alg;
        MRQ_MINLPProb *prob;
        //MRQ_GlobalCutGenerator *cutGen; //deprecated
        MRQ_NewGlobalCutGenerator *cutGenerator;
        
        //branchAndBound::BBL_UserCallbacks *bblCallBacks;
        MRQ_BBLCallbacks *bblCallBacks;
        
        clock_t clockStart;
        double timeStart;
        
        
    protected:
        
        
        bool tryUpdateBestSolution(const int threadNumber, const int n, double* solution, const double fsolution, const long unsigned int iter);
        
        int getBestSolutionCopy(const int n, double *solution, double &fsolution);
        
        double getLowerBound() const;
        
        double getUpperBound() const;
        
        //(only for nonlinear BB)
        int BB_getNumberOfOpenNodes(long unsigned int &nnodes) const;
        
        /*(only for nonlinear BB): get a number of branch-and-bound open nodes. User takes the ownership of the nodes. Note, nodes goten are take of from open nodes list. It is user responsibility explore this nodes...
        */ 
        long unsigned int BB_getOpenNodes(long unsigned int numberOfNodes, MRQ_NewBBNode* &nodes, bool sortNodesByBound);
        
        // that function set the cut (only for nonlinear BB):
        //  lb <=  a'x + 0.5x'Qx <= ub
        //only put the lower triangle of Q
        int addGlobalQuadraticCut(const int nza, const int *acols, const double *avalues, const int nzQ, const int *Qrows, const int* Qcols, const double *Qvalues, const double lb, const double ub);
        
        void BB_printOpenNodesList();
        
        /*only for BB*/
        MRQ_Preprocessor* BB_getPreprocessorPointer(unsigned int threadNumber);
        
        /*only for BB*/
        minlpproblem:: MIP_ConstraintsByColumnsStorager * BB_getConstraintsByColumnsStoragerPointer();
        
        
        friend class MRQ_Algorithm;
        friend class MRQ_BranchAndBound;
        friend class MRQ_BBLCallbacks;
        
        
    public:
        
        /*
        * That method is called in the beggining of the algorithm execution:
        * 
        * algCode: code of the algorithm being applied
        * 
        * numberOfThreads: number os threads being used by the algorithm
        */
        
        virtual int beforeAll(const MRQ_ALG_CODE algCode, const unsigned int numberOfThreads)
        {
            return 0;
        }
        
        
        /* that method is to generate the root node on branch-and-bound. It can be useful if you are generating nodes from your own class derived from MRQ_BBNode. If that is not the case, you do not have to do anything and Muriqui code will generate the rootNode from MRQ_BBNode class. All you have to do is generate a object from a class derived from MRQ_BBNode on variable rootNode. You do not need set any attribute or call any method of MRQ_BBNode. If you let rootNode == NULL, Muriqui will generate the rootNode
        * 
        */
        
        virtual int BB_generateRootNode( MRQ_NewBBNode* &rootNode )
        {
            return 0;
        }
        
        
        
        /*
        * that method is called, for each thread, just once time by run, before respective thread enter in the branch-and-bound loop. It can be useful to perform some initialization for each thread before they explore nodes.
        * 
        * threadNumber: number of the thread that is calling the callback function
        * 
        * lb: current best lower bound 
        * 
        * ub: current best upper bound
        * 
        * nlpSolver: interface object for nonlinear solver
        * 
        */ 
        virtual int beforeBBLoop(const unsigned int threadNumber, const double lb, const double ub, MRQ_NLPSolver &nlpSolver)
        {
            return 0;
        }
        
        
        /*
        * that method is called, for each thread, just once time by run, after respective thread exits from the branch-and-bound loop. It can be useful to perform some finitialization for each thread after they explore nodes.
        * 
        * threadNumber: number of the thread that is calling the callback function
        * 
        * lb: current best lower bound 
        * 
        * ub: current best upper bound
        * 
        * nlpSolver: interface object for nonlinear solver
        * 
        * threadReturnCode: return code of the respective thread
        * 
        */ 
        virtual void afterBBLoop(const unsigned int threadNumber, const double lb, const double ub, MRQ_NLPSolver &nlpSolver, const int threadReturnCode)
        {
        }
        
        
        
        /*method called before nonlinear Branch-and-Bound solve the relaxation in the current node of Branch-and-Bound algorithm. That method is called if the parameter in_call_before_solving_relax_callback is set to true.
        *
        * nlx and nux points to current bounds to variables in the node and will be used in solving process. So, you can change those bounds, for example, to fix aditional variables. Take care to do it to do not damaging the BB correctness. The changes of bounds in a node will not be propagated to its node descendants (If you want do it, you need set again the respective bounds in node descendants).
        * 
        * nlpSolver is an interface object representing the nlp solver used to solve the relaxation. You can add aditional cuts. Do not forget. Cuts will be valid for all nodes solved by the same thread until you remove the cuts. Have in mind BB can work on several threads of execution and each thread has its own MRQ_nlpSolver object. Each node is exploited by only one unique thread. You can set the initial solution for the solver also.
        * 
        * 
        * pruneNode is a flag that allow the user prune the current node before solving the relaxation. Make sure set that flag appropriately.
        * 
        * 
        */
        
        virtual int BB_beforeSolvingRelaxation(const unsigned int threadNumber, MRQ_NewBBNode &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, MRQ_NLPSolver &nlpSolver, bool &pruneNode)
        {
            pruneNode = false;
            return 0;
        }
        
        
        /* That function is called after BB solve the relaxation in the current node. That method is called if the parameter in_call_after_solving_relax_callback is set to true.
        * 
        * threadNumber: number of the thread that is calling the callback function
        * 
        * iter: number of the current iteration
        * 
        * lb: current lower bound to problem
        * 
        * nlx: array with lower bounds in the current node
        * 
        * nux: array with upper bounds in the current node
        * 
        * nlpSOlver: object representing the nlp solver used to solve the relaxation. You can add aditional cuts. Do not forget. Cuts will be valid for all nodes solved by the same thread until you remove the cuts. Have in mind BB can work on several threads of execution and each thread has its own MRQ_nlpSolver object. Each node is exploited by only one unique thread. You can set the initial solution for the solver also.
        * 
        * status: MRQ_RETRUN_CODE of the relaxation solving process
        * 
        * solution: solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * 
        * objFSolution: objective function on solution
        * 
        * dualObjFSolution: dual objective function on solution returned by solver (if solver does not provide this value, it is value as NaN)
        * 
        * nodeLowerBound: node lower bound. You can update that value and improve the lower bound for the current node. If your lower bound is better than the current lower bound, it will be used in node descendants.
        * 
        * solveAgain: (output) flag to indicate if you would like resolve the relaxation in that same node (it is usefull if you add some cut). That callback method will be called again after the resolution.
        * 
        * pruneNode: (output) flag to prune that node.
        * 
        * branchEvenIntegerSol: (output) flag to sinalize to the current thread if the branching should be done even if the soluiton is integer. That flag is initialized using parameter in_branch_even_integer sol. If you do not wanna change the flag, you do not need do anything.
        */
        
        virtual int BB_afterSolvingRelaxation(const unsigned int threadNumber, MRQ_NewBBNode &node, const long unsigned int iter, const double lb, const double ub, const double *nlx, const double *nux, MRQ_NLPSolver &nlpSolver, const int status, double &objFSolution, double &dualObjFSolution, double *sol, double *constrs, double *dualSolC, double *dualSolV, bool &pruneNode, bool &solveAgain, bool &branchEvenIntegerSol)
        {
            pruneNode = false;
            solveAgain = false;
            return 0;
        }
        
        
        /* That method can be used to choose a index to perform branching if the optimal solution of the relaxation is not integer. That method only will be called if the branching strategy is set to MRQ_BB_BS_USER_INDEX_CHOICE
        * 
        * 
        * threadNumber: number of the thread that is calling the callback function
        * 
        * iter: number of the current iteration
        * 
        * lb: current lower bound to problem
        * 
        * nlx: array with lower bounds in the current node
        * 
        * nux: array with upper bounds in the current node
        * 
        * status: MRQ_RETRUN_CODE of the relaxation solving process
        * 
        * fsolution: objective function on solution
        * 
        * sol: solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * 
        * constrs: constraint values in the solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * dualSol: dual solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * 
        * sizeIndices: (output) number of indices choosen to branch (> 0)
        * 
        * indices: (output) array containing the indices to perform branch.
        * 
        * 
        */
        
        virtual int BB_chooseIndexToBranch(const unsigned int threadNumber, MRQ_NewBBNode &node, const long unsigned int iter, const double lb, const double ub, const double *nlx, const double *nux, const int status, double &fsolution, const double *sol, const double *constrs, const double *dualSolC, const double *dualSolV, int &sizeIndices, unsigned int *indices)
        {
            return 0;
        }
        
        
        
        /*
        * That method can be used to users generate new nodes in branching procedure by themselves. In this way, an user can apply his own strategy to branching the space. That method only will be called if the branching strategy is set to MRQ_BB_BS_USER_NODE_GENERATION
        *
        *
        * threadNumber: number of the thread that is calling the callback function
        * 
        * iter: number of the current iteration
        * 
        * lb: current lower bound to problem
        * 
        * ub: current upper bound to problem
        * 
        * nlx: array with lower bounds in the current node
        * 
        * nux: array with upper bounds in the current node
        * 
        * status: OPT_RETRUN_CODE of the relaxation solving process
        * 
        * fsolution: objective function on solution
        * 
        * sol: solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * 
        * constrs: constraint values in the solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * 
        * dualSol: dual solution gotten by the relaxation (if it is available, i.e., if the relaxation is feasible and it was solved successfully)
        * 
        * userNodeGenerator: object that allow user generate new nodes. User must call the method generateNode to do it. Atention: nodeBounds of new node should be ordered by index of variables in a ascendet way.
        * 
        */
        
        virtual int BB_generateNodes(const unsigned int threadNumber, MRQ_NewBBNode &node, const long unsigned int iter, const double lb, const double ub, const double *nlx, const double *nux, const int status, const double fsolution, const double *sol, const double *constrs, const double *dualSolC, const double *dualSolV, MRQ_NewUserNodeGenerator2 &userNodeGenerator )
        {
            return 0;
        }
        
        
        
        /*
        * That method is called at the end of each iteration in all algorithms.
        * 
        * algCode: code of the algorithm being applied
        * 
        * threadNumber: number of the thread that is calling the callback function
        * 
        * iter: current iteration
        * 
        * cpuTime: cpu time spent by the algorithm
        * 
        * wallTime: wall time (clock) spent by the algoritm
        * 
        * lb: best lower bound for the MINLP problem
        * 
        * ub: best upper bound for the MINLP problem
        * 
        */
        
        virtual int endOfIteration(const MRQ_ALG_CODE algCode, const unsigned int threadNumber, const long unsigned int iter, const double cpuTime, const double wallTime, const double lb, const double ub)
        {
            return 0;
        }
        
        
        /*
        * That method is called imediatelly before algorithms try update the best solution found. You can change the objective value or some variable values under your own responsibility. 
        * 
        * After that callback, algorihms will perform the test on objValue to update or no the best solution found.
        * 
        * sol: solution being used to try update best solution
        * objValue: objective value in sol. Note you can change that value. If you do not change, the original objective value is used.
        *
        * ub: current upper bound. Note when that callback is called, the respective algorithm still did not try update the solution.
        * 
        * iter: iteration which solution was found.
        * 
        */
        virtual int updatingBestSolution(const MRQ_ALG_CODE algCode, const unsigned int threadNumber, double* solution, double &objValue, const double ub, const long unsigned int iter)
        {
            return 0;
        }
        
        
        
        /*
        * That method is called when the best solution is effectvelly updated. Note, it is called after the updating.
        *
        * threadNumber: number of the thread that is calling the callback function
        * 
        * newSol: new best solution
        * 
        * oldBestObj: old best objective value
        * 
        * newBestObj: new best objective value, i.e., the objective value of the new solution
        * 
        */ 
        virtual void newBestSolution( const unsigned int threadNumber, const double *newSol, const double oldBestObj, const double newBestObj, const long unsigned int iter )
        {
        }
        
        
        /*
        * That method is called in linear approximation based branch-and-bounds like LP/NLP-BB and LP-BB before the MILP solver solves a node in its branch-and-bound tree. This method is called if the parameter in_call_before_solve_callback_in_milp_bb is set as true. Warning: This method is not called in nonlinear branch-and-bound!
        * 
        * algCode (input): code of the algorithm being applied
        * 
        * threadNumber (input): number of the thread that is calling the callback function
        * 
        * milpSolverCallbackInterface (input): an object to perform operations on the MILP solver inside callback
        * 
        * solsToBuildLazyConstraints (output): an object to (optionally) pass solutions to algorithm. Lazy constraints with linearizations over problem's nonlinear functions will me built over the constraints stored in this object.  (ATTENTION: THIS FUNCTIONALITY IS NOT WORKING BY NOW (it is not my fault. It is due to cplex limitations))
        * 
        */
        virtual int linearApp_beforeSolveInMILPBB(const MRQ_ALG_CODE algCode, const unsigned int threadNumber, MRQ_MILPSolverCallbackInterface &milpSolverCallbackInterface, MRQ_SolutionStorer &solsToBuildLazyConstraints )
        {
            return 0;
        }
        
        
        /*
        * That method is called in linear approximation based branch-and-bounds like LP/NLP-BB and LP-BB before the MILP solver branching a node in its branch-and-bound tree. This method is called if the parameter in_call_branching_callback_in_milp_bb is set as true. Warning: This method is not called in nonlinear branch-and-bound! By now, user have to personalize branching only by means of MRQ_MILPSolverCallbackInterface object
        * 
        * algCode (input): code of the algorithm being applied
        * 
        * threadNumber (input): number of the thread that is calling the callback function
        * 
        * milpSolverCallbackInterface (input): an object to perform operations on the MILP solver inside callback
        * 
        * 
        */
        virtual int linearApp_branchingInMILPBB(const MRQ_ALG_CODE algCode, const unsigned int threadNumber, MRQ_MILPSolverCallbackInterface &milpSolverCallbackInterface)
        {
            return 0;
        }
        
        
        /*
        * That method is called at the end of the Algorithms.
        *
        * algCode: code of the algorithm being applied.
        * 
        * iters: number of iterations
        * 
        * cpuTime: cpu time spent by the algorithm
        * 
        * wallTime: wall time (clock) spent by the algoritm
        * 
        * lb: best lower bound for the MINLP problem
        * 
        * ub: best upper bound for the MINLP problem
        * 
        */
        virtual void afterAll(const MRQ_ALG_CODE algCode, const long unsigned int iters, const double cpuTime, const double wallTime, const double lb, const double ub)
        {
        }
        
    };
    
    
    //function to call method insideRun. It is usefull, for example, to run an algorithm in a MINLPProblem adopting different variable bounds
    int MRQ_insideRun( MRQ_Algorithm *alg, MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, const int thnumber,  const double insideSolverMaxTime, double *nlx, double *nux );
    
    
}
    

#endif 	//MRQ_ADVANCED_HPP_
