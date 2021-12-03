/*
 * MRQ_algClasses.hpp
 *
 *  Created on: 27/08/2013
 *      Author: yo
 */

#ifndef MRQ_ALGCLASSES_HPP_
#define MRQ_ALGCLASSES_HPP_


/*
 * Classes for algorithms in the new object oriented implementation of Muriqui solver
 *
 *
 * Author: Wendel Melo
 * Date: 16-July-2013
 *
 *
 */

#include <ctime>
#include <vector>

#include "MRQ_constants.hpp"
#include "MRQ_dataStructures.hpp"
#include "BBL_branchAndBound.hpp"



namespace muriqui
{

    typedef minlpproblem::MIP_IntegratedHessian MRQ_IntegratedHessian;
    typedef minlpproblem::MIP_Preprocessing MRQ_Preprocessor;

    typedef optsolvers::OPT_LPSolver  MRQ_LPSolver;

    typedef branchAndBound::BBL_Mutex MRQ_Mutex;


    class MRQ_UserCallbacks;
    class MRQ_DynConstrSetUnity;


    class MRQ_Algorithm
    {
        
    public:
        
        bool in_assume_convexity;
        bool in_call_end_of_iteration_callback;
        bool in_call_new_best_solution_callback;
        bool in_call_update_best_sol_callback;
        bool in_fix_int_vars_from_nlp_relax_sol; //fix integer vars based on theorem from Marcia Fampa's Thesis.
        
        bool in_preprocess_lin_constr;
        bool in_preprocess_quad_constrs;
        bool in_preprocess_obj_function;
        bool in_print_parameters_values;
        bool in_set_special_nlp_solver_params;
        bool in_store_history_solutions;
        bool in_use_initial_solution;
        bool in_use_dynamic_constraint_set;
        
        unsigned long int in_max_iterations;
        unsigned int in_number_of_threads;
        unsigned int in_printing_frequency;
        int in_print_level;

        MRQ_MILP_SOLVER in_milp_solver;
        MRQ_NLP_SOLVER  in_nlp_solver;
        
        double in_absolute_convergence_tol;
        double in_absolute_feasibility_tol;
        double in_integer_tol;
        double in_lower_bound;
        double in_max_time;
        double in_max_cpu_time;
        double in_relative_convergence_tol;
        double in_relative_feasibility_tol;
        double in_upper_bound;
        
        MRQ_UserCallbacks *in_user_callbacks;
        
        
        bool out_feasible_solution;
        unsigned long int out_number_of_feas_sols;
        unsigned long int out_number_of_iterations;
        unsigned long int out_number_of_iterations_to_first_feas_sol;
        unsigned long int out_number_of_iterations_to_best_sol;
        unsigned int out_number_of_threads;
        int out_user_callback_error_code;
        MRQ_RETURN_CODE out_return_code;
        MRQ_ALG_CODE out_algorithm;
        
        double out_clock_time;
        double out_cpu_time;
        double out_clock_time_to_first_feas_sol;
        double out_cpu_time_to_first_feas_sol;
        
        double out_clock_time_to_best_sol;
        double out_cpu_time_to_best_sol;
        
        
        double out_lower_bound;
        double out_upper_bound;
        double out_obj_opt_at_continuous_relax;
        
        MRQ_SolutionHistory out_sol_hist;
        
        double *out_best_sol; //best solution
        double out_best_obj;
        
        
        MRQ_Algorithm();

        virtual ~MRQ_Algorithm();
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL);
        
        void copyParametersFrom(const MRQ_Algorithm &source);
        
        virtual void desallocateInitialSolution();
        
        virtual std::string getAlgorithmInitials() const;
        
        virtual std::string getAlgorithmName() const;
        
        virtual int getBestSolutionCopy(const int n, double *solution, double &fsolution) const;
        
        //size or dcs can be NULL; We just return over arguments havinga valid pointers 
        virtual void getDCS0Array(int *size, MRQ_DynConstrSetUnity *dcs) const;
        
        //size or dcs can be NULL; We just return over arguments havinga valid pointers
        virtual void getDCS1Array(int *size, MRQ_DynConstrSetUnity *dcs) const;
        
        virtual bool isLinearApproximationAlgorithm() const;
        
        virtual void printParameters(std::ostream &out = std::cout) const;
        
        int readParametersFromFile(const char* fileName, const bool printErrorMsgs = true);
        
        //each parameter should be preceded by its type in the file, ex: int in_print_level 4
        int readParametersWithTypeFromFile(const char* fileName, const bool printErrorMsgs = true, const bool printFileOpenError = false);
        
        virtual void resetOutput();
        
        virtual void resetParameters();
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) = 0;
        
        virtual int setInitialSolution(const int n, const double *xI);
        
        virtual int setIntegerParameter(const char *name, const long int value);
        
        virtual int setDoubleParameter(const char *name, const double value);
        
        int setParameters(const MRQ_GeneralSolverParams &params);
        
        //that method should be used to set enumeration parameters
        virtual int setStringParameter(const char *name, const char *value);
        
        

    protected:

        bool run_by_inside;
        int thnumber;
        //int insideNumberOfThreads; //it is not used by branch and bound, only to set number of threads in solvers objects by other algorithms...
        double zl, zu;
        double *xInit;
        double insideSolverMaxTime;
        
        double *nlx, *nux; //we can apply the algorithm on a subproblem, i. e., using new bounds...
        
        //MRQ_MINLPProb *prob;
        
        int ndcs0;
        int ndcs1;
        MRQ_DynConstrSetUnity *dcs0;
        MRQ_DynConstrSetUnity *dcs1;
        
        

        //virtual int addSolutionToHistory(const int n, const double *x, const double objF, const int iter, const double time);

        int allocateBestSol(const unsigned int n);
        
        int allocateDCS0(const unsigned int size);
        
        int allocateDCS1(const unsigned int size);
        
        int allocateInitialSolution(const unsigned int n);
        
        //preproclc and preprocuc can be NULL. In this case, we do not store new bounds for constraints. Anyway, in general, bounds for constranints do not change and minlpproblem::preprocessing, except when a constraint is detected redundant
        virtual MRQ_RETURN_CODE algorithmInitialization(const int nthreads, const bool allocBoundsCopy, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, MRQ_MINLPProb& prob, double* &lx, double* &ux, MRQ_Preprocessor* preprocessor, bool *updtSomeConstrBounds = NULL, double** preproclc = NULL, double** preprocuc = NULL, bool callUserCallbackBeforeAll = true);
        
        virtual int _algorithmInitialization(); //initialization that does not need input parameters 
        
        virtual void algorithmFinalization( const int nthreads, MRQ_MINLPProb& prob, double* &lx, double* &ux );
        
        virtual bool checkTerminationCriterions(const int threadNumber, const double zl, const double zu, long unsigned int iter, const double timeStart, const clock_t& clockStart, MRQ_RETURN_CODE& retCode) const;

        void desallocateBestSol();
        
        void desallocateDCS0();
        
        void desallocateDCS1();
        
        //method to try get a good feasible solution. Return  MRQ_HEURISTIC_SUCCESS or MRQ_OPTIMAL_SOLUTION if a feasible solution is found and false otherwise
        //virtual int getFeasibleSolution( MRQ_MINLPProb& prob, double* lx, double* ux, const double* xInitial, const double maxTime, const bool stopOnFirst, double* sol, double& obj);
        
        virtual int insideRun(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, const int thnumber, const double insideSolverMaxTime, double *nlx, double *nux);
        
        
        void printSubSolvers(const bool printMilp, const bool printNlp, const bool printGlobal);
        
        virtual bool tryUpdateBestSolution(const unsigned int threadNumber, const int n, double* sol, double objValue, const long unsigned int iter, const clock_t& clockStart, const double timeStart, const bool storeOnHistory);
        
        virtual bool tryUpdateBestSolutionWithLock(const unsigned int threadNumber, const unsigned int nthreads, MRQ_Mutex *mutex, const int n, double* sol, double objValue, const long unsigned int iter, const clock_t &clockStart, const double timeStart, const bool storeOnHistory);
        
        
        friend class MRQ_UserCallbacks;
        friend class MRQ_BBLCallbacks;
        
        friend int MRQ_insideRun( MRQ_Algorithm *alg, MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, const int thnumber, const double insideSolverMaxTime, double *nlx, double *nux );
    };





    class MRQ_LAAPointsStoring;
    class MRQ_GradientsEvaluation;
    class MRQ_BBLCallbacks;
    class MRQ_BonminHybrid;
    class MRQ_MILPSolverCallbackInterface;
    class MRQ_MILPCallbackData;
    class MRQ_MasterMILPProb;

    /*#if OPT_HAVE_CPLEX
        int CPXPUBLIC MRQ_labbLazyConstraintCplexCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
        
        int CPXPUBLIC MRQ_labbCplexBeforeSolveCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
        
        int CPXPUBLIC MRQ_labbCplexBranchingCallback(CPXCENVptr env, void *cbdata, int          wherefrom, void *cbhandle, int brtype, int sos, int nodecnt, int bdcnt, const int *nodebeg, const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p);
        
        void CPXPUBLIC MRQ_labbCplexDeleteNodeCallback(CPXCENVptr env, int wherefrom, void *cbhandle, int seqnum,       void *handle);
    #endif

    #if OPT_HAVE_GUROBI
        int __stdcall MRQ_bboaGurobiCallback( GRBmodel *model, void *cbdata, int where, void *usrdata);
    #endif

    #if OPT_HAVE_GLPK
        void MRQ_bboaGlpkCallback(glp_tree *T, void *info);
    #endif */


    class MRQ_LinearApproxAlgorithm : public MRQ_Algorithm
    {
        
    public:
        
        bool in_call_before_solve_callback_in_milp_bb; //just work for cplex
        bool in_call_branching_callback_in_milp_bb; //just work for cplex
        bool in_measure_nlp_time;
        bool in_set_obj_lower_bound_on_master_problem;
        bool in_set_quadratics_in_master_problem;
        
        
        MRQ_OBJ_LIN_STRATEGY  in_obj_linearisation_strategy;
        MRQ_CONSTR_LIN_STRATEGY in_constr_linearisation_strategy;
        MRQ_QUAD_APP_MASTER_STRATEGY in_quad_app_master_strategy;
        
        double in_eps_to_active_constr_to_linearisation;  //relative tolerance to decide if a constraint will be linearized
        
        
        /**** parameters to pseudo pruning (only to lazy constranints based algorithms by now) ****/
        
        bool in_try_pseudo_pruning_before_solving_relaxations; //warning: milp solvers can have no support for this. If you set and get some crashes, segfault etc, disable this option. 
        
        //unsigned int in_min_number_of_brachings_per_var_before_pseudo_pruning;
        
        MRQ_BB_PSEUDO_PRUNING_STRATEGY in_pseudo_pruning_strategy;
        
        double in_alpha_to_balance_estimative_in_pseudo_pruning; //in the interval [0 1]. Here, since we cannot have the lowest increase, we just multiply this parameter by the pseudo costs for each variable estimative.
        
        double in_absolute_convergence_tol_for_pseudo_pruning;
        
        double in_relative_convergence_tol_for_pseudo_pruning;
        
        double in_absolute_upper_bound_slack_for_pseudo_pruning;
        
        double in_relative_upper_bound_slack_factor_for_pseudo_pruning;
        
        /**** end of parameters to pseudo pruning (only to lazy constranints based algorithms by now) ****/
        
        
        //bool in_consider_bounds_to_aux_vars_in_oa_problem; //to consider lower and upper bounds to auxiliary variables in oa (milp) approximate problem. I was thinking it could help milp solver to go faster but sometimes the opposite occours.
        
        unsigned int out_number_of_obj_linears_saved_by_zu;
        unsigned int out_number_of_obj_linears_saved_by_infeasibility;
        unsigned int out_number_of_constr_linears_saved;
        
        long unsigned int out_number_of_nlp_probs_solved;
        long unsigned int out_number_of_milp_solver_iters;
        
        double out_cpu_time_on_box_to_constr_calculation;
        
        double out_clock_time_of_nlp_solving;
        double out_cpu_time_of_nlp_solving;
        
        
        /**** output to pseudo prunning (only to lazy constranints based algorithms by now) ****/
        
        long unsigned int out_number_of_pseudo_prunes_on_solving;
        long unsigned int out_number_of_pseudo_prunes_on_branching;
        
        /**** end of output to pseudo prunning (only to lazy constranints based algorithms by now) ****/
        
        
        
        MRQ_LinearApproxAlgorithm();
        
        virtual ~MRQ_LinearApproxAlgorithm();
        
        
        int addPointsToLinearisation(const int npoints, const int dim,  double** points);
        
        int addPointsToLinearisation(const int dim, MRQ_SolutionHistory& points);
        
        int addPointToLinearisation( const int dim,  const double* point);
        
        void copyParametersFrom(const MRQ_LinearApproxAlgorithm &source);
        
        int deletePointsToLinearisation(void);
        
        virtual bool isLinearApproximationAlgorithm() const override;
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetOutput() override;
        
        virtual void resetParameters() override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        virtual int setDoubleParameter(const char *name, const double value) override;
        
        //that method should be used to set enumeration parameters
        virtual int setStringParameter(const char *name, const char *value) override;
        
        
        
    protected:
        
        /********* for lazy constraints based algorithms *********/
        unsigned int nthreads;
        unsigned int nthreads_lazy;
        /********* end for lazy constraints based algorithms *********/
        
        /********** for pseudo pruning *******/
        
        int pp_nvars;
        int pp_nI;
        
        unsigned int pp_nthreads;
        
        int *pp_intVars;
        
        long unsigned int *pp_nBranchs; //number of branching for each variable
        
        long unsigned int *pp_counter1; //pseudo pruning counter for before solve relacation
        long unsigned int *pp_counter2; //pseudo pruning counter
        
        double *pp_lowestBoundPseudoPruned; //to store the lowest lower bound from pseudo pruned nodes
        
        double **pp_lupcosts; //to get lower and upper pseudo costs (we need an array per thread)
        double **pp_sols; //to get solution values
        
        int allocateMemoryForPseudoPruning(unsigned int nthreads, unsigned int nvars, unsigned int nI);
        
        void deallocateMemoryForPseudoPruning();
        
        /********** end for pseudo pruning *******/
        
        
        int nPoints;  //Number of user points to add linearisations.
        double **points;	//User points to add linearisations
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        virtual int _algorithmInitialization() override;
        
        virtual void algorithmFinalization( const int nthreads, MRQ_MINLPProb& prob, double* &lx, double* &ux ) override;
        
        
        int setLazyConstraintCallbackOnMilpSolver( MRQ_LPSolver *master, void *callbackData);
        
        int setBeforeSolveCallbackOnMilpSolver( MRQ_LPSolver *master, void *callbackData, bool enforceOriginalIndices = true);
        
        int setBranchingCallbackOnMilpSolver( MRQ_LPSolver *master, void *callbackData, bool enforceOriginalIndices = true);
        
        
        int setDeleteNodeCallbackOnMilpSolver( MRQ_LPSolver *master, void *callbackData, bool enforceOriginalIndices = true);
        
        
        
        int setMasterProblemBase( MRQ_MasterMILPProb *masterMilp, const int thnumber,  MRQ_MINLPProb& prob, const int solver, const bool setLinearObj, const bool setQuad, const bool setVarTypes, const double* lx, const double* ux, const int nauxvars, MRQ_GeneralSolverParams* params, bool *auxConstrEval, MRQ_LAAPointsStoring *laps, const int nthreads );
        
        //method to be implemented by BB based on OA and ECP adopting lazy constraints
        virtual int solverCallbackLazyConstraints( MRQ_MILPCallbackData &data);
        
        //method to be implemented by BB based on OA and ECP toenable a before solve callback
        virtual int solverCallbackBeforeSolve( MRQ_MILPCallbackData &data);
        
        virtual int solverCallbackBranching( MRQ_MILPCallbackData &data);
        
        int addLazyConstraintsLinearizationOnSolution( const unsigned int thnumber, MRQ_MILPSolverCallbackInterface *callbackSolver, MRQ_MINLPProb &prob, MRQ_GradientsEvaluation &gradEval, const bool incQuadsInMaster, const bool linearizeObj, const bool *auxConstrEval, const double *psol, const double *objSol, const double *pconstr, const double *masterConstr, const int *indices, const double *plc, const double *puc, double *auxVars, bool *auxConstrEval2);
        
        
        friend class MRQ_BBLCallbacks;
        friend class MRQ_BonminHybrid;
        
    public:
        
        int getNPoints(); //return number of points to add linearisations
        
        
        #if OPT_HAVE_CPLEX
            friend int CPXPUBLIC MRQ_labbLazyConstraintCplexCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
            
            friend int CPXPUBLIC MRQ_labbCplexBeforeSolveCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
            
            friend int CPXPUBLIC MRQ_labbCplexBranchingCallback(CPXCENVptr env, void *cbdata, int          wherefrom, void *cbhandle, int brtype, int sos, int nodecnt, int bdcnt, const int *nodebeg, const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p);
            
            friend void CPXPUBLIC MRQ_labbCplexDeleteNodeCallback(CPXCENVptr env,
           int wherefrom, void *cbhandle, int seqnum,       void *handle);
            
        #endif
        
        #if OPT_HAVE_GUROBI
            friend int __stdcall MRQ_bboaGurobiCallback( GRBmodel *model, void *cbdata, int where, void *usrdata);
        #endif
        
        #if OPT_HAVE_GLPK
            friend void MRQ_bboaGlpkCallback(glp_tree *T, void *info);
        #endif
        
    };

    


    class MRQ_ExtCutPlan : public MRQ_LinearApproxAlgorithm
    {
        
    public:
        //bool in_binary_cut_when_int_fixed_master_infeasible; // when we fix integer variables in master problem and problem gets infeasible, use the binary cut. Otherwise, we try solve a kind of feasibility problem...
        //bool in_fix_int_vars_to_try_improve_cut; //fix integer variables in master problem to try improve the solution and cuts...
        
        bool in_delete_intermediate_linearizations_in_each_iteration; //only keep linearization on first and last point in each iteration
        
        bool in_refine_final_solution_using_nlp;
        
        unsigned int in_max_number_of_fixed_relax_master_problem_solved_per_iteration;
        
        unsigned int in_subiter_printing_frequency;
        
        
        bool out_all_subproblems_py_completely_solved; //true if all subproblems Py (where we fix integer variable) were completely solved
        long unsigned int out_number_of_subproblems_py_completely_solved; //number of subproblems Py (where we fix integer variable) completely solved;
        
        long unsigned int out_number_of_master_relaxation_solved;
        
        #if MRQ_DEBUG_MODE
            long unsigned int out_number_of_integer_solution_repetitions; //number of repetions of integer solutions (the single repetion of integer soltuion marked as optimal is disconsidered)
            MRQ_SolutionHistory out_integer_solution_history; //this history is to save integer solution to see if integer solutions are repeting. Here, we count even optimal integer solution
        #endif
        
        
        MRQ_ExtCutPlan();
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        virtual void resetOutput() override;
        
        virtual void resetParameters() override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;

        virtual ~MRQ_ExtCutPlan();
    };
    
    
    class MRQ_OuterApp : public MRQ_LinearApproxAlgorithm
    {
        
    public:
        
        bool in_binarie_cut_when_nlp_infeasible;
        bool in_round_first_nlp_relaxation_solution;
        bool in_use_first_nlp_relaxation;
        
        
        MRQ_OuterApp();
        
        virtual ~MRQ_OuterApp();
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        //virtual void resetOutput() override;
        
        virtual void resetParameters() override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        
    protected:
        
        void updateBestSol( const int n, double* newsol, const double* dualSol, const double objValue, optsolvers::OPT_LPSolver* master, muriqui::MRQ_LAAPointsStoring* laps, const bool linearizeObj, muriqui::MRQ_QUAD_APP_MASTER_STRATEGY quadAppStrategy, muriqui::MRQ_IntegratedHessian* intHess, double* auxVars, const clock_t& clockStart, const int timeStart, const long unsigned int iter);
        
    };






    #if OPT_HAVE_CPLEX
        int MRQ_bblpecpCplexCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
    #endif

    #if OPT_HAVE_GUROBI
        int MRQ_bblpecpGurobiCallback( GRBmodel *model, void *cbdata, int where, void *usrdata);
    #endif

    #if OPT_HAVE_GLPK
        void MRQ_bblpecpGlpkCallback(glp_tree *T, void *info);
    #endif


    class MRQ_LPBBExtCutPlan : public MRQ_ExtCutPlan
    {
    public:
        
        //bool in_fix_int_vars_to_try_improve_cut; //fix integer variables in master problem to try improve the solution and cuts...
        
        
        MRQ_LPBBExtCutPlan();
        
        virtual ~MRQ_LPBBExtCutPlan();
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        //virtual void resetOutput() override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        
    protected:
        
        MRQ_Mutex SEMAPH_history;
        
    private:
        
        MRQ_Mutex SEMAPH_updtSol;
        MRQ_Mutex SEMAPH_updtOut;
        
        virtual int solverCallbackLazyConstraints( MRQ_MILPCallbackData &data) override;
        
    };




    class MRQ_ExtSupHypPlane : public MRQ_LinearApproxAlgorithm
    {
        
    public:
        
        //bool in_delete_linearizations_from_int_vars_fixing;
        bool in_fix_int_vars_to_try_improve_cut; //fix integer variables in master problem to try improve the solution and cuts..
        bool in_linearize_on_interior_sol;
        bool in_try_solve_interior_problem_if_cont_relax_fail;
        
        unsigned int in_max_lp_subiters;
        unsigned int in_max_lp_subiter_to_improve_obj_app;
        
        MRQ_ESHP_INTERIOR_POINT_STRATEGY in_interior_point_strategy;
        MRQ_CONSTR_LIN_STRATEGY in_lp_constr_linearisation_strategy;
        
        
        double in_absolute_tol_to_check_previous_sol;
        double in_cont_relax_absolute_convergence_tol;
        double in_cont_relax_relative_convergence_tol;
        
        double in_delta_to_inc_eps_to_active_constraints_to_linearization;
        
        //double in_eps_lp;
        //double in_eps_milp;
        double in_eps_to_interior_sol; //eps to consider the interior solution valid e.g. 0.01
        double in_eps_to_line_search;
        
        double in_eps_to_enforce_interior_sol_on_cont_relax_sol;
        
        double in_relative_tol_to_check_previous_sol;
        
        
        unsigned int out_number_of_lp_iterations;
        
        
        MRQ_ExtSupHypPlane();
        
        virtual ~MRQ_ExtSupHypPlane();
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetOutput() override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        virtual int setDoubleParameter(const char *name, const double value) override;
        
        //that method should be used to set enumeration parameters
        virtual int setStringParameter(const char *name, const char *value) override;
        
    protected:
        
        //bool updateBestSol( const int n, double objValue, double* sol, const clock_t& clockStart, const int timeStart, const long unsigned int iter );
    };



    class MRQ_LPBBExtSupHypPlane : public MRQ_ExtSupHypPlane
    {
    public:
        
        long unsigned int out_number_of_master_relaxation_solved;
        
        MRQ_LPBBExtSupHypPlane();
        
        virtual ~MRQ_LPBBExtSupHypPlane();
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        virtual void resetOutput() override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
    protected:
        
        MRQ_Mutex SEMAPH_history;
        
    private:
        
        MRQ_Mutex SEMAPH_updtSol;
        
        
        virtual int solverCallbackLazyConstraints( MRQ_MILPCallbackData &data) override;
    };


    class MRQ_LPNLPBBOuterApp :public MRQ_LinearApproxAlgorithm
    {
    public:
        
        
        bool in_use_first_nlp_relaxation;
        bool in_binarie_cut_when_nlp_infeasible;
        bool in_linearize_obj_in_nl_feas_solutions;
        
        
        MRQ_LPNLPBBOuterApp();
        
        virtual ~MRQ_LPNLPBBOuterApp();
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
    protected:
        
        MRQ_Mutex SEMAPH_history;
        
    private:
        
        MRQ_Mutex SEMAPH_updtSol;
        
        
        
        virtual int solverCallbackLazyConstraints( MRQ_MILPCallbackData &data) override;
        
        
    };


    class MRQ_BonminHybrid : public MRQ_LPNLPBBOuterApp
    {
    public:
        
        long unsigned int in_out_app_max_iterations;
        unsigned int in_nlp_relaxation_solving_frequence;
        double in_outer_app_max_time;
        double in_outer_app_max_cpu_time;
        
        
        MRQ_OuterApp in_outer_app;
        
        
        
        MRQ_BonminHybrid();
        
        virtual ~MRQ_BonminHybrid();
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        virtual int setDoubleParameter(const char *name, const double value) override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        
    private:
        
        //MRQ_Mutex SEMAPH_incIterCount;
        
        //long unsigned int iterCount; //to count iterations to apply continuous realaxtion
        
        int solverCallbackCutToBonminHyb(MRQ_MILPCallbackData &data);
        
        
        #if OPT_HAVE_CPLEX
            friend int CPXPUBLIC MRQ_bonminHybCplexCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);
        #endif
        
    };



    typedef branchAndBound::BBL_MTNodeListManager MRQ_MTNodeListManager;
    //typedef branchAndBound::BBL_PruneCounter MRQ_PruneCounter;

    class MRQ_PruneCounter : public branchAndBound::BBL_PruneCounter
    {
    public:
        long unsigned int pseudopruning1; //to count pseudo pruning when branching variable(s) is already choosen
        long unsigned int pseudopruning2; //to count pseudo pruning before choosing braching variable(s)
        
        MRQ_PruneCounter():branchAndBound::BBL_PruneCounter(){
                reset();
        }
        
        void reset()
        {
            branchAndBound::BBL_PruneCounter::reset();
            pseudopruning1 = 0;
            pseudopruning2 = 0;
        }
    };



    class MRQ_ThreadData;
    class MRQ_Points;
    class MRQ_PseudoCostCalc;
    //class MRQ_VarBoundsUpdater;
    class MRQ_GlobalCutList;
    class MRQ_NodeGenerator;
    class MRQ_UserCallbacks;



    class MRQ_BranchAndBound : public MRQ_Algorithm
    {
        
    public:
        
        //bool in_branch_even_integer_sol; //that only will work if user set a callback to generate node or choose index to branch...
        
        bool in_calculate_pseudo_cost_average_above_error_estimative;
        
        bool in_call_after_bb_loop_callback;
        bool in_call_before_bb_loop_callback;
        bool in_call_before_solving_relax_callback;
        bool in_call_after_solving_relax_callback;
        bool in_call_generate_root_node_callback;
        bool in_count_total_prunes;
        bool in_consider_relax_infeas_if_solver_fail;
        bool in_consider_relax_infeas_if_solver_get_unbounded;
        
        bool in_only_apply_pseudo_pruning_on_fixed_integer_vars;
        
        bool in_repeat_strong_branching_if_get_bounds_updating;
        
        bool in_stop_multibranch_after_first_bound_prune; //strategy like constraint branching and multibranchng will be avoided after the first prune by bound
        bool in_use_dual_obj_to_bound_prunning;
        bool in_use_early_branching;
        
        bool in_use_feas_heuristic_diving;
        bool in_use_feas_heuristic_fp;
        bool in_use_feas_heuristic_oafp;
        bool in_use_feas_heuristic_rens;
        
        
        bool in_use_outer_app;
        bool in_use_outer_app_as_heuristic;
        
        //bool in_use_parent_primal_sol_as_initial_sol;
        //bool in_use_parent_dual_sol_as_initial_sol;
        //bool in_use_round_heuristic;
        
        
        int in_bounds_updt_frequency;
        int in_feas_heuristic_frequency;
        int in_igma2_frequency;
        int in_lists_reorganization_frequency;
        unsigned int in_max_number_of_branchings_in_constraint_branching;
        unsigned int in_max_tree_level_to_count_prunes_by_level;
        unsigned int in_min_number_of_bound_prunes_per_var_before_pseudo_pruning;
        int in_number_of_branching_vars; //if <= 0, we will set it automatically
        int in_number_of_node_sublists;
        unsigned int in_number_of_pos_iters_to_update_var_bounds; //when a feasible solution is found, we update var bounds after this number of iterations (if that strategy is enable, of course). It avoids try update var bounds too many times if several feasible solutions are found in a small interval of iterations
        //int in_number_of_pseudo_costs_rely;
        
        
        int in_outer_app_frequence; //Number of iterations beetwen two calling to outer approximation algorithm applied to original problem
        int in_outer_app_subprob_frequence; //Number of iterations beetwen two calling to outer approximation algorithm applied to subproblems in the tree
        int in_rounding_heuristic_call_iter_frequence;
        
        long int in_seed_to_random_numbers;
        
        MRQ_BB_BOUND_LIN_UPDT_STRATEGY in_bounds_linear_updt_strategy;
        MRQ_BB_BRANCH_STRATEGY in_branching_strategy;
        MRQ_BB_BRANCH_STRATEGY in_variable_choosing_strategy_to_constraint_branching;
        MRQ_BB_CONSTRAINT_BRANCH_STRATEGY in_constr_branching_strategy;
        MRQ_BB_EXP_STRATEGY in_exp_strategy;
        MRQ_BB_INT_HEURISTICS_STRATEGY in_igma2_strategy;
        MRQ_BB_INT_HEURISTICS_STRATEGY in_int_feas_heurs_strategy;
        MRQ_LP_SOLVER in_linear_bounds_updt_solver;
        
        MRQ_BB_PARENT_NODE_BOUNDS_STORAGE_STRATEGY in_parent_node_bounds_storage_strategy;
        MRQ_BB_PARENT_SOL_STORING_STRATEGY in_parent_sol_storing_strategy;
        MRQ_ROUNDING_STRATEGY in_rounding_heuristic_strategy;
        
        MRQ_BB_PSEUDO_PRUNING_STRATEGY in_pseudo_pruning_strategy;
        
        
        double in_feas_heuristic_max_time;
        double in_bounds_updt_factor_to_subgroup; //only if in_linear_bounds_updt_strategy is MRQ_LBUS_SUBGROUP
        double in_nlp_relative_opt_tol_on_early_branching;
        double in_nlp_relative_opt_tol_to_int_sol_on_early_branching;
        
        double in_nlp_relative_primal_tol_on_early_branching;
        double in_nlp_relative_primal_tol_to_int_sol_on_early_branching;
        
        double in_nlp_relative_dual_tol_on_early_branching;
        double in_nlp_relative_dual_tol_to_int_sol_on_early_branching;
        
        double in_absolute_convergence_tol_for_pseudo_pruning;
        double in_relative_convergence_tol_for_pseudo_pruning;
        
        double in_absolute_upper_bound_slack_for_pseudo_pruning;
        double in_relative_upper_bound_slack_factor_for_pseudo_pruning;
        
        double in_alpha_to_balance_estimative_in_pseudo_pruning; // alpha to balance the estimative to objective increasing in pseudo pruning [0 1). If alpha is 1, we choose as estimative the statistical increasing mean value in bound pruning. If alpha is 0, we choose as estimative the minimum value already seen in bound pruning.
        
        double in_outer_app_subprob_time;
        double in_outer_app_time;
        double in_pseudo_cost_mu;
        double in_pseudo_cost_mu_to_variable_choosing_on_constraint_branching;
        
        double in_integer_neighborhood_factor_to_ssr;
        
        /********* begin of igma2 parameters **********/
        
        bool in_igma2_set_max_dist_constr_on_bin_vars;
        bool in_igma2_set_special_gap_min_solver_params;
        bool in_igma2_solve_local_search_problem_even_on_non_int_sol;
        
        
        //MRQ_NLP_SOLVER in_igma2_gap_min_solver;
        MRQ_IGMA2_NEIGHBORHOOD_STRATEGY in_igma2_neighborhood_strategy;
        
        double in_igma2_factor_to_max_dist_constr;
        double in_igma2_factor_to_max_dist_constr_on_bin_vars;
        double in_igma2_percentual_to_rectangular_neighborhood;
        
        /********* end of igma2 parameters **********/
        
        
        
        MRQ_GeneralSolverParams *in_OAmilpParams;
        MRQ_GeneralSolverParams *in_OAnlpParams;
        
        bool out_nlp_failure_in_some_relaxation;
        
        int out_seed_to_random_numbers;
        unsigned int out_number_of_strong_branching_calculations_to_pseudo_costs;
        long unsigned int out_number_of_constraint_branchings;
        long unsigned int out_number_of_open_nodes;
        long unsigned int out_number_of_iters_having_wrong_lower_bounds;
        
        
        double out_pseudo_cost_average_above_error_estimative;  //a estimative of standard deviation in objective increasing by pseudo-costs.
        
        MRQ_PruneCounter out_prune_counter;
        MRQ_PruneCounter *out_prune_counter_by_level;
        
        
        
        
        
        MRQ_BranchAndBound();
        
        ~MRQ_BranchAndBound();
        
        
        virtual void resetOutput() override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        
        void desallocatePruneCountersByLevel(); //user can free the prune counter by level
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        //warning: overwritte array already stored. DCS0: constraints that will be disable when a variable is fix in a value 
        int setDCS0Array(const int size, const MRQ_DynConstrSetUnity *dcs);
        
        //warning: overwritte array already stored. DCS1: constraints that will be enable when a variable is fix in a value 
        int setDCS1Array(const int size, const MRQ_DynConstrSetUnity *dcs);
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        virtual int setDoubleParameter(const char *name, const double value) override;
        
        //that method should be used to set enumeration parameters
        virtual int setStringParameter(const char *name, const char *value) override;
        
        
    protected:
        
        friend class MRQ_UserCallbacks;
        friend class MRQ_BBLCallbacks;
        
        MRQ_BBLCallbacks *_bblCallbacks;
        
        
        
        int allocatePruneCountersByLevel(const int nlevels);
        
        
        
        virtual bool tryUpdateBestSolution(const unsigned int threadNumber, const int n, double* sol, double objValue, const long unsigned int iter, const clock_t& clockStart, const double timeStart, const bool storeOnHistory) override;
        
        virtual bool tryUpdateBestSolution(const unsigned int threadNumber, const int n, double* sol, double objValue, const long unsigned int iter, const clock_t& clockStart, const double timeStart, const bool storeOnHistory, const bool addSolToOALin); 
        
        virtual int getBestSolutionCopy(const int n, double *solution, double &fsolution);
        
    };





    class MRQ_ContinuousRelax : public MRQ_Algorithm
    {
        
    public:
        
        bool in_set_integer_vars_as_integers;
        
        bool out_nlp_feasible_solution;
        int out_original_solver_return_code;
        double *out_constraint_values;
        double *out_dual_sol;
        
        
        MRQ_ContinuousRelax();
        
        virtual ~MRQ_ContinuousRelax();
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetParameters() override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        virtual void resetOutput() override;
        
        
    protected:
        
        
        int allocateDualAndConstraintsArrays( const int ndual, const int m);
        
        void desallocateDualAndConstraintsArrays();
    };



    class MRQ_IGMA0 : public MRQ_Algorithm
    {
    public:
        
        bool in_declare_relax_infeas_if_solver_fail;
        bool in_use_integers_vars_on_gap_min_problem;
        bool in_stop_heuristcs_on_first_solution;  //try heuristics only until the first feasible solution be found
        bool in_use_heuristcs;
        
        int in_number_of_threads; //to global solver used
        
        
        MRQ_IGMA0_EPS_STRATEGY in_eps_updt_strategy;
        
        MRQ_GLOBAL_SOLVER in_global_solver;
        
        
        double in_absolute_obj_cut_epsilon;
        double in_decrease_fator_integer_tol;
        double in_relative_obj_cut_epsilon;
        double in_max_time_on_heuristcs;
        double in_alpha; //parameter to balance objective cut on MRQ_IGMA1_EPS_BIN_SEARCH strategy
        
        
        
        
        MRQ_IGMA0();
        
        virtual ~MRQ_IGMA0();
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* globalSolverParams = NULL, MRQ_GeneralSolverParams* subGlobalSolverParams= NULL) override;
        
        virtual int setDoubleParameter(const char *name, const double value) override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        virtual int setStringParameter(const char *name, const char *value) override;
    };


    class MRQ_IGMA1BBCallbacks;

    class MRQ_IGMA1 : public MRQ_Algorithm
    {
        
    public:
        
        bool in_adopt_obj_cut; //only for heuristic mode...
        bool in_consider_relax_infeas_if_solver_fail;
        bool in_enable_gap_min_solver_premature_stoping;
        bool in_heuristic_mode;
        bool in_reorganize_lists;
        bool in_set_special_gap_min_solver_params;
        bool in_set_gap_exp_on_constr; //set gap expressio in constraints instead of objective function. Note, a auxiliary varialbe will be added..
        bool in_set_gap_ubound_constr;
        //bool in_stop_on_first_sol;
        bool in_use_random_initial_sols;
        bool in_use_general_feas_heuristics;
        
        
        int in_lists_reorganization_frequency;
        unsigned int in_max_improvments_of_best_sol;
        unsigned int in_max_nonimprovment_integer_sols; //only for in_heuristic_mode = true and in_adopt_obj_cut = false
        unsigned int in_min_number_of_iters_on_gap_min_premature_stop;
        
        
        long int in_seed_to_random_numbers;
        
        
        MRQ_BB_CONSTRAINT_BRANCH_STRATEGY in_constr_branching_strategy;
        MRQ_BB_EXP_STRATEGY in_exp_strategy;
        
        MRQ_NLP_SOLVER in_gap_min_solver;
        MRQ_IGMA_GAP_MIN_OBJ_STRATEGY in_gap_min_obj_strategy;
        
        
        
        double in_absolute_obj_cut_eps;
        double in_relative_obj_cut_eps;
        
        
        double in_factor_to_increase_abs_obj_cut_eps;
        double in_factor_to_increase_rel_obj_cut_eps;
        
        double in_factor_to_increase_abs_obj_cut_eps_on_zero;
        
        double in_absolute_obj_tol_to_active_obj_cut; //absolute tolerance to check if objective cut is active (heuristic mode)
        double in_relative_obj_tol_to_active_obj_cut; //relative tolerance to check if objective cut is active (heuristic mode)
        
        double in_factor_to_min_gap_obj_by_avg_gap;
        //double in_integer_tol_on_gap_min_premature_stop;
        
        double in_lower_bound_to_random_sol; //just to replace -Infinity or a lower value to a variable lower bound
        double in_upper_bound_to_random_sol; //just to replace Infinity or a upper value to a variable upper bound
        
        double in_min_improv_to_obj_on_gap_min_premature_stop;
        double in_min_improv_to_infeas_on_gap_min_premature_stop;
        
        
        
        
        
        unsigned long int out_number_of_open_nodes;
        unsigned long int out_number_of_inner_iterations;
        unsigned long int out_number_of_prunes_by_obj_cut_active;
        
        
        MRQ_IGMA1();
        
        virtual ~MRQ_IGMA1();
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetOutput() override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* gapMinSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        virtual int setDoubleParameter(const char *name, const double value) override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        virtual int setStringParameter(const char *name, const char *value) override;
        
        
    protected:
        
        
        friend class MRQ_IGMA1BBCallbacks;
    };



    class MRQ_IGMA2 : public MRQ_BranchAndBound
    {
        
    public:
        
        bool in_set_max_dist_constr_on_bin_vars;
        bool in_set_special_gap_min_solver_params;
        bool in_solve_local_search_problem_even_on_non_int_sol;
        
        
        MRQ_IGMA2_NEIGHBORHOOD_STRATEGY in_neighborhood_strategy;
        
        double in_factor_to_max_dist_constr;
        double in_factor_to_max_dist_constr_on_bin_vars;
        double in_percentual_to_rectangular_neighborhood;
        
        MRQ_UserCallbacks *in_user_callbacks; //some day we could implement some callback interfaces on this algorithm. The original in_user_callbacks form MRQ_Algorithm is used for us to implement igma2.
        
        
        
        bool out_feas_sol_on_gap_min_problem;
        bool out_feas_sol_found_by_igma2_procedure; //if is false, the feasible solution was found by some b&b procedure out of igma2
        
        
        MRQ_IGMA2();
        
        virtual ~MRQ_IGMA2();
        
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetOutput() override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* gapMinSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        virtual int setDoubleParameter(const char *name, const double value) override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        virtual int setStringParameter(const char *name, const char *value) override;
    };




    class MRQ_Heuristic : public MRQ_Algorithm
    {
        
    public:
        
        bool in_solve_nlp_as_local_search_at_end;
        
        bool in_use_random_seed_to_random_numbers;
        
        long int in_seed_to_random_numbers;
        
        long int out_seed_to_random_numbers;
        
        MRQ_Heuristic();
        
        virtual ~MRQ_Heuristic();
        
        
        
        void copyParametersFrom(const MRQ_Heuristic &source);
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetParameters() override;
        
        virtual void resetOutput() override;
        
        //virtual bool checkTerminationCriterions(const int threadNumber, const double zl, const double zu, long unsigned int iter, const double timeStart, const clock_t clockStart, int& retCode);
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
    };



    class MRQ_FeasibilityPump : public MRQ_Heuristic
    {
        
    public:
        
        bool in_set_linear_obj_term_on_bin_vars_at_nlp;
        bool in_set_norm1_on_nlp;
        
        unsigned int in_last_iters_considered_to_cycles;
        unsigned int in_max_cycle_subiters;
        
        double in_lower_bound_to_pi;
        double in_upper_bound_to_pi;
        
        
        MRQ_FeasibilityPump();
        
        virtual ~MRQ_FeasibilityPump();
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        virtual int setDoubleParameter(const char *name, const double value) override;
    };



    class MRQ_OAFeasibilityPump : public MRQ_LinearApproxAlgorithm//, public MRQ_Heuristic
    {
    public:
        
        bool in_set_norm1_on_nlp; //if false, we use an approximation of norm2
        bool in_set_linear_obj_term_on_bin_vars_at_nlp; //if true, nlp consider a linear term for binary variables. Otherwise, a quadratic term will be used. The objective measures the distance from integer solution
        bool in_solve_nlp_as_local_search_at_end;
        //int in_max_number_of_gap_equals;
        
        
        MRQ_OAFeasibilityPump();
        
        virtual ~MRQ_OAFeasibilityPump();
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
    };



    class MRQ_Diving : public MRQ_Heuristic
    {
    public:
        
        bool in_consider_relax_infeas_if_solver_fail;
        
        MRQ_DIVE_SELECT_STRATEGY in_dive_selec_strategy;
        
        double in_percentual_of_add_var_fixing;
        
        
        
        MRQ_Diving();
        
        virtual ~MRQ_Diving();
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetParameters() override;
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        virtual int setDoubleParameter(const char *name, const double value) override;
        
        //that method should be used to set enumeration parameters
        virtual int setStringParameter(const char *name, const char *value) override;
        
        
    private:
        
        
        MRQ_RETURN_CODE selectFracVariable(muriqui::MRQ_MINLPProb& prob, const int* intVars, const double* sol, const int* nrowsVars, double* auxVars, bool& intSol, unsigned int& index, bool& roundUp);
    };

    class MRQ_HeuristicExecutor;

    class MRQ_RENS : public MRQ_Heuristic 
    {
    public:
        
        bool in_apply_only_heuristics_on_subproblem; //apply several heuristcs
        bool in_solve_continuous_relaxation_on_subproblem; //solve continuous relaxation for the subproblem, considering the reduced feasible region. Since we solve continuous relaxation for the original problem, we coud skip this step, but maybe the new constraints can help preprocessor to reduce bounds and get a better solution
        
        bool in_stop_on_first_sol; //stop execution in the first feasible solution founded
        
        double in_continuous_neighborhood_factor;
        double in_integer_neighborhood_factor; //factor to define local branching constraint (only if you do not use MRQ_RENS_NS_ORIGINAl as in_neighborhood_strategy)
        
        MRQ_ALG_CODE in_algorithm_to_solve_subproblem;
        
        MRQ_SUBPROBLEM_NEIGHBORHOOD_STRATEGY in_neighborhood_strategy;
        
        MRQ_Algorithm *in_algorithm_object_to_solve_subproblem; //put here a pointer to algorithm that you would like to be used to solve subproblem 
        
        MRQ_HeuristicExecutor *in_heuristic_exec_object_to_solve_subproblem; //put here a pointer to heuristic executor that you would like to be used to solve subproblem  
        
        
        MRQ_ALG_CODE out_algorithm_to_solve_subproblem;
        
        
        MRQ_RENS();
        
        virtual ~MRQ_RENS();
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetParameters() override;
        
        virtual void resetOutput() override;
        
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, MRQ_GeneralSolverParams* algParams);
        
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        virtual int setDoubleParameter(const char *name, const double value) override;
        
        //that method should be used to set enumeration parameters
        virtual int setStringParameter(const char *name, const char *value) override;
        
        
    };
    
    
    class MRQ_StructuredStochasticRounding : public MRQ_Heuristic
    {
        
    public:
        
        //bool in_preprocess_after_handling_constraints;
        bool in_preprocess_after_variable_rounding;
        bool in_random_order_to_threat_classes;
        bool in_random_order_to_threat_constraints_in_each_class;
        bool in_solve_continuous_relaxation;
        bool in_solve_minlp_as_local_search;
        bool in_stop_local_search_solving_on_first_improvment_solution;
        
        unsigned int in_max_number_of_improvments;
        
        
        double in_integer_neighborhood_factor; //factor to define local branching constraint (only if you do not use MRQ_SNS_ORIGINAL as in_neighborhood_strategy)
        double in_continuous_neighborhood_factor; //factor to define local branching constraint (only if you use MRQ_SNS_EUCLIDEAN_NEIGHBORHOOD as in_neighborhood_strategy). Radius of the neighborhood will be defined line in_continuous_neighborhood_factor*nC*alpha, where nC is the number of continuous variables and alpha is the greatest absolute value in continuous initial solution
        
        double in_relative_convergence_tol_to_local_search;
        
        double in_min_probability_to_round;
        
        MRQ_SUBPROBLEM_NEIGHBORHOOD_STRATEGY in_neighborhood_strategy;
        MRQ_SSR_PREPROCESSING_STRATEGY in_preprocessing_point_strategy;
        MRQ_SSR_VARIABLES_ADDITIONAL_FIXING_STRATEGY in_additional_vars_fixing_strategy;
        
        MRQ_SSR_VARIABLE_BOUNDS_UPDT_STRATEGY in_rounding_var_bounds_updt_strategy;
        
        MRQ_SSR_CONTINUOUS_RELAXATION_STRATEGY_TO_STOCHASTIC_ROUNDING in_cont_relax_strategy_to_stoch_rounding;
        
        MRQ_ALG_CODE in_local_search_algorithm;
        
        
        MRQ_GeneralSolverParams* in_alg_params;
        
        MRQ_Algorithm *in_alg; //if user put a pointer here, we use to aply algorithm. User has still the ownership of the pointer
        
        
        bool out_vars_fixed_by_stoch_rounding;
        
        double out_cpu_time_at_nlp_integer_fixed_sol;
        double out_obj_at_nlp_integer_fixed_sol;
        double out_obj_at_first_sol;
        
        MRQ_ALG_CODE out_local_search_alg_code;
        
        
        
        MRQ_StructuredStochasticRounding();
        
        virtual ~MRQ_StructuredStochasticRounding();
        
        virtual int checkAlgorithmRequirements(MRQ_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL) override;
        
        virtual void printParameters(std::ostream &out = std::cout) const override;
        
        virtual void resetParameters() override;
        
        virtual void resetOutput() override;
        
        
        virtual int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams = NULL, MRQ_GeneralSolverParams* nlpSolverParams = NULL) override;
        
        int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, MRQ_GeneralSolverParams* algParams);
        
        
        virtual int setIntegerParameter(const char *name, const long int value) override;
        
        virtual int setDoubleParameter(const char *name, const double value) override;
        
        //that method should be used to set enumeration parameters
        virtual int setStringParameter(const char *name, const char *value) override;
        
    };

    
    class MRQ_RensExecutor;

    class MRQ_HeuristicExecutor
    {
        
    private:
        
        int ai; //algorithm index iterator...
        
        
        void __resetMyParameters();
        
        
        int __run( MRQ_MINLPProb& prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, double zu, double &obj, double *sol, bool cycleRunning, MRQ_Algorithm **successAlg, const bool insideRun, const int thnumber, const int nThreads, const double insideSolverMaxTime, double* nlx, double* nux);
        
        
        
        
    protected:
        
        //successAlg is an output argument and points to algorithm that gives the solution. It can be NULL
        int insideRun(MRQ_MINLPProb& prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams, double zu, double &obj, double *sol, bool cycleRunning, MRQ_Algorithm **successAlg, const int thnumber, const int nThreads, const double insideSolverMaxTime, double* nlx, double* nux);
        
        
        
        
    public:
        
        //const static int nAlgs = 6; //number of algorithms available
        
        
        MRQ_OAFeasibilityPump oafp;
        MRQ_Diving diving;
        MRQ_FeasibilityPump fp;
        MRQ_IGMA1 igma1;
        MRQ_IGMA2 igma2;
        MRQ_RENS rens;
        MRQ_StructuredStochasticRounding ssr;
        
        const std::vector< std::pair< bool* const, MRQ_Algorithm* const> > algs = { {&in_use_oafp, &oafp}, {&in_use_fp, &fp},  {&in_use_diving, &diving},  {&in_use_igma2, &igma2}, {&in_use_igma1, &igma1}, {&in_use_rens, &rens}, {&in_use_ssr, &ssr}  };
        
        
    public:
        
        bool in_stop_on_first_sol;
        
        bool in_use_oafp;// = use_flags[2]; 
        bool in_use_fp;// = use_flags[0];
        bool in_use_diving;// = use_flags[3];
        bool in_use_igma1;// = use_flags[1];
        bool in_use_igma2;
        bool in_use_rens;
        bool in_use_ssr;
        
        double in_max_total_cpu_time;
        double in_max_total_clock_time; //note, that attribute does not affect each algorithm's maximum cpu time, just avoid new algorithms be applied if the maximum time is reached.
        
        
        
        MRQ_HeuristicExecutor();
        
        
        //virtual void printParameters(std::ostream &out = std::cout) const;
        void printParameters(std::ostream &out = std::cout) const;
        
        void resetParameters();
        
        
        
        bool* getAlgUseFlagPointer( const int number );
        
        
        //number ranges from 0 to nAlgs-1. That method is usefull to add new heuristics classes in the future performing minimal changes in the class code. You only have to change nAlgs and assigned a number to new heuristics
        MRQ_Algorithm* getAlgPointer(int number);
        
        int getNumberOfAlgs() const;
        
        //return True if the pointer alg is one of the objects in this 
        bool hasAlgPointer(const MRQ_Algorithm *alg);
        
        //set initial soultion for all algorithms. If setParameter is true, this methdo set in_use_initial_solution;
        int setInitialSolution(const int n, const double *sol);
        
        
        //set all heuristics to use solvers
        void setSolvers(MRQ_MILP_SOLVER milpSolver,  MRQ_NLP_SOLVER nlpSolver);
        
        
        void setMaxCPUTimes( const double maxTime);
        
        
        //set all heuristics to use the same maximum cpu time
        void setMaxTimes(const double maxTime);
        
        
        void setMaxIters(const long unsigned int maxIters);
        
        void setNumberOfThreads(const unsigned int numberOfThreads);
        
        void setPrintLevels(const int level);
        
        
        //if cycleRunning is true, run start from the next algorithm no used in the last calling...
        int run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams *milpSolverParams, MRQ_GeneralSolverParams *nlpSolverParams, double zu, double &obj, double *sol = NULL, bool cycleRunning = true, MRQ_Algorithm **successAlg = NULL );
        
        
        int setIntegerParameter(const char *name, const long int value);
        
        int setDoubleParameter(const char *name, const double value);
        
        
        
        int setIntegerParameterToAlgs(const char *name, const long int value);
        
        int setDoubleParameterToAlgs(const char *name, const double value);
        
        int setStringParameterToAlgs(const char *name, const char * value);
        
        int setParametersToAlgs(const MRQ_GeneralSolverParams &params);
        
        
        friend class MRQ_BBLCallbacks;
        friend class MRQ_IGMA0;
        friend class MRQ_IGMA1;
        friend class MRQ_RensExecutor;
    };



    class MRQ_AlgorithmParameterSetter : public optsolvers::OPT_ParameterSetter
    {
    public:
        MRQ_Algorithm *algorithm;
            
        MRQ_AlgorithmParameterSetter(MRQ_Algorithm *algorithm)
        {
            this->algorithm = algorithm;
        }
        
        virtual ~MRQ_AlgorithmParameterSetter(){}
        
        virtual int setDoubleParameter(const char *name, const double value)
        {
            return algorithm->setDoubleParameter(name, value);
        }
        
        virtual int setIntegerParameter(const char *name, const long int value)
        {
            return algorithm->setIntegerParameter(name, value);
        }
        
        virtual int setStringParameter(const char *name, const char *value)
        {
            return algorithm->setStringParameter(name, value);
        }
    };



    MRQ_Algorithm* MRQ_newAlgorithm(int algCode, int numberOfNLEqualityConstraints);


    int MRQ_testRun(MRQ_MINLPProb &prob, const char *probName, const bool printAlgParameters);
}










#endif /* MRQ_ALGCLASSES_HPP_ */






