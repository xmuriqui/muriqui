

#ifndef __MRQ_SOLVERS_HPP__
#define __MRQ_SOLVERS_HPP__

#include <cassert>
#include <vector>

#include "MRQ_constants.hpp"
#include "MRQ_dataStructures.hpp"
#include "MRQ_tools.hpp"

#include "OPT_solvers.hpp"



namespace muriqui
{
    
    typedef optsolvers::OPT_NLPSolver  MRQ_NLPSolver;
    
    typedef optsolvers::OPT_LPSolver  MRQ_LPSolver;
    
    
    class MRQ_Norm1ConstrHandler
    {
    public:
        
        int addAuxVariables(const int n, const int naux, MRQ_LPSolver* master);
        
        //we assume auxilary variables are added after the n first variables from minp problem. Sol can be NULL. In this case, RHS will not be set
        int addNorm1constr(const int n, const int nIntVars, const int* intVars, const double* sol, MRQ_LPSolver* master, const bool roundSol, int& begNormConstrs);
        
        int changeNorm1constr( const int begNormConstrs, const int nIntVars, const int* intVars, const double* sol, MRQ_LPSolver* master, const bool roundSol);
    };
    
    
    
    
    
    class MRQ_MasterMILPProb
    {
        int thnumber;
        int begNormConstrs;
        bool *auxConstrEval;
        int *auxCols;
        double *auxConstr, *auxVals;
        MRQ_GradientsEvaluation *gradEval;
        MRQ_MINLPProb *prob;
        
        int allocateMemory(const int sizeInds, const int sizeConstrs, const bool allocateGradEval, MRQ_MINLPProb *prob);
        
        
    public:
        
        
        optsolvers::OPT_LPSolver *master;
        
        
        MRQ_MasterMILPProb();
        
        ~MRQ_MasterMILPProb();
        
        
        void desallocateMemory();
        
        void initialize();
        
        //set the variables, linear constraints, variables bounds, variable types and linear part of objective function...
        int setProblemBase( const int thnumber, MRQ_MINLPProb& prob, const int solver, const bool setLinearObj, const bool setQuad, const bool setVarTypes, const double* lx, const double* ux, const int nauxvars, MRQ_GeneralSolverParams* params, const int nthreads );
        
        int addBinaryCut(const int nI, const int* intVars, const double* x, double* auxVals);
        
        //auxConstrEval2, auxCols, auxVals and constrValues are optional output parameters. If you do not pass, internal structures will be used.
        int addConstraintLinearisationPointsToMILP( const double epsActiveConstrToLinearisation, unsigned int* nConstrLinearsSaved, const int nPoints, double** points, const bool incQuadsInMaster, const int constrLinStrat, bool* constrEval, bool* auxConstrEval2 = 0, int* auxCols = 0, double* auxVals = 0, double* constrValues = 0);
        
        int addNorm1constr(const int n, const int nIntVars, const int* intVars, const double *sol);
        
        int addObjLinearisationPointsToMILP( const int nPoints, double **points, const double zu, const bool incQuadsInMaster, const int objLinStrat, int *auxCols = NULL, double *auxVals = NULL, MRQ_LAAPointsStoring *laps = NULL );
        
        
        //auxCols and auxVars are output arguments and they will have the representation of the contraint added.
        int addLinearizedObjFunction( bool newx, const double* x, const bool incQuadsInMaster, const int in_obj_linearisation_strategy, const double zu, MRQ_LAAPointsStoring *laps = NULL, int* auxCols = NULL, double* auxVals = NULL, double* objValue = NULL);
        
        //auxCols and auxVars are output arguments and they will have the representation of the contraint added.
        int addLinearizedObjFunction( bool newx, const double* x, const bool incQuadsInMaster, int* auxCols = NULL, double* auxVals = NULL, double* objValue = NULL);
        
        int addLinearizedNLConstraints( bool newx, const double* x, const bool incQuadsInMaster, const bool* constrEval, int* auxCols, double* auxVals, double* constrValues = NULL, bool constrCalculated = false);
        
        int addLinearizedNLConstraintsByStrategy( const double epsActiveConstrToLinearisation, unsigned int* nConstrLinearsSaved, const bool newx, const double* x, const bool incQuadsInMaster, int constrLinStrat, const bool* constrEval, bool* auxConstrEval2 = 0, int* auxCols = 0, double* auxVals = 0, double* constrValues = 0, const bool constrCalculated = false, const double* masterConstrValues = NULL, const double *newlc = NULL, const double *newuc = NULL  );
        
        
        int changeNorm1constr(const int nIntVars, const int* intVars, const double *sol);
        
        bool minlpProbhasNLConstraints();
        
        
        int setVarBounds(const int n, const double* lb, const double* ub);
        
    };
    
    
    
    class MRQ_NlpFeasBuilder
    {
    public:
        
        //set auxilary variables for constraints and set the objective functio also. Note, we assume solver has already the relaxation builded. nAuxVars is an output argument and returns the number of auxiliary varas added
        int setAuxVariables( optsolvers::OPT_LPSolver *solver, int *nAuxVars = NULL );
    };
    
    
    typedef optsolvers::OPT_NLPNonObjEval MRQ_NLPNonObjEval;
    
    
    //nlp feasibility problem
    class MRQ_NLPFeasProb
    {
    public:
        
        MRQ_NLPNonObjEval *eval;
        optsolvers::OPT_LPSolver *solver;
        
        
        MRQ_NLPFeasProb();
        
        ~MRQ_NLPFeasProb();
        
        void desallocate();
        
        void initialize();
        
        int setProblem(const int solverCode, MRQ_MINLPProb &prob, const double *lx, const double *ux, MRQ_GeneralSolverParams *params, const int thnumber, const bool setSpecialParams, const int nthreads, const double maxCpuTime, const double maxTime);
    };
    
    
    
    //nlp interior point problem
    class MRQ_NLPIPProblem
    {
    public:
        
        MRQ_NLPNonObjEval *eval;
        optsolvers::OPT_LPSolver *solver;
        
        
        MRQ_NLPIPProblem();
        
        ~MRQ_NLPIPProblem();
        
        void desallocate();
        
        void initialize();
        
        int setProblem(const int solverCode, MRQ_MINLPProb &prob, const double *lx, const double *ux, MRQ_GeneralSolverParams *params, const int thnumber, const bool setSpecialParams, const bool considerQuadToInterior, const int nthreads, const double maxCpuTime, const double maxTime);
    };
    
    
    
    
    class MRQ_NLPFeasPumpProb
    {
        int begNorm1Constrs;
        int nn1vars;
        int *n1varIndex; //variable indixes to norm 1
        
        int allocaten1varIndex(const int size);
        
        
    public:
        
        bool setLinearObjTermOnBinVars;
        bool useNorm1;
        
        MRQ_NLPNonObjEval *eval; //taking advantage the classe builded 
        optsolvers::OPT_LPSolver *solver;
        
        
        MRQ_NLPFeasPumpProb();
        
        ~MRQ_NLPFeasPumpProb();
        
        void desallocate();
        
        void initialize();
        
        
        
        int setProblemBase(const int solverCode, MRQ_MINLPProb& prob, const double *lx, const double *ux, const double *lc, const double *uc, MRQ_GeneralSolverParams* params, const int thnumber, const bool setSpecialParams, const bool setLinearObjTermOnBinVars, const bool useNorm1, const int nthreads, const double maxCpuTime, const double maxTime);
        
        
        //if setLinearObjTermOnBinVars is true, we set the terms in objective function related to bin vars as linear terms. Otherwise, we use a quadratic terms when setLinearObjTermOnBinVars is false, or for non binary variables
        int setObjective( const int nI, const int* intVars, const double* lx, const double* ux, const double* sol );
        
    };
    
    
    
    //here, we consider already the objective cut. NOTE: WE CONSIDER OBJECT CUT CONSTRAINT HAS THE INDEX m, WHERE M is the numberof original constraints in the problem
    class MRQ_NoObjWithObjCutNLPEval : public optsolvers::OPT_NonLinearEval
    {
        bool mynewx;
        bool nlConstrs;
        int origm;
        minlpproblem::MIP_NonLinearEval *oeval;
        
        double objFactor;
        
        
    public:
        
        //double *auxVals;
        
        MRQ_NoObjWithObjCutNLPEval(const bool nlConstrs, minlpproblem::MIP_NonLinearEval* originalEval, int origm, double objFactor);
        
        virtual ~MRQ_NoObjWithObjCutNLPEval();
        
        void initialize(const bool nlConstrs, minlpproblem::MIP_NonLinearEval* originalEval, int origm, double objFactor);
        
        virtual int eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value);
        
        virtual int eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *constrEval, const double *x, double *values);
        
        virtual int eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *values);
        
        virtual int eval_grad_nl_constrs_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *constrEval, const double *x, minlpproblem::MIP_SparseMatrix& jacobian);
        
        virtual int eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double objFactor, const double *lambda, minlpproblem::MIP_SparseMatrix& hessian);
    };
    
    
    
    
    class MRQ_GapMinProb
    {
        int *auxCols;
        double *auxVals; //nao esqueca 
        
        int allocateAuxStructures(const int n);
        
        int setObjCutConstr(MRQ_MINLPProb& prob, const double zu);
        
    public:
        
        //bool objCut;
        int auxObjVarIndex;
        int intGapConstrIndex;
        int objCutIndex;
        minlpproblem::MIP_NonLinearEval *eval;
        optsolvers::OPT_QPSolver *solver;
        
        MRQ_GapMinProb();
        
        ~MRQ_GapMinProb();
        
        void desallocate();
        
        void initialize();
        
        int setProblem(const int solverCode, MRQ_MINLPProb& prob, const double* lx, const double* ux, MRQ_GeneralSolverParams* params, const int thnumber, const bool setSpecialParams, const bool setGapExpOnConstr, const bool setGapUpperBound, const int nthreads, const double maxCpuTime, const double maxTime, const int naddvars, const int naddconstrs);
        
        
        //set or update objective cut... 
        int updateObjCutConstr(MRQ_MINLPProb& prob, const double zu);
    };
    

    
    class MRQ_Rounding
    {
        
    public:
        
        virtual ~MRQ_Rounding();
        
        bool roundSolution( const int thnumber, MRQ_MINLPProb &prob, const int nI, const int *intVars, const double absFeasTol, const double relFeasTol, const double* x, const double zu, const bool localSearchIfFeas, const bool localSearchIfInfeas, const double* nlx, const double* nux, double *auxConstr, double *auxVars,  MRQ_NLPSolver& nlp, double* sol, double& fsol );
        
        
        
    private:
        
        virtual bool roundSolution( const int thnumber, muriqui::MRQ_MINLPProb& prob, const int nI, const int* intVars, const double absFeasTol, const double relFeasTol, const double* x, const double* nlx, const double* nux, double* auxConstr, double* sol, double& fsol );
        
        bool localSearchOverNLPSOlver(MRQ_MINLPProb &prob, const int nI, const int *intVars, const double* nlx, const double* nux, double *auxVars, MRQ_NLPSolver& nlp, double* sol, double& fsol );
    };
    
    
    class MRQ_RandomRounding : public MRQ_Rounding
    {
    public:
        
        MRQ_Random random;
        
        MRQ_RandomRounding(long int seed);
        
        virtual ~MRQ_RandomRounding();
        
    private:
        
        virtual bool roundSolution( const int thnumber, muriqui::MRQ_MINLPProb& prob, const int nI, const int* intVars, const double absFeasTol, const double relFeasTol, const double* x, const double* nlx, const double* nux, double* auxConstr, double* sol, double& fsol ) override;
    };
    
    
    
    
    class MRQ_DynConstrSetSetter
    {
        
        inline int __setDyncConstrSet(const int ndcs, const MRQ_DynConstrSetUnity *dcs, const double *nlx, const double *nux, const bool nlFlagValue, MRQ_NLPSolver *nlp)
        {
            for(int i = 0; i < ndcs; i++)
            {
                const int cind = dcs[i].constrIndex;
                const int vind = dcs[i].varIndex;
                
                if( nlx[vind] == nux[vind] && nlx[vind] == dcs[i].varValue )
                    nlp->setConstrNLFlag(cind, nlFlagValue ); //we check indices before start BB
            }
            
            return 0;
        }
        
        
    public:
        
        
        inline int checkIndices(const int n, const int m, const int ndcs, const MRQ_DynConstrSetUnity *dcs)
        {
            for(int i = 0; i < ndcs; i++)
            {
                const int cind = dcs[i].constrIndex;
                const int vind = dcs[i].varIndex;
                
                if( cind < 0 || cind >= m || vind < 0 || vind >= n )
                    return MRQ_INDEX_FAULT;
            }
            
            return 0;
        }
        
        
        inline int resetOriginalNLConstrFlags(const int ndcs, const MRQ_DynConstrSetUnity *dcs, const bool *originalNLConstrFlags, MRQ_NLPSolver *nlp)
        {
            for(int i = 0; i < ndcs; i++)
            {
                const int cind = dcs[i].constrIndex;
                
                nlp->setConstrNLFlag(cind, originalNLConstrFlags[cind] ); //we check indices before start BB
            }
            
            return 0;
        }
        
        
        inline int setDyncConstrSet( const int ndcs0, const MRQ_DynConstrSetUnity *dcs0, const int ndcs1, const MRQ_DynConstrSetUnity *dcs1, const double *nlx, const double *nux, const bool *originalNLConstrFlags, MRQ_NLPSolver *nlp )
        {
            
            resetOriginalNLConstrFlags(ndcs0, dcs0, originalNLConstrFlags, nlp);
            resetOriginalNLConstrFlags(ndcs1, dcs1, originalNLConstrFlags, nlp);
            
            __setDyncConstrSet(ndcs0, dcs0, nlx, nux, false, nlp);
            __setDyncConstrSet(ndcs1, dcs1, nlx, nux, true, nlp);
            
            return 0;
        }
    };
    
    
    
    //DEPRECATED: use MRQ_BoundsUpdaterSolver
    class MRQ_LPboundsUpdater
    {
        
        inline int __optOnVariable(const int index, double &value);
        
    public:
        
        MRQ_LPSolver *lp;
        
        
        MRQ_LPboundsUpdater();
        
        ~MRQ_LPboundsUpdater();
        
        int buildProblem(const int solverCode, MRQ_MINLPProb &prob, const double zu);
        
        int calcMinValueOnVariable(const int index, double &value);
        
        int calcMaxValueOnVariable(const int index, double &value);
        
        void deallocate();
        
        int setVariablesBounds( const int n, const double *lx, const double *ux );
        
        int updateObjectiveCut(MRQ_MINLPProb &prob, const double zu);
    };
    
    
    class MRQ_BoundsUpdaterSolver
    {
        
        inline int __optOnVariable(const int index, double &value);
        
    public:
        
        unsigned int objCutConstrIndex;
        MRQ_LPSolver *lp;
        MRQ_NoObjWithObjCutNLPEval *objCutEval;
        
        MRQ_BoundsUpdaterSolver();
        
        ~MRQ_BoundsUpdaterSolver();
        
        int buildProblem(const int solverCode, const MRQ_MINLPProb &prob, const bool onlyLinearPart, const bool setObjCut, const double zu, const unsigned thnumber, const unsigned int nThreads, MRQ_GeneralSolverParams *solverParams, const double *lc, const double *uc, const bool setSpecialNlpSolverParams  );
        
        int calcMinValueOnVariable(const int index, double &value);
        
        int calcMaxValueOnVariable(const int index, double &value);
        
        void deallocate();
        
        int setVariablesBounds( const int n, const double *lx, const double *ux );
        
        int updateObjectiveCut(const MRQ_MINLPProb &prob, const double zu);
        
        //calculates new variable bounds to nonfixed variables. Note, lx and ux are input and output arguments. lx[indVars[i] ] has the lower bound of indVars[i]. If allIntVars is true, we assume all variables in indVars are supposed to be integers
        int calculateNewVarBounds( const unsigned int n, const int nIndVars, const int *indVars, double *lx, double *ux, const bool allIntVars, bool &infeasible );
    };
    
    
    
    
    
    class MRQ_NewBBNode;
    
    //only to binary problems
    class MRQ_BinVarsOptNlpRelaxSolFixer 
    {
        
    public:
        
        //if node is not NULL, store the new bounds inside node
        int fixBinVarsFromNlpRelaxSol(const int nI, const int *intVars, const double zu, const double objRelaxSol, const double *duallx, const double *dualux, int &nFixed, double *lx, double *ux, MRQ_NewBBNode *node = NULL);
        
    };
    
    
    
    
    static inline MRQ_Rounding* MRQ_newRounding(int roundingStrategy, int seed)
    {
        switch(roundingStrategy)
        {
            case MRQ_RS_NEAREST_INTEGER:
                return new (std::nothrow) MRQ_Rounding;
                break;
            case MRQ_RS_PROBABILISTIC:
                return new (std::nothrow) MRQ_RandomRounding(seed);
                break;
            default:
                #if MRQ_DEBUG_MODE
                    MRQ_PRINTERRORMSGP("Invalid code to set rounding object:  ", roundingStrategy );
                #endif
                return NULL;
        }
    }
    
    
    static inline int MRQ_setVarBounds( optsolvers::OPT_LPSolver *solver, const int n, const double *lb, const double *ub )
    {
        int r;
        
        r = solver->setnVariablesBounds(n, lb, ub);
        
        return r == 0 ? 0 : MRQ_MILP_SOLVER_ERROR;
    }
    
    
    static inline int MRQ_setSpecialParameters(optsolvers::OPT_Solver *solver)
    {
        const int scode = solver->getSolverCode();
        int r;
        
        
        if( scode == optsolvers::OPT_MOSEK )
        {
            
            r = solver->setDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_REL_GAP", 1.0e-4 ); //default: 1.0e-6
            
            r += solver->setDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_MU_RED", 1.0e-6 ); //default: 1.0e-12
            
            r += solver->setDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_PFEAS", 1.0e-6 ); //default: 1.0e-8
            
            r += solver->setDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_DFEAS", 1.0e-6 ); //default: 1.0e-8
            
            r += solver->setDoubleParameter( "MSK_DPAR_INTPNT_TOL_INFEAS", 1.0e-6 ); //default: 1.0e-8
            
            
            #if MRQ_DEBUG_MODE
                if( r != 0 )
                {
                    MRQ_PRINTERRORNUMBER(r);
                    return MRQ_NLP_SOLVER_ERROR;
                }
            #endif
        }
        else if( scode == optsolvers::OPT_IPOPT )
        {
            r = 0;
            
            r += solver->setDoubleParameter("tol", 1.0e-5);
            r += solver->setDoubleParameter("acceptable_tol", 1.0e-3);
            r += solver->setDoubleParameter("mu_target", 1.0e-8);
            
            
            #if MRQ_DEBUG_MODE
                if( r != 0 )
                {
                    MRQ_PRINTERRORNUMBER(r);
                    return MRQ_NLP_SOLVER_ERROR;
                }
            #endif
        }
        else if( scode == optsolvers::OPT_NLP_KNITRO )
        {
            //tuner: parameter to tune knitro
            //ftol: parameter to stop knitro if no changes in objective are detected
            r = solver->setDoubleParameter("opttol", 1.0e-4); //defaut: 1e-6
            r += solver->setDoubleParameter("feastol", 1.0e-4); //defaut: 1e-6
            
            #if MRQ_DEBUG_MODE
                if( r != 0 )
                {
                    MRQ_PRINTERRORNUMBER(r);
                    return MRQ_NLP_SOLVER_ERROR;
                }
            #endif
        }
        
        
        return 0;
    }
    
    //nlx, nux, nlc, nuc are new bounds for variables and constraints, respectivelly. They can be NULL...
    inline int MRQ_setNLPRelaxProb( const MRQ_MINLPProb &prob, const double *nlx, const double *nux, const double *nlc, const double *nuc, optsolvers::OPT_Solver *solver, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int thnumber, const bool setSpecialParameters, optsolvers::OPT_GeneralSolverParams *userParams, const int nthreads, const double maxCpuTime, const double maxTime, const int naddvars, const int naddconstrs)
    {
        const int scode = solver->getSolverCode();
        int r;
        
        #if MRQ_DEBUG_MODE
            //by now, we only accpet True as setVarbOunds.
            assert( setVarBounds ); //TODO: remove this assertion after test all algorithms...
        #endif
        
        
        r = optsolvers::OPT_setMINLPProblem( prob, solver, setObj, setConstrs, false, setVarType, naddvars, naddconstrs);
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        if( setVarBounds )
        {
            const auto n = prob.n;
            const double *plx = nlx ? nlx : prob.lx;
            const double *pux = nux ? nux : prob.ux;
            
            optsolvers::OPT_LPSolver *lp = (optsolvers::OPT_LPSolver*) solver;
            
            r = lp->setnVariablesBounds(n, plx, pux);
            MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        }
        
        if( nlc )
        {
            const decltype(prob.m) m = prob.m;
            const double *plc = nlc ? nlc : prob.lc;
            const double *puc = nlc ? nuc : prob.uc;
            
            optsolvers::OPT_LPSolver *lp = (optsolvers::OPT_LPSolver*) solver;
            
            for(decltype(prob.m) i = 0; i < m; i++)
            {
                r = lp->setConstraintBounds(i, plc[i], puc[i]);
                MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
            }
        }
        
        
        solver->setThreadNumber(thnumber);
        
        solver->setNumberOfThreads(nthreads);
        
        
        if( !std::isinf(maxCpuTime) )
        {
            r = solver->setMaxCPUTime(maxCpuTime);
            #if MRQ_DEBUG_MODE
                if( r != 0 )
                    MRQ_PRINTERRORNUMBER(r);
            #endif
        }
        
        if( !std::isinf(maxTime) )
        {
            r = solver->setMaxTime(maxTime);
            #if MRQ_DEBUG_MODE
                if( r != 0 )
                    MRQ_PRINTERRORNUMBER(r);
            #endif
        }
        
        if( scode == optsolvers::OPT_MOSEK  )
        {
            r = solver->setIntegerParameter( "MSK_IPAR_CHECK_CONVEXITY", 0); //turn off convexity checker
            
            r = solver->setDoubleParameter( "MSK_DPAR_CHECK_CONVEXITY_REL_TOL", 1e-4); //increasying tolerance for convex checking. Even we disabling the convex cheking, mosek is cheking anyway...
            
            if( prob.hasNLConstraints() || prob.hasObjNLTerm() )
            {
                #if OPT_HAVE_MOSEK
                    r += solver->setIntegerParameter( "MSK_IPAR_INTPNT_STARTING_POINT", MSK_STARTING_POINT_SATISFY_BOUNDS );
                #endif
            }
            
        }
        else if( scode == optsolvers::OPT_IPOPT )
        {
            r = solver->setStringParameter("expect_infeasible_problem", "yes");
            #if MRQ_DEBUG_MODE
                if( r != 0 )
                    MRQ_PRINTERRORNUMBER(r);
            #endif
            
            
            r = solver->setStringParameter( "honor_original_bounds", "yes" );
            #if MRQ_DEBUG_MODE
                if( r != 0 )
                    MRQ_PRINTERRORNUMBER(r);
            #endif
        }
        
        
        if( setSpecialParameters )
            MRQ_setSpecialParameters(solver);
        
        
        if(userParams)
        {
            int r = solver->setParameters( *userParams );
            
            #if MRQ_DEBUG_MODE
                if( r != 0 )
                    MRQ_PRINTERRORNUMBER(r);
            #endif
        }
        
        return 0;
    }
    
    
    inline int MRQ_setInteriorPointProbClosedToNLPRelaxProb( const MRQ_MINLPProb &prob, const double *lx, const double *ux, const double *lc, const double *uc, optsolvers::OPT_LPSolver *solver, const int thnumber, const bool setSpecialParameters, optsolvers::OPT_GeneralSolverParams *userParams, const int nthreads, const double maxCpuTime, const double maxTime, const bool considerQuadConstr, const double epsToEnforceInteriorSol)
    {
        const int m = prob.m;
        
        int r;
        
        r = MRQ_setNLPRelaxProb(prob, lx, ux, lc, uc, solver, true, true, true, false, thnumber, setSpecialParameters, userParams, nthreads, maxCpuTime, maxTime, 0, 0);
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        
        const bool *nlConstr = prob.nlConstr;
        const double *plc = lc ? lc : prob.lc, *puc = lc ? uc : prob.uc;
        const MRQ_SparseMatrix *QC = prob.QC;
        
        for( int i = 0; i < m; i++ )
        {
            if( nlConstr[i] || ( QC[i].getNumberOfElements() > 0 && considerQuadConstr )  )
            {
                double nlci, nuci;
                
                if(plc[i] > -MIP_INFINITY)
                {
                    const double eps = MRQ_abs(plc[i]) > 1.0 ? MRQ_abs(plc[i])* epsToEnforceInteriorSol : epsToEnforceInteriorSol ;
                    
                    nlci = plc[i] + eps;
                }
                else
                    nlci = -OPT_INFINITY;
                
                if(puc[i] < MIP_INFINITY)
                {
                    const double eps = MRQ_abs(puc[i]) > 1.0 ? MRQ_abs(puc[i])* epsToEnforceInteriorSol : epsToEnforceInteriorSol ;
                    
                    nuci = puc[i] - eps;
                }
                else
                    nuci = OPT_INFINITY;
                
                r = solver->setConstraintBounds(i, nlci, nuci);
                MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
            }
        }
        
        
        return 0;
    }
    
    
    inline double MRQ_calculateDeltaRHSToConstraintLinearisation(const double constrValue, const int n, const double *grad, const double *sol)
    {
        double deltaRHS = 0.0;//we do that count separeted to avoid numerical errors...
        
        for(int j = 0; j < n; j++)
        {
            deltaRHS += grad[j]*sol[j];
        }
        
        deltaRHS += -constrValue;
        return deltaRHS;
    }
    
    
    /*
    * Calculates binary cut to linearisation based algorithms. cols and values describes columns indices and coeficients of this constraints. nz has the number of nonzeros in cols and values. This cosntraints is less than rhs.
    * 
    * We assume all variables are binaries.
    */ 
    void MRQ_calculateBinCut(const int nI, const int* intVars, const double* sol, double* values, double& rhs);
    
    
    /*Here, it is implicit cols and values has prob.n+1 elements*/
    int MRQ_calculateObjLinearizedConstraint(MRQ_MINLPProb &prob, const int thnumber, bool newx, const double* x, const bool incQuadsInMaster, const double* objValue, int* cols, double* values, double &rhs);
    
    
    void MRQ_calculateConstraintsToBeLinearizedByStrategy( MRQ_MINLPProb &prob, const int constrLinStrat, const double epsActiveConstrToLinearisation, const bool* constrEval, const double* constrValues, const double *masterConstrValues, const double *newlc, const double *newuc, unsigned int *nConstrLinearsSaved, bool* outputConstrEval);
    
    
    
    

    
    inline enum MRQ_LP_SOLVER MRQ_getDefaultLPSolverCode()
    {
        const std::vector<MRQ_LP_SOLVER> solversCode = {MRQ_LP_CPLEX, MRQ_LP_GUROBI, MRQ_LP_XPRESS, MRQ_LP_MOSEK, MRQ_LP_KNITRO, MRQ_LP_CBC, MRQ_LP_GLPK, MRQ_LP_IPOPT, MRQ_LP_ALGENCAN, MRQ_LP_WORHP};
        const unsigned int nSolversCode = solversCode.size();
        
        for(unsigned int i = 0; i < nSolversCode; i++)
        {
            if( optsolvers::OPT_isSolverAvailable(solversCode[i]) )
                return solversCode[i];
        }
        
        return MRQ_UNDEFINED_LP;
    }
    
    
    inline MRQ_MILP_SOLVER MRQ_getDefaultMILPSolverCode()
    {
        const std::vector<MRQ_MILP_SOLVER> solversCode = {MRQ_CPLEX, MRQ_GUROBI, MRQ_XPRESS, MRQ_MILP_KNITRO, MRQ_MILP_MOSEK, MRQ_CBC, MRQ_GLPK};
        const unsigned int nSolversCode = solversCode.size();
        
        for(unsigned int i = 0; i < nSolversCode; i++)
        {
            if( optsolvers::OPT_isSolverAvailable(solversCode[i]) )
                return solversCode[i];
        }
        
        return MRQ_UNDEFINED_MILP;
    }
    
    
    inline MRQ_NLP_SOLVER MRQ_getDefaultNLPSolverCode()
    {
        const std::vector<MRQ_NLP_SOLVER> solversCode = {MRQ_NLP_KNITRO, MRQ_IPOPT, MRQ_NLP_MOSEK, MRQ_ALGENCAN, MRQ_WORHP}; //Mosek is not the first solver because it only address convex problems
        const unsigned int nSolversCode = solversCode.size();
        
        
        
        for(unsigned int i = 0; i < nSolversCode; i++)
        {
            if( optsolvers::OPT_isSolverAvailable(solversCode[i]) )
                return solversCode[i];
        }
        
        return MRQ_UNDEFINED_NLP;
    }
    
    
    inline MRQ_NLP_SOLVER MRQ_getDefaultMinGapNLPSolverCode()
    {
        const std::vector<MRQ_NLP_SOLVER> solversCode = {MRQ_NLP_KNITRO, MRQ_IPOPT, MRQ_ALGENCAN, MRQ_WORHP}; //Mosek cannot be used here because does not address nonconvex problems
        const unsigned int nSolversCode = solversCode.size();
        
        
        for(unsigned int i = 0; i < nSolversCode; i++)
        {
            if( optsolvers::OPT_isSolverAvailable(solversCode[i]) )
                return solversCode[i];
        }
        
        return MRQ_UNDEFINED_NLP;
    }
    
    
}




#endif
