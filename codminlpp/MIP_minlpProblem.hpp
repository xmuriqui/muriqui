/*
* That library implements a representation of Mixed Integer Nonlinear Programming Problems:
* 
* Class for store the MINLP problem data
*
* Min  c'x + 0.5x'Qx + f(x) + d
* subject to:
*
* 	l_c <= a_ix + 0.5x'Q_ix + g_i(x) <= u_c
*
*  l_x <= x <= u_x
*
*  x_i is integer for i \in I
* 
* 
* 
* Author: Wendel Melo
* 
* Date 25-Feb-2014
* 
*/


#ifndef MIP_MINLPPROBLEM_HPP
#define MIP_MINLPPROBLEM_HPP



#include "SPM_NewSparseMatrix.hpp"
#include "MIP_constants.hpp"


//using namespace spm;

namespace minlpproblem
{
    
    typedef newspm::SPM_NewSparseMatrix<int, double> MIP_SparseMatrix;
    typedef MIP_SparseMatrix::Iterator MIP_SparseMatrixIterator;
    typedef MIP_SparseMatrix::RowIndexIterator MIP_SparseMatrixRowIndexIterator;
    //typedef spm::SPM_SparseRow<double> MIP_SparseRow;
    //typedef spm::SPM_SparseElement<double> MIP_SparseElement;
    
    
    
    
    
    
    inline bool MIP_isLinearProblemType( const int value )
    {
        return value == MIP_PT_LP ||value == MIP_PT_MILP;
    }
    
    
    /*inline bool MIP_isQuadraticProblemType(const int value)
    {
        return value == MIP_PT_LP || value == MIP_PT_MILP || value == MIP_PT_QP || value == MIP_PT_MIQP || value == MIP_PT_QCP || value == MIP_PT_MIQCP;
    } */
    
    
    inline bool MIP_isQuadOnlyProblemType( const int value)
    {
        return value == MIP_PT_LP || value == MIP_PT_MILP || value == MIP_PT_QP || value == MIP_PT_MIQP;
    }
    
    
    
    inline bool MIP_isQuadConstrProblemType( const int value )
    {
        return value == MIP_PT_QCP || value == MIP_PT_MIQCP;
    }
    
    
    inline bool MIP_isNonlinearProblemType(const int value)
    {
        return value == MIP_PT_NLP || value == MIP_PT_MINLP;
    }
    
    
    inline bool MIP_isValidVarType(const int value)
    {
        return MIP_VARTYPE_BEGIN <= value && value <= MIP_VARTYPE_END;
    }
    
    
    inline bool MIP_isIntegerType(const int value)
    {
        return value == MIP_VT_INTEGER;
    }
    
    
    inline bool MIP_isIntegerTypeToBranch(const int value)
    {
        return MIP_isIntegerType(value) || value == MIP_VT_CONTINGER;
    }
    
    //l and u are bound to variable
    inline bool MIP_isBinaryVar(const int varType, const double l, const double u)
    {
        return MIP_isIntegerType(varType) && l > -1.0 && u < 2.0;
    }
    
    
    class MIP_NonLinearEval
    {
    public:

        //that method is called before other functions evaluations
        virtual int initialize(int nthreads, int n, int m, int nzNLJac, int nzNLLagHess){ return 0; }

        virtual int eval_nl_obj_part(int threadnumber, int n, bool newx, const double *x, double &value) = 0;
        
        //if constrEval is null, this method should evaluate all nonlinear constraints parts
        virtual int eval_nl_constrs_part(int threadnumber, int n, int m, bool newx, const bool *constrEval, const double *x, double *values) = 0;

        virtual int eval_grad_nl_obj_part(int threadnumber, int n, bool newx, const double *x, double *values) = 0;
        
        //if constrEval is null, this method should evaluate all nonlinear constraints parts
        virtual int eval_grad_nl_constrs_part(int threadnumber, int n, int m, int nz, bool newx, const bool *constrEval, const double *x, MIP_SparseMatrix& jacobian) = 0;
        
        //if lambda is null, method must assume all lambdas are zero
        virtual int eval_hessian_nl_lagran_part(int threadnumber, int n, int m, int nz, bool newx, const double *x, double objFactor, const double *lambda, MIP_SparseMatrix& hessian) = 0;

        virtual void finalize(int nthreads, int n, int m, int nzNLJac, int nzNLLagHess){  }

        virtual ~MIP_NonLinearEval(){};
    };
    
    
    
    
    
    
    /*Class for store the MINLP problem data
    *
    * Min  c'x + 0.5x'Qx + f(x) + d
    * subject to:
    *
    * 	l_c_i <= a_ix + 0.5x'Q_ix + g_i(x) <= u_c_i
    *
    *  l_x <= x <= u_x
    *
    *  x_i is integer for i \in I
    *
    */
    class MIP_MINLPProb{
        
    protected:
        
        bool objLinearPart; //flag to having nonzero linear coefficients at objective fucntion
        int nI; //total number of integer variables
        
        
        void updateNlConstrsFlag();
        
        void updateNumberOfIntegerVars();
        
        
    public:
        
        bool hasNlObj; //flag to indicate if objective function has a general nonlinear term
        
        bool hasNlConstrs; 
        bool *nlConstr; //flags to indicate if constraints have general nonlinear terms
        
        int n; //total number of variables and 
        int m; //number of constraints
        
        int *xtype; //type of variables (continuous, binary or general integer)
        
        int *xprior; //priorities to each variable (usefull in branching process )

        double *lx, *ux; //lower and upper bound to variables
        double *c;  //vector of linear coefficient in objective function
        double d;	//constant in the objective function
        
        
        double *lc, *uc; //lower and upper bound to constraints

        double *x;
        //double *xInit;
        double objFactor, objValue;

        MIP_SparseMatrix Q; //quadratic terms in objective function
        
        MIP_SparseMatrix *QC; //quadratic terms in constraints
        
        MIP_SparseMatrix A; //linear terms in constraints
        
        MIP_SparseMatrix J; //jacobian of constraints (only general nonlinear parts)
        
        MIP_SparseMatrix lagH; //hessian of lagrangian
        
        
        MIP_NonLinearEval *nlEval;
        
        
        
        MIP_MINLPProb();
        
        ~MIP_MINLPProb();
        
        
        int addConstraints(const int ncons);
        
        int addVariables(const int nvars);
        
        //check a triple sparse structure. It can be usefull to check a sparse matrix structure before set some matrix in the problem...
        
        
        /*check obj and Constr Derivatives. step is a general step to calculate finite diferences. Morevoer, you can can specify a particular step for each variavble in steps (it can be NULL also)
        */ 
        int checkFisrtDerivatives(bool checkObj, bool checkConstr, const double *x, const double step, const double *steps, const double tolerance, bool &answer) const;
        
        
        //check obj and Constr Derivatives. step is a general step to calculate finite diferences. Morevoer, you can can specify a particular step for each variavble in steps (it can be NULL also). If lambda is NULL, we set lambda
        int checkSecondDerivatives(bool checkObj, bool checkConstr, const double* x, const double objFactor, const double* lambda, const double step, const double* steps, const double tolerance, bool& answer) const;
        
        
        //if constrEval is NULL, we evaluate all constraints
        int constraintsEval(const int threadnumber, const bool newx, const bool *constrEval, const double *x, double *values) const;
        
        
        //eval only linear part of constraints. evaluation of constraint constrIndices[i] will be put on values[i]. So, values should have at least space for nIndices elements
        int constraintLinearPartEvaluation(const int nIndices, const int *constrIndices, const double *x, double *values );
        
        
        int copyProblemFrom(const minlpproblem::MIP_MINLPProb& other );
        
        
        
        //that function just call the Nonlinear evaluation object to calculate the Jacobian. Maybe we do not need it, but in the filter will be easier do alterations...
        int nlJacobianEval(const int threadnumber, const bool newx, const bool *constrEval, const double *x, MIP_SparseMatrix& jacobian) const;
        
        
        
        //if constrEval is NULL, we evaluate all constraints having nonlinear terms...
        int nlObjAndConstraintsEval(const bool evalObj, const bool evalConstrs, const int threadnumber, bool newx, const bool* constrEval, const double* x, double& objValue, double* constrValues) const;
        
        int objEval(const int threadnumber, const bool newx, const double *x, double &value,  double aditionalNlObjFactor = 1.0) const;
        
        
        //complete gradient of objective function
        int objGradEval(const int threadnumber, const bool newx, const double *x, double *values, double aditionalNlObjFactor = 1.0) const;
        
        
        //that function just call the Nonlinear evaluation object to calculate the hessian of the lagrangian. Maybe we do not need it, but in the future will be easier do alterations...
        int nlpHessianEval(const int threadnumber, const bool newx, const double *x, const double objFactor, const double *lambda, MIP_SparseMatrix& hessian) const;
        
        
        
        void deallocateMatrices();
        
        int deleteConstraints(const int ncons, const int *cons);
        
        void deleteJacobianStructure();
        
        void deleteJacobianStructureLine(const int line);
        
        void deleteLagrangianStructure();
        
        void deleteLagrangianStructureLine(const int line);
        
        int getConstraintBounds( const int constrIndex, double &lb, double &ub ) const;
        
        void getConstraintLowerBounds( double *lbs ) const;
        
        void getConstraintUpperBounds( double *ubs ) const;
        
        //function get the number of linear, quadratic and nonlinear constraints
        void getConstraintStatistcs(int* ml, int* mq, int* mnl) const;
        
        int getConstraintLinearCoef( const int row, const int col, double& value, bool *inStructure = NULL ) const;
        
        int getConstraintLinearPart( const int row, int* nzs, int* cols, double* values ) const;
        
        int getConstraintQuadCoef( const int constrIndex, const int row, const int col, double& value, bool *inStructure = NULL) const;
        
        int getConstraintQuadCoefMatrix( const int constrIndex, int *nzs, int *rows, int *cols, double *values ) const;
        
        int getConstraintQuadCoefMatrix( const int constrIndex, int *rowStart, int *cols, double *values ) const;
        
        //just lower triangle of the row...
        int getConstraintQuadCoefMatrixRow( const int constrIndex, const int row, int* nzs, int* cols, double* values ) const;
        
        //get linear terms by triple sparse
        void getConstraintsLinearPart(int* nzs, int* rows, int* cols, double* values) const;
        
        //get linear terms by compressed row
        void getConstraintsLinearPart(int* rowStart, int* cols, double* values) const;
        
        int getConstraintsNonLinearTermFlag(const int constrIndex, bool &flag) const;
        
        //that function return the number of continuous variables and fill inds with the continuous variables indices
        int getContinuousIndices(int *inds) const;
        
        //that function return the number of integer variables...
        int getIntegerIndices(int* inds ) const;
        
        double getIntGap(const double *sol) const;
        
        int getJacobianRowStructure(const int row, int* nzs, int* cols) const;
        
        int getJacobianStructure(int *nzs, int* rows, int* cols) const;
        
        //get by compressed row format 
        int getJacobianStructure(int* rowStart, int* cols) const;
        
        int getLagrangianHessianRowStructure( const int row, int *nzs, int *cols) const;
        
        int getLagrangianHessianStructure(int *nzs, int* rows, int* cols ) const;
        
        int getLagrangianHessianStructure(int *rowStart, int *cols) const;
        
        MIP_NonLinearEval* getNonLinearEvaluationObject(void) const;
        
        int getNumberOfBinaryVars() const;
        
        int getNumberOfConstraints() const;
        
        int getNumberOfConstraintQuadCoefMatrixTerms( const int constrIndex, int& nzs) const;
        
        int getNumberOfIntegerVars() const;
        
        //return the number of columns in jacobian structure
        int getNumberOfJacobianNonZeros() const;
        
        int getNumberOfJacobianNonZerosAtRow(const int row, int &nzs) const;
        
        //return the number of columns in jacobian structure
        int getNumberOfLagrangianHessianNonZeros(void) const;
        
        int getNumberOfLagrangianHessianNonZerosAtRow( const int row, int &nzs) const;
        
        int getNumberOfLinearCoefs(void) const;
        
        int getNumberOfLinearCoefsInConstr(const int constrIndex, int &nzs) const;
        
        int getNumberOfLinearConstraints() const;
        
        int getNumbersOfNLConstraints( int &nNLEqualityConstraints, int &nNLInequalityConstraints, int &nNLFreeConstraints) const;
        
        int getNumberOfNLConstraints() const;
        
        int getNumberOfNLEqualityConstraints(void) const;
        
        int getNumberOfObjQuadTerms() const;
        
        int getNumberOfObjQuadCoefMatrixRowTerms(int row, int &nzs) const;
        
        int getNumberOfQuadCoefsInConstr(const int constrIndex, int &nzs) const;
        
        int getNumberOfQuadConstraints() const;
        
        int getNumberOfQuadEqualityConstraints() const;
        
        int getNumberOfQuadMatricesInConstrs() const;
        
        int getNumberOfVars() const;
        
        double getObjConstant() const;
        
        double getObjFactor() const;
        
        int getObjLinearCoef( const int varIndex, double &value ) const;
        
        void getObjLinearCoefs( double *value ) const;
        
        int getObjQuadCoef( const int row, const int col, double &value, bool *inStructure = NULL ) const;
        
        //just lower triangle of the row...
        int getObjQuadCoefsMatrixRow( const int row, int *nzs, int *cols, double *values ) const;
        
        void getObjQuadCoefsMatrix( int *nzs, int *rows, int *cols, double *values ) const;
        
        //get by compressed row format
        void getObjQuadCoefsMatrix(int *rowStart, int *cols, double *values) const;
        
        MIP_PROBLEMTYPE getProblemType(void) const;
        
        
        //return a list with constraints having  quadratic coefficient matrix, even if nonlinear term is present. This method return number of indices.
        int getQuadMatrixConstraintInds( int *indices ) const;
        
        //that function build the reverse map of function getIntegerIndices return the number of integer variables... Size of inds should be at least n
        int getReverseIntegerIndices(int* inds ) const;
        
        int getVariableBounds( const int varIndex, double &lb, double &ub ) const;
        
        void getVariableLowerBounds(double *lbs) const;
        
        void getVariableUpperBounds(double *ubs) const;
        
        int getVariableType(const int index, MIP_VARTYPE &type) const;
        
        bool hasConstraintNLTerm(int index) const;
        
        bool hasLinCoefObj(void) const;
        
        bool hasNLConstraints(void) const;
        
        bool hasNLTerm(void) const;
        
        bool hasObjNLTerm(void) const;
        
        bool hasQuadMatrixInSomeConstraint(void) const;
        
        void initialize(void);
        
        bool isBinaryProblem(void) const;
        
        //recieve the constraints already evaluated to check if they are feasible
        bool isConstrValuesFeasible( const double absTol, const double relTol, const double *constrValues ) const;
        
        /* 
        * constrValues can be NULL. In this case, the method allocate a interval array to do the evaluations. So, it is faster pass the array constrValues...
        * constrEval can be NULL also. In this case, all constraints should be evaluated
        * */
        int isFeasibleToConstraints(const int thnumber, const double* x, const bool newx, const bool* constrEval, const double absTol, const double relTol, bool& answer, double* constrValues = NULL) const;
        
        bool isIntegerSolution(const double *sol, const double int_tol) const;
        
        int isLinearConstraint(const int index, bool &response) const;
        
        int isQuadraticConstraint(const int index, bool &response) const;
        
        void print(std::ostream& out = std::cout, const bool printDerivativeStructures = false) const;
        
        int readMIQCPModelInFile(const char* fileName);
        
        
        int removeConstraints( const int nconstrs, const int *indexes );
        
        int removeVars( const int nvars, const int *indexes );
        
        //warning: that method just set a coefficient already in the structure. If coefficient is not in he structure, an error code is retruned...
        int setConstraintLinearCoefInStructure( const int constrIndex, const int varIndex, const double value );
        
        
        //warning: that method is really inefficient. It is better use methods setConstraintLinearPart or setConstraintsAllLinearPart
        int setConstraintsLinearColumn( const int col,  const int nzs, const int* rows, const double* values );
        
        
        //if ncols = 0, we will use ncols = n
        int setConstraintLinearPart(const int row, const double* a, const int ncols = 0, const double zeroTol = 0.0);
        
        //overwrite the entire row...
        int setConstraintLinearPart(const int row, const int nz, const int* cols, const double* vals);
        
        //int setConstraintLinearPart(const int row, MIP_SparseRow &a);
        
        //warning: that method overwrites all linear elements alredy defined.
        int setConstraintsLinearPart(MIP_SparseMatrix &A);
        
        //warning: that method overwrites all linear elements alredy defined.
        int setConstraintsLinearPart(const int nzs, const int *rows, const int *cols, const double *values);
        
        /*
        set by compressed row format
        warning: that method overwrites all linear elements alredy defined. */
        int setConstraintsLinearPart(const int *rowStart, const int *cols, const double *values);
        
        int setConstraintNonLinearTermFlag(const int index, const bool answer);
        
        int setConstraintsNonLinearTermFlag(const bool *answer);
        
        int setConstraintLowerBound(const int index, const double lb);
        
        int setConstraintLowerBounds(const double *lb);
        
        int setConstraintQuadCoefsMatrix(const int index, const int nzs, const int* rows, const int* cols, const double* vals);
        
        //set by compressed row format (you must set only the lower triangle)
        int setConstraintQuadCoefsMatrix(const int index, const int* rowStart, const int* cols, const double* vals);
        
        //M is a full matrix ordered by rows... we just look to lower triangle in M
        int setConstraintQuadCoefsMatrix(const int index, const double* M, const double zeroTol = 0.0, const int nrows = 0, const int ncols = 0);
        
        int setConstraintQuadCoefsMatrixRow( const int constrIndex, const int row, const int nzs, const int *cols, const double *values );
        
        int setConstraintQuadCoefsMatrix(const int index, const MIP_SparseMatrix& M);
        
        int setConstraintUpperBound(const int index, const double ub);
        
        int setConstraintUpperBounds(const double *ub);
        
        //int setInitialSolutions(const double *xI);
        
        //int setInitialSolution(const int index, const double value);
        
        int setJacobianStructure(const int nzs, const int* rows, const int* cols);
        
        //set by compressed row format
        int setJacobianStructure(const int *rowStart, const int* cols);
        
        //that method copy the structure from input matrix
        int setJacobianStructure( MIP_SparseMatrix &M);
        
        //warning: that method store the jacobian of a nonlinear constraint. If that line have already been defined, it will be overwritten
        int setJacobianRowStructure( const int row, const int nzs, const int* cols);
        
        //int setJacobianRowStructure( const int index, MIP_SparseRow &row );
        
        int setLagrangianHessianStructure(const int nzs, const int* rows, const int* cols);
        
        //set by compressed row format (you must set just the lower triangle)
        int setLagrangianHessianStructure(const int *rowStart, const int* cols);
        
        //just lower triangle...
        int setLagrangianHessianStructure( MIP_SparseMatrix &M);
        
        //just lower triangle...
        int setLagrangianHessianRowStructure(const int row, const int nzs, const int* cols);
        
        
        void setNonLinearEvaluationObject(MIP_NonLinearEval *nl);
        
        void setObjConstant(const double constant);
        
        void setObjFactor(const double factor);
        
        int setObjLinearCoefficient(const int index, const double c);
        
        int setObjLinearCoefficients(const double *c);
        
        int setObjLinearCoefficients(const int nzs, const int *cols, const double *c);
        
        void setObjNonLinearTermFlag(const bool answer);
        
        int setObjQuadCoefsMatrix(const int nzs, const int *rows, const int *cols, const double *vals);
        
        //set by compressed row format (you must set he lower triangle)
        int setObjQuadCoefsMatrix(const int *rowStart, const int *cols, const double *vals);
        
        int setObjQuadCoefsMatrix(const double *M, const double zeroTol = 0.0, const int nrows = 0, const int ncols = 0);
        
        int setObjQuadCoefsMatrix(const MIP_SparseMatrix &M);
        
        //just lower triangle of the row...
        int setObjQuadCoefsMatrixRow( const int row, const int nzs, const int *cols, const double *values );
        
        int setParametersAndAllocate(const int n, const int m);
        
        //set first n lower bounds
        int setVariableLowerBounds(const int n, const double *lb);
        
        
        int setVariableLowerBounds(const int ninds, const int *inds, const double *lb);
        
        
        int setVariableLowerBound(const int index, const double lb);
        
        //set variable priority. It can be consider in B&B algorithms to choosing an index to perform branch. Higher is priority value, higher is the variable priority over others... default value is zero.
        int setVariablePriority(const int index, const int priority);
        
        //set first n upper bounds
        int setVariableUpperBounds(const int n, const double *ub);
        
        int setVariableUpperBounds(const int ninds, const int *inds, const double *ub);
        
        
        int setVariableType(const int index, const MIP_VARTYPE type);
        
        int setVariableTypes(const MIP_VARTYPE *types);
        
        int setVariableUpperBound(const int index, const double ub);
        
        int writeMIQCPModelInAMPLFile(const char* outFileName, const char* solverOption = 0);
        
        int writeMIQCPModelInGAMSFile(const char* outFileName, const char* probName, const char* solverOption = NULL, const bool optfile = false, const double optca = 1.0e-4, const double optcr = 1.0e-5);
        
        int writeMIQCPModelInFile(const char* fileName);
        
        int writeMIQCPModelInLPFile(const char* fileName);
        
        
    protected:
        
        //int allocateMatrices();
        
    };
    
    
    
    class MIP_IntegratedHessian
    {
        int allocate(const unsigned int n);
        
        void setSPMonHessianValues( const int n, const minlpproblem::MIP_SparseMatrix& M, const double factor, double* values );
        
    public:
        
        MIP_MINLPProb *prob;
        
        MIP_SparseMatrix lagH;
        
        unsigned int nzs;
        
        int *rows;
        int *cols;
        double *values;
        
        
        unsigned int *hessIndex; //matrix for getting indexes for triple sparse in hessian...
        
        
        MIP_IntegratedHessian(MIP_MINLPProb *prob);
        
        ~MIP_IntegratedHessian();
        
        
        
        
        int allocateTripleSparseArrays(const unsigned int size );
        
        int buildStructures( const bool setLagHessianMatrixCopy, const bool setTripleSparArrays );
        
        void desallocate();
        
        void desallocateTripleSparseArrays();
        
        int evalCompleteHessian( const int thnumber, const bool newx, const double *x, const double objFactor, const double *lambda, MIP_SparseMatrix *lagH = NULL, int *rows = NULL, int *cols = NULL, double *values = NULL );
        
        void initialize(MIP_MINLPProb *prob);
        
        void setTripleSparseArrays( int *rows = NULL, int *cols = NULL );
        
    };
    
    
    /* Class to storage the linear constraits where each variable apears 
     */
    class MIP_ConstraintsByColumnsStorager
    {
        
    public:
        
        unsigned int *constrColsOffset;
        unsigned int *constrCols; //here, we save the constraints indexes where each variable appears (only to linear constraints). So, for variable i, its constraint indices start from constrCols[ constrColsOffset[i] ] until constrCols[ constrColsOffset[i+1]-1 ]
        
        
        MIP_ConstraintsByColumnsStorager();
        
        ~MIP_ConstraintsByColumnsStorager();
        
        //save structures to allowe we know ehat are the constraints indices where each variable appears (only to linear constraints)
        int storageConstraintsByColumns(const MIP_MINLPProb &prob, const bool considerQuadConstrs  );
        
        void deallocate();
        
    };
    
    
    
    class MIP_Preprocessing
    {
    public:
        
        double in_abs_feas_tol;
        double in_rel_feas_tol;
        
        double in_abs_bound_tol;
        double in_rel_bound_tol;
        
        char *auxFlags; //to store redundant constraints...
        double *auxValues;
        const MIP_MINLPProb *prob;
        
        
        
        MIP_Preprocessing(const MIP_MINLPProb *prob = NULL );
        
        ~MIP_Preprocessing();
        
        
        int allocateMemory( const int nvars, const int nconstrs );
        
        void deallocateMemory();
        
        void initialize(const MIP_MINLPProb *prob);
        
        void resetParameters();
        
        
        //note, lx and ux are input and output arguments. They will be overwritte with new bounds... if a constraint is detected redundant, we put (-)MIP_INFINITY in it rspective bound...
        int preprocess(const bool preprocQuadConstrs, const bool preprocObj, const double zu, double* lx, double* ux, bool& varBoundsUpdt, bool &constrBoundsUpdt, const double* inlc = NULL, const double* inuc = NULL, double* outlc = NULL, double* outuc = NULL );
        
        
        /*
         * Here, we use constraint by column storager to speedup the preprocessing. In this sittuation, we can receive a list having variables which we will use to evaluate their respective constraints. The idea is using this method in cases where user changes some variable bounds by itself and it desires evaluate the effect in the preprocessing. 
         * 
         * note, lx and ux are input and output arguments. They will be overwritte with new bounds... if a constraint is detected redundant, we put (-)MIP_INFINITY in it rspective bound... 
         * 
         */
        int preprocess(const unsigned int nTargetVariables, const int *targetVariables, const MIP_ConstraintsByColumnsStorager &ccstorager, const bool preprocQuadConstrs, const bool preprocObj, const double zu, double* lx, double* ux, bool& varBoundsUpdt, bool &constrBoundsUpdt, const double* inlc = NULL, const double* inuc = NULL, double* outlc = NULL, double* outuc = NULL );
        
        
        /*just preprocess if obj is linear or quadratic 
         * if varsUpdt is not NULL, we set varsUpdt[i] as true if variable i has its bounds changed. WE DO NOT INITIALIZE varsUpdt!!!!!!!
         */
        int preprocessObjF( const double zu, double *lx, double *ux, bool &updtBounds, bool &intVarBoundsUpdt, bool *varsUpdt = nullptr );
        
        
    private:
        
        
        
        //Calculate the minimum value that xQx can assume. If we have a multiply like x*y and x or y is fixed, "turns" it in a linear term in array coefs. WARNING: coef will be initialized with zeros only if there is some new coefficient, i.e., nNewCoefs > 0
        inline int __calcQminValue(const MIP_SparseMatrix &M, const double factor, const double *lx, const double *ux, double *coefs, double &constValue, unsigned int &nNewCoefs);
        
        /*if leqConstr is true, we assume constraint is <= b. Otherwise, we assume is >= b. That method return true if some bound changes.
         * if varsUpdt is not NULL, we set varsUpdt[i] as true if variable i has its bounds changed. WE DO NOT INITIALIZE varsUpdt!!!!!!!
        */
        bool __preprocLinConstr( const unsigned int nzsRow, const int *rowCols, const double *rowValues, double b, const bool leqConstr, const int* varTypes, double* lx, double* ux, bool& infeasible, bool& redundant, bool &intVarUpdt, bool *varsUpdt = nullptr );
        
        /* That method return true if some bound changes
         * if varsUpdt is not NULL, we set varsUpdt[i] as true if variable i has its bounds changed. WE DO NOT INITIALIZE varsUpdt!!!!!!!
         */
        bool __preprocFullLinConstr( const int n, const double *row, double b, const int* varTypes, double* lx, double* ux, bool& infeasible, bool& redundant, bool &intVarBoundsUpdt, bool* varsUpdt = nullptr);
        
        inline bool __tryUpdateBound( const double rhsMinusSum, const double coef, const int varType,  double &l, double &u);
        
        //if varsUpdt is not NULL, we set varsUpdt[i] as true if variable i has its bounds changed. WE DO NOT INITIALIZE varsUpdt!!!!!!!
        int __preprocQuadConstr( const MIP_SparseMatrix& Q, const unsigned int nzsRow, const int *rowCols, const double *rowValues, double b, const bool leqConstr, const int* varTypes, double* lx, double* ux, bool &updtBounds, bool& redundant, bool &intVarBoundsUpdt, bool *varsUpdt = nullptr ); 
        
        
        /*
         * if varsUpdt is not NULL, we set varsUpdt[i] as true if variable i has its bounds changed. WE DO NOT INITIALIZE varsUpdt!!!!!!!
         */
        inline int __preprocessGeneralConstr(const bool preprocQuadConstrs, const MIP_SparseMatrix &QC, const unsigned int ainzs, const int* aiCols, const double* aiValues, const double &lc, const double &uc, char &flagRedunConstr, const int *xtype, double *lx, double *ux, double *outlci, double *outuci, bool &varBoundsUpdt, bool &constrBoundsUpdt, bool &intVarBoundsUpdt, bool *varsUpdt = nullptr );
        
    };
    
    
    
    
    
    
    class MIP_Preprocessing_old
    {
    public:
        
        double in_abs_feas_tol;
        double in_rel_feas_tol;
        
        double in_abs_bound_tol;
        double in_rel_bound_tol;
        
        char *auxFlags; //to store redundant constraints...
        double *auxValues;
        const MIP_MINLPProb *prob;
        
        
        MIP_Preprocessing_old(const MIP_MINLPProb *prob = NULL );
        
        ~MIP_Preprocessing_old();
        
        
        int allocateMemory( const int nvars, const int nconstrs );
        
        void desallocateMemory();
        
        void initialize(const MIP_MINLPProb *prob);
        
        void resetParameters();
        
        
        //note, lx and ux are input and output arguments. They will be overwritte with new bounds... if a constraint is detected redundant, we put (-)MIP_INFINITY in it rspective bound...
        int preprocess(const bool preprocQuadConstrs, const bool preprocObj, const double zu, double* lx, double* ux, bool& updtVarBounds, bool &updtConstrBounds, const double* inlc = NULL, const double* inuc = NULL, double* outlc = NULL, double* outuc = NULL );
        
        
        //just preprocess if obj is linear or quadratic
        int preprocessObjF( const double zu, double *lx, double *ux, bool &updtBounds );
        
        
        
        
    private:
        
        inline int __calcQminValue(const MIP_SparseMatrix &M, const double factor, const double *lx, const double *ux, double *coefs, double &constValue);
        
        //if leqConstr is true, we assume constraint is <= b. Otherwise, we assume is >= b. That method return true if some bound changes 
        bool __preprocLinConstr( const unsigned int nzsRow, const int *rowCols, const double *rowValues, double b, const bool leqConstr, const int* varTypes, double* lx, double* ux, bool& infeasible, bool& redundant);  //bool __preprocLinConstr( const MIP_SparseRow& row, double b, const bool leqConstr, const int* varTypes, double* lx, double* ux, bool& infeasible, bool& redundant);
        
        //That method return true if some bound changes 
        bool __preprocFullLinConstr( const int n, const double *row, double b, const int* varTypes, double* lx, double* ux, bool& infeasible, bool& redundant);
        
        inline bool __tryUpdateBound( const double rhsMinusSum, const double coef, const int varType,  double &l, double &u);
        
        
        int __preprocQuadConstr( const MIP_SparseMatrix& Q, const unsigned int nzsRow, const int *rowCols, const double *rowValues, double b, const bool leqConstr, const int* varTypes, double* lx, double* ux, bool &updtBounds, bool& redundant ); //int __preprocQuadConstr( const MIP_SparseMatrix& Q, const MIP_SparseRow &a, double b, const bool leqConstr, const int* varTypes, double* lx, double* ux, bool &updtBounds, bool& redundant );
        
    };
    
    
    
    
    
    
    
    
    
    
    //that class is sueful when we define a new MINLP prob from a previous problem, but we do not change  do not alter any general nonlinear ter, we just change number os variables or linear or quadratic constraints. So, we just need this to encapsulate the original callback to perform nonlinear evaluations. In this way, the user callbacks will receive the correct numbers of variables and constraints.
    class MIP_EncapsulatedNonLinearEval : public MIP_NonLinearEval
    {
        int on, om;
        
    public:
        
        //MIP_MINLPProb *oprob;
        MIP_NonLinearEval *originalEval;
        
        
        MIP_EncapsulatedNonLinearEval(MIP_MINLPProb *originalProb, MIP_NonLinearEval *originalEval) : MIP_NonLinearEval()
        {
            //this->oprob = originalProb;
            this->originalEval = originalEval;
            
            on = originalProb->n;
            om = originalProb->m;
        }
        
        
        //that method is called before other functions evaluations
        virtual int initialize(int nthreads, int n, int m, int nzNLJac, int nzNLLagHess) override
        { 
            return originalEval->initialize(nthreads, on, om, nzNLJac, nzNLLagHess ); 
        }
        
        virtual int eval_nl_obj_part(int threadnumber, int n, bool newx, const double *x, double &value) override
        {
            return originalEval->eval_nl_obj_part( threadnumber, on, newx, x, value);
        }
        
        virtual int eval_nl_constrs_part(int threadnumber,  int n, const int m, bool newx, const bool *constrEval, const double *x, double *values) override
        {
            return originalEval->eval_nl_constrs_part(threadnumber, on, om, newx, constrEval, x, values);
        }
        
        virtual int eval_grad_nl_obj_part(int threadnumber, int n, bool newx, const double *x, double *values) override
        {
            return originalEval->eval_grad_nl_obj_part(threadnumber, on, newx, x, values);
        }
        
        virtual int eval_grad_nl_constrs_part(int threadnumber, int n, int m, int nz, bool newx, const bool *constrEval, const double *x, MIP_SparseMatrix& jacobian) override
        {
            return originalEval->eval_grad_nl_constrs_part(threadnumber, on, om, nz, newx, constrEval, x, jacobian);
        }

        //virtual int eval_hessian_nl_obj_part(const int threadnumber, const int n, const int nz, const bool newx, const double *x, MIP_SparseMatrix& hessian) = 0;

        virtual int eval_hessian_nl_lagran_part(int threadnumber, int n,  int m, int nz, bool newx, const double *x, double objFactor, const double *lambda, MIP_SparseMatrix& hessian) override
        {
            return originalEval->eval_hessian_nl_lagran_part(threadnumber, on, om, nz, newx, x, objFactor, lambda, hessian);
        }
        
        virtual void finalize(int nthreads, int n, int m, int nzNLJac, int nzNLLagHess) override
        {
            return originalEval->finalize(nthreads, on, om, nzNLJac, nzNLLagHess);
        }

        virtual ~MIP_EncapsulatedNonLinearEval(){};
    };
    
    
    /*class to store indices of special classes of
     * constraints on binary variables
     */
    class MIP_BinSumConstrsIndsByClass
    {
    private:
        
        unsigned int nI;
        
        void initialize();
        
        //if classes[classNumber] is not allocated, we alloc for maxNumberOfIndices indices.
        int addIndexToClassArray(const unsigned int classNumber, const unsigned int index, const unsigned maxNumberOfIndices );
        
        int addIndexToKnapsackIndsArray( const unsigned int intVarIndex, const unsigned int knapsackConstrIndex, const unsigned maxNumberOfIndices);
        
        
        void adjustClassesArrays();
        
        void adjustKnapsacksArrays();
        
    public:
        
        const static unsigned int nBinSumConstrClasses = 9;
        
        //number of constraints in each one of MIP_NUMBEROFBINSUMCONSTRCLASSES classes
        unsigned int nClasses[nBinSumConstrClasses];
        
        
        /*CLASS 0: indices of constraints in the form: 
        * x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} = b
        * 
        * b > 0, all x are binary. Here, we sort the constraints by b value.
        */
        //unsigned int *class1; //covered
        
        /*CLASS 1: indices of constraints in the form: 
        * x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} = x_{n_1} + b
        * or
        * x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} = b
        * all x are binary.
        */
        //unsigned int *class2; //covered
        
        /*CLASS 2: indices of constraints in the form: 
        * x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} <= x_{n_1} + b
        * or
        * x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} <= b
        * all x are binary. Here, we have only one variable having coefficient -1.
        */
        //unsigned int *class3; //covered
        
        /*CLASS 3: indices of constraints in the form: 
        * x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} <= b
        * 
        * all x are binary. Here, we can sort the constraints by b value.
        * all x are binary.
        */
        //unsigned int *class4; //covered
        
        /*CLASS 4: indices of constraints in the form 
         * x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} >= b
         * 
         * all x are binary, b > 0
         */
        //unsigned int *class5; //covered
        
        
        /*CLASS 5: indices of constraints in the form (flow constraints): 
        * x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} - x_{n_2} - x_{n_3} - ... - x_{n_q} = b  
        * all x are binary. Here, we have more than one variable having coefficient -1.
        */
        //unsigned int *class6; //covered
        
        
        /*CLASS 6: indices os constraints in the form:
         * x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} <=  x_{n_1} + x_{n_2} + x_{n_3} + ... + x_{n_q} + b
         * or
         * x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} - x_{n_2} - x_{n_3} - ... - x_{n_q} <= b
         */ 
        //unsigned int *class7; //
        
        
        /*CLASS 7: indices of constraints in the form (0-1 knapsack constraints): 
        * a_{p_1}*x_{p_1} + a_{p_2}*x_{p_2} + a_{p_3}*x_{p_3} + ... + a_{p_k}*x_{p_k} <= b
        * all x are binary, at least one a_{p_j} is different of 0 and 1.
        */
        //unsigned int *class8; //covered
        
        
        /*CLASS 8: indices of constraints in the form: 
        * a_{p_1}*x_{p_1} + a_{p_2}*x_{p_2} + a_{p_3}*x_{p_3} + ... + a_{p_k}*x_{p_k} >= b
        * all x are binary, at least one a_{p_j} is different of 0 and 1.
        */
        
        /*TODO: a class to reverse the class 6 (class 6 is a generalziation of class 2, but I think we do not need reverse the class 2 also)
         * 
         * x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} - x_{n_2} - x_{n_3} - ... - x_{n_q} >= b
         * 
        */
        
        //unsigned int *class9; //covered
        
        unsigned int* classes[nBinSumConstrClasses]; //classes[i] has a pointer to an array having indices of constraints belong to class i
        
        
        
        
        /*knapsackIndicesByVar[i] has the knapsack constraints indices that includes variable reverseIntVars[i] (get reverseIntVars from MIP_MINLPProb) */
        unsigned int **knapsackInds;
        unsigned int *nKnapsackInds; //nKnapsackInds[i] has number of knapsack indices in knapsackInds[i]
        
        
        MIP_BinSumConstrsIndsByClass();
        
        ~MIP_BinSumConstrsIndsByClass();
        
        /*lx, ux,lc and uc, can be NULL. In this case, we will use values from prob.
         * reverseIntVars are used only if calculateKnapsakIndices is set to true. You can pass NULL pointers. In this case, this arrays will be get from prob just in time. 
         */
        int calculateIndices( const MIP_MINLPProb &prob, const double *lx = NULL, const double *ux = NULL, const double *lc = NULL, const double *uc = NULL, int *reverseIntVars = NULL, bool calculateKnapsakIndices = false, bool  sortClass0And3IndicesByb = false, bool substituteFixVars = false );
        
        void deallocate();
        
        
        
        
    };
    
    
    
    
    
    int MIP_checkRowStructure(const unsigned int nrows, const unsigned int ncols, const bool symmetric, unsigned int row, unsigned int nzs, int *cols, double *values = NULL);
    
    
    
    int MIP_checkTripleSparseStructure(const int nrows, const int ncols, const bool symmetric, const unsigned int nzs, const int *rows, const int *cols, const double *values = NULL);
    
    
    
    //Jacobian must be already calculated
    void MIP_constrCompleteGrad(const MIP_MINLPProb &prob, const MIP_SparseMatrix &Jac, const int constr, const double *x, double *grad, const bool initializeWithZeros);
    
    
    
    //lagrangian hessian must be already calculated. (calculate only lower triangle part). Remember: prob.objfactor is already considered inside this function
    void MIP_completeLagHessianRow(const MIP_MINLPProb &prob, const int mquad, const int *quadIndex, const MIP_SparseMatrix &lagH, const double objFactor, const double *lambda, const int rowIndex, double *rowValues, const bool initializeWithZeros);
    
    
    
};




#endif




