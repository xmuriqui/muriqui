


#ifndef __MRQ_TOOLS_HPP__
#define __MRQ_TOOLS_HPP__


//#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <cmath>

#include <random>

#include "MRQ_constants.hpp"
#include "MRQ_dataStructures.hpp"


#if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    #include <iostream>
    #include <fstream>
    #include <sstream> 
    #include <map>
    #include <thread>
#endif





#define MRQ_DYN_UPDT_BNDS_ON_STRONG_BRANCH_CALCS 0



#define MRQ_STR(att)	#att
//flag below is to expand other defines before stringficate
#define MRQ_EXPSTR(att)		MRQ_STR(att)

#define MRQ_STRPARINTVALUE(att)  MRQ_STR(att) << ": (" << att << ")"

#define MRQ_STRATT(att)   #att, att
#define MRQ_ATTSTR(att)   att, #att


#define MRQ_STRFFATT(att)  #att ": " << att



#ifdef __FILE__
    #ifdef __LINE__
        #define MRQ_DEF_GETFILELINE 1
    #endif
#endif


#ifdef MRQ_DEF_GETFILELINE
        
    #define MRQ_GETFILELINE  \
        " on file: '" __FILE__  "' line: '" << __LINE__ <<  "' function: " << __func__
#else
    #define MRQ_GETFILELINE ""
#endif


# define MRQ_PREPRINT "muriqui: "


#define MRQ_PRINTMEMERROR std::cerr << MRQ_PREPRINT "Memory error" << MRQ_GETFILELINE << "\n"


#define MRQ_PRINTERROR std::cerr << MRQ_PREPRINT "Error" << MRQ_GETFILELINE << "\n"


#define MRQ_PRINTERRORNUMBER(number) std::cerr << MRQ_PREPRINT "Error " << number << MRQ_GETFILELINE << "\n"


#define MRQ_PRINTCALLBACKERRORNUMBER(number) std::cerr << MRQ_PREPRINT "Callback function error " << number << MRQ_GETFILELINE << "\n"

#define MRQ_PRINTCALLBACKERRORNUMBERWITHNAME(number, name) std::cerr << MRQ_PREPRINT "Callback function error " << number << " at function: " << name << MRQ_GETFILELINE << "\n"


#define MRQ_PRINTERRORMSG(msg) std::cerr << MRQ_PREPRINT << msg << MRQ_GETFILELINE << "\n"

#define MRQ_PRINTERRORMSGP(msg, arg) std::cerr << MRQ_PREPRINT << msg << arg << MRQ_GETFILELINE << "\n"

#define MRQ_PRINTMSG(msg) std::cout << MRQ_PREPRINT << msg;

#define MRQ_getchar()  MRQ_PRINTERRORMSG("stopped in a getchar"), getchar()



#define MRQ_IFERRORRETURN(error, retCode) \
    if(error){ MRQ_PRINTERRORNUMBER(error); return retCode;}

#define MRQ_IFMEMERRORRETURN(errorExp) \
    if(errorExp){ MRQ_PRINTMEMERROR; return MRQ_MEMORY_ERROR;}

#define MRQ_IFERRORGOTOLABEL(error, varToGetRetCode, retCode, destinLabel) \
    if(error){ MRQ_PRINTERRORNUMBER(error); varToGetRetCode = retCode; goto destinLabel;  }

#define MRQ_IFMEMERRORGOTOLABEL(errorExp, varToGetRetCode, destinLabel) \
    if(errorExp){ MRQ_PRINTMEMERROR; varToGetRetCode = MRQ_MEMORY_ERROR; goto destinLabel;  }

#define MRQ_IFERRORSETVAR(error, varToGetRetCode, retCode) \
    if(error){ MRQ_PRINTERRORNUMBER(error); varToGetRetCode = retCode; }

#define MRQ_IFCALLBACKERRORSETVAR(error, varToGetRetCode) \
    if(error){ MRQ_PRINTCALLBACKERRORNUMBER(error); varToGetRetCode = MRQ_CALLBACK_FUNCTION_ERROR; }

#define MRQ_IFCALLBACKERRORGOTOLABEL(error, varToGetRetCode, destinLabel) \
    if(error){ MRQ_PRINTCALLBACKERRORNUMBER(error); varToGetRetCode = MRQ_CALLBACK_FUNCTION_ERROR; goto destinLabel;  }



namespace muriqui
{
    
    int MRQ_getNumCores(void);
    
    double MRQ_getTime(void);
    
    void MRQ_helloOnceTime(std::ostream &out = std::cout);
    
    int MRQ_tryReadAlgChoiceFile(const char *algChoiceFileName, int printLevel, MRQ_ALG_CODE &readAlgCode);
    
    
    int MRQ_writeProblemParameters(const char *probName, MRQ_MINLPProb &prob, FILE *outFile );
    
    template <class myClass>
    static inline void MRQ_secFree(myClass* &p)
    {
        if(p)
        {
            free(p);
            p = NULL;
        }
    }
    
    template <class myClass>
    void MRQ_malloc(myClass* &a, size_t nElements)
    {
        a = (myClass*) malloc( nElements * sizeof(myClass) );
    }
    
    
    template <class myClass>
    void MRQ_calloc(myClass* &a, size_t nElements)
    {
        a = (myClass*) calloc(nElements, sizeof(myClass));
    }
    
    template <class myClass>
    static inline int MRQ_realloc(myClass* &a, size_t numberOfElements)
    {
        if(numberOfElements == 0)
        {
            MRQ_secFree(a);
        }
        else
        {
            void *p = realloc(a, numberOfElements * sizeof(myClass) );
            MRQ_IFMEMERRORRETURN(!p);
            
            a = (myClass*) p;
        }
        
        return 0;
    }
    
    
    template <class myClass>
    static inline void MRQ_secDelete(myClass* &p)
    {
        if(p)
        {
            delete p;
            p = NULL;
        }
    }
    
    
    template <class myClass>
    static inline void MRQ_secDeleteArray(myClass* &p)
    {
        if(p)
        {
            delete[] p;
            p = NULL;
        }
    }
    


    template <class myClass1, class myClass2>
    static inline void MRQ_copyArray(unsigned int size, const myClass1 *origin, myClass2 *destin)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(decltype(size) i = 0; i < size; i++)
            destin[i] = origin[i];
    }
    
    
    template <class sizeClass, class myClass1, class myClass2>
    static inline void MRQ_copyArrayAnySize(sizeClass size, const myClass1 *origin, myClass2 *destin)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(decltype(size) i = 0; i < size; i++)
            destin[i] = origin[i];
    }
    
    
    template <class myClass>
    static inline void MRQ_setAllArray( unsigned int size, myClass *a, myClass value )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            a[i] = value;
    }
    
    
    template <class myClass>
    static inline void MRQ_sumAllArray( unsigned int size, myClass *a, myClass value )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            a[i] += value;
    }
    
    template <class myClass>
    static inline void MRQ_multiplyAllArray(unsigned int size, myClass *a, myClass value )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            a[i] *= value;
    }
    
    //set array sequentially in the interval [begin end].  Here, we supose difference between begin and end can be calculated by an unsigned int.
    template <class myClass>
    static inline void MRQ_setSequentialValuesArray_(myClass *a, myClass begin, myClass end)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(myClass i = begin; i <= end; i += 1) //we put i += 1 to let more generic
            a[i - begin] = i;
    }
    
    
    template <class myClass>
    static inline void MRQ_swap(myClass &a, myClass &b)
    {
        myClass c = a;
        a = b;
        b = c;
    }
    
    
    template <class myClass>
    static inline myClass MRQ_max(const myClass a, const myClass b ) 
    {
        return a >= b ? a : b;
    }
    
    template <class myClass>
    static inline myClass MRQ_min(const myClass a, const myClass b)
    {
        return a <= b ? a : b;
    }
    
    template <class myClass>
    static inline myClass MRQ_abs(const myClass a)
    {
        return a >= 0 ? a : -a;
    }
    
    
    //set item if itemName matches with userName
    template <class myClass>
    static inline int MRQ_setEnum( const char *itemName, const myClass itemValue, const char *userName, myClass &item )
    {
        if( strcmp(itemName, userName) == 0 )
        {
            item = itemValue;
            return 0;
        }
        
        return MRQ_NAME_ERROR;
    }
    
    
    //set value with itemValueStr if itemValue matches with userValue
    template <class myClass>
    static inline int MRQ_setStrByEnum( const myClass itemValue, const char *itemValueStr, const myClass userValue, char *value)
    {
        if( itemValue == userValue )
        {
            strcpy( value, itemValueStr );
            return 0;
        }
        
        return MRQ_VALUE_ERROR;
    }
    
    
    
    //set att to userValue if attName matches with userName
    template <class myClass>
    static inline int MRQ_setAtt( const char *attName, myClass &att, const char *userName, myClass  userValue  )
    {
        if( strcmp(attName, userName) == 0)
        {
            att = userValue;
            return 0;
        }
        
        return MRQ_NAME_ERROR;
    }


    template <class myClass>
    static inline int MRQ_setStrAtt( const char *attName, myClass &att, const char *userName, const char *userValue )
    {
        if( strcmp(attName, userName) == 0 )
        {
            const int r = MRQ_strToEnum(userValue, att);
            
            return r == 0 ? 0 : MRQ_VALUE_ERROR; //if we got a error in the value, we return a positive value to diferentiate of the case where attribute name is wrong
        }
        
        return MRQ_NAME_ERROR;
    }
    
    
    template <class myClass>
    static inline int MRQ_reallocateArray(const unsigned int size, myClass* &a)
    {
        if(size > 0)
        {
            myClass* p = (myClass*) realloc(a, size*sizeof(myClass) );
            if(!p)
            {
                #if MRQ_DEBUG_MODE
                    MRQ_PRINTMEMERROR;
                #endif
                
                return MRQ_MEMORY_ERROR;
            }
            
            a = p;
        }
        else
        {
            MRQ_secFree(a);
        }
        
        return 0;
    }
    
    
    inline double MRQ_calcCPUTtime(const clock_t& clockStart, const clock_t& clockEnd)
    {
        return (double( clockEnd - clockStart ))/CLOCKS_PER_SEC;
    }
    
    
    inline double MRQ_calcCPUTtime(const clock_t& clockStart)
    {
        return MRQ_calcCPUTtime(clockStart, clock());
    }
    
    
    inline double MRQ_gap(const double x)
    {
        return MRQ_abs( x - round(x) );
    }
    
    //gap of binary variable...
    inline double MRQ_binGap(const double x)
    {
        return MRQ_min( x, 1-x );
    }
    
    inline bool MRQ_isIntegerDouble(const double value)
    {
        return floor(value) == value;
    }
    
    
    inline bool MRQ_isValidVarType(const int value)
    {
        return value == MRQ_VT_CONTINUOUS || value == MRQ_VT_INTEGER;
    }
    
    
    inline bool MRQ_isIntegerType(const int value)
    {
        return value == MRQ_VT_INTEGER;
    }
    
    inline int MRQ_isBinary(const double l, const double u)
    {
        return l > -1.0 && l <= 1.0  &&  u >= 0.0 && u < 2.0;
    }
    
    
    inline double MRQ_zuWithTol(const double zu, const double abs_tol, const double rel_tol )
    {
        double r = zu; //we put it to avoid problems with multiple threads	
        return r - MRQ_max(abs_tol, MRQ_abs(r)*rel_tol );
    }
    
    
    
    /*inline void MRQ_constrCompleteGrad(MRQ_MINLPProb &prob, MRQ_SparseMatrix &Jac, const int constr, const double *x, double *grad)
    {
        Jac.getFullRow(constr, grad);
        
        if( prob.QC[constr].getNumberOfElements() > 0 )
            prob.QC[constr].quadraticGradientEvaluation(x, grad, true);
        
        prob.A.getFullRowAccumulation(constr, grad);
    } */
    
    
    inline bool MRQ_isBinarieProblemAtRegion(MRQ_MINLPProb &prob, const double *lx, const double *ux)
    {
        const int n = prob.n;
        int i;
        
        
        for(i = 0; i < n; i++)
        {
            if( MRQ_isIntegerType( prob.xtype[i] ) && lx[i] != ux[i] && ( lx[i] <= -1.0 || lx[i] > 1.0 || ux[i] < 0 || ux[i] >= 2.0 )  )
            {
                return false;
            }
        }
        
        return true;
    }
    
    
    inline bool MRQ_isIntegerSol(const int nI, const int *intVars, const double *sol, const double intTol  )
    {
        
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            //if( MRQ_abs( sol[ind] -round(sol[ind]) ) > intTol )
            if( MRQ_gap( sol[ind] ) > intTol )
                return false;
        }
        
        return true;
    }
    
    
    
    
    inline void MRQ_fixIntVarsOnSol( const int n, const int *xtype, const double *sol, optsolvers::OPT_LPSolver &solver )
    {
        for(int i = 0; i < n; i++)
        {
            if( minlpproblem::MIP_isIntegerType( xtype[i] ) )
            {
                const double fix = round( sol[i] );
                            
                const int r = solver.setVariableBounds( i, fix, fix );
                
                #if MRQ_DEBUG_MODE
                    if( r != 0 )
                    {
                        std::cerr << MRQ_PREPRINT << "Error " << MRQ_GETFILELINE << std::endl;
                        
                        assert(false);
                    }
                #endif
            }
        }
        
    }
    
    
    inline void MRQ_unfixIntegerVars( const int n, const int *xtype, const double *lx, const double *ux, optsolvers::OPT_LPSolver &solver )
    {
        for(int i = 0; i < n; i++)
        {
            if( minlpproblem::MIP_isIntegerType( xtype[i] ) )
            {
                int r = solver.setVariableBounds(i, lx[i], ux[i] )  ;
                
                #if MRQ_DEBUG_MODE
                    if( r != 0 )
                    {
                        std::cerr << MRQ_PREPRINT << "Error " << MRQ_GETFILELINE << std::endl;
                        
                        assert(false);
                    }
                #endif
            }
        }
        
    }
    
    
    
    inline int MRQ_fixIntVarsOnSolByList( const int nI, const int *intVars, const double* sol, optsolvers::OPT_LPSolver& solver )
    {
        int retCode = 0;
        
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            const double fix = round(sol[ind]);
            
            const int r = solver.setVariableBounds( ind, fix, fix );
            
            
            if( r != 0 )
            {
                #if MRQ_DEBUG_MODE
                    MRQ_PRINTERRORNUMBER(r);
                    //std::cout << "i: " << i << " ind: " << ind << " sol: " << sol[ind] << " fix: " << fix << "\n";
                    //assert(false);
                    
                    //MRQ_getchar();
                #endif
                
                retCode = MRQ_NLP_SOLVER_ERROR;
            }
            
        }
        
        return retCode;
    }
    
    
    inline int MRQ_unfixIntegerVarsByList( const int nI, const int *intVars, const double *lx, const double *ux, optsolvers::OPT_LPSolver& solver )
    {
        int retCode = 0;
        
        for( int i = 0; i < nI; i++ )
        {
            const int ind = intVars[i];
            
            const int r = solver.setVariableBounds(ind, lx[ind], ux[ind]);
            
            #if MRQ_DEBUG_MODE
                if( r != 0 )
                {
                    MRQ_PRINTERRORNUMBER(r);// std::cerr << MRQ_PREPRINT << "Error " << MRQ_GETFILELINE << std::endl;
                    
                    retCode = MRQ_NLP_SOLVER_ERROR;
                    
                    assert(false);
                }
            #endif
        }
        
        return retCode;
    }
    
    
    #if MRQ_BB_SUPER_THREAD_DEBUG_MODE
        extern std::map< std::thread::id, std::ostream *>  MRQ_thsOut; //defined in tools.cpp
        /*extern MRQ_Mutex SEMAPH_thsOut; //defined in tools.cpp */
        
        void MRQ_createFileToThread( const std::thread::id &tid );
        
        
        void MRQ_closeFileToThread( const std::thread::id &tid );
        
    #endif
    
    
    class MRQ_Random
    {
        #if MRQ_HAVE_CPP_2011
            std::mt19937_64 gen;
        #endif
    
        long int seed;
        
    public:
        
        MRQ_Random();
        
        MRQ_Random(const long int seed);
        
        //~IQD_Random();
        
        long int getSeed();
        
        //if NULL, current time is used like seed...
        long int setSeed(const long int* seed = NULL);
        
        //generates true with probability prob
        bool randBool(const double prob);
        
        //generates a random integer in a interval [begin  end]
        int randInt(const int begin, const int end);
        
        //generates a random unsigned integer in a interval [begin  end]
        unsigned int randUInt(const unsigned int begin, const unsigned int end);
        
        //generates a uniform random real in the interval [0 1)
        double random();
        
        //generates a normal random real in the interval [0 1)
        double randomNormal();
        
        //generates a random real in the interval [begin end)
        double random(const double begin, const double end);
        
        //generates a random real in the interval [begin end)
        double randomNormal(const double begin, const double end);
        
    };
    
    //put array alements in a random order
    template <class arrayClass>
    void MRQ_shuffleArray(const unsigned int size, arrayClass *array, MRQ_Random &random)
    {
        if( size > 0)
        {
            for(unsigned int i = 0; i < size-1; i++)
            {
                const unsigned int targetIndex = random.randUInt(i, size-1); //we change element i with some other element k >= i. Note, the last element cannot be changed with anyone after it, but it can be changed in previous iterations
                MRQ_swap( array[i], array[targetIndex] );
            }
        }
    }
    
    
    
    
    //class to calculate nolinear function gradients (it is usefull to lienarize constraints in a faster way...)
    
    class MRQ_GradientsEvaluation
    {
        
    public:
        unsigned int thnumber;
        MRQ_SparseMatrix J;
        MRQ_MINLPProb *prob;
        
        //MRQ_GradientsEvaluation(const int thnumber, MRQ_MINLPProb *prob);
        
        void constraintCompleteGradient(const int constr, const double *x, double *grad);
        
        void desallocate();
        
        int evaluateJacobian(const bool newx, const bool *constrEval, const double *x);
        
        int initialize(const unsigned int thnumber, muriqui::MRQ_MINLPProb* prob);
    };
    
    
    
    // class to store points and obj linear approximations. It is useful for new strategies for linearization in linear apporximation algorithms...
    
    class MRQ_LAAPointsStoring
    {
        unsigned int nallocated;
        
        int allocateMorePoints( const unsigned int nnewpoints );
        
        int checkByNonObjCut( optsolvers::OPT_LPSolver &master, const double zu, const int n, const double *objLin, const double rhs, int &indMaster, const double *point);
        
        //Point value is the approximation of objective function...
        int _checkByNonObjCut(optsolvers::OPT_LPSolver &master, const double zu, const int n, int &indMaster, const double pointValue);
        
    public:
        
        unsigned int blockSize; //size of block to allocation new points.
        int npoints;
        
        unsigned int nObjLinRem; //total number of linearization of objective removed
        
        //bool *used; //flag to marck if point[i] is used to linearize objective function
        MRQ_SparseMatrix points;
        //MRQ_SparseMatrix objAp; //objAp[i] has the linear approximation of the objective function in points[i]. In each line, we store the coefficients of variables and rhs in the last position
        
        int *indMaster; //indMilp[i] has the index of objCut on point[i] in Milp masters problem
        
        
        
        
        MRQ_LAAPointsStoring(const unsigned int n = 0);
        
        ~MRQ_LAAPointsStoring();
        
        void desallocate();
        
        void initialize(const unsigned int n = 0);
        
        //n is the dimention of point. Note objApp is the linear approximation of objective fucntion in point, and it should have dimension n+1, and the last position is the RHS of the linear approximation
        int addPoint( unsigned int n, const double *point);
        
        
        
        //that method should be called when we have an update in zu value. Here, we check if we can remove linearization on some points using the linearization on other points...
        int updateObjLinearizationsByNonObjCutPoints2( optsolvers::OPT_LPSolver &master, const double zu, double *auxVals );
        
        
        //this method remove from milp the linearizations of objective function cut by the criterial of MRQ_OLS_NON_OBJ_CUT_POINTS. Actually, only the coefficients of the constraints are removed, and the constraint get empty. ;) RHS of linearization is in the last position of objLin
        int updateObjLinearizationsByNonObjCutPointsByNewPoint( optsolvers::OPT_LPSolver& master, const double zu, const double* objLin);
        
        
    };
    
    
    
    //class to dynamic constraint set. I am not sure if it better define it here or in the algorithms. The idea is help solvers defining variable to enable or disable a constraint to complemenet a big-M strategy. For a constraint, it will be enable/disbale if a variable be fixed in a one or zero value. If we detect a constraint should be disabled, we set the nonlinear flag to constraint to false. In this way, we avoid calculate values and derivatives for this respective constraint. Note that strategy can be used togheter big-M. Original nonlinear flag will be kept while respective associated variable is not fixed... (Note you can use that to enable nonlinear constraints that are originally disable). Note also linear and quadratic part of cosntarints will kept on in any situation... Note if you use that in a algorithm that never fix variables, it will not have any affects. That was developed to be used in algorithms like branch and bound. Use your own risk. 
    
    class MRQ_DynConstrSetUnity
    {
    public:
        int constrIndex;
        int varIndex;
        double varValue;
    };
    
    
    
    
    
    /*
    * we perform a line search betwen two solutions: startSol and endSol. One of those sols should be an infesible solution while the other one should be feasible.
    * 
    * Here, we work on a lambda in [0 1]:
    * 
    * sol = (lambda)*startSol + (1 - lambda)*endSol
    * 
    * We calculate lambda which that sol is in the boundary of feasilbe solution
    * 
    * 
    */
    MRQ_RETURN_CODE MRQ_lineSearch2( const int thnumber, muriqui::MRQ_MINLPProb& prob, const double epsToLambda, const double* startSol, const double* endSol, const double absFeasTol, const double relFeasTol, const bool* cEval, double* auxConstr, double* outSol, double *outLambda = NULL );
    
    
    
    
}



#endif




