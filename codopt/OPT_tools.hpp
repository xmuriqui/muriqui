


#ifndef OPT_TOOLS_HPP
#define OPT_TOOLS_HPP


#include <cstdlib>
#include <ctime>
#include <typeinfo>

#include <set>
#include <map>

#include "OPT_solvers.hpp"


#if OPT_CPP_MULTITHREADING
    #include <mutex>
#endif


#if OPT_SUPER_THREAD_DEBUG_MODE //to debug on cases havultiple threads, we save files with some debug messages
    #include <iostream>
    #include <fstream>
    #include <sstream> 
    #include <map>
    #include <thread>
#endif


#ifdef __FILE__
    #ifdef __LINE__
        #define OPT_DEF_GETFILELINE 1
    #endif
#endif


#ifdef OPT_DEF_GETFILELINE

    #define OPT_GETFILELINE  \
        " on file: '" __FILE__  "' function: '" << __func__ <<  "' line: " << __LINE__
#else
    #define OPT_GETFILELINE ""
#endif


#define OPT_STR(att)	#att
//flag below is to expand other defines before stringficate
#define OPT_EXPSTR(att)		OPT_STR(att)


#define OPT_PREPRINT "optsolvers: "


#define OPT_PRINTNONSUPMSG std::cerr << OPT_PREPRINT "method: " << __func__ << " is not supported by " << typeid(*this).name() << "\n"


#define OPT_PRINTMEMERROR std::cerr << OPT_PREPRINT << "Memory error" << OPT_GETFILELINE << "\n"


#define OPT_PRINTERROR std::cerr << OPT_PREPRINT << "Error" << OPT_GETFILELINE << "\n"


#define OPT_PRINTERRORNUMBER(number) std::cerr << OPT_PREPRINT << "Error " << number << OPT_GETFILELINE << "\n"


#define OPT_PRINTERRORMSG(msg) std::cerr << OPT_PREPRINT << msg << OPT_GETFILELINE << "\n"

#define OPT_PRINTERRORMSGP(msg, arg) std::cerr << OPT_PREPRINT << msg << arg << OPT_GETFILELINE << "\n"

#define OPT_PRINTCALLBACKERRORNUMBER(number) std::cerr << OPT_PREPRINT << "Callback function error " << number << OPT_GETFILELINE << "\n"


#define OPT_getchar()  OPT_PRINTERRORMSG("stopped in a getchar"), getchar() 


#define OPT_LIBNOTAVAILABLERET(solverCode)  std::cerr << OPT_PREPRINT << "Solver '" << OPT_getSolverName(solverCode) << "' is not avaiable in this compilation. " << OPT_GETFILELINE << "\n" ; return OPT_LIBRARY_NOT_AVAILABLE


#define OPT_OPERATIONNOTIMPLEMENTEDRET(solverCode)  std::cerr << OPT_PREPRINT << "Operation " << __func__ << " is not implemented yet for '" << OPT_getSolverName(solverCode)  << OPT_GETFILELINE << "\n" ; return OPT_OPERATION_NOT_IMPLEMENTED


#define OPT_OPERATIONNOTSUPORTEDRET(solverCode)  std::cerr << OPT_PREPRINT << "Operation " << __func__ << " is not supported yet for '" << OPT_getSolverName(solverCode)  << OPT_GETFILELINE << "\n" ; return OPT_OPERATION_NOT_SUPPORTED




#define OPT_IFERRORRETURN(error, retCode) \
    if(error){ OPT_PRINTERRORNUMBER(error); return retCode;}

#define OPT_IFMEMERRORRETURN(errorExp) \
    if(errorExp){ OPT_PRINTMEMERROR; return OPT_MEMORY_ERROR;}

#define OPT_IFERRORGOTOLABEL(error, varToGetRetCode, retCode, destinLabel) \
    if(error){ OPT_PRINTERRORNUMBER(error); varToGetRetCode = retCode; goto destinLabel;  }

#define OPT_IFMEMERRORGOTOLABEL(errorExp, varToGetRetCode, destinLabel) \
    if(errorExp){ OPT_PRINTMEMERROR; varToGetRetCode = OPT_MEMORY_ERROR; goto destinLabel;  }

#define OPT_IFERRORSETVAR(error, varToGetRetCode, retCode) \
    if(error){ OPT_PRINTERRORNUMBER(error); varToGetRetCode = retCode; }

#define OPT_IFCALLBACKERRORSETVAR(error, varToGetRetCode) \
    if(error){ OPT_PRINTCALLBACKERRORNUMBER(error); varToGetRetCode = OPT_CALLBACK_FUNCTION_ERROR; }

#define OPT_IFCALLBACKERRORGOTOLABEL(error, varToGetRetCode, destinLabel) \
    if(error){ OPT_PRINTCALLBACKERRORNUMBER(error); varToGetRetCode = OPT_CALLBACK_FUNCTION_ERROR; goto destinLabel;  }





namespace optsolvers{
    
    
    //mutex to implement mutual exclusion on multithread proceduring. In this way, the rest of code do not need worry about haw library to multithreading is being used...
    class OPT_Mutex
    {
        
    public:
        
        #if OPT_CPP_MULTITHREADING
            std::mutex mymutex;
        #else
            #if OPT_OMP_MULTITHREADING
                omp_lock_t mutex;
            #endif
        #endif
        
        
        OPT_Mutex();
        
        void initialize();
        
        int lock( );
        
        int tryLock( );
        
        int unlock( );
        
        void destroy();
        
        ~OPT_Mutex();
    };
    
    
    
    static inline double OPT_getTime(void)
    {
        #if OPT_HAVE_CLOCK_GETTIME
            timespec Time;
            
            clock_gettime(CLOCK_REALTIME, &Time);
            
            return Time.tv_sec + Time.tv_nsec/1.0e9;
        #else
            return (double) time(NULL);
        #endif
    }
    
    
    inline double OPT_calcCPUTtime(const clock_t& clockStart, const clock_t& clockEnd)
    {
        return (double( clockEnd - clockStart ))/CLOCKS_PER_SEC;
    }
    
    
    inline double OPT_calcCPUTtime(const clock_t& clockStart)
    {
        return OPT_calcCPUTtime(clockStart, clock());
    }
    
    
    template <class myClass>
    void OPT_malloc(myClass* &a, const unsigned int nElements)
    {
        a = (myClass*) malloc(nElements* sizeof(myClass));
    }
    
    template <class myClass>
    void OPT_calloc(myClass* &a, const unsigned int nElements)
    {
        a = (myClass*) calloc(nElements, sizeof(myClass));
    }
    
    template <class myClass>
    void OPT_secFree( myClass* &p )
    {
        if( p )
        {
            free(p);
            p = NULL;
        }
    }
    
    
    template <class myClass>
    static inline int OPT_realloc(myClass* &a, const unsigned int nElements)
    {
        if(nElements == 0)
        {
            OPT_secFree(a);
        }
        else
        {
            void *p = realloc(a, nElements * sizeof(myClass) );
            OPT_IFMEMERRORRETURN(!p);
            
            a = (myClass*) p;
        }
        
        return 0;
    }
    
    
    template <class myClass>
    void OPT_secDelete( myClass* &p )
    {
        if( p )
        {
            delete p;
            p = NULL;
        }
    }
    
    
    template <class myClass>
    void OPT_secDeleteArray( myClass* &p )
    {
        if( p )
        {
            delete[] p;
            p = NULL;
        }
    }
    
    
    template <class myClass1, class myClass2>
    void OPT_swap( myClass1 &obj1, myClass2 &obj2 )
    {
        myClass2 temp = obj2;
        
        obj2 = obj1; 
        obj1 = temp;
    }
    
    
    template <class myClass>
    myClass OPT_min(myClass obj1, myClass obj2 )
    {
        return obj1 <= obj2 ? obj1 : obj2;
    }
    
    
    template <class myClass>
    myClass OPT_max( myClass obj1, myClass obj2 )
    {
        return obj1 >= obj2 ? obj1 : obj2;
    }
    
    
    template <class myClass>
    myClass OPT_abs( myClass obj )
    {
        return obj >= 0 ? obj : -obj;
    }
    
    
    template <class myClass>
    inline myClass* OPT_getPointerToConstArray(const myClass *array)
    {
        return  *( (myClass **) ((void **) &array) );
    }
    
    
    template <class myClass>
    inline void OPT_setAllArray(const unsigned int size, myClass *array, const myClass value)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            array[i] = value;
    }
    
    
    template <class myClass1, class myClass2>
    inline void OPT_copyArray(unsigned int size, const myClass1* source, myClass2 *destin)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            destin[i] = source[i];
    }
    
    
    template <class myClass0, class myClass1, class myClass2>
    inline void OPT_copyArrayTimesFactor(unsigned int size, const myClass0 factor, const myClass1* source, myClass2 *destin)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            destin[i] = factor*source[i];
    }
    
    
    
    template <class myClass1>
    inline void OPT_copyShiftArray( unsigned int size, const myClass1 shift, const myClass1 *source, myClass1 *destin )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for( unsigned int i = 0; i < size; i++ )
            destin[i] = source[i] + shift;
    }
    
    
    template <class myClass1, class myClass2>
    inline void OPT_multiplyAllArray( unsigned int size, const myClass1 factor, myClass2 *array )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for( unsigned int i = 0; i < size; i++ )
            array[i] *= factor;
    }
    
    
    //set array sequentially in the interval [begin end].  Here, we supose difference between begin and end can be calculated by an unsigned int.
    template <class myClass>
    static inline void OPT_setSequentialValuesArray(myClass *a, myClass begin, myClass end)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(myClass i = begin; i <= end; i += 1) //we put i += 1 to let more generic
            a[i - begin] = i;
    }
    
    
    //inner product between two arrays x (1Xn) and y (1Xn), i.e., x*y' 
    template <class myClass>
    static inline myClass OPT_innerProduct(const unsigned int n, const myClass *x, const myClass *y)
    {
        myClass r = 0;
        
        for(unsigned int i = 0; i < n; i++)
            r += x[i]*y[i];
        
        return r;
    }
    
    
    template <class myClass>
    inline myClass* OPT_getPointerFromConstPointer(const myClass* p)
    {
        return *((myClass**) ((void**) &p));
    }
    
    
    template <class valueClass, class containerClass>
    inline int OPT_addValuesToContainer( unsigned int nvalues, const valueClass *values, containerClass &container)
    {
        try
        {
            for(unsigned int i = 0; i < nvalues; i++)
                container.insert(values[i]);
        }
        catch(std::bad_alloc& ba)
        {
            return OPT_MEMORY_ERROR;
        }
        
        return 0;
    }
    
    
    template <class containerClass>
    inline int OPT_addNzRowToContainer(const OPT_SparseMatrix &M, containerClass &container)
    {
        //for(auto it = M.beginRowIndex(); it != M.endRowIndex(); ++it)
            //container.insert(*it);
        
        return OPT_addValuesToContainer( M.nNzRowIndices, M.nzRowIndex, container);
    }
    
    
    
    
    inline void OPT_setInicialSolution(const unsigned int n, const unsigned int m, const double* xSource, const double* dualConstrsSource, const double* dualVarsSource, double *xTarget, double *dualConstrsTarget, double *dualVarsTarget)
    {
        if( xSource )
        {
            OPT_copyArray(n, xSource, xTarget);
        }
        else
        {
            if( n >= 1 )
                xTarget[0] = NAN;
        }
        
        
        if( dualConstrsSource )
        {
            OPT_copyArray(m, dualConstrsSource, dualConstrsTarget);
        }
        else
        {
            if( m >= 1 )
                dualConstrsTarget[0] = NAN;
        }
        
        
        if( dualVarsSource )
        {
            OPT_copyArray(n, dualVarsSource, dualVarsTarget);
        }
        else
        {
            if( n >= 1 )
                dualVarsTarget[0] = NAN;
        }
    }



    
    
    
    
    /* This function receives two sets of indices and values and unifies them in a single set. Values in the second set can  overwite values in the first set for indices appearing in both sets. There is no problem if output array is one of input array. Function will work even in this case!
    */ 
    inline int OPT_unifyArraysOfIndicesAndValues(const unsigned int inputSize1, const int *inputIndex1, const double *inputValues1, const unsigned int inputSize2, const int *inputIndex2, const double * inputValues2, unsigned int &outputSize, int *outputIndex, double *outputValues )
    {
        std::map<unsigned int, double> coefs; //maps are ordered in c++
        
        outputSize = 0;
        
        try
        {
            for(unsigned int i = 0; i < inputSize1; i++)
                coefs[ inputIndex1[i] ] = inputValues1[i];
            
            for(unsigned int i = 0; i < inputSize2; i++)
                coefs[ inputIndex2[i] ] = inputValues2[i]; //auxValues2[ cols[i] ] = values[i];
        }
        catch(std::bad_alloc& ba)
        {
            return OPT_MEMORY_ERROR;
        }
        
        for( auto &it : coefs )
        {
            outputIndex[ outputSize ] = it.first;
            outputValues[ outputSize ] = it.second;
            outputSize++;
        }
        
        return 0;
    }
    
    
    
    /*inline void OPT_accumulateRowInArray( minlpproblem::MIP_SparseMatrix &M, const int rowIndex, const double factor, double *a )
    {
        minlpproblem::MIP_SparseRow &row = M[rowIndex];
        const unsigned int nel = row.getNumberOfElements();
        
        for(unsigned int j = 0; j < nel; j++)
            a[ row[j].getColumn() ] += factor * row[j].getValue();
    }*/
    
    
    //calculate nqcons and indqcons 
    inline void OPT_calcNqconstQaudIndex(const int n, const int mquad, const int *quadIndex, minlpproblem::MIP_SparseMatrix *QC, int *nqcons, int *indqcons)
    {
        for(int i = 0; i < n; i++)
        {
            nqcons[i] = 0;
            for( int j = 0; j < mquad; j++ )
            {
                if( QC[ quadIndex[j] ].getNumberOfElementsAtRow(i) > 0 )//if( QC[ quadIndex[j] ][i].getNumberOfElements() > 0 )
                {
                    nqcons[i]++;
                    indqcons[i] = quadIndex[j];
                }
            }
        }
    }
    
    //set in flag[i] as true if row i has some element in M. flags is not initialized
    static inline void OPT_setNonzeroRows(const minlpproblem::MIP_SparseMatrix &M, bool *flags)
    {
        const unsigned int nNzRowIndices = M.nNzRowIndices;
        const unsigned int *nzRowIndex = M.nzRowIndex;
        
        for( unsigned int i = 0; i < nNzRowIndices; i++ )
            flags[ nzRowIndex[i] ]= true;
    }
    
    
    
    static inline int OPT_fillNzRowsLagH(const OPT_MINLPProb &prob, const bool *auxCEval, const int mquad, const int *quadIndex, std::unordered_set<int> &nzRowsLagH)
    {
        const OPT_SparseMatrix *QC = prob.QC;
        int r;
        
        r = OPT_addNzRowToContainer(prob.lagH, nzRowsLagH);
        OPT_IFERRORRETURN(r, r);
        
        r = OPT_addNzRowToContainer(prob.Q, nzRowsLagH);
        OPT_IFERRORRETURN(r, r);
        
        for(int i = 0; i < mquad; i++)
        {
            int ind = quadIndex[i];
            
            if(auxCEval[ind])
            {
                r = OPT_addNzRowToContainer( QC[ind], nzRowsLagH);
                OPT_IFERRORRETURN(r, r);
            }
        }
        
        return 0;
    }
    
    
    static inline int OPT_fillColsNzRowJac(const OPT_MINLPProb &prob, const bool *constrEval, int **colsNzRowJac, int *sizeColsNzRowJac )
    {
        const int m = prob.m;
        
        const OPT_SparseMatrix &A = prob.A;
        const OPT_SparseMatrix &J = prob.J;
        const OPT_SparseMatrix *QC = prob.QC;
        
        for(int i = 0; i < m; i++)
        {
            if(!constrEval[i])
            {
                sizeColsNzRowJac[i] = 0;
                continue;
            }
            
            std::set<int> indexSet;
            
            
            OPT_addValuesToContainer( A.getNumberOfElementsAtRow(i), A[i], indexSet );
            
            OPT_addValuesToContainer( J.getNumberOfElementsAtRow(i), J[i], indexSet );
            
            
            if(QC[i].getNumberOfElements() > 0)
            {
                //we add all index in QC[i] (both rows and columns)
                
                OPT_addValuesToContainer( QC[i].getNumberOfElements(), QC[i][0], indexSet );
                
                OPT_addValuesToContainer( QC[i].nNzRowIndices, QC[i].nzRowIndex, indexSet );
            }
            
            sizeColsNzRowJac[i] = indexSet.size();
            
            if( sizeColsNzRowJac[i] > 0 )
            {
                colsNzRowJac[i] = (int*) malloc( sizeColsNzRowJac[i] * sizeof(int) );
                OPT_IFMEMERRORRETURN(!colsNzRowJac[i]);
                
                int *pindex = colsNzRowJac[i];
                int j = 0;
                for(auto index : indexSet)
                {
                    pindex[j] = index;
                    j++;
                }
            }
            
        }
        
        
        return 0;
    }
    
    
    static inline void OPT_fillColsNzRowHess( const OPT_MINLPProb &prob, const bool *constrEval, const int mquad, const int *quadIndex, const std::unordered_set<int> &nzRowsLagH,  std::unordered_set<int> *colsNzRowHess )
    {
        const OPT_SparseMatrix &Q = prob.Q;
        const OPT_SparseMatrix &lagH = prob.lagH;
        const OPT_SparseMatrix *QC = prob.QC;
        
        int k = 0;
        
        //note indices are not ordered, but there is no problem about that
        
        for(auto it: nzRowsLagH)
        {
            const int i = it;
            
            OPT_addValuesToContainer( lagH.getNumberOfElementsAtRow(i),  lagH[i], colsNzRowHess[k]);
            
            OPT_addValuesToContainer( Q.getNumberOfElementsAtRow(i), Q[i], colsNzRowHess[k]);
            
            for(int p = 0; p < mquad; p++)
            {
                const int j = quadIndex[p];
                
                if(constrEval[j] && QC[j].getNumberOfElements() > 0)
                    OPT_addValuesToContainer( QC[j].getNumberOfElementsAtRow(i), QC[j][i], colsNzRowHess[k]);
            }
            
            k++;
        }
    }
    
    
    static inline int OPT_fillColsNzRowHess( const OPT_MINLPProb &prob, const bool *constrEval, const int mquad, const int *quadIndex, const int nNzRowsLagH, const int* nzRowsLagH, int **colsNzRowHess, int *sizeColsNzRowHess )
    {
        const OPT_SparseMatrix &Q = prob.Q;
        const OPT_SparseMatrix &lagH = prob.lagH;
        const OPT_SparseMatrix *QC = prob.QC;
        int r;
        
        
        //note indices are not ordered, but there is no problem about that
        
        //for(auto it: nzRowsLagH)
        for(int k = 0; k < nNzRowsLagH; k++)
        {
            const int i = nzRowsLagH[k];
            std::set<int> indexSet; // we choose keeping the order of indices
            
            
            r = OPT_addValuesToContainer( lagH.getNumberOfElementsAtRow(i), lagH[i], indexSet);
            OPT_IFERRORRETURN(r, r);
            
            r = OPT_addValuesToContainer( Q.getNumberOfElementsAtRow(i), Q[i], indexSet);
            OPT_IFERRORRETURN(r, r);
            
            for(int p = 0; p < mquad; p++)
            {
                const int j = quadIndex[p];
                
                if(constrEval[j] && QC[j].getNumberOfElements() > 0)
                {
                    r = OPT_addValuesToContainer( QC[j].getNumberOfElementsAtRow(i), QC[j][i], indexSet);
                    OPT_IFERRORRETURN(r, r);
                }
            }
            
            sizeColsNzRowHess[k] = indexSet.size();
            if(sizeColsNzRowHess[k] > 0)
            {
                colsNzRowHess[k] = (int *) malloc( sizeColsNzRowHess[k]* sizeof(int) );
                OPT_IFMEMERRORRETURN(!colsNzRowHess[k]);
                
                int *pinds = colsNzRowHess[k];
                
                int j = 0;
                for(auto index: indexSet)
                {
                    pinds[j] = index;
                    j++;
                }
                
                #if OPT_DEBUG_MODE
                    assert(sizeColsNzRowHess[k] == j);
                #endif
            }
            
        }
        
        return 0;
    }
    
    
    
    //that function get a pointer to a matriz if the jacobian is exclusivelly composed by this matriz (J or A). Otherwise, return NULL
    static inline const OPT_SparseMatrix* OPT_getSingleJacobianPointer(const OPT_MINLPProb &prob, const int mquad)
    {
        const OPT_SparseMatrix *singleJ = NULL;
        
        if(mquad == 0)
        {
            if(prob.A.getNumberOfElements() == 0)
                singleJ = &prob.J;
            else if( prob.J.getNumberOfElements() == 0)
                singleJ = &prob.A;
        }
        
        return singleJ;
    }
    
    
    
    static inline const OPT_SparseMatrix* OPT_getSingleLagHPointer(const OPT_MINLPProb &prob, const int mquad, const int *quadIndex)
    {
        const OPT_SparseMatrix *singleH = NULL;
        
        const int hasH = prob.hasNlConstrs || prob.hasNlConstrs ? 1 : 0;
        
        const int hasQ = prob.Q.getNumberOfElements() > 0 ? 1 : 0;
        
        
        if(hasQ + mquad + hasH == 1)
        {
            if(hasH)
                singleH = &prob.lagH;
            
            if(hasQ)
                singleH = &prob.Q;
            
            if(mquad)
                singleH = &prob.QC[ quadIndex[0] ];
        }
        
        
        return singleH;
    }
    
    
    #if OPT_SUPER_THREAD_DEBUG_MODE
        extern std::map< std::thread::id, std::ostream* >  OPT_thsOut; //defined in tools.cpp
        
        void OPT_createFileToThread( const std::thread::id &tid );
        
        void OPT_closeFileToThread( const std::thread::id &tid );
    #endif
    
    
    
    
    
    //If Jacobian is composed by a single matrix, set singleJ to pointer to this matrix. Otherwise, this fucntion will allocate jacArrays with jac indices.
    // This fucntion returns the number of nonzeros in jacobian
    int OPT_getSingleJacPointerOrCalcJacIndices(const OPT_MINLPProb &prob, const int mquad, const bool hasFreeConstrs, const bool *constrsEval, const OPT_SparseMatrix* &singleJ, int &nnz_jac_g, int* &sizeColsNzRowJac, int** &colsNzRowJac);
    
    
    
    //If lagrangian hessian is composed by a single matrix, set singleLagH to pointer to this matrix. Otherwise, this fucntion will allocate lag hessArrays with jac indices.
    // This fucntion returns the number of nonzeros in lagrangian hessian
    int OPT_getSingleLagHPointerOrCalcLagHIndices(const OPT_MINLPProb &prob, const int mquad, const int *quadIndex, const bool *constrsEval, const OPT_SparseMatrix* &singleLagH, int &nnz_h_lag, int &nNzRowsLagH, int* &nzRowsLagH, int* &sizeColsNzLagH, int** &colsNzRowLagH);
    
    
    //this method fills array values with values from compete jacobian and return the size of this array.
    int OPT_getValuesFromCompleteJacobian(const OPT_MINLPProb& prob, const OPT_SparseMatrix& J, const int* sizeColsNzRowJac, int** colsNzRowJac, const bool* auxCEval, const int newm, const double* x, double* auxVars, double* values);
    
    
    
    //this method fills array values with values from compete jacobian and return the size of this array.
    int OPT_getValuesFromCompleteLagHessian(const OPT_MINLPProb& prob, const OPT_SparseMatrix& lagH, const int mquad, const int *quadIndex, const double obj_factor, const double *lambda, const int nNzRowsLagH, const int *nzRowsLagH, const int *sizeColsNzLagH, int **colsNzRowLagH, double *auxVars, double *values);
    
    
    
    int OPT_evalCompleteLagrangianHessian(const int thnumber, const bool newx, const double *x, const OPT_MINLPProb& prob, OPT_SparseMatrix& lagH, const int mquad, const int *quadIndex, const double obj_factor, const double aditional_nl_obj_factor, const double *lambda, const bool* auxCEval, const int newm, const int nNzRowsLagH, const int *nzRowsLagH, const int *sizeColsNzLagH, int **colsNzRowLagH, double *auxConstrsVars, double *auxConstrs, int &sizeValues, double *values);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
};



#endif

