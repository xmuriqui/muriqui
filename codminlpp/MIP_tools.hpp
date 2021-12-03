


#ifndef MIP_TOOLS_HPP
#define MIP_TOOLS_HPP

#include <cmath>
#include <cstring>
#include <iostream>


#if MIP_CPP_MULTITHREADING
    #include <mutex>
#endif


#include "MIP_constants.hpp"

#define MIP_PREPRINT "minlpproblem: "


#ifdef __FILE__
    #ifdef __LINE__
        #define MIP_DEF_GETFILELINE 1
    #endif
#endif


#ifdef MIP_DEF_GETFILELINE

    #define MIP_GETFILELINE  \
        " on file: '" __FILE__  "' function: '" << __func__ <<  "' line: " << __LINE__
#else
    #define MIP_GETFILELINE ""
#endif




#define MIP_PRINTMEMERROR std::cerr << MIP_PREPRINT << "Memory error" << MIP_GETFILELINE << "\n"


#define MIP_PRINTERROR std::cerr << MIP_PREPRINT << "Error" << MIP_GETFILELINE << "\n"


#define MIP_PRINTERRORNUMBER(number) std::cerr << MIP_PREPRINT << "Error " << number << MIP_GETFILELINE << "\n"

#define MIP_PRINTCALLBACKERRORNUMBER(number) std::cerr << MIP_PREPRINT << "Callback function error " << number << MIP_GETFILELINE << "\n"


#define MIP_PRINTERRORMSG(msg) std::cerr << MIP_PREPRINT << msg << MIP_GETFILELINE << "\n"

#define MIP_PRINTERRORMSGP(msg, arg) std::cerr << MIP_PREPRINT << msg << arg << MIP_GETFILELINE << "\n"

#define MIP_PRINTMSG(msg) std::cout << MIP_PREPRINT << msg;



#define MIP_getchar()  MIP_PRINTERRORMSG("stopped in a getchar"), getchar()



#define MIP_IFERRORRETURN(error, retCode) \
    if(error){ MIP_PRINTERRORNUMBER(error); return retCode;}

#define MIP_IFMEMERRORRETURN(errorExp) \
    if(errorExp){ MIP_PRINTMEMERROR; return MIP_MEMORY_ERROR;}

#define MIP_IFERRORGOTOLABEL(error, varToGetRetCode, retCode, destinLabel) \
    if(error){ MIP_PRINTERRORNUMBER(error); varToGetRetCode = retCode; goto destinLabel;  }

#define MIP_IFMEMERRORGOTOLABEL(errorExp, varToGetRetCode, destinLabel) \
    if(errorExp){ MIP_PRINTMEMERROR; varToGetRetCode = MIP_MEMORY_ERROR; goto destinLabel;  }

#define MIP_IFERRORSETVAR(error, varToGetRetCode, retCode) \
    if(error){ MIP_PRINTERRORNUMBER(error); varToGetRetCode = retCode; }

#define MIP_IFCALLBACKERRORSETVAR(error, varToGetRetCode) \
    if(error){ MIP_PRINTCALLBACKERRORNUMBER(error); varToGetRetCode = MIP_CALLBACK_FUNCTION_ERROR; }

#define MIP_IFCALLBACKERRORRETURN(error) \
    if(error){ MIP_PRINTCALLBACKERRORNUMBER(error); return MIP_CALLBACK_FUNCTION_ERROR; }


#define MIP_IFCALLBACKERRORGOTOLABEL(error, varToGetRetCode, destinLabel) \
    if(error){ MIP_PRINTCALLBACKERRORNUMBER(error); varToGetRetCode = MIP_CALLBACK_FUNCTION_ERROR; goto destinLabel;  }

//to GAMS
#define MIP_IFERRORGOTOLABELANDPRINTBUFFER(errorExp, varToGetRetCode, retCode, buffer, destinLabel) \
    if(errorExp){ MIP_PRINTERRORMSG(buffer); varToGetRetCode = retCode; goto destinLabel;  }




namespace minlpproblem
{
    
    //class to write files using a max number of bytes in a line
    class MIP_StreamMBLine
    {
        
    public:
        
        unsigned int maxBytes; //maxbytes in a line
        unsigned int nBytes; //number of bytes in the current line
        std::ostream *out;
        
        
        
        MIP_StreamMBLine( std::ostream *outputstream, unsigned int maxBytes = 255 )
        {
            out = outputstream;
            nBytes = 0;
            this->maxBytes = maxBytes;
        }
        
        
        void endline()
        {
            *out << std::endl;
            nBytes = 0;
        }
        
        
        //if size is negative, function calculate the size
        void print( const char *s, int size = -1 )
        {
            if( size < 0)
                size = strlen(s);
            
            //we deserve two bytes because in some systems, end of line is made by CR + LF
            if( size + nBytes >= maxBytes -1 )
                endline();
            
            
            *out << s;
            
            nBytes += size;
        }
    };
    
    
    
    //mutex to implement mutual exclusion on multithread proceduring. In this way, the rest of code do not need worry about haw library to multithreading is being used...
    class MIP_Mutex
    {
        
    public:
        
        #if MIP_CPP_MULTITHREADING
            std::mutex mymutex;
        #else
            #if MIP_OMP_MULTITHREADING
                omp_lock_t mutex;
            #endif
        #endif
        
        
        MIP_Mutex();
        
        void initialize();
        
        int lock( );
        
        int tryLock( );
        
        int unlock( );
        
        void destroy();
        
        ~MIP_Mutex();
    };
    
    
    
    
    
    
    template <class myClass>
    inline myClass MIP_max(const myClass &a, const myClass &b ) 
    {
        return a >= b ? a : b;
    }

    template <class myClass>
    inline myClass MIP_min(const myClass &a, const myClass &b)
    {
        return a <= b ? a : b;
    }
    
    
    template <class myClass>
    inline void MIP_setAllArray( unsigned int size, myClass *a, myClass value )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for( unsigned int i = 0; i < size; i++ )
            a[i] = value;
    }
    
    
    //function to shift array in marked positions...
    template <class myClass>
    inline unsigned int MIP_shitfValuesInArray( const unsigned int size, myClass *a, bool *remIndexes )
    {
        unsigned int nsize = size;
        
        for( unsigned int i = size-1; i >= 1; i-- )
        {
            if( remIndexes[i] )
            {
                for(unsigned int j =i+1; j < nsize; j++)
                    a[j-1] = a[j];
                
                nsize--;
            }
            
        }
        
        //first index
        if( remIndexes[0] )
        {
            //note we "should" decrease nsize after the for above, but since we did it before, we use j < nsize instead of j < nsize-1
            for(unsigned int j = 1; j < nsize; j++)
                a[j-1] = a[j];
            
            nsize--;
        }
        
        return nsize;
    }
    
    
    template < class myClass>
    void MIP_printLowerTriangleMatrix( unsigned int n, myClass *M )
    {
        for(unsigned int i = 0, k = 0; i < n; i++)
        {
            for(unsigned int j = 0; j <= i; j++, k++)
                std::cout << M[k] << " ";
            
            std::cout << std::endl;
        }
    }
    
    
    
    

    template <class myClass>
    inline void MIP_secFree(myClass* &p)
    {
        if(p)
        {
            free(p);
            p = NULL;
        }
    }


    template <class myClass1, class myClass2>
    inline void MIP_copyArray(const unsigned int size, const myClass1 *origin, myClass2 *destin)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            destin[i] = origin[i];
    }
    
    
    template <class myClass1, class myClass2>
    inline void MIP_accumulateArray(const unsigned int size, const myClass1 *origin, myClass2 *destin)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            destin[i] += origin[i];
    }
    
    
    template <class myClass>
    static inline void MIP_multiplyAllArray(unsigned int size, myClass *a, myClass value )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            a[i] *= value;
    }
    
    
    template <class myClass>
    inline myClass MIP_abs(const myClass &a)
    {
        return a >= 0 ? a : -a;
    }
    
    
    inline bool MIP_isIntegerDouble(const double value)
    {
        return floor(value) == value;
    }
    
    
    template <class myClass>
    void MIP_malloc(myClass* &a, const unsigned int nElements)
    {
        a = (myClass*) malloc( nElements * sizeof(myClass) );
    }
    
    
    template <class myClass>
    void MIP_calloc(myClass* &a, const unsigned int nElements)
    {
        a = (myClass*) calloc(nElements, sizeof(myClass));
    }
    
    template <class myClass>
    static inline int MIP_realloc(myClass* &a, const unsigned int size)
    {
        if(size == 0)
        {
            MIP_secFree(a);
        }
        else
        {
            void *p = realloc(a, size * sizeof(myClass) );
            MIP_IFMEMERRORRETURN(!p);
            
            a = (myClass*) p;
        }
        
        return 0;
    }
    
    
    //allocte a (m by n) matrix
    template <class myClass>
    inline int MIP_allocateMatrix( myClass** &M, const unsigned int m, const unsigned int n )
    {
        int code;
        unsigned int i;
        
        
        MIP_malloc(M, m); //M = (myClass **) malloc( m * sizeof(myClass *) );
        
        if( !M )
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTMEMERROR;
            #endif
            code = MIP_MEMORY_ERROR;
            goto termination;
        }
        
        MIP_malloc(M[0], m*n); //M[0] = (myClass *) malloc( m*n * sizeof(myClass) );
        
        if( !M[0] )
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTMEMERROR;
            #endif
            code = MIP_MEMORY_ERROR;
            goto termination;
        }
        
        for(i = 1; i < m; i++)
            M[i] = &M[0][i*n];
        
        code = 0;
        
    termination:
        
        return code;
    }
    
    
    
    
    
    //allocte a (m by n) zeros matrix
    template <class myClass>
    inline int MIP_allocateZeroMatrix( myClass** &M, const unsigned int m, const unsigned int n )
    {
        int code;
        unsigned int i;
        
        
        MIP_malloc(M, m); //M = (myClass **) malloc( m * sizeof(myClass *) );
        if( !M )
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTMEMERROR;
            #endif
            code = MIP_MEMORY_ERROR;
            goto termination;
        }
        
        MIP_calloc(M[0], m*n); //M[0] = (myClass *) calloc( m*n , sizeof(myClass) );
        if( !M[0] )
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTMEMERROR;
            #endif
            code = MIP_MEMORY_ERROR;
            goto termination;
        }
        
        for(i = 1; i < m; i++)
            M[i] = &M[0][i*n];
        
        code = 0;
        
    termination:
        
        return code;
    }


    template <class myClass>
    inline void MIP_freeMatrix( myClass** &M )
    {
        if( M )
        {
            if( M[0] )
                free( M[0] );
            free(M);
        }
    }
    
    
    inline double MIP_gap(const double x)
    {
        return MIP_abs(x - round(x));
    }
    
    
    template <class myClass>
    inline void MIP_swap(myClass &v1, myClass &v2)
    {
        myClass temp = v1;
        
        v1 = v2;
        v2 = temp;
    }
    
    //sort by selection sort
    template <class indexClass, class weightClass>
    inline void MIP_sortIndicesByWeight(unsigned int size, indexClass *indices, const weightClass *weights  )
    {
        
        for(decltype(size) i = 0; i < size; i++)
        {
            indexClass lowestInd = i;
            weightClass lowest = weights[ indices[i] ];
            
            for(decltype(size) j = i+1; j < size; j++)
            {
                auto indj = indices[j];
                
                if( weights[indj] < lowest )
                {
                    lowest = weights[indj];
                    lowestInd = j;
                }
            }
            
            if( lowestInd != i )
                MIP_swap(indices[i], indices[lowestInd]);
        }
        
    }
    
    
    
    /*inline bool MIP_isValidVarType(const int value)
    {
        return value == MIP_VT_CONTINUOUS || value == MIP_VT_INTEGER;
    } */
    
    
};


#endif

