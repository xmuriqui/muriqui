
#ifndef _BBL_TOOLS_HPP
#define _BBL_TOOLS_HPP


#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <climits>

#include <ostream>

#include "BBL_branchAndBound.hpp"



#define BBL_PREPRINT "branch-and-bound: "

#ifdef __FILE__
    #ifdef __LINE__
        #define BBL_DEF_GETFILELINE 1
    #endif
#endif


#ifdef BBL_DEF_GETFILELINE

    #define BBL_GETFILELINE  \
        " on file: " << __FILE__ << " line: " << __LINE__
#else
    #define BBL_GETFILELINE ""
#endif


#define BBL_PRINTMEMERROR std::cerr << BBL_PREPRINT << "Memory error" << BBL_GETFILELINE << "\n"


#define BBL_PRINTERROR std::cerr << BBL_PREPRINT << "Error" << BBL_GETFILELINE << "\n"


#define BBL_PRINTERRORNUMBER(number) std::cerr << BBL_PREPRINT << "Error " << number << BBL_GETFILELINE << "\n"


#define BBL_PRINTERRORMSG(msg) std::cerr << BBL_PREPRINT << msg << BBL_GETFILELINE << "\n"


#define BBL_getchar()  BBL_PRINTERRORMSG("stopped in a getchar"), getchar() 



#define BBL_IFERRORRETURN(error, retCode) \
    if(error){ BBL_PRINTERRORNUMBER(error); return retCode;}

#define BBL_IFMEMERRORRETURN(errorExp) \
    if(errorExp){ BBL_PRINTMEMERROR; return BBL_MEMORY_ERROR;}

#define BBL_IFERRORGOTOLABEL(error, varToGetRetCode, retCode, destinLabel) \
    if(error){ BBL_PRINTERRORNUMBER(error); varToGetRetCode = retCode; goto destinLabel;  }

#define BBL_IFMEMERRORGOTOLABEL(errorExp, varToGetRetCode, destinLabel) \
    if(errorExp){ BBL_PRINTMEMERROR; varToGetRetCode = BBL_MEMORY_ERROR; goto destinLabel;  }

#define BBL_IFERRORSETVAR(error, varToGetRetCode, retCode) \
    if(error){ BBL_PRINTERRORNUMBER(error); varToGetRetCode = retCode; }

#define BBL_IFCALLBACKERRORSETVAR(error, varToGetRetCode) \
    if(error){ BBL_PRINTCALLBACKERRORNUMBER(error); varToGetRetCode = BBL_CALLBACK_FUNCTION_ERROR; }

#define BBL_IFCALLBACKERRORGOTOLABEL(error, varToGetRetCode, destinLabel) \
    if(error){ BBL_PRINTCALLBACKERRORNUMBER(error); varToGetRetCode = BBL_CALLBACK_FUNCTION_ERROR; goto destinLabel;  }






namespace branchAndBound{
    
    
    
    
    
    typedef unsigned int BBL_uint;
    
    
    
    template <class myClass>
    inline myClass BBL_max(const myClass &a, const myClass &b ) 
    {
        return a >= b ? a : b;
    }
    
    
    template <class myClass>
    inline myClass BBL_min(const myClass &a, const myClass &b)
    {
        return a <= b ? a : b;
    }
    
    
    template <class myClass>
    inline myClass BBL_abs( const myClass &a )
    {
        return a >= 0 ? a : -a;
    }
    
    
    
    template <class myClass>
    inline void BBL_secFree(myClass* &p)
    {
        if(p)
        {
            free(p);
            p = NULL;
        }
    }


    template <class myClass>
    inline void BBL_secDelete(myClass* &p)
    {
        if(p)
        {
            delete p;
            p = NULL;
        }
    }
    
    template <class myClass>
    inline void BBL_secDeleteArray(myClass* &p)
    {
        if(p)
        {
            delete[] p;
            p = NULL;
        }
    }
    
    
    template <class myClass1, class myClass2>
    inline void BBL_swap( myClass1 &v1, myClass2 &v2 )
    {
        myClass1 vt = v1;
        v1 = v2;
        v2 = vt;
    }
    
    
    template <class myClass1, class myClass2>
    inline void BBL_copySequence(const unsigned int size, myClass1 *source, myClass2 *destin)
    {
        
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            destin[i] = source[i];
        
        /*for(i = 0; i < size; i++)
        {
            *destin = *source;
            destin++;
            source++;
        }*/
        
    }
    
    
    template <class myClass>
    void BBL_malloc(myClass* &a, const unsigned int nElements)
    {
        a = (myClass*) malloc( nElements * sizeof(myClass) );
    }
    
    
    template <class myClass>
    void BBL_calloc(myClass* &a, const unsigned int nElements)
    {
        a = (myClass*) calloc(nElements, sizeof(myClass));
    }
    
    template <class myClass>
    static inline int BBL_realloc(myClass* &a, const unsigned int size)
    {
        if(size == 0)
        {
            BBL_secFree(a);
        }
        else
        {
            void *p = realloc(a, size * sizeof(myClass) );
            BBL_IFMEMERRORRETURN(!p);
            
            a = (myClass*) p;
        }
        
        return 0;
    }
    
    
    
    
    inline double BBL_getTime(void)
    {
        #if BBL_HAVE_CLOCK_GETTIME
            timespec Time;
            
            clock_gettime(CLOCK_REALTIME, &Time);
            
            return Time.tv_sec + Time.tv_nsec/1.0e9;
        #else
            return (double) time(NULL);
        #endif
    }
    
    
    inline unsigned int BBL_getNumCores()
    {
        unsigned int r;
        
        #if BBL_CPP_MULTITHREADING
            r = std::thread::hardware_concurrency();
        #elif BBL_OMP_MULTITHREADING
            r = omp_get_num_procs();
        #else
            r = 1;
        #endif
        
        return r;
    }
    
    
    inline double BBL_zuWithTol(const double zu, const double abs_tol, const double rel_tol )
    {
        double r = zu; //we put it to avoid problems with multiple threads
        return r - BBL_max(abs_tol, BBL_abs(r)*rel_tol );
    }
    
    
    
    inline BBL_EXT_EXP_STRATEGY BBL_int2extExpStrategy( const int number )
    {
        switch( number )
        {
            case BBL_EES_DEPTH:
                return BBL_EES_DEPTH;
            case BBL_EES_WIDTH:
                return BBL_EES_WIDTH;
            case BBL_EES_BEST_LIMIT:
                return BBL_EES_BEST_LIMIT;
            case BBL_EES_DEPTH_BEST_LIMIT:
            default:
                return BBL_EES_DEPTH_BEST_LIMIT;
            
        }
    }
    
    
    /* This functions returns the greater sequential integer representable in 32 bits float. Since we have 23 bits in mantissa and there is na implicit bit in the representation. So the number is 2^24
    * source: https://stackoverflow.com/questions/3793838/which-is-the-first-integer-that-an-ieee-754-float-is-incapable-of-representing-e
    */ 
    inline int BBL_greatestSequenceIntInaFloat()
    {
        return 1 << 24;
    }
    
    inline int BBL_lowestSequenceIntInaFloat()
    {
        return -(1 << 24);
    }
    
    
    inline int64_t BBL_greatestSequenceIntInaDouble() //we use int64 instead of long int to work on mingw
    {
        const int64_t v = 1;
        
        return v << 53;
    }
    
    inline int64_t BBL_lowestSequenceIntInaDouble()  //we use int64 instead of long int to work on mingw
    {
        const int64_t v = 1;
        return -(v << 53);
    }
    
    
    //we use int64 instead of long int to work on mingw
    inline int64_t BBL_lowestSequence(BBL_UNION_TYPES_NODE_BOUNDS_POINTER strategy)
    {
        if( strategy == BBL_UTNBP_SCHAR || strategy == BBL_UTNBP_SHORT_SCHAR )
            return SCHAR_MIN;
        else if( strategy == BBL_UTNBP_SHORT_INT || strategy == BBL_UTNBP_SHORT_SHORT_INT )
            return SHRT_MIN;
        else if( strategy == BBL_UTNBP_FLOAT || strategy == BBL_UTNBP_SHORT_FLOAT )
            return BBL_lowestSequenceIntInaFloat();
        else if( strategy == BBL_UTNBP_DOUBLE || strategy == BBL_UTNBP_SHORT_DOUBLE )
            return BBL_lowestSequenceIntInaDouble();
        else
            assert(false); //value was not coded here
    }
    
    
    inline long int BBL_lowestSequence(BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY parentNodeBoundsStrategy)
    {
        return  BBL_lowestSequence( BBL_pnbs2usfdnbp(parentNodeBoundsStrategy, 0) );
    }
    
    
    
    inline long int BBL_greatestSequence(BBL_UNION_TYPES_NODE_BOUNDS_POINTER strategy)
    {
        if( strategy == BBL_UTNBP_SCHAR || strategy == BBL_UTNBP_SHORT_SCHAR)
            return SCHAR_MAX;
        else if( strategy == BBL_UTNBP_SHORT_INT || strategy == BBL_UTNBP_SHORT_SHORT_INT)
            return SHRT_MAX;
        else if( strategy == BBL_UTNBP_FLOAT || strategy == BBL_UTNBP_SHORT_FLOAT)
            return BBL_greatestSequenceIntInaFloat();
        else if( strategy == BBL_UTNBP_DOUBLE || strategy == BBL_UTNBP_SHORT_DOUBLE)
            return BBL_greatestSequenceIntInaDouble();
        else
            assert(false); //value was not coded here
    }
    
    
    inline long int BBL_greatestSequence(BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY parentNodeBoundsStrategy)
    {
        return BBL_greatestSequence( BBL_pnbs2usfdnbp(parentNodeBoundsStrategy, 0) );
    }
    
    
    
    //here, we assume there is space for the new element in array
    template <class indType, class numberType>
    inline void BBL_addElementOnSortedNodeBoundsArray(const unsigned int size, BBL_TNodeBoundsSol<indType, numberType> *array, const BBL_TNodeBoundsSol<indType, numberType> &element)
    {
        if( size == 0 )
        {
            array[0] = element;
        }
        else  if( array[size-1].ind < element.ind )
        {
            //we put in the end
            array[size] = element;
        }
        else
        {
            unsigned int j;
            
            for(j = 0; j < size; j++)
            {
                if( array[j].ind > element.ind  )
                {
                    //shifting positions
                                    
                    //do not use copyarray here because it uses vectorization
                    for(unsigned int k = size; k > j; k--)
                        array[k] = array[k-1];
                    
                    array[j] = element;
                    break;
                }
            }
            
            #if BBL_DEBUG_MODE
                assert(j < size); //we have to enter in the break above
            #endif 
        }
        
    }
    
    
    template <class indType, class numberType>
    inline void BBL_addElementOnSortedNodeBoundsArray(const unsigned int size, BBL_TNodeBoundsSol<indType, numberType> *a, unsigned int index, const double li, const double ui, const double soli)
    {
        BBL_TNodeBoundsSol<indType, numberType> el;
        
        el.ind = index;
        el.l = li;
        el.u = ui;
        el.sol = soli;
        
        BBL_addElementOnSortedNodeBoundsArray(size, a, el);
    }
    
    
    
    
    template <class indClass, class myClass>
    inline void BBL_getVarBoundsOnNodeFromNodeBounds( const unsigned int nElements, const BBL_TNodeBounds<indClass, myClass> *bounds, double *nlx, double *nux )
    {
        for(unsigned int i = 0; i < nElements; i++)
        {
            auto ind = bounds[i].ind;
            
            nlx[ind] = bounds[i].l;
            nux[ind] = bounds[i].u;
        }
    }
    
    
    template <class indClass, class myClass>
    inline void BBL_getVarIndicesOnNodeFromNodeBounds( const unsigned int nElements, const BBL_TNodeBounds<indClass, myClass> *bounds, unsigned int *indices )
    {
        for(unsigned int i = 0; i < nElements; i++)
        {
            indices[i] = bounds[i].ind;
        }
    }
    
    
    template <class indClass, class myClass>
    inline void BBL_getVarBoundsOnNodeFromArray( const unsigned int nElements, const BBL_Array< BBL_TNodeBounds<indClass, myClass> > &array, double *nlx, double *nux )
    {
        BBL_getVarBoundsOnNodeFromNodeBounds(nElements, array.a, nlx, nux);
    }
    
    
    template <class indClass, class myClass>
    inline void BBL_getVarBoundsOnNodeFromArraySize( const BBL_ArraySize< BBL_TNodeBounds<indClass, myClass> > &arraySize, double *nlx, double *nux )
    {
        BBL_getVarBoundsOnNodeFromNodeBounds( arraySize.size, arraySize.a, nlx, nux );
    }
    
    
    template <class indClass, class myClass>
    inline void BBL_getVarBoundsOnNodeFromNodeBoundsSol( const unsigned int nElements, const BBL_TNodeBoundsSol<indClass, myClass> *bounds, double *nlx, double *nux )
    {
        for(unsigned int i = 0; i < nElements; i++)
        {
            auto ind = bounds[i].ind;
            
            nlx[ind] = bounds[i].l;
            nux[ind] = bounds[i].u;
        }
    }
    
    
    
    
    
    template <class indClass, class myClass>
    inline void BBL_printNodeBounds(const unsigned int nElements, const BBL_TNodeBounds<indClass, myClass> *array, std::ostream &out )
    {
        for(unsigned i = 0; i < nElements; i++)
        {
            out << "\tvar: " << array[i].ind << " l: " << (double) array[i].l << " u: " << (double) array[i].u;
        }
    }
    
    
    template <class indClass, class myClass>
    inline void BBL_printNodeBoundsSol(const unsigned int nElements, const BBL_TNodeBoundsSol<indClass, myClass> *array, std::ostream &out )
    {
        for(unsigned i = 0; i < nElements; i++)
        {
            out << "\tvar: " << array[i].ind << " l: " << (double) array[i].l << " u: " << (double) array[i].u << " sol: " << array[i].sol;
        }
    }
    
    
    
    
    template <class indClass, class myClass>
    inline void BBL_printArraySizeNodeBounds( const BBL_ArraySize< BBL_TNodeBounds<indClass, myClass> > &arraySize, std::ostream &out )
    {
        BBL_printNodeBounds(arraySize.size, arraySize.a, out);
    }
    
    
    
    template <class pClass, class dClass>
    static inline void BBL_writeAndShift( pClass* &pointer, dClass data )
    {
        dClass *pd = (dClass*) pointer;
        pd[0] = data;
        pointer = (pClass*) (pd+1);
    }
    
    
    #if 0
    //this function write in the buffer and shift it!
    //variable index are written like 32-bits usnigned int and variable bounds like 64 bits double
    template <class indClass, class myClass>
    inline long unsigned int BBL_writeVarBoundsInaBufferArrayFromNodeBoundsAndShift(const unsigned int nElements, const BBL_TNodeBounds<indClass, myClass> *bounds, void* &buffer)
    {
        uint32_t ind;
        double l, u;
        long unsigned int bytesWritten = 0;
        
        
        for(unsigned int i = 0; i < nElements; i++)
        {
            ind = bounds[i].ind;
            l = bounds[i].l;
            u = bounds[i].u;
            
            BBL_writeAndShift(buffer, ind);
            bytesWritten += sizeof(ind);
            
            BBL_writeAndShift(buffer, l);
            bytesWritten += sizeof(l);
            
            BBL_writeAndShift(buffer, u);
            bytesWritten += sizeof(u);
        }
        
        return bytesWritten;
    }
    
    
    //variable index are written like 32-bits usnigned int and variable bounds like 64 bits double
    template <class indClass, class myClass>
    inline long unsigned int BBL_writeVarBoundsInaBufferArrayFromArrayAndShift(const unsigned int nElements, BBL_Array< BBL_TNodeBounds<indClass, myClass> > &array, void* &buffer)
    {
        return BBL_writeVarBoundsInaBufferArrayFromNodeBoundsAndShift(nElements, array.a, buffer);
    }
    
    
    //variable index are written like 32-bits usnigned int and variable bounds like 64 bits double
    template <class indClass, class myClass>
    inline long unsigned int BBL_writeVarBoundsInaBufferArrayFromArraySizeAndShift(BBL_ArraySize< BBL_TNodeBounds<indClass, myClass> > &array, void* &buffer)
    {
        return BBL_writeVarBoundsInaBufferArrayFromNodeBoundsAndShift(array.size, array.a, buffer);
    }
    
    #endif
    
    
    template <class indClass, class myClass>
    inline void BBL_buildInitialSolutionFromNodeBoundsSol( unsigned int nElements, const BBL_TNodeBoundsSol<indClass, myClass> *bounds, double *sol )
    {
        for(unsigned int i = 0; i < nElements; i++)
        {
            const unsigned int ind = bounds[i].ind;
            
            if( sol[ind] > bounds[i].u )
                sol[ind] = bounds[i].u;
            else if( sol[ind] < bounds[i].l )
                sol[ind] = bounds[i].l;
        }
    }
    
    
    /*This will work only for BBL_Node. Anyway, I want declare a template... :<) */
    template <class nodeClass>
    inline void BBL_appendToEncadeateList( nodeClass* &head, nodeClass *tail )
    {
        nodeClass *pnode;
        
        if(head == NULL)
        {
            head = tail;
            return;
        }
        
        if(tail == NULL)
            return;
        
        for( pnode = head;  pnode->next != NULL ; pnode = pnode->next);
        
        pnode->next = tail;
        tail->previous = pnode;
    }
    
    
}






#endif






























