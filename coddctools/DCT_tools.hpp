
#ifndef DCT_TOOLS_HPP
#define DCT_TOOLS_HPP


#include <cstdlib>
#include <cstddef>
#include <cassert>

#include <iostream>
#include <thread>

#include "DCT_dctools.hpp"

namespace dctools
{
    
    
    #ifdef __FILE__
        #ifdef __LINE__
            #define DCT_DEF_GETFILELINE 1
        #endif
    #endif


    #ifdef DCT_DEF_GETFILELINE
            
        #define DCT_GETFILELINE  \
            " on file: '" __FILE__  "' line: '" << __LINE__ <<  "' function: " << __func__
    #else
        #define DCT_GETFILELINE ""
    #endif


    # define DCT_PREPRINT "dctools: "

    
    #define DCT_OSPRINTMEMERROR(out) out << DCT_PREPRINT "Memory error" << DCT_GETFILELINE << "\n"
    
    #define DCT_PRINTMEMERROR DCT_OSPRINTMEMERROR(std::cerr)

    
    #define DCT_OSPRINTERROR(out) out << DCT_PREPRINT "Error" << DCT_GETFILELINE << "\n"
    
    #define DCT_PRINTERROR DCT_OSPRINTERROR(std::cerr)
    
    
    
    #define DCT_OSPRINTERRORNUMBER(out, number) out << DCT_PREPRINT "Error " << number << DCT_GETFILELINE << "\n"
    
    #define DCT_PRINTERRORNUMBER(number) DCT_OSPRINTERRORNUMBER(std::cerr, number)
    
    
    #define DCT_OSPRINTCALLBACKERRORNUMBER(out, number) out << DCT_PREPRINT << "Callback function error " << number << DCT_GETFILELINE << "\n"
    
    #define DCT_PRINTCALLBACKERRORNUMBER(number) DCT_OSPRINTCALLBACKERRORNUMBER(std::cerr, number)
    
    
    #define DCT_OSPRINTERRORMSG(out, msg) out << DCT_PREPRINT << msg << DCT_GETFILELINE << "\n"
    
    #define DCT_PRINTERRORMSG(msg) DCT_OSPRINTERRORMSG(std::cerr, msg)

    
    #define DCT_OSPRINTERRORMSGP(out, msg, arg) out << DCT_PREPRINT << msg << arg << DCT_GETFILELINE << "\n"
    
    #define DCT_PRINTERRORMSGP(msg, arg) DCT_OSPRINTERRORMSGP(std::cerr, msg, arg)

    
    #define DCT_OSPRINTMSG(out, msg) out << DCT_PREPRINT << msg
    
    #define DCT_PRINTMSG(msg) DCT_OSPRINTMSG(std::cout, msg)


    #define DCT_getchar()  DCT_PRINTERRORMSG("stopped in a getchar"), getchar()
    
    
    
    
    #define DCT_OSIFERRORRETURN(out, error, retCode) \
    if(error){ DCT_OSPRINTERRORNUMBER(out, error); return retCode;}
    
    #define DCT_OSIFMEMERRORRETURN(out, errorExp) \
        if(errorExp){ DCT_OSPRINTMEMERROR(out); return DCT_RC_MEMORY_ERROR;}

    #define DCT_OSIFERRORGOTOLABEL(out, error, varToGetRetCode, retCode, destinLabel) \
        if(error){ DCT_OSPRINTERRORNUMBER(out, error); varToGetRetCode = retCode; goto destinLabel;  }

    #define DCT_OSIFMEMERRORGOTOLABEL(out, errorExp, varToGetRetCode, destinLabel) \
        if(errorExp){ DCT_OSPRINTMEMERROR(out); varToGetRetCode = DCT_RC_MEMORY_ERROR; goto destinLabel;  }

    #define DCT_OSIFERRORSETVAR(out, error, varToGetRetCode, retCode) \
        if(error){ DCT_OSPRINTERRORNUMBER(out, error); varToGetRetCode = retCode; }

    #define DCT_OSIFCALLBACKERRORSETVAR(out, error, varToGetRetCode) \
        if(error){ DCT_OSPRINTCALLBACKERRORNUMBER(out, error); varToGetRetCode = DCT_OSCALLBACK_FUNCTION_ERROR; }

    #define DCT_OSIFCALLBACKERRORGOTOLABEL(out, error, varToGetRetCode, destinLabel) \
        if(error){ DCT_OSPRINTCALLBACKERRORNUMBER(out, error); varToGetRetCode = DCT_OSCALLBACK_FUNCTION_ERROR; goto destinLabel;  }
    
    
    #define DCT_IFERRORRETURN(error, retCode) \
        DCT_OSIFERRORRETURN(std::cerr, error, retCode)
                    
    #define DCT_IFMEMERRORRETURN(errorExp) \
        DCT_OSIFMEMERRORRETURN(std::cerr, errorExp)

    #define DCT_IFERRORGOTOLABEL(error, varToGetRetCode, retCode, destinLabel) \
        DCT_OSIFERRORGOTOLABEL(std::cerr, error, varToGetRetCode, retCode, destinLabel)

    #define DCT_IFMEMERRORGOTOLABEL(errorExp, varToGetRetCode, destinLabel) \
        DCT_OSIFMEMERRORGOTOLABEL(std::cerr, errorExp, varToGetRetCode, destinLabel)

    #define DCT_IFERRORSETVAR(error, varToGetRetCode, retCode) \
        DCT_OSIFERRORSETVAR(std::cerr, error, varToGetRetCode, retCode)

    #define DCT_IFCALLBACKERRORSETVAR(error, varToGetRetCode) \
        DCT_OSIFCALLBACKERRORSETVAR(std::cerr, error, varToGetRetCode)

    #define DCT_IFCALLBACKERRORGOTOLABEL(error, varToGetRetCode, destinLabel) \
        DCT_OSIFCALLBACKERRORGOTOLABEL(std::cerr, error, varToGetRetCode, destinLabel)
        
        
    
    
    
    
    template <class myClass>
    void DCT_malloc(myClass* &a, size_t size)
    {
        a = (myClass*) malloc( size * sizeof(myClass) );
    }
    
    
    template <class myClass>
    void DCT_calloc(myClass* &a, size_t size)
    {
        a = (myClass*) calloc(size, sizeof(myClass));
    }
    
    template <class myClass>
    static inline void DCT_secFree(myClass* &p)
    {
        if(p)
        {
            free(p);
            p = NULL;
        }
    }
    
    template <class myClass>
    static inline int DCT_realloc(myClass* &a, size_t size)
    {
        if(size == 0)
        {
            DCT_secFree(a);
        }
        else
        {
            void *p = realloc(a, size * sizeof(myClass) );
            if(!p)
            {
                DCT_PRINTMEMERROR;
                return DCT_RC_MEMORY_ERROR;
            }
            a = (myClass*) p;
        }
        
        return 0;
    }
    
    
    
    template <class myClass>
    static inline void DCT_secDelete(myClass* &p)
    {
        if(p)
        {
            delete p;
            p = NULL;
        }
    }
    
    
    template <class myClass>
    static inline void DCT_secDeleteArray(myClass* &p)
    {
        if(p)
        {
            delete[] p;
            p = NULL;
        }
    }
    


    template <class myClass1, class myClass2>
    static inline void DCT_copyArray(const unsigned int size, const myClass1 *origin, myClass2 *destin)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            destin[i] = origin[i];
    }
    
    
    template <class myClass>
    static inline void DCT_setAllArray( const unsigned int size, myClass *a, myClass value )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            a[i] = value;
    }
    
    
    template <class myClass>
    static inline void DCT_sumAllArray( const unsigned int size, myClass *a, myClass value )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            a[i] += value;
    }
    
    template <class myClass>
    static inline void DCT_multiplyAllArray(const unsigned int size, myClass *a, myClass value )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            a[i] *= value;
    }
    
    
    template <class myClass>
    static inline myClass DCT_swap(myClass &a, myClass &b)
    {
        myClass c = a;
        a = b;
        b = c;
    }
    
    
    template <class myClass>
    static inline myClass DCT_max(const myClass a, const myClass b ) 
    {
        return a >= b ? a : b;
    }
    
    template <class myClass>
    static inline myClass DCT_min(const myClass a, const myClass b)
    {
        return a <= b ? a : b;
    }
    
    template <class myClass>
    static inline myClass DCT_abs(const myClass a)
    {
        return a >= 0 ? a : -a;
    }
    
    
    //write data in the memory position pointed by pointer and shift the pointer for the memory position after data
    template <class pClass, class dClass>
    static inline void DCT_writeAndShift( pClass* &pointer, dClass data )
    {
        dClass *pd = (dClass*) pointer;
        pd[0] = data;
        pointer = (pClass*) (pd+1);
    }
    
    template <class pClass, class dClass>
    static inline void DCT_readAndShift( pClass* &pointer, dClass &data )
    {
        dClass *pd = (dClass*) pointer;
        data = pd[0];
        pointer = (pClass*) (pd+1);
    }
    
    
    /*inline unsigned int DCT_sizeofOpenNodeRep( unsigned int nVarBounds )
    {
        //connect code, lower bound, nbounds, and nbounds tuples of <index, lb, ub>
        return sizeof(DCT_Int32) + sizeof(double) + sizeof(DCT_UInt32) + nVarBounds*( sizeof(DCT_UInt32) + 2*sizeof(double) );
    }*/
    
    
    //count the size of string discarding white spaces in the beggining and in the end. \0 is discarded also
    inline unsigned int DCT_trimSize(const char *text)
    {
        const char *otext = text;
        unsigned int size = 0, size2 = 0;
        
        
        //skiping white spaces in the beggining 
        while(*text == ' ' || *text == '\t' || *text == '\n')
        {
            text++;
        }
        
        
        while(*text)
        {
            size++;
            text++;
        }
        
        
        //now, text point to the end of string
        text--;
        while(*text == ' ' || *text == '\t' || *text == '\n')
        {
            size2++;
            
            if(text == otext)
                break;
            
            text--;
        }
        
        #if DCT_DEBUG_MODE
            assert(size >= size2);
        #endif
        
        size -= size2;
        
        return size;
    }
    
    //copy string discarding white spaces in the beggining and in the end
    inline void DCT_copyTrimString(const char *source, char *destin)
    {
        const char *start = source;
        const char *end;
        
        if(source)
        {
            while(*start == ' ' || *start == '\t' || *start == '\n')
                start++;
            
            end = start;
            
            while(*end)
                end++;
            
            //now, end points to \0
            end--;
            while(*end == ' ' || *end == '\t' || *end == '\n')
            {
                if(end == start)
                    break;
                end--;
            }
            
            
            for(const char *p = start; p <= end; p++, destin++ )
                *destin = *p;
        }
        
        *destin='\0';
    }
    
    
    
    
    inline unsigned int DCT_getNumCores()
    {
        unsigned int r;
        
        #if DCT_CPP_MULTITHREADING
            r = std::thread::hardware_concurrency();
        #elif DCT_OMP_MULTITHREADING
            r = omp_get_num_procs();
        #else
            r = 1;
        #endif
        
        return r;
    }
    
    
    
    inline double DCT_getTime(void)
    {
        #if DCT_HAVE_CLOCK_GETTIME
            timespec Time;
            
            clock_gettime(CLOCK_REALTIME, &Time);
            
            return Time.tv_sec + Time.tv_nsec/1.0e9;
        #else
            return (double) time(NULL);
        #endif
    }
    
    
    
    
    
    
}



#endif

