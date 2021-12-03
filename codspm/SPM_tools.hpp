#ifndef SPM_TOOLS_HPP
#define SPM_TOOLS_HPP



#ifdef __FILE__
	#ifdef __LINE__
	    #define SPM_DEF_GETFILELINE 1
	#endif
#endif


#ifdef SPM_DEF_GETFILELINE

    #define SPM_GETFILELINE  \
	     " on file: '" __FILE__  "' function: '" << __func__ <<  "' line: " << __LINE__
#else
    #define SPM_GETFILELINE ""
#endif


#define SPM_PREPRINT "sparsematrix: "


#define SPM_PRINTNONSUPMSG std::cerr << SPM_PREPRINT "method: " << __func__ << " is not supported by " << typeid(*this).name() << "\n"


#define SPM_PRINTMEMERROR std::cerr << SPM_PREPRINT << "Memory error" << SPM_GETFILELINE << "\n"


#define SPM_PRINTERROR std::cerr << SPM_PREPRINT << "Error" << SPM_GETFILELINE << "\n"


#define SPM_PRINTERRORNUMBER(number) std::cerr << SPM_PREPRINT << "Error " << number << SPM_GETFILELINE << "\n"


#define SPM_PRINTERRORMSG(msg) std::cerr << SPM_PREPRINT << msg << SPM_GETFILELINE << "\n"

#define SPM_PRINTERRORMSGP(msg, arg) std::cerr << SPM_PREPRINT << msg << arg << SPM_GETFILELINE << "\n"

#define SPM_PRINTCALLBACKERRORNUMBER(number) std::cerr << SPM_PREPRINT << "Callback function error " << number << SPM_GETFILELINE << "\n"


#define SPM_PRINTINDEXERROR std::cerr << SPM_PREPRINT << "Index error" << SPM_GETFILELINE << std::endl

#define SPM_PRINTINDEXERRORP(arg) std::cerr << SPM_PREPRINT << "Index error " << arg << SPM_GETFILELINE << std::endl


#define SPM_getchar()  SPM_PRINTERRORMSG("stopped in a getchar"), getchar() 


#define SPM_LIBNOTAVAILABLERET(solverCode)  std::cerr << SPM_PREPRINT << "Solver '" << SPM_getSolverName(solverCode) << "' is not avaiable in this compilation. " << SPM_GETFILELINE << "\n" ; return SPM_LIBRARY_NOT_AVAILABLE


#define SPM_OPERATIONNOTIMPLEMENTEDRET(solverCode)  std::cerr << SPM_PREPRINT << "Operation " << __func__ << " is not implemented yet for '" << SPM_getSolverName(solverCode)  << SPM_GETFILELINE << "\n" ; return SPM_OPERATION_NOT_IMPLEMENTED


#define SPM_OPERATIONNOTSUPORTEDRET(solverCode)  std::cerr << SPM_PREPRINT << "Operation " << __func__ << " is not supported yet for '" << SPM_getSolverName(solverCode)  << SPM_GETFILELINE << "\n" ; return SPM_OPERATION_NOT_SUPPORTED




#define SPM_IFERRORRETURN(error, retCode) \
	if(error){ SPM_PRINTERRORNUMBER(error); return retCode;}

#define SPM_IFMEMERRORRETURN(errorExp) \
	if(errorExp){ SPM_PRINTMEMERROR; return SPM_MEMORY_ERROR;}

#define SPM_IFERRORGOTOLABEL(error, varToGetRetCode, retCode, destinLabel) \
	if(error){ SPM_PRINTERRORNUMBER(error); varToGetRetCode = retCode; goto destinLabel;  }

#define SPM_IFMEMERRORGOTOLABEL(errorExp, varToGetRetCode, destinLabel) \
	if(errorExp){ SPM_PRINTMEMERROR; varToGetRetCode = SPM_MEMORY_ERROR; goto destinLabel;  }

#define SPM_IFERRORSETVAR(error, varToGetRetCode, retCode) \
	if(error){ SPM_PRINTERRORNUMBER(error); varToGetRetCode = retCode; }

#define SPM_IFCALLBACKERRORSETVAR(error, varToGetRetCode) \
	if(error){ SPM_PRINTCALLBACKERRORNUMBER(error); varToGetRetCode = SPM_CALLBACK_FUNCTION_ERROR; }

#define SPM_IFCALLBACKERRORGOTOLABEL(error, varToGetRetCode, destinLabel) \
	if(error){ SPM_PRINTCALLBACKERRORNUMBER(error); varToGetRetCode = SPM_CALLBACK_FUNCTION_ERROR; goto destinLabel;  }

	

#define SPM_MY_ZERO_TO_EVAL 1e-10




namespace newspm{

#define SPM_DEBUG_MODE WAXM_DEBUG_MODE
#define SPM_EMAIL "wendelmelo@ufu.br"

enum SPM_RETURN_CODES
{
	SPM_MEMORY_ERROR 		= -101,
	SPM_ELEMENT_NOT_PRESENT	= -102,
	SPM_BAD_DEFINITIONS 	= -103,
	SPM_REPETEAD_ELEMENT 	= -104,
	SPM_BAD_VALUE 			= -105,
	SPM_UPPER_TRIANGLE 		= -106,
	SPM_INDEX_FAULT 		= -107,
	SPM_INTERNAL_ERROR		= -108
};


template <class myClass>
inline myClass SPM_abs( const myClass num )
{
    return ( num >= 0 ? num : -num );
}

template <class myClass>
inline myClass SPM_max( const myClass a, const myClass b )
{
    return (a >= b ? a : b);
}

template <class myClass>
inline myClass SPM_min( const myClass a, const myClass b )
{
    return (a <= b ? a : b);
}






template <class myClass1, class myClass2>
inline void SPM_swap( myClass1 &a, myClass2 &b )
{
	myClass1 ca = a;
	
	a = b;
	b = ca;
}


template <class myClass>
inline void SPM_secFree(myClass* &p)
{
	if(p)
	{
		free(p);
		p = NULL;
	}
}


template <class myClass>
inline void SPM_setAllArray( const unsigned int size, myClass *a, const myClass value )
{
	#pragma ivdep
	#pragma GCC ivdep
	for( unsigned int i = 0; i < size; i++ )
		a[i] = value;
}


template <class myClass>
inline void SPM_addToAllArray( const unsigned int size, myClass *a, const myClass value )
{
	#pragma ivdep
	#pragma GCC ivdep
	for( unsigned int i = 0; i < size; i++ )
		a[i] += value;
}


template <class myClass>
inline void SPM_multiplyAllArray( const unsigned int size, myClass *a, const myClass value )
{
	#pragma ivdep
	#pragma GCC ivdep
	for( unsigned int i = 0; i < size; i++ )
		a[i] *= value;
}


template <class type1, class type2>
inline void SPM_copyArray(const unsigned int size, const type1* source, type2 *destin)
{
	#pragma ivdep
	#pragma GCC ivdep
	for( unsigned int i = 0; i < size; i++)
		destin[i] = source[i];
}


//shift left the elements in array by "shift positions". Note shift cannot be negative, otherwise, you wanna a right shift
template <class myClass>
inline void SPM_shiftLeftArray(const unsigned int size, myClass *a, const unsigned int shift = 1)
{
	//we cannot use use unsigned int because the index will not be negative
	//for(int i = 0; i < size; i++)
		//a[i-shift] = a[i];
	
	myClass* const pa = &a[ -((int) shift) ];
	//do not use vectorization because pa is builder on a. So, we can have error here...
	for(unsigned int i = 0; i < size; i++)
		pa[i] = a[i];
}


//shift Right the elements in array by "shift positions". Note shift cannot be negative, otherwise, you wanna a left shift
template <class myClass>
inline void SPM_shiftRightArray(const unsigned int size, myClass *a, const unsigned int shift = 1)
{
	if( size > 0 )
	{
		myClass* const pa = &a[shift];
		//do not use vectorization because pa is builder on a. So, we can have error here...
		for(unsigned int i = size-1u; ; i--) //right shift. We start from end
		{
			pa[i] = a[i];
			if(i == 0) //i can be unsigned int
				break;
		}
	}
}


};





#endif
