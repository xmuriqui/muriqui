


#ifndef MIP_CONSTANTS_HPP
#define MIP_CONSTANTS_HPP

#include "MIP_config.hpp"



namespace minlpproblem
{
    enum MIP_RETURN_CODE
    {
        MIP_SUCESS 					= 0,
        MIP_BAD_DEFINITIONS			= -1,
        MIP_BAD_VALUE				= -2,
        MIP_INDEX_FAULT				= -3,
        MIP_MEMORY_ERROR			= -4,
        MIP_REPETEAD_INDEXES		= -5,
        MIP_UPPER_TRIANGLE_INDEX	= -6,
        MIP_UNDEFINED_ERROR			= -7,
        MIP_CALLBACK_FUNCTION_ERROR	= -8,
        MIP_LIBRARY_NOT_AVAILABLE   = -9,
        MIP_INFEASIBILITY			= -10,
        MIP_NOT_APPLICABLE			= -11,
        MIP_FAILURE					= -12,
        
        MIP_RETURN_CODE_BEGIN 		= MIP_FAILURE,
        MIP_RETURN_CODE_END			= MIP_SUCESS
    };
    
    
    
    
    enum MIP_VARTYPE
    {
        MIP_VT_CONTINUOUS 	= 201,
        MIP_VT_INTEGER 		= 202,
        MIP_VT_CONTINGER	= 203, //special kind of variable to internal use
        //MIP_VT_BINARY = 	204
        
        MIP_VARTYPE_BEGIN 	= 201,
        MIP_VARTYPE_END		= 203
    };
    
    
    enum MIP_PROBLEMTYPE
    {
        MIP_PT_LP		= 140,
        MIP_PT_MILP,
        MIP_PT_QP,
        MIP_PT_MIQP,
        MIP_PT_QCP,
        MIP_PT_MIQCP,
        MIP_PT_NLP,
        MIP_PT_MINLP
    };
}


#endif
