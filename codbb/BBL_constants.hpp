#ifndef _BBL_CONSTANTS_HPP
#define _BBL_CONSTANTS_HPP


#include "BBL_config.hpp"
#include <cstdint>
#include <climits>
#include <cassert>


#if !BBL_CPP_MULTITHREADING
	#define BBL_OMP_MULTITHREADING 1
#endif

#if BBL_CPP_MULTITHREADING || BBL_OMP_MULTITHREADING
	#define BBL_MULTITHREADING 1
#endif





#define BBL_INFINITY  WAXM_INFINITY




namespace branchAndBound{
	
	enum BBL_RETURN_CODES
	{
		BBL_FEASIBLE_SOLUTION			= 2,
		BBL_OPTIMAL_SOLUTION			= 1,
		BBL_INFEASIBLE_PROBLEM			= -1,
		BBL_NO_FEASIBLE_SOLUTION_FOUND	= -2,
		BBL_MEMORY_ERROR 				= -3,
		BBL_UNDEFINED		 			= -4,
		BBL_UNDEFINED_ERROR 			= -5,
		BBL_MAX_TIME_STOP				= -6,
		BBL_MAX_ITERATIONS_STOP			= -7,
		//BBL_CALLBACK_FUNCTION_ERROR		= -8,
		BBL_STOP_REQUIRED_BY_USER		= -9,
		BBL_SOLVER_ERROR				= -10,
		BBL_BAD_DEFINITIONS				= -11,
		BBL_USER_CALLBACK_NOT_IMPLEMENTED = -12
	};
	
	
	enum BBL_EXP_STRATEGY
	{
		BBL_ES_DEPTH	= 30001,
		BBL_ES_WIDTH	= 30002,
		BBL_ES_BEST_LIMIT	= 30003
	};
	
	
	enum BBL_EXT_EXP_STRATEGY
	{
		BBL_EES_DEPTH	= BBL_ES_DEPTH,
		BBL_EES_WIDTH	= BBL_ES_WIDTH,
		BBL_EES_BEST_LIMIT	= BBL_ES_BEST_LIMIT,
		BBL_EES_DEPTH_BEST_LIMIT = BBL_ES_BEST_LIMIT + 1
	};
	
	
	enum BBL_BRANCH_STRATEGY
	{
		BBL_BS_USER_INDEX_CHOICE = 		40001,
		BBL_BS_USER_NODE_GENERATION = 	40002
	};
	
	
	enum BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY
	{
		BBL_PNBSS_SCHAR				= 50001,
		BBL_PNBSS_SHORT_INT,
		BBL_PNBSS_FLOAT,
		BBL_PNBSS_DOUBLE,
		
		
		BBL_PNBSS_SCHAR_SPECIAL_STRUCTS_TO_INT_VARS,		//not implemented yet
		BBL_PNBSS_SHORT_INT_SPECIAL_STRUCTS_TO_INT_VARS,	//not implemented yet
		BBL_PNBSS_FLOAT_SPECIAL_STRUCTS_TO_INT_VARS,		//not implemented yet
		BBL_PNBSS_DOUBLE_SPECIAL_STRUCTS_TO_INT_VARS,		//not implemented yet
		
		
		BBL_PNBSS_SINGLE_MIN = BBL_PNBSS_SCHAR,
		BBL_PNBSS_SINGLE_MAX = BBL_PNBSS_DOUBLE,
		
		BBL_PNBSS_ESP_MIN = BBL_PNBSS_SCHAR_SPECIAL_STRUCTS_TO_INT_VARS,
		BBL_PNBSS_ESP_MAX = BBL_PNBSS_DOUBLE_SPECIAL_STRUCTS_TO_INT_VARS,
	};
	
	
	enum BBL_UNION_TYPES_NODE_BOUNDS_POINTER
	{ /*Note: this enumeration should be stored in a char. So, make sure to attribute values that can be converted in 8 bits!!!!*/
		BBL_UTNBP_SHORT_SCHAR 		= 0,	//index: unsigned short int - bounds: signed char 
		BBL_UTNBP_SHORT_SHORT_INT 	= 1,	//index: unsigned short int - bounds: short int
		BBL_UTNBP_SHORT_FLOAT 		= 2,	//index: unsigned short int - bounds: float
		BBL_UTNBP_SHORT_DOUBLE 		= 3,	//index: unsigned short int - bounds: double
		
		BBL_UTNBP_SCHAR 			= 4,	//index: unsigned int - bounds: signed char
		BBL_UTNBP_SHORT_INT 		= 5,	//index: unsigned int - bounds: short int
		BBL_UTNBP_FLOAT 			= 6,	//index: unsigned int - bounds: float
		BBL_UTNBP_DOUBLE 			= 7,	//index: unsigned int - bounds: double
	};
	
	
	#if 0
	enum BBL_UNION_FLOAT_OR_DOUBLE_NODE_BOUNDS_POINTER
	{ /*Note: this enumeration should be stored in a char. So, make sure to attribute values that can be converted in 8 bits!!!!*/
		BBL_UFDNBP_FLOAT = 1,
		BBL_UFDNBP_DOUBLE = 2,
	};
	
	
	inline BBL_UNION_FLOAT_OR_DOUBLE_NODE_BOUNDS_POINTER BBL_pnbs2ufdnbp(BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY strategy)
	{
		if( strategy == BBL_PNBSS_FLOAT || strategy == BBL_PNBSS_FLOAT_SPECIAL_STRUCTS_TO_INT_VARS )
			return BBL_UFDNBP_FLOAT;
		else
			return BBL_UFDNBP_DOUBLE;
	}
	#endif
	
	
	inline BBL_UNION_TYPES_NODE_BOUNDS_POINTER BBL_pnbs2usfdnbp(BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY strategy, const unsigned int maxVars)
	{
		if( strategy == BBL_PNBSS_SCHAR || strategy == BBL_PNBSS_SCHAR_SPECIAL_STRUCTS_TO_INT_VARS )
		{
			if( maxVars <= USHRT_MAX + 1 )
				return BBL_UTNBP_SHORT_SCHAR;
			else
				return BBL_UTNBP_SCHAR;
		}
		else if( strategy == BBL_PNBSS_SHORT_INT || strategy == BBL_PNBSS_SHORT_INT_SPECIAL_STRUCTS_TO_INT_VARS )
		{
			if( maxVars <= USHRT_MAX + 1 )
				return BBL_UTNBP_SHORT_SHORT_INT;
			else
				return BBL_UTNBP_SHORT_INT;
		}
		else if( strategy == BBL_PNBSS_FLOAT || strategy == BBL_PNBSS_FLOAT_SPECIAL_STRUCTS_TO_INT_VARS )
		{
			if( maxVars <= USHRT_MAX + 1 )
				return BBL_UTNBP_SHORT_FLOAT;
			else
				return BBL_UTNBP_FLOAT;
		}
		else
		{
			#if BBL_DEBUG_MODE
				assert(strategy == BBL_PNBSS_DOUBLE || strategy == BBL_PNBSS_DOUBLE_SPECIAL_STRUCTS_TO_INT_VARS);
			#endif
			
			if( maxVars <= USHRT_MAX + 1 )
				return BBL_UTNBP_SHORT_DOUBLE;
			else
				return BBL_UTNBP_DOUBLE;
		}
	}
	
	
}




#endif
