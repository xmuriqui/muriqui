
#ifndef BBL_CONFIG_H_
#define BBL_CONFIG_H_

#include "../WAXM_config.h"

#define BBL_HAVE_CLOCK_GETTIME 	WAXM_HAVE_CLOCK_GETTIME
#define BBL_CPP_MULTITHREADING 	WAXM_CPP_MULTITHREADING

#define BBL_DEBUG_MODE WAXM_DEBUG_MODE


#define BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS 	1	//setting this option as 1 will save 4 or 8 bytes per node. It can be usefull to address large problems. However, you could not be able to create your own version of branch-and-bound node


#define BBL_USE_SHORT_INTEGER_TO_DEPTH_IN_BB_NODE 1 //setting this option as 1 will enable usage of short integers to store node deep. In general, short integers are enough to do the work in the most part of the case, but if you have a very large problem maybe you colud have a BB tree deeper than a short int can represent. Supposing 16 bits to short int, we could have in the maximum 65536 levels. Anyway, even in the case where your problem has more variables than it, all nodes can be pruned before reach this point. Moreover, usage of short ints only would affect the calculation of number of prunes by level (when user enable it). I suppose explorarton strategies would not be affect because we simulate data strutures to perform it (stack for depth and queue for width).



#define BBL_USE_SHORT_INTEGER_TO_NUMBER_OF_MY_BOUNDS_IN_BB_NODE 1	//setting this option as 1 will enable usage of short integers to and number of proper bounds (number of bounds set in the last branching) in the nodes (note, number of parent bounds still is a normal integer). In general, short integers are enough to do the work, but we can have some exception, for example, if you ahave the insane case where you are branching over more than 65536 variables in a single branching (supposin short integers having 16 bits).



#endif
