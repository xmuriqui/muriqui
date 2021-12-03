


#include <math.h> //we include math.h to use INFINITY

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <climits>
#include <cassert>
#include <cstring>


#include <iostream>
#include <new>

#if BBL_CPP_MULTITHREADING
	#include <chrono> //only in c++ 2011
#endif



#include "BBL_branchAndBound.hpp"
#include "BBL_tools.hpp"


using namespace std;
using namespace branchAndBound;



#define BBL_SLEEP_THREADS_MILISECONDS 10



#define BBL_LOW_PRINTING 1
#define BBL_MEDIAN_PRINTING 4
#define BBL_HIGH_PRINTING 7





template <class myClass>
inline void BBL_printSequence( const unsigned int size, myClass &s )
{
	unsigned int i;
	
	for(i = 0; i < size; i++)
	{
		std::cout << *s;
		s++;
	}
}










BBL_HistorySolution::BBL_HistorySolution()
{
	nvars = iter = -1;
	time = cputime = -1.0;
	sol = NULL;
}


BBL_HistorySolution::BBL_HistorySolution(const int n, const long unsigned int iter, const double time, const double cputime, const double *sol, const double objF)
{
	int i;
	
	i = allocateSol(n);
	if(i != 0)
	{
		return;
	}
	
	this->nvars = n;
	this->iter = iter;
	this->time = time;
	this->cputime = cputime;
	
	for(i = 0; i < n; i++)
		this->sol[i] = sol[i];
	
	this->objF = objF;
	
}


BBL_HistorySolution::~BBL_HistorySolution()
{
	freeSolution();
}


int BBL_HistorySolution::allocateSol(const int n)
{
	sol = (double *) malloc( n * sizeof(double) );
	if( !sol )
	{
		#if BBL_DEBUG_MODE
			BBL_PRINTMEMERROR;
		#endif
		return BBL_MEMORY_ERROR;
	}
	
	return 0;
}


void BBL_HistorySolution::freeSolution()
{
	BBL_secFree(sol);
}


/*int BBL_HistorySolution::getnvars()
{
	return nvars;
}

int BBL_HistorySolution::getiter()
{
	return iter;
}

double BBL_HistorySolution::gettime()
{
	return time;
} */

int BBL_HistorySolution::getsolution(double *solution) const
{
	int i;
	
	if( sol == NULL )
	{
	    return -1;
	}
	else
	{
	    for(i = 0; i < nvars; i++)
			solution[i] = sol[i];
		
	    return 0;
	}
}






BBL_SolutionHistory::BBL_SolutionHistory()
{
	nsols = 0;
	hsols = NULL;
}


BBL_SolutionHistory::~BBL_SolutionHistory()
{
	desallocate();
}


void BBL_SolutionHistory::desallocate()
{
	int i;
	
	for(i = 0; i < nsols; i++)
	{
		if( hsols[i] )
			delete hsols[i];
	}
	
	BBL_secFree(hsols);
	
	nsols = 0;
}


/*int BBL_SolutionHistory::getnsols()
{
	return nsols;
} */


int BBL_SolutionHistory::addSolution(const int n, const long unsigned int iter, const double time, const double cputime, const double* sol, const double objF)
{
	BBL_HistorySolution **haux;
	
	haux = (BBL_HistorySolution **) realloc( hsols, (nsols + 1) * sizeof(BBL_HistorySolution *) );
	
	if(haux)
	{
		hsols = haux;
		hsols[nsols] = new (nothrow) BBL_HistorySolution(n, iter, time, cputime, sol, objF);
		
		if( !hsols[nsols] )
		{
			#if BBL_DEBUG_MODE
				BBL_PRINTMEMERROR;
			#endif
			return BBL_MEMORY_ERROR;
		}
		
		nsols++;
		return 0;
	}
	else
	{
		#if BBL_DEBUG_MODE
			BBL_PRINTMEMERROR;
		#endif
		return BBL_MEMORY_ERROR;
	}
}


//it is only a pointer, not a copy. Do not free. If you want a copy of solution use the method getsol of MRQ_HistorySolution pointer
BBL_HistorySolution * BBL_SolutionHistory::getHistSolPointer(const int index) const
{
	return hsols[index];
}





/*
int BBL_UserNodeGenerator::setx( BBL_Array< double >*& x, const double* sol, const double* dualSol) const
{
	const unsigned int size = (usePrimalSol ? n : 0) + ( useDualSol ? ndual : 0 );
	unsigned int ind;
	int code;
	
	
	x = new (nothrow) BBL_Array<double>;
		
	if(!x)
	{
		#if BBL_DEBUG_MODE
			BBL_PRINTMEMERROR;
		#endif
		code = BBL_MEMORY_ERROR;
		goto termination;
	}
	
	
	
	code = x->allocate( size );
	
	if(code != 0)
	{
		delete x;
		x = NULL;
		
		#if BBL_DEBUG_MODE
			BBL_PRINTMEMERROR;
		#endif
		code = BBL_MEMORY_ERROR;
		goto termination;
	}
	
	//x->incPointerCounter();
	
	if( usePrimalSol )
	{
		BBL_copySequence(n, sol, x->a);
		ind = n;
	}
	else
	{
		ind = 0;
	}
	
	if(useDualSol)
		BBL_copySequence(ndual, dualSol, &(x->a[ind]) );
	
	
	code = 0;
	
termination:

	return code;
}
*/


inline int BBL_UserNodeGenerator::insertNode( BBL_Node *node)
{
	int code;
	
	#if 0
	const bool userInitSol = (initSol || initDual);
	
	if(usePrimalSol || useDualSol)
	{
		if( userInitSol )
		{
			const double *ps = initSol ? initSol : sol;
			const double *ds = initDual ? initDual : dualSol;
			/*BBL_Array<double> *myx ;
			
			//if user provide a different primal or dual solution than optimal solution from the current node, we need set x for this son.
			
			code = setx( myx, ps, ds );
			
			if( code != 0 )
				goto termination; */
			
			r = node->setParentSol(n, ps, ndual, ds);
			BBL_IFERRORGOTOLABEL(r, code, r, termination);
		}
		else
		{
			/* //this node use the standard initial solution (relaxation's solution)
			
			if( !x )
			{
				code = setx(x, sol, dualSol);
				
				if( code != 0 )
					goto termination;
			}
			
			node->setxParentPointer(x); */
			
			r = node->setParentSol(n, sol, ndual, dualSol);
			BBL_IFERRORGOTOLABEL(r, code, r, termination);
		}
	}
	#endif
	
	
	
	
	//node->myBounds = boundsArray;
	//node->nMyBounds = nBounds;
	
	
	
	if(nnodes == 0)
	{
		nodes = node;
	}
	else
	{
		tail->next = node;
		node->previous = tail;
	}
	
	tail = node;
	
	nnodes++;
	
	code = 0;
	
	#if BBL_DEBUG_MODE
		assert( node->next == NULL);
	#endif
	
	
	
//termination:
	
	
	return code;
}





BBL_UserNodeGenerator::BBL_UserNodeGenerator(const unsigned int n, const unsigned int ndual, const bool use_parent_primal_sol, const bool use_parent_dual_sol, const BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY nodesType)//: parentBounds( BBL_pnbs2usfdnbp(nodesType) )
{
	this->n = n;
	this->ndual = ndual;
	
	this->usePrimalSol = use_parent_primal_sol;
	this->useDualSol = use_parent_dual_sol;
	
	this->nodesType = nodesType;
	
	this->parentBounds = NULL;
	this->parentBoundsNoInherit = NULL;
	
	auxlx = NULL;
	
}


BBL_UserNodeGenerator::~BBL_UserNodeGenerator()
{
	deallocate();
}




inline int BBL_UserNodeGenerator::generateBoundsArray( const double *nodelx, const double *nodeux, BBL_Node &node) const
{
	unsigned int i, j, nbounds;
	int r, code;
	
	
	for(i = 0, nbounds = 0; i < n; i++)
	{
		if(nodelx[i] != olx[i] || nodeux[i] != oux[i])
			nbounds++;
	}
	
	//boundsArray = (BBL_NodeBoundsSol *) malloc( nbounds * sizeof(BBL_NodeBoundsSol) );
	
	r = node.allocateNodeBounds(nbounds);
	BBL_IFERRORGOTOLABEL(r, code, r, termination);
	
	
	for(i = 0, j = 0; i < n; i++)
	{
		if( nodelx[i] != olx[i] || nodeux[i] != oux[i] )
		{
			/*boundsArray[j].ind = i;
			boundsArray[j].l = nodelx[i];
			boundsArray[j].u = nodeux[i];
			boundsArray[j].sol = sol[i];*/
			
			node.myBounds.setArrayElement(j, &i, &nodelx[i], &nodeux[i], &sol[i]);
			j++;
		}
	}
	
	
	code = 0;
	
termination:
	
	return code;
}



int BBL_UserNodeGenerator::generateNode( const double *nodelx, const double *nodeux, BBL_Node *newNode, double parentLowerBound, const double *initSol , const double *initDual )
{
	//BBL_ArraySize< BBL_NodeBounds > *boundsArray = NULL;
	//BBL_NodeBoundsSol *boundsArray = NULL;
	int r, code;
	//unsigned int nbounds;
	
	if(lb > parentLowerBound)
		parentLowerBound = lb;
	
	
	if( newNode == NULL )
	{
		newNode = new (std::nothrow) BBL_Node( BBL_pnbs2usfdnbp( nodesType, n) );
		BBL_IFMEMERRORGOTOLABEL(!newNode, code, desallocate_memory);
	}
	
	
	/* boundsArray = new (std::nothrow) BBL_ArraySize<BBL_NodeBounds>;
	if(!boundsArray)
	{
		#if BBL_DEBUG_MODE
			BBL_PRINTMEMERROR;
		#endif
		code = BBL_MEMORY_ERROR;
		goto desallocate_memory;
	} */
	
	if( parentBoundsNoInherit == NULL )
	{
		const unsigned int parentDepth = parent->getDepth();
		
		parentBoundsNoInherit = new (std::nothrow) BBL_BasePointer<BBL_ParentNodeInfo>( BBL_pnbs2usfdnbp( nodesType, n) );
		BBL_IFMEMERRORGOTOLABEL(!parentBoundsNoInherit, code, desallocate_memory);
		
		parentBoundsNoInherit->incPointerCounter();
		
		if( parentDepth < parentBoundsNoInherit->a.getMaxDepth() )
		{
			parentBoundsNoInherit->a.depth = parentDepth;
		}
		else
		{ //we reach the maximum depth representation
			parentBoundsNoInherit->a.depth = parentDepth - 1;
		}
	}
	
	if(usePrimalSol || useDualSol)
	{
		if(parentBoundsNoInherit->a.xParent == NULL )
		{
			const double *ps = initSol ? initSol : sol;
			const double *ds = initDual ? initDual : dualSol;
		
			r = parentBoundsNoInherit->a.setParentSol(n, ps, useDualSol ? ndual : 0, ds);
			BBL_IFERRORGOTOLABEL(r, code, r, desallocate_memory);
		}
	}
	
	if( parentBoundsNoInherit->a.lb < parentLowerBound )
		parentBoundsNoInherit->a.lb = parentLowerBound;
	
	
	
	newNode->setParentInfoPointer(parentBoundsNoInherit);
	
	
	r = generateBoundsArray( nodelx, nodeux, *newNode);
	BBL_IFERRORGOTOLABEL(r, code, r, desallocate_memory);
	
	
	r = insertNode(newNode);
	BBL_IFERRORGOTOLABEL(r, code, r, desallocate_memory);
	
	
	
	code = 0;
	
desallocate_memory:
	
	return code;
}




int BBL_UserNodeGenerator::generateNode( const unsigned int nNewBounds, const BBL_NodeBoundsSol *newBounds, const bool inheritParentBounds, const bool isNewBoundsAscOrdered, BBL_Node *newNode, double nodeLowerBound, const double *parentSol, const double *parentDualSol )
{
    //const unsigned int npMyBounds = parent->nMyBounds;
    //BBL_NodeBoundsSol *pMyBounds = parent->myBounds;

    bool nodeAlocated = false;
    //unsigned int i, j, k;
    int r, code;
    //unsigned int nbounds; //total number of bounds
    double realNodeLowerBound;
    BBL_NodeBoundsSol *boundsArray = NULL;

    //double *pl, *pu;

    //printf("generateNode. nNewBounds: %u\n", nNewBounds);

    if( lb > nodeLowerBound )
        nodeLowerBound = lb;

    realNodeLowerBound = BBL_max( nodeLowerBound, parent->getLowerBound() );

    if( newNode == NULL )
    {
        newNode = BBL_generateNewNode(nodesType, n); 
        BBL_IFMEMERRORGOTOLABEL(!newNode, code, desallocate_memory);
        
        nodeAlocated = true;
    }
    #if BBL_DEBUG_MODE
    else
    {
        assert( !newNode->previous && !newNode->next );
    }
    #endif





    if( this->parentBounds == NULL )
    {
        #if 0
        const BBL_ParentNodeBounds *parentBounds = &parent->getParentBounds()->a; //const BBL_NodeBounds* const parentBounds = parent->parentBounds->a;
        const unsigned int nParentBounds = parentBounds->getSize(); //const unsigned int nParentBounds = parent->parentBounds->size;
        
        int r;
        //BBL_NodeBounds *boundsArray = NULL;
        
        //i runs on parentBounds and j runs newBounds
        for(i = j = nbounds = 0; i < nParentBounds && j < npMyBounds; nbounds++)
        {
            unsigned int pbind;
            
            parentBounds->getArrayElement(i, &pbind, NULL, NULL);
            
            
            if(pbind < pMyBounds[j].ind) //if(parentBounds[i].ind < pMyBounds[j].ind)
            {
                i++;
            }
            else
            {
                if( pbind == pMyBounds[j].ind ) //if( parentBounds[i].ind == pMyBounds[j].ind )
                    i++;
                j++;
            }
        }
        
        if( i < nParentBounds )
        {
            nbounds += nParentBounds - i;
        }
        else //only i or j is lower than their respective limits
        {
            if( j < npMyBounds )
                nbounds += npMyBounds - j;
        }
        
        #if MRQ_DEBUG_MODE
            assert( nbounds <= nNewBounds + nDadBounds );
        #endif
        
        
        //boundsArray = new (std::nothrow) BBL_NodeBounds[nbounds];
        //BBL_IFMEMERRORGOTOLABEL( !boundsArray, code, desallocate_memory );
        
        //r = this->parentBounds.pointer.allocateArraySize();
        //BBL_IFERRORGOTOLABEL(r, code, r, desallocate_memory);
        
        this->parentBounds = new (std::nothrow) BBL_BasePointer<BBL_ParentNodeBounds>( BBL_pnbs2usfdnbp(nodesType) );
        
        r = this->parentBounds->a.allocateElements(nbounds);
        BBL_IFERRORGOTOLABEL(r, code, r, desallocate_memory);
        
        this->parentBounds->incPointerCounter();
        
        
        for(i = j = k = 0; i < nParentBounds && j < npMyBounds; k++)
        {
            unsigned int pbind;
            double pbl, pbu;
            
            parentBounds->getArrayElement(i, &pbind, &pbl, &pbu);
            
            
            if( pbind < pMyBounds[j].ind )
            {
                //boundsArray[k].ind = pbind;//parentBounds[i].ind;
                //boundsArray[k].l = pbl;//parentBounds[i].l;
                //boundsArray[k].u = pbu;//parentBounds[i].u;
                
                this->parentBounds->a.setArrayElement(k, &pbind, &pbl, &pbu);
                
                i++;
            }
            else
            {
                //boundsArray[k].ind = pMyBounds[j].ind;
                //boundsArray[k].l = pMyBounds[j].l;
                //boundsArray[k].u = pMyBounds[j].u;
                
                const unsigned int &pbind2 = pMyBounds[j].ind;
                const double &pbl2 = pMyBounds[j].l, &pbu2 = pMyBounds[j].u;
                
                this->parentBounds->a.setArrayElement(k, &pbind2, &pbl2, &pbu2);
                
                if(pbind == pMyBounds[j].ind)//if(parentBounds[i].ind == pMyBounds[j].ind)
                    i++;
                j++;
            }
        }
        
        
        if( i < nParentBounds )
        {
            j = nParentBounds;
            for( ; i < j; i++, k++ )
            {
                unsigned int pbind;
                double pbl, pbu;
                
                parentBounds->getArrayElement(i, &pbind, &pbl, &pbu);
                
                
                //boundsArray[k].ind = pbind;//parentBounds[i].ind;
                //boundsArray[k].l = pbl;//parentBounds[i].l;
                //boundsArray[k].u = pbu;//parentBounds[i].u;
                
                this->parentBounds->a.setArrayElement(k, &pbind, &pbl, &pbu);
            }
        }
        else
        {
            for( ; j < npMyBounds; j++, k++ )
            {
                unsigned int pbind = pMyBounds[j].ind;
                double pbl = pMyBounds[j].l, pbu = pMyBounds[j].u;
                
                //boundsArray[k].ind = pMyBounds[j].ind;
                //boundsArray[k].l = pMyBounds[j].l;
                //boundsArray[k].u = pMyBounds[j].u;
                
                this->parentBounds->a.setArrayElement(k, &pbind, &pbl, &pbu);
            }
        }
        
        
        //for( unsigned int w = 0; w < nbounds; w++ )
            //boundsArray[w].sol = sol[ boundsArray[w].ind ];
        
        /*this->parentBounds = new (std::nothrow) BBL_ArraySize<BBL_NodeBounds>;
        
        if( !this->parentBounds )
        {
            #if BBL_DEBUG_MODE
                BBL_PRINTMEMERROR;
            #endif
            code = BBL_MEMORY_ERROR;
            goto desallocate_memory;
        }
        
        this->parentBounds->setArray(boundsArray);
        this->parentBounds->size = nbounds;*/
        
        #endif
        
        //BBL_getchar();
        
        if( inheritParentBounds )
        {
            r = parent->generateParentInfoForChilds( BBL_pnbs2usfdnbp(nodesType, n), this->parentBounds);
            BBL_IFERRORGOTOLABEL(r, code, r, desallocate_memory);
        }
        else
        {
            const unsigned int parentDepth = parent->getDepth(); //note getDepth return the parent depth plus 1. Do not replace for parentInfo->a.depth
            
            parentBounds = new (std::nothrow) BBL_BasePointer<BBL_ParentNodeInfo>( BBL_pnbs2usfdnbp(nodesType, n) );
            BBL_IFMEMERRORGOTOLABEL(!parentBounds, code, desallocate_memory);
            
            parentBounds->incPointerCounter();
            
            if( parentDepth < parentBounds->a.getMaxDepth() )
            {
                parentBounds->a.depth = parentDepth;
            }
            else
            { //we reach the maximum depth representation. We just repeat the depth of the parent. Note parent->getDepth() returnd parent->parentInfo->a.depth + 1, so we have to decrease 1
                parentBounds->a.depth = parentDepth - 1;
            }
        }
    }


    if(usePrimalSol || useDualSol)
    {
        
        if( parentBounds->a.xParent == NULL )
        {
            const double *ps = parentSol ? parentSol : sol;
            const double *ds = parentDualSol ? parentDualSol : dualSol;
            
            r = parentBounds->a.setParentSol(n, ps, useDualSol ? ndual : 0, ds );
            BBL_IFERRORGOTOLABEL(r, code, r, desallocate_memory);
        }
    }



    if( realNodeLowerBound > parentBounds->a.lb )
        parentBounds->a.lb = realNodeLowerBound;


    newNode->setParentInfoPointer( this->parentBounds );


    

    /*boundsArray = (BBL_NodeBoundsSol *) malloc( nNewBounds * sizeof(BBL_NodeBoundsSol) );
    BBL_copySequence( nNewBounds, newBounds, boundsArray ); */

    
    
        
        
    if( !isNewBoundsAscOrdered )
    {
        //performing a buble sort to order indices...
        bool change;
        
        BBL_malloc(boundsArray, nNewBounds);
        BBL_IFMEMERRORGOTOLABEL( !boundsArray, code, desallocate_memory );
        
        BBL_copySequence( nNewBounds, newBounds, boundsArray );
        
        auto myNNewBounds = nNewBounds;
        
        do
        {
            change = false;
            
            for(unsigned int i = 1; i < myNNewBounds; i++)
            {
                if( boundsArray[i-1].ind > boundsArray[i].ind )
                {
                    BBL_swap( boundsArray[i-1], boundsArray[i] );
                    change = true;
                }
            }
            
            myNNewBounds--; //in each iteration k, the k-th greatest element is the right position in the array. So, we cannot check the last k positions more.
            
        }while(change);
    }
    
    
    {
        const BBL_NodeBoundsSol *pnewBounds = isNewBoundsAscOrdered ? newBounds : boundsArray;
        
        
        r = newNode->allocateNodeBounds( nNewBounds );
        BBL_IFERRORGOTOLABEL(r, code, r, desallocate_memory);

        for(unsigned int i = 0; i < nNewBounds; i++)
        {
            newNode->myBounds.setArrayElement(i, &pnewBounds[i].ind, &pnewBounds[i].l, &pnewBounds[i].u, &pnewBounds[i].sol );
        }
    }



    #if 0
    if( isNewBoundsAscOrdered )
    {
        
        /*if(inheritParentBounds)
        {
            
            //i runs on parentBounds and j runs newBounds
            for(i = j = nbounds = 0; i < nParentBounds && j < npMyBounds; nbounds++)
            {
                if(parentBounds[i].ind < pMyBounds[j].ind)
                {
                    i++;
                }
                else
                {
                    if( parentBounds[i].ind == pMyBounds[j].ind )
                        i++;
                    j++;
                }
            }
            
            if( i < nParentBounds )
            {
                nbounds += nParentBounds - i;
            }
            else //only i or j is lower than their respective limits
            {
                if( j < npMyBounds )
                    nbounds += npMyBounds - j;
            }
            
            #if MRQ_DEBUG_MODE
                assert( nbounds <= nNewBounds + nDadBounds );
            #endif
            
        }
        else
        {
            nbounds = npMyBounds;
        }
        
        
        boundsArray = (BBL_NodeBoundsSol *) malloc( nbounds * sizeof(BBL_NodeBoundsSol) );

        if( !boundsArray )
        {
            #if BBL_DEBUG_MODE
                BBL_PRINTMEMERROR;
            #endif
            code = BBL_MEMORY_ERROR;
            goto desallocate_memory;
        }
        
        
        for(i = j = k = 0; i < nParentBounds && j < npMyBounds; k++)
        {
            if( parentBounds[i].ind < pMyBounds[j].ind )
            {
                boundsArray[k].ind = parentBounds[i].ind;
                boundsArray[k].l = parentBounds[i].l;
                boundsArray[k].u = parentBounds[i].u;
                
                i++;
            }
            else
            {
                boundsArray[k] = pMyBounds[j];
                
                if(parentBounds[i].ind == pMyBounds[j].ind)
                    i++;
                j++;
            }
        }
        
        
        if( i < nParentBounds )
        {
            j = nParentBounds;
            for( ; i < j; i++, k++ )
            {
                boundsArray[k].ind = parentBounds[i].ind;
                boundsArray[k].l = parentBounds[i].l;
                boundsArray[k].u = parentBounds[i].u;
            }
        }
        else
        {
            for( ; j < npMyBounds; j++, k++ )
                boundsArray[k] = pMyBounds[j];
        }
        
        
        for( unsigned int w = 0; w < nbounds; w++ )
            boundsArray[w].sol = sol[ boundsArray[w].ind ];
        
        */
        
        
        
    }
    else
    {
        if( inheritParentBounds )
        {
            pl = nlx;
            pu = nux;
            
            BBL_copySequence(n, pl, auxlx);
            BBL_copySequence(n, pu, auxux);
            
            
            for(i = 0; i < npMyBounds; i++)
            {
                k = pMyBounds[i].ind;
                auxlx[k] = pMyBounds[i].l;
                auxux[k] = pMyBounds[i].u;
            }
            
            
            
            code = generateBoundsArray(auxlx, auxux, boundsArray, nbounds);
            
            if( code != 0 )
            {
                #if BBL_DEBUG_MODE
                    BBL_PRINTERRORNUMBER(code);
                #endif
                code = BBL_MEMORY_ERROR;
                goto desallocate_memory;
            }
            
        }
        else
        {
            //pl = olx;
            //pu = oux;
            
            
            nbounds = npMyBounds;
            
            boundsArray = (BBL_NodeBoundsSol *) malloc( nbounds * sizeof(BBL_NodeBoundsSol) );
            
            if( !boundsArray )
            {
                #if BBL_DEBUG_MODE
                    BBL_PRINTMEMERROR;
                #endif
                code = BBL_MEMORY_ERROR;
                goto desallocate_memory;
            }
            
            #pragma GCC ivdep
            for(unsigned int i = 0; i < nbounds; i++)
                boundsArray[i] = pMyBounds[i];
            
            
            #pragma GCC ivdep
            for(unsigned int i = 0; i < nbounds; i++)
                boundsArray[i].sol = sol[ boundsArray[i].ind ];
            
        }
        
        
    }
    #endif


    code = insertNode(newNode);

    if( code != 0 )
    {
        #if BBL_DEBUG_MODE
            BBL_PRINTERRORNUMBER(code);
        #endif
        goto desallocate_memory;
    }


    code = 0;

    desallocate_memory:


    if( code != 0 )
    {
        if( nodeAlocated )
        {
            delete newNode; 
        }
    }
    
    if(boundsArray)     free( boundsArray );


    return code;
}



unsigned int BBL_UserNodeGenerator:: getNumberOfNodes() const
{
	return nnodes;
}



BBL_Node* BBL_UserNodeGenerator::getNodesPointer() const
{
	return nodes;
}



int BBL_UserNodeGenerator::allocate(const unsigned int n)
{
	auxlx = (double *) malloc( 2*n*sizeof(double) );
	if(!auxlx)
	{
		#if BBL_DEBUG_MODE
			BBL_PRINTMEMERROR;
		#endif
		return BBL_MEMORY_ERROR;
	}
	
	auxux = &auxlx[n];
	
	return 0;
}



void BBL_UserNodeGenerator::deallocate()
{
	BBL_secFree(auxlx);
	
	if( parentBounds )
	{
		parentBounds->decPointerCounter();
		parentBounds = NULL;
	}
	
	if( parentBoundsNoInherit )
	{
		parentBoundsNoInherit->decPointerCounter();
		parentBoundsNoInherit = NULL;
	}
}



void BBL_UserNodeGenerator::initialize( BBL_Node *parentNode, double *olx, double *oux, double *nlx, double *nux, double lb, double *sol, double *dualSol)
{
	parent = parentNode;
	this->olx = olx;
	this->oux = oux;
	this->nlx = nlx;
	this->nux = nux;
	this->lb = lb;
	this->sol = sol;
	this->dualSol = dualSol;
	
	nnodes = 0;
	x = NULL;
	
	if(parentBounds)
	{
		parentBounds->decPointerCounter();
		parentBounds = NULL;
	}
	
	if(parentBoundsNoInherit)
	{
		parentBoundsNoInherit->decPointerCounter();
		parentBoundsNoInherit = NULL;
	}
	
	//parentBounds.setNodeBoundsPointer(NULL, NULL);
}









bool BBL_UserCallbacks::tryUpdateBestSolution(const int threadNumber, double* solution, const double fsolution, const long unsigned int iter)
{
	return bb->tryUpdateBestSolution(threadNumber, iter, solution, fsolution, *this);
}


bool BBL_UserCallbacks::tryUpdateLowerBound(const double lb)
{
	return bb->tryUpdateLowerBound(lb);
}


double BBL_UserCallbacks::getLowerBound() const
{
	return bb->zl;
}


long unsigned int BBL_UserCallbacks::getNodes(long unsigned int numberOfNodes, BBL_Node* &nodes, bool sortNodesByBound)
{
	return bb->getNodes(numberOfNodes, nodes, sortNodesByBound);
}


long unsigned int BBL_UserCallbacks::getNumberOfOpenNodes() const
{
	return bb->nodes->getNumberOfNodes();
}


bool BBL_UserCallbacks::getSomeBoundPrune() const
{
	return bb->someBoundPrune;
}


bool BBL_UserCallbacks::getSomeInfeasPrune() const
{
	return bb->someInfeasPrune;
}


bool BBL_UserCallbacks::getSomeOptPrune() const
{
	return bb->someOptPrune;
}


bool BBL_UserCallbacks::getSomeUserPrune() const
{
	return bb->someUserPrune;
}


double BBL_UserCallbacks::getUpperBound() const
{
	return bb->zu;
}


int BBL_UserCallbacks::getBestSolutionCopy(double *solution, double &fsolution) const
{
	return bb->getBestSolutionCopy(solution, fsolution);
}


double BBL_UserCallbacks::getBestSolutionObj() const
{
	return bb->getBestSolutionObj();
}


void BBL_UserCallbacks::getVariableBoundsArrayPointers( double*& lx, double*& ux ) const
{
	lx = bb->lx;
	ux = bb->ux;
}


bool BBL_UserCallbacks::hasFeasibleSolution() const
{
	return bb->out_best_obj < bb->in_infinity_value;
}


void BBL_UserCallbacks::printOpenNodesList() const
{
	if(bb)
	{
		if(bb->nodes)
			bb->nodes->print();
	}
}


void BBL_UserCallbacks::lockAuxiliaryThreadsBeforeBBLoop()
{
	bb->lockAuxThreadsBeforeLoop();
}


void BBL_UserCallbacks::unlockAllAuxiliaryThreads()
{
	bb->unlockAuxThreads();
}


BBL_BranchAndBound::BBL_BranchAndBound()
{
	out_best_sol = NULL;
	out_prune_counter_by_level = NULL;
	
	initializeProtectedData();
	
	resetParameters();
	resetOutput();
}


BBL_BranchAndBound::~BBL_BranchAndBound()
{
	desallocateThreadData();
	desallocateBestSol();
	desallocatePruneCountersByLevel();
	desallocateBoundArrays(); //not necessary, but ok
}


void BBL_BranchAndBound::resetOutput()
{
	desallocateBestSol();
	out_best_obj = INFINITY;
    out_first_sol_iter = -1; //highest possible value
	out_best_sol_iter = -1; //highest possible value
	
	out_feasible_sol = false;
	out_number_of_feas_sols = 0;
	out_number_of_iterations = 0;
	out_number_of_open_nodes = 0;
	out_number_of_threads = 1;
	out_return_code = BBL_UNDEFINED_ERROR;
	out_return_subcode = BBL_UNDEFINED;
	
	out_clock_time = out_cpu_time = 0.0;
    out_cpu_time_to_fisrt_sol = -1.0;
    out_clock_time_to_fisrt_sol = -1.0;
    out_cpu_time_to_best_sol = -1.0;
    out_clock_time_to_best_sol = -1.0;
	out_lower_bound = -INFINITY;
	out_upper_bound = INFINITY;
	out_obj_opt_at_root_relax = NAN;
	
	out_sol_hist.desallocate();
}


void BBL_BranchAndBound::resetParameters()
{
	in_call_after_bb_loop_callback = false;
	in_call_before_bb_loop_callback = false;
	in_call_before_solving_relax_callback = false;
	in_call_end_of_iteration_callback = false;
	in_call_new_best_solution_callback = false;
	in_call_updating_best_solution_callback = false;
	
	in_count_total_prunes = false;
	in_consider_relax_infeas_if_solver_fail = true;
	in_prune_nodes_by_bound = true;
	in_reorganize_lists = true;
	in_store_history_solutions = false;
	in_store_parent_dual_solution_on_nodes = false;
	in_store_parent_primal_solution_on_nodes = false;
	in_use_dual_obj_to_bound_prunning = false;
	
	in_lists_reorganization_frequency = 10000;
	in_max_tree_level_to_count_prunes_by_level = 0;
	in_number_of_node_sublists = 100;
	in_number_of_threads = 0;
	in_print_level = BBL_MEDIAN_PRINTING;
	in_printing_frequency = 10000;
	
	in_branching_strategy = BBL_BS_USER_INDEX_CHOICE;
	in_exp_strategy = BBL_EES_DEPTH_BEST_LIMIT;
	in_parent_node_bounds_strategy = BBL_PNBSS_DOUBLE;
	
	in_max_iterations = ULONG_MAX;
	
	in_absolute_convergence_tol = 1.0e-6;
	in_infinity_value = BBL_INFINITY;
	in_lower_bound = -INFINITY;
	in_upper_bound = INFINITY;
	in_max_time = INFINITY;
	in_max_cpu_time = INFINITY;
	
	in_relative_convergence_tol = 1.0e-4;
}


//that function only can be called by thread 0
void BBL_BranchAndBound::lockAuxThreadsBeforeLoop()
{
	if( __lockAuxThreadsBeforeLoop == false )
	{
		__lockAuxThreadsBeforeLoop = true;
		SEMAPH_lockAuxThs.lock(nthreads);
	}
}


void BBL_BranchAndBound::unlockAuxThreads()
{
	if( __lockAuxThreadsBeforeLoop )
	{
		__lockAuxThreadsBeforeLoop = false;
		SEMAPH_lockAuxThs.unlock(nthreads);
	}
}



int BBL_BranchAndBound::algorithmInitialization(const unsigned int n, const unsigned int m, const unsigned int ndual, const unsigned int nthreads, const double *lx, const double *ux)
{
	int ret;
	
	__lockAuxThreadsBeforeLoop = false;
	
	nodes = NULL;
	
	someBoundPrune = false;
	someInfeasPrune = false;
	someOptPrune = false;
	someUserPrune = false;
	
	endThreads = false;
	
	thRunning = NULL;
	thlbCurrentNode = NULL;
	thReturnCodes = NULL;
	thReturnSubCodes = NULL;
	
	this->n = n;
	this->m = m;
	this->ndual = ndual;
	
	zl = BBL_max( in_lower_bound, -in_infinity_value );
	zu = BBL_min( in_upper_bound, in_infinity_value );
	
	iter = 0;
	
	ret = allocateThreadData(nthreads);
	if(ret != 0)
	{
		#if BBL_DEBUG_MODE
			BBL_PRINTERRORNUMBER(ret);
		#endif
		return ret;
	}
	
	out_prune_counter.reset();
	
	
	if(in_max_tree_level_to_count_prunes_by_level > 0)
	{
		ret = allocatePruneCountersByLevel( in_max_tree_level_to_count_prunes_by_level );
		
		if(ret != 0)
		{
			#if BBL_DEBUG_MODE
				BBL_PRINTERRORNUMBER(ret);
			#endif
			return ret;
		}
	}
	
	
	ret = allocateBoundArrays(n);
	if(ret != 0)
		return ret;
	
	
	BBL_copySequence(n, lx, this->lx);
	BBL_copySequence(n, ux, this->ux);
	
	
	return allocateBestSol(n, m);
}



void BBL_BranchAndBound::algorithmFinalization( )
{
	
	out_feasible_sol = out_best_obj < in_infinity_value;
	out_number_of_iterations = iter;
	out_number_of_threads = nthreads;
	out_lower_bound = zl;
	out_upper_bound = zu;
	out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
	out_clock_time = BBL_getTime() - timeStart;

}


int BBL_BranchAndBound::allocateBestSol(const unsigned int n, const unsigned int m)
{
	BBL_malloc(out_best_sol, n);//out_best_sol = (double *) malloc( n * sizeof(double) );
	
	if( !out_best_obj )
	{
		#if BBL_DEBUG_MODE
			BBL_PRINTMEMERROR;
		#endif
		return BBL_MEMORY_ERROR;
	}
	
	for(unsigned int i = 0; i < n; i++)
		out_best_sol[i] = NAN;
	
	return 0;
}


int BBL_BranchAndBound::allocateBoundArrays(const unsigned int n)
{
	BBL_malloc(lx, n); //lx = (double *) malloc( n * sizeof(double) );
	BBL_malloc(ux, n);//ux = (double *) malloc( n * sizeof(double) );
	
	BBL_IFMEMERRORRETURN(!lx || !ux);
	
	return 0;
}


int BBL_BranchAndBound::allocatePruneCountersByLevel(const int nlevels)
{
	desallocatePruneCountersByLevel();
	
	out_prune_counter_by_level = new (nothrow) BBL_PruneCounter[nlevels];
	BBL_IFMEMERRORRETURN(!out_prune_counter_by_level);
	
	return 0;
}


int BBL_BranchAndBound::allocateThreadData(const unsigned int nthreads)
{
	BBL_calloc(thRunning, nthreads); //thRunning = (signed char *) calloc( nthreads, sizeof(thRunning) );
	BBL_malloc(thReturnCodes, nthreads); //thReturnCodes = (int *) malloc( nthreads * sizeof(int) );
	BBL_malloc(thReturnSubCodes, nthreads); //thReturnSubCodes = (int *) malloc( nthreads * sizeof(int) );
	BBL_malloc(thlbCurrentNode, nthreads); //thlbCurrentNode = (double *) malloc( nthreads * sizeof(double) );
	
	BBL_IFMEMERRORRETURN(!thRunning || !thReturnCodes || !thReturnSubCodes || !thlbCurrentNode);
	
	return 0;
}


bool BBL_BranchAndBound::checkTerminationCriterions(const int threadNumber, const double zl, const double zu, long unsigned int iter, int& retCode)
{
	//we already have a feasible solution...
	if(zu - zl < in_absolute_convergence_tol || zu - zl < BBL_abs(zu)*in_relative_convergence_tol)
	{
		if(out_best_obj < in_infinity_value)
		{
			retCode = BBL_OPTIMAL_SOLUTION;
			//return true;
		}
		else
		{
			retCode = BBL_NO_FEASIBLE_SOLUTION_FOUND; //BBL_INFEASIBLE_PROBLEM;
			//return true;
		}
		
		return true;
	}
	
	
	if( iter >= in_max_iterations )
	{
		retCode = BBL_MAX_ITERATIONS_STOP;
		return true;
	}
	
	
	if( in_max_cpu_time < INFINITY )
	{
		const double cpuTime = (double(clock()-clockStart) )/CLOCKS_PER_SEC;
		
		if( cpuTime >= in_max_cpu_time )
		{
			retCode = BBL_MAX_TIME_STOP;
			return true;
		}
	}
	
	
	if( in_max_time < INFINITY )
	{
		const double wallTime = BBL_getTime() - timeStart;
		
		if( wallTime >= in_max_time )
		{
			retCode = BBL_MAX_TIME_STOP;
			return true;
		}
	}
	
	
	return false;
}


void BBL_BranchAndBound::desallocateBestSol()
{
	BBL_secFree(out_best_sol);
}


void BBL_BranchAndBound::desallocateBoundArrays()
{
	BBL_secFree(lx);
	BBL_secFree(ux);
}


void BBL_BranchAndBound::desallocatePruneCountersByLevel()
{
	BBL_secDeleteArray( out_prune_counter_by_level );
}


void BBL_BranchAndBound::desallocateThreadData()
{
	BBL_secFree(thRunning);
	BBL_secFree(thReturnCodes);
	BBL_secFree(thReturnSubCodes);
	BBL_secFree(thlbCurrentNode);
}





int BBL_BranchAndBound::generateNodes(const unsigned int nBranchVars, unsigned int *branchVars, const double *breakValues1, const double *breakValues2,  BBL_Node *parent, const double *lxp, const double *uxp, const double nodelb, const double *solp, const double *dualSolp, const bool use_parent_primal_sol, const bool use_parent_dual_sol, BBL_Node* &nodes )
{
	const unsigned int totalNodes = 1 << nBranchVars; //(int) pow(2, nBranchVars) 
	//const BBL_BasePointer<BBL_ParentNodeBounds>* grandBounds = parent->getParentBounds();
	//const unsigned int nGrandBounds = grandBounds->a.getSize();
	
	bool repeat;
	int r, code;
	unsigned int i, ind;
	double vsol, vl, vu;
	BBL_Node *auxNode;
	
	BBL_BasePointer<BBL_ParentNodeInfo> *newParentInfo = NULL; 
	
	
	
	if( nodes == NULL )
	{
		nodes = BBL_generateNewNode( in_parent_node_bounds_strategy, n);
		BBL_IFMEMERRORGOTOLABEL(!nodes, code, termination);
	}
	
	
	//we check the quantity of bounds
	
	//sorting the indices of variable branching... sorry by buble sort...
	
	do
	{
		repeat = false;
		for(i = 1; i < nBranchVars; i++)
		{
			if(branchVars[i] < branchVars[i-1])
			{
				BBL_swap( branchVars[i], branchVars[i-1] );
				repeat = true;
			}
		}
		
	}while( repeat );
	
	
	#if 0
	
	newParentBounds = new (std::nothrow) BBL_BasePointer<BBL_ParentNodeBounds> ( BBL_pnbs2usfdnbp( in_parent_node_bounds_strategy) );
	BBL_IFMEMERRORGOTOLABEL(!newParentBounds, code, termination);
	
	newParentBounds->incPointerCounter();
	
	//i runs on parentBounds and j runs newBounds
	for(i = j = nb = 0; i < nGrandBounds && j < nParentBounds; nb++)
	{
		unsigned int ind;
		
		grandBounds->a.getArrayElement(i, &ind, NULL, NULL);
		
		if( ind < parentBounds[j].ind ) //if( grandBounds[i].ind < parentBounds[j].ind )
		{
			i++;
		}
		else
		{
			if( ind == parentBounds[j].ind ) //if( grandBounds[i].ind == parentBounds[j].ind )
				i++;
			j++;
		}
	}
	
	if( i < nGrandBounds )
	{
		nb += nGrandBounds - i;
		
		#if BBL_DEBUG_MODE
			assert( j == nParentBounds );
		#endif
	}
	else //only i or j is lower than their respective limits
	{
		#if BBL_DEBUG_MODE
			assert( i == nGrandBounds );
		#endif
		
		if( j < nParentBounds )
			nb += nParentBounds - j;
	}
	
	
	#if BBL_DEBUG_MODE
		assert( nb <= nGrandBounds + nParentBounds);
	#endif
	
	
	r = newParentBounds->a.allocateElements(nb);
	BBL_IFERRORGOTOLABEL(r, code, r, termination);
	
	
	//code = boundsArray->allocate(nb);
	//BBL_IFERRORGOTOLABEL(code, code, code, termination);
	
	
	//filling boundsArray...
	for(i = j = k = 0; i < nGrandBounds && j < nParentBounds; k++)
	{
		unsigned int gbind;
		
		grandBounds->a.getArrayElement(i, &gbind, NULL, NULL);
		
		//printf("parentBounds[%d].ind: %d  branchVars[%d]: %d     ", i, parentBounds[i].ind, j, branchVars[j]);
		
		if( gbind < parentBounds[j].ind ) //if( grandBounds[i].ind < parentBounds[j].ind )
		{
			ind = gbind; //grandBounds[i].ind;
			i++;
		}
		else
		{
			ind = parentBounds[j].ind;
			if( gbind == ind ) //if( grandBounds[i].ind == ind )
				i++;
			
			j++;
		}
		
		//printf("ind: %d lxDad[%d]: %f uxDad[%d]: %f\n", ind, ind, lxDad[ind], ind, uxDad[ind] );
		
		//boundsArray->a[k].ind = ind;
		//boundsArray->a[k].l = lxp[ind];
		//boundsArray->a[k].u = uxp[ind];
		
		newParentBounds->a.setArrayElement(k, &ind, &lxp[ind], &uxp[ind]);
	}
	
	
	if( i < nGrandBounds )
	{
		j = nGrandBounds;
		for(  ; i < j; i++, k++)
		{
			//ind = grandBounds[i].ind;
			
			grandBounds->a.getArrayElement(i, &ind, NULL, NULL);
			
			//boundsArray->a[k].ind = ind;
			//boundsArray->a[k].l = lxp[ind];
			//boundsArray->a[k].u = uxp[ind];
			
			newParentBounds->a.setArrayElement(k, &ind, &lxp[ind], &uxp[ind]);
		}
	}
	else
	{
		
		for( ; j < nParentBounds; j++, k++)
		{
			ind = parentBounds[j].ind;
			newParentBounds->a.setArrayElement(k, &ind, &lxp[ind], &uxp[ind]);
		}
		
	}
		
	//boundsArray->size = nb;
	#endif
	
	
	r = parent->generateParentInfoForChilds( BBL_pnbs2usfdnbp( in_parent_node_bounds_strategy, n), newParentInfo );
	BBL_IFERRORGOTOLABEL(r, code, r, termination);
	
	
	if( use_parent_primal_sol )
	{
		r = newParentInfo->a.setParentSol( n, solp, use_parent_dual_sol ? ndual : 0, dualSolp );
		BBL_IFERRORGOTOLABEL(r, code, r, termination);
	}
	
	if( nodelb > newParentInfo->a.lb )
		newParentInfo->a.lb = nodelb;
	
	/*if( use_parent_primal_sol )
	{
		BBL_copySequence(n, solp, x->a);
		ind = n;
	}
	else
	{
		ind = 0;
	}
	if( use_parent_dual_sol )
		BBL_copySequence(ndual, dualSolp, &(x->a[n]) ); */
	
	
	//MRQ_copyArray(nBranchVars, branchVars, lastBranch->a);
	//lastBranch->size = nBranchVars;
	
	
	
	
	
	//now, we allocate the nodes... 
	
	for(i = 0, auxNode = nodes;  ; )
	{
		r = auxNode->allocateNodeBounds( nBranchVars );
		BBL_IFERRORGOTOLABEL(r, code, r, termination);
		
		
		auto &myBounds = auxNode->myBounds;
		
		for(unsigned int j = 0; j < nBranchVars; j++)
		{
			ind = branchVars[j];
			vsol = solp[ind];
			
			//myBoundsp[j].ind = ind;
			//myBoundsp[j].sol = solp[ind];
			
			//bitwise operations
			//i has the number of the node
			//1 << j is the same: (int) pow(2, j)
			//test if the j-th bit in i is zero
			if( (i & (1 << j) ) == 0 )
			{
				//left branch...
				//myBoundsp[j].l = lxp[ind]; 
				//myBoundsp[j].u = breakValues1[j];
				
				vl = lxp[ind];
				vu = breakValues1[j];
			}
			else
			{
				//myBoundsp[j].l = breakValues2[j];
				//myBoundsp[j].u = uxp[ind];
				
				vl = breakValues2[j];
				vu = uxp[ind];
			}
			
			myBounds.setArrayElement(j, &ind, &vl, &vu, &vsol );
		}
		
		//auxNode->nMyBounds = nBranchVars;
		
		//auxNode->depth = newdepth;
		//auxNode->lb = lb;
		
		//if(use_parent_primal_sol || use_parent_dual_sol)
			//auxNode->setxParentPointer(x);
		
		//auxNode->setParentBoundsPointer(boundsArray);
		//auxNode->dadBounds->size = nb;
		auxNode->setParentInfoPointer(newParentInfo);
		
		
		i++;
		if( i >= totalNodes ) //could be ==
		{
			break;
		}
		
		
		if( auxNode->next == NULL )
		{
			auxNode->next = BBL_generateNewNode( in_parent_node_bounds_strategy, n ); //auxNode->next = new (nothrow) BBL_Node;
			BBL_IFMEMERRORGOTOLABEL(!(auxNode->next), code, termination);
			
			auxNode->next->previous = auxNode;
		}
		
		auxNode = (BBL_Node *) auxNode->next;
		
	}
	
	
	//user can have generated more nodes than necessary in nodes. So, we delete this unnecessary nodes
	
	
	for( BBL_Node *waste = auxNode->next ; waste ; )
	{
		BBL_Node *next = waste->next;
		delete waste;
		waste = next;
	}
	
	
	auxNode->next = NULL;
	
	
	code = 0;
	
termination:
	
	
	if( code != 0 )
	{
		auxNode = nodes;
		
		if( auxNode )
		{
			while(true)
			{
				if( auxNode->next )
				{
					auxNode = (BBL_Node *) auxNode->next;
					delete (BBL_Node *) auxNode->previous;
				}
				else
				{
					//end of the allocated list
					delete auxNode;
					break;
				}
			}
		}
		
		nodes = NULL;
	}
	
	if(newParentInfo)
		      newParentInfo->decPointerCounter();
	
	return code;
}


double  BBL_BranchAndBound::getBestSolutionObj()
{
	return out_best_obj;
}


long unsigned int BBL_BranchAndBound::getNodes(long unsigned int numberOfNodes, BBL_Node* &nodes, bool sortNodesByBound)
{
	return this->nodes->getNodes(numberOfNodes, nodes, sortNodesByBound);
}


int BBL_BranchAndBound::getBestSolutionCopy(double *solution, double &fsolution)
{
	if(solution)
	{
		const unsigned int n = this->n;
		
		SEMAPH_SolBounds.lock(nthreads);
			fsolution = out_best_obj;
			BBL_copySequence(n, out_best_sol, solution);
		SEMAPH_SolBounds.unlock(nthreads);
	}
	else
	{
		fsolution = getBestSolutionObj();
	}
	
	return 0;
}



inline void BBL_BranchAndBound:: incBoundPruneCounters( BBL_PruneCounter &pcounter, BBL_PruneCounter *pcounterlevel, const unsigned int depth )
{
	someBoundPrune = true;
	
	if(in_count_total_prunes)
		pcounter.bound++;
	
	if(in_max_tree_level_to_count_prunes_by_level > depth)
		pcounterlevel[depth].bound++;
}



inline void BBL_BranchAndBound:: incInfeasPruneCounters( BBL_PruneCounter &pcounter, BBL_PruneCounter *pcounterlevel, const unsigned int depth )
{
	someInfeasPrune = true;
	
	if(in_count_total_prunes)
		pcounter.infeas++;
	
	if(in_max_tree_level_to_count_prunes_by_level > depth)
		pcounterlevel[depth].infeas++;
}



inline void BBL_BranchAndBound:: incOptPruneCounters( BBL_PruneCounter &pcounter, BBL_PruneCounter *pcounterlevel, const unsigned int depth )
{
	someOptPrune = true;
	
	if(in_count_total_prunes)
		pcounter.opt++;
	
	if(in_max_tree_level_to_count_prunes_by_level > depth)
		pcounterlevel[depth].opt++;
}



inline void BBL_BranchAndBound:: incUserPruneCounters( BBL_PruneCounter &pcounter, BBL_PruneCounter *pcounterlevel, const unsigned int depth )
{
	someUserPrune = true;
	
	if(in_count_total_prunes)
		pcounter.user++;
	
	if(in_max_tree_level_to_count_prunes_by_level > depth)
		pcounterlevel[depth].user++;
}






void BBL_BranchAndBound::initializeProtectedData()
{
	ux = lx = NULL;
	//it is not necessaty, but we initialize thread data even so
	thRunning = NULL;
	thReturnCodes = NULL;
	thReturnSubCodes = NULL;
	thlbCurrentNode = NULL;
	
	nodes = NULL; //not really necessary also...
}



bool BBL_BranchAndBound::tryUpdateBestSolution( const unsigned int threadNumber, const long unsigned int iter, double* sol, double objValue, BBL_UserCallbacks& user_calbacks)
{
	bool updt = true, change = false;
	unsigned int i;
    decltype(out_number_of_feas_sols) nfeas;
	double obj = objValue; //user can change obj value. In this way, we preserve original objValue 
	double oldBestObj, zlaux;
    double clocktime, cputime;
	
	
	if( in_call_updating_best_solution_callback )
	{
		updt = user_calbacks.updatingBestSolution( threadNumber, sol, obj, zu, iter) == 0;
		
		//user can change the objecive value. I think, in general, user do not decrease obj value, only increase it. Anyway, we adopt the minimum value between old and new objective to try update zu.
		//if(obj < objValue)
        objValue = obj;
	}
	
	
	if(updt && objValue < zu )  //note, the strict correct is to compare objValue with zu inside the mutual exclusion zone. However, zu never increases its value. So, if objValue is already greater than zu, so, we do not need set semaphore to check in a strict correct way.
	{
        cputime = (double(clock()-clockStart) )/CLOCKS_PER_SEC;
        clocktime = BBL_getTime() - timeStart;
        
		SEMAPH_SolBounds.lock(nthreads);
		{
			//user can change objValue. We only put a possible new objective value in out_best_obj, and we let zu with the minimum value between original obj value and new obj value.
			//We need adopt the minimum value to guarantee correct bound prune in branch-and-bound.
			if(objValue < zu)
				zu = objValue;
			
			//obj can be greater than objValue (it can be have been corrected by user). So, we just use the new value to store out_best_obj.
			if( obj < out_best_obj )
			{
				BBL_copySequence(n, sol, out_best_sol);
				oldBestObj = out_best_obj;
				
				out_best_obj = obj;
				//zu = BBL_min(obj, objValue); 
				
				out_number_of_feas_sols++;
                
                out_best_sol_iter = iter;
                out_cpu_time_to_best_sol = cputime;
                out_clock_time_to_best_sol = clocktime;
                
                change = true;
				nfeas = out_number_of_feas_sols;
			}
		}
		SEMAPH_SolBounds.unlock(nthreads);
	}
	
	
	if(change)
	{
        if( nfeas == 1 ) //first feasible solution
        {
            out_first_sol_iter = iter;
            out_cpu_time_to_fisrt_sol = cputime;
            out_clock_time_to_fisrt_sol = clocktime;
        }
        
        
		if( in_store_history_solutions )
		{
			out_sol_hist.addSolution(n, iter, clocktime, cputime, sol, obj);
		}
		
		
		if( in_call_new_best_solution_callback )
			user_calbacks.newBestSolution( threadNumber, sol, oldBestObj, out_best_obj, iter); //we pass sol instead of out_bets_sol because out_best_sol can be changing by other thread
		
		
		zlaux = zl;
		
		if( in_prune_nodes_by_bound )
		{
			SEMAPH_DelNodes.lock(nthreads);
				i = nodes->pruneNodesByBound(zu, zlaux, in_max_tree_level_to_count_prunes_by_level);
			SEMAPH_DelNodes.unlock(nthreads);
			
			if(i > 0)
			{
				someBoundPrune = true;
				if(in_count_total_prunes)
				{
					SEMAPH_Attrib.lock(nthreads);
						out_prune_counter.bound += i;
					SEMAPH_Attrib.unlock(nthreads);
				}
			}
		}
		
		
		for(unsigned i = 0; i < nthreads; i++)
		{
		    if( thRunning[i] && thlbCurrentNode[i] < zlaux )
				zlaux = thlbCurrentNode[i];
		}
		
		
		SEMAPH_SolBounds.lock(nthreads);
			if( zlaux > zl )
				zl = zlaux;
		SEMAPH_SolBounds.unlock(nthreads);
		
		
		if(in_exp_strategy == BBL_EES_DEPTH_BEST_LIMIT && expStrategy == BBL_EES_DEPTH)
		{
			nodes->reorganizeToBestLimit(zlaux, 0.5*(zlaux + zu) );
			
			expStrategy = BBL_EES_BEST_LIMIT;
		}
	}
	
	
	
	return change;
}



bool BBL_BranchAndBound::tryUpdateLowerBound(const double lb)
{
	bool r = false;
	
	SEMAPH_SolBounds.lock(nthreads);
		if( lb > zl )
		{
			zl = lb;
			r = true;
		}
	SEMAPH_SolBounds.unlock(nthreads);
	
	return r;
}




//we need that function because c++11 threads only works on functions, not on methods... 
//that is ridiculus but you cannot pass arguments for thread as reference. Only copy or pointer...
int branchAndBound::BBL_bbLoop( BBL_UserCallbacks* userCalbacks, BBL_BranchAndBound* bb, const unsigned int number )
{
	return bb->bbLoop( *userCalbacks, number );
}





int BBL_BranchAndBound::run(const unsigned int nPrimalVars, const unsigned int nCons, const unsigned int nDualVars, const double *lx, const double *ux, branchAndBound::BBL_UserCallbacks& userCalbacks)
{
    int ret;
    unsigned int nroots = 1;
    BBL_Node *root = NULL;

    #if BBL_CPP_MULTITHREADING
        std::thread *myThreads = NULL;
    #endif




    timeStart = BBL_getTime();
    clockStart = clock();


    if( in_number_of_threads < 1 )
        nthreads = BBL_getNumCores();
    else
        nthreads = in_number_of_threads;

    //
    ret = algorithmInitialization(nPrimalVars, nCons, nDualVars, nthreads, lx, ux);
    if( ret != 0 )
    {
        BBL_PRINTERRORNUMBER(ret);
        
        out_return_code = ret;
        goto termination;
    }


    //if( in_print_level > 1 )
        //std::cout << std::endl << BBL_PREPRINT "Starting Branch-And-Bound algorithm" << std::endl << std::endl;

    nodes = new (std::nothrow) BBL_MTNodeListManager;
    if( !nodes )
    {
        if( in_print_level >= BBL_LOW_PRINTING )
            BBL_PRINTMEMERROR;
        
        out_return_code = BBL_MEMORY_ERROR;
        goto termination;
    }

    if( in_exp_strategy == BBL_EES_DEPTH_BEST_LIMIT )
    {
        if( zu < in_infinity_value )
            expStrategy = BBL_ES_BEST_LIMIT;
        else
            expStrategy = BBL_ES_DEPTH;
    }
    else
    {
        expStrategy = in_exp_strategy;
    }



    ret = nodes->allocateLists( expStrategy, nthreads, in_number_of_node_sublists, zl, zu, in_max_tree_level_to_count_prunes_by_level );

    if( ret != 0 )
    {
        if( in_print_level >= BBL_LOW_PRINTING )
            BBL_PRINTMEMERROR;
        
        out_return_code = BBL_MEMORY_ERROR;
        goto termination;
    }






    userCalbacks.bb = this;

    ret = userCalbacks.beforeAll( nthreads, this->lx, this->ux );

    if( ret != 0 )
    {
        if( in_print_level >= BBL_LOW_PRINTING )
            std::cerr << BBL_PREPRINT "Callback function userCalbacks.beforeAll returned error code " << ret << BBL_GETFILELINE  << "\n";
        
        out_return_code = BBL_STOP_REQUIRED_BY_USER;
        out_return_subcode = ret;
        
        goto termination;
    }



    //generating the root node(s)

    ret = userCalbacks.generateRootNode(root);
    if( ret != 0 )
    {
        if( in_print_level >= BBL_LOW_PRINTING )
            std::cerr << BBL_PREPRINT "Callback function userCalbacks.beforeAll returned error code " << ret << BBL_GETFILELINE << "\n";
        
        out_return_code = BBL_STOP_REQUIRED_BY_USER;
        out_return_subcode = ret;
        
        goto termination;
    }


    if( root )
    {
        for( BBL_Node *auxNode = root->next ; auxNode ; auxNode = auxNode->next )
            nroots++;
    }
    else
    {
        //generating the root node;
        root = BBL_generateNewNode(in_parent_node_bounds_strategy, n); //new (std::nothrow) BBL_Node;
    }


    if( !root )
    {
        if( in_print_level >= BBL_LOW_PRINTING )
            BBL_PRINTMEMERROR;
        
        out_return_code = BBL_MEMORY_ERROR;
        goto termination;
    }




    for( BBL_Node *node = root; node; node = node->next )
    {
        
        if( node->isParentInfoNull() )
        {
            BBL_BasePointer<BBL_ParentNodeInfo> *parentBounds = new (std::nothrow) BBL_BasePointer<BBL_ParentNodeInfo>( BBL_pnbs2usfdnbp(in_parent_node_bounds_strategy, n) );
            
            BBL_IFMEMERRORGOTOLABEL(!parentBounds, out_return_code, termination);
            
            //parentBounds.allocateArraySize();
            //parentBounds->incPointerCounter();
            
            //we consider root nodes as being non sibling; They do not share lastBranch...
            node->setParentInfoPointer(parentBounds);
            
            
            //we consider root nodes as being non sibling; They do not share lastBranch...
            /*node->parentBounds = new (std::nothrow) BBL_ArraySize <BBL_NodeBounds>;
            
            if( !node->parentBounds )
            {
                if( in_print_level >= BBL_LOW_PRINTING )
                    BBL_PRINTMEMERROR;
                
                out_return_code = BBL_MEMORY_ERROR;
                goto termination;
            }
            
            node->parentBounds->incPointerCounter();*/
        }
        
    }


    nodes->insertNodes(0, nroots, root);

    root = NULL;

    if( in_print_level > BBL_LOW_PRINTING ) //on low printing, we just print basic information
        //cout << "iter lb  ub  gap nodes at open" << endl;
        printf( BBL_PREPRINT "%-10s  %-14s  %-14s  %-14s  %-10s  %s\n", "iter", "lower bound", "upper bound", "gap", "open nodes", "thread");



    for(unsigned int i = 0; i < nthreads; i++) 
        thReturnCodes[i] = BBL_UNDEFINED; //initializing the thread return codes...

    for(unsigned int i = 0; i < nthreads; i++) 
        thReturnSubCodes[i] = BBL_UNDEFINED; //initializing the thread return codes...

    for(unsigned int i = 0; i < nthreads; i++)
        thlbCurrentNode[i] = -INFINITY;

    endThreads = false;




    #if BBL_CPP_MULTITHREADING
        
        if( nthreads > 1 )
        {
            myThreads = new (nothrow) thread[nthreads-1];
            
            if( !myThreads )
            {
                if(in_print_level >= BBL_LOW_PRINTING)
                    BBL_PRINTMEMERROR;
                
                out_return_code = BBL_MEMORY_ERROR;
                goto termination;
            }
            
            
            for(BBL_uint i = 1; i < nthreads; i++)
            {
                //that is ridiculus, but you cannot pass arguments for thread as reference. Only copy or pointer...
                myThreads[i-1] = std::thread( BBL_bbLoop, &userCalbacks, this, i);
            }
        }
        
        
        //using the main thread also
        bbLoop(userCalbacks, 0);
        
        for(BBL_uint i = 1; i < nthreads; i++)
            myThreads[i-1].join();
        
        
    #elif BBL_OMP_MULTITHREADING
        
        omp_set_num_threads( nthreads );
        
        #pragma omp parallel
        {
            bbLoop( userCalbacks, omp_get_thread_num() );
        }
        
    #else

        bbLoop( userCalbacks, 0);
        
    #endif


    out_return_code = BBL_UNDEFINED;
    for(unsigned i = 0; i < nthreads; i++)
    {
        if( thReturnCodes[i] != BBL_UNDEFINED )
        {
            //we adopt the first non BBL_UNDEFINED_ERROR return code 
            out_return_code = thReturnCodes[i];
            out_return_subcode = thReturnSubCodes[i];
            break;
        }
    }



    termination:


    for( BBL_Node *node = root; node; )
    {
        BBL_Node *next = node->next;
        
        delete node;
        node = next;
    }


    if( nodes )
    {
        if( in_max_tree_level_to_count_prunes_by_level > 0 )
            nodes->acumulatePruneLevelCounters( in_max_tree_level_to_count_prunes_by_level, out_prune_counter_by_level);
        
        
        out_number_of_open_nodes = nodes->getNumberOfNodes();
        
        delete nodes;
    }


    #if BBL_CPP_MULTITHREADING
        if( myThreads )
            delete[] myThreads;
    #endif



    desallocateThreadData();

    desallocateBoundArrays();

    algorithmFinalization();


    userCalbacks.afterAll( out_number_of_iterations, out_cpu_time, out_clock_time, out_lower_bound, out_upper_bound);


    if(in_print_level > BBL_LOW_PRINTING)
        cout << "cpu time: " << out_cpu_time << endl;


    return out_return_code;
}





int BBL_BranchAndBound::bbLoop( BBL_UserCallbacks& userCallbacks, const unsigned int thnumber)
{
    bool prune;
    unsigned int ui, uj;
    BBL_BRANCH_STRATEGY branchStrat = in_branching_strategy;
    int ret;
    unsigned long int  myiter, threadIter = 0;
    double zlaux, zuaux;
    double objValue, dualObjValue, nodelb;

    unsigned int *indices = NULL;

    double *nlx = NULL, *nux = NULL;
    double *sol = NULL, *dualSol = NULL;

    double *values1 = NULL, *values2 = NULL;

    BBL_RETURN_CODES retCode;
    BBL_Node *auxBBLNodeBase;
    BBL_Node *currNode = NULL; //we put NULL because sometimes we abort BB procedure if user request before desallocate current Node. So, we need test at end...
    BBL_UserNodeGenerator *userNodeGen = NULL;

    BBL_PruneCounter  *pruneCounterLevel = NULL, pruneCounter;


    //std::cout << "Thread " << thnumber << " entrou no bb.loop" << endl;

    if( __lockAuxThreadsBeforeLoop && thnumber > 0 )
    {
        SEMAPH_lockAuxThs.lock(nthreads); //we just put it here to stop all threads except thread 0. In some point,  thread 0 will unlock this semaphore and each thread will unlock the next...
        SEMAPH_lockAuxThs.unlock(nthreads);
    }

    if( endThreads )
    { // we need it here because some thread could be locked by SEMAPH_lockAuxThs until the end of execution... when thread is unlock execution can have been finished.
        thReturnCodes[ thnumber ] = BBL_UNDEFINED;
        goto termination;
    }


    //std::cout << "Thread " << thnumber << " passou do semaforo " << endl;

    BBL_malloc( nlx, 2*n );
    BBL_malloc( sol,  n + ndual );
    BBL_IFMEMERRORGOTOLABEL(!nlx || !sol, thReturnCodes[thnumber], termination);
    

    nux = &nlx[n];
    dualSol = &sol[n];


    //now, we allocate userNodeGenerator because user can change branching strategy during the BB procedure
    //if( in_branching_strategy == BBL_BS_USER_NODE_GENERATION)
    {
        userNodeGen = new (std::nothrow) BBL_UserNodeGenerator( n, ndual, in_store_parent_primal_solution_on_nodes, in_store_parent_dual_solution_on_nodes, in_parent_node_bounds_strategy );
        BBL_IFMEMERRORGOTOLABEL(!userNodeGen, thReturnCodes[thnumber], termination);
        
        ret = userNodeGen->allocate(n);
        BBL_IFERRORGOTOLABEL(ret, thReturnCodes[thnumber], ret, termination);
    }
    //else
    {
        BBL_malloc(indices, n); //indices = (unsigned int *) malloc( n * sizeof(unsigned int) );
        BBL_malloc(values1, n); //values1 = (double *) malloc( n * sizeof(double) );
        BBL_malloc(values2, n); //values2 = (double *) malloc( n * sizeof(double) );
        
        BBL_IFMEMERRORGOTOLABEL(!indices || !values1 || !values2, thReturnCodes[thnumber], termination);
    }
    
        

    if(in_max_tree_level_to_count_prunes_by_level > 0)
    {
        pruneCounterLevel = new (nothrow) BBL_PruneCounter[ in_max_tree_level_to_count_prunes_by_level ];
        
        BBL_IFMEMERRORGOTOLABEL(!pruneCounterLevel, thReturnCodes[thnumber], termination);
    }


    if( in_call_before_bb_loop_callback )
    {
        int ret = userCallbacks.beforeBBLoop(thnumber, zl, zu);
        
        if( ret != 0 )
        {
            if( in_print_level >= BBL_LOW_PRINTING )
                std::cerr << BBL_PREPRINT "Callback function userCalbacks.beforeBBLoop returned error code " << ret << BBL_GETFILELINE << "\n";
            
            thReturnCodes[thnumber] = BBL_STOP_REQUIRED_BY_USER;
            thReturnSubCodes[thnumber] = ret;
            
            goto termination;
        }
    }



    //starting bb loop
    while( true )
    {
        
        //printf("Oi 1 Thread: %d\n", thnumber);
        //fflush(stdout);
        
        //getting a node to exploit
        while( true )
        {
            double value;
            
            if( endThreads )
            {
                thReturnCodes[ thnumber ] = BBL_UNDEFINED;
                goto termination;
            }
            
            if( in_prune_nodes_by_bound )
            {
                zuaux = zu;
                value = BBL_zuWithTol(zuaux, in_absolute_convergence_tol, in_relative_convergence_tol);
            }
            else
            {
                value = INFINITY;
            }
            
            
            //trying to get a node to exploit
            long unsigned int longuint = nodes->getNodePointer( value, auxBBLNodeBase, in_max_tree_level_to_count_prunes_by_level );
            
            if(longuint > 0)
            {
                someBoundPrune = true;
                if(in_count_total_prunes)
                    pruneCounter.bound += longuint;
            }
            
            if( auxBBLNodeBase )
            {
                break;
            }
            
            //we failed to get a node to exploit
            
            if( thnumber == 0 && __lockAuxThreadsBeforeLoop )
            { //so, we do have no more nodes to exploit and other threads are blocked before the bb loop. In this case, we must stop the algorithm because there are no more possible nodes to exploit. Note, only thread 0 can change flag __lockAuxThreadsBeforeLoop. So, I think we must not have problems about mutual exclusion on __lockAuxThreadsBeforeLoop here.
                goto endLoop;
            }
            
            
            //thRunning[thnumberOrder] can be 0 ou -1 in this case. We only give up the threads is all thRunning[thnumberOrder] are -1. In this way, all threads should pass by here at least 2 times before end. In this way, we garantee an incorrect giving up in a strange situation...
            
            thlbCurrentNode[thnumber] = in_infinity_value;
            
            
            if( thRunning[thnumber] == 1 )
            {
                if( nthreads == 1 )
                    goto endLoop;
                
                thRunning[thnumber] = 0;
                uj = 1;
            }
            else
            { //thRunning[thnumber] == 0 or == -1
            
                //checking if some thread is still running
                bool allNonPositive = true;
                for( ui = 0; ui < nthreads; ui++)
                {
                    //we compare if thRunning is non negative. In this way, we guarantee all threads passed by here at least twice.
                    
                    if( thRunning[ui] > 0 )
                    {
                        allNonPositive = false;
                        break;
                    }
                }
                
                //only if all threads has thRunning as 0 or -1, we put this thRunning as -1.
                
                if( allNonPositive )
                {
                    thRunning[thnumber] = -1; //we just put -1 if all other threads have thRunning as 0 or -1.
                    
                    ////only if all threads have thRunning as -1 we end B&B
                    
                    uj = 0;  
                    for( ui = 0; ui < nthreads; ui++)
                    {
                        if( thRunning[ui] > -1 )
                        {
                            uj = 1;
                            break;
                        }
                    }
                }
                else
                {
                    uj = 1;
                }
            }
            
            
            
            if(uj == 1)
            {//we still have theads running. So we sleep to wait new nodes in the list
                
                #if BBL_CPP_MULTITHREADING
                    
                    //std::this_thread::sleep_for(  std::chrono::milliseconds( BBL_SLEEP_THREADS_MILISECONDS)  );
                    
                #elif BBL_OMP_MULTITHREADING
                    //nanosleep(&tim, NULL); I would like to use nanosleep, but it only works on POSIX systems
                #else
                    goto endLoop; //we have only 1 thread. If there is no nodes more to exploit, we end branch-and-bound
                #endif
            }
            else
            {
                goto endLoop;
            }
            
            
        }
        
        //printf("Oi 2 Thread: %d\n", thnumber);
        //fflush(stdout);
        
        thRunning[ thnumber ] = 1; //we let it here to avoid a possible incorret giving up of threads...
        
        currNode = (BBL_Node *) auxBBLNodeBase;
        
        thlbCurrentNode[ thnumber ] = currNode->getLowerBound(); //do not use currNode->heurlb because the nodes are not ordered by this value. Otherwise, the algorithm will set currNode->heurlb as the global lower bound incorrectly...
        
        
        if( expStrategy == BBL_ES_BEST_LIMIT )
        {
            zlaux = thlbCurrentNode[0];
            
            for(BBL_uint i = 1; i < nthreads; i++)
            {
                if( thlbCurrentNode[i] < zlaux )
                    zlaux = thlbCurrentNode[i];
            }
            
            if( zlaux > zl )
            {
                SEMAPH_SolBounds.lock(nthreads);
                    if( zlaux > zl )
                        zl = zlaux;
                SEMAPH_SolBounds.unlock(nthreads);
            }
        }
        
        
        
        SEMAPH_Attrib.lock(nthreads);
            iter++;
            myiter = iter;
        SEMAPH_Attrib.unlock(nthreads);
        
        threadIter++;
        
        
        BBL_copySequence(n, lx, nlx);
        BBL_copySequence(n, ux, nux);
        currNode->getVarBoundsOnNode(nlx, nux);
        
        
        prune = false;
        
        if( in_call_before_solving_relax_callback )
        {
            ret = userCallbacks.beforeSolvingRelaxation( thnumber, *currNode, myiter, zl, zu, nlx, nux, prune);
            
            if(prune)
                incUserPruneCounters(pruneCounter, pruneCounterLevel, currNode->getDepth());
            
            if(ret != 0)
            {
                if( in_print_level >= BBL_LOW_PRINTING )
                    std::cerr << BBL_PREPRINT "Callback function userCalbacks.beforeSolvingRelaxation returned error code " << ret << BBL_GETFILELINE << endl;
                
                thReturnCodes[thnumber] = BBL_STOP_REQUIRED_BY_USER;
                thReturnSubCodes[thnumber] = ret;
                
                goto termination;
            }
        }
        
        //printf("Oi 3 Thread: %d\n", thnumber);
        //fflush(stdout);
        
        if( !prune )
        {
            bool genFeasSol = false;
            retCode = BBL_UNDEFINED;
            objValue = dualObjValue = NAN;
            nodelb = currNode->getLowerBound();
            
            ret = userCallbacks.solveSubProblem(thnumber, *currNode, myiter, zl, zu, nlx, nux, retCode, objValue, dualObjValue, sol, dualSol, genFeasSol, prune, nodelb, branchStrat);
            
            //cout << "1.1 prune: " << prune << endl;
            
            if(ret != 0)
            {
                if( in_print_level >= BBL_LOW_PRINTING )
                    std::cerr << BBL_PREPRINT "Callback function userCalbacks.solveSubProblem returned error code " << ret << BBL_GETFILELINE << "\n";
                
                thReturnCodes[thnumber] = BBL_STOP_REQUIRED_BY_USER;
                thReturnSubCodes[thnumber] = ret;
                
                goto termination;
            }
            
            //for(unsigned int i = 0; i < n; i++)
                //cout << "sol[" << i << "]: " << sol[i] << endl;
            
            if( myiter == 1 )
            {
                out_obj_opt_at_root_relax = in_use_dual_obj_to_bound_prunning ? dualObjValue : objValue;
                
                if( retCode == BBL_OPTIMAL_SOLUTION )
                {
                    if( nodes->getNumberOfNodes() == 0 )
                    {
                        tryUpdateLowerBound( out_obj_opt_at_root_relax );
                    }
                }
            }
            
            if( prune )
            {
                if( in_max_tree_level_to_count_prunes_by_level > currNode->getDepth() )
                    incUserPruneCounters(pruneCounter, pruneCounterLevel, currNode->getDepth());
            }
            else
            {
                nodelb = BBL_max( nodelb, in_use_dual_obj_to_bound_prunning || retCode != BBL_OPTIMAL_SOLUTION ? dualObjValue : objValue );
                
                zuaux = zu;
                zuaux = BBL_zuWithTol( zuaux, in_absolute_convergence_tol, in_relative_convergence_tol );
                
                //cout << "1.3 prune: " << prune << " nodelb: " <<  nodelb << " zu: " << zu << endl ;
                
                if( in_prune_nodes_by_bound && nodelb >= zuaux )
                {
                    prune = true;
                    incBoundPruneCounters( pruneCounter, pruneCounterLevel, currNode->getDepth() );
                }
                else
                {
                    
                    if( genFeasSol )
                        tryUpdateBestSolution(thnumber, myiter, sol, objValue, userCallbacks);
                    
                    
                    if( retCode == BBL_OPTIMAL_SOLUTION  )
                    {
                        
                        if( genFeasSol )
                        {
                            //that is not ideal, zu can have been updated by other thread and this prune would be by bound, not by optimallity, but that is the best that we can do now, and, in practical terms, there is no difference about the kind of prune...
                            zuaux = zu;
                            zuaux = BBL_zuWithTol( zuaux, in_absolute_convergence_tol, in_relative_convergence_tol );
                            
                            
                            if( nodelb >= zuaux ) //we need to do this test because spatial B&B can have feasible solution in a node, but it is not the optimal solution...
                            {
                                prune = true;
                                
                                incOptPruneCounters( pruneCounter, pruneCounterLevel, currNode->getDepth() );
                            }
                        }
                        
                    }
                    else if( retCode == BBL_FEASIBLE_SOLUTION )
                    {
                        //we do nothing, but we put it here to avoid prune by infeasiblity in the next if
                    }
                    else if( retCode == BBL_INFEASIBLE_PROBLEM || in_consider_relax_infeas_if_solver_fail )
                    {
                        prune = true;
                        
                        incInfeasPruneCounters( pruneCounter, pruneCounterLevel, currNode->getDepth() );
                    
                    }
                    else
                    {
                        //we have an error at branchAndBound solving... we let user generates new nodes...
                    }
                    
                } //end of if( nodelb >= zuaux )
                
                //cout << "1.4 prune: " << prune << endl;
            }
            
        }
        
        //printf("Oi 4 Thread: %d\n", thnumber);
        //fflush(stdout);
        
        
        if( !prune )
        {
            //branching...
            
            
            if( branchStrat == BBL_BS_USER_NODE_GENERATION  )
            {
                userNodeGen->initialize(currNode, lx, ux, nlx, nux, nodelb, sol, dualSol);
                
                
                ret = userCallbacks.generateNodes( thnumber, *currNode, myiter, zl, zu, nlx, nux, retCode, objValue, sol, dualSol, *userNodeGen );
                
                if( ret != 0 )
                {
                    if( in_print_level >= BBL_LOW_PRINTING )
                        std::cerr << BBL_PREPRINT "Callback function userCalbacks.generateNodes returned error code " << ret << BBL_GETFILELINE << "\n";
                    
                    thReturnCodes[thnumber] = BBL_STOP_REQUIRED_BY_USER;
                    thReturnSubCodes[thnumber] = ret;
                    
                    goto termination;
                }
                
                
                
                const unsigned int nnodes = userNodeGen->getNumberOfNodes();
                
                
                
                if( nnodes > 0 )
                {
                    nodes->insertNodesDifsLowerBounds( thnumber, nnodes, userNodeGen->getNodesPointer()  );
                }
                else
                {
                    if( in_print_level > BBL_MEDIAN_PRINTING )
                    {
                        BBL_PRINTERRORMSG("Warning: user callback method generateNodes did not generate any node to perform branching on current iteration!\n\n");
                    }
                    
                }
                
            }
            else //in_branching_strategy == BBL_BS_USER_INDEX_CHOICE
            {
                BBL_Node *newNodes = NULL;
                
                unsigned int ninds = 0;
                
                
                ret = userCallbacks.chooseIndexToBranch( thnumber, *currNode, myiter, zl, zu, nlx, nux, retCode, objValue, dualObjValue, sol, dualSol, ninds, indices, values1, values2, newNodes);
                
                if(ret != 0)
                {
                    if( in_print_level >= BBL_LOW_PRINTING )
                        std::cerr << BBL_PREPRINT "Callback function userCalbacks.chooseIndexToBranch returned error code " << ret << BBL_GETFILELINE << endl;
                    
                    thReturnCodes[thnumber] = BBL_STOP_REQUIRED_BY_USER;
                    thReturnSubCodes[thnumber] = ret;
                    
                    goto termination;
                }
                
                
                
                if( ninds > 0 )
                {
                    ret = generateNodes( ninds, indices, values1, values2, currNode, nlx, nux, nodelb, sol, dualSol, in_store_parent_primal_solution_on_nodes, in_store_parent_dual_solution_on_nodes, newNodes );
                    
                    if(ret != 0)
                    {
                        if( in_print_level >= BBL_LOW_PRINTING )
                            std::cerr << BBL_PREPRINT "generateNodes returned error code " << ret << BBL_GETFILELINE << endl;
                        
                        thReturnCodes[thnumber] = ret;
                        
                        goto termination;
                    }
                    
                    nodes->insertNodes(thnumber, 1 << ninds, newNodes);
                }
                else
                {
                    if( in_print_level > BBL_MEDIAN_PRINTING )
                    {
                        BBL_PRINTERRORMSG("Warning: Callback function userCalbacks.chooseIndexToBranch did not choose any index to perform branch");
                    }
                }
                
            }
            
        }
        
        //printf("Oi 5 Thread: %d\n", thnumber);
        //fflush(stdout);
        
        if(in_call_end_of_iteration_callback)
        {
            ret = userCallbacks.endOfIteration( thnumber, myiter, (double(clock()-clockStart) )/CLOCKS_PER_SEC, BBL_getTime() - timeStart, zl, zu, *currNode, nlx, nux );
            
            if( ret != 0 )
            {
                if(in_print_level >= BBL_LOW_PRINTING)
                    std::cerr <<  BBL_PREPRINT "Callback function userCalbacks.endOfIteration returned error code " << ret << BBL_GETFILELINE << "\n";
                    
                thReturnCodes[thnumber] = BBL_STOP_REQUIRED_BY_USER;
                thReturnSubCodes[thnumber] = ret;
                
                goto termination;
            }
        }
        
        
        if( in_print_level > BBL_LOW_PRINTING && ( (myiter - 1) % in_printing_frequency == 0 ) )
        {
            char line[200];
            zlaux = zl;
            zuaux = zu;
            
            const double gap = zuaux - zlaux;
            long unsigned int nnodes = nodes->getNumberOfNodes();
            
            sprintf( line, BBL_PREPRINT "%-10ld  %+-14e  %+-14e  %+-14e  %-10ld  %u\n", myiter, zlaux, zuaux, gap, nnodes, thnumber );
            
            SEMAPH_Print.lock(nthreads);
            printf("%s", line );
            //cout << myiter << "  "<< zlaux << " " << zuaux << " " << gap << " " << nnodes << " " << thnumber << "\n" ;
            SEMAPH_Print.unlock(nthreads);
        }
        
        
        if( currNode->hasSharedData() )
        {
            //we just set the semaphore if the node has data being shared. In this way, we avoid semaphore usage and blocking other threads...
            
            SEMAPH_DelNodes.lock(nthreads);
                currNode->deallocateSharedData(); //we only call this methos in this zone. That is better than call the destructor directly
            SEMAPH_DelNodes.unlock(nthreads);
        }
        
        delete currNode;
        currNode = NULL;
        
        
        //if( checkTerminationCriterions(thnumber, zl, zu, myiter, (double(clock()-clockStart) )/CLOCKS_PER_SEC, BBL_getTime() - timeStart,  thReturnCodes[thnumber] )  )
        if( checkTerminationCriterions(thnumber, zl, zu, myiter, thReturnCodes[thnumber] )  )
        {
            goto termination;
        }
        
        //printf("Oi 6 Thread: %d\n", thnumber);
        //fflush(stdout);
        
        if( in_reorganize_lists && expStrategy == BBL_ES_BEST_LIMIT && threadIter % in_lists_reorganization_frequency == 0)
        {
            nodes->reorganizeNodeLists(thnumber);
        }
        
        
        
    }



    endLoop:	

    if( out_best_obj < in_infinity_value )
        thReturnCodes[thnumber] = BBL_OPTIMAL_SOLUTION;
    else
        thReturnCodes[thnumber] = BBL_NO_FEASIBLE_SOLUTION_FOUND;



    termination:


    if( in_call_after_bb_loop_callback )
    {
        userCallbacks.afterBBLoop(thnumber, zl, zu, thReturnCodes[thnumber]);
    }
        

    endThreads = true;
    thRunning[thnumber] = -1;


    if(thnumber == 0) 
    {
        //just to guarantee other threads are not loocked
        if( __lockAuxThreadsBeforeLoop )
            unlockAuxThreads(); 
    }


    if( in_count_total_prunes || in_max_tree_level_to_count_prunes_by_level > 0 )
    {
        SEMAPH_Attrib.lock(nthreads);
        {
            if( in_count_total_prunes )
                out_prune_counter.accumulate( pruneCounter);
            
            for(unsigned int i = 0; i < in_max_tree_level_to_count_prunes_by_level; i++)
                out_prune_counter_by_level[i].accumulate( pruneCounterLevel[i] );
        }
        SEMAPH_Attrib.unlock(nthreads);
    }

    if(nlx)		free(nlx);
    //if(nux)		free(nux);
    if(sol)		free(sol);
    //if(dualSol)	free(dualSol);

    if( userNodeGen )	delete userNodeGen;

    if(pruneCounterLevel)	delete[] pruneCounterLevel;

    if(indices)		free(indices);
    if(values1)		free(values1);
    if(values2)		free(values2);

    if(currNode)	delete currNode;


    //printf("Oi 7 Thread: %d\n", thnumber);
    //fflush(stdout);

    std::cout << BBL_PREPRINT "Thread " << thnumber << " finishing\n";

    return thReturnCodes[ thnumber ];
}






























