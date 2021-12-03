
#include <cstdlib>
#include <cstdio>

#include <cassert>
#include <math.h>
#include <cstring>



#include "BBL_branchAndBound.hpp"
#include "BBL_tools.hpp"




using namespace std;

using namespace branchAndBound;





BBL_PruneCounter::BBL_PruneCounter()
{
	reset();
}


void BBL_PruneCounter::reset()
{
	bound = infeas = opt = user = 0;
}


void BBL_PruneCounter::accumulate(BBL_PruneCounter& other)
{
	bound += other.bound;
	infeas += other.infeas;
	opt += other.opt;
	user += other.user;
}



BBL_Mutex::BBL_Mutex()
{
	initialize();
}


void BBL_Mutex::initialize()
{
	#if BBL_CPP_MULTITHREADING
		
	#else	
		#if BBL_OMP_MULTITHREADING
			omp_init_lock(&mutex);
		#endif
	#endif
}


int BBL_Mutex::lock( const unsigned int nthreads )
{
	
	#if BBL_CPP_MULTITHREADING
		if( nthreads > 1 )
			mymutex.lock();
	#else
		#if BBL_OMP_MULTITHREADING
			if( nthreads > 1 )
				omp_set_lock(&mutex);
		#endif
	#endif
	
	return 0;
}


int BBL_Mutex::tryLock( const unsigned int nthreads )
{
	#if BBL_CPP_MULTITHREADING
		if( nthreads > 1 )
			return mymutex.try_lock() == true ? 0 : BBL_UNDEFINED_ERROR;
	#else
		#if BBL_OMP_MULTITHREADING
			if( nthreads > 1 )
				return omp_test_lock(&mutex) == 1 ? 0 : BBL_UNDEFINED_ERROR;
		#endif
	#endif
	
	return 0;
}


int BBL_Mutex::unlock( const unsigned int nthreads )
{
	#if BBL_CPP_MULTITHREADING
		if( nthreads > 1 )
			mymutex.unlock();
	#else
		#if BBL_OMP_MULTITHREADING
			if( nthreads > 1 )
				omp_unset_lock(&mutex);
		#endif
	#endif
	
	return 0;
}


void BBL_Mutex::destroy()
{
	#if BBL_OMP_MULTITHREADING
		omp_destroy_lock(&mutex);
	#endif
}


BBL_Mutex::~BBL_Mutex()
{
	destroy();
}





BBL_NodeList:: BBL_NodeList(const int strategy)
{
    initialize(strategy);
}


void BBL_NodeList::initialize(const int strategy)
{
	head = tail = iter = NULL;
	firstlb = INFINITY;
    nNodes = 0;
    expStrat = strategy;
}



//if maxLevelPruneCounter > 0, that method accumulates in pruneLevelCounter the number of prunes by level
long unsigned int BBL_NodeList::pruneNodesByBound(const double zu, double &zl, const unsigned int maxLevelPruneCounter, unsigned int *pruneLevelCounter)
{
	long unsigned int nodesPruned = 0;
    bool chgIter = false;
    BBL_Node *p;
	
	zl = -INFINITY;
	
    if(head == NULL)
		return 0;
    
	
    if(expStrat == BBL_ES_BEST_LIMIT)
    {
		
		if(iter->getLowerBound() >= zu)
		    chgIter = true; //iter will be prune...
		
		//we run starting from the end.
		for(  ; (tail != NULL) && tail->getLowerBound() >= zu ; tail = p )
		{
		    p = tail->previous;
			
			if( tail->getDepth() < maxLevelPruneCounter )
				pruneLevelCounter[ tail->getDepth() ]++;
			
		    delete tail;
		    nodesPruned++; //nNodes--;
		}
		
		if(tail == NULL)
		{//we pruned all nodes
		    head = iter = NULL;
			firstlb = INFINITY;
		}
		else
		{
		    //The nodes are ordered by lower bound
		    zl = head->getLowerBound();
		    
		    tail->next = NULL;
		    if(chgIter)
				iter = tail;
		}
    }
    else
    {
		if(head == NULL)
		    return 0;
		
		zl = BBL_min(head->getLowerBound(), tail->getLowerBound());
	
	
		if(head == tail)
		{
		    if(head->getLowerBound() >= zu)
		    {
				delete head;
				nodesPruned = 1; //nNodes--;
				nNodes -= 1;
				head = tail = iter = NULL;
				firstlb = INFINITY;
		    }
		    
		    return nodesPruned;
		}
	
	
		for(iter = head->next; iter != tail; iter = p)
		{
		    p = iter->next;
			
			if( iter->getLowerBound() < zl )
				zl = iter->getLowerBound();
		    
		    if(iter->getLowerBound() >= zu)
		    {
				//head and tail are not considered. So, we have a previus and a next
				iter->previous->next = iter->next;
				iter->next->previous = iter->previous;
				
				if( iter->getDepth() < maxLevelPruneCounter)
					pruneLevelCounter[iter->getDepth()]++;
				
				delete iter;
				nodesPruned++; //nNodes--;
		    }
		    
		}
	
		//now, we check the head and the tail
		if( head->getLowerBound() >= zu )
		{
		    p = head;
		    head = head->next; //the next exist because tail is diferent to head
		    firstlb = head->getLowerBound();
		    
		    if( p->getDepth() < maxLevelPruneCounter )
				pruneLevelCounter[p->getDepth()]++;
			
		    delete p;
		    nodesPruned++; //nNodes--;
		    head->previous = NULL;
		}
	
		if( tail->getLowerBound() >= zu )
		{
		    p = tail;
		    
		    if(head == tail)
		    {
				head = tail = iter = NULL;
				firstlb = INFINITY;
		    }
		    else
		    {
				//if(tail->previous)
				    //tail->previous = NULL;
				
				if(tail == iter)
				    iter = tail->previous;
				
				tail = tail->previous;
				tail->next = NULL;
		    }
		    
		    if( p->getDepth() < maxLevelPruneCounter )
				pruneLevelCounter[p->getDepth()]++;
		    
		    delete p;
		    nodesPruned++; //nNodes--;
		}
		
    }
    
    nNodes -= nodesPruned;
    return nodesPruned;
}





long unsigned int BBL_NodeList::getNodePointer(const double zu, BBL_Node* &node, const unsigned int maxLevelPruneCounter, unsigned int* pruneLevelCounter)
{
    //we return the first node better than zu...
    long unsigned int nodesPruned = 0;
    //bool fiter = false;
    
    for( ; ((head != NULL) && head->getLowerBound() >= zu) ; head = node )
    {
		//if(head == iter)
		    //fiter = true;
		
		node = head->next;
		//if(*node) //we do not need do it because the head will be removed...
			//(*node)->previous = NULL;
		
		if( head->getDepth() < maxLevelPruneCounter )
			pruneLevelCounter[head->getDepth()]++;
		
		delete head;
		nodesPruned++; //nNodes--;
    }
    
    
    node = head; //that node is NULL or is lower than zu
    if(head == NULL)
    {
		tail = iter = NULL;
		firstlb = INFINITY; //head could have been deleted
    }
    else
    {
		head = head->next;
		if(head)
		{
		    head->previous = NULL;
			firstlb = head->getLowerBound();
		}
		else
		    tail = NULL;
		
		/*if(fiter) //the node pointed by iter was deleted
		    iter = head;
		else  */
			if( node == iter) 
			    iter = head; //iter is only used in best limit strategy. So, we only need check if iter points to old head, i.e., *node. Note that if head->getLowerBound() is >= zu, so all nodes will be prune in a best limit strategy...
		
		node->previous = NULL;
		node->next = NULL;
		
		nNodes--;
    }
    
    nNodes -= nodesPruned;
    return nodesPruned;
}




void BBL_NodeList::reorganizeToBestLimit(void )
{
	#if BBL_DEBUG_MODE
		unsigned int onNodes = nNodes;
	#endif
	
	unsigned long int n;
	double lb;
	BBL_Node *oiter, *nodes;
	
	
	if( expStrat == BBL_ES_BEST_LIMIT )
		return;
	
	//backing-up the list
	oiter = head;
	//otail = tail;
	//oiter = iter;
	
	//reseting the list
	initialize(BBL_ES_BEST_LIMIT);
	
	//now, we reinsert the nodes. The nodes will be inserted in the new strategy. We only can insert a list of nodes having the same lower bound...
	
	
	while( oiter )
	{
		nodes = oiter;
		lb = oiter->getLowerBound();
		n = 1;
		
		for( oiter = oiter->next; oiter ; oiter = oiter->next )
		{
			if( oiter->getLowerBound() != lb )
			{
				oiter->previous->next = NULL;
				oiter->previous = NULL;
				break;
			}
			
			n++;
		}
		
		insertNodes(n, nodes);
	}
	
	#if BBL_DEBUG_MODE
		assert( nNodes == onNodes );
	#endif
}



/*
void BBL_NodeList::reorganizeTofDad(void )
{
	#if BBL_DEBUG_MODE
		unsigned int onNodes = nNodes;
	#endif
	
	unsigned long int n;
	double fDad;
	BBL_Node *oiter, *nodes;
	
	
	if( expStrat == BBL_ES_BEST_DAD_RELAX )
		return ;
	
	//backing up the list
	oiter = head;
	
	//reseting the list
	initialize(BBL_ES_BEST_DAD_RELAX);
	
	
	//now, we reinsert the nodes. The nodes will be inserted in the new strategy. We only can insert a list of nodes having the same lower bound...
	
	while( oiter )
	{
		nodes = oiter;
		fDad = oiter->fDad;
		n = 1;
		
		for( oiter = oiter->next; oiter ; oiter = oiter->next )
		{
			if( oiter->fDad != fDad )
			{
				oiter->previous->next = NULL;
				oiter->previous = NULL;
				break;
			}
			
			n++;
		}
		
		insertNodes(n, nodes);
	}
	
	#if BBL_DEBUG_MODE
		assert( nNodes == onNodes );
	#endif
}

*/





//we assume that the nodes have the same parent 
void BBL_NodeList::insertNodes(const unsigned int numberOfNodes, BBL_Node* nodes)
{
	#if BBL_DEBUG_MODE
	    unsigned int n = numberOfNodes;
	#endif
		
    BBL_Node *endn;
    
    nNodes += numberOfNodes;
    
    for(endn = nodes; endn->next != NULL; endn = endn->next)
    {
		#if BBL_DEBUG_MODE
		    n--;
		#endif
    }
    
    #if BBL_DEBUG_MODE
		//printf("n: %d\n", n);
		assert(n == 1);
    #endif
    
    /*if(expStrat == BBL_ES_BEST_DAD_RELAX)
	{
		if(head == NULL)
		{
			head = nodes;
			tail = endn;
			iter = head;
		}
		else
		{
			if( nodes->fDad >= tail->fDad ) // we are sure tail is not null
			{
				//we must insert after the tail
				tail->next = nodes;
				nodes->previous = tail;
				tail = endn;
				return;
			}
			
			if( nodes->fDad <= head->fDad )
			{//we must insert before the head
				endn->next = head;
				head->previous = endn;
				head = nodes;
				return;
			}
			
			if( nodes->fDad > iter->fDad )
			{
				
				if( nodes->fDad <= ( iter->fDad + tail->fDad )/2.0 )
				{
					//we start from iter. We are sure that nodes must be before tail
					for( iter = iter->next ; iter->fDad < nodes->fDad ; iter = iter->next );
					
					//nodes must be before iter
					
					iter->previous->next = nodes;
					nodes->previous = iter->previous;
					
					iter->previous = endn;
					endn->next = iter;
				}
				else
				{
					//we start from tail. We are sure that nodes must be after (current) iter
					for(iter = tail->previous; iter->fDad > nodes->fDad; iter = iter->previous);
					
					//nodes must be after iter
					iter->next->previous = endn;
					endn->next = iter->next;
					
					iter->next = nodes;
					nodes->previous = iter;
				}
				
				return;
			}
			else //nodes->fDad <= iter->fDad
			{
				
				if( nodes->fDad > (head->fDad + iter->fDad)/2.0 )
				{
					//we start from iter. We are sure nodes must be after head.
					for( iter = iter->previous; iter->fDad > nodes->fDad; iter = iter->previous );
					
					//nodes must be after iter
					iter->next->previous = endn;
					endn->next = iter->next;
					
					iter->next = nodes;
					nodes->previous = iter;
				}
				else //nodes->fDad <= (head->fDad + iter->fDad)/2.0
				{
					//we start from head. We are shure nodes must be before iter.
					for( iter = head->next; iter->fDad < nodes->fDad; iter = iter->next);
					
					//nodes must be before iter
					iter->previous->next = nodes;
					nodes->previous = iter->previous;
					
					iter->previous = endn;
					endn->next = iter;
				}
				
				return;
			}
			
		}
		
		
	}
    else */
	if(expStrat == BBL_ES_BEST_LIMIT)
    {
		if(head == NULL)
		{
		    head = nodes;
			firstlb = head->getLowerBound();
		    tail =  endn;
			iter = head;
			
		}
		else
		{
		    
		    if( nodes->getLowerBound() >= tail->getLowerBound() ) // we are sure tail is not null
		    {
				//we must insert after the tail
				tail->next = nodes;
				nodes->previous = tail;
				tail = endn;
				return;
		    }
		    
		    if( nodes->getLowerBound() <= head->getLowerBound() )
		    {//we must insert before the head
				endn->next = head;
				head->previous = endn;
				head = nodes;
				firstlb = head->getLowerBound();
				return;
		    }
		    
		    if( nodes->getLowerBound() > iter->getLowerBound()  )
		    {
				
				if( nodes->getLowerBound() <= ( iter->getLowerBound() + tail->getLowerBound() )/2.0 )
				{
					//we start from iter. We are sure that nodes must be before tail
					for( iter = iter->next ; iter->getLowerBound() < nodes->getLowerBound() ; iter = iter->next );
					
					//nodes must be before iter
					iter->previous->next = nodes;
					nodes->previous = iter->previous;
					
					iter->previous = endn;
					endn->next = iter;
				}
				else
				{
					//we start from tail. We are sure that nodes must be after (current) iter
					for(iter = tail->previous; iter->getLowerBound() > nodes->getLowerBound(); iter = iter->previous );
					
					//nodes must be after iter
					iter->next->previous = endn;
					endn->next = iter->next;
					
					iter->next = nodes;
					nodes->previous = iter;
				}
			
				return;
		    }
		    else //nodes->getLowerBound() <= iter->getLowerBound()
		    {
				
				if( nodes->getLowerBound() > ( head->getLowerBound() + iter->getLowerBound() )/2.0 )
				{
					//we start from iter. We are sure nodes must be after head.
					for( iter = iter->previous; iter->getLowerBound() > nodes->getLowerBound(); iter = iter->previous );
					
					//nodes must be after iter
					iter->next->previous = endn;
					endn->next = iter->next;
					
					iter->next = nodes;
					nodes->previous = iter;
				}
				else // nodes->getLowerBound() <= ( head->getLowerBound() + iter->getLowerBound() )/2.0
				{
					//we start from head. We are shure nodes must be before iter.
					for( iter = head->next; iter->getLowerBound() < nodes->getLowerBound(); iter = iter->next );
					
					//nodes must be before iter
					iter->previous->next = nodes;
					nodes->previous = iter->previous;
					
					iter->previous = endn;
					endn->next = iter;
				}
				
				return;
		    }
		     
		}
		
    }
    else if(expStrat == BBL_ES_DEPTH)
    {
		//we should insert in the begin of encadeate list...
		//nodes->previous = NULL;
		
		//p->next = head;
		
		if(head == NULL)
		    tail = endn;
		else
		{
		    head->previous = endn;
		    endn->next = head;
		}
		
		head = nodes;
		firstlb = head->getLowerBound();
    }
    else
    {
		//width exploration... we should insert the node in the tail...
		if(tail == NULL)
		{
		    head = nodes;
			firstlb = head->getLowerBound();
		}
		else
		{
		    tail->next = nodes;
		    nodes->previous = tail;
		}
		tail = endn;
    }
    
    return;
}




unsigned long int BBL_NodeList::deallocateAllNodes(void)
{
    unsigned long int n = 0;
    BBL_Node *p;
    
    for(iter = head; iter != NULL; iter = p, n++ )
    {
		p = iter->next;
		delete iter;
    }
    
    #if BBL_DEBUG_MODE
	    //printf("n: %lu nNodes: %lu\n", n, nNodes);
		assert( n == nNodes );
    #endif
    
    nNodes = 0;
    head = NULL;
	firstlb = INFINITY;
	tail = NULL;
    
    return n;
}


void BBL_NodeList::print(void ) const
{
    int i = 0;
    BBL_Node *p;

    printf("\t********************************************************************************\n");
    printf("\tnNodes: %lu expStrat: %d\n", nNodes, expStrat);
    for(p = head; p != NULL; p = p->next, i++ )
    {
		printf("\tnode: %d depth: %u lb: %f\n", i, p->getDepth(), p->getLowerBound());
		//printf("\t-----------------------------------------------------\n");
		p->print();
		//for(int j = 0; j < (int) p->nBounds; j++)
			//printf("\tvar: %u l: %f u: %f \n", p->bounds[j].ind, p->bounds[j].l, p->bounds[j].u);
		
		printf("\t.....................................................\n");
    }
    printf("\t********************************************************************************\n");
    
}


unsigned long int BBL_NodeList::countNodes(void ) const
{
    unsigned long int i = 0;
    BBL_Node *p;
	
    for(p = head; p; p = p->next, i++);
    
    return i;
}


int BBL_NodeList::getFirstNodeLowerBound(double& lb) const
{
	if(head)
	{
		//lb = head->getLowerBound();
		lb = firstlb;
		return 0;
	}
	else
	{
		return -1;
	}
}


/*int BBL_NodeList::getFirstNodefDad(double& fdad)
{
	if(head)
	{
		fdad = head->fDad;
		return 0;
	}
	else
	{
		return -1;
	}
} */



int BBL_NodeList::getFirstNodeDepth(const double zu, unsigned int& depth, long unsigned int &nNodesPruned)
{
	BBL_Node *nodeAux;
	
	nNodesPruned = 0;
	
	
	for( ; head; head = nodeAux)
	{
		
		nodeAux = head->next;
		
		if( head->getLowerBound() >= zu )
		{
			if(nodeAux)
			{
				nodeAux->previous = NULL;
				firstlb = nodeAux->getLowerBound();
			}
			else
			{
				firstlb = INFINITY;
			}
			
			delete head;
			nNodesPruned++;
			
		}
		else
		{
			break;
		}
		
	}
	
	nNodes -= nNodesPruned;
	
	if(head)
	{
		depth = head->getDepth();
		return 0;
	}
	else
	{
		tail = NULL;
		return -1;
	}
}


long unsigned int BBL_NodeList::getFirstNodes(const long unsigned int numberOfNodes, BBL_Node* &nodes)
{
	long unsigned int i, aux;
	
	nodes = head;
	
	if( numberOfNodes >= nNodes )
	{
		head = tail = iter = NULL;
		firstlb = INFINITY;
		
		i = nNodes;
		nNodes = 0;
		
		return i;
	}
	else
	{
		if( nNodes >= 2*numberOfNodes )
		{
			if(numberOfNodes == 0)
			{
				nodes = NULL;
				return 0;
			}
			
			//we start from the head
			
			aux = numberOfNodes - 1;
			for(i = 0, iter = head; i < aux; i++, iter = iter->next);
			
			head = iter = iter->next;
			firstlb = head->getLowerBound();
		}
		else
		{
			//we start from the tail
			aux = nNodes - numberOfNodes - 1; //we are sure numberOfNodes < nNodes
			
			for(i = 0, iter = tail;  i < aux; i++, iter = iter->previous);
			
			head = iter;
			firstlb = head->getLowerBound();
		}
		
		head->previous->next = NULL;
		head->previous = NULL;
		
		
		nNodes -= numberOfNodes;
		return numberOfNodes;
	}
	
}

long unsigned int BBL_NodeList::getLastNodes(const long unsigned int numberOfNodes, BBL_Node* &nodes)
{
	long unsigned int i, aux;
	
	if( numberOfNodes >= nNodes)
	{
		nodes = head;
		
		head = tail = iter = NULL;
		firstlb = INFINITY;
		
		i = nNodes;
		nNodes = 0;
		
		return i;
	}
	else
	{
		if( nNodes <= 2*numberOfNodes )
		{ //if nNodes and numberOfNodes are both 0, we had enter in the if above
			//we start from the head;
			
			aux = nNodes - numberOfNodes - 1;
			for(i = 0, iter = head; i < aux; i++, iter = iter->next);
			
			tail = iter;
		}
		else
		{
			if(numberOfNodes == 0)
			{
				nodes = NULL;
				return 0;
			}
			
			//we start from the tail
			aux = numberOfNodes - 1;
			for(i = 0, iter = tail; i < aux; i++, iter = iter->previous);
			
			tail = iter = iter->previous;
		}
		
		nodes = tail->next;
		(nodes)->previous = NULL;
		
		tail->next = NULL;
		
		nNodes -= numberOfNodes;
		return numberOfNodes;
	}
	
}


/*long unsigned int BBL_NodeList::getNumberOfNodes()
{
	return nNodes;
}*/




BBL_NodeList::~BBL_NodeList()
{
    deallocateAllNodes();
}




BBL_NodeListManager::BBL_NodeListManager()
{
	initialize();
}

void BBL_NodeListManager::initialize()
{
	expStrat = BBL_ES_BEST_LIMIT;
	nNodeLists = 1000;
	nodeLists = NULL;
	llb = NULL;
	pruneLevelCounter = NULL;
	nNodes = 0;
	reorgFactorToLastEmptyList = 0.05;
}



void BBL_NodeListManager::calculateLBintervalsForLists(const double begLBInterval, const double endLBInterval)
{
	unsigned int i;
	double step;
	
	llb[0] = -INFINITY;
	
	if( nNodeLists > 1 )
	{
		if( isinf(endLBInterval) )
		{
			for(i = 1; i < nNodeLists; i++)
				llb[i] = INFINITY; //all nodes will be in the first node list...
		}
		else
		{
			step = (endLBInterval - begLBInterval)/(nNodeLists - 1);
			
			for(i = 1; i < nNodeLists; i++)
				llb[i] = begLBInterval + step*i;
		}
	}
}



void BBL_NodeListManager::acumulatePruneLevelCounters(const unsigned int maxLevelPruneCounter, BBL_PruneCounter* counter) const
{
	unsigned int i;
	
	for(i = 0; i < maxLevelPruneCounter; i++)
		counter[i].bound += pruneLevelCounter[i];
}



int BBL_NodeListManager::allocateLists(const int strategy, const unsigned int nNodeLists, const double begLBInterval, const double endLBInterval, const unsigned int maxLevelPruneCounter)
{
	unsigned int i;
	int code;
	
	expStrat = strategy;
	
	/*if( strategy == BBL_ES_BEST_LIMIT )
		this->nNodeLists = nNodeLists;
	else
		this->nNodeLists = 1; */
	
	this->nNodeLists = nNodeLists; //we let allocate all list requested although depth and width only use one list. We do it because user can change to best limi strategy after a feasible solution, for example...
	
	if(maxLevelPruneCounter > 0)
	{
		pruneLevelCounter = (unsigned int *) calloc( maxLevelPruneCounter, sizeof(unsigned int) );
		if(!pruneLevelCounter)
		{
			code = BBL_MEMORY_ERROR;
			goto desallocate_memory;
		}
	}
	
	
	nodeLists = new(nothrow) BBL_NodeList[nNodeLists];
	llb = (double *) malloc( nNodeLists * sizeof(double) );
	
	if(!nodeLists || !llb)
	{
		code = BBL_MEMORY_ERROR;
		goto desallocate_memory;
	}
	
	for(i = 0; i < nNodeLists; i++)
		nodeLists[i].initialize(expStrat);
	
	
	calculateLBintervalsForLists(begLBInterval, endLBInterval);
	
	code = 0;
	
desallocate_memory:
	
	return code;
	
}



void BBL_NodeListManager::deallocateLists()
{
	if(nodeLists)
	{
		delete[] nodeLists;
		nodeLists = NULL;
	}
	
	if(llb)
	{
		free(llb);
		llb = NULL;
	}
	
	if(pruneLevelCounter)
	{
		free(pruneLevelCounter);
		pruneLevelCounter = NULL;
	}
}


long unsigned int BBL_NodeListManager::pruneNodesByBound(const double zu, double& zl, const unsigned int maxLevelPruneCounter)
{
	int i; //do not use unsigend here..
	long unsigned int nodesPruned = 0;
	double zl_lists = INFINITY, aux;
	
	
	for(i = nNodeLists - 1; i >= 0; i--)
	{
		
		if( nodeLists[i].getNumberOfNodes() > 0 )
		{
			
			nodesPruned +=  nodeLists[i].pruneNodesByBound(zu, aux, maxLevelPruneCounter, pruneLevelCounter);
			
			//printf("Podei a lista %u. Nos podados: %d\n", i, auxInt);
			
			if( !isinf(aux) && aux < zl_lists )
				zl_lists = aux;
		}
		
		if(llb[i] < zu )
		{
			break; //list before that do not have nodes to prune.
		}
	}
	
	if( isinf(zl_lists) )
		zl = -INFINITY; //zl_list has the value INFINITY
	else
		zl = zl_lists;
	
	
	nNodes -= nodesPruned;
	
	return nodesPruned;
}



long unsigned int BBL_NodeListManager::getNodePointer(const double zu, BBL_Node* &p, const unsigned int maxLevelPruneCounter)
{
	long unsigned int n = 0;
	p = NULL;
	
	for(decltype(nNodeLists) i = 0; i < nNodeLists; i++)
	{
		if( nodeLists[i].getNumberOfNodes() > 0 )
		{
			n += nodeLists[i].getNodePointer(zu, p, maxLevelPruneCounter, pruneLevelCounter);
			
			if(p)
			{
				nNodes--;
				break;
			}
		}
	}
	
	nNodes -= n;
	return n;
}



void BBL_NodeListManager::insertNodes(const unsigned int numberOfNodes, BBL_Node* nodes)
{
	unsigned int i, list = 0, beg, end;
	const double lb = nodes->getLowerBound();
	
	nNodes += numberOfNodes;
	
	if( expStrat == BBL_ES_BEST_LIMIT )
	{
		// we need find the correct list to put the nodes
		if( lb >= llb[nNodeLists - 1] ) //we put <= instead of < because we can get a strange case where all sublists has the same llb
		{
			list = nNodeLists - 1; //we must insert in the last one
		}
		else
		{
			beg = 0;
			end = nNodeLists - 2; //if we have only one list, we will enter in the if above
			
			while( true )
			{
				i = (beg + end)/2;
				
				if( lb < llb[i] )
				{
					end = i - 1;
				}
				else 
				{
					if( lb < llb[i + 1] ) //lb >= llb[i]
					{
						list = i; //we find the place
						break;
					}
					else
					{
						beg = i + 1;
					}
				}
			}
		}
	}
	
	nodeLists[list].insertNodes(numberOfNodes, nodes);
}




unsigned long int BBL_NodeListManager::deallocateAllNodes(void )
{
	unsigned int i;
	unsigned long int n = 0;
	
	for(i = 0; i < nNodeLists; i++)
	{
		n += nodeLists[i].deallocateAllNodes();
	}
	
	#if BBL_DEBUG_MODE
		//printf("n: %d nNodes: %d\n", n, nNodes);
		assert( n == nNodes );
    #endif
	
	nNodes = 0;
	
	return n;
}


void BBL_NodeListManager::print(void ) const
{
	unsigned int i;
	
	printf("total nNodes: %lu\n", nNodes);
	
	for(i = 0; i < nNodeLists; i++)
	{
		if( nodeLists[i].getNumberOfNodes() > 0 )
		{
			printf("\n\nnode list %u llb[%u]: %f\n", i, i, llb[i]);
			nodeLists[i].print();
		}
	}
}


unsigned long int BBL_NodeListManager::countNodes(void ) const
{
	long unsigned int nodes = 0;
	unsigned int i;
	
	for(i = 0; i < nNodeLists; i++)
	{
		nodes += nodeLists[i].countNodes();
	}
	
	return nodes;
}



void BBL_NodeListManager::reorganizeToBestLimit(const double begLBInterval, const double endLBInterval)
{
	unsigned int i;
	
	if( expStrat == BBL_ES_BEST_LIMIT )
		return;
	
	//in either depth or width strategies, we only use the first node list. Anyway, we need change the expStrat in other lists
	
	calculateLBintervalsForLists(begLBInterval, endLBInterval);
	
	for(i = 0; i < nNodeLists; i++)
		nodeLists[i].reorganizeToBestLimit();
	
	
	expStrat = BBL_ES_BEST_LIMIT;
}





//that method balance the nodes in the lists. All list will have the same quantity of nodes in general (except maybe the last one). It only makes sense to use with BEST_LIMIT strategy.
void BBL_NodeListManager::reorganizeNodeLists()
{
	unsigned int i, j;
	unsigned int n, aux, aux2;
	const long unsigned int nNodesAtEachSubList = ceil( (double)nNodes/nNodeLists );
	BBL_Node *nodes;
	const double lbToUseIfPreviousIsZero = 0.01; //I do not know what to do about this value...
	
	
	if(nNodes == 0 )//if( expStrat != BBL_ES_BEST_LIMIT || nNodes == 0 )
		return;
	
	
	for(i = 0; i < nNodeLists - 1; i++) //the last list will be unbalanced, but ok
	{
		n = nodeLists[i].getNumberOfNodes();
		
		if( i*nNodesAtEachSubList < nNodes )
		{
			if( n > nNodesAtEachSubList )
			{
				aux = n - nNodesAtEachSubList;
				const auto nnodesret = nodeLists[i].getLastNodes(aux, nodes);
				
				#if BBL_DEBUG_MODE
					assert( nnodesret == aux );
				#endif
				
				nodeLists[i+1].insertNodes(aux, nodes); //that insertion will be fast because the nodes will be before the head in that new list
			}
			else if( n < nNodesAtEachSubList )
			{
				aux = nNodesAtEachSubList - n;
				
				for(j = i + 1; j < nNodeLists; j++)
				{
					aux2 = nodeLists[j].getNumberOfNodes();
					if( aux2 > 0 )
					{
						if( aux2 >= aux )
						{
							const unsigned int nnodesret = nodeLists[j].getFirstNodes(aux, nodes);
							
							#if BBL_DEBUG_MODE
								assert( nnodesret == aux );
							#endif
							
							nodeLists[i].insertNodes(aux, nodes); //that insertion will be fast because the nodes will be after the tail in that new list.
							break;
						}
						else
						{
							const unsigned int nnodesret = nodeLists[j].getFirstNodes(aux2, nodes);
							
							#if BBL_DEBUG_MODE
								assert( nnodesret == aux2 );
							#endif
							
							nodeLists[i].insertNodes(aux2, nodes); //that insertion will be fast because the nodes will be after the tail in that new list.
							aux = aux - aux2;
						}
					}
				}
			}
			
			//updating llb[i]
			nodeLists[i].getFirstNodeLowerBound( llb[i] ); 
		}
		else
		{
			//we have already allocated all nodes. It will happen if the number of lists is greater than the number of nodes.
			
			if( llb[i - 1] > 0 )
				llb[i] = (1 + reorgFactorToLastEmptyList) * llb[i - 1];
			else if( llb[i - 1] < 0)
				llb[i] = (1 - reorgFactorToLastEmptyList) * llb[i - 1];
			else // llb[i - 1] == 0.0 ???? 
				llb[i] = lbToUseIfPreviousIsZero; //I do not konw waht to do...
		}
		
	}
	
	//updating the llb of the last list
	nodeLists[nNodeLists - 1].getFirstNodeLowerBound( llb[nNodeLists - 1] ); 
	
	llb[0] = -INFINITY; //to avoid problems, we mantain the first llb at -INFINITY
}


/*long unsigned int BBL_NodeListManager::getNumberOfNodes() const
{
	return nNodes;
}*/



long unsigned int BBL_NodeListManager::getFirstNodes(long unsigned int numberOfNodes, BBL_Node* &nodes)
{
	long unsigned int mynnodes = numberOfNodes, nNodesGotten;
	BBL_Node *myNodes, *pLastNode;
	
	nodes = NULL;
	
	for(decltype(nNodeLists) i = 0; i < nNodeLists; i++)
	{
		const auto nNodesAtList = nodeLists[i].getNumberOfNodes();
		
		if( nNodesAtList == 0 )
			continue;
		
		decltype(nNodesAtList) auxNodesToGet = BBL_min(mynnodes, nNodesAtList);
		
		
		const auto nodesGotten = nodeLists[i].getFirstNodes(auxNodesToGet, myNodes);
		
		#if BBL_DEBUG_MODE
			assert( nodesGotten == auxNodesToGet );
		#endif
		
		if( nodes == NULL )
		{
			nodes = pLastNode = myNodes;
		}
		else
		{
			#if BBL_DEBUG_MODE
				assert(pLastNode->next == NULL);
				assert(myNodes->previous == NULL);
			#endif
			
			pLastNode->next = myNodes;
			myNodes->previous = pLastNode;
		}
		
		mynnodes -= auxNodesToGet;
		
		if( mynnodes <= 0 ) //we aready got the requested nodes. we do not need updhate pLastNode
			break;
		
		//updating pLastNode
		for( ;   pLastNode->next != NULL   ;   pLastNode = pLastNode->next );
	}
	
	#if BBL_DEBUG_MODE
		assert( numberOfNodes >= mynnodes );
	#endif
	
	nNodesGotten = numberOfNodes - mynnodes;
	
	#if BBL_DEBUG_MODE
		assert( this->nNodes >= nNodesGotten );
	#endif
	
	this->nNodes -= nNodesGotten;
	
	return nNodesGotten;
}




long unsigned int BBL_NodeListManager::getLastNodes(long unsigned int numberOfNodes, BBL_Node* &nodes)
{
	long unsigned int mynnodes = numberOfNodes, nNodesGotten;
	BBL_Node *myNodes;
	
	nodes = NULL;
	
	for(decltype(nNodeLists) i = nNodeLists; i > 0; i--)
	{
		decltype(i) indList = i-1;
		const auto nNodesAtList = nodeLists[indList].getNumberOfNodes();
		
		if( nNodesAtList == 0 )
			continue;
		
		decltype(nNodesAtList) auxNodesToGet = BBL_min(mynnodes, nNodesAtList);
		
		const auto nodesGotten = nodeLists[indList].getLastNodes(auxNodesToGet, myNodes);
		
		#if BBL_DEBUG_MODE
			assert( nodesGotten == auxNodesToGet );
		#endif
		
		BBL_appendToEncadeateList(myNodes, nodes);
		nodes = myNodes;
		
		mynnodes -= auxNodesToGet;
		if( mynnodes <= 0 )
			break;
	}
	
	#if BBL_DEBUG_MODE
		assert( numberOfNodes >= mynnodes );
	#endif
	
	nNodesGotten = numberOfNodes - mynnodes;
	
	#if BBL_DEBUG_MODE
		assert( this->nNodes >= nNodesGotten );
	#endif
	
	this->nNodes -= nNodesGotten;
	
	return nNodesGotten;
}



int BBL_NodeListManager::getFirstNodeDepth(const double zu, unsigned int& depth, long unsigned int &nNodesPruned)
{
	unsigned int i;
	int aux;
	long unsigned int nNodesPrunedList;
	
	nNodesPruned = 0;
	
	for(i = 0; i < nNodeLists; i++)
	{
		aux = nodeLists[i].getFirstNodeDepth(zu, depth, nNodesPrunedList);
		
		nNodesPruned += nNodesPrunedList;
		
		if(aux == 0)
			break; //we already have the depth of first node
	}
	
	nNodes -= nNodesPruned;
	return aux; //we return the return value of the last list...
}


int BBL_NodeListManager::getFirstNodeLowerBound(double& lb) const
{
	unsigned int i;
	int aux;
	
	for(i = 0; i < nNodeLists; i++)
	{
		aux = nodeLists[i].getFirstNodeLowerBound(lb);
		if(aux == 0)
			return aux;
	}
	
	return aux; //we return the return value of the last list...
}




BBL_NodeListManager::~BBL_NodeListManager()
{
	   deallocateLists();
}


BBL_MTNodeListManager::BBL_MTNodeListManager()
{
	initialize();
}


void BBL_MTNodeListManager::initialize()
{
	expStrat = BBL_ES_BEST_LIMIT;
	nThreads = 1;
	nNodeLists = 1;
	nodeLists = NULL;
	SEMAPH_lists = SEMAPH_listQueries = NULL;
}


void BBL_MTNodeListManager::acumulatePruneLevelCounters(const unsigned int maxLevelPruneCounter, BBL_PruneCounter* counter) const
{
	unsigned int i;
	
	for(i = 0; i < nThreads; i++)
		nodeLists[i].acumulatePruneLevelCounters(maxLevelPruneCounter, counter);
}




int BBL_MTNodeListManager::allocateLists(const int strategy, const unsigned int numberOfThreads, const unsigned int nNodeSubLists, const double begLBInterval, const double endLBInterval, const unsigned int maxLevelPruneCounter)
{
	int code, aux;
	unsigned int i;
	expStrat = strategy;
	nThreads = nNodeLists = numberOfThreads;
	
	
	
	nodeLists = new (nothrow) BBL_NodeListManager[nNodeLists];
	SEMAPH_lists = new (nothrow) BBL_Mutex[nNodeLists];
	SEMAPH_listQueries = new (nothrow) BBL_Mutex[nNodeLists];
	if( !nodeLists || !SEMAPH_lists || !SEMAPH_listQueries )
	{
		code = BBL_MEMORY_ERROR;
		goto desallocate_memory;
	}
		
	
	
	for(i = 0; i < nNodeLists; i++)
	{
		aux = nodeLists[i].allocateLists(expStrat, nNodeSubLists, begLBInterval, endLBInterval, maxLevelPruneCounter);
		if(aux != 0)
		{
			code = aux;
			goto desallocate_memory;
		}
	}
	
	
	code = 0;
	
desallocate_memory:

	return code;
}


long unsigned int BBL_MTNodeListManager::countNodes(void ) const
{
	long unsigned int n = 0;
	unsigned int i;
	
	for(i = 0; i < nNodeLists; i++)
		n += nodeLists[i].countNodes();
	
	return n;
}


long unsigned int BBL_MTNodeListManager::desallocateAllNodes(void )
{
	long unsigned int n = 0;
	unsigned int i;
	
	for(i = 0; i < nNodeLists; i++)
		n += nodeLists[i].deallocateAllNodes();
	
	/*#if BBL_DEBUG_MODE
		assert(n == nNodes);
	#endif */
	
	//nNodes = 0;
	
	
	return n;
}


void BBL_MTNodeListManager::deallocate()
{
	BBL_secDeleteArray( nodeLists );
	BBL_secDeleteArray( SEMAPH_lists );
	BBL_secDeleteArray( SEMAPH_listQueries );
}


long unsigned int BBL_MTNodeListManager::getNodePointer(const double zu, BBL_Node* &p, const unsigned int maxLevelPruneCounter)
{
	bool bestFound = false;
	unsigned int i;
	long unsigned int n = 0, nodesPrunedList;
	unsigned int depth, bestDepth = 0;
	unsigned int indexBest;
	int aux;
	double lb, bestlb = INFINITY;
	
	
	p = NULL;
	
	if( expStrat == BBL_ES_BEST_LIMIT )
	{
		for(i = 0; i < nNodeLists; i++)
		{
			SEMAPH_listQueries[i].lock(nThreads);
				aux = nodeLists[i].getFirstNodeLowerBound(lb);
				
				if(aux == 0 && lb <= bestlb)
				{
					if(bestFound)
						SEMAPH_listQueries[indexBest].unlock(nThreads); //we do not need current list indexBestlb more
					
					bestlb = lb;
					indexBest = i;
					bestFound = true;
					//we maintain list i locked because by now it is the list whose we will get the node
				}
				else
				{
			SEMAPH_listQueries[i].unlock(nThreads);
				}
		}
		
		
		if(bestFound)
		{
			SEMAPH_lists[indexBest].lock(nThreads);// BBL_lockMutex(nThreads, &listSemphs[indexBest]);
				n += nodeLists[indexBest].getNodePointer(zu, p, maxLevelPruneCounter);
			SEMAPH_lists[indexBest].unlock(nThreads);// BBL_unlockMutex(nThreads, &listSemphs[indexBest]);
			
			
			SEMAPH_listQueries[indexBest].unlock(nThreads); // BBL_unlockMutex(nThreads, &listQuerySemphs[indexBest]);
		}
		
		
		
	}
	else
	{
		//if we use the mutexex listQuerySemphs here, we can got a error at exploration order..
		for(i = 0; i < nNodeLists; i++)
		{
			SEMAPH_lists[i].lock(nThreads);
			
				aux = nodeLists[i].getFirstNodeDepth(zu, depth, nodesPrunedList);
				
				n += nodesPrunedList;
				
				if(aux == 0)
				{
					if( expStrat == BBL_ES_DEPTH) 
					{
						if(depth >= bestDepth )
						{
							if(bestFound)
								SEMAPH_lists[indexBest].unlock(nThreads); // BBL_unlockMutex(nThreads, &listSemphs[indexBest] ); //we do not need current list indexBest more
							
							bestDepth = depth;
							indexBest = i;
							bestFound = true;
						}
						else
							SEMAPH_lists[i].unlock(nThreads); // BBL_unlockMutex(nThreads, &listSemphs[i]);
					}
					else //if(expStrat == BBL_ES_WIDTH)
					{
						if(bestFound)
						{
							if(depth <= bestDepth)
							{
								SEMAPH_lists[indexBest].unlock(nThreads); // BBL_unlockMutex(nThreads, &listSemphs[indexBest] );
								bestDepth = depth;
								indexBest = i;
							}
							else
								SEMAPH_lists[i].unlock(nThreads); // BBL_unlockMutex(nThreads, &listSemphs[i]);
						}
						else
						{
							bestDepth = depth;
							indexBest = i;
							bestFound = true; //we maintain tlist i locked because by now it is the list whose we will get the node
						}
					}
					
				}
				else
				{
					//printf("cai no else! aux: %d\n", aux);
			SEMAPH_lists[i].unlock(nThreads); 
				}
		}
		
		if(bestFound)
		{
				n += nodeLists[indexBest].getNodePointer(zu, p);
			SEMAPH_lists[indexBest].unlock(nThreads); // BBL_unlockMutex(nThreads, &listSemphs[indexBest]);
		}
	}
	
	
	
	/*if(*p)
	{
		nNodes--;
	}
	
	nNodes -= n; */
	return n;
}



long unsigned int BBL_MTNodeListManager::getNodes(long unsigned int numberOfNodes, BBL_Node* &nodes, bool sortNodesByBound)
{
	decltype(numberOfNodes) mynnodes = numberOfNodes;
	BBL_Node *pLastNode = NULL;
	
	
	nodes = NULL;
	
	for( decltype(nNodeLists) i = 0; i < nNodeLists; i++ )
	{
		BBL_Node *myNodes = NULL;
		long unsigned int nodesGotten = 0, auxNodesToGet = 0;
		long unsigned int nNodesAtList;
		
		SEMAPH_listQueries[i].lock(nThreads); 
		{
			SEMAPH_lists[i].lock(nThreads);
			{
				nNodesAtList = nodeLists[i].getNumberOfNodes();
				
				if( nNodesAtList > 0 )
				{
					auxNodesToGet = BBL_min(mynnodes, nNodesAtList);
					
					nodesGotten =nodeLists[i].getFirstNodes(auxNodesToGet, myNodes);
					//we insert the new nodes ot of mutual exclusion zone...
				}
			}
			SEMAPH_lists[i].unlock(nThreads);
		}
		SEMAPH_listQueries[i].unlock(nThreads);
		
		
		#if BBL_DEBUG_MODE
			//std::cout << "nNodesAtList: " << nNodesAtList << " nodesGotten: " << nodesGotten << " auxNodesToGet: " << auxNodesToGet << "\n";
			assert( nodesGotten == auxNodesToGet );
		#endif
		
		
		if( myNodes )
		{
			if(nodes == NULL)
			{
				nodes = pLastNode = myNodes;
			}
			else
			{
				#if BBL_DEBUG_MODE
					assert(pLastNode->next == NULL);
					assert(myNodes->previous == NULL);
				#endif
				
				pLastNode->next = myNodes;
				myNodes->previous = pLastNode;
			}
			
			mynnodes -= auxNodesToGet;
		}
		
		
		if( mynnodes <= 0 ) //we aready got the requested nodes. we do not need updhate pLastNode
			break;
		
		//updating pLastNode
		if(pLastNode)
		{
			for( ;   pLastNode->next != NULL   ;   pLastNode = pLastNode->next );
		}
	}
	
	#if BBL_DEBUG_MODE
		assert( numberOfNodes >= mynnodes );
	#endif
	
	if(sortNodesByBound)
	{ // to sort nodes, we use a BBL_NodeList just to perform insertion algorithm. Since BBL_NodeList use a internal pointer iterator, I believe the insertion algorthim will not be so expensive
		
		BBL_NodeList auxNodeList(BBL_ES_BEST_LIMIT);
		BBL_Node *p1, *p2;
		
		for( p1 = nodes ; p1 != NULL ;  )
		{
			p2 = p1->next;
			
			p1->previous = NULL;
			p1->next = NULL;
			
			auxNodeList.insertNodes(1, p1);
			
			p1 = p2;
		}
		
		//now, we get the pointer to list's head
		nodes = auxNodeList.head;
		
		auxNodeList.initialize(); //we do not want deallocate our nodes; So, we
	}
	
	
	return (numberOfNodes - mynnodes);
}


long unsigned int BBL_MTNodeListManager::getNumberOfNodes(void ) const
{
	long unsigned int n = 0;
	unsigned int i;
	
	for(i = 0; i < nNodeLists; i++)
		n += nodeLists[i].getNumberOfNodes();
	
	return n;
}


//function to insert nodes that can have diferent lower bounds
void BBL_MTNodeListManager::insertNodesDifsLowerBounds(const unsigned int listNumber, const unsigned int numberOfNodes, BBL_Node* nodes)
{
	unsigned int i = 1;
	BBL_Node *auxNodes = nodes;
	BBL_Node *next = auxNodes->next;
	
	#if BBL_DEBUG_MODE
		unsigned int nnodes = 0;
	#endif
	
	
	
	//we insert nodes having same lower bound togheter. If the lower bound is differente, we have to separte them..
	
	while(true)
	{
		if( next != NULL )
		{
			if( auxNodes->getLowerBound() == next->getLowerBound() )
			{
				i++;
			}
			else
			{
				next->previous->next = NULL;
				next->previous = NULL;
				
				insertNodes(listNumber, i, auxNodes);
				
				#if BBL_DEBUG_MODE
					nnodes += i;
				#endif
				
				i = 1;
				auxNodes = next;
			}
			
			next = next->next;
		}
		else
		{
			break;
		}
	}
	
	#if BBL_DEBUG_MODE
		nnodes += i;
		//std::cout << "nnodes: " << nnodes << " numberOfNodes: " << numberOfNodes << "\n";
		assert( nnodes == numberOfNodes );
	#endif
	
	insertNodes(listNumber, i, auxNodes);
}


void BBL_MTNodeListManager::insertNodes(const unsigned int listNumber, const unsigned int numberOfNodes, BBL_Node* nodes)
{
	SEMAPH_lists[listNumber].lock(nThreads); //BBL_lockMutex(nThreads, &listSemphs[listNumber]);
	
		nodeLists[listNumber].insertNodes(numberOfNodes, nodes);
	
	SEMAPH_lists[listNumber].unlock(nThreads); //BBL_unlockMutex(nThreads, &listSemphs[listNumber]);
	
	//nNodes += numberOfNodes;
}


void BBL_MTNodeListManager::print(const int listNumber) const
{
	unsigned int i;
	
	if(listNumber >= 0)
	{
		nodeLists[listNumber].print();
	}
	else
	{
		for(i = 0; i < nNodeLists; i++)
		{
			
			if( nodeLists[i].getNumberOfNodes() > 0 )
			{
				printf("##########################################################################\n");
				printf("Node list %u:\n", i);
				nodeLists[i].print();
				printf("##########################################################################\n");
			}
		}
	}
}



unsigned int BBL_MTNodeListManager::pruneNodesByBound(const double zu, double& zl, const unsigned int maxLevelPruneCounter)
{
	unsigned int i, nodesPruned = 0;
	double zl_lists = INFINITY, myzl;
	
	
	for(i = 0; i < nNodeLists; i++)
	{
		SEMAPH_lists[i].lock(nThreads);// BBL_lockMutex(nThreads, &listSemphs[i]);
			nodesPruned += nodeLists[i].pruneNodesByBound(zu, myzl, maxLevelPruneCounter);
		SEMAPH_lists[i].unlock(nThreads); // BBL_unlockMutex(nThreads, &listSemphs[i]);
		
		if( !isinf(myzl) && myzl < zl_lists )
			zl_lists = myzl;
	}
	
	if( isinf(zl_lists) )
		zl = -INFINITY; //zl_list has the value INFINITY
	else
		zl = zl_lists;
	
	
	//nNodes -= nodesPruned;
	
	return nodesPruned;
}


void BBL_MTNodeListManager::reorganizeToBestLimit(const double begLBInterval, const double endLBInterval )
{
	if( expStrat == BBL_ES_BEST_LIMIT)
		return;
	
	
	for( decltype(nNodeLists) i = 0; i < nNodeLists; i++ )
	{
		SEMAPH_listQueries[i].lock(nThreads); //BBL_lockMutex(nThreads, &listQuerySemphs[i]);
			SEMAPH_lists[i].lock(nThreads); //BBL_lockMutex(nThreads, &listSemphs[i]);
		
				nodeLists[i].reorganizeToBestLimit(begLBInterval, endLBInterval);
				nodeLists[i].reorganizeNodeLists();
			
			SEMAPH_lists[i].unlock(nThreads); //BBL_unlockMutex(nThreads, &listSemphs[i]);
		SEMAPH_listQueries[i].unlock(nThreads); //BBL_unlockMutex(nThreads, &listQuerySemphs[i]);
	}
	
	expStrat = BBL_ES_BEST_LIMIT;
}




void BBL_MTNodeListManager::reorganizeNodeLists(const unsigned int listNumber)
{
	SEMAPH_lists[listNumber].lock(nThreads); //BBL_lockMutex(nThreads, &listSemphs[listNumber]);
		nodeLists[listNumber].reorganizeNodeLists();
	SEMAPH_lists[listNumber].unlock(nThreads); //BBL_unlockMutex(nThreads, &listSemphs[listNumber]);
}


BBL_MTNodeListManager::~BBL_MTNodeListManager()
{
	   deallocate();
}







































