
#include <cmath> //we include math.h to use INFINITY


#include <cstdlib>
#include <cstdio>
#include <climits>
#include <cassert>

#include <cstdint>



#include "BBL_branchAndBound.hpp"
#include "BBL_tools.hpp"


using namespace std;
using namespace branchAndBound;



#include <new>



#include "BBL_node.hpp"
#include "BBL_tools.hpp"


using namespace std;
using namespace branchAndBound;




#if 0
BBL_NodeBase::BBL_NodeBase()
{
    //nBounds = 0;
    //nBranchVars = 0;
    //bounds = NULL;
    //lb = -INFINITY;
    next = previous = NULL;
    //lastBranch = NULL;
    //xParent = NULL;
    //dualDad = NULL;
}

void BBL_NodeBase::print(std::ostream& out) const
{
    /* unsigned int i;
    for(i = 0; i < nBounds; i++)
    {
        printf("\tvar: %u l: %f u: %f  ", bounds[i].ind, bounds[i].l, bounds[i].u );
    }
    printf("\n"); */
}

int BBL_NodeBase::getParentSolution(const unsigned int n, double* sol) const
{
    if( !xParent )
        return BBL_BAD_DEFINITIONS;
    
    
    BBL_copySequence(n, xParent->a, sol);
    
    return 0;
}

void BBL_NodeBase::setxParentPointer( BBL_Array<double> *p)
{
    xParent = p;
    p->incPointerCounter();
}

void BBL_NodeBase::deallocate(void)
{
    deallocateSharedData();
}

void BBL_NodeBase::deallocateSharedData()
{
    if(xParent)
    {
        xParent->decPointerCounter();
        xParent = NULL;
    }
}

bool BBL_NodeBase::hasSharedData() const
{
    if( xParent )
    {
        if( xParent->getnPointers( ) > 1 )
            return true;
    }
    
    return false;
}

BBL_NodeBase::~BBL_NodeBase()
{
    deallocate();
}
#endif



#if 0
BBL_FloatOrDoubleNodeBoundsPointer::BBL_FloatOrDoubleNodeBoundsPointer (BBL_UNION_FLOAT_OR_DOUBLE_NODE_BOUNDS_POINTER type)
{
    this->type = type;
    value.floatNodeBounds = NULL; //to initialize the union
}


BBL_FloatOrDoubleNodeBoundsPointer:: ~BBL_FloatOrDoubleNodeBoundsPointer()
{
    if(type == BBL_UFDNBP_FLOAT)
    {
        if(value.floatNodeBounds)
            value.floatNodeBounds->decPointerCounter();
    }
    else
    {
        if(value.doubleNodeBounds)
            value.doubleNodeBounds->decPointerCounter();
    }
}


void BBL_FloatOrDoubleNodeBoundsPointer::getVarBoundsOnNode(double *nlx, double *nux) const
{
    if( type == BBL_UFDNBP_FLOAT )
    {
        BBL_getVarBoundsOnNodeFromArraySize( *value.floatNodeBounds, nlx, nux );
    }
    else
    {
        assert(value.doubleNodeBounds);
        BBL_getVarBoundsOnNodeFromArraySize( *value.doubleNodeBounds, nlx, nux );
    }
}


void BBL_FloatOrDoubleNodeBoundsPointer:: print( std::ostream &out ) const
{
    if( type == BBL_UFDNBP_FLOAT )
    {
        BBL_printArraySizeNodeBounds(*value.floatNodeBounds, out);
    }
    else
    {
        BBL_printArraySizeNodeBounds(*value.doubleNodeBounds, out);
    }
}


int BBL_FloatOrDoubleNodeBoundsPointer:: allocateArraySize()
{
    void *p;
    
    if(type == BBL_UFDNBP_FLOAT)
    {
        p = value.floatNodeBounds = new (std::nothrow) BBL_ArraySize<BBL_FloatNodeBounds>;
    }
    else
    {
        p = value.doubleNodeBounds = new (std::nothrow) BBL_ArraySize<BBL_DoubleNodeBounds>;
    }
    
    BBL_IFMEMERRORRETURN(!p);
    
    return 0;
}


int BBL_FloatOrDoubleNodeBoundsPointer:: allocateElements(const unsigned int size)
{
    int r;
    
    if(type == BBL_UFDNBP_FLOAT)
    {
        r = value.floatNodeBounds->allocate(size);
        BBL_IFERRORRETURN(r, r);
        
        value.floatNodeBounds->size = size;
    }
    else
    {
        r = value.doubleNodeBounds->allocate(size);
        BBL_IFERRORRETURN(r, r);
        
        value.doubleNodeBounds->size = size;
    }
    
    return 0;
}


int BBL_FloatOrDoubleNodeBoundsPointer:: getArrayElement( const unsigned int arrayIndex, unsigned int *ind, double *l, double *u) const
{
    if(type == BBL_UFDNBP_FLOAT)
    {
        const auto &elem = value.floatNodeBounds->a[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
    }
    else
    {
        const auto &elem = value.doubleNodeBounds->a[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
    }
    
    return 0;
}


int BBL_FloatOrDoubleNodeBoundsPointer:: setArrayElement( const unsigned int arrayIndex, const unsigned int *ind, const double *l, const double *u)
{
    if(type == BBL_UFDNBP_FLOAT)
    {
        auto &elem = value.floatNodeBounds->a[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
    }
    else
    {
        auto &elem = value.doubleNodeBounds->a[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
    }
    
    return 0;
}

#endif










BBL_ClassUnionNodeBoundsPointer:: BBL_ClassUnionNodeBoundsPointer (BBL_UNION_TYPES_NODE_BOUNDS_POINTER type)
{
    this->type = type;
    value.floatNodeBounds = NULL; //to initialize the union
}


BBL_ClassUnionNodeBoundsPointer:: ~BBL_ClassUnionNodeBoundsPointer()
{
    deallocate();
}


void BBL_ClassUnionNodeBoundsPointer::getVarBoundsOnNode( const unsigned nodeBoundsSize, double *nlx, double *nux) const
{
    if( type == BBL_UTNBP_SHORT_SCHAR )
    {
        BBL_getVarBoundsOnNodeFromNodeBounds( nodeBoundsSize, value.shortScharNodeBounds, nlx, nux );
    }
    else if( type == BBL_UTNBP_SHORT_SHORT_INT )
    {
        BBL_getVarBoundsOnNodeFromNodeBounds( nodeBoundsSize, value.shortShortIntNodeBounds, nlx, nux );
    }
    else if( type == BBL_UTNBP_SHORT_FLOAT )
    {
        BBL_getVarBoundsOnNodeFromNodeBounds( nodeBoundsSize, value.shortFloatNodeBounds, nlx, nux );
    }
    else if( type == BBL_UTNBP_SHORT_DOUBLE)
    {
        BBL_getVarBoundsOnNodeFromNodeBounds( nodeBoundsSize, value.shortDoubleNodeBounds, nlx, nux );
    }
    
    
    else if( type == BBL_UTNBP_SCHAR )
    {
        BBL_getVarBoundsOnNodeFromNodeBounds( nodeBoundsSize, value.scharNodeBounds, nlx, nux );
    }
    else if( type == BBL_UTNBP_SHORT_INT )
    {
        BBL_getVarBoundsOnNodeFromNodeBounds( nodeBoundsSize, value.shortIntNodeBounds, nlx, nux );
    }
    else if( type == BBL_UTNBP_FLOAT )
    {
        BBL_getVarBoundsOnNodeFromNodeBounds( nodeBoundsSize, value.floatNodeBounds, nlx, nux );
    }
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        BBL_getVarBoundsOnNodeFromNodeBounds( nodeBoundsSize, value.doubleNodeBounds, nlx, nux );
    }
}


void BBL_ClassUnionNodeBoundsPointer:: print( const unsigned nodeBoundsSize, std::ostream &out ) const
{
    if( type == BBL_UTNBP_SHORT_SCHAR )
        BBL_printNodeBounds(nodeBoundsSize, value.shortScharNodeBounds, out);
    else if( type == BBL_UTNBP_SHORT_SHORT_INT )
        BBL_printNodeBounds(nodeBoundsSize, value.shortShortIntNodeBounds, out);
    else if( type == BBL_UTNBP_SHORT_FLOAT )
        BBL_printNodeBounds(nodeBoundsSize, value.shortFloatNodeBounds, out);
    else if( type == BBL_UTNBP_SHORT_DOUBLE )
        BBL_printNodeBounds(nodeBoundsSize, value.shortDoubleNodeBounds, out);
    
    
    else if( type == BBL_UTNBP_SCHAR )
        BBL_printNodeBounds(nodeBoundsSize, value.scharNodeBounds, out);
    else if( type == BBL_UTNBP_SHORT_INT )
        BBL_printNodeBounds(nodeBoundsSize, value.shortIntNodeBounds, out);
    else if( type == BBL_UTNBP_FLOAT )
        BBL_printNodeBounds(nodeBoundsSize, value.floatNodeBounds, out);
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        BBL_printNodeBounds(nodeBoundsSize, value.doubleNodeBounds, out);
    }
}



int BBL_ClassUnionNodeBoundsPointer:: allocateElements(const unsigned int size)
{
    void *p;
    
    if( type == BBL_UTNBP_SHORT_SCHAR )
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        value.shortScharNodeBounds = (decltype(value.shortScharNodeBounds) ) malloc( size * sizeof(*value.shortScharNodeBounds) ) ;
        p = value.shortScharNodeBounds;
    }
    else if( type == BBL_UTNBP_SHORT_SHORT_INT )
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        value.shortShortIntNodeBounds = (decltype(value.shortShortIntNodeBounds) ) malloc( size * sizeof(*value.shortShortIntNodeBounds) );
        p = value.shortShortIntNodeBounds;
    }
    else if(type == BBL_UTNBP_SHORT_FLOAT)
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        value.shortFloatNodeBounds = (decltype(value.shortFloatNodeBounds) ) malloc( size * sizeof(*value.shortFloatNodeBounds) );
        p = value.shortFloatNodeBounds;
    }
    else if( type == BBL_UTNBP_SHORT_DOUBLE )
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        value.shortDoubleNodeBounds = (decltype(value.shortDoubleNodeBounds)) malloc( size * sizeof(*value.shortDoubleNodeBounds) );
        p = value.shortDoubleNodeBounds;
    }
    
    
    else if( type == BBL_UTNBP_SCHAR )
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        value.scharNodeBounds = (decltype(value.scharNodeBounds) ) malloc( size * sizeof(*value.scharNodeBounds) ) ;
        p = value.scharNodeBounds;
    }
    else if( type == BBL_UTNBP_SHORT_INT )
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        value.shortIntNodeBounds = (decltype(value.shortIntNodeBounds) ) malloc( size * sizeof(*value.shortIntNodeBounds) );
        p = value.shortIntNodeBounds;
    }
    else if(type == BBL_UTNBP_FLOAT)
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        value.floatNodeBounds = (decltype(value.floatNodeBounds) ) malloc( size * sizeof(*value.floatNodeBounds) );
        p = value.floatNodeBounds;
    }
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        value.doubleNodeBounds = (decltype(value.doubleNodeBounds)) malloc( size * sizeof(*value.doubleNodeBounds) );
        p = value.doubleNodeBounds;
    }
    
    BBL_IFMEMERRORRETURN(!p);
    
    return 0;
}


int BBL_ClassUnionNodeBoundsPointer:: getArrayElement( const unsigned int arrayIndex, unsigned int *ind, double *l, double *u) const
{
    if( type == BBL_UTNBP_SHORT_SCHAR )
    {
        const auto &elem = value.shortScharNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
    }
    else if( type == BBL_UTNBP_SHORT_SHORT_INT )
    {
        const auto &elem = value.shortShortIntNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
    }
    else if(type == BBL_UTNBP_SHORT_FLOAT)
    {
        const auto &elem = value.shortFloatNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
    }
    else if(type == BBL_UTNBP_SHORT_DOUBLE )
    {
        const auto &elem = value.shortDoubleNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
    }
    
    
    else if( type == BBL_UTNBP_SCHAR )
    {
        const auto &elem = value.scharNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
    }
    else if( type == BBL_UTNBP_SHORT_INT )
    {
        const auto &elem = value.shortIntNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
    }
    else if(type == BBL_UTNBP_FLOAT)
    {
        const auto &elem = value.floatNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
    }
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        const auto &elem = value.doubleNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
    }
    
    return 0;
}


int BBL_ClassUnionNodeBoundsPointer:: setArrayElement( const unsigned int arrayIndex, const unsigned int *ind, const double *l, const double *u)
{
    if( type == BBL_UTNBP_SHORT_SCHAR )
    {
        auto &elem = value.shortScharNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
    }
    else if(type == BBL_UTNBP_SHORT_SHORT_INT)
    {
        auto &elem = value.shortShortIntNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
    }
    else if(type == BBL_UTNBP_SHORT_FLOAT)
    {
        auto &elem = value.shortFloatNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
    }
    else if(type == BBL_UTNBP_SHORT_DOUBLE)
    {
        auto &elem = value.shortDoubleNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
    }
    
    
    else if( type == BBL_UTNBP_SCHAR )
    {
        auto &elem = value.scharNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
    }
    else if(type == BBL_UTNBP_SHORT_INT)
    {
        auto &elem = value.shortIntNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
    }
    else if(type == BBL_UTNBP_FLOAT)
    {
        auto &elem = value.floatNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
    }
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        
        auto &elem = value.doubleNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
    }
    
    return 0;
}








BBL_ClassUnionNodeBoundsSolPointer:: BBL_ClassUnionNodeBoundsSolPointer (BBL_UNION_TYPES_NODE_BOUNDS_POINTER type)
{
    this->type = type;
    value.floatNodeBounds = NULL; //to initialize the union
}


BBL_ClassUnionNodeBoundsSolPointer:: ~BBL_ClassUnionNodeBoundsSolPointer()
{
    deallocate();
}


void BBL_ClassUnionNodeBoundsSolPointer::getVarBoundsOnNode( const unsigned nodeBoundsSize, double *nlx, double *nux) const
{
    if( type == BBL_UTNBP_SHORT_SCHAR )
    {
        BBL_getVarBoundsOnNodeFromNodeBoundsSol( nodeBoundsSize, value.shortScharNodeBounds, nlx, nux );
    }
    else if( type == BBL_UTNBP_SHORT_SHORT_INT )
    {
        BBL_getVarBoundsOnNodeFromNodeBoundsSol( nodeBoundsSize, value.shortShortIntNodeBounds, nlx, nux );
    }
    else if( type == BBL_UTNBP_SHORT_FLOAT )
    {
        BBL_getVarBoundsOnNodeFromNodeBoundsSol( nodeBoundsSize, value.shortFloatNodeBounds, nlx, nux );
    }
    else if( type == BBL_UTNBP_SHORT_DOUBLE)
    {
        BBL_getVarBoundsOnNodeFromNodeBoundsSol( nodeBoundsSize, value.shortDoubleNodeBounds, nlx, nux );
    }
    
    
    else if( type == BBL_UTNBP_SCHAR )
    {
        BBL_getVarBoundsOnNodeFromNodeBoundsSol( nodeBoundsSize, value.scharNodeBounds, nlx, nux );
    }
    else if( type == BBL_UTNBP_SHORT_INT )
    {
        BBL_getVarBoundsOnNodeFromNodeBoundsSol( nodeBoundsSize, value.shortIntNodeBounds, nlx, nux );
    }
    else if( type == BBL_UTNBP_FLOAT )
    {
        BBL_getVarBoundsOnNodeFromNodeBoundsSol( nodeBoundsSize, value.floatNodeBounds, nlx, nux );
    }
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        BBL_getVarBoundsOnNodeFromNodeBoundsSol( nodeBoundsSize, value.doubleNodeBounds, nlx, nux );
    }
}


void BBL_ClassUnionNodeBoundsSolPointer:: print( const unsigned nodeBoundsSize, std::ostream &out ) const
{
    if( type == BBL_UTNBP_SHORT_SCHAR )
        BBL_printNodeBoundsSol(nodeBoundsSize, value.shortScharNodeBounds, out);
    else if( type == BBL_UTNBP_SHORT_SHORT_INT )
        BBL_printNodeBoundsSol(nodeBoundsSize, value.shortShortIntNodeBounds, out);
    else if( type == BBL_UTNBP_SHORT_FLOAT )
        BBL_printNodeBoundsSol(nodeBoundsSize, value.shortFloatNodeBounds, out);
    else if( type == BBL_UTNBP_SHORT_DOUBLE )
        BBL_printNodeBoundsSol(nodeBoundsSize, value.shortDoubleNodeBounds, out);
    
    
    else if( type == BBL_UTNBP_SCHAR )
        BBL_printNodeBoundsSol(nodeBoundsSize, value.scharNodeBounds, out);
    else if( type == BBL_UTNBP_SHORT_INT )
        BBL_printNodeBoundsSol(nodeBoundsSize, value.shortIntNodeBounds, out);
    else if( type == BBL_UTNBP_FLOAT )
        BBL_printNodeBoundsSol(nodeBoundsSize, value.floatNodeBounds, out);
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        BBL_printNodeBoundsSol(nodeBoundsSize, value.doubleNodeBounds, out);
    }
}



int BBL_ClassUnionNodeBoundsSolPointer:: allocateElements(const unsigned int size)
{
    void *p;
    
    if(size == 0)
    {
        deallocate();
        return 0;
    }
    
    if( type == BBL_UTNBP_SHORT_SCHAR )
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        p = realloc( value.shortScharNodeBounds, size * sizeof(*value.shortScharNodeBounds) );
        BBL_IFMEMERRORRETURN(!p);
        value.shortScharNodeBounds = (decltype(value.shortScharNodeBounds) ) p;
    }
    else if( type == BBL_UTNBP_SHORT_SHORT_INT )
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        p = realloc( value.shortShortIntNodeBounds, size * sizeof(*value.shortShortIntNodeBounds) );
        BBL_IFMEMERRORRETURN(!p);
        value.shortShortIntNodeBounds = (decltype(value.shortShortIntNodeBounds) ) p;
    }
    else if(type == BBL_UTNBP_SHORT_FLOAT)
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        p = realloc( value.shortFloatNodeBounds, size * sizeof(*value.shortFloatNodeBounds) );
        BBL_IFMEMERRORRETURN(!p);
        value.shortFloatNodeBounds = (decltype(value.shortFloatNodeBounds) ) p;
    }
    else if( type == BBL_UTNBP_SHORT_DOUBLE )
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        p = realloc( value.shortDoubleNodeBounds, size * sizeof(*value.shortDoubleNodeBounds) );
        BBL_IFMEMERRORRETURN(!p);
        value.shortDoubleNodeBounds = (decltype(value.shortDoubleNodeBounds) ) p;
    }
    
    
    else if( type == BBL_UTNBP_SCHAR )
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        p = realloc( value.scharNodeBounds, size * sizeof(*value.scharNodeBounds) );
        BBL_IFMEMERRORRETURN(!p);
        value.scharNodeBounds = (decltype(value.scharNodeBounds) ) p;
    }
    else if( type == BBL_UTNBP_SHORT_INT )
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        p = realloc( value.shortIntNodeBounds, size * sizeof(*value.shortIntNodeBounds) );
        BBL_IFMEMERRORRETURN(!p);
        value.shortIntNodeBounds = (decltype(value.shortIntNodeBounds) ) p;
    }
    else if(type == BBL_UTNBP_FLOAT)
    {
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        p = realloc( value.floatNodeBounds, size * sizeof(*value.floatNodeBounds) );
        BBL_IFMEMERRORRETURN(!p);
        value.floatNodeBounds = (decltype(value.floatNodeBounds) ) p;
    }
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        //we cannot use BBL_malloc because we cannot create a reference of packed structured member...
        p = realloc( value.doubleNodeBounds, size * sizeof(*value.doubleNodeBounds) );
        BBL_IFMEMERRORRETURN(!p);
        value.doubleNodeBounds = (decltype(value.doubleNodeBounds) ) p;
    }
    
    
    return 0;
}


int BBL_ClassUnionNodeBoundsSolPointer:: getArrayElement( const unsigned int arrayIndex, unsigned int *ind, double *l, double *u, double *sol) const
{
    if( type == BBL_UTNBP_SHORT_SCHAR )
    {
        const auto &elem = value.shortScharNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
        if(sol)
            *sol = elem.sol;
    }
    else if( type == BBL_UTNBP_SHORT_SHORT_INT )
    {
        const auto &elem = value.shortShortIntNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
        if(sol)
            *sol = elem.sol;
    }
    else if(type == BBL_UTNBP_SHORT_FLOAT)
    {
        const auto &elem = value.shortFloatNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
        if(sol)
            *sol = elem.sol;
    }
    else if(type == BBL_UTNBP_SHORT_DOUBLE )
    {
        const auto &elem = value.shortDoubleNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
        if(sol)
            *sol = elem.sol;
    }
    
    
    else if( type == BBL_UTNBP_SCHAR )
    {
        const auto &elem = value.scharNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
        if(sol)
            *sol = elem.sol;
    }
    else if( type == BBL_UTNBP_SHORT_INT )
    {
        const auto &elem = value.shortIntNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
        if(sol)
            *sol = elem.sol;
    }
    else if(type == BBL_UTNBP_FLOAT)
    {
        const auto &elem = value.floatNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
        if(sol)
            *sol = elem.sol;
    }
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        const auto &elem = value.doubleNodeBounds[arrayIndex];
        if(ind)
            *ind = elem.ind;
        if(l)
            *l = elem.l;
        if(u)
            *u = elem.u;
        if(sol)
            *sol = elem.sol;
    }
    
    return 0;
}


int BBL_ClassUnionNodeBoundsSolPointer:: setArrayElement( const unsigned int arrayIndex, const unsigned int *ind, const double *l, const double *u, const double *sol)
{
    if( type == BBL_UTNBP_SHORT_SCHAR )
    {
        auto &elem = value.shortScharNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
        if(sol)
            elem.sol = *sol;
    }
    else if(type == BBL_UTNBP_SHORT_SHORT_INT)
    {
        auto &elem = value.shortShortIntNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
        if(sol)
            elem.sol = *sol;
    }
    else if(type == BBL_UTNBP_SHORT_FLOAT)
    {
        auto &elem = value.shortFloatNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
        if(sol)
            elem.sol = *sol;
    }
    else if(type == BBL_UTNBP_SHORT_DOUBLE)
    {
        auto &elem = value.shortDoubleNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
        if(sol)
            elem.sol = *sol;
    }
    
    
    else if( type == BBL_UTNBP_SCHAR )
    {
        auto &elem = value.scharNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
        if(sol)
            elem.sol = *sol;
    }
    else if(type == BBL_UTNBP_SHORT_INT)
    {
        auto &elem = value.shortIntNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
        if(sol)
            elem.sol = *sol;
    }
    else if(type == BBL_UTNBP_FLOAT)
    {
        auto &elem = value.floatNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
        if(sol)
            elem.sol = *sol;
    }
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        
        auto &elem = value.doubleNodeBounds[arrayIndex];
        if(ind)
            elem.ind = *ind;
        if(l)
            elem.l = *l;
        if(u)
            elem.u = *u;
        if(sol)
            elem.sol = *sol;
    }
    
    return 0;
}


void BBL_ClassUnionNodeBoundsSolPointer::buildInitialSolution( const unsigned nodeBoundsSize, double *sol) const
{
    if( type == BBL_UTNBP_SHORT_SCHAR )
        BBL_buildInitialSolutionFromNodeBoundsSol(nodeBoundsSize, value.shortScharNodeBounds, sol);
    else if( type == BBL_UTNBP_SHORT_SHORT_INT )
        BBL_buildInitialSolutionFromNodeBoundsSol(nodeBoundsSize, value.shortShortIntNodeBounds, sol);
    else if( type == BBL_UTNBP_SHORT_FLOAT )
        BBL_buildInitialSolutionFromNodeBoundsSol(nodeBoundsSize, value.shortFloatNodeBounds, sol);
    else if( type == BBL_UTNBP_SHORT_DOUBLE )
        BBL_buildInitialSolutionFromNodeBoundsSol(nodeBoundsSize, value.shortDoubleNodeBounds, sol);
    
    
    else if( type == BBL_UTNBP_SCHAR )
        BBL_buildInitialSolutionFromNodeBoundsSol(nodeBoundsSize, value.scharNodeBounds, sol);
    else if( type == BBL_UTNBP_SHORT_INT )
        BBL_buildInitialSolutionFromNodeBoundsSol(nodeBoundsSize, value.shortIntNodeBounds, sol);
    else if( type == BBL_UTNBP_FLOAT )
        BBL_buildInitialSolutionFromNodeBoundsSol(nodeBoundsSize, value.floatNodeBounds, sol);
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        BBL_buildInitialSolutionFromNodeBoundsSol(nodeBoundsSize, value.doubleNodeBounds, sol);
    }
}




void BBL_ClassUnionNodeBoundsSolPointer::addElementOnSortedArray(const unsigned int nodeBoundsSize, unsigned int index, const double l, const double u, const double sol)
{
    if( type == BBL_UTNBP_SHORT_SCHAR )
        BBL_addElementOnSortedNodeBoundsArray(nodeBoundsSize, value.shortScharNodeBounds, index, l, u, sol);
    else if( type == BBL_UTNBP_SHORT_SHORT_INT )
        BBL_addElementOnSortedNodeBoundsArray(nodeBoundsSize, value.shortShortIntNodeBounds, index, l, u, sol);
    else if( type == BBL_UTNBP_SHORT_FLOAT )
        BBL_addElementOnSortedNodeBoundsArray(nodeBoundsSize, value.shortFloatNodeBounds, index, l, u, sol);
    else if( type == BBL_UTNBP_SHORT_DOUBLE )
        BBL_addElementOnSortedNodeBoundsArray(nodeBoundsSize, value.shortDoubleNodeBounds, index, l, u, sol);
    
    
    else if( type == BBL_UTNBP_SCHAR )
        BBL_addElementOnSortedNodeBoundsArray(nodeBoundsSize, value.scharNodeBounds, index, l, u, sol);
    else if( type == BBL_UTNBP_SHORT_INT )
        BBL_addElementOnSortedNodeBoundsArray(nodeBoundsSize, value.shortIntNodeBounds, index, l, u, sol);
    else if( type == BBL_UTNBP_FLOAT )
        BBL_addElementOnSortedNodeBoundsArray(nodeBoundsSize, value.floatNodeBounds, index, l, u, sol);
    else
    {
        #if BBL_DEBUG_MODE
            assert(type == BBL_UTNBP_DOUBLE);
        #endif
        BBL_addElementOnSortedNodeBoundsArray(nodeBoundsSize, value.doubleNodeBounds, index, l, u, sol);
    }
}






BBL_ParentNodeInfo::BBL_ParentNodeInfo(BBL_UNION_TYPES_NODE_BOUNDS_POINTER type) : parentBounds(type)
{
    depth = -1;
    lb = -INFINITY;
    xParent = NULL;
    
    nParentBounds = 0;
}

BBL_ParentNodeInfo::~BBL_ParentNodeInfo()
{
    deallocate();
}

unsigned int BBL_ParentNodeInfo::getNumberOfBounds() const 
{
    return nParentBounds; 
}


unsigned int BBL_ParentNodeInfo::getMaxDepth() const
{
    decltype(depth) aux = -1;
    return aux;
}


int BBL_ParentNodeInfo::getBoundArrayElement(const unsigned int arrayIndex, unsigned int *ind, double *l, double *u) const
{
    return parentBounds.getArrayElement(arrayIndex, ind, l, u);
}

int BBL_ParentNodeInfo::setBoundArrayElement(const unsigned int arrayIndex, const unsigned int *ind, const double *l, const double *u)
{
    return parentBounds.setArrayElement(arrayIndex, ind, l, u);
}


int BBL_ParentNodeInfo::allocateElements(const unsigned int nElements)
{
    int r = parentBounds.allocateElements(nElements);
    if(r)
        return r;
    
    nParentBounds = nElements;
    
    return r;
}

void BBL_ParentNodeInfo::deallocate()
{
    if(xParent)
    { //we cannot use BBL_secFree because we cannot create a reference of packed structured member...
        free(xParent);
        xParent = NULL;
    }
    
    parentBounds.deallocate();
}


int BBL_ParentNodeInfo::setParentSol(const unsigned int sizePrimal, const double *primalSol, const unsigned int sizeDual, const double *dualSol)
{
    if(primalSol)
    {
        const unsigned int size = dualSol ? sizePrimal + sizeDual : sizePrimal ;
        
        //BBL_realloc(xParent, size);
        //we cannot use BBL_realloc because we cannot create a reference of packed structured member...
        if( size == 0 )
        {
            free(xParent);
            xParent = NULL;
        }
        else
        {
            void *p = realloc(xParent, size * sizeof(*xParent) );
            BBL_IFMEMERRORRETURN(!p);
            
            xParent = (decltype(xParent)) p;
        }
        
        
        BBL_IFMEMERRORRETURN(!xParent);
        
        BBL_copySequence(sizePrimal, primalSol, xParent);
        
        if( dualSol )
            BBL_copySequence(sizeDual, dualSol, &xParent[sizePrimal]);
    }
    
    return 0;
}











BBL_Node::BBL_Node(BBL_UNION_TYPES_NODE_BOUNDS_POINTER type) : myBounds(type)
{
    nMyBounds = 0;
    //myBounds = NULL;
    
    parentInfo = NULL;
    //lb = -BBL_INFINITY;
    
    previous = NULL;
    next = NULL;
}



int BBL_Node::allocateNodeBounds( const unsigned int size )
{
    /*myBounds = ( BBL_NodeBoundsSol* ) malloc( size * sizeof(BBL_NodeBoundsSol) );
    //bounds = new (nothrow) BBL_NodeBounds[size];
    if(!myBounds)
    {
        #if BBL_DEBUG_MODE
            BBL_PRINTMEMERROR;
        #endif
        return BBL_MEMORY_ERROR;
    }*/
    
    if( size > getMaxNMyBounds() )
    {
        BBL_PRINTERRORMSG("Number of bounds allocation requested is greater than representation! requested size: " << size << " max representation: " << getMaxNMyBounds() );
        return BBL_BAD_DEFINITIONS;
    }
    
    int r = myBounds.allocateElements(size);
    BBL_IFERRORRETURN(r, r);
    
    nMyBounds = size;
    return 0;
}



int BBL_Node::buildInitialSolution(const unsigned int nprimal, double* sol) const
{
    int ret = getParentPrimalSolution(nprimal, sol);
    BBL_IFERRORRETURN(ret, ret);
    
    //we do not look to parentBounds because if we store the solution of the parent node, this solution already satisfy the parentBounds
    
    myBounds.buildInitialSolution(nMyBounds, sol);
    
    return 0;
}



int BBL_Node::getParentPrimalSolution( const unsigned int nprimal, double* sol ) const
{
    const double *psol = parentInfo->a.xParent;
    
    if( !psol )
        return BBL_BAD_DEFINITIONS;
    
    BBL_copySequence( nprimal, psol, sol );
    
    return 0;
}


//warning: nprimal is the number of primal variables stored at node. If you are not storing primal variables in the node, set nprimal as 0!
int BBL_Node::getParentDualSolution( const unsigned int nprimal, const unsigned int ndual, double *sol ) const
{
    const double *psol = parentInfo->a.xParent;
    
    if( !psol )
        return BBL_BAD_DEFINITIONS;
    
    BBL_copySequence( ndual, &(psol[nprimal]), sol );
    
    return 0;
}



void BBL_Node::getMyVarBoundsOnNode(double *nlx, double *nux) const
{
    /*for(unsigned int i = 0; i < nMyBounds; i++)
    {
        const unsigned int ind = myBounds[i].ind;
        
        nlx[ind] = myBounds[i].l;
        nux[ind] = myBounds[i].u;
    } */
    myBounds.getVarBoundsOnNode(nMyBounds, nlx, nux);
}


//that method set olny positions correspondets to new bounds
void BBL_Node::getVarBoundsOnNode(double *nlx, double *nux) const
{
    if(parentInfo)
    {
        const BBL_ParentNodeInfo &pnb = parentInfo->a;
        pnb.parentBounds.getVarBoundsOnNode(pnb.nParentBounds, nlx, nux);
    }
    
    getMyVarBoundsOnNode(nlx, nux);
}


/*void BBL_Node::setParentBoundsPointer( BBL_ArraySize< BBL_FNodeBounds >* pfloat, BBL_ArraySize< BBL_NodeBounds >* pdouble )
{
    parentBounds.setNodeBoundsPointer(pfloat, pdouble);
}*/


void BBL_Node:: setParentInfoPointer( BBL_BasePointer<BBL_ParentNodeInfo> *parentBounds)
{
    if(this->parentInfo)
        this->parentInfo->decPointerCounter();
    
    this->parentInfo = parentBounds;
    parentBounds->incPointerCounter();
}


void BBL_Node::deallocate()
{
    /*if( myBounds )
    {
        //we cannot use BBL_secFree because we cannot create a reference of packed structured member...
        free(myBounds);
        myBounds = NULL;
    }*/
    myBounds.deallocate();
    nMyBounds = 0;
    deallocateSharedData();
}



void BBL_Node::print( ostream &out ) const
{
    if(parentInfo)
    {
        const BBL_ParentNodeInfo &pnb = parentInfo->a;
        
        out << "\tLower bound: " << pnb.lb << "\n";
    }
    
    out << "\tParent bounds: " << "\n";
    if(parentInfo)
    {
        const BBL_ParentNodeInfo &pnb = parentInfo->a;
        
        pnb.parentBounds.print(pnb.nParentBounds, out);
        //std::cout << "pnb.nRealBounds: " << pnb.nRealBounds << "\n";
    }
    
    out << "\n\tMy bounds: " << nMyBounds << "\n";
    if(nMyBounds > 0)
    {
        //for(i = 0; i < nMyBounds; i++)
            //out << "\tvar: " <<  myBounds[i].ind <<  " l: " << myBounds[i].l << " u: " << myBounds[i].u ;
        myBounds.print(nMyBounds, out);
    }
    out << "\n" ;
}


unsigned int BBL_Node::getNumberOfTotalBounds() const
{
    return (parentInfo ? parentInfo->a.nParentBounds : 0 ) + nMyBounds;
}


long unsigned int BBL_Node::writeVarBoundsInaBufferArray(void *buffer) const
{
    uint32_t ind;
    double l, u;
    long unsigned int bytesWriten = 0;
    
    
    if(parentInfo)
    {
        const unsigned int nParentBounds = parentInfo->a.nParentBounds;
        const BBL_ClassUnionNodeBoundsPointer &pnbp = parentInfo->a.parentBounds;
        
        for(unsigned int i = 0; i < nParentBounds; i++)
        {
            pnbp.getArrayElement(i, &ind, &l, &u);
            
            BBL_writeAndShift(buffer, ind);
            bytesWriten += sizeof(ind);
            
            BBL_writeAndShift(buffer, l);
            bytesWriten += sizeof(l);
            
            BBL_writeAndShift(buffer, u);
            bytesWriten += sizeof(u);
        }
    }
    
    /*for(unsigned int i = 0; i < nMyBounds; i++)
    {
        ind = myBounds[i].ind;
        l = myBounds[i].l;
        u = myBounds[i].u;
        
        BBL_writeAndShift(buffer, ind);
        bytesWriten += sizeof(ind);
        
        BBL_writeAndShift(buffer, l);
        bytesWriten += sizeof(l);
        
        BBL_writeAndShift(buffer, u);
        bytesWriten += sizeof(u);
    }*/
    
    //bytesWriten += BBL_writeVarBoundsInaBufferArrayFromNodeBoundsAndShift(nMyBounds, myBounds, buffer);
    
    
    for(unsigned int i = 0; i < nMyBounds; i++)
    {
        myBounds.getArrayElement(i, &ind, &l, &u, NULL);
        
        BBL_writeAndShift(buffer, ind);
        bytesWriten += sizeof(ind);
        
        BBL_writeAndShift(buffer, l);
        bytesWriten += sizeof(l);
        
        BBL_writeAndShift(buffer, u);
        bytesWriten += sizeof(u);
    }
    
    return bytesWriten;
}



bool BBL_Node::isParentInfoNull()
{
    return parentInfo == NULL;
}


const BBL_BasePointer<BBL_ParentNodeInfo>* BBL_Node::getParentInfo() const
{
    return parentInfo;
}


double BBL_Node::getLowerBound() const
{
    if( parentInfo )
        return parentInfo->a.lb;
    else
        return -INFINITY;
}


unsigned int BBL_Node::getMaxNMyBounds() const
{
    decltype(nMyBounds) aux = -1; //since nMyBounds is unsigned, we hope putting -1 give us the maximum value of the type;
    
    return aux;
}


unsigned int BBL_Node::getDepth() const
{
    decltype(BBL_ParentNodeInfo::depth) depth; //depth can be a short unsigned int. We have to declare in this way since we take advantage a possible overflow (BBL_ParentNodeInfo initializes depth with -1 for root node).
    
    if(parentInfo)
        depth = parentInfo->a.depth + 1;
    else
        depth = 0;
    
    return depth;
}


bool BBL_Node::hasSharedData() const
{
    if( parentInfo )
    {
        if( parentInfo->getnPointers() > 1 )
            return true;
    }
    
    return false;
}


void BBL_Node::deallocateSharedData()
{
    if( parentInfo )
    {
        parentInfo->decPointerCounter();
        parentInfo = NULL;
    }
}


int BBL_Node::setParentSol(const unsigned int sizePrimal, const double *primalSol, const unsigned int sizeDual, const double *dualSol)
{
    return parentInfo->a.setParentSol(sizePrimal, primalSol, sizeDual, dualSol);
}



int BBL_Node::generateParentInfoForChilds( BBL_UNION_TYPES_NODE_BOUNDS_POINTER type,  BBL_BasePointer<BBL_ParentNodeInfo>* &newParentBounds) const
{
    unsigned int i, j, k;
    int r, code;


    newParentBounds = new (std::nothrow) BBL_BasePointer<BBL_ParentNodeInfo>(type);
    BBL_IFMEMERRORGOTOLABEL(!newParentBounds, code, termination);

    newParentBounds->incPointerCounter();

    if( parentInfo )
    {
        const BBL_ParentNodeInfo &parBounds= parentInfo->a;
        const unsigned int nParentBounds = parBounds.getNumberOfBounds();
        unsigned int nbounds;
        const unsigned int nodeDepth = getDepth(); //note getDepth return the parent depth plus 1. Do not replace for parentInfo->a.depth
        
        BBL_ParentNodeInfo &newParBounds = newParentBounds->a;
        
        
        newParBounds.lb = parBounds.lb ;
        
        
        if( nodeDepth < newParBounds.getMaxDepth() ) 
        {
            newParBounds.depth = nodeDepth; //remember nodeDepth is parBounds.depth + 1
        }
        else //getDepth() == newParBounds.getMaxDepth()
        { //we reach the maximum representation of deep. we just reply the depth
            newParBounds.depth = parBounds.depth;
        }
        
        
        //i runs on parentBounds and j runs newBounds
        for(i = j = nbounds = 0; i < nParentBounds && j < nMyBounds; nbounds++)
        {
            unsigned int pbind, mbind;
            
            parBounds.getBoundArrayElement(i, &pbind, NULL, NULL);
            
            myBounds.getArrayElement(j, &mbind, NULL, NULL, NULL);
            
            if(pbind < mbind) //if(pbind < myBounds[j].ind) //if(parentBounds[i].ind < pMyBounds[j].ind)
            {
                i++;
            }
            else
            {
                if( pbind == mbind ) //if( pbind == myBounds[j].ind ) //if( parentBounds[i].ind == pMyBounds[j].ind )
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
            if( j < nMyBounds )
                nbounds += nMyBounds - j;
        }
        
        #if MRQ_DEBUG_MODE
            assert( nbounds <= nNewBounds + nDadBounds );
        #endif
        
        r = newParBounds.allocateElements(nbounds);
        BBL_IFERRORGOTOLABEL(r, code, r, termination);
        
        //i runs on parentBounds and j runs newBounds
        for(i = j = k = 0; i < nParentBounds && j < nMyBounds; k++)
        {
            unsigned int pbind, mbind;
            double pbl, pbu, mbl, mbu;
            
            parBounds.getBoundArrayElement(i, &pbind, &pbl, &pbu);
            myBounds.getArrayElement(j, &mbind, &mbl, &mbu, NULL);
            
            if( pbind < mbind )//if( pbind < myBounds[j].ind )
            {
                newParBounds.setBoundArrayElement(k, &pbind, &pbl, &pbu);
                i++;
            }
            else
            {
                const unsigned int &pbind2 = mbind; //const unsigned int &pbind2 = myBounds[j].ind;
                const double &pbl2 = mbl, &pbu2 = mbu; //const double &pbl2 = myBounds[j].l, &pbu2 = myBounds[j].u;
                
                newParBounds.setBoundArrayElement(k, &pbind2, &pbl2, &pbu2);
                
                if(pbind == mbind)//if(pbind == myBounds[j].ind)//if(parentBounds[i].ind == pMyBounds[j].ind)
                    i++;
                j++;
            }
        }
        
        
        if( i < nParentBounds )
        {
            unsigned int pbind;
            double pbl, pbu;
            
            j = nParentBounds;
            for( ; i < j; i++, k++ )
            {
                parBounds.getBoundArrayElement(i, &pbind, &pbl, &pbu);
                
                newParBounds.setBoundArrayElement(k, &pbind, &pbl, &pbu);
            }
        }
        else
        {
            unsigned int pbind;
            double pbl, pbu;
            
            for( ; j < nMyBounds; j++, k++ )
            {
                //unsigned int pbind = myBounds[j].ind;
                //double pbl = myBounds[j].l, pbu = myBounds[j].u;
                
                myBounds.getArrayElement(j, &pbind, &pbl, &pbu, NULL);
                
                newParBounds.setBoundArrayElement(k, &pbind, &pbl, &pbu);
            }
        }
        
        
    }
    else
    {
        //If we have no parent info, we assume that is the root node. So, its child will have deep 1. Since object newParBounds refers to parent bound, we set depth as 0
        newParentBounds->a.depth = 0;
    }


    code = 0;

    termination:


    return code;
}



BBL_Node::~BBL_Node()
{
    deallocate();
}



#if 0

BBL_Node2::BBL_Node2( BBL_UNION_FLOAT_OR_DOUBLE_NODE_BOUNDS_POINTER type):BBL_Node(type)
{ 
    byteParentFixed = NULL;
    sintParentBounds = NULL;
}


//that method set olny positions correspondets to new bounds
void BBL_Node2::getVarBoundsOnNode(double *nlx, double *nux) const
{
    /*Order to set is:
    * parentBounds
    * sintParentBounds
    * byteParentFixed
    * myBounds
    */
    
    parentBounds.getVarBoundsOnNode(nlx, nux);
    
    if(sintParentBounds)
        BBL_getVarBoundsOnNodeFromArraySize(*sintParentBounds, nlx, nux);
    
    if(byteParentFixed)
    {
        const unsigned int size = byteParentFixed->size;
        const BBL_ByteVarFixed *a = byteParentFixed->a;
        
        for(unsigned int i = 0; i < size; i++)
        {
            auto ind = a[i].index;
            nlx[ind] = nux[ind] = a[i].value;
        }
    }
    
    getMyVarBoundsOnNode(nlx, nux);
}



void BBL_Node2::print( ostream &out ) const
{
    out << "\tParent bounds: \n";
    parentBounds.print(out);
    
    if(sintParentBounds)
    {
        out << "\tshort int Parent bounds: \n";
        BBL_printArraySizeNodeBounds(*sintParentBounds, out);
    }
    
    
    if(byteParentFixed)
    {
        out << "fixed Parent bounds: \n";
        const unsigned int nfix = byteParentFixed->size;
        auto a = byteParentFixed->a;
        
        for(unsigned int i = 0; i < nfix; i++)
        {
            out << "\tvar: " << a[i].index << " value: " << a[i].value;
        }
    }
    
    out << endl << "\tMy bounds: " << nMyBounds << "\n";
    
    for(unsigned int i = 0; i < nMyBounds; i++)
    {
        out << "\tvar: " <<  myBounds[i].ind <<  " l: " << myBounds[i].l << " u: " << myBounds[i].u ;
    }
    out << "\n" ;
}



unsigned int BBL_Node2::getNumberOfTotalBounds() const
{
    return parentBounds.getSize() + (sintParentBounds ? sintParentBounds->size : 0) + (byteParentFixed ? byteParentFixed->size : 0) + nMyBounds;
}



long unsigned int BBL_Node2::writeVarBoundsInaBufferArray(void *buffer) const
{
    const unsigned int nParentBounds = parentBounds.getSize();
    uint32_t ind;
    double l, u;
    long unsigned int bytesWriten = 0;
    
    for(unsigned int i = 0; i < nParentBounds; i++)
    {
        parentBounds.getArrayElement(i, &ind, &l, &u);
        
        BBL_writeAndShift(buffer, ind);
        bytesWriten += sizeof(ind);
        
        BBL_writeAndShift(buffer, l);
        bytesWriten += sizeof(l);
        
        BBL_writeAndShift(buffer, u);
        bytesWriten += sizeof(u);
    }
    
    if(sintParentBounds)
        bytesWriten += BBL_writeVarBoundsInaBufferArrayFromArraySizeAndShift( *sintParentBounds, buffer);
    
    
    if(byteParentFixed)
    {
        const unsigned int bsyze = byteParentFixed->size;
        const BBL_ByteVarFixed *a = byteParentFixed->a;
        
        for(unsigned int i = 0; i < bsyze; i++)
        {
            ind = a[i].index;
            l = a[i].value; //variable is fixed! Lower and upper bound are the same value
            
            BBL_writeAndShift(buffer, ind);
            bytesWriten += sizeof(ind);
            
            BBL_writeAndShift(buffer, l);
            BBL_writeAndShift(buffer, l);
            bytesWriten += 2*sizeof(l);
        }
    }
    
    
    bytesWriten += BBL_writeVarBoundsInaBufferArrayFromNodeBoundsAndShift(nMyBounds, myBounds, buffer);
    
    
    return bytesWriten;
}


#endif





















