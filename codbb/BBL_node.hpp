
#ifndef _BBL_NODE_HPP
#define _BBL_NODE_HPP

#include <cassert>
#include <cstdint>

#include "BBL_constants.hpp"


namespace branchAndBound{
    
    
    #if 0 
    {
        //we avoid create a base class to NodeBOunds because virtual class turn objects larger in the memory. It is necessary more bytes to store an obect of a virtual class
        class BBL_BaseNodeBounds
        {
        public:
            
            inline virtual unsigned int getIndex() const = 0;
            
            inline virtual double getLowerBound() const = 0; //we make this method return double because that is the dominant type. However, some derived class can get this attribute like float, or even int, for example.
            
            inline virtual double getUpperBound() const = 0; //we make this method return double because that is the dominant type. However, some derived class can get this attribute like float, or even int, for example.
            
            inline virtual void setIndex(const unsigned int value) = 0;
            
            inline virtual void setLowerBound(const double value) = 0; //we make this method receive double because that is the dominant type. However, some derived class can get this attribute like float, or even int, for example.
            
            inline virtual void setUpperBound(const double value) = 0; //we make this method receive double because that is the dominant type. However, some derived class can get this attribute like float, or even int, for example.
        };
    }
    #endif
    
    
    #ifdef _MSC_VER
        //for microsoft compiler
        #define BBL_PRE_PACK   __pragma( pack(push, 1) )
        #define BBL_POS_PACK    __pragma( pack(pop))
    #else
        //g++ and icpc compilers
        #define BBL_PRE_PACK  
        #define BBL_POS_PACK   __attribute__((__packed__))
    #endif
    
    
    
    template <class indexType, class numberType>
    BBL_PRE_PACK class BBL_TNodeBounds //: BBL_BaseNodeBounds
    {
    public:
        //bool lastBranch;
        indexType ind; //index of variable
        numberType l, u;  //lower and upper bound to variable
        
        
        /* we do not declarate methods like virtual because virtual classes spends more bytes in the memory
        inline virtual unsigned int getIndex() const override{ return ind;}
        
        inline virtual double getLowerBound() const override{ return l; }
        
        inline virtual double getUpperBound() const override{ return u; }
        
        inline virtual void setIndex(const unsigned int value) override{ ind = value; }
        
        inline virtual void setLowerBound(const double value) override{ l = value; }
        
        inline virtual void setUpperBound(const double value) override{ u = value; }*/
    }BBL_POS_PACK;  //we really need save memory in this structure. Due to packing, if numberType is a char, compiler will waste 2 bytes due to pack. The most part of memory spent by branch-and-bound is due to memory. We can have literally millions of nodes in the memory, and, each node, an array having several BBL_TNodeBounds. So, we re telling to compiler to do not use padding in this structure to do not waste bytes
    
    
    typedef BBL_TNodeBounds<short unsigned int, double> BBL_ShortDoubleNodeBounds;
    typedef BBL_TNodeBounds<short unsigned int, float> BBL_ShortFloatNodeBounds;
    typedef BBL_TNodeBounds<short unsigned int, short int> BBL_ShortShortIntNodeBounds;
    typedef BBL_TNodeBounds<short unsigned int, int8_t> BBL_ShortSCharNodeBounds;
    
    
    typedef BBL_TNodeBounds<unsigned int, double> BBL_DoubleNodeBounds;
    typedef BBL_TNodeBounds<unsigned int, float> BBL_FloatNodeBounds;
    typedef BBL_TNodeBounds<unsigned int, short int> BBL_ShortIntNodeBounds;
    typedef BBL_TNodeBounds<unsigned int, int8_t> BBL_SCharNodeBounds;
    
    typedef BBL_DoubleNodeBounds BBL_NodeBounds;
    
    
    template <class indexType, class numberType>
    BBL_PRE_PACK class BBL_TNodeBoundsSol : public BBL_TNodeBounds<indexType, numberType>
    {
    public:
        double sol; //solution in the parent node (for pseudo cost)
        
        /*BBL_NodeBoundsSol& operator = (const BBL_NodeBounds &other)
        {
            ind = other.ind;
            l = other.l;
            u = other.u;
            
            return *this;
        } */
        
    }BBL_POS_PACK ;
    
    
    typedef BBL_TNodeBoundsSol<short unsigned int, double> BBL_ShortDoubleNodeBoundsSol;
    typedef BBL_TNodeBoundsSol<short unsigned int, float> BBL_ShortFloatNodeBoundsSol;
    typedef BBL_TNodeBoundsSol<short unsigned int, short int> BBL_ShortShortIntNodeBoundsSol;
    typedef BBL_TNodeBoundsSol<short unsigned int, int8_t> BBL_ShortSCharNodeBoundsSol;
    
    
    typedef BBL_TNodeBoundsSol<unsigned int, double> BBL_DoubleNodeBoundsSol;
    typedef BBL_TNodeBoundsSol<unsigned int, float> BBL_FloatNodeBoundsSol;
    typedef BBL_TNodeBoundsSol<unsigned int, short int> BBL_ShortIntNodeBoundsSol;
    typedef BBL_TNodeBoundsSol<unsigned int, int8_t> BBL_SCharNodeBoundsSol;
    
    
    typedef BBL_DoubleNodeBoundsSol BBL_NodeBoundsSol;
    
    
    
    
    #if 0
    { //we do not unify more BBL_Array to save bytes. Virtual classes spend more bytes to be stored
        /*
        * This class is defined to can unify all classes BBL_Array in a same pointer type
        */ 
        class BBL_BaseArray
        {
        public:
            virtual inline void* getPointerToElement(const unsigned int index) const = 0;
            
            virtual inline const unsigned int* getPointerToArraySize()
            {
                return NULL; //BBL_Array inherits from this classe and has no information abaout array size. For this case, we return NULL to flag this information is not available
            }
        };
    }
    #endif
    
    
    /*That class store an array of yourClass and it has a counter to number of appointments.
    *It is like a smart pointer...
    **/
    template <class yourClass>
    class BBL_Array //: public BBL_BaseArray
    {
        unsigned int nPointers;
        
    public:
        yourClass *a;
        
        BBL_Array()
        {
            a = NULL;
            nPointers = 0;
        }
        
        inline void setArray(yourClass *array)
        {
            a = array;
            nPointers = 0;
        }
        
        inline void decPointerCounter(void)
        {
            if( nPointers <= 1 )
                delete this; //suicide! Nobody points to that pointer more...
            else
                nPointers--;
        }
        
        inline void incPointerCounter(void)
        {
            nPointers++;
        }
        
        inline unsigned int getnPointers(void) const
        {
            return nPointers;
        }
        
        
        inline int allocate(const unsigned int size)
        {
            a = new (std::nothrow) yourClass[size];
            
            if(!a)
                return BBL_MEMORY_ERROR;
            
            nPointers = 0;
            return 0;
        }
        
        
        inline void deallocate(void)
        {
            if(a)
            {
                delete[] a;
                a = NULL;
            }
        }
        
        ~BBL_Array()
        {
            deallocate();
        }
    };
    
    
    //in that class, we store the size of the array also...
    template <class yourClass>
    class BBL_ArraySize : public BBL_Array<yourClass>
    {
        
    public:
        
        unsigned int size;
        
        BBL_ArraySize():BBL_Array<yourClass>()
        {
            size = 0;
        }
        
        //we do not desalloc any memory here. So, we do not declare destructors as virtual because virtual classes spends more bytes and we really need save bytes to open nodes list
        ~BBL_ArraySize()
        {
        }
    };
    
    
    /*template <class yourClass>
    class BBL_NodeBoundsArraySize : public BBL_ArraySize< BBL_TNodeBounds<yourClass> > 
    {
    public:
        
        inline BBL_BaseNodeBounds* getIndexPointer(const unsigned int index)
        {
            return &a[index];
        }
    }; */
    
    
    /*I am not sure if that is the best decision, but i decided create a union. The purpose is because in some branch-and-bound executions, we can really have millions of nodes and, in this way, overcome the memory. SO, to try save memory, we are trying to give the possibility of store parent bounds using floats for lower and upper bounds to variables, instead of doubles. So, we create a union to put in the node. User will decide if it wants use the float array or double array. I hope God helps me here! Ok, I could use derived classes from BBL_NodeBase, one class to float array, and another to double one, but, it would bring some problems to create class BBL_Node2 that must derived from Node. If we had 2 classe for BB_Node, we would need creat two classes of BBL_Node2, one for each class in the ramification*/
    
    #if 0
    union BBL_UnionFloatOrDoubleNodeBoundsPointer
    {
        BBL_ArraySize<BBL_FloatNodeBounds> *floatNodeBounds;
        BBL_ArraySize<BBL_DoubleNodeBounds> *doubleNodeBounds;
    };
    
    
    /*creating a class to our union be descrtivie*/
    BBL_PRE_PACK class BBL_FloatOrDoubleNodeBoundsPointer
    {
        char type; //flag to point what is the type used inside union. We use: 0 for float and 1 for double
        
    public:
        
        union BBL_UnionFloatOrDoubleNodeBoundsPointer value;
        
        
        BBL_FloatOrDoubleNodeBoundsPointer(BBL_UNION_FLOAT_OR_DOUBLE_NODE_BOUNDS_POINTER type);
        
        ~BBL_FloatOrDoubleNodeBoundsPointer();
        
        
        inline void getNodeBoundsPointer( BBL_ArraySize<BBL_FloatNodeBounds>* &floatPointer, BBL_ArraySize<BBL_DoubleNodeBounds>* &doublePointer)
        {
            if( type == BBL_UFDNBP_FLOAT )
            {
                floatPointer = value.floatNodeBounds;
                doublePointer = NULL;
            }
            else
            {
                floatPointer = NULL;
                doublePointer = value.doubleNodeBounds;
            }
        }
        
        
        inline void setNodeBoundsPointer( BBL_ArraySize<BBL_FloatNodeBounds>* floatPointer, BBL_ArraySize<BBL_DoubleNodeBounds>* doublePointer)
        {
            if( type == BBL_UFDNBP_FLOAT )
            {
                if(value.floatNodeBounds)
                    value.floatNodeBounds->decPointerCounter();
                
                if(floatPointer)
                    floatPointer->incPointerCounter();
                
                value.floatNodeBounds = floatPointer;
            }
            else
            {
                if(value.doubleNodeBounds)
                    value.doubleNodeBounds->decPointerCounter();
                
                if(doublePointer)
                    doublePointer->incPointerCounter();
                
                value.doubleNodeBounds = doublePointer;
            }
        }
        
        inline void setNodeBoundsPointer( BBL_FloatOrDoubleNodeBoundsPointer &nodeBounds)
        {
            #if BBL_DEBUG_MODE
                assert(type == nodeBounds.type);
            #endif
            
            setNodeBoundsPointer( nodeBounds.value.floatNodeBounds, nodeBounds.value.doubleNodeBounds );
        }
        
        
        inline void incPointerCounter()
        {
            if( type == BBL_UFDNBP_FLOAT )
                value.floatNodeBounds->incPointerCounter();
            else
                value.doubleNodeBounds->incPointerCounter();
        }
        
        
        /*inline BBL_BaseArray* getBaseArrayVarBounds()
        {
            if(type == BBL_UFDNBP_FLOAT)
                return value.floatNodeBounds;
            else
                return value.doubleNodeBounds;
        }*/
        
        
        inline bool isNodeBoundsNull()
        {
            if(type == BBL_UFDNBP_FLOAT)
                return value.floatNodeBounds == NULL;
            else
                return value.doubleNodeBounds == NULL;
        }
        
        
        inline unsigned int getSize() const
        {
            if(type == BBL_UFDNBP_FLOAT)
            {
                return value.floatNodeBounds ? value.floatNodeBounds->size : 0;
            }
            else
            {
                return value.doubleNodeBounds ? value.doubleNodeBounds->size : 0;
            }
        }
        
        
        void getVarBoundsOnNode(double *nlx, double *nux) const;
        
        void print( std::ostream &out = std::cout) const;
        
        int allocateArraySize();
        
        int allocateElements(const unsigned int size);
        
        int getArrayElement(const unsigned int arrayIndex, unsigned int *ind, double *l, double *u) const;
        
        int setArrayElement(const unsigned int arrayIndex, const unsigned int *ind, const double *l, const double *u);
        
    }BBL_POS_PACK ;
    
    #endif
    
    
    
    union BBL_UnionNodeBoundsPointer
    {
        BBL_ShortSCharNodeBounds 	*shortScharNodeBounds;
        BBL_ShortShortIntNodeBounds *shortShortIntNodeBounds;
        BBL_ShortFloatNodeBounds 	*shortFloatNodeBounds;
        BBL_ShortDoubleNodeBounds 	*shortDoubleNodeBounds;
        
        BBL_SCharNodeBounds 	*scharNodeBounds;
        BBL_ShortIntNodeBounds 	*shortIntNodeBounds;
        BBL_FloatNodeBounds 	*floatNodeBounds;
        BBL_DoubleNodeBounds 	*doubleNodeBounds;
    };
    
    
    
    /*creating a class to our union be descritivie*/
    BBL_PRE_PACK class BBL_ClassUnionNodeBoundsPointer
    {
        int8_t type; //flag to point what is the type used inside union. We use: 0 for float and 1 for double
        
    public:
        
        union BBL_UnionNodeBoundsPointer value;
        
        
        BBL_ClassUnionNodeBoundsPointer(BBL_UNION_TYPES_NODE_BOUNDS_POINTER type);
        
        ~BBL_ClassUnionNodeBoundsPointer();
        
        
        inline void setNodeBoundsPointer( const BBL_ClassUnionNodeBoundsPointer &nodeBounds)
        {
            #if BBL_DEBUG_MODE
                assert(type == nodeBounds.type);
            #endif
            
            value = nodeBounds.value;
            
            //setNodeBoundsPointer( nodeBounds.value.floatNodeBounds, nodeBounds.value.doubleNodeBounds );
        }
        
        
        inline bool isNodeBoundsNull()
        {
            return value.floatNodeBounds == NULL; //we just hav to check if the pointer is NULL, it does not care what type is
        }
        
        inline void deallocate()
        {
            if(value.floatNodeBounds)
            {
                free(value.floatNodeBounds);
                value.floatNodeBounds = NULL;
            }
        }
        
        void getVarBoundsOnNode(const unsigned nodeBoundsSize, double *nlx, double *nux) const;
        
        void print( const unsigned nodeBoundsSize, std::ostream &out = std::cout) const;
        
        int allocateElements(const unsigned int size);
        
        int getArrayElement(const unsigned int arrayIndex, unsigned int *ind, double *l, double *u) const;
        
        int setArrayElement(const unsigned int arrayIndex, const unsigned int *ind, const double *l, const double *u);
        
    }BBL_POS_PACK ;
    
    
    union BBL_UnionNodeBoundsSolPointer
    {
        BBL_ShortSCharNodeBoundsSol 	*shortScharNodeBounds;
        BBL_ShortShortIntNodeBoundsSol 	*shortShortIntNodeBounds;
        BBL_ShortFloatNodeBoundsSol 	*shortFloatNodeBounds;
        BBL_ShortDoubleNodeBoundsSol 	*shortDoubleNodeBounds;
        
        BBL_SCharNodeBoundsSol 		*scharNodeBounds;
        BBL_ShortIntNodeBoundsSol 	*shortIntNodeBounds;
        BBL_FloatNodeBoundsSol 		*floatNodeBounds;
        BBL_DoubleNodeBoundsSol 	*doubleNodeBounds;
    };
    
    /*creating a class to our union be descritivie*/
    BBL_PRE_PACK class BBL_ClassUnionNodeBoundsSolPointer
    {
    public:
        union BBL_UnionNodeBoundsSolPointer value; //we declare "value" before "type" to try making advantage the align in the memory
        
    private:
        
        int8_t type; //flag to point what is the type used inside union. We use: 0 for float and 1 for double
        
    public:
        
        BBL_ClassUnionNodeBoundsSolPointer(BBL_UNION_TYPES_NODE_BOUNDS_POINTER type);
        
        ~BBL_ClassUnionNodeBoundsSolPointer();
        
        
        inline void setNodeBoundsPointer( const BBL_ClassUnionNodeBoundsSolPointer &nodeBounds)
        {
            #if BBL_DEBUG_MODE
                assert(type == nodeBounds.type);
            #endif
            
            value = nodeBounds.value;
            
            //setNodeBoundsPointer( nodeBounds.value.floatNodeBounds, nodeBounds.value.doubleNodeBounds );
        }
        
        
        inline bool isNodeBoundsNull()
        {
            return value.floatNodeBounds == NULL; //we just hav to check if the pointer is NULL, it does not care what type is
        }
        
        inline void deallocate()
        {
            if(value.floatNodeBounds)
            {
                free(value.floatNodeBounds);
                value.floatNodeBounds = NULL;
            }
        }
        
        void getVarBoundsOnNode(const unsigned nodeBoundsSize, double *nlx, double *nux) const;
        
        void print( const unsigned nodeBoundsSize, std::ostream &out = std::cout) const;
        
        int allocateElements(const unsigned int size);
        
        int getArrayElement(const unsigned int arrayIndex, unsigned int *ind, double *l, double *u, double *sol) const;
        
        int setArrayElement(const unsigned int arrayIndex, const unsigned int *ind, const double *l, const double *u, const double *sol);
        
        void buildInitialSolution(const unsigned nodeBoundsSize, double *sol) const;
        
        void addElementOnSortedArray(const unsigned int nodeBoundsSize, unsigned int index, const double l, const double u, const double sol);
        
    }BBL_POS_PACK ;
    
    
    
    
    
    
    
    
    
    
    BBL_PRE_PACK class BBL_ParentNodeInfo
    {
    public:
        
        double lb;				// lower Bound to parent node
        
        double *xParent;		//Dad's x and dual solution.
        
        BBL_ClassUnionNodeBoundsPointer parentBounds; //pointer to parent bound array bound
        
        unsigned int nParentBounds; 	//number of parent node bounds
        
        #if BBL_USE_SHORT_INTEGER_TO_DEPTH_IN_BB_NODE
            unsigned short int
        #else
            unsigned int
        #endif
        depth;		// parent node's depth
        
        
        BBL_ParentNodeInfo(BBL_UNION_TYPES_NODE_BOUNDS_POINTER type);
        
        ~BBL_ParentNodeInfo();
        
        int allocateElements(const unsigned int nElements);
        
        void deallocate();
        
        int getBoundArrayElement(const unsigned int arrayIndex, unsigned int *ind, double *l, double *u) const;
        
        unsigned int getMaxDepth() const;
        
        unsigned int getNumberOfBounds() const;
        
        int setBoundArrayElement(const unsigned int arrayIndex, const unsigned int *ind, const double *l, const double *u);
        
        
        int setParentSol(const unsigned int sizePrimal, const double *primalSol, const unsigned int sizeDual = 0, const double *dualSol = NULL);
        
    }BBL_POS_PACK ;
    
    
    
    
    /*That class store an array of yourClass and it has a counter to number of appointments.
    *It is like a smart pointer...
    **/
    template <class yourClass>
    BBL_PRE_PACK class BBL_BasePointer 
    {
        unsigned short int nPointers;
        
    public:
        yourClass a;  //yourClass should have a deallocate method
        
        BBL_BasePointer()
        {
            nPointers = 0;
        }
        
        inline void decPointerCounter(void)
        {
            if( nPointers <= 1 )
                //a.deallocate(); deallocate will be calleb when this obejct will be deleted since a is not a pointer
                delete this; //suicide! Nobody points to that pointer more...
            else
                nPointers--;
        }
        
        inline void incPointerCounter(void)
        {
            nPointers++;
        }
        
        inline unsigned int getnPointers(void) const
        {
            return nPointers;
        }
        
        
        ~BBL_BasePointer()
        {
            //deallocate();
        }
    }BBL_POS_PACK ;
    
    
    
    
    template <>
    BBL_PRE_PACK class BBL_BasePointer< BBL_ParentNodeInfo>
    {
        unsigned short int nPointers;
        
    public:
        
        BBL_ParentNodeInfo a;
        
        BBL_BasePointer(BBL_UNION_TYPES_NODE_BOUNDS_POINTER type):a(type)
        {
            nPointers = 0;
        }
        
        inline void decPointerCounter(void)
        {
            if( nPointers <= 1 )
            {
                //a.deallocate(); deallocate will be calleb when this obejct will be deleted since a is not a pointer
                delete this; //suicide! Nobody points to that pointer more...
            }
            else
                nPointers--;
        }
        
        inline void incPointerCounter(void)
        {
            nPointers++;
        }
        
        inline unsigned int getnPointers(void) const
        {
            return nPointers;
        }
        
        
        ~BBL_BasePointer()
        {
            //deallocate();
        }
    }BBL_POS_PACK ;
    
    
    
    
    
    #if 0
    class BBL_NodeBase
    {
        
    public:
        
        //unsigned int depth;		//node's depth
        //double lb;			//Lower Bound to node
        
        BBL_NodeBase *previous, *next;	//to encadeate list of nodes
        
        
        
        
        BBL_NodeBase();
        
        virtual ~BBL_NodeBase();
        
        
        virtual int buildInitialSolution(const unsigned int nprimal, double* sol) const = 0;
        
        virtual void print(std::ostream &out = std::cout) const;
        
        int getParentSolution(const unsigned int n, double* sol) const;
        
        void setxParentPointer( BBL_Array<double> *p);
        
        //void setdualDadPointer( BBL_Array<double> *p);
        
        virtual void deallocateSharedData();
        
        virtual void deallocate(void);
        
        
        //return true if some data is shared between other objects by means of some BBL_Array
        virtual bool hasSharedData() const;
        
        //void setlastBranchPointer( BBL_ArraySize<unsigned int> *p);
        
        //that method set olny positions correspondets to new bounds
        virtual void getVarBoundsOnNode(double *nlx, double *nux) const = 0;
        
        
    protected:
        
        //virtual int allocateNodeBounds( const unsigned int size ) = 0;
        
    };
    #endif
    
    
    
    
    BBL_PRE_PACK class BBL_Node
    {
    
    public:
        
        BBL_Node *previous, *next;	//to encadeate list of nodes (declare before myBounds to keep this data aligned in the memory and optimize performance t manage open node list)
        
        
        #if BBL_USE_SHORT_INTEGER_TO_NUMBER_OF_MY_BOUNDS_IN_BB_NODE
            unsigned short int
        #else
            unsigned int
        #endif
        nMyBounds;
        
        BBL_ClassUnionNodeBoundsSolPointer myBounds; //Bounds to variables...
        
        
    protected:
        
        BBL_BasePointer<BBL_ParentNodeInfo> *parentInfo;
        
        
    public:
        
        
        BBL_Node(BBL_UNION_TYPES_NODE_BOUNDS_POINTER type);
        
        
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        ~BBL_Node();
        
        //allocate memory only for the bounds for this node. This method does not allocate memory for bounds inherits from parent. 
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        int allocateNodeBounds( const unsigned int size );
        
        
        //that method build a initial soultion using parent solution stored at node. The solution is updated to respect new bounds in this node. If parent solution is not stored, return BBL_BAD_DEFINITIONS
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        int buildInitialSolution(const unsigned int nprimal, double* sol) const;
        
        
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        void deallocate();
        
        
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        void deallocateSharedData();
        
        
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        int generateParentInfoForChilds( BBL_UNION_TYPES_NODE_BOUNDS_POINTER type,  BBL_BasePointer<BBL_ParentNodeInfo>* &newParentBounds ) const;
        
        
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        unsigned int getDepth() const;
        
        
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        double getLowerBound() const;
        
        
        unsigned int getMaxNMyBounds() const;
        
        
        //get the number of indices having bounds changed in this nodes. Note: we can have repetiions here and repetions are counted also
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        unsigned int getNumberOfTotalBounds() const;
        
        
        //warning: nprimal is the number of primal variables stored at node. If you are not storing primal variables in the node, set nprimal as 0!
        int getParentDualSolution( const unsigned int nprimal, const unsigned int ndual, double *sol ) const;
        
        
        const BBL_BasePointer<BBL_ParentNodeInfo>* getParentInfo() const;
        
        
        int getParentPrimalSolution( const unsigned int nprimal, double *sol ) const;
        
        //that method set olny positions correspondets to new bounds
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        void getVarBoundsOnNode(double *nlx, double *nux) const;
        
        
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        bool hasSharedData() const;
        
        
        
        bool isParentInfoNull();
        
        
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        void print(std::ostream& out = std::cout) const;
        
        //void setParentBoundsPointer( BBL_ArraySize< BBL_FNodeBounds >* pfloat, BBL_ArraySize< BBL_NodeBounds >* pdouble );
        
        
        
        void setParentInfoPointer(BBL_BasePointer<BBL_ParentNodeInfo> *parentBounds); 
        
        
        
        int setParentSol(const unsigned int sizePrimal, const double *primalSol, const unsigned int sizeDual = 0, const double *dualSol = NULL);
        
        
        /* this method is due to dctools and it writes all bounds stored in the node in buffer array. Order is:
        * 1 - Parent inherit Node Bounds
        * 2 - Own Node bounds
        * 
        * Let nb be the total number of bounds. So, we will write in buffer a nb-tuple having:
        * < variable index (uint32), lower bound (double), upper bound (double) >
        * 
        * Note the same variable can have bounds sets in Parent bounds and node bounds. So, its bounds will apear two times in the buffer
        */
        #if !BBL_DO_NOT_SET_BBL_NODE_AS_VIRTUAL_CLASS
        virtual
        #endif 
        long unsigned int writeVarBoundsInaBufferArray(void *buffer) const;
        
        
        
    protected:
        
        void getMyVarBoundsOnNode(double *nlx, double *nux) const;
        
        
    }BBL_POS_PACK ;
    
    
    
    /*//class to box constraints that can be represented using short integer
    class BBL_SIntNodeBounds
    {
    public:
        unsigned int index;
        short int l, u;
    }; */
    
    
    #if 0
    /* This class is design to store nodes in integer problems. BBL_Node could be used for this proposit also, but BBL_Node2 tries optimize the storage spending less bytes. If your are implementing a pure spatial B&B having no integer variables, BB_Node is more appropriated to you.
    * The objective here is to save bytes on representation, enabling BBL to be used in large problems by more cpu time whitout overcome the system memory.
    */
    
    //class to represent variables being fixed at integer values between 0 and 255.
    BBL_PRE_PACK class BBL_ByteVarFixed
    {
    public:
        unsigned int index;
        unsigned char value;
    }BBL_POS_PACK ;
    
    
    class BBL_Node2 : public BBL_Node
    {
    public:
        
        /*Order to set is:
        * parentBounds
        * sintParentBounds
        * byteParentFixed
        * myBounds
        */
        
        BBL_ArraySize<BBL_ShortIntNodeBounds> *sintParentBounds;
        
        BBL_ArraySize<BBL_ByteVarFixed> *byteParentFixed; //indices of variables in the parent node fixed at some 8 bits value. 
        
        
        
        BBL_Node2( BBL_UNION_FLOAT_OR_DOUBLE_NODE_BOUNDS_POINTER type);
        
        //we do not need the method buildInitialSolution because this method only overwits using myBounds since we assume the parent sol already satisfy the parent bounds
        //virtual int buildInitialSolution(const unsigned int nprimal, double* sol) const override;
        
        //that method set olny positions correspondents to new bounds
        virtual void getVarBoundsOnNode(double *nlx, double *nux) const override;
        
        virtual void print(std::ostream& out = std::cout) const override;
        
        //get the number of indices having bounds changed in this nodes. Note: we can have repetiions here and repetions are counted also
        virtual unsigned int getNumberOfTotalBounds() const override;
        
        /*this method is due to dctools and writes all bounds stored in the node in buffer array. Order is:
        * 1 - Parent inherit Node Bounds
        * 2 - Own Node bounds
        * 
        * Let nb be the total number of bounds. So, we will write in buffer a nb-tuple having:
        * < variable index (uint32), lower bound (double), upper bound (double) >
        * 
        * Note the same variable can have bounds sets in Parent bounds and node bounds. So, its bounds will apear two times in the buffer
        */
        virtual long unsigned int writeVarBoundsInaBufferArray(void *buffer) const override;
        
        virtual void deallocate() override {}
    };
    
    #endif
    
    inline BBL_Node* BBL_generateNewNode(BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY strategy, const unsigned int maxVars) 
    {
        BBL_UNION_TYPES_NODE_BOUNDS_POINTER unionType = BBL_pnbs2usfdnbp(strategy, maxVars);
        BBL_Node* node;
        
        
        /*if( strategy == BBL_PNBSS_FLOAT_SPECIAL_STRUCTS_TO_INT_VARS || strategy == BBL_PNBSS_DOUBLE_SPECIAL_STRUCTS_TO_INT_VARS)
        {
            //node = new (std::nothrow) BBL_Node2(unionType); //BBL_Node is optimized to store
            assert(false);
        }
        else */
        {
            node = new (std::nothrow) BBL_Node(unionType);
        }
        
        return node;
    }
    
    
    
    
}



#endif
