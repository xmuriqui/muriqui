/* Definitions to build a Muriqui Server on distributed computation
* 
* 
* 
* 
* 
*/

#ifndef MRQ_SERVER_HPP_
#define MRQ_SERVER_HPP_



#include "muriqui.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_advanced.hpp"
#include "DCT_bbserver.hpp"



namespace muriqui
{
    
    
    #define MRQ_CURRENT_NODE_REPRESENTATION_VERSION 1
    
    
    
    class MRQ_SegFaultException : public std::exception { };
    
    
    class MRQ_AbortException : public std::exception {} ;


    
    
    
    
    
    class MRQ_ServerServiceCoreGenerator : public dctools::DCT_ServiceCoreGenerator
    {
    public:
        
        virtual int generateServiceCore(dctools::DCT_ServiceCore* &serviceCore) override;
    };
    
    
    
    
    class MRQ_ServerServiceCore : public dctools::DCT_ServiceCore
    {
        
    public:
        
        MRQ_MINLPProb *prob;
        double *olx, *oux;
        
        
        MRQ_ServerServiceCore();
        
        virtual ~MRQ_ServerServiceCore();
        
        
        void deallocate();
        
        void initialize();
        
        int __run(dctools::DCT_BBServer *server, dctools::DCT_ServerConnection *connection, const long unsigned int serviceNumber, const unsigned int nthreads, dctools::DCT_Int64 nBasicInputParameters, dctools::DCT_Byte *basicInputParameters, dctools::DCT_FileNames &inputFiles, dctools::DCT_AllGeneralParams &allGeneralParams, dctools::DCT_Int64 sizeOfOpenNodeRep, dctools::DCT_Byte *openNodeRep, double lowerBound, double upperBound, dctools::DCT_UInt32 nRequestedNodes, bool &stopService, dctools::DCT_SERVER_CLOSING_CODE &responseCodeToClientIfStop );
        
        
        virtual int run(dctools::DCT_BBServer *server, dctools::DCT_ServerConnection *connection, const long unsigned int serviceNumber, const unsigned int nthreads, dctools::DCT_Int64 nBasicInputParameters, dctools::DCT_Byte *basicInputParameters, dctools::DCT_FileNames &inputFiles, dctools::DCT_AllGeneralParams &allGeneralParams, dctools::DCT_Int64 sizeOfOpenNodeRep, dctools::DCT_Byte *openNodeRep, double lowerBound, double upperBound, dctools::DCT_UInt32 nRequestedNodes, bool &stopService, dctools::DCT_SERVER_CLOSING_CODE &responseCodeToClientIfStop ) override;
        
        
        
        //int __run(dctools::DCT_BBServer *server, dctools::DCT_ServerConnection *connection, const long unsigned int serviceNumber, const unsigned int nthreads, int32_t nBasicInputParameters, int32_t *basicInputParameters, dctools::DCT_FileNames &inputFiles, dctools::DCT_AllGeneralParams &allGeneralParams, dctools::DCT_VarBounds &varsBounds, double lowerBound, double upperBound, dctools::DCT_UInt32 nRequestedNodes, bool &stopService, dctools::DCT_SERVER_CLOSING_CODE &responseCodeToClientIfStop );
        
        
        //virtual int run(dctools::DCT_BBServer *server, dctools::DCT_ServerConnection *connection, const long unsigned int serviceNumber, const unsigned int nthreads, int32_t nBasicInputParameters, int32_t *basicInputParameters, dctools::DCT_FileNames &inputFiles, dctools::DCT_AllGeneralParams &allGeneralParams, dctools::DCT_VarBounds &varsBounds, double lowerBound, double upperBound, dctools::DCT_UInt32 nRequestedNodes, bool &stopService, dctools::DCT_SERVER_CLOSING_CODE &responseCodeToClientIfStop ) override;
        
    };
    
    
    
    //class to represent several open nodes to de sent to client. We use that for MILP solver B&B where we cannot get several nodes in a single call. So we use this class to store several open nodes, and so, we send the several nodes to client.
    class MRQ_OpenNodesStorer
    {
    public:
        
        unsigned int nNodes;
        unsigned int maxNodes;
        dctools::DCT_Byte **nodesRep; 
        dctools::DCT_UInt64 *nBytesInNodeRep; //nBytesInNodeRep[i] has the number os bytes in nodesRep[i]
        
        MRQ_OpenNodesStorer();
        
        ~MRQ_OpenNodesStorer();
        
        void initialize();
        
        int allocate(unsigned int maxNodes);
        
        void deallocateAll();
        
        void deallocateNodes();
        
    };
    
    

    class MRQ_OpenNodeWriter
    {
    public:
        dctools::DCT_UInt64 bufferSize;
        dctools::DCT_Byte *buffer;
        
        MRQ_OpenNodeWriter();
        
        ~MRQ_OpenNodeWriter();
        
        void desallocate();
        
        void initialize();
        
        int reallocateBuffer(long unsigned int newSize);
        
        int write( dctools::DCT_BaseSocket &socket, dctools::DCT_CONECTION_CODE ccode, const MRQ_NewBBNode *node, const dctools::DCT_VarBounds *rootVarBounds = NULL);
        
        //olx and oux are original bounds for variables from the problem definition. nlx and nux are bounds for variables in the respective node.
        int write( dctools::DCT_BaseSocket &socket, dctools::DCT_CONECTION_CODE ccode, const double lb, const double *olx, const double *oux, const double *nlx, const double *nux, const unsigned int nI, const int *intVars );
        
        
        //if necessary, output array will realloc to store the node. So, maxSizeOutputArray is input/output argument. If reallocation is necessary and reallocateWithSlack is true, we reallocate the dobule space necessary to avoid future reallocations. If reallocateWithSlack is false, we justa reallocate the strict necessary space.
        int write( dctools::DCT_Byte* &output, dctools::DCT_UInt64 &maxSizeOutputArray, const bool reallocateWithSlack, const bool writeHeader, dctools::DCT_CONECTION_CODE ccode, const double lb, const double *olx, const double *oux, const double *nlx, const double *nux, const unsigned int nI, const int *intVars, dctools::DCT_UInt64 &finalSizeOutput );
        
        int write( dctools::DCT_BaseSocket &socket, dctools::DCT_CONECTION_CODE ccode, const MRQ_OpenNodesStorer &openNodesStorer);
    };
    
    
    
    class MRQ_ServerBBCallbacks : public MRQ_UserCallbacks
    {
        bool stopAlg;
        bool accNodesToSend; //flag to know if we have to accumulate open nodes to send
        long unsigned int nextIterToSendLowerBound;
        long unsigned int myiter;
        double myzl;
        
        //DCT_VarBoundsWriter varBoundsWriter;
        
    public:
        
        bool requestFirstBestSol;
        bool setLazyConstraintsOnReceivedSol;
        bool stopService;
        unsigned int n;
        unsigned int nI;  //number of integer vars 
        unsigned int nthreads;
        dctools::DCT_UInt32 nRequestedNodes;
        dctools::DCT_Byte *rootOpenNodeRep; //just for the case when server received more than one open node to explore
        int *intVars;	//indices of integer vars
        double *auxSentSol; //n+3 array
        double *auxVars;
        double *nodeBounds; //for get nlx and nux in linear approximation algorithms (only for thread responsing the client)
        double *receivedObjAndSol; //just one thread receives solution
        double *olx, *oux; //original bounds for variables (just pointers)
        dctools::DCT_ServerConnection *connection;
        dctools::DCT_VarBounds *rootVarBounds;
        MRQ_OpenNodeWriter nodeWriter;
        
        unsigned int minNumberOfNodestoSendNode;
        unsigned int frequencyToSendLowerBound;
        unsigned int maxNumberOfNodesToSend; //maximum number of nodes sent to another server in a single call
        double minRelativeGaptoUsePseudoCosts;
        
        MRQ_OpenNodesStorer openNodesStorer;
        
        MRQ_Mutex SEMAPH_connectionWrite;
        
        MRQ_Mutex SEMAPH_startOthersThreadsAfterSendingNodesToOtherServes; //This semaphore is to block threads (except threda 0) to enter in the B&B loop before thread 0 sending the nodes to other servers. We adopt this strategy, to try genertae a balanced division of the space between all servers. So, remaining threads only start after thread 0 divider space between all servers.
        
        MRQ_MINLPProb *prob;
        MRQ_Algorithm *myAlg;
        
        MRQ_ServerBBCallbacks(dctools::DCT_ServerConnection* connection);
        
        virtual ~MRQ_ServerBBCallbacks();
        
        
        virtual int beforeAll(const MRQ_ALG_CODE algCode, const unsigned int numberOfThreads) override;
        
        virtual int BB_generateRootNode( MRQ_NewBBNode* &rootNode ) override;
        
        virtual int beforeBBLoop(const unsigned int threadNumber, const double lb, const double ub, MRQ_NLPSolver &nlpSolver) override;
        
        virtual void afterBBLoop(const unsigned int threadNumber, const double lb, const double ub, MRQ_NLPSolver &nlpSolver, const int threadReturnCode) override;
        
        virtual int BB_beforeSolvingRelaxation(const unsigned int threadNumber, MRQ_NewBBNode &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, MRQ_NLPSolver &nlpSolver, bool &pruneNode) override;
        
        virtual int BB_afterSolvingRelaxation(const unsigned int threadNumber, MRQ_NewBBNode &node, const long unsigned int iter, const double lb, const double ub, const double *nlx, const double *nux, MRQ_NLPSolver &nlpSolver, const int status, double &objFSolution, double &dualObjFSolution, double *sol, double *constrs, double *dualSolC, double *dualSolV, bool &pruneNode, bool &solveAgain, bool &branchEvenIntegerSol) override;
        
        virtual void newBestSolution( const unsigned int threadNumber, const double *newSol, const double oldBestObj, const double newBestObj, const long unsigned int iter ) override;
        
        virtual int linearApp_beforeSolveInMILPBB(const MRQ_ALG_CODE algCode, const unsigned int threadNumber, MRQ_MILPSolverCallbackInterface &milpSolverCallbackInterface, MRQ_SolutionStorer &solsToBuildLazyConstraints ) override; //deprecated! we are not using this function more... TODO: remove this method
        
        
        virtual int linearApp_branchingInMILPBB(const MRQ_ALG_CODE algCode, const unsigned int threadNumber, MRQ_MILPSolverCallbackInterface &milpSolverCallbackInterface) override;
        
        virtual int endOfIteration(const MRQ_ALG_CODE algCode, const unsigned int threadNumber, const long unsigned int iter, const double cpuTime, const double wallTime, const double lb, const double ub) override;
        
        virtual void afterAll(const MRQ_ALG_CODE algCode, const long unsigned int iters, const double cpuTime, const double wallTime, const double lb, const double ub);
        
        
        int allocateStructures(MRQ_Algorithm *alg, MRQ_MINLPProb *prob);
        
        int general_beforeSolvingRelaxation(const unsigned int threadNumber, const MRQ_ALG_CODE algCode, MRQ_MILPSolverCallbackInterface *milpSolverCallbackInterface, MRQ_NewBBNode *node, const long unsigned int iter, const double lb, const double nodelb, double *nlx, double *nux,  bool &solReceived, bool &receivedSolUpdtBestSol, double *receivedObjAndSol, bool &pruneNode);
        
        void deallocate();
    };









    template < class key, class value >
    int MRQ_copyContainerMap(const std::map<key, value> &source, std::map<key, value> &destination)
    {
        try
        {
            for( const auto &pairKeyValue : source )
            {
                destination[ pairKeyValue.first ] = pairKeyValue.second;
            }
        }
        catch( std::bad_alloc &ba )
        {
            MRQ_PRINTMEMERROR;
            return MRQ_MEMORY_ERROR;
        }
        
        return 0;
    }
    
    
    static inline int MRQ_buildParamsFromDCTParams(MRQ_GeneralSolverParams &mrqParams, const dctools::DCT_GeneralParams &dctParams)
    {
        //unfortunatelly, perform an attribution generates errors in some compilers like gcc, so we put our hands to work!
        int r1, r2, r3;
        
        r1 = MRQ_copyContainerMap(dctParams.dblParams, mrqParams.dblParams);
        if(r1 != 0)
        {
            MRQ_PRINTERRORNUMBER(r1);
        }
        
        r2 = MRQ_copyContainerMap(dctParams.intParams, mrqParams.intParams);
        if(r2 != 0)
        {
            MRQ_PRINTERRORNUMBER(r2);
        }
        
        r3 = MRQ_copyContainerMap(dctParams.strParams, mrqParams.strParams);
        if(r3 != 0)
        {
            MRQ_PRINTERRORNUMBER(r3);
        }
        
        /*try
        {
            mrqParams.dblParams = dctParams.dblParams;
            mrqParams.intParams = dctParams.intParams;
            mrqParams.strParams = dctParams.strParams;
        }
        catch(std::bad_alloc &ba)
        {
            MRQ_PRINTMEMERROR;
            return MRQ_MEMORY_ERROR;
        }*/
        
        return (r1 | r2 | r3 ) ? MRQ_UNDEFINED_ERROR : 0;
    }
    
    
    
    
    
    
    inline dctools::DCT_OPTIMIZATION_RETURN_CODE MRQ_retCode2DCT_optRetCode( int mrqRetCode )
    {
        const static DCT_Dict<int, dctools::DCT_OPTIMIZATION_RETURN_CODE> dict = {
            {MRQ_OPTIMAL_SOLUTION , dctools::DCT_ORC_OPTIMAL_SOLUTION},
            {MRQ_INFEASIBLE_PROBLEM , dctools::DCT_ORC_INFEASIBLE_PROBLEM},
            {MRQ_UNBOUNDED_PROBLEM , dctools::DCT_ORC_UNBOUNDED_PROBLEM},
            {MRQ_CALLBACK_FUNCTION_ERROR , dctools::DCT_ORC_EVALUATION_ERROR},
            {MRQ_LIBRARY_NOT_AVAILABLE ,  dctools::DCT_ORC_LIBRARY_NOT_AVAILABLE},
            {MRQ_MAX_ITERATIONS_STOP , dctools::DCT_ORC_MAX_ITERATIONS},
            {MRQ_MAX_TIME_STOP , dctools::DCT_ORC_MAX_TIME},
            {MRQ_MILP_SOLVER_ERROR , dctools::DCT_ORC_MILP_SOLVER_ERROR},
            {MRQ_NLP_SOLVER_ERROR , dctools::DCT_ORC_NLP_SOLVER_ERROR},
            {MRQ_MEMORY_ERROR , dctools::DCT_ORC_MEMORY_ERROR}
        };
        
        if( dict.count(mrqRetCode) == 0 )
            return dctools::DCT_ORC_UNDEFINED_ERROR;
        else
            return dict.at(mrqRetCode);
    }
    
    
    
    inline bool MRQ_isServerAlgorithm(const int algCode)
    {
        return algCode == MRQ_BB_ALG || algCode == MRQ_LP_BB_ECP_BASED_ALG || algCode == MRQ_LP_NLP_BB_OA_BASED_ALG || algCode == MRQ_LP_BB_ESH_BASED_ALG;
    }
    
    
    
    //nTotalVarBounds should consider varaible bounds for all nodes being sent
    //this function returns the total number of bytes
    inline unsigned int MRQ_newSizeofOneOpenNodeRep(unsigned int nnodes, unsigned int nTotalVarBounds, bool considerHeader = true )
    {
        /* 1 - connect code (int32),
        *  2 - size (in bytes of the array having nodes) (uint64),
        *  3 - a number identifing the version of node representation used (uint32),
        *  4 - number of nodes (uint32),
        * 
        *  For each node, we have:
        *  		5 - lower bound (double), 
        *  		6 - nbounds (uint32), 
        *  		7 - and nbounds tuples of <index, lb, ub>
        *
        * Header includes itens 1-4
        */
        
        unsigned int sizeOfHeader = 0;
        
        if(considerHeader)
        {
            sizeOfHeader = sizeof(dctools::DCT_Int32) + sizeof(dctools::DCT_UInt64) + sizeof(dctools::DCT_UInt32) + sizeof(dctools::DCT_UInt32);
        }
        
        return  sizeOfHeader + nnodes*( sizeof(double) + sizeof(dctools::DCT_UInt32) ) + nTotalVarBounds*( sizeof(dctools::DCT_UInt32) + 2*sizeof(double) );
    }
    
    
    
    
    inline int MRQ_readNodeFromOpenNodeBufferAndShift(dctools::DCT_Byte* &buffer, double &nodelb, dctools::DCT_VarBounds &varBounds)
    {
        dctools::DCT_UInt32 nbounds = 0;
        int retCode = 0;
        
        dctools::DCT_readAndShift(buffer, nodelb);
        dctools::DCT_readAndShift(buffer, nbounds);
        
        for( decltype(nbounds) i = 0; i < nbounds; i++ )
        {
            dctools::DCT_UInt32 index;
            double lb, ub;
            dctools::DCT_Bounds bounds;
            
            dctools::DCT_readAndShift(buffer, index);
            dctools::DCT_readAndShift(buffer, lb);
            dctools::DCT_readAndShift(buffer, ub);
            
            bounds.lb = lb;
            bounds.ub = ub;
            
            try
            {
                varBounds[index] = bounds;
            }
            catch (std::bad_alloc &e)
            {
                MRQ_PRINTMEMERROR;
                retCode = dctools::DCT_RC_MEMORY_ERROR;
                goto termination;
            }
        }
        
        
    termination:
        
        return retCode;
    }
    
    
}


#endif /* MRQ_SERVER_HPP_ */
