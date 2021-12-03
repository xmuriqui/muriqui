

#include <cstdio>
#include <cstring>
#include <climits>
#include <csignal>
#include <cmath>

#include <iostream>
#include <new>
#include <exception>

#include "MRQ_server.hpp"
#include "MRQ_ampl.hpp"  //we have to include MRQ_ampl by last to avoid problems about macro expansion on ASL defines (I hate ASL)
#include "MIP_gams.hpp"



using namespace dctools;
using namespace branchAndBound;
using namespace muriqui;





MRQ_OpenNodeWriter::MRQ_OpenNodeWriter()
{
    initialize();
}

MRQ_OpenNodeWriter::~MRQ_OpenNodeWriter()
{
    desallocate();
}

void MRQ_OpenNodeWriter::desallocate()
{
    MRQ_secFree(buffer);
    bufferSize = 0;
}

void MRQ_OpenNodeWriter::initialize()
{
    bufferSize = 0;
    buffer = NULL;
}

int MRQ_OpenNodeWriter::reallocateBuffer(long unsigned int newSize)
{
    //we do not use MRQ_realloc because realloc copy bytes allocated, and we do not need this
    if(buffer)
    {
        free(buffer);
        buffer = NULL;
            bufferSize = 0;
    }
    
    if(newSize > 0)
    {
        MRQ_malloc(buffer, newSize);
        MRQ_IFMEMERRORRETURN(!buffer);
            bufferSize = newSize;
    }
    
    
    /*int r = MRQ_realloc(buffer, newSize);
    if(r != 0)
    {
        MRQ_PRINTERRORNUMBER(r);
        return r;
    } 
    
    sizeBuffer = newSize; */
    return 0;
}





//here, we receive a connection code to flag if we are sending a node by DCT_CC_OPEN_NODE_RESPONSE or DCT_CC_INITIAL_OPEN_NODE_RESPONSE
int MRQ_OpenNodeWriter::write( DCT_BaseSocket &socket, DCT_CONECTION_CODE ccode, const MRQ_NewBBNode *node, const DCT_VarBounds *rootVarBounds )
{
    /*
    * Here, we write:
    *  1 - The connection code (DCT_OPEN_NODE_RESPONSE or DCT_CC_INITIAL_OPEN_NODE_RESPONSE) code (int32)
    *  2 - The total number of bytes representing nodes being sent (uint64)
    *  3 - A number identifying the version of node representation. (This can be useful in the future, if we decide change the node representation and we sttill wish keep compatibility with old versions) (uint32)
    *  4 - The number k of nodes being sent (uint32)
    * 	For each one of the k nodes:
    *  5 - The node lower bound (double)
    *  6 - The number n of tuple of bounds transfered (uint32)
    *  7 - n tuples of < index (uint32), lb (double), ub (double)>
    *
    * 
    * we try send all data in just one call to speed up node transfer between servers and clients
    */
    
    #if DCT_DEBUG_MODE
        assert( ccode == DCT_CC_OPEN_NODE_RESPONSE || ccode == DCT_CC_INITIAL_OPEN_NODE_RESPONSE );
    #endif
    
    const DCT_Int32 int32ccode = ccode; //we are not sure if a numeration will alwasy have 32 bit in all machines...
    DCT_UInt32 numberOfNodes = 0; //number of nodes to be sent
    const DCT_UInt32 nodeRepVersion = MRQ_CURRENT_NODE_REPRESENTATION_VERSION;
    const unsigned int nRootBounds = rootVarBounds ? rootVarBounds->size() : 0;
    
    unsigned int totalBounds = 0;//(rootVarBounds ? rootVarBounds->size() : 0) + node->getNumberOfTotalBounds(); //const unsigned int totalBounds = (rootVarBounds ? rootVarBounds->size() : 0) + (parentBounds ? parentBounds->getSize() : 0) + nNodeBounds;
    
    long unsigned int totalBytes = 0;
    int r;
    unsigned int myTotalBounds = 0;
    long unsigned int bytesWritten;
    DCT_UInt64 nNodeBytes;
    DCT_Byte *pBuffer;
    
    
    
    //calculating total number of nodes to be sent and number of bounds
    for( decltype(node) pNode = node ; pNode ; pNode = (MRQ_NewBBNode*) pNode->next )
    {
        numberOfNodes++;
        totalBounds += pNode->getNumberOfTotalBounds();
    }
    
    totalBounds += (numberOfNodes*nRootBounds);
    
    totalBytes = MRQ_newSizeofOneOpenNodeRep(numberOfNodes, totalBounds);
    
    
    
    if( totalBytes > bufferSize )
    {
        int r = reallocateBuffer(2*totalBytes); //we take advantage the realloc to allocate a little more memory
        if(r != 0)
        {
            MRQ_PRINTERRORNUMBER(r);
            return DCT_RC_MEMORY_ERROR; //do not return a Muriqui code here, to do have confusion about the return codes. 
        }
    }
    
    nNodeBytes = totalBytes - sizeof(int32ccode) - sizeof(nNodeBytes); //the number of bytes being sent representing the open nodes (ie.e this excludes the connect code and this number of bytes being sent).
    
    pBuffer = buffer;
    
    DCT_writeAndShift(pBuffer, int32ccode);
    DCT_writeAndShift(pBuffer, nNodeBytes); //writting the number of bytes being sent representing the open nodes
    DCT_writeAndShift(pBuffer, nodeRepVersion);
    DCT_writeAndShift(pBuffer, numberOfNodes); //wrritting the number of nodes being sent
    
    //sending specific node information
    for( decltype(node) pNode = node ; pNode ; pNode = (MRQ_NewBBNode*) pNode->next )
    {
        double nodelb = pNode->getBestLowerBound();
        DCT_UInt32 nodeTotalBounds = nRootBounds + pNode->getNumberOfTotalBounds();
        
        DCT_writeAndShift(pBuffer, nodelb);
        DCT_writeAndShift(pBuffer, nodeTotalBounds);
        
        if(nRootBounds > 0)
        {
            for(auto &pairVBounds : *rootVarBounds)
            {
                const auto &ind = pairVBounds.first;
                const DCT_Bounds &bounds = pairVBounds.second;
                
                DCT_writeAndShift(pBuffer, ind);
                DCT_writeAndShift(pBuffer, bounds.lb);
                DCT_writeAndShift(pBuffer, bounds.ub);
                myTotalBounds++;
            }
        }
        
        bytesWritten = pNode->writeVarBoundsInaBufferArray(pBuffer);
        
        pBuffer += bytesWritten;
    }
    
    
    #if MRQ_DEBUG_MODE
        //std::cout << "pBuffer: " << (long int) pBuffer << " &buffer[totalBytes-1]: " << (long int) &buffer[totalBytes-1] << "\n";
        assert( pBuffer-1 == &buffer[totalBytes-1] );
    #endif
    
    
    
    #if 0
    //note, some indices can repeat here because they are already in rootVarBounds. But I think there is no problem about that
    if( parentBounds )
    {
        const unsigned int size = parentBounds->getSize();
        unsigned int ind;
        double l, u;
        
        for(unsigned int i = 0; i < size; i++)
        {
            parentBounds->getArrayElement(i, &ind, &l, &u);
            
            DCT_writeAndShift(pBuffer, ind);
            DCT_writeAndShift(pBuffer, l);
            DCT_writeAndShift(pBuffer, u);
            myTotalBounds++;
        }
    }
    
    //note, some indices can repeat here because they are already in parent bounds and\or rootVarBounds. But I think there is no problem about that
    for(unsigned int i = 0; i < nNodeBounds; i++)
    {
        DCT_writeAndShift(pBuffer, nodeBounds[i].ind);
        DCT_writeAndShift(pBuffer, nodeBounds[i].l);
        DCT_writeAndShift(pBuffer, nodeBounds[i].u);
        myTotalBounds++;
    }
    #endif
    
    //#if MRQ_DEBUG_MODE
        //assert(myTotalBounds == totalBounds);
        //assert( &buffer[totalBytes-1] == pBuffer-1 );
    //#endif
    
    
    
    //#if MRQ_DEBUG_MODE
        //assert( pBuffer[bytesWritten-1] == buffer[totalBytes-1] );
    //#endif
    
    
    r = socket.writeData(buffer, totalBytes);
    if(r != 0)
    {
        MRQ_PRINTERRORNUMBER(r);
        return r;
    }
    
    return 0;
}



#if 0
//here, we receive a connection code to flag if we are sending a node by DCT_CC_OPEN_NODE_RESPONSE or DCT_CC_INITIAL_OPEN_NODE_RESPONSE
int MRQ_OpenNodeWriter::write( dctools::DCT_BaseSocket &socket, dctools::DCT_CONECTION_CODE ccode, const double lb, const double *olx, const double *oux, const double *nlx, const double *nux, const unsigned int nI, const int *intVars )
{
    /*
    * Here, we write:
    *  1 - The connection code (DCT_OPEN_NODE_RESPONSE or DCT_CC_INITIAL_OPEN_NODE_RESPONSE) code (int32)
    *  2 - The total number of bytes representing nodes being sent (uint64)
    *  3 - The number k of nodes being sent (uint32) (in this case, just one node)
    * 	For each one of the k nodes:
    *  4 - The node lower bound (double)
    *  5 - The number n of tuple of bounds transfered (uint32)
    *  6 - n tuples of < index (uint32), lb (double), ub (double)>
    *
    * 
    * we try send all data in just one call to speed up node transfer between servers and clients
    */
    
    #if DCT_DEBUG_MODE
        assert( ccode == DCT_CC_OPEN_NODE_RESPONSE || ccode == DCT_CC_INITIAL_OPEN_NODE_RESPONSE );
    #endif
    
    int r;
    const DCT_UInt32 numberOfNodes = 1; //number of nodes to be sent
    const DCT_UInt32 nodeRepVersion = MRQ_CURRENT_NODE_REPRESENTATION_VERSION;
    DCT_UInt32 totalBounds = 0, auxTotalBounds = 0;
    long unsigned int totalBytes; 
    DCT_UInt64 nNodeBytes;
    DCT_Byte *pBuffer;
    
    
    
    for(unsigned int i = 0; i < nI; i++)
    {
        const auto ind = intVars[i];
        
        if( nlx[ind] != olx[ind] || nux[ind] != oux[ind] )
            totalBounds++;
    }
    
    
    totalBytes = MRQ_newSizeofOneOpenNodeRep(numberOfNodes, totalBounds); //DCT_sizeofOpenNodeRep(totalBounds);
    
    if( totalBytes > sizeBuffer )
    {
        int r = reallocateBuffer(2*totalBytes); //we take advantage the realloc to allocate a little more memory
        if(r != 0)
        {
            MRQ_PRINTERRORNUMBER(r);
            return DCT_RC_MEMORY_ERROR; //do not return a Muriqui code here, to do have confusion about the return codes. 
        }
    }
    
    nNodeBytes = totalBytes - sizeof(ccode) - sizeof(nNodeBytes); //the number of bytes being sent representing the open nodes (i.e. this excludes the connect code and this number of bytes being sent).
    
    pBuffer = buffer;
    
    DCT_writeAndShift(pBuffer, ccode);   //writing the connect code
    DCT_writeAndShift(pBuffer, nNodeBytes); //writing the number of bytes being sent representing the open nodes
    DCT_writeAndShift(pBuffer, nodeRepVersion);
    DCT_writeAndShift(pBuffer, numberOfNodes); //writing the number of nodes being sent
    DCT_writeAndShift(pBuffer, lb);
    DCT_writeAndShift(pBuffer, totalBounds);
    
    for(unsigned int i = 0; i < nI; i++)
    {
        const DCT_UInt32 ind = intVars[i];
        
        if( nlx[ind] != olx[ind] || nux[ind] != oux[ind] )
        {
            const double lbind = nlx[ind];
            const double ubind = nux[ind];
            
            DCT_writeAndShift(pBuffer, ind);
            DCT_writeAndShift(pBuffer, lbind);
            DCT_writeAndShift(pBuffer, ubind);
            
            #if MRQ_DEBUG_MODE
                auxTotalBounds++;
            #endif
        }
    }
    
    std::cout << "totalBounds: " << totalBounds << "\n";
    
    #if DCT_DEBUG_MODE
        assert(auxTotalBounds == totalBounds);
    #endif
    
    r = socket.writeData(buffer, totalBytes);
    MRQ_IFERRORRETURN(r, r);
    
    return 0;
}
#endif


//here, we receive a connection code to flag if we are sending a node by DCT_CC_OPEN_NODE_RESPONSE or DCT_CC_INITIAL_OPEN_NODE_RESPONSE
int MRQ_OpenNodeWriter::write( dctools::DCT_BaseSocket &socket, dctools::DCT_CONECTION_CODE ccode, const double lb, const double *olx, const double *oux, const double *nlx, const double *nux, const unsigned int nI, const int *intVars )
{
    const bool reallocateWithSlack = true, writeHeader = true;
    int r;
    dctools::DCT_UInt64 totalBytes; 
    
    #if MRQ_DEBUG_MODE
        totalBytes = -1;
    #endif
    
    
    r = write(buffer, bufferSize, reallocateWithSlack, writeHeader, ccode, lb, olx, oux, nlx, nux, nI, intVars, totalBytes);
    MRQ_IFERRORRETURN(r, r);
    
    
    r = socket.writeData(buffer, totalBytes);
    MRQ_IFERRORRETURN(r, r);
    
    return 0;
}



int MRQ_OpenNodeWriter::write( dctools::DCT_Byte* &output, dctools::DCT_UInt64 &maxSizeOutputArray, const bool reallocateWithSlack, const bool writeHeader, dctools::DCT_CONECTION_CODE ccode, const double lb, const double *olx, const double *oux, const double *nlx, const double *nux, const unsigned int nI, const int *intVars, dctools::DCT_UInt64 &finalSizeOutput )
{
    /*
    * Here, we write:
    *  1 - The connection code (DCT_OPEN_NODE_RESPONSE or DCT_CC_INITIAL_OPEN_NODE_RESPONSE) code (int32)
    *  2 - The total number of bytes representing nodes being sent (uint64)
    *  3 - A number identifying the version of node representation. (This can be useful in the future, if we decide change the node representation and we sttill wish keep compatibility with old versions) (uint32)
    *  4 - The number k of nodes being sent (uint32) (in this case, just one node)
    * 	For each one of the k nodes:
    *  5 - The node lower bound (double)
    *  6 - The number n of tuple of bounds transfered (uint32)
    *  7 - n tuples of < index (uint32), lb (double), ub (double)>
    *
    * 
    * we try send all data in just one call to speed up node transfer between servers and clients
    */
    
    #if MRQ_DEBUG_MODE
        if(writeHeader)
        {
            assert( ccode == DCT_CC_OPEN_NODE_RESPONSE || ccode == DCT_CC_INITIAL_OPEN_NODE_RESPONSE );
        }
    #endif
    
    int r;
    const DCT_UInt32 numberOfNodes = 1; //number of nodes to be sent
    const DCT_UInt32 nodeRepVersion = MRQ_CURRENT_NODE_REPRESENTATION_VERSION;
    DCT_UInt32 totalBounds = 0, auxTotalBounds = 0;
    long unsigned int totalBytes; 
    DCT_UInt64 nNodeBytes;
    DCT_Byte *pBuffer;
    
    
    
    finalSizeOutput = 0;
    
    for(unsigned int i = 0; i < nI; i++)
    {
        const auto ind = intVars[i];
        
        if( nlx[ind] != olx[ind] || nux[ind] != oux[ind] )
        {
            //std::cout << "ind: " << ind << " l: " << nlx[ind] << " u: " << nux[ind] << "\n";
            totalBounds++;
        }
    }
    
    
    totalBytes = MRQ_newSizeofOneOpenNodeRep(numberOfNodes, totalBounds, writeHeader);
    
    if( totalBytes > maxSizeOutputArray )
    {
        DCT_UInt64 newSize = reallocateWithSlack ? 2*totalBytes : totalBytes;
        
        MRQ_secFree(output);
        MRQ_malloc(output, newSize);
        //int r = reallocateBuffer(newSize); //we take advantage the realloc to allocate a little more memory
        if(!output)
        {
            maxSizeOutputArray = 0;
            MRQ_PRINTMEMERROR;
            return DCT_RC_MEMORY_ERROR; //do not return a Muriqui code here, to do have confusion about the return codes. 
        }
        
        maxSizeOutputArray = newSize;
    }
    
    pBuffer = output;
    
    
    if(writeHeader)
    {
        const DCT_Int32 int32ccode = ccode; //we are not sure if a numeration will alwasy have 32 bit in all machines...
        
        nNodeBytes = totalBytes - sizeof(int32ccode) - sizeof(nNodeBytes); //the number of bytes being sent representing the open nodes (i.e. this excludes the connect code and this number of bytes being sent).
        
        DCT_writeAndShift(pBuffer, int32ccode);   //writing the connect code
        DCT_writeAndShift(pBuffer, nNodeBytes); //writing the number of bytes being sent representing the open nodes
        DCT_writeAndShift(pBuffer, nodeRepVersion);
        DCT_writeAndShift(pBuffer, numberOfNodes); //writing the number of nodes being sent
    }
    
    
    DCT_writeAndShift(pBuffer, lb);
    DCT_writeAndShift(pBuffer, totalBounds);
    
    for(unsigned int i = 0; i < nI; i++)
    {
        const DCT_UInt32 ind = intVars[i];
        
        if( nlx[ind] != olx[ind] || nux[ind] != oux[ind] )
        {
            const double lbind = nlx[ind];
            const double ubind = nux[ind];
            
            DCT_writeAndShift(pBuffer, ind);
            DCT_writeAndShift(pBuffer, lbind);
            DCT_writeAndShift(pBuffer, ubind);
            
            #if MRQ_DEBUG_MODE
                auxTotalBounds++;
            #endif
        }
    }
    
    //std::cout << "totalBounds: " << totalBounds << "\n";
    
    #if DCT_DEBUG_MODE
        assert(auxTotalBounds == totalBounds);
        assert(output + totalBytes == pBuffer);
    #endif
        
    finalSizeOutput = totalBytes;
    
    
    //r = socket.writeData(buffer, totalBytes);
    //MRQ_IFERRORRETURN(r, r);
    
    return 0;
}



int MRQ_OpenNodeWriter::write( dctools::DCT_BaseSocket &socket, dctools::DCT_CONECTION_CODE ccode, const MRQ_OpenNodesStorer &openNodesStorer)
{
    /*
    * Here, we write:
    *  1 - The connection code (DCT_OPEN_NODE_RESPONSE or DCT_CC_INITIAL_OPEN_NODE_RESPONSE) code (int32)
    *  2 - The total number of bytes representing nodes being sent (uint64)
    *  3 - A number identifying the version of node representation. (This can be useful in the future, if we decide change the node representation and we sttill wish keep compatibility with old versions) (uint32)
    *  4 - The number k of nodes being sent (uint32) (in this case, just one node)
    * 	For each one of the k nodes:
    *  5 - The node lower bound (double)
    *  6 - The number n of tuple of bounds transfered (uint32)
    *  7 - n tuples of < index (uint32), lb (double), ub (double)>
    *
    * 
    * we try send all data in just one call to speed up node transfer between servers and client
    */
    
    #if MRQ_DEBUG_MODE
        assert( ccode == DCT_CC_OPEN_NODE_RESPONSE || ccode == DCT_CC_INITIAL_OPEN_NODE_RESPONSE );
    #endif
    
    int r;
    const DCT_Int32 int32ccode = ccode;
    const DCT_UInt32 numberOfNodes = openNodesStorer.nNodes; //number of nodes to be sent
    const DCT_UInt32 nodeRepVersion = MRQ_CURRENT_NODE_REPRESENTATION_VERSION;
    long unsigned int totalBytes; 
    DCT_UInt64 nNodeBytes = 0;
    
    
    const auto nBytesInNodeRep = openNodesStorer.nBytesInNodeRep;
    DCT_Byte ** openNodesRep = openNodesStorer.nodesRep; 
    
    DCT_Byte *pBuffer;
    
    
    for( DCT_UInt32 i = 0; i < numberOfNodes; i++ )
    {
        nNodeBytes += nBytesInNodeRep[i];
    }
    
    nNodeBytes += sizeof(nodeRepVersion) + sizeof(numberOfNodes);
    
    totalBytes = nNodeBytes + sizeof(int32ccode) + sizeof(nNodeBytes);
    
    if( totalBytes > bufferSize )
    {
        int r = reallocateBuffer(totalBytes);
        MRQ_IFERRORRETURN(r, r);
    }
    
    pBuffer = buffer;
    
    DCT_writeAndShift(pBuffer, int32ccode); //writing the connect code
    DCT_writeAndShift(pBuffer, nNodeBytes); //writing the number of bytes being sent representing the open nodes
    DCT_writeAndShift(pBuffer, nodeRepVersion);//writing the node representation version
    DCT_writeAndShift(pBuffer, numberOfNodes); //writing the number of nodes being sent
    
    //writting the nodes
    for( DCT_UInt32 i = 0; i < numberOfNodes; i++ )
    {
        auto nBytesForThisNode = nBytesInNodeRep[i];
        
        MRQ_copyArrayAnySize( nBytesForThisNode, openNodesRep[i], pBuffer);
        
        pBuffer += nBytesForThisNode;
    }
    
    #if MRQ_DEBUG_MODE
        assert( pBuffer == buffer + totalBytes ); //making sure we set the correct number of butes
    #endif
    
    r = socket.writeData(buffer, totalBytes);
    MRQ_IFERRORRETURN(r, r);
    
    return 0;
}




MRQ_OpenNodesStorer::MRQ_OpenNodesStorer()
{
    initialize();
}


int MRQ_OpenNodesStorer::allocate(unsigned int maxNodes)
{
    deallocateAll();
    
    MRQ_calloc(nodesRep, maxNodes);
    MRQ_calloc(nBytesInNodeRep, maxNodes);
    MRQ_IFMEMERRORRETURN(!nodesRep || !nBytesInNodeRep);
    
    this->maxNodes = maxNodes;
    
    return 0;
}


void MRQ_OpenNodesStorer::deallocateAll()
{
    if(nodesRep)
    {
        deallocateNodes();
        
        free(nodesRep);
        nodesRep = NULL;
    }
    
    MRQ_secFree(nBytesInNodeRep);
    
    maxNodes = 0;
}


//here, we just deallocate nodes representation, but we stiil keep nodesRep
void MRQ_OpenNodesStorer::deallocateNodes()
{
    if( nodesRep )
    {
        for(decltype(maxNodes) i = 0; i < maxNodes; i++)
        {
            if( nodesRep[i] )
                free(nodesRep[i]);
        }
    }
    
    nNodes = 0;
}


void MRQ_OpenNodesStorer::initialize()
{
    nNodes = 0;
    maxNodes = 0;
    nodesRep = NULL;
    nBytesInNodeRep = NULL;
}


MRQ_OpenNodesStorer::~MRQ_OpenNodesStorer()
{
    deallocateAll();
}



MRQ_ServerBBCallbacks::MRQ_ServerBBCallbacks(DCT_ServerConnection  *connection) : MRQ_UserCallbacks()
{
    this->connection = connection;
    rootOpenNodeRep = NULL;
    intVars = NULL;
    auxSentSol = NULL;
    auxVars = NULL;
    nodeBounds = NULL;
    receivedObjAndSol = NULL;
    requestFirstBestSol = true;
    setLazyConstraintsOnReceivedSol = false; //cplex does not work with lazy constraints on solve callback...
    stopAlg = false;
    stopService = false;
    nI = 0;
    nRequestedNodes = 0;
    myiter = 0;
    minNumberOfNodestoSendNode = 2;
    maxNumberOfNodesToSend = 1000;
    
    myzl = -INFINITY;
    frequencyToSendLowerBound = 20000;
    minRelativeGaptoUsePseudoCosts = 0.0; //so, we will use pseudocost always
}


MRQ_ServerBBCallbacks::~MRQ_ServerBBCallbacks()
{
    deallocate();
}


int MRQ_ServerBBCallbacks::allocateStructures(MRQ_Algorithm *alg, MRQ_MINLPProb *prob)
{
    const unsigned int nvars = prob->n;
    
    
    this->prob = prob;
    
    nI = prob->getNumberOfIntegerVars();
    
    MRQ_malloc(intVars, nI);
    MRQ_malloc(auxSentSol, nvars + 3);
    MRQ_malloc(auxVars, 2*(nvars+1) );
    MRQ_malloc(nodeBounds, 2*(nvars+1) );
    MRQ_malloc(receivedObjAndSol, nvars+1);
    
    MRQ_IFMEMERRORRETURN(!intVars || !auxSentSol || !auxVars || !nodeBounds || !receivedObjAndSol);
    
    prob->getIntegerIndices(intVars);
    
    return 0;
}


void MRQ_ServerBBCallbacks::deallocate()
{
    MRQ_secFree(intVars);
    MRQ_secFree(auxSentSol);
    MRQ_secFree(auxVars);
    MRQ_secFree(nodeBounds);
    MRQ_secFree(receivedObjAndSol);
}


int MRQ_ServerBBCallbacks::beforeAll(const MRQ_ALG_CODE algCode, const unsigned int numberOfThreads)
{
    //Note, beforeAll can be called by others algorithms insinde B&B code, like, for example, some heuristics. So, we need to restrict this lock to B&B
    
    if( MRQ_isServerAlgorithm(algCode) ) 
    {
        stopAlg = false;
        nthreads = numberOfThreads;
        
        
        if( algCode == MRQ_BB_ALG ) 
        {
            nextIterToSendLowerBound = 2;
            
            if(nRequestedNodes > 0)
            {
                SEMAPH_startOthersThreadsAfterSendingNodesToOtherServes.lock(nthreads); //we lock this semaphore to allow just thread 0 start B&B exploration. Ony after thread 0 sent initial nodes to other servers, we let other trhreads start exploration. // Note, if nRequestedNodes is 0, so, this server should not send any initial node to other servers, and we can skip this lock.
            }
        }
        else
        {
            nextIterToSendLowerBound = 1000; //to avoid problems where we have several open nodes to run, we put 1000. So, we just send lower bound if a node requires more than 1000 iterations. In this way, we avoid another server request to this server an open node and never get response, since this server is exploring several open nodes and any node do not requires more than 1 or two iterations...
            
            myiter = 0;
            myzl = -INFINITY;
            
            //WE ASSUME THAT IS A LINEAR APPROXIMATION ALGORITHM
            int r = openNodesStorer.allocate( maxNumberOfNodesToSend );
            MRQ_IFERRORRETURN(r, r);
            
        }
        
        accNodesToSend = false;
    }
    
    return 0;
}


//by now, we prefer set the bounds in the original problem to reduce the size of open nodes. Also, it would be needed change bb code to call BB_generateRootNode...
#if 0
/*we generate our own root node to put the varBounds received from client.
* It is better do it because if we change the bounds in the original problem,
* we can get some inconsistencies when a server pass a node to another
*/ 
int MRQ_ServerBBCallbacks::BB_generateRootNode( MRQ_NewBBNode* &rootNode )
{
    const unsigned int nbounds = rootVarBounds->size();
    unsigned int k = 0;
    int code, r;
    BBL_NodeBounds *nodeBounds;
    
    MRQ_getchar();
    
    if(nbounds == 0)
    {
        std::cout << "budega!\n";
        MRQ_getchar();
        code = 0;
        goto termination;
    }
    
    
    rootNode = new (std::nothrow) MRQ_NewBBNode;
    
    if(!rootNode)
    {
        MRQ_PRINTMEMERROR;
        code = MRQ_MEMORY_ERROR;
        goto termination;
    }
    
    
    rootNode->parentBounds = new (std::nothrow) BBL_ArraySize <BBL_NodeBounds>;
    
    if(!rootNode->parentBounds)
    {
        MRQ_PRINTMEMERROR;
        code = MRQ_MEMORY_ERROR;
        goto termination;
    }
    
    
    r = rootNode->parentBounds->allocate(nbounds);
    if(r != 0)
    {
        MRQ_PRINTERRORNUMBER(r);
        code = MRQ_MEMORY_ERROR;
        goto termination;
    }
    
    
    nodeBounds = rootNode->parentBounds->a;
    
    for( auto &pairVBounds : *rootVarBounds )
    {
        const DCT_Bounds &bounds = pairVBounds.second;
        
        nodeBounds[k].ind = pairVBounds.first;
        nodeBounds[k].l = bounds.lb;
        nodeBounds[k].u = bounds.ub;
        
        k++;
    }
    
    #if MRQ_DEBUG_MODE
        assert(k == nbounds);
    #endif
    
        std::cout << "no recebido\n";
    rootNode->print();
    MRQ_getchar();
    
    
    code = 0;
    
termination:
    
    if(code != 0)
        MRQ_secDelete(rootNode);
    
    
    return code;
}
#endif



int MRQ_ServerBBCallbacks::BB_generateRootNode( MRQ_NewBBNode* &rootNode )
{
    int retCode;
    unsigned int n = prob->n;
    const branchAndBound::BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY parentNodeBoundsStrategy = (branchAndBound::BBL_PARENT_NODE_BOUNDS_STORAGE_STRATEGY) ( (MRQ_BranchAndBound*) myAlg)->in_parent_node_bounds_storage_strategy;
    
    DCT_UInt32 numberOfNodes;
    DCT_Byte *paux = rootOpenNodeRep;
    double *lx = prob->lx; // *ux = prob->ux;
    MRQ_NewBBNode *p;
    DCT_VarBounds nodeVarsBounds;
    
    
    rootNode = NULL; //initializing rootNode. It is important for the case if we have some meory error, we have to delete nodes allocated
    
    DCT_readAndShift(paux, numberOfNodes);
    
    #if MRQ_DEBUG_MODE
        assert(numberOfNodes > 1);
        assert( ((MRQ_BranchAndBound*) myAlg)->in_branching_strategy != MRQ_BB_BS_STBRANCH_PSEUDO_COSTS); //we cannot use pseudo-costs here because we have several disconected root nodes
    #endif
    
    
    for(decltype(numberOfNodes) k = 0; k < numberOfNodes; k++)
    {
        double nodelb;
        unsigned int j = 0;
        
        //we are readinge the node variable bounds to a map (DCT_VarBounds) because sender server can specify two or more time bounds to the same variable sine node bounds can be set directy in minlp problem...
        int r = MRQ_readNodeFromOpenNodeBufferAndShift(paux, nodelb, nodeVarsBounds);
        MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        
        MRQ_NewBBNode *node = new (std::nothrow) MRQ_NewBBNode( parentNodeBoundsStrategy, n);
        MRQ_IFMEMERRORGOTOLABEL(!node, retCode, termination);
        
        
        if( k == 0 )
        {
            rootNode = node;
        }
        else
        {
            node->previous = p;
            p->next = node;
        }
        
        p = node;
        
        #if MRQ_SET_LBHEUR_ON_BB_NODE
            node->heurlb = nodelb;
        #endif
        
        
        r = node->allocateNodeBounds( nodeVarsBounds.size() );
        
        //std::cout << "###################################################################\n";
        //std::cout << "Node: " << k << "\n";
        for(auto &pairvbounds : nodeVarsBounds)
        {
            const auto &ind = pairvbounds.first;
            const DCT_Bounds &bounds = pairvbounds.second;
            double vlb = bounds.lb, vub = bounds.ub;
            double auxSolValue; //auxSolValue is just to simulate the optimal solution of parent node, since we do not know this solution. We will simulate a solution value consider a gap of 0.5
            
            if( lx[ind] == vlb )
                auxSolValue = vub + 0.5;
            else
                auxSolValue = vlb - 0.5;
            
            
            int r = node->myBounds.setArrayElement(j, &ind, &vlb, &vub, &auxSolValue);
            MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            //std::cout << "ind: " << ind << " l: " << bounds.lb << " u: " <<  bounds.ub << "\t";
            
            j++;
        }
        
        //std::cout << "\n";
        
        //node->print();
        
        nodeVarsBounds.clear(); //we need a new nodeVarsBounds for each iteration
    }
    
    //MRQ_getchar();
    
    retCode = 0;
    
termination:
    
    if(retCode != 0) 
    { //deleting nodes already allocated
        MRQ_NewBBNode *p1, *p2;
        
        for(p1 = rootNode ; p1; p1 = p2)
        {
            p2 = (decltype(p2)) p1->next;
            delete p1;
        }
    }
    
    return retCode;
}



int MRQ_ServerBBCallbacks::beforeBBLoop(const unsigned int threadNumber, const double lb, const double ub, MRQ_NLPSolver &nlpSolver)
{
    //std::cout << DCT_PREPRINT "Thread " << threadNumber << " entered at MRQ_ServerBBCallbacks::beforeBBLoop \n";
    //DCT_getchar();
    
    if( threadNumber != 0 && nRequestedNodes > 0 )
    {
        SEMAPH_startOthersThreadsAfterSendingNodesToOtherServes.lock(nthreads); //thread 0 locked this semaphore in MRQ_ServerBBCallbacks::beforeAll. So, all other threads will be blocked here.  We lock this semaphore to allow just thread 0 start B&B exploration. Ony after thread 0 sent initial nodes to other servers, we let other trhreads start exploration.
        SEMAPH_startOthersThreadsAfterSendingNodesToOtherServes.unlock(nthreads); //open the semaphore to other threads
    }
    
    return 0;
}



void MRQ_ServerBBCallbacks::afterBBLoop(const unsigned int threadNumber, const double lb, const double ub, MRQ_NLPSolver &nlpSolver, const int threadReturnCode)
{
    //std::cout << DCT_PREPRINT "Thread " << threadNumber << " entered at MRQ_ServerBBCallbacks::afterBBLoop \n";
    //DCT_getchar();
    
    if( threadNumber == 0 && nRequestedNodes > 0 )
        SEMAPH_startOthersThreadsAfterSendingNodesToOtherServes.unlock(nthreads); // we unlock just to guarentee all threads will be unlock if some error coccours befoe thread 0 sent nodes to other servers.
}



int MRQ_ServerBBCallbacks:: BB_beforeSolvingRelaxation( const unsigned int threadNumber, MRQ_NewBBNode &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, MRQ_NLPSolver &nlpSolver, bool &pruneNode)
{
    bool solReceieved, receivedSolUpdtBestSol;
    
    int r = general_beforeSolvingRelaxation(threadNumber, MRQ_BB_ALG, NULL, &node, iter, lb, node.getBestLowerBound(), nlx, nux, solReceieved, receivedSolUpdtBestSol, receivedObjAndSol, pruneNode);
    MRQ_IFERRORRETURN(r, r);
    
    return 0;
}



int MRQ_ServerBBCallbacks::general_beforeSolvingRelaxation(const unsigned int threadNumber, const MRQ_ALG_CODE algCode, MRQ_MILPSolverCallbackInterface *milpSolverCallbackInterface, MRQ_NewBBNode *node, const long unsigned int iter, const double lb, const double nodelb, double *nlx, double *nux, bool &solReceived, bool &receivedSolUpdtBestSol, double *receivedObjAndSol, bool &pruneNode)
{
    bool linAppAlg = milpSolverCallbackInterface != NULL;
    const unsigned int indThreadReadConnection = algCode == MRQ_BB_ALG ? nthreads/2u : 0;
    
    
    #if MRQ_DEBUG_MODE
        if(linAppAlg)
        {
            assert(node == NULL);
            assert(nlx == NULL);
            assert(nux == NULL);
        }
        else
        {
            assert(node != NULL);
            assert(nlx != NULL);
            assert(nux != NULL);
        }
    #endif
    
    solReceived = false;
    receivedSolUpdtBestSol = false;
    pruneNode = false;
    
    
    //before threat other connection codes, first we send the initial open nodes. We do it because otherwise, we would have to set SEMAPH_connectionWrite to threat the connection codes received and we are avoiding it since in some moment, the server can receive several and several requestiong for open nodes and setting semaphore all time spending computation
    if( nRequestedNodes > 0 )
    {
        int r = 0;
        long unsigned int nOpenNodes;
        
        if( milpSolverCallbackInterface )
        {
            r = milpSolverCallbackInterface->getNumberOfOpenNodes(nOpenNodes);
            MRQ_IFERRORRETURN(r, r);
            
            //r = milpSolverCallbackInterface->getCurrentNodeLowerBound(nodelb);
            //MRQ_IFERRORRETURN(r, r);
            
        }
        else
        {
            r = BB_getNumberOfOpenNodes(nOpenNodes);
            MRQ_IFERRORRETURN(r, r);
            
            //nodelb = node->getBestLowerBound();
        }
        
        
        SEMAPH_connectionWrite.lock(nthreads); //we should set some semaphore here. It could be another semaphore, but we take advantage of SEMAPH_connectionWrite because we would have to lock it anyway to send the open node.
        {
            if(nRequestedNodes > 0) //here, we check again because the correct is check nRequestedNodes inside exclusion zone. We just put the first if to avoid set the semaphore unnecessarily when nRequestedNodes is zero.
            {
                //std::cout << "nOpenNodes: " << nOpenNodes << " nRequestedNodes: " << nRequestedNodes << "\n";
                //MRQ_getchar();
                
                if( nOpenNodes >= nRequestedNodes )
                { //note we have, at least nRequestedNodes + 1 to exploring(coutning with the node in this thread stored at variable node). So, we send the current node
                    
                    if(linAppAlg)
                    {
                        /*we put the nlx and nux calculation inside muatl exclucusion zone because auxVars é shared between all threads*/
                        
                        nlx = auxVars;
                        nux = &auxVars[n+1];
                        
                        r = milpSolverCallbackInterface->getVarsLowerBoundsOnNode(n, nlx);
                        MRQ_IFERRORRETURN(r, r);
                    
                        r = milpSolverCallbackInterface->getVarsUpperBoundsOnNode(n, nux);
                        MRQ_IFERRORRETURN(r, r);
                        
                        
                        r = nodeWriter.write( *connection, DCT_CC_INITIAL_OPEN_NODE_RESPONSE, nodelb, olx, oux, nlx, nux, nI, intVars );
                    }
                    else
                    {
                        r = nodeWriter.write( *connection, DCT_CC_INITIAL_OPEN_NODE_RESPONSE, node, rootVarBounds );
                    }
                    
                    //MRQ_getchar();
                    
                    if(r == 0)
                    {
                        nRequestedNodes--;
                        pruneNode = true; //we are sending this node to client, so, it does not belong us
                    }
                }
            }
        }
        SEMAPH_connectionWrite.unlock(nthreads);
        
        
        
        if( algCode == MRQ_BB_ALG && nRequestedNodes <= 0 )
        { //we sent all initial nodes to other servers. So, we can unlock other threads to they could help in the B&B exploration
            SEMAPH_startOthersThreadsAfterSendingNodesToOtherServes.unlock(nthreads);
        }
        
        
        if(r != 0)
        {
            MRQ_PRINTERRORNUMBER(r);
            return r;// //here, we do not threat error
        }
        
        
        if( pruneNode )
        {
            std::cout << MRQ_PREPRINT "Sending open node\n";
            
            if(linAppAlg)
            {
                for(unsigned int i = 0; i < nI; i++)
                {
                    const auto ind  = intVars[i];
                    
                    if( olx[ind] != nlx[ind] || oux[ind] != nux[ind] )
                    {
                        const auto l = nlx[ind];
                        const auto u = nux[ind];
                        
                        std::cout << "var: " << ind << " l: " << l << " u: " << u << "\t";
                    }
                }
                std::cout << "\n";
            }
            else
            {
                node->print();
            }
            
            //MRQ_getchar();
        }
    }
    //we just put this test to only allow a specific thread be responsible for all comunication with the client. So, other threads will not be blocked by the mutex and could explore nodes normally.
    else if(threadNumber == indThreadReadConnection)
    {
        //printf("[");
        //fflush(stdout);
        
        bool sendLowerBound = iter > nextIterToSendLowerBound;
        bool sendNoNodeResp = false;
        bool sendNode = false;
        bool requestBestSol = requestFirstBestSol;
        
        DCT_Int32 request;
        
        int readReturnCode = -1, rr = 0;
        
        long unsigned int nOpenNodes;
        
        
        if(stopAlg)
        {
            MRQ_PRINTERRORMSG("An error was gotten. Giving up.");
            return DCT_RC_UNDEFINED_ERROR;
        }
        
        
        if(sendLowerBound)
            nextIterToSendLowerBound += frequencyToSendLowerBound; //do not move this from here
        
        
        //connection here is non blocking to read operations. So, it is possible we have nothing to read. In this case, we do not need treat anything .
        //since now, just a single thread can execute this block of code, we only need mutex to write operations (any thread can still send a new best solution).
        
        
        if( accNodesToSend )
        {
            int r; // do not use the variable r delared outside
            
            /* that is for linear approximation algorithms. Since client request open nodes, we do not readsocket until we collect all opne nodes to send to client.
            */
            
            #if MRQ_DEBUG_MODE
                //std::cout << "openNodesStorer.nNodes: " << openNodesStorer.nNodes << " openNodesStorer.maxNodes: " << openNodesStorer.maxNodes << "\n";
                assert(openNodesStorer.nNodes < openNodesStorer.maxNodes ); //we are storing a new node. So, we have to have space for it!
            #endif
            
            DCT_Byte* &nodeRep = openNodesStorer.nodesRep[ openNodesStorer.nNodes ];
            DCT_UInt64 &sizeNodeRep = openNodesStorer.nBytesInNodeRep[ openNodesStorer.nNodes ];
            
            
            nlx = nodeBounds; //only this thread use that 
            nux = &nlx[n+1];
                
            r = milpSolverCallbackInterface->getVarsLowerBoundsOnNode(n, nlx);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            r = milpSolverCallbackInterface->getVarsUpperBoundsOnNode(n, nux);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            
            
            int rr = nodeWriter.write( nodeRep, sizeNodeRep, false, false, DCT_CC_OPEN_NODE_RESPONSE, nodelb, olx, oux, nlx, nux, nI, intVars, sizeNodeRep );
            MRQ_IFERRORRETURN(rr, rr);
            
            openNodesStorer.nNodes++;
            
            
            if( openNodesStorer.nNodes == maxNumberOfNodesToSend ) 
            {
                sendNode = true;
            }
            else
            {
                #if MRQ_DEBUG_MODE
                    assert(openNodesStorer.nNodes <= maxNumberOfNodesToSend); //number of nodes stored cannot be greather than maxNumberOfNodesToSend
                #endif
            }
            
            pruneNode = true;
            
            
            /*std::cout << MRQ_PREPRINT "Storing open node:\n";
            
            for(unsigned int i = 0; i < nI; i++)
            {
                const auto ind  = intVars[i];
                
                if( olx[ind] != nlx[ind] || oux[ind] != nux[ind] )
                {
                    const auto l = nlx[ind];
                    const auto u = nux[ind];
                    
                    std::cout << "var: " << ind << " l: " << l << " u: " << u << "\t";
                }
            }
            std::cout << "\n"; */
            
            
            //pruneNode = true; //we are sending this node to client, so, it does not belong us
        }
        
        
    
        //SEMAPH_connection.lock(nthreads);
        {
            if( !accNodesToSend ) 
            {
                readReturnCode = connection->readData(&request, sizeof(request));
                
                if(readReturnCode == 0)
                {
                    //data were read with success
                    if( request == DCT_CC_LOWER_BOUND_REQUEST  )
                    {
                        sendLowerBound = true;
                    }
                    else if( request == DCT_CC_OPEN_NODE_REQUEST )
                    {
                        /* Here, if we have more than one node to exploit (beyond this current node), we send this current node to client. */ 
                        
                        int r;
                        
                        if(linAppAlg)
                            r = milpSolverCallbackInterface->getNumberOfOpenNodes( nOpenNodes);
                        else
                            r = BB_getNumberOfOpenNodes(nOpenNodes);
                            
                        MRQ_IFERRORRETURN(r, r);
                        
                        
                        //std::cout << "Recebi pedido de no em aberto. ";
                        
                        if( nOpenNodes < minNumberOfNodestoSendNode )
                        {
                            sendNoNodeResp = true;
                        }
                        else
                        {
                            if(linAppAlg)
                                accNodesToSend = true;
                            else
                                sendNode = true;
                        }
                        
                    }
                    else if( request == DCT_CC_NEW_SOLUTION_FOUND )
                    {
                        /*
                        * Here, we read:
                        * 		1 - The number of variables (uint32)
                        * 		2 - The objective value (double)
                        * 		3 - solution (n doubles)
                        * 
                        */
                        
                        DCT_UInt32 nvars;
                        
                        //here, we need read more data, so, we have to change the reading proccess to block, otherwise, maybe we read before data arrives.
                        
                        connection->enableBlock();
                        
                        rr = connection->readData( &nvars, sizeof(nvars) );
                        if(rr != 0)
                        {
                            MRQ_PRINTERRORNUMBER(rr);
                        }
                        
                        rr = connection->readData( receivedObjAndSol, (n+1)*sizeof(double) ); //we just read the solution here. The updation, we perform outside this exclusion zone
                        if(rr != 0)
                        {
                            MRQ_PRINTERRORNUMBER(rr);
                        }
                        
                        solReceived = true;
                        
                        //now, we restore the nonblocking proccess
                        connection->disableBlock();
                        
                    }
                    else if( request == DCT_CC_NO_SOLUTION_AVAILABLE)
                    { //we do nothing about this
                    }
                    else if( request == DCT_CC_STOP_SERVICE )
                    {
                        stopService = true;
                        return DCT_SCC_REQUEST_BY_CLIENT;
                    }
                    else
                    {
                        MRQ_PRINTERRORMSGP("Unknow connect code recieved from client: ", request);
                        MRQ_getchar();
                    }
                }
                
                
            }
            
            
            
            if( sendLowerBound || sendNoNodeResp || sendNode || requestBestSol )
            {
                int rr2 = rr = 0;
                double zl = lb;
                //double nlb;
                
                unsigned int nNodesToSend = 0; 
                MRQ_NewBBNode *otherNodesToSend = NULL; //we set this pointer only for BB alg to set other nodes from open node list to send to other servers
                
                
                if( linAppAlg )
                { 
                    
                    if( sendLowerBound )
                    { //if we are in linAppAlg, lb variable should be set as infinity. So, we read zl from milpSolverCallbackInterface
                        int r = milpSolverCallbackInterface->getMIPDualBound(zl); 
                        MRQ_IFERRORRETURN(r, r);
                        
                        if( zl > myzl )
                            myzl = zl;
                    }
                    
                }
                
                
                SEMAPH_connectionWrite.lock(nthreads);
                {
                    if(sendLowerBound)
                    {
                        const DCT_UInt32 sizelbresponse = sizeof(DCT_Int32) + sizeof(double);
                        DCT_Int32 lbresponse[5];
                        double *pd = (double*) &lbresponse[1];
                        
                        /*
                        * Here, we write: 
                        * 	1 -  DCT_LOWER_BOUND_RESPONSE code (32 bit)
                        * 	2 -  lower bound (double)
                        */
                        
                        //we build lb response out of exclusion region...
                        lbresponse[0] = DCT_CC_LOWER_BOUND_RESPONSE;
                        *pd = zl;
                        
                        rr2 = connection->writeData( lbresponse, sizelbresponse );
                        
                        if(rr2 != 0)
                        {
                            MRQ_PRINTERRORNUMBER(rr2);
                            //here, we do not threat error
                        }
                        
                        std::cout << MRQ_PREPRINT "sending lower bound to client: " << zl << "\n";
                    }
                    
                    if( sendNoNodeResp )
                    {
                        const DCT_Int32 noNodeResp = DCT_CC_NO_OPEN_NODE_AVAILABLE;
                        
                        rr = connection->writeData( &noNodeResp, sizeof(noNodeResp) );
                        if(rr != 0)
                        {
                            MRQ_PRINTERRORNUMBER(rr);
                            //here, we do not threat error
                        }
                        
                        std::cout << MRQ_PREPRINT "no open node available!\n";
                    }
                    else if( sendNode )
                    {
                        
                        if( linAppAlg )
                        {
                            int r;
                            
                            /*//we put nlx and nux calculations inside mutual exclusion zone because all threads share auxVars
                            nlx = auxVars;
                            nux = &auxVars[n+1];
                                
                            r = milpSolverCallbackInterface->getVarsLowerBoundsOnNode(n, nlx);
                            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
                            
                            r = milpSolverCallbackInterface->getVarsUpperBoundsOnNode(n, nux);
                            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR); */
                            
                            
                            //rr = nodeWriter.write( *connection, DCT_CC_OPEN_NODE_RESPONSE, nodelb, olx, oux, nlx, nux, nI, intVars );
                            
                            rr = nodeWriter.write( *connection, DCT_CC_OPEN_NODE_RESPONSE, openNodesStorer );
                            
                            accNodesToSend = false;
                            openNodesStorer.nNodes = 0;
                            
                            //MRQ_getchar();
                        }
                        else
                        {
                            #if MRQ_DEBUG_MODE
                                assert( algCode == MRQ_BB_ALG );
                            #endif
                            
                            nNodesToSend =  MRQ_min<long unsigned int>(nOpenNodes/2, maxNumberOfNodesToSend);
                            
                            
                            
                            if( nNodesToSend <= 0)
                                nNodesToSend = 1;
                            
                            
                            nNodesToSend = BB_getOpenNodes(nNodesToSend, otherNodesToSend, true);
                            
                            if(nNodesToSend == 0)
                            {
                                MRQ_PRINTERRORMSG("Error: no BB node to send in the open list. Sending current node");
                                otherNodesToSend = node;
                                nNodesToSend = 1; 
                                pruneNode = true; //we are sending this node to client, so, it does not belong us
                            }
                            //nodesToSend = node;
                            
                            rr = nodeWriter.write( *connection, DCT_CC_OPEN_NODE_RESPONSE, otherNodesToSend, rootVarBounds );
                            if(rr != 0)
                            {
                                MRQ_PRINTERRORNUMBER(rr);
                                //here, we do not threat error
                            }
                            
                            
                            {
                                std::cout << MRQ_PREPRINT "Sending " << nNodesToSend << " open nodes:\n";
                                
                                /*unsigned int j = 0;
                                for( auto p = otherNodesToSend; p; p = (decltype(p)) p->next )
                                {
                                    std::cout << "###################################################################\n" 
                                        "Node " << j << ":\n";
                                    p->print();
                                    j++;
                                }
                                
                                MRQ_getchar(); */
                            }
                            
                        }
                        
                        if(rr != 0)
                        {
                            MRQ_PRINTERRORNUMBER(rr);
                            //here, we do not threat error
                        }
                        
                        
                        
                        /*if(linAppAlg)
                        {
                            std::cout << MRQ_PREPRINT "Sending open node:\n";
                            
                            for(unsigned int i = 0; i < nI; i++)
                            {
                                const auto ind  = intVars[i];
                                
                                if( olx[ind] != nlx[ind] || oux[ind] != nux[ind] )
                                {
                                    const auto l = nlx[ind];
                                    const auto u = nux[ind];
                                    
                                    std::cout << "var: " << ind << " l: " << l << " u: " << u << "\t";
                                }
                            }
                            std::cout << "\n";
                            
                            
                            //pruneNode = true; //we are sending this node to client, so, it does not belong us
                        } */
                        
                        
                        //MRQ_getchar();
                        
                        
                    }
                    
                    
                    if(requestBestSol)
                    { //it is possible some other server found a solution before this server starts. So, we request the solution
                        const DCT_Int32 ccode = DCT_CC_BEST_SOLUTION_REQUEST;
                        
                        int r = connection->writeData(&ccode, sizeof(ccode));
                        if(r != 0)
                        {
                            //by now, we do not care if we got some error
                            MRQ_PRINTERRORNUMBER(r);
                        }
                        
                        requestFirstBestSol = false;
                    }
                    
                    
                }SEMAPH_connectionWrite.unlock(nthreads);
                
                
                if( otherNodesToSend && otherNodesToSend != node )
                { //deleting nodes sent
                    MRQ_NewBBNode *p1, *p2;
                    
                    #if MRQ_DEBUG_MODE
                        unsigned int k = 0;
                    #endif
                    
                    for(p1 = otherNodesToSend; p1; p1 = p2)
                    {
                        p2 = (decltype(p2)) p1->next;
                        delete p1;
                        #if MRQ_DEBUG_MODE
                            k++;
                        #endif
                    }
                    
                    #if MRQ_DEBUG_MODE
                        //std::cout << "nNodesToSend: " << nNodesToSend << " k: " << k << "\n";
                        assert(k == nNodesToSend);
                    #endif
                }
                
                
                
                if(rr2 == DCT_RC_CLOSED_CONNECTION || rr == DCT_RC_CLOSED_CONNECTION)
                { //connection with the client was closed. So, we finish the execution
                    return DCT_RC_CLOSED_CONNECTION;
                }
                
                
            }
            
            
        }
        //SEMAPH_connection.unlock(nthreads);
    
        
        //if r == 0, a data was read from socket with success
        if(readReturnCode == 0)
        {
            if(request == DCT_CC_NEW_SOLUTION_FOUND)
            {
                if(rr == 0)
                {
                    MRQ_PRINTMSG("Trying updating solution received from client.\n");
                    receivedSolUpdtBestSol = tryUpdateBestSolution(threadNumber, n, &receivedObjAndSol[1], receivedObjAndSol[0], iter);
                }
                else
                {
                    MRQ_PRINTERRORMSG("Error at reading solution from socket");
                    //MRQ_getchar();
                }
            }
            /*else if(request == DCT_STOP_SERVICE)
            {
                //we just return a nonzero code to stop the branch-and-bound
                return DCT_STOP_SERVICE;
            }*/
            
        }
        
        //printf("] ");
        //fflush(stdout);
    }
    
    return 0;
}



int MRQ_ServerBBCallbacks::linearApp_beforeSolveInMILPBB(const MRQ_ALG_CODE algCode, const unsigned int threadNumber, MRQ_MILPSolverCallbackInterface &milpSolverCallbackInterface, MRQ_SolutionStorer &solsToBuildLazyConstraints )
{
    bool pruneNode = false;
    int r;
    double nodelb;
    const double zu = getUpperBound();
    
    
    r = milpSolverCallbackInterface.getCurrentNodeLowerBound(nodelb);
    MRQ_IFERRORRETURN(r, r);
    
    //std::cout << "nodelb: " << nodelb << " zu: " << zu << " ";
    
    if( nodelb >= zu ) //TODO: consider here absolute and relative convergence tolerance to check this prune
    {
        pruneNode = true;
    }
    else
    {
        bool solReceieved, receivedSolUpdtBestSol;
        #if 0
        long int iter;
        
        r = milpSolverCallbackInterface.getNumberOfIterations(iter);
        if(r != 0)
        {
            MRQ_PRINTERRORNUMBER(r);
            
            myiter++; //so we use our own count. The ideal would be use a emaphore here, but it would be so expensive, and I think it is good enough in this way...
            iter = myiter;
        }
        #endif
        
        myiter++;
        
        if( myiter % myAlg->in_printing_frequency == 1 )
            std::cout << MRQ_PREPRINT << myiter << "-  lb: " << myzl << " ub: " << zu << " current node lb: " << nodelb << "\n";
        
        r = general_beforeSolvingRelaxation(threadNumber, algCode, &milpSolverCallbackInterface, NULL, myiter, -INFINITY, nodelb, NULL, NULL, solReceieved, receivedSolUpdtBestSol, receivedObjAndSol, pruneNode);
        MRQ_IFERRORRETURN(r, r);
        
        if( setLazyConstraintsOnReceivedSol && receivedSolUpdtBestSol ) //so, we have received a solution and this solution updates best solution found. Note: just one thread can receive a solution
        {	
            MRQ_PRINTMSG("Received solution and updated best solution. Adding lazzy constraints\n");
            
            r = solsToBuildLazyConstraints.addSolution(n, receivedObjAndSol[0], &receivedObjAndSol[1]);
            MRQ_IFERRORRETURN(r, r);
        }
        
    }
    
    
    if( pruneNode ) 
    {
        //MRQ_PRINTMSG("pruning node\n");
        r = milpSolverCallbackInterface.pruneCurrentNode();
        MRQ_IFERRORRETURN(r, r);
    }
    
    return 0;
}



int MRQ_ServerBBCallbacks::linearApp_branchingInMILPBB(const MRQ_ALG_CODE algCode, const unsigned int threadNumber, MRQ_MILPSolverCallbackInterface &milpSolverCallbackInterface)
{
    bool pruneNode = false;
    int r;
    double nodelb;
    const double zu = getUpperBound();
    
    //std::cout << "threadNumber: " << threadNumber << "\n";
    
    r = milpSolverCallbackInterface.getCurrentNodeLowerBound(nodelb);
    MRQ_IFERRORRETURN(r, r);
    
    if( nodelb >= zu ) //TODO: consider here absolute and relative convergence tolerance to check this prune
    {
        pruneNode = true;
    }
    else
    {
        bool solReceieved, receivedSolUpdtBestSol;
        
        myiter++;
        
        if( myiter % myAlg->in_printing_frequency == 1 )
            std::cout << MRQ_PREPRINT << myiter << "-  lb: " << myzl << " ub: " << zu << " current node lb: " << nodelb << "\n";
        
        r = general_beforeSolvingRelaxation(threadNumber, algCode, &milpSolverCallbackInterface, NULL, myiter, -INFINITY, nodelb, NULL, NULL, solReceieved, receivedSolUpdtBestSol, receivedObjAndSol, pruneNode);
        MRQ_IFERRORRETURN(r, r);
    }
    
    
    if( pruneNode ) 
    {
        //MRQ_PRINTMSG("pruning node\n");
        r = milpSolverCallbackInterface.pruneCurrentNode();
        MRQ_IFERRORRETURN(r, r);
    }
    
    return 0;
    
    
}



int MRQ_ServerBBCallbacks::endOfIteration( const MRQ_ALG_CODE algCode, const unsigned int threadNumber, const long unsigned int iter, const double cpuTime, const double wallTime, const double lb, const double ub)
{
    BB_printOpenNodesList();
    MRQ_getchar();
    return 0;
}



int MRQ_ServerBBCallbacks:: BB_afterSolvingRelaxation(const unsigned int threadNumber, MRQ_NewBBNode &node, const long unsigned int iter, const double lb, const double ub, const double *nlx, const double *nux, MRQ_NLPSolver &nlpSolver, const int status, double &objFSolution, double &dualObjFSolution, double *sol, double *constrs, double *dualSolC, double *dualSolV, bool &pruneNode, bool &solveAgain, bool &branchEvenIntegerSol)
{
    #if MRQ_DEBUG_MODE
        if( nlpSolver.retCode == optsolvers::OPT_OPTIMAL_SOLUTION && nlpSolver.objValue < lb - 0.05 )
        {
            std::cout << MRQ_PREPRINT "warning: node relaxation lower than lower bound. nlp.objValue: " << nlpSolver.objValue << " lb: " << lb << " nlpSolver.retCode: " << nlpSolver.retCode << " nlpSolver.origSolverRetCode: " << nlpSolver.origSolverRetCode <<  " objFSolution: " << objFSolution << "\n";
            //MRQ_getchar();
        }
    #endif
    
    pruneNode = false;
    solveAgain = false;
    return 0;
}


void MRQ_ServerBBCallbacks:: newBestSolution( const unsigned int threadNumber, const double *newSol, const double oldBestObj, const double newBestObj, const long unsigned int iter )
{
    /*
    * Since we write in the socket to send the new solution, we take advantage to send lower bound also.
    *
    * Here, we write:
    * 		1 - The conncetion code DCT_NEW_SOLUTION_FOUND (int32)
    * 		2 - The number n of variables in the problem (ATTENTION: to turn easier to client, we write this value as DOUBLE)
    * 		3 - The lower bound of the problem (double)
    * 		4 - The objective value in the new solution (double)
    * 		5 - n double values specifing the solution (double)
    * 
    */
    
    int r;
    const double lb = getLowerBound(); //linear approximation based BB could return -MRQ_INFINITY for this value. But, there is no problem about that (I think)...
    
        
    
    DCT_Int32 ccode = DCT_CC_NEW_SOLUTION_FOUND;
    
    //MRQ_getchar();
    
    std::cout << MRQ_PREPRINT "New solution found. Obj: " << newBestObj << " sending to client." << std::endl;
    
    SEMAPH_connectionWrite.lock(nthreads);
    {
        r = connection->writeData(&ccode, sizeof(ccode));
        
        if(r != 0)
        {
            MRQ_PRINTERRORNUMBER(r);
            stopAlg = true; //we cannot stop algorithm from here
        }
        
        auxSentSol[0] = n;
        auxSentSol[1] = lb;
        auxSentSol[2] = newBestObj;
        MRQ_copyArray(n, newSol, &auxSentSol[3]);
        
        r = connection->writeData(auxSentSol, (n+3)*sizeof(auxSentSol[0]));
    }
    SEMAPH_connectionWrite.unlock(nthreads);
    
    if(r != 0)
    {
        MRQ_PRINTERRORNUMBER(r);
        stopAlg = true; //we cannot stop algorithm from here
    }
}



void MRQ_ServerBBCallbacks::afterAll(const MRQ_ALG_CODE algCode, const long unsigned int iters, const double cpuTime, const double wallTime, const double lb, const double ub)
{
    int r;
    
    if( openNodesStorer.nNodes > 0 )
    {
        #if MRQ_DEBUG_MODE
            assert( accNodesToSend );
            assert( algCode != MRQ_BB_ALG );
        #endif
        
        MRQ_PRINTMSG("sending remainder nodes to client");
        
        SEMAPH_connectionWrite.lock(nthreads);
        {
            r = nodeWriter.write( *connection, DCT_CC_OPEN_NODE_RESPONSE, openNodesStorer );
        }
        SEMAPH_connectionWrite.unlock(nthreads);
        
        accNodesToSend = false;
        openNodesStorer.nNodes = 0;
        
        //I am not sure what I should do if we get some error in the socket
    }
    
}




int MRQ_ServerServiceCoreGenerator:: generateServiceCore(DCT_ServiceCore* &serviceCore)
{
    serviceCore = new (std::nothrow) MRQ_ServerServiceCore;
    MRQ_IFMEMERRORRETURN(!serviceCore);
    
    return 0;
}



MRQ_ServerServiceCore::MRQ_ServerServiceCore():dctools::DCT_ServiceCore()
{
    initialize();
}


MRQ_ServerServiceCore::~MRQ_ServerServiceCore()
{
    deallocate();
}


void MRQ_ServerServiceCore::deallocate()
{
    MRQ_secDelete(prob);
    MRQ_secFree(olx);
    MRQ_secFree(oux);
}



void MRQ_ServerServiceCore::initialize()
{
    prob = NULL;
    olx = NULL;
    oux = NULL;
}



int MRQ_ServerServiceCore::run(DCT_BBServer *server, DCT_ServerConnection *connection, const long unsigned int serviceNumber, const unsigned int nthreads, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_FileNames &inputFiles, DCT_AllGeneralParams &allGeneralParams, DCT_Int64 sizeOfOpenNodeRep, DCT_Byte *openNodeRep, double lowerBound, double upperBound, DCT_UInt32 nRequestedNodes, bool &stopService, DCT_SERVER_CLOSING_CODE &responseCodeToClientIfStop )
{
    int retCode = DCT_RC_UNDEFINED_ERROR;
    
    try
    {
        retCode = __run(server, connection, serviceNumber, nthreads, nBasicInputParameters, basicInputParameters, inputFiles, allGeneralParams, sizeOfOpenNodeRep, openNodeRep, lowerBound, upperBound, nRequestedNodes, stopService, responseCodeToClientIfStop);
    }
    catch( MRQ_SegFaultException &e )
    {
        DCT_FinalResults fr; //we create a final results object to call sendResume and finalize this execution. Note the optRetCode will be undefined error by default in DCT_FinalResults
        
        std::cerr << MRQ_PREPRINT "Service " << serviceNumber << " generated a segmentation fault!\n";
        
        
        sendResume(connection, fr, 0, NULL);
    }
    catch( MRQ_AbortException &e )
    {
        DCT_FinalResults fr; //we create a final results object to call sendResume and finalize this execution. Note the optRetCode will be undefined error by default in DCT_FinalResults
        
        std::cerr << MRQ_PREPRINT "Service " << serviceNumber << " generated an abort signal (maybe some assert?)!\n";
        
        sendResume(connection, fr, 0, NULL);
    }
    catch(...)
    {
        DCT_FinalResults fr; //we create a final results object to call sendResume and finalize this execution. Note the optRetCode will be undefined error by default in DCT_FinalResults
        
        //we cath any exception to avoid our server be aborted. In this way, only the serviceNumber will be interrupted
        std::cerr << MRQ_PREPRINT "Service " << serviceNumber << " generated some exception!\n";
        
        sendResume(connection, fr, 0, NULL);
    }
    
    return retCode;
}



int MRQ_ServerServiceCore::__run( DCT_BBServer *server, DCT_ServerConnection *connection, const long unsigned int serviceNumber, const unsigned int nthreads, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_FileNames &inputFiles, DCT_AllGeneralParams &allGeneralParams, DCT_Int64 sizeOfOpenNodeRep, DCT_Byte *openNodeRep, double lowerBound, double upperBound, DCT_UInt32 nRequestedNodes, bool &stopService, DCT_SERVER_CLOSING_CODE &responseCodeToClientIfStop )
{
    int retCode = 0, r;
    int algCode = MRQ_UNDEFINED_ALG;
    
    DCT_UInt32 numberOfReceivedNodes = 0;
    
    std::string *inputFileName;
    std::string newInputFileName;
    
    minlpproblem::MIP_NonLinearEval *myEval = NULL;
    
    MRQ_Algorithm *alg = NULL;
    
    MRQ_ServerBBCallbacks bbcallbacks(connection);
    
    MRQ_GeneralSolverParams  muriquiParams, milpParams, nlpParams, globalParams;
    
    DCT_FinalResults finalResults;
    DCT_VarBounds *rootVarsBounds = NULL;
    
    const DCT_UInt32 nIntBasicInputParameters = nBasicInputParameters/sizeof(DCT_Int32);
    
    
    if( nIntBasicInputParameters > 0 )
    {
        algCode = ((DCT_Int32*) basicInputParameters)[0];
    }
    
    
    if( !MRQ_isServerAlgorithm(algCode) )
        algCode = MRQ_LP_BB_ECP_BASED_ALG; //our default algorithm. We do not allow server run remainder algorithms
    
    
    alg = MRQ_newAlgorithm(algCode, -1);
    MRQ_IFMEMERRORGOTOLABEL(!alg, retCode, termination);
    
    
    
    //std::cout << "prob: " << prob << "\n";
    
    if( prob == NULL )
    {
        if( inputFiles.size() == 0 )
        {
            MRQ_PRINTERRORMSG("Input file is missing");
            retCode = MRQ_BAD_DEFINITIONS;
            goto termination;
        }
        
        prob = new (std::nothrow) MRQ_MINLPProb;
        MRQ_IFMEMERRORGOTOLABEL(!prob, retCode, termination);
        
        
        if( inputFiles.count(MRQ_IFT_AMPL_MODEL_FILE) > 0 )
        {
            MRQ_ReadAmplModel reader;
            
            inputFileName = &inputFiles[MRQ_IFT_AMPL_MODEL_FILE];
            
            
            std::cout << "inputFileName: " << *inputFileName << "\n";
            
            //checking if file has extension .nl
            if( !( inputFileName->size() < 3 || inputFileName->substr( inputFileName->size() -3, std::string::npos) == ".nl") )
            {
                newInputFileName = *inputFileName + ".nl";
                
                //ampl only reads the file if it has the extension .nl . I HATE AMPL!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                r = rename(inputFileName->c_str(), newInputFileName.c_str() );
                if( r == 0 )
                {
                    inputFiles[MRQ_IFT_AMPL_MODEL_FILE] = newInputFileName;
                }
            }
            
            r = reader.readProblem( (char*) inputFiles[MRQ_IFT_AMPL_MODEL_FILE].c_str() , *prob, myEval);
            MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_UNDEFINED_ERROR, termination);
        }
        else if( inputFiles.count(MRQ_IFT_GAMS_MODEL_FILE) > 0 )
        { //put here code to read a model from gams. It was not tested yet!!!!!!!
            minlpproblem::MIP_GamsModelReader reader;
            
            inputFileName = &inputFiles[MRQ_IFT_GAMS_MODEL_FILE];
            
            
            DCT_PRINTMSG("Reading GAMS model!");
            
            r = reader.readProblem( inputFileName->c_str(), *prob, myEval );
            MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_UNDEFINED_ERROR, termination);
        }
        else
        {
            MRQ_PRINTERRORMSG("Input file type not supported");
            return MRQ_BAD_DEFINITIONS;
        }
        
        MRQ_malloc(olx, prob->n);
        MRQ_malloc(oux, prob->n);
        MRQ_IFMEMERRORGOTOLABEL(!olx || !oux, retCode, termination);
        
        prob->getVariableLowerBounds(olx);
        prob->getVariableUpperBounds(oux);
    }
    else
    { //so, we already have the MINLP prob from a previous execution. We only restore the original bounds before owerwrite with the node bounds
        prob->setVariableLowerBounds(prob->n, olx);
        prob->setVariableUpperBounds(prob->n, oux);
        
        if( inputFiles.count(MRQ_IFT_AMPL_MODEL_FILE) > 0 )
        {
            inputFileName = &inputFiles[MRQ_IFT_AMPL_MODEL_FILE];
        }
        else if( inputFiles.count(MRQ_IFT_GAMS_MODEL_FILE) > 0 )
        {
            inputFileName = &inputFiles[MRQ_IFT_GAMS_MODEL_FILE];
        }
        else
        {
            MRQ_PRINTERRORMSG("Input file type not supported");
            return MRQ_BAD_DEFINITIONS;
        }
    }
    
    
    
    
    //reading the node for this problem
    if( sizeOfOpenNodeRep > 0 )
    {
        DCT_Int32 nodeRepVersion;
        DCT_Byte *paux = openNodeRep;
        
        DCT_readAndShift(paux, nodeRepVersion);
        
        #if DCT_DEBUG_MODE
            assert(nodeRepVersion == MRQ_CURRENT_NODE_REPRESENTATION_VERSION); //when we have more than one version of node representation, remove that
        #endif
        
        bbcallbacks.rootOpenNodeRep = paux;
        DCT_readAndShift(paux, numberOfReceivedNodes);
        
        //std::cout << "numberOfNodes: " << numberOfReceivedNodes << "\n";
        //MRQ_getchar();
        
        if(numberOfReceivedNodes == 1)
        {
            double nodelb;
            
            bbcallbacks.rootOpenNodeRep = NULL; //we have only one node to start exploration. So, we will replace problem bounds directly and will not generate root node.
            
            rootVarsBounds = new (std::nothrow) DCT_VarBounds;
            MRQ_IFMEMERRORGOTOLABEL(!rootVarsBounds, retCode, termination);
            
            
            int r = MRQ_readNodeFromOpenNodeBufferAndShift(paux, nodelb, *rootVarsBounds);
            MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            #if MRQ_DEBUG_MODE
                assert( &paux[-1] == &openNodeRep[sizeOfOpenNodeRep-1] );
            #endif
        }
        
        
    }
    
    
    
    
    //setting new variable bounds
    std::cout <<MRQ_PREPRINT "lower bound: " << lowerBound << " upperBound: " << upperBound << "\n";
    
    if(rootVarsBounds)
    {
        std::cout << MRQ_PREPRINT "My bounds to this problem: \n";
        
        for(auto &pairvbounds : *rootVarsBounds)
        {
            const auto &ind = pairvbounds.first;
            const DCT_Bounds &bounds = pairvbounds.second;
            
            int r = prob->setVariableLowerBound(ind, bounds.lb);
            MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_BAD_DEFINITIONS, termination);
            
            r = prob->setVariableUpperBound(ind, bounds.ub);
            MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_BAD_DEFINITIONS, termination);
            
            std::cout << "ind: " << ind << " l: " << bounds.lb << " u: " <<  bounds.ub << "\t";
        }
    }
    std::cout << "\n";
    //MRQ_getchar();
    
    
    //if(prob->n <= 30)
        //prob->print();
    
    
    
    if( nIntBasicInputParameters > 0 )
    {
        DCT_Int32 *intBasicInputParameters = (DCT_Int32*) basicInputParameters;
        //std::cout << "basicInputParameters[0]: " << basicInputParameters[0] << "\n";
        
        if( nIntBasicInputParameters >= 2 && intBasicInputParameters[1] != MRQ_UNDEFINED_MILP )
            alg->in_milp_solver = MRQ_intToMILPSolver( intBasicInputParameters[1]);
        
        if( nIntBasicInputParameters >= 3 && intBasicInputParameters[2] != MRQ_UNDEFINED_NLP )
            alg->in_nlp_solver = MRQ_intToNLPSolver( intBasicInputParameters[2]);
        
        if( nIntBasicInputParameters >= 5)
            bbcallbacks.maxNumberOfNodesToSend = intBasicInputParameters[4];
        
        if( nIntBasicInputParameters >= 6)
            bbcallbacks.frequencyToSendLowerBound = intBasicInputParameters[5];
        
        //std::cout << "intBasicInputParameters[4]: " << intBasicInputParameters[4] << "\n";
        //MRQ_getchar();
    }
    
    
    #if 1
    //setting muriqui parameters
    {
        r = MRQ_buildParamsFromDCTParams(muriquiParams, *allGeneralParams[MRQ_DC_GPT_MURIQUI]);
        MRQ_IFERRORRETURN(r, r);
        
        r = MRQ_buildParamsFromDCTParams(milpParams, *allGeneralParams[MRQ_DC_GPT_MILP_SOLVER]);
        MRQ_IFERRORRETURN(r, r);
        
        r = MRQ_buildParamsFromDCTParams(nlpParams, *allGeneralParams[MRQ_DC_GPT_NLP_SOLVER]);
        MRQ_IFERRORRETURN(r, r);
        
        r = MRQ_buildParamsFromDCTParams(globalParams, *allGeneralParams[MRQ_DC_GPT_GLOBAL_SOLVER]);
        MRQ_IFERRORRETURN(r, r);
        
        /*std::cout << "\nmuriqui parameters: \n";
        muriquiParams.print();
        
        std::cout << "\nmilp parameters: \n";
        milpParams.print();
        
        std::cout << "\nnlp parameters: \n";
        nlpParams.print();
        
        std::cout << "\nglobal parameters: \n";
        globalParams.print(); */
        
    }
    #endif
    
    
    r = alg->setParameters(muriquiParams);
    if(r != 0)
    {
        MRQ_PRINTERRORMSG("An error was gotten to set Muriqui parameters");
    }
    
    
    //std::cout <<MRQ_PREPRINT "lowerBound: " << lowerBound << " upperBound: " << upperBound << "\n";
    //DCT_getchar();
    
    
    if(upperBound < MRQ_INFINITY )
    {
        if( lowerBound > -MRQ_INFINITY && alg->out_algorithm == MRQ_BB_ALG )
        {
            if( (upperBound - lowerBound)/(MRQ_abs(upperBound) + 0.0001) < bbcallbacks.minRelativeGaptoUsePseudoCosts )
                ( (MRQ_BranchAndBound*) alg)->in_branching_strategy = MRQ_BB_BS_BIN_FIRST_HIGHEST_INT_GAP; //if the gap is lowest than x%, we disable pseudo-cost because probably we spent few iterations
        }
        
        
        alg->in_upper_bound = upperBound + MRQ_abs(upperBound)*1e-8; //we just put that because on the first ietrations, we will ask the current best solution to client, and we desire the solution sent be used to update the best solution. Due to it, we add the tolerance on upperBound. So, the solution could be useful, for example, as linearization point to outter approximation.
    }
    
    if( lowerBound > -MRQ_INFINITY )
        alg->in_lower_bound = lowerBound;
    
    alg->in_number_of_threads = nthreads;
    
    bbcallbacks.n = prob->getNumberOfVars();
    bbcallbacks.nRequestedNodes = nRequestedNodes;
    bbcallbacks.rootVarBounds = rootVarsBounds;
    bbcallbacks.olx = olx;
    bbcallbacks.oux = oux;
    bbcallbacks.myAlg = alg;
    bbcallbacks.requestFirstBestSol = true;
    
    
    r = bbcallbacks.allocateStructures(alg, prob);
    MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_MEMORY_ERROR, termination);
    
    
    alg->in_user_callbacks = &bbcallbacks;
    
    alg->in_call_new_best_solution_callback = true;
    
    if( alg->out_algorithm == MRQ_BB_ALG )
    {
        MRQ_BranchAndBound *bb = (MRQ_BranchAndBound*) alg;
        
        bb->in_call_after_bb_loop_callback = true;
        bb->in_call_before_bb_loop_callback = true;
        bb->in_call_before_solving_relax_callback = true;
        #if MRQ_DEBUG_MODE
            bb->in_call_after_solving_relax_callback = true;
        #endif
        
        if( numberOfReceivedNodes > 1)
        {
            bb->in_call_generate_root_node_callback = true; //we have to put the received nodes
            bb->in_branching_strategy = MRQ_BB_BS_BIN_FIRST_HIGHEST_INT_GAP; //we cannot use pseudo-costs here because we have several disconected root nodes
            
            //bb->in_call_end_of_iteration_callback = true;  //just to debug
        }
    }
    else if( alg->isLinearApproximationAlgorithm() )
    {
        MRQ_LinearApproxAlgorithm *la = (MRQ_LinearApproxAlgorithm *) alg;
        
        //la->in_call_before_solve_callback_in_milp_bb = true;
        la->in_call_branching_callback_in_milp_bb = true;
        
        
        if( numberOfReceivedNodes > 1 )
        {
            DCT_Int32 nodeRepVersion;
            DCT_Byte *paux = openNodeRep;
            DCT_UInt32 nnodes;
            DCT_VarBounds varBounds;
            
            double nodelb;
            
            
            DCT_readAndShift(paux, nodeRepVersion);
            #if DCT_DEBUG_MODE
                assert(nodeRepVersion == MRQ_CURRENT_NODE_REPRESENTATION_VERSION); //when we have more than one version of node representation, remove that
            #endif
            
            DCT_readAndShift(paux, nnodes);
            #if DCT_DEBUG_MODE
                assert(nnodes == numberOfReceivedNodes);
            #endif
            
            std::cout << MRQ_PREPRINT <<  "My receved nodes: ";
            
            
            finalResults.optCode = DCT_ORC_OPTIMAL_SOLUTION ;
            finalResults.lowerBound = INFINITY;
            
            //since we cannot intruce nodes received in MILP solver B&B, we will make a single run for each node
            for(decltype(nnodes) i = 0; i < nnodes; i++ )
            {
                int r = MRQ_readNodeFromOpenNodeBufferAndShift(paux, nodelb, varBounds);
                MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
                
                
                /*std::cout << MRQ_PREPRINT <<  "#########################################################\n";
                std::cout << "number of bounds: " << varBounds.size() <<  "\n";
                for(auto &pairvbounds : varBounds)
                {
                    const auto &ind = pairvbounds.first;
                    const DCT_Bounds &bounds = pairvbounds.second;
                    
                    std::cout << "ind: " << ind << " l: " << bounds.lb << " u: " <<  bounds.ub << "\t";
                }
                std::cout << "\n"; */
                
                
                
                for(auto &pairvbounds : varBounds)
                {
                    const auto &ind = pairvbounds.first;
                    const DCT_Bounds &bounds = pairvbounds.second;
                    
                    int r = prob->setVariableLowerBound(ind, bounds.lb);
                    MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_BAD_DEFINITIONS, termination);
                    
                    r = prob->setVariableUpperBound(ind, bounds.ub);
                    MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_BAD_DEFINITIONS, termination);
                }
                
                
                
                bbcallbacks.rootVarBounds = &varBounds;
                
                alg->run(*prob, &milpParams, &nlpParams);
                
                printf("____________________________________________________________________________\n");
                printf("Muriqui Server - algorithm: %s (%d) Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld\n", alg->getAlgorithmName().c_str(), alg->out_algorithm,  inputFileName->c_str(), alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
                printf("____________________________________________________________________________\n");
                
                alg->in_print_parameters_values = false; //we do not want print parameters several times
                bbcallbacks.requestFirstBestSol = false; //we do not need request first best sol in the next node run.
                
                //if( alg->out_return_code != MRQ_INFEASIBLE_PROBLEM )
                    //MRQ_getchar();
                
                if(bbcallbacks.stopService)
                {
                    stopService = true;
                    responseCodeToClientIfStop = DCT_SCC_REQUEST_BY_CLIENT;
                }
                
                if( alg->out_upper_bound < alg->in_upper_bound )
                    alg->in_upper_bound = alg->out_upper_bound; //passing a possible better upper bound to next iterations
                
                
                if( alg->out_return_code != MRQ_OPTIMAL_SOLUTION && alg->out_return_code != MRQ_INFEASIBLE_PROBLEM )
                    MRQ_retCode2DCT_optRetCode(alg->out_return_code);
                
                if( alg->out_number_of_threads > finalResults.nthreads )
                    finalResults.nthreads = alg->out_number_of_threads;
                
                finalResults.nIters += alg->out_number_of_iterations;
                
                if( alg->out_lower_bound < finalResults.lowerBound )
                    finalResults.lowerBound = alg->out_lower_bound;
                
                if( alg->out_upper_bound < finalResults.upperBound )
                    finalResults.upperBound = alg->out_upper_bound;
                
                finalResults.cpuTime += alg->out_cpu_time;
                finalResults.clockTime += alg->out_clock_time;
                
                
                //restoring the original bounds
                r = prob->setVariableLowerBounds(prob->n, olx);
                MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_UNDEFINED_ERROR, termination);
                
                r = prob->setVariableUpperBounds(prob->n, oux);
                MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_UNDEFINED_ERROR, termination);
                
                varBounds.clear();
            }
            
            finalResults.feasSolFound = false; //since we have several executions, we always set feasSolFound as false. I believe the is no problem about that since feasible soluteions are sent to client as soon as they era found.
            finalResults.algorithm = alg->out_algorithm;
            
            
            
            //MRQ_getchar();
        }
            
        
    }
    
    
    if( alg->out_algorithm == MRQ_BB_ALG || numberOfReceivedNodes <= 1   )
    {
        alg->run(*prob, &milpParams, &nlpParams);
        
        printf("____________________________________________________________________________\n");
        printf("Muriqui Server - algorithm: %s (%d) Problem: %s Return code: %d obj function: %0.10f time: %0.2f cpu time: %0.2f iters: %ld\n", alg->getAlgorithmName().c_str(), alg->out_algorithm,  inputFileName->c_str(), alg->out_return_code, alg->out_best_obj, alg->out_clock_time, alg->out_cpu_time, alg->out_number_of_iterations);
        printf("____________________________________________________________________________\n");
        
        if(bbcallbacks.stopService) //( alg->out_return_code == MRQ_CALLBACK_FUNCTION_ERROR && alg->out_user_callback_error_code == DCT_SCC_REQUEST_BY_CLIENT )
        {
            stopService = true;
            responseCodeToClientIfStop = DCT_SCC_REQUEST_BY_CLIENT;
        }
        
        finalResults.feasSolFound = alg->out_feasible_solution;
        finalResults.algorithm = alg->out_algorithm;
        finalResults.optCode = MRQ_retCode2DCT_optRetCode(alg->out_return_code);
        finalResults.nthreads = alg->out_number_of_threads;
        finalResults.nIters = alg->out_number_of_iterations;
        
        finalResults.nServerCalls = 0; //some day, each server can call others servers (subservers) to perform the job
            
        finalResults.lowerBound = alg->out_return_code == MRQ_INFEASIBLE_PROBLEM ? MRQ_INFINITY : alg->out_lower_bound;
        
        finalResults.upperBound = alg->out_upper_bound;
        if(alg->out_feasible_solution)
            finalResults.objBestSol = alg->out_best_obj;
        finalResults.objAtFirstRelax = alg->out_obj_opt_at_continuous_relax;
        finalResults.cpuTime = alg->out_cpu_time;
        finalResults.clockTime= alg->out_clock_time;
    
    }
    
    
    
    r = sendResume(connection, finalResults, alg->out_feasible_solution ? prob->n: 0, alg->out_best_sol );
    MRQ_IFERRORGOTOLABEL(r, retCode, MRQ_UNDEFINED_ERROR, termination);
    
    
termination:
    
    #if MRQ_DEBUG_MODE
        assert( bbcallbacks.openNodesStorer.nNodes == 0 ); //we cannot have remainder nodes in openNodesStorer
    #endif
    
    
    if(alg)		delete alg;
    if(rootVarsBounds)	delete rootVarsBounds;
    
    
    return retCode;
}

