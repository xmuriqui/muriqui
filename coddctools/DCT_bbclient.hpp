

#ifndef DCT_BBSERVER_HPP
#define DCT_BBSERVER_HPP



#include "DCT_dctools.hpp"
#include "DCT_sockets.hpp"

#include <vector>
#include <queue>
#include <ostream>
#include <condition_variable>


namespace dctools
{
	
	#define DCT_CLIENT_DEBUG_MODE 1
	
	#define DCT_SECONDS_TO_SLEEP_BEFORE_CHECK_AGAIN_ALL_SERVER_END 1
	
	
	//version of client is defined by MAJOR_VERSION.MINOR_VERSION ex: 0.3, 1.6, 2.1 etc 
	
	
	class DCT_ServerError
	{
	public:
		unsigned int server;
		int error;
	};
	
	/*
	 * Class to manage the queue of DCT_VarBounds. This class have already a internal mutex to perform a producer/consumer multithrading task. Actually, this class works like a receptor of DCT_VarBound and this purpose is receive node bounds wich are requested by servers.
	 */ 
	class DCT_VarBoundsReceptorQueue
	{
		bool memoryError;
		bool able;
		
	public:
		
		std::mutex mtx;
		std::condition_variable cv;
		std::queue<DCT_Byte*> queue;
		DCT_Mutex SEMAPH_queue;
		DCT_Mutex SEMAPH_readCheck;
		
		unsigned int secondsToSleep;
		
		DCT_VarBoundsReceptorQueue();
		
		~DCT_VarBoundsReceptorQueue();
		
		
		void disable();
		
		void enable();
		
		bool isEmpty();
		
		//push a DCT_VarBounds pointer in the queue and take its ownership
		int push(DCT_Byte* vb);
		
		int pop(DCT_Byte* &vb);
		
	};
	
	
	class DCT_BBClient
	{
		bool server0FinishedSomeExecution;
		bool someOptimalRetCode;
		bool someServerEnded;
		DCT_Byte *serverRunning;
		DCT_Byte **serverPVarsBounds; //we store the node bouns in each server. So if server go down or disconect by some reazon, we can still atribute the node for other server
		int *threadReturnCodes;
		double *serverLowerBounds;
		
		
		DCT_Mutex SEMAPH_bestSol;
		DCT_Mutex SEMAPH_serverRunning;
		DCT_Mutex SEMPAH_FinalResult;
		DCT_Mutex *SEMAPH_socketsWrite;
		DCT_Mutex SEMAPH_checkOpenPVarBounds;
		DCT_Mutex SEMAPH_serverError;
		DCT_Mutex SEMAPH_initialOpenNodes;
		
		std::vector<DCT_Byte*> initialOpenNodes; //vector to the first server put the initial open nodes to other servers
		
		DCT_VarBoundsReceptorQueue *varBoundsRecepQueue;
		std::vector<DCT_ClientSocket> sockets;
		
		
		
		
		int clientThread(const unsigned int serverNumber, DCT_ClientSocket &socket, DCT_FileNames *inputFiles, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_AllGeneralParams *allGeneralParams, DCT_UInt32 nInitialRequestedNodes, DCT_Byte *pvarsBounds);
		
		int getService( dctools::DCT_ClientSocket& socket, const unsigned int maxThreads, unsigned int& nthreads);
		
		void finishServer(const unsigned int serverNumber, const int clientThreadReturnCode);
		
		int sendBasicDataToServer(DCT_ClientSocket &socket, DCT_FileNames *inputFiles, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_AllGeneralParams *allGeneralParams, DCT_UInt32 numberOfOpenNodesToSendToOtherServers);
		
		//if objValue is NULL, we assume objective value is in the first index in solution
		int sendSolution(DCT_ClientSocket &socket, DCT_UInt32 nvars, const double *solution, const double *objValue = NULL, std::ostream &out = std::cerr);
		
		/* send a message with the best solution found for all servers, except the server indicated by excluded variable. If objValue is NULL, we assume objValue is in the first position of solution */
		int sendSolutionToServers(DCT_UInt32 nvars, const double *objValue, const double *solution, unsigned int excluded, std::ostream &out);
		
		int startService(DCT_ClientSocket &socket, DCT_Byte *openNodeRep, double lb);
		
		int treatServerRequests( DCT_Int32 connectCode, unsigned int serverNumber, double* &auxSol, DCT_FinalResults *finalResults, double* &finalSol, DCT_UInt32 &nvarsFinalSol, std::ostream &out);
		
		int tryUpdateBestSol(const unsigned int serverNumber, const unsigned int n, double objValue, const double *solution, bool *updated = NULL, std::ostream &out = std::cerr);
		
		
		
	public:
		
		//TODO: erase this two attributes below
		//DCT_UInt32 solutionSize; //size of solution (in bytes)
		//void *solution;
		
		unsigned int nConnectedServers;
		unsigned int nvars;
		unsigned int secondsToSleepWaitingANode; //seconds to a thread sleep waiting to receive a node from a server
		unsigned int maxSimultaneuousServersToSendBasicData; //number of maximum simultaneous servers that we send basic input data
		double bestObjValue;
		double *bestSol;
		double *arrayBestObjAndSol; //this array accumulates space for objective value and solution
		
		DCT_FinalResults finalResults;
		
		std::vector<DCT_ServerError> serverError;
		
		
		DCT_BBClient();
		
		~DCT_BBClient();
		
		void deallocate();
		
		void initialize();
		
		
		
		int run(DCT_Components *servers, 	DCT_SeveralNumberOfComponents *allNumberOfServers,  DCT_FileNames *inputFiles, DCT_UInt64 nBasicInputParameters = 0, DCT_Byte *basicInputParameters = NULL,  DCT_AllGeneralParams *allGeneralParams = NULL);
		
		
		//this function enable the usage of a handler to threat the SIGINT and SIGTERM signals (generated when user press CTRL+C and kill comand respectivelly). In this case, application will treat the signal sending messages to serves stop their services and close 
		int setSigIntHandler();
		
		
		//function to iterrupt and stop the services in progress: useful if the program receives a signal to STOP
		void stopOperations();
		
		
		//this function disable the usage of handle to treat SIGINT signal
		int unsetSigIntHandler();
		
		
		friend int DCT_clientThread(DCT_BBClient *bbclient, const unsigned int serverNumber, DCT_ClientSocket *socket, DCT_FileNames *inputFiles, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_AllGeneralParams *allGeneralParams, DCT_UInt32 nInitialRequestedNodes, DCT_Byte *pvarsBounds);
		
		friend int DCT_semaphSendBasicDataToServer(DCT_BBClient *bbclient, unsigned int serverNumber, DCT_Mutex *SEMAPH_socketWrite, DCT_ClientSocket *socket, DCT_FileNames *inputFiles, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_AllGeneralParams *allGeneralParams, DCT_UInt32 numberOfOpenNodesToSendToOtherServers);
	};
	
	
	int DCT_clientThread(DCT_BBClient *bbclient, const unsigned int serverNumber, DCT_ClientSocket *socket, DCT_FileNames *inputFiles, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_AllGeneralParams *allGeneralParams, DCT_UInt32 nInitialRequestedNodes, DCT_Byte *pvarsBounds);
	
	
	int DCT_semaphSendBasicDataToServer(DCT_BBClient *bbclient, unsigned int serverNumber, DCT_Mutex *SEMAPH_socketWrite, DCT_ClientSocket *socket, DCT_FileNames *inputFiles, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_AllGeneralParams *allGeneralParams, DCT_UInt32 numberOfOpenNodesToSendToOtherServers);
	
}




#endif
