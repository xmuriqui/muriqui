
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cstring>
#include <csignal>

#include <chrono>
#include <new>
#include <vector>
#include <exception>
#include <set>

#include <iostream>
#include <fstream>

#include "DCT_bbclient.hpp"
#include "DCT_tools.hpp"


#include <thread>



using namespace dctools;


template <class myClass>
int DCT_writeParametersFromMap( dctools::DCT_BaseSocket &socket, const std::map<std::string, myClass> &map )
{
	DCT_UInt32 k = 0, size = map.size();
	int r;
	//size_t bytesWritten;
	
	
	r = socket.writeData(&size, sizeof(size));
	if(r != 0)
	{
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORNUMBER(r);
		#endif
		return DCT_RC_SERVER_ERROR;
	}
	
	
	for( const auto &pairParam : map )
	{
		const std::string &name = pairParam.first;
		const myClass &value = pairParam.second;
		
		//writes, size of string and string
		r = DCT_writeStringWithSizeAnd0(socket, name.c_str() );
		
		if(r != 0)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORNUMBER(r);
			#endif
			return DCT_RC_SERVER_ERROR;
		}
		
		r = socket.writeData(&value, sizeof(value) );
		
		if(r != 0)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORNUMBER(r);
			#endif
			return DCT_RC_SERVER_ERROR;
		}
		
		k++;
	}
	
	#if DCT_DEBUG_MODE
		assert(size == k);
	#endif
	
	return 0;
}


template <>
int DCT_writeParametersFromMap<std::string>( dctools::DCT_BaseSocket &socket, const std::map<std::string, std::string> &map )
{
	DCT_UInt32 k = 0, size = map.size();
	int r;
	//size_t bytesWritten;
	
	
	r = socket.writeData(&size, sizeof(size));
	if(r != 0)
	{
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORNUMBER(r);
		#endif
		return DCT_RC_SERVER_ERROR;
	}
	
	
	
	for( const auto &pairParam : map )
	{
		const std::string &name = pairParam.first;
		const std::string &value = pairParam.second;
		
		
		//writes, size of string and string
		r = DCT_writeStringWithSizeAnd0(socket, name.c_str() );
		
		if(r != 0)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORNUMBER(r);
			#endif
			return DCT_RC_SERVER_ERROR;
		}
		
		
		r = DCT_writeStringWithSizeAnd0(socket, value.c_str() );
		
		if(r != 0)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORNUMBER(r);
			#endif
			return DCT_RC_SERVER_ERROR;
		}
		
		k++;
	}
	
	#if DCT_DEBUG_MODE
		assert(size == k);
	#endif
	
	return 0;
}




DCT_VarBoundsReceptorQueue::DCT_VarBoundsReceptorQueue()
{
	secondsToSleep = DCT_SECONDS_TO_SLEEP_WAITING_FOR_BOUNDS;
	memoryError = false;
	enable();
}

DCT_VarBoundsReceptorQueue::~DCT_VarBoundsReceptorQueue()
{
	std::unique_lock<std::mutex> lock(mtx);
	
	while( !queue.empty() )
	{
		delete queue.front();
		queue.pop();
	}
}


void DCT_VarBoundsReceptorQueue::disable()
{
	able = false;
	cv.notify_all();
}


void DCT_VarBoundsReceptorQueue::enable()
{
	able = true; //when this object is disable, we have no threads locked boy its mutex
}


bool DCT_VarBoundsReceptorQueue::isEmpty()
{
	return queue.empty();
}


//push a DCT_VarBounds pointer in the queue and take its ownership
int DCT_VarBoundsReceptorQueue::push(DCT_Byte* vb)
{
	int code = 0;
	//std::unique_lock<std::mutex> lock(mtx); //unique lock object call the lock method of mutex. On its destruction the method unlock will be called
	
	//even if vb is NULL, we must put it in the queue anyway, because cv.notify_one can be called before the other thread call cv.wait. So, it is necessary put the NULL in the queue to thread 
	
	SEMAPH_queue.lock();
	{
		try
		{
			queue.push(vb);
		}
		catch(std::bad_alloc &ba)
		{
			code = DCT_RC_MEMORY_ERROR;
			memoryError = true;
		}
	}
	SEMAPH_queue.unlock();
	
	
	cv.notify_one(); //cv.notify_all();
	
	
	return code;
}

int DCT_VarBoundsReceptorQueue::pop(DCT_Byte* &vb)
{
	int code = 0;
	vb = NULL;
	
	if(!able && queue.empty() )
	{
		return 0; 
	}
	
	
	
	SEMAPH_readCheck.lock(); //to guarantee just one thread in cheking the queue
	{
		std::unique_lock<std::mutex> lock(mtx); //the constructor unique_lock DOES call mtx.lock. I do not why, but condtition_variables needes that. //The destructor of unique_lock DOES call mtx.unlock.  //Do not put this declaration before SEMAPH_readCheck.lock()
		
		
		while( queue.empty() && able ) //we have to put while because C++ documentation of conditional variable says about "spurious wake-up calls without notify be called"
		{
			
			
			//std::cout << "Vou dar o wait" << std::endl;
			
			cv.wait_for(lock, std::chrono::seconds(secondsToSleep) ); //wait method will put this thread to sleep until it be notified. Before sleep, the mutex.unlock method will be called to unblock other threads. //we have to use a time also to wake up because we can have a very strange situation where this thread check is the queue is empty and, after that,  other thread put a varBounds in the queue and calls cv.notify_one before this thread calls cv.wait
			
			//this thread only will be notifyed when we call cv.notify_one ou cv.notify_all
			
			//std::cout << "Sai do wait" << std::endl;
			
		}
		
		SEMAPH_queue.lock();
		{
			if( !queue.empty() )
			{
				vb = queue.front();
				queue.pop();
			}
		}
		SEMAPH_queue.unlock();
	}
	SEMAPH_readCheck.unlock();
	
	
	if(memoryError)
		code = DCT_RC_MEMORY_ERROR;
	
	return code;
}






DCT_BBClient::DCT_BBClient()
{
	initialize();
}


DCT_BBClient::~DCT_BBClient()
{
	unsetSigIntHandler();
	deallocate();
}


void DCT_BBClient::deallocate()
{
	//solutionSize = 0;
	//DCT_secFree(solution);
	
	if(serverPVarsBounds)
	{
		for(unsigned int i = 0; i < nConnectedServers; i++)
		{
			if(serverPVarsBounds[i])
				free(serverPVarsBounds[i]);
		}
		free(serverPVarsBounds);
		serverPVarsBounds = NULL;
	}
	
	DCT_secFree(serverRunning);
	DCT_secFree(threadReturnCodes);
	DCT_secFree(serverLowerBounds);
	nvars = 0;
	DCT_secFree(arrayBestObjAndSol);
	bestSol = NULL;
	DCT_secDeleteArray(SEMAPH_socketsWrite);
	//DCT_secDeleteArray(SEMAPH_varBoundsQueue);
	
	sockets.clear();
	serverError.clear();
	
	initialOpenNodes.clear();
	
	DCT_secDeleteArray(varBoundsRecepQueue);
	
	/*if(serverPVarsBounds)
	{
		for(unsigned int i = 0; i < nConnectedServers; i++)
		{
			if(serverPVarsBounds)
				delete serverPVarsBounds;
		}
		
		free(serverPVarsBounds);
		      serverPVarsBounds = NULL;
	}*/
}


void DCT_BBClient::initialize()
{
	//solutionSize= 0;
	
	server0FinishedSomeExecution = false;
	someServerEnded = false;
	someOptimalRetCode = false;
	
	nConnectedServers = 0;
	secondsToSleepWaitingANode = DCT_SECONDS_TO_SLEEP_WAITING_FOR_BOUNDS;
	
	maxSimultaneuousServersToSendBasicData = 8;
	
	//solution = NULL;
	serverRunning = NULL;
	serverPVarsBounds = NULL;
	threadReturnCodes = NULL;
	serverLowerBounds = NULL;
	
	nvars = 0;
	bestObjValue = INFINITY;
	bestSol = NULL;
	arrayBestObjAndSol = NULL;
	
	SEMAPH_socketsWrite = NULL;
	//SEMAPH_varBoundsQueue = NULL;
	
	varBoundsRecepQueue = NULL;
	serverPVarsBounds = NULL;
	
	finalResults.initialize();
}


int DCT_BBClient::getService( DCT_ClientSocket &socket, const unsigned int maxThreads, unsigned int &nthreads )
{
	const unsigned int ndatareq = 2;
	const unsigned int ndataresp = 2;
	int r; 
	size_t bytesTransfered;
	DCT_UInt32 data[2] = { DCT_CC_SERVICE_REQUEST, maxThreads };
	
	nthreads = 0;
	
	
	
	r = socket.writeData(data, ndatareq * sizeof(data[0]), &bytesTransfered );
			
	if( bytesTransfered != ndatareq * sizeof(data[0]) )
	{
		DCT_PRINTERRORMSG("Error to write data in socket\n");
		
		DCT_getchar();
		return DCT_RC_UNDEFINED_ERROR;
	}
	
	
	bytesTransfered = 0;
	do
	{
		size_t bytest;
		
		r = socket.readData(data, ndataresp * sizeof(data[0]), &bytest );
		if( r != 0 )
		{
			DCT_PRINTERRORNUMBER(r);
			return DCT_RC_UNDEFINED_ERROR;
		}
		
		bytesTransfered += bytest;
		
	} while( bytesTransfered < ndataresp * sizeof(data[0]) );
	
	{
		std::string str;
		
		DCT_enumToStr((DCT_CONECTION_CODE) data[0], str);
		
		std::cout << DCT_PREPRINT "server response: " << data[0] << " (" << str << ") reserved threads: " << data[1] << "\n";
	}
	
	if( data[0] == DCT_CC_SERVICE_ACCEPTED )
	{
		nthreads = data[1];
		return 0;
	}
	else
	{
		return DCT_RC_REFUSED_SEVICE;
	}
	
}


void DCT_BBClient::finishServer(const unsigned int serverNumber, const int clientThreadReturnCode)
{
	someServerEnded = true;
	serverRunning[serverNumber] = -1;
	
	//we set mutex because we can have some thread asking for a open node to serverNumber. So, it is better be sure there is no thread asking anything...
	SEMAPH_socketsWrite[serverNumber].lock();
	{
		sockets[serverNumber].closeSocket();
	}
	SEMAPH_socketsWrite[serverNumber].unlock();
	
	varBoundsRecepQueue[serverNumber].disable(); //it is better disable queue after close the socket. So, we are sure this queue will not more filled
	
	if( clientThreadReturnCode != 0)
	{
		DCT_ServerError se;
		
		se.server = serverNumber;
		se.error = DCT_ORC_SERVER_CONNECTION_ERROR;
		
		SEMAPH_serverError.lock();
		{
			serverError.push_back( se );
		}
		SEMAPH_serverError.unlock();
	}
	
	threadReturnCodes[serverNumber] = clientThreadReturnCode;
}


//sizeSolExpected must be in bytes
inline int DCT_readFinalResultsAndSol(DCT_BaseSocket &socket, DCT_FinalResults &finalResults, DCT_UInt32 nvarsExpected, double* &sol, DCT_UInt32 *nvarsReceived = NULL, std::ostream &out = std::cerr)
{
	int r, code;
	DCT_UInt32 sizeOfServerFinalResults, mynvarsReceived = 0;
	DCT_Byte *myBufferToReadFinalResults = NULL;
	
	
	/*
	 * Here, we read:
	 *
	 * 1 - the size of DCT_FinalResults (uint32)
	 * 
	 * 2 - the DCT_FinalResults data
	 * 
	 * 3 - Number of variables in solution (uint32)
	 *
	 * 4 - The solution
	 */ 
	
	
	r = socket.readData( &sizeOfServerFinalResults, sizeof(sizeOfServerFinalResults) );
	DCT_OSIFERRORGOTOLABEL(out, r, code, DCT_RC_SERVER_ERROR, termination);
	
	if( sizeOfServerFinalResults <= sizeof(finalResults) )
	{
		myBufferToReadFinalResults = (DCT_Byte*) &finalResults;
		
		if( sizeOfServerFinalResults < sizeof(finalResults) )
			finalResults.initialize();//we have bytes in finalResults that will not be overwritten, so we initialize it
	}
	else
	{ //size of server's Final Results is greater than client's one. It can happen if version of dctools in the server is more recent. So, we need allocate a buffer
		DCT_malloc(myBufferToReadFinalResults, sizeOfServerFinalResults);
		DCT_IFMEMERRORGOTOLABEL(!myBufferToReadFinalResults, code, termination);
	}
	
	
	r = socket.readData( myBufferToReadFinalResults, sizeOfServerFinalResults, NULL );
	DCT_OSIFERRORGOTOLABEL(out, r, code, DCT_RC_SERVER_ERROR, termination);
	
	if( myBufferToReadFinalResults != (DCT_Byte*) &finalResults )
		DCT_copyArray<DCT_Byte>( sizeof(finalResults), myBufferToReadFinalResults, (DCT_Byte*) &finalResults );
	
	
	r = socket.readData( &mynvarsReceived, sizeof(mynvarsReceived), NULL );
	DCT_OSIFERRORGOTOLABEL(out, r, code, DCT_RC_SERVER_ERROR, termination);
	
	
	if(mynvarsReceived > 0)
	{
		if( mynvarsReceived > nvarsExpected )
		{
			if( sol == NULL )
			{
				DCT_malloc(sol, mynvarsReceived);
				DCT_OSIFMEMERRORGOTOLABEL(out, !sol, code, termination);
			}
			else
			{
				out << DCT_PREPRINT "Warning: Size of solution received is greater than size of solution expected! size of solution received: " << mynvarsReceived << " size of solution expected: " << nvarsExpected << DCT_GETFILELINE << "\n";
				//code = DCT_VALUE_ERROR;
				//goto termination;
				
				assert(false);
			}
		}
		
		r = socket.readData( sol, mynvarsReceived*sizeof(sol[0]) );
		DCT_OSIFERRORGOTOLABEL(out, r, code, DCT_RC_SERVER_ERROR, termination);
		
		if(nvarsReceived)
			*nvarsReceived = mynvarsReceived; //do not put it on termination
	}
	
	
	code = 0;
	
termination:
	
	if( myBufferToReadFinalResults != NULL && myBufferToReadFinalResults != (DCT_Byte*) &finalResults )
		free(myBufferToReadFinalResults);
	
	return code;
}



int dctools::DCT_semaphSendBasicDataToServer(DCT_BBClient *bbclient, unsigned int serverNumber, DCT_Mutex *SEMAPH_socketWrite, DCT_ClientSocket *socket, DCT_FileNames *inputFiles, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_AllGeneralParams *allGeneralParams, DCT_UInt32 numberOfOpenNodesToSendToOtherServers)
{
	int rcode;
	
	SEMAPH_socketWrite->lock();
		rcode = bbclient->sendBasicDataToServer(*socket, inputFiles, nBasicInputParameters, basicInputParameters, allGeneralParams, numberOfOpenNodesToSendToOtherServers);
	SEMAPH_socketWrite->unlock();
	
	if(rcode != 0)
	{
		std::ostream *serr = &std::cerr;
		
		DCT_OSPRINTERRORNUMBER(*serr, rcode);
		bbclient->finishServer(serverNumber, rcode);
	}
	
	return rcode;
}



int  DCT_BBClient::sendBasicDataToServer( DCT_ClientSocket &socket, DCT_FileNames *inputFiles, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_AllGeneralParams *allGeneralParams, DCT_UInt32 numberOfOpenNodesToSendToOtherServers)
{
	//when we enter in this function, socket is already open and a service was accepted by the server
	
	int r, code;
	size_t bytesTransfered;
	const unsigned int maxRequest = 10;
	DCT_Int32 request[maxRequest];
	FILE *file = NULL;
	
	
	//sending general parameters
	if( allGeneralParams && allGeneralParams->size() > 0)
	{
		const unsigned int nrequest = 2;
		
		for( auto &pairGeneralParams : *allGeneralParams )
		{
			const DCT_Int32 key = pairGeneralParams.first;
			DCT_GeneralParams *params = pairGeneralParams.second;
			
			
			request[0] = DCT_CC_SOLVER_PARAMETERS;
			request[1] = key;
			
			
			r = socket.writeData(request, nrequest*sizeof(request[0]) );
			DCT_IFERRORGOTOLABEL(r, code, DCT_RC_SERVER_ERROR, termination);
			
			r = DCT_writeParametersFromMap(socket, params->dblParams);
			DCT_IFERRORGOTOLABEL(r, code, DCT_RC_SERVER_ERROR, termination);
			
			r = DCT_writeParametersFromMap(socket, params->intParams);
			DCT_IFERRORGOTOLABEL(r, code, DCT_RC_SERVER_ERROR, termination);
			
			r = DCT_writeParametersFromMap(socket, params->strParams);
			DCT_IFERRORGOTOLABEL(r, code, DCT_RC_SERVER_ERROR, termination);
			
		}
	}
	
	
	//sending the basicInputParameters
	if( nBasicInputParameters > 0 )
	{
		//const unsigned int nrequest = 2;
		DCT_Byte *paux = (DCT_Byte*) request;
		DCT_Int32 ccode = DCT_CC_BASIC_INPUT_PARAMETERS;
		
		
		DCT_writeAndShift(paux, ccode);
		DCT_writeAndShift(paux, nBasicInputParameters);
		
		//request[0] = DCT_CC_BASIC_INPUT_PARAMETERS;
		//request[1] = nBasicInputParameters;
		
		int r = socket.writeData(request, sizeof(ccode) + sizeof(nBasicInputParameters) );
		DCT_IFERRORGOTOLABEL(r, code, DCT_RC_SERVER_ERROR, termination);
		
		r = socket.writeData(basicInputParameters, nBasicInputParameters* sizeof(basicInputParameters[0]) );
		DCT_IFERRORGOTOLABEL(r, code, DCT_RC_SERVER_ERROR, termination);
	}
	
	
	#if 0
	//sending variable bounds
	if( varsBounds && varsBounds->size() > 0)
	{
		r = DCT_writeVarBounds(*socket, *varsBounds);
		
		if(r != 0)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORNUMBER(r);
			#endif
			code = r;
			goto termination;
		}
	}
	#endif
	
	//sendig the input files
	if(inputFiles)
	{
		
		//checking if input file exists
		for(auto &pairInputFile : *inputFiles)
		{
			std::string &inputFile = pairInputFile.second;
			
			const unsigned int nrequest = 2;
			
			DCT_Int32 fileType = pairInputFile.first;
			DCT_UInt64 fileSize; 
			
			
			request[0] = DCT_CC_INPUT_PROBLEM_FILE;
			request[1] = fileType;
			
			
			
			file = fopen(inputFile.c_str(), "rb");
			
			if(!file)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORMSGP( "error to open input file: ", inputFile);
				#endif
				
				code = DCT_RC_NAME_ERROR;
				goto termination;
			}
			
			
			fileSize = DCT_getFileSize(file);
			
			//std::cout << "fileSize: " << fileSize << "\n";
			
			//first we send the request code and fileType
			socket.writeData(request, nrequest * sizeof(request[0]), &bytesTransfered );
			
			if( bytesTransfered != nrequest * sizeof(request[0]) )
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORMSG("error to send data to server");
				#endif
				
				code = DCT_RC_SERVER_ERROR;
				goto termination;
			}
			
			
			//now, we send the file size
			socket.writeData(&fileSize, sizeof(fileSize), &bytesTransfered);
			
			if( bytesTransfered != sizeof(fileSize) )
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORMSG("error to send data to server");
				#endif
				
				code = DCT_RC_SERVER_ERROR;
				goto termination;
			}
			
			
			//finally, we send the file
			r = socket.writeFile(file);
			
			fclose(file);
			file = NULL;
			
			DCT_IFERRORGOTOLABEL(r, code, DCT_RC_SERVER_ERROR, termination);
		}
	}
	
	
	//here, we send the information about how many nodes the server should separate in the beggining of exploration to send to other servers
	if(numberOfOpenNodesToSendToOtherServers > 0)
	{
		const unsigned int nrequest = 2;
		
		request[0] = DCT_CC_NUMBER_OF_INITIAL_OPEN_NODES_REQUEST;
		request[1] = numberOfOpenNodesToSendToOtherServers;
		
		int r = socket.writeData(request, nrequest* sizeof(request[0]) );
		DCT_IFERRORGOTOLABEL(r, code, DCT_RC_SERVER_ERROR, termination);
	}
	
	
	code = 0;
	
termination:
	
	if(file)	fclose(file);
	
	return code;
}


int DCT_BBClient::sendSolution(DCT_ClientSocket &socket, DCT_UInt32 nvars, const double *solution, const double *objValue, std::ostream &out)
{
	const DCT_UInt32 request[2] = {DCT_CC_NEW_SOLUTION_FOUND, nvars};
	int r;
	
	/*
	 * Here, we write:
	 * 		1 - The connection code DCT_NEW_SOLUTION_FOUND (int32)
	 * 		2 - The number n of variables (uint32)
	 * 		3 - The objective value (double)
	 * 		4 - The values of n variables (double)
	 * 
	 */ 
	
	
	r = socket.writeData(&request, 2*sizeof(request[0]) );
	if( r != 0 )
	{
		DCT_OSPRINTERRORNUMBER(out, r);
		return r;
	}
	
	
	if( objValue == NULL )
	{
		//we assume solution has the objective value in the first position
		r = socket.writeData(solution, (nvars+1)* sizeof(solution[0]) );
		DCT_OSIFERRORRETURN(out, r, r);
		
	}
	else
	{
		r = socket.writeData(&objValue, sizeof(objValue) );
		DCT_OSIFERRORRETURN(out, r, r);

		r = socket.writeData(solution, nvars * sizeof(solution[0]) );
		DCT_OSIFERRORRETURN(out, r, r);
	}
	
	return 0;
}


//if objValue is NULL, we assume objValue is in the first position of solution
int DCT_BBClient::sendSolutionToServers(DCT_UInt32 nvars, const double *objValue, const double *solution, unsigned int excluded, std::ostream &out)
{
	int r;
	
	for(unsigned int i = 0; i < nConnectedServers; i++)
	{
		if(i == excluded)
			continue;
		
		if(serverRunning[i] <= 0)
			continue;
			
		DCT_ClientSocket &socket = sockets[i];
		
		SEMAPH_socketsWrite[i].lock();
		{
			r = sendSolution(socket, nvars, solution, objValue, out);
		}
		SEMAPH_socketsWrite[i].unlock();
		
		if(r != 0)
		{ //here, we do not abort if some socket get error
			DCT_OSPRINTERRORNUMBER(out, r);
		}
	}
	
	return 0;
}



int DCT_BBClient::startService(DCT_ClientSocket &socket, DCT_Byte *openNodeRep, double lb)
{
	int r;
	const unsigned int maxRequest = 5;
	DCT_Int32 request[maxRequest];
	double nodelb = -INFINITY;
	
	
	request[0] = DCT_CC_CLEAR_VAR_BOUNDS;
	r = socket.writeData(request, sizeof(request[0]) );
	if(r != 0)
	{
		DCT_PRINTERRORNUMBER(r);
		return r;
	}
	
	
	//sending variable bounds
	/*if( varsBounds && varsBounds->size() > 0)
	{
		//r = DCT_writeVarBounds(socket, *varsBounds);
		r = varBoundsWriter->write(socket, *varsBounds);
		if(r != 0)
		{
			DCT_PRINTERRORNUMBER(r);
			return r;
		}
	}*/
	
	if( openNodeRep )
	{
		DCT_Int32 ccode; //just for read the connect code
		DCT_UInt64 size;
		
		DCT_Byte *auxp = openNodeRep;
		
		DCT_readAndShift(auxp, ccode);
		DCT_readAndShift(auxp, size); //after the connect code DCT_CC_OPEN_NODE_RESPONSE, we have the size (in bytes) of the array representing the nodes being sent
		
		#if DCT_DEBUG_MODE
			assert(ccode == DCT_CC_OPEN_NODE_RESPONSE);
		#endif
		
		size += sizeof(ccode) + sizeof(size); //the total size of array openNodeRep is the size of array representing open nodes plus the size o connect code and the size of number representing the size :)
		
		
		#if 0
		DCT_UInt32 nbounds;
		DCT_readAndShift(auxp, nodelb);
		DCT_readAndShift(auxp, nbounds);
		
		
		//std::cout << "sending node to server. node lb: " << nodelb << " server lb: " << lb << "\n";
			
		if( nodelb < lb )
		{
			std::cerr << "Warning. node lb is lower than server lb. node lb: " << nodelb << " server lb: " << lb << "\n";
			//assert(false);
			//DCT_getchar();
		}
		
		size = DCT_sizeofOpenNodeRep(nbounds);
		
		#endif
		
		
		r = socket.writeData( openNodeRep, size );
		if(r != 0)
		{
			DCT_PRINTERRORNUMBER(r);
			return r;
		}
	}
	
	
	
	//sending lower and upper bound
	if( !std::isinf(lb) && lb > nodelb )
	{
		/*
		 * Here, we write:
		 * 		1 - The connection code DCT_LOWER_BOUND_RESPONSE (int32)
		 * 		2 - The lower bound value (double)
		 */ 
		request[0] = DCT_CC_LOWER_BOUND_RESPONSE;
		double *pd = (double*) &request[1];
		
		pd[0] = lb;
		
		//std::cout << "Enviando lower bound response!" << std::endl;
		r = socket.writeData(request, sizeof(DCT_Int32) + sizeof(double) );
		if(r != 0)
		{
			DCT_PRINTERRORNUMBER(r);
			return r;
		}
	}
	
	if( !std::isinf(bestObjValue) )
	{
		/*
		 * Here, we write:
		 * 		1 - The connection code DCT_UPPER_BOUND_RESPONSE (int32)
		 * 		2 - The upper bound value (double)
		 */ 
		request[0] = DCT_CC_UPPER_BOUND_RESPONSE;
		double *pd = (double*) &request[1];
		
		pd[0] = bestObjValue;
		
		//std::cout << "Enviando upper bound!" << std::endl;
		r = socket.writeData(request, sizeof(DCT_Int32) + sizeof(double) );
		if(r != 0)
		{
			DCT_PRINTERRORNUMBER(r);
			return r;
		}
	}
	
	request[0] = DCT_CC_START_SERVICE;
	
	//starting the service
	//std::cout << "Iniciando o serviço!" << std::endl;
	r = socket.writeData(request, sizeof(request[0]) );
	if(r != 0)
	{
		DCT_PRINTERRORNUMBER(r);
		return r;
	}
	
	
	
	//DCT_getchar();
	
	#if 0
	//server will send to us two int32 numbers
	
	r = socket->readData(request, 2*sizeof(request[0]), NULL);
	if(r != 0)
	{
		DCT_PRINTERRORNUMBER(r);
		return r;
	}
	#endif
	
	
	return 0;
}



int dctools::DCT_clientThread(DCT_BBClient *bbclient, const unsigned int serverNumber, DCT_ClientSocket *socket, DCT_FileNames *inputFiles, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_AllGeneralParams *allGeneralParams, DCT_UInt32 nInitialRequestedNodes, DCT_Byte *pvarsBounds)
{
	int code;
	
	try 
	{
		code = bbclient->clientThread(serverNumber, *socket, inputFiles, nBasicInputParameters, basicInputParameters, allGeneralParams, nInitialRequestedNodes, pvarsBounds);
	}
	catch(...)
	{ //some exception have been raised. Here, we catch all exceptions to avoid our client stops.
		
		std::cerr << DCT_PREPRINT "Some exception was raised about the management of server " << serverNumber << std::endl;
		
		code = DCT_RC_UNDEFINED_ERROR;
		bbclient->finishServer(serverNumber, code);
	}
	
	return code;
}




//finalSolSize is an input/output argument. On input, you must pass the size (in bytes) of finalSol array (It can be zero). On output, the size in bytes if finalSol need be allocated. 
int DCT_BBClient::treatServerRequests( DCT_Int32 connectCode, unsigned int serverNumber, double* &auxSol, DCT_FinalResults *finalResults, double* &finalSol, DCT_UInt32 &nvarsFinalSol, std::ostream &out)
{
	double lbread = INFINITY;
	DCT_ClientSocket &socket = sockets[serverNumber];
	
	
	if(connectCode == DCT_CC_LOWER_BOUND_RESPONSE)
	{
		
		const int r = socket.readData(&lbread, sizeof(lbread));
		DCT_OSIFERRORRETURN(out, r, r);
		
	}
	else if(connectCode == DCT_CC_NEW_SOLUTION_FOUND)
	{
		int r;
		unsigned int n;
		double nvars;
		double objValue;
		double *sol;
		
		/* Here, we read:
		 * 		1 - The number n of variables in the problem (ATTENTION: to turn easier to client, we read this value as DOUBLE)
		 * 		2 - The lower bound of the problem (double)
		 * 		3 - The objective value in the new solution (double)
		 * 		4 - n double values specifing the solution (double)
		 */
		
		r = socket.readData( &nvars, sizeof(nvars) );
		DCT_OSIFERRORRETURN(out, r, r);
		
		n = nvars;
		if(auxSol == NULL)
		{
			DCT_malloc(auxSol, n+2);
			DCT_OSIFMEMERRORRETURN(out, !auxSol);
		}
		
		/*if(arrayBestObjAndSol == NULL)
		{
			this->nvars = n;
			DCT_malloc(arrayBestObjAndSol, n+1);
			DCT_OSIFMEMERRORRETURN(out, !arrayBestObjAndSol);
			
			bestSol = &arrayBestObjAndSol[1];
		}*/
		
		
		r = socket.readData( auxSol, (n+2)*sizeof(auxSol[0]) );
		DCT_OSIFERRORRETURN(out, r, r);
		
		
		lbread = auxSol[0];
		objValue = auxSol[1];
		sol = &auxSol[2];
		
		r = tryUpdateBestSol(serverNumber, n, objValue, sol, NULL, out);
		DCT_OSIFERRORRETURN(out, r, r);
		
		
	}
	else if(connectCode == DCT_CC_BEST_SOLUTION_REQUEST)
	{
		const bool hasSol = !std::isinf(bestObjValue);
		int r;
		const DCT_Int32 ccode = DCT_CC_NO_SOLUTION_AVAILABLE;
		
		
		SEMAPH_socketsWrite[serverNumber].lock();
		{
			if( hasSol )
			{
				SEMAPH_bestSol.lock();
				{
					arrayBestObjAndSol[0] = bestObjValue;
					
					#if DCT_DEBUG_MODE
						assert( bestSol == &arrayBestObjAndSol[1] );
					#endif
					
					r = sendSolution(socket, nvars, arrayBestObjAndSol);
				}
				SEMAPH_bestSol.unlock();
			}
			else
				r = socket.writeData(&ccode, sizeof(ccode));
		}
		SEMAPH_socketsWrite[serverNumber].unlock();
		
		DCT_OSIFERRORRETURN(out, r, r);
		
	}
	else if(connectCode == DCT_CC_FINAL_RESULTS)
	{
		int r = DCT_readFinalResultsAndSol(socket, *finalResults, nvarsFinalSol, finalSol, &nvarsFinalSol, out);
		DCT_OSIFERRORRETURN(out, r, r);
		
	}
	else
	{
		int r;
		size_t bytesRead = 0;
		const int sizeData = 10000;
		DCT_Byte data[sizeData];
		
		out << DCT_PREPRINT  "An unexpected connection code was read: " << connectCode << DCT_GETFILELINE << "\n";
		
		socket.disableBlock();
		
		r = socket.readData(data, sizeData, &bytesRead);
		
		out << "r: " << r << " bytesRead: " << bytesRead << "\n";
		
		socket.enableBlock();
		
		return DCT_RC_VALUE_ERROR;
	}
	
	
	if(!std::isinf(lbread))
		serverLowerBounds[serverNumber] = lbread;
	
	
	return 0;
}



int DCT_BBClient::run(DCT_Components *servers, 	DCT_SeveralNumberOfComponents *allNumberOfServers, DCT_FileNames *inputFiles, DCT_UInt64 nBasicInputParameters, DCT_Byte *basicInputParameters,  DCT_AllGeneralParams *allGeneralParams )
{
	const double timeStart = DCT_getTime();
	const unsigned int maxnservers = servers->size() + allNumberOfServers->getMaxNumberOfServers();
	
	int code;
	//unsigned int MAX_SIMULTANEOUS_SERVERS_TO_SEND_BASIC_DATA = 1;
	
	
	double *auxSol = NULL;
	std::thread *threads = NULL;
	std::thread *threadsToSendBasicData = NULL;
	//DCT_FinalResults *serverFinalResults = NULL;
	
	
	
	deallocate();
	initialize();
	
	
	nConnectedServers = 0;
	
	if(inputFiles)
	{
		//checking if input files exist
		for(auto &pairInputFile : *inputFiles)
		{
			std::string &inputFile = pairInputFile.second;
			
			FILE *read = fopen( inputFile.c_str(), "r" );
			DCT_IFERRORGOTOLABEL(!read, code, DCT_RC_NAME_ERROR, termination);
			
			fclose(read);
		}
	}
	
	
	
	sockets.reserve(maxnservers);
	
	//first, we just conect with "regular" list of servers, after we think how connect with allNumberOfServers
	
	{
		
		for(auto server: *servers)
		{
			unsigned int nthreads;
			DCT_ClientSocket &socket = sockets[ nConnectedServers ];
			
			std::cout << DCT_PREPRINT "Trying to connect to server " << server->ipaddress << " on port " << server->port << std::endl; //we use endl to perform a flush
			
			int r = socket.connectToAServer( server->port, server->ipaddress );
			
			if(r != 0)
			{
				DCT_PRINTERRORNUMBER(r);
				
				std::cerr << DCT_PREPRINT << "Error to connect with server " << server->ipaddress << " on port " << server->port << std::endl; //we use endl to perform a flush
				
				continue;
			}
			
			std::cout << DCT_PREPRINT "Requesting service to server " << server->ipaddress << " on port " << server->port << std::endl;
			
			r = getService(socket, server->maxThreads, nthreads);
			if(r != 0)
			{
				DCT_PRINTERRORNUMBER(r);
				continue;
			}
			
			nConnectedServers++;
		}
		
		//TODO: treat here allNumberOfServers
	}
	
	
	
	if( nConnectedServers == 0 )
	{
		DCT_PRINTERRORMSG("Sorry! No available server to run your application");
		
		code = DCT_RC_NO_SERVER_AVAILABLE;
		goto termination;
	}
	
	
	//sockets.resize(connectedServers); resize does not work if you do not have an explicity copy constructor 
	
	
	threadsToSendBasicData = new (std::nothrow) std::thread[nConnectedServers-1];
	threads = new (std::nothrow) std::thread[nConnectedServers];
	//serverFinalResults = new (std::nothrow) DCT_FinalResults[nConnectedServers];
	SEMAPH_socketsWrite = new (std::nothrow) DCT_Mutex[nConnectedServers];
	//SEMAPH_varBoundsQueue = new (std::nothrow) DCT_Mutex[nConnectedServers];
	varBoundsRecepQueue = new (std::nothrow) DCT_VarBoundsReceptorQueue[nConnectedServers];
	
	DCT_calloc(serverRunning, nConnectedServers);
	DCT_calloc(threadReturnCodes, nConnectedServers);
	DCT_malloc(serverLowerBounds, nConnectedServers);
	DCT_calloc(serverPVarsBounds, nConnectedServers);
	
	
	DCT_IFMEMERRORGOTOLABEL(!threadsToSendBasicData || !threads || !SEMAPH_socketsWrite || !varBoundsRecepQueue || !serverRunning || !serverLowerBounds || !serverPVarsBounds, code, termination );
	
	
	try
	{
		initialOpenNodes.reserve(nConnectedServers);
	}
	catch(std::bad_alloc& ba)
	{
		DCT_PRINTMEMERROR;
		code = DCT_RC_MEMORY_ERROR;
		goto termination;
	}
	
	
	DCT_setAllArray<double>(nConnectedServers, serverLowerBounds, INFINITY);
	
	for(unsigned int i = 0; i < nConnectedServers; i++)
		varBoundsRecepQueue[i].secondsToSleep = secondsToSleepWaitingANode;
	
	
	finalResults.initialize();
	
	
	//first, we run the first server, to generate the first partitions
	serverRunning[0] = 1;
	someServerEnded = false;
	
	
	//threadsToSendBasicData[0] = std::thread(DCT_semaphSendBasicDataToServer, this, 0, &SEMAPH_socketsWrite[0], &sockets[0], inputFiles, nBasicInputParameters, basicInputParameters, allGeneralParams, nConnectedServers - 1);
	
	{
		//sending basic data to first server
		int r = DCT_semaphSendBasicDataToServer(this, 0, &SEMAPH_socketsWrite[0], &sockets[0], inputFiles, nBasicInputParameters, basicInputParameters, allGeneralParams, nConnectedServers - 1);
		
		DCT_IFERRORGOTOLABEL(r, code, r, termination);
	}
	
	
	//we start the first thread to start explorartion and send the first nodes
	threads[0] = std::thread(DCT_clientThread, this, 0, &sockets[0], inputFiles, nBasicInputParameters, basicInputParameters, allGeneralParams, nConnectedServers - 1, (DCT_Byte *) NULL );
	
	
	//sending basic data to other servers. Note, we can send the datas before we get the initial nodes to other serves
	for(decltype(nConnectedServers) i = 1; i < nConnectedServers; i++)
	{
		//to avoid have problems about connections, we just send data to 10 servers by time
		if( (i-1) % maxSimultaneuousServersToSendBasicData == 0 )
		{
			for( decltype(i) j = 1; j < i; j++)
			{
				if( threadsToSendBasicData[j-1].joinable() ) //std::thread can be joined just one time. So, if we already have call join for this thread, we cannot call join again, otherwise, we get an execution error. 
					threadsToSendBasicData[j-1].join();
			}
		}
		
		threadsToSendBasicData[i-1] = std::thread(DCT_semaphSendBasicDataToServer, this, i, &SEMAPH_socketsWrite[i], &sockets[i], inputFiles, nBasicInputParameters, basicInputParameters, allGeneralParams, 0);
		
	}
	
	
	for(unsigned int i = 1; i < nConnectedServers; i++)
	{
		bool continueWait = true;
		DCT_Byte *initialOpenNodeRep = NULL;
		
		
		if( threadsToSendBasicData[i-1].joinable() ) //std::thread can be joined just one time. So, if we already have call join for this thread, we cannot call join again, otherwise, we get an execution error. If joinable returns false, is because we already call join for this thread, so, we are sure the basic data were already sent
			threadsToSendBasicData[i-1].join(); //making sure we already sent the basic data to this server
		
		
		if( threadReturnCodes[i] != 0 )
		{
			//DCT_semaphSendBasicDataToServer got some failure to send basic data. So, we cannot start service in this server
			continue;
		}
		
		
		while(continueWait)
		{
			//we must check if server 0 finished before check if initialOpenNodeRep has some node to this server. So, we guarantee that we will get a possible initialNode stored to this server. If we revert and check initialOpenNodes before check if server 0 Finished, we can have a bad luck situation where server 0 put a open node to this server after we check initialOpenNodes and, immediately finished before we check server0FinishedSomeExecution
			if( server0FinishedSomeExecution )
			{
				continueWait = false;//so server 0 stop its execution without provide all initial open nodes requested. We let the remaining server start and request open nodes to some other server. //Do not break here because before stop, server 0 can have put some open node
			}
			
			SEMAPH_initialOpenNodes.lock();
			{
				if( initialOpenNodes.size() >= i )
				{
					initialOpenNodeRep = initialOpenNodes[i-1];
					continueWait = false;
				}
			}
			SEMAPH_initialOpenNodes.unlock();
			
		}
		
		//std::cout << "initialOpenNodeRep: " << (!initialOpenNodeRep ? "not null" : "NULL") << "\n";
		//DCT_getchar();
		
		if(initialOpenNodeRep)
		{
			serverRunning[i] = 1; //if serverRunning[i] is 1, client wo not request an open node to another server to send to server i, and server i will explorer the represented by openNodeRep
		}
		
		threads[i] = std::thread(DCT_clientThread, this, i, &sockets[i], inputFiles, nBasicInputParameters, basicInputParameters, allGeneralParams, 0, initialOpenNodeRep );
	}
	
	
	DCT_secDeleteArray(threadsToSendBasicData);
	initialOpenNodes.clear(); //we can erase vector initialOpenNodes
	
	for(unsigned int i = 0; i < nConnectedServers; i++)
		threads[i].join();
	
	
	if( bestSol != NULL )
	{
		unsigned int serverErrorSize = serverError.size();
		
		finalResults.feasSolFound = true;
		finalResults.upperBound = bestObjValue;
		finalResults.objBestSol = bestObjValue;
		
		//check if we have a false infeasible problem status in serverError. It can happen if some server finish a execution before the first feasible solution be found. The server will flag infeasible problem and this flag will be put in the server error if we do not have feasible solution yet. But, after, if some other service found a feasible solution, we have to remove infeasible problem status from server error
		
		for(unsigned int i = serverErrorSize; i > 0; i--)
		{
			unsigned int ind = i-1;
			
			if( serverError[ind].error == DCT_ORC_INFEASIBLE_PROBLEM)
			{ //so, we shift the serverError array to erase this error, since we have a feasible solution and the ret code is infeasible problem.
				for(unsigned int k = ind; k < serverErrorSize-1; k++)
					serverError[k] = serverError[k+1];
				
				serverErrorSize--;
			}
		}
		
		if( serverErrorSize < serverError.size() )
			serverError.resize(serverErrorSize);
	}
	
	
	
	{
		//check if some node got unexplored due to some error in the connection with server (servers can go down or be disconnected) or some other error in the client
		
		//check if some node got unexplored in serverPVarsBounds, i.e., some connection error with the server
		for(unsigned int i = 0; i < nConnectedServers; i++)
		{
			if( serverPVarsBounds[i] != NULL )
			{ //the node gots unexplored
				DCT_ServerError se;
				
				se.server = i;
				se.error = DCT_ORC_SERVER_CONNECTION_ERROR;
				
				serverError.push_back( se ); //here, we only have one thread. So, we do not need semaphore
				
				std::cerr << DCT_PREPRINT << "Error: An open node got unexplored in the server " << i << "\n";
			}
		}
		
		//check if some node got unexplored in the server queues, i.e., client error
		for(unsigned int i = 0; i < nConnectedServers; i++)
		{
			bool error = false;
			
			std::queue<DCT_Byte*> &queue = varBoundsRecepQueue[i].queue;
			
			while( !queue.empty() )
			{
				auto p = queue.front();
				queue.pop();
				if(p != NULL)
				{ //so, we have some node unexplored
					delete p;
					error = true;
				}
			}
			
			if(error)
			{
				DCT_ServerError se;
				
				se.server = i;
				se.error = DCT_ORC_CLIENT_ERROR;
				
				serverError.push_back( se ); //here, we only have one thread. So, we do not need semaphore
				
				std::cerr << DCT_PREPRINT "Error: Open nodes got unexplored in the server queue " << i << "\n";
			}
		}
	}
	
	
	
	
	if( serverError.size() == 0  )
	{
		//we did not error in any server
		finalResults.optCode = DCT_ORC_OPTIMAL_SOLUTION;
	}
	else
	{
		//if all server error are equals, we set if
		const unsigned int serverErrorSize = serverError.size();
		bool difError = false;
		
		for(unsigned int i = 1; i < serverErrorSize; i++)
		{
			if( serverError[i].error != serverError[0].error )
			{
				difError = true;
				break;
			}
		}
		
		if( !difError)
		{
			finalResults.optCode = serverError[0].error;
		}
	}
	
	
	//TODO: some day we need to change that when we really work with more servers
	//finalResults = serverFinalResults[0];
	
	
	
	code = 0;
	
termination:
	
	if(auxSol)		free(auxSol);
	if(threadsToSendBasicData)		delete[] threadsToSendBasicData;
	if(threads)		delete[] threads;
	//if(serverFinalResults)	delete[] serverFinalResults;
	//if(serverVarBounds)		delete[] serverVarBounds;
	
	
	DCT_secFree(serverRunning);
	DCT_secFree(threadReturnCodes);
	DCT_secFree(serverLowerBounds);
	DCT_secDeleteArray(SEMAPH_socketsWrite);
	//DCT_secDeleteArray(SEMAPH_varBoundsQueue);
	DCT_secDeleteArray(varBoundsRecepQueue);
	
	if(serverPVarsBounds)
	{
		for(unsigned int i = 0; i < nConnectedServers; i++)
		{
			if(serverPVarsBounds[i])
				delete serverPVarsBounds[i];
		}
		
		free(serverPVarsBounds);
		      serverPVarsBounds = NULL;
	}
	
	sockets.clear();
	
	finalResults.clockTime = DCT_getTime() - timeStart;
	
	return code;
}



DCT_Mutex SEMAPH_clientSetToSigIntHandler;
std::set<DCT_BBClient *> clientSetToSigIntHandler;


static void DCT_clientSigIntHandler(int signal)
{
	std::cerr << DCT_PREPRINT "Signal " << signal << " was received\n";
	
	SEMAPH_clientSetToSigIntHandler.lock();
	{
		for( auto &client : clientSetToSigIntHandler )
			client->stopOperations();
	}
	SEMAPH_clientSetToSigIntHandler.unlock();
	
	exit(-1);
}


int DCT_BBClient::setSigIntHandler()
{
	int code = 0;
	
	//we protect using a Semaphore because we can have a situation where a user sets several threads having DCT_BBClient differents. So, we protect the access
	SEMAPH_clientSetToSigIntHandler.lock();
	{
		#ifdef SIGINT
			signal(SIGINT, DCT_clientSigIntHandler);
		#endif
		#ifdef SIGQUIT
			signal(SIGQUIT, DCT_clientSigIntHandler);
		#endif
		#ifdef SIGTERM
			signal(SIGTERM, DCT_clientSigIntHandler);
		#endif
		#ifdef SIGABRT
			signal(SIGABRT, DCT_clientSigIntHandler);
		#endif
		#ifdef SIGHUP
			signal(SIGHUP, DCT_clientSigIntHandler);
		#endif
		
			
		try
		{
			clientSetToSigIntHandler.insert(this);
		}
		catch(std::bad_alloc &ba)
		{
			code = DCT_RC_MEMORY_ERROR;
		}
	}
	SEMAPH_clientSetToSigIntHandler.unlock();
	
	return code;
}


//function to iterrupt and stop the services in progress
void DCT_BBClient::stopOperations()
{
	DCT_Int32 request = DCT_CC_STOP_SERVICE;
	int r;
	
	for(unsigned int i = 0; i < nConnectedServers; i++)
	{
		std::cout << "Encerrando servidor " << i << "\n";
		
		SEMAPH_socketsWrite[i].lock();
		{
			r = sockets[i].writeData(&request, sizeof(request));
		}
		SEMAPH_socketsWrite[i].unlock();
		if(r != 0)
			DCT_PRINTERRORNUMBER(r);
			
		sockets[i].closeSocket();
	}
}


int DCT_BBClient::unsetSigIntHandler()
{
	SEMAPH_clientSetToSigIntHandler.lock();
	{
		if( clientSetToSigIntHandler.count(this) > 0 )
			clientSetToSigIntHandler.erase(this);
		
		if( clientSetToSigIntHandler.size() == 0 )
		{
			//we unset our handlers
			#ifdef SIGINT
				signal(SIGINT, SIG_DFL);
			#endif
			#ifdef SIGQUIT
				signal(SIGQUIT, SIG_DFL);
			#endif
			#ifdef SIGTERM
				signal(SIGTERM, SIG_DFL);
			#endif
			#ifdef SIGABRT
				signal(SIGABRT, SIG_DFL);
			#endif
			#ifdef SIGHUP
				signal(SIGHUP, SIG_DFL);
			#endif
		}
	}
	SEMAPH_clientSetToSigIntHandler.unlock();
	
	return 0;
}


int DCT_BBClient::clientThread(const unsigned int serverNumber, DCT_ClientSocket &socket, DCT_FileNames *inputFiles, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_AllGeneralParams *allGeneralParams, DCT_UInt32 nInitialRequestedNodes, DCT_Byte *pvarsBounds)
{
	const unsigned int nservers = nConnectedServers;
	bool someServerRunning = true; //we need start this as true to work with first server that already begin with its serverRunning set as true
	int r, code = 0;
	unsigned int nvarsFinalSol = 0;
	unsigned int serverllb = 0; //we initialize for the first server
	DCT_Int32 responseCode;
	
	double *auxSol = NULL;
	double *finalSol = NULL;
	
	std::ofstream myfile;
	
	std::ostream *sout = &std::cout;
	std::ostream *serr = &std::cerr;
	
	
	//DCT_Byte *pvarsBounds = NULL;
	//DCT_VarBoundsWriter varBoundsWriter;
	DCT_FinalResults myFinalResults;
	
	
	#if DCT_SAVE_CLIENT_SERVER_LOGS
	{
		char fname[50];
		
		sprintf(fname, DCT_CLIENT_SERVER_LOGS_FILE_NAME_PREFIX "%03u.log", serverNumber);
		
		try
		{
			myfile.open(fname);
		}
		catch(std::exception &e)
		{
		}
		
		if( myfile.is_open() )
		{
			sout = &myfile;
			serr = &myfile;
		}
		else
		{
			DCT_OSPRINTERRORMSGP(*serr, "Failure to open file ", fname);
		}
	}
	#endif
	
	
	{
		//now, to save time, we sent the basic data before enter in this function, i. e., before we have the open initial node to server
		
		/*SEMAPH_socketsWrite[serverNumber].lock();
		{
			r = sendBasicDataToServer(socket, inputFiles, nBasicInputParameters, basicInputParameters, allGeneralParams, nInitialRequestedNodes);
		}
		SEMAPH_socketsWrite[serverNumber].unlock(); */
		
		
		/*r = DCT_semaphSendBasicDataToServer(this, serverNumber, &SEMAPH_socketsWrite[serverNumber], &socket, inputFiles, nBasicInputParameters, basicInputParameters, allGeneralParams, nInitialRequestedNodes);
		
		DCT_OSIFERRORGOTOLABEL(*serr, r, code, r, termination); */
	} 
	
	while(true)
	{
		if( serverRunning[serverNumber] <= 0 )
		{
			//is the server is not running, we must request a open node from some other server to make this server run.
			DCT_ClientSocket *socketllb;
			const DCT_Int32 requestCode = DCT_CC_OPEN_NODE_REQUEST;
			
			pvarsBounds = NULL;
			
			//first, we check if some server go down (for example, it was disconnected) and let  some open node
			if( someServerEnded )
			{
				SEMAPH_checkOpenPVarBounds.lock();
				{
					for(unsigned int i = 0; i < nservers; i++)
					{
						if( serverRunning[i] < 0 && serverPVarsBounds[i] != NULL ) //so, server i already finished, but it has let some open node 
						{
							pvarsBounds = serverPVarsBounds[i];
							serverPVarsBounds[i] = NULL;
							#if DCT_DEBUG_MODE
								assert(i != serverNumber);
								assert(someServerRunning == true);
							#endif
							break;
						}
					}
				}
				SEMAPH_checkOpenPVarBounds.unlock();
			}
			
			
			while( pvarsBounds == NULL )
			{
				#if DCT_DEBUG_MODE
					*sout << DCT_PREPRINT "server: " << serverNumber << " serverRunning: ";
					for(unsigned int i = 0; i < nservers; i++)
						*sout << (int) serverRunning[i] << " ";
					*sout << "\n";
				#endif
				
				
				//so, we must check all other servers to find a server running and request an open node to it. Note, only the first server will start this function with serverRunning == true
				
				
				//we choose a running server with the lowest lower bound to request a open node.
				
				serverllb = -1;
				double llb = INFINITY;
				
				someServerRunning = false;
				
				for(unsigned int i = 0; i < nservers; i++)
				{
					if( serverRunning[i] > 0 )
					{
						serverRunning[serverNumber] = 0; //by safe, we put zero here to avoid all threads give up since serverRunning[serverNumber] can be -1 (when all serverRunnings are -1, all proccess end)
						
						someServerRunning = true;
						if( serverLowerBounds[i] <= llb )
						{
							serverllb = i;
							llb = serverLowerBounds[i];
						}
					}
				}
				
				if(someServerRunning == false)
				{
					//we just abort if all serverRunning are -1
					
					if( serverRunning[serverNumber] == 0 )
					{
						serverRunning[serverNumber] = -1;
						//to decrease the chance of a bizarr error, we put this thread to sleep before check again if all servers are not runnings. So, maybe some server have time to change it status to run. I know it a terrible solution (actually it nos even a solution), but that is the best that we can do by now. It decreases client charge of proccess also because thread will not be check all time until the last server finish...
						
						std::this_thread::sleep_for(std::chrono::seconds( DCT_SECONDS_TO_SLEEP_BEFORE_CHECK_AGAIN_ALL_SERVER_END));
						
						continue;
					}
					else
					{
						#if DCT_DEBUG_MODE
							assert( serverRunning[serverNumber] == -1 );
						#endif
						
						varBoundsRecepQueue[serverNumber].disable();
						
						//we check to see if all serverRunning are -1
						bool allEnd = true;
						
						for(unsigned int i = 0; i < nservers; i++)
						{
							if( serverRunning[i] >= 0 )
							{
								allEnd = false;
								break;
							}
						}
						
						if(allEnd)
							break; //so, the work is finished
						else
							continue;
					}
				}
				
				
				//we have at least one server running. So, we request a open node
				
				#if DCT_DEBUG_MODE
					assert(serverllb >= 0 && serverllb < nservers);
					assert(serverllb != serverNumber);
				#endif
				
				socketllb = &sockets[serverllb];
				
				while( pvarsBounds == NULL )
				{
					if( serverRunning[serverllb] <= 0 )
					{
						//serverllb is not running more
						break;
					}
					
					*sout << DCT_PREPRINT "server " << serverNumber << " has requested a open node to server " << serverllb << std::endl;
					
					SEMAPH_socketsWrite[serverllb].lock();
					{
						//we make a request to serverllb, but only thread serverllb can read the data received. So, we will be monitoring the queue of varBounds of that thread to see our varBounds list
						r = socketllb->writeData( &requestCode, sizeof(requestCode) );
					}
					SEMAPH_socketsWrite[serverllb].unlock();
					if( r!= 0 )
					{//we prefer do not terminate this thread if some error had occurred because other thread could raise some exception and have finished. So, we just print the error and pray for thigns go fine anyway :)
						DCT_OSPRINTERRORNUMBER(*serr, r);
					}
					//DCT_OSIFERRORGOTOLABEL(*serr, r, code, r, termination);
					
					
					r = varBoundsRecepQueue[serverllb].pop(pvarsBounds); //the pop method will put this thread to sleep until another thread put some data in the queue to be popped. This is curious, because this thread makes the request for a node, but, other thread will receive the node (the thread wich is managing the the server serverllb)
					DCT_OSIFERRORGOTOLABEL(*serr, r, code, r, termination);
				} //end of while( pvarsBounds == NULL )
				
			} //end of while( pvarsBounds == NULL )
			
			
		} //end of if( serverRunning[serverNumber] <= 0 )
		#if DCT_DEBUG_MODE
		else
		{
			assert( pvarsBounds != NULL || serverNumber == 0 );
		}
		#endif
		
		if(someServerRunning == false)
		{
			//so, the work is finished because we do not have any server running
			break;
		}
		
		*sout << DCT_PREPRINT "receiving open node from server " << serverllb << " to send to server " << serverNumber << std::endl; //we use endl to perform a flush
				
				
		//we must set the server is running here to avoid some other thread stop the procedure because all serverRunnings are marked as false
				
		serverRunning[serverNumber] = 1;
		varBoundsRecepQueue[serverNumber].enable();
		
		
		
		serverPVarsBounds[serverNumber] = pvarsBounds;
		
		
		SEMAPH_socketsWrite[serverNumber].lock();
		{
			r = startService(socket, pvarsBounds, serverLowerBounds[serverllb] );
		}
		SEMAPH_socketsWrite[serverNumber].unlock();
		
		if(r != 0)
		{
			DCT_OSPRINTERRORNUMBER(*serr, r);
			code = r;
			break;
		}
		
		varBoundsRecepQueue[serverNumber].enable();
		//serverLowerBounds[serverNumber] = serverLowerBounds[serverllb]; if we put this server with the same lower bound than serverllb, we can get a bad sittuation where servers request nodes to serverNumber instead of serverllb. So, to avoid that, we let serverLowerBounds[serverNumber] as INFINITY until servNumber sent us its lower bound. 
		
		//to treat messages from the server
		do
		{
			r = socket.readData( &responseCode, sizeof(responseCode) );
			if(r != 0)
			{
				DCT_OSPRINTERRORNUMBER(*serr, r);
				std::cerr << DCT_PREPRINT "server: " << serverNumber << " failure to read from sokect. Error: " << r << "\n";
				code = DCT_RC_SERVER_ERROR; //r;
				goto termination;
			}
			
			{
				std::string strCode;
				
				DCT_enumToStr( (DCT_CONECTION_CODE) responseCode, strCode );
				
				std::cout << DCT_PREPRINT "server: " << serverNumber << " lb: " << serverLowerBounds[serverNumber] << " ub: " << bestObjValue << " response code: " << responseCode << " (" << strCode << ")" << " node queue size: " << varBoundsRecepQueue[serverNumber].queue.size() << "\n";
				
			}
			
			//if( responseCode == DCT_VAR_BOUNDS || responseCode == DCT_NO_OPEN_NODE_AVAILABLE || responseCode == DCT_SERVICE_NOT_STARTED  )
			if( responseCode == DCT_CC_OPEN_NODE_RESPONSE || responseCode == DCT_CC_NO_OPEN_NODE_AVAILABLE || responseCode == DCT_CC_NO_OPEN_NODE_SERVICE_NOT_STARTED || responseCode == DCT_CC_INITIAL_OPEN_NODE_RESPONSE )
			{
				//we just treat DCT_SERVICE_NOT_STARTED because another thread could have asked a open node when this server was not running
				
				DCT_Byte *openNodeRep = NULL;
				
				if( responseCode == DCT_CC_OPEN_NODE_RESPONSE || responseCode == DCT_CC_INITIAL_OPEN_NODE_RESPONSE )
				{
					/* Here, we read:
					 * 	1 - The size m of open nodes being sent (uint64)
					 *  2 - m bytes representing the nodes. We just pass this to another server and let this server responsible to interpret this m bytes.
					*/
					DCT_Int32 ccode = DCT_CC_OPEN_NODE_RESPONSE;
					DCT_UInt64 nBytesOfNodesSent;
					DCT_UInt64 sizeOpenNodeRep;
					DCT_Byte *paux;
					
					r = socket.readData( &nBytesOfNodesSent, sizeof(nBytesOfNodesSent) );
					if(r != 0 )
					{
						DCT_OSPRINTERRORNUMBER(*serr, r);
						code = r;
						goto termination;
					}
					
					
					sizeOpenNodeRep = sizeof(DCT_Int32) + sizeof(DCT_UInt64) + nBytesOfNodesSent;
					
					DCT_malloc(openNodeRep, sizeOpenNodeRep);
					if( !openNodeRep )
					{
						DCT_OSPRINTMEMERROR(*serr);
						code = DCT_RC_MEMORY_ERROR;
						goto termination;
					}
					
					paux = openNodeRep;
					
					DCT_writeAndShift(paux, ccode);
					DCT_writeAndShift(paux, nBytesOfNodesSent);
					
					r = socket.readData(paux, nBytesOfNodesSent);
					if(r != 0 )
					{
						DCT_OSPRINTERRORNUMBER(*serr, r);
						code = r;
						goto termination;
					}
					
					#if DCT_DEBUG_MODE
						assert(&openNodeRep[sizeOpenNodeRep-1] == &paux[nBytesOfNodesSent-1]);
					#endif
					
					
					/*double lb;
					unsigned int sizeOpenNodeRep, sizeBytesBounds;
					DCT_UInt32 nbounds;
					DCT_Byte *paux;
					
					r = socket.readData(&lb, sizeof(lb) );
					if(r != 0 )
					{
						DCT_OSPRINTERRORNUMBER(*serr, r);
						code = r;
						goto termination;
					}
					
					r = socket.readData(&nbounds, sizeof(nbounds) );
					if(r != 0 )
					{
						DCT_OSPRINTERRORNUMBER(*serr, r);
						code = r;
						goto termination;
					}
					
					
					sizeBytesBounds = nbounds*( sizeof(DCT_UInt32) + 2*sizeof(double) );
					sizeOpenNodeRep = DCT_sizeofOpenNodeRep(nbounds);
					
					DCT_malloc(openNodeRep, sizeOpenNodeRep);
					if( !openNodeRep )
					{
						DCT_OSPRINTMEMERROR(*serr);
						code = DCT_RC_MEMORY_ERROR;
						goto termination;
					}
					
					paux = openNodeRep;
					
					//writing the open node representation
					DCT_writeAndShift(paux, DCT_CC_OPEN_NODE_RESPONSE);
					DCT_writeAndShift(paux, lb);
					DCT_writeAndShift(paux, nbounds);
					
					r = socket.readData(paux, sizeBytesBounds );
					if(r != 0 )
					{
						free(openNodeRep);
						
						DCT_OSPRINTERRORNUMBER(*serr, r);
						code = r;
						goto termination;
					}
					
					std::cout << "nbounds: " << nbounds << "\n";
					
					#if DCT_DEBUG_MODE
						if(sizeBytesBounds > 0)
							assert( &paux[sizeBytesBounds-1] == &openNodeRep[sizeOpenNodeRep-1]); 
					#endif */
					
				}
				
				/*r = 0;
				SEMAPH_varBoundsQueue[serverNumber].lock();
				{
					try
					{
						varBoundsQueue[serverNumber].push(vb);
					}
					catch(std::bad_alloc &ba)
					{
						r = DCT_MEMORY_ERROR;
					}
				}
				SEMAPH_varBoundsQueue[serverNumber].unlock(); */
				
				//std::cout << "Vou mexer no receptor de nos" << std::endl; 
				
				if( responseCode == DCT_CC_INITIAL_OPEN_NODE_RESPONSE )
				{
					#if DCT_DEBUG_MODE
						assert(serverNumber == 0);
					#endif
					
					SEMAPH_initialOpenNodes.lock();
					{
						initialOpenNodes.push_back(openNodeRep);
					}
					SEMAPH_initialOpenNodes.unlock();
				}
				else
				{
					r = varBoundsRecepQueue[serverNumber].push(openNodeRep);
				}
				//std::cout << "Mexi no receptor de nos" << std::endl;
				
				if(r != 0)
				{
					if(openNodeRep)  free(openNodeRep);
					DCT_OSPRINTMEMERROR(*serr);
					code = r;
					goto termination;
				}
			}
			else
			{
				r = treatServerRequests( responseCode, serverNumber, auxSol, &myFinalResults, finalSol, nvarsFinalSol, *sout );
				
				DCT_OSIFERRORGOTOLABEL(*serr, r, code, r, termination);
				
			}
			
			
		}while(responseCode != DCT_CC_FINAL_RESULTS);
		
		if(serverNumber == 0)
			server0FinishedSomeExecution = true;
		
		//we should not disable varBoundsRecepQueue here because we can have a open request to open node to be proccessed by the serverNumber. In this case, disable the varBoundsRecepQueue[serverNumber] will generate a bug since this thread will put a pointer in the queue, but, since the queue was disable, the waiting thread will not more hoping get a node from this queue. So, ONLY DISABLE varBoundsRecepQueue[serverNumber] when you are sure the server is not working more.
		//varBoundsRecepQueue[serverNumber].disable();
		
		serverLowerBounds[serverNumber] = INFINITY; //we put infinity here to avoid others servers request lower bound to this node
		
		//serves finalizes exploration of node. So, we can free its varBounds
		DCT_secFree(serverPVarsBounds[serverNumber]); //note: we put NULL or serverPVarsBounds to sinalize which that node is explored. So, if the server ends execution and serverPVarsBounds is not NULL. it is because we had some problem about the server and we lost it, probabily the server was disconect or go down. So, we have to give the respective pvarsBounds to other server
		
		
		if( myFinalResults.feasSolFound )
		{
			tryUpdateBestSol(serverNumber, nvarsFinalSol, myFinalResults.objBestSol, finalSol, NULL, *sout);
		}
		
		
		if( myFinalResults.optCode == DCT_ORC_OPTIMAL_SOLUTION )
			someOptimalRetCode = true;
		
		//in this point, we reach the DCT_FINAL_RESULTS connection code
		if( (myFinalResults.optCode != DCT_ORC_OPTIMAL_SOLUTION && myFinalResults.optCode != DCT_ORC_INFEASIBLE_PROBLEM) || (myFinalResults.optCode == DCT_ORC_INFEASIBLE_PROBLEM && someOptimalRetCode == false ) ) //actually even if  optCode is DCT_ORC_INFEASIBLE_PROBLEM and someOptimalRetCode is false, it is not necessarelly an error, since another server could still find a feasible solution. Even so, we flag the infeasibility and check after all threads have finished.
		{
			DCT_ServerError se;
			
			se.server = serverNumber;
			se.error = myFinalResults.optCode;
			
			SEMAPH_serverError.lock();
			{
				serverError.push_back( se );
			}
			SEMAPH_serverError.unlock();
		}
		
		*sout << DCT_PREPRINT "server: " << serverNumber << " myFinalResults.algorithm: " << myFinalResults.algorithm << "\n";
		*sout << DCT_PREPRINT "server: " << serverNumber << " myFinalResults.clockTime: " << myFinalResults.clockTime << "\n";
		*sout << DCT_PREPRINT "server: " << serverNumber << " myFinalResults.cpuTime: " << myFinalResults.cpuTime << "\n";
		*sout << DCT_PREPRINT "server: " << serverNumber << " myFinalResults.nIters: " << myFinalResults.nIters << "\n";
		*sout << DCT_PREPRINT "server: " << serverNumber << " myFinalResults.feasSolFound: " << (int) myFinalResults.feasSolFound << "\n";
		*sout << DCT_PREPRINT "server: " << serverNumber << " myFinalResults.lowerBound: " << myFinalResults.lowerBound << std::endl;
		//std::cout << ": " <<  << "\n";
		//DCT_getchar();
		
		
		SEMPAH_FinalResult.lock();
		{
			finalResults.algorithm = myFinalResults.algorithm;
			//finalResults.clockTime += myFinalResults.clockTime; //we set our finalResults.clockTime by ourselves
			finalResults.cpuTime += myFinalResults.cpuTime;
			finalResults.nIters += myFinalResults.nIters;
			finalResults.feasSolFound = finalResults.feasSolFound || myFinalResults.feasSolFound;
			finalResults.nServerCalls += myFinalResults.nServerCalls + 1;
			if( std::isinf(finalResults.lowerBound) )
			{
				finalResults.lowerBound = myFinalResults.lowerBound;
			}
			else
			{
				if(myFinalResults.lowerBound < finalResults.lowerBound)
					finalResults.lowerBound = myFinalResults.lowerBound; //we adopt the lowest lower bound as the general lower bound
			}
		}
		SEMPAH_FinalResult.unlock();
		
		
		serverRunning[serverNumber] = 0;
		
	}
	
	
	
	
	
	
termination:
	
	
	finishServer(serverNumber, code);
	
	*sout << DCT_PREPRINT << "server " << serverNumber << " has finished works with the code " << code << std::endl;
	
	if(auxSol)		free(auxSol);
	if(finalSol)	free(finalSol);
	
	if(myfile.is_open())
		myfile.close();
	
	return code;
}





int DCT_BBClient::tryUpdateBestSol(const unsigned int serverNumber, const unsigned int n, double objValue, const double *solution, bool *updated, std::ostream &out)
{
	bool myUpdate = false;
	int code;
	
	SEMAPH_bestSol.lock();
	{
		if(objValue < bestObjValue )
		{
			
			if(arrayBestObjAndSol == NULL)
			{
				this->nvars = n;
				DCT_malloc(arrayBestObjAndSol, n+1);
				
				if( !arrayBestObjAndSol )
				{
					SEMAPH_bestSol.unlock();
					DCT_OSIFMEMERRORGOTOLABEL( out, true, code, termination );
				}
				
				bestSol = &arrayBestObjAndSol[1];
			}
			
			bestObjValue = objValue;
			arrayBestObjAndSol[0] = objValue;
			DCT_copyArray(n, solution, bestSol);
			myUpdate = true;
		}
		
	}
	SEMAPH_bestSol.unlock();
	
	if(myUpdate)
	{
		int r = sendSolutionToServers(n, NULL, arrayBestObjAndSol, serverNumber, out);
		DCT_OSIFERRORGOTOLABEL(out, r, code, r, termination);
	}
	
	
	code = 0;
	
termination:

	if(updated)
		*updated = myUpdate;
	
	return code;
}

























