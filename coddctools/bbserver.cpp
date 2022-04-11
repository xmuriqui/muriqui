
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cassert>
#include <cmath>

#include <iostream>
#include <new>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "DCT_bbserver.hpp"
#include "DCT_tools.hpp"


#include <thread>




using namespace dctools;


class DCT_RunningThread
{
	
public:
	
	bool running;
	std::thread *pthread;
	
	DCT_RunningThread()
	{
		running = false;
		pthread = NULL;
	}
};





int DCT_ServiceCore::sendResume(DCT_ServerConnection *connection, DCT_FinalResults &finalResults, DCT_UInt32 nvars, double *solution)
{
	const unsigned int nresponse = 2;
	const DCT_Int32 response[2] = {DCT_CC_FINAL_RESULTS, sizeof(finalResults)};
	int r;
	
	/*
	 * Here, we write:
	 * 
	 * 1 - the connection code DCT_FINAL_RESULTS (int32)
	 * 
	 * 2 - the bytes size of DCT_FinalResults (uint32). (ok, In a first view, it would not be necessary, but may, in future versions, we can have a DCT_FinalResults having more fields, and so, more bytes, in this case, we have to keep compatibility with previous versions)
	 * 
	 * 3 - the DCT_FinalResults data 
	 * 
	 * 4 - Number of variables in solution (uint32) (if user do not want pass a solution, just put zero here)
	 * 
	 * 5 - The solution
	 * 
	 * 
	 */ 
	
	r = connection->writeData( response, nresponse * sizeof(response[0]) );
	DCT_IFERRORRETURN(r, DCT_RC_SERVER_ERROR);
	
	
	r = connection->writeData( &finalResults, sizeof(finalResults) );
	DCT_IFERRORRETURN(r, DCT_RC_SERVER_ERROR);
	
	
	r = connection->writeData( &nvars, sizeof(nvars) );
	DCT_IFERRORRETURN(r, DCT_RC_SERVER_ERROR);
	
	
	if(nvars > 0)
	{
		r = connection->writeData( solution, nvars *sizeof(solution[0]) );
		DCT_IFERRORRETURN(r, DCT_RC_SERVER_ERROR);
	}
	
	return 0;
}



DCT_BBServer::DCT_BBServer(DCT_ServiceCoreGenerator &serviceCoreGenerator)
{
	maxThreads = DCT_getNumCores();
	nthreads = 0; //threads being used
	maxThreadsByService = maxThreads;
	
	nservices = 0; //services being executed
	maxServices = UINT_MAX;
	
	this->serviceCoreGenerator = &serviceCoreGenerator;
}



DCT_BBServer::~DCT_BBServer()
{
}


//a negative value to maxThreads indicates that server must use the number of cores in the hardware minus the absolute value of maxThreads
void DCT_BBServer::setParams(int maxThreads, unsigned int maxThreadsByService, unsigned int maxSimultaneuousServices)
{
	const unsigned int ncores = DCT_getNumCores();
	
	
	if( maxThreads < 0 )  //a negative value to maxThreads indicates that server must use the number of cores in the hardware minus the absolute value of maxThreads
	{
		if( (unsigned int) -maxThreads >= ncores )
			this->maxThreads = 1;
		else
			this->maxThreads = ncores + maxThreads; //remember: maxThreads is negative
	}
	else
	{
		#if DCT_DEBUG_MODE
			assert( maxThreads >= 0 && ncores > 0);
		#endif
		
		this->maxThreads = DCT_min<unsigned int>(maxThreads, ncores);
		
		if(this->maxThreads == 0)
			this->maxThreads = ncores;
	}
	
	
	
	this->maxThreadsByService = DCT_min(maxThreadsByService, this->maxThreads);
	
	if(this->maxThreadsByService <= 0)
		this->maxThreadsByService = this->maxThreads;
	
	
	this->maxServices = maxSimultaneuousServices;
}


int DCT_BBServer::reserveThreads( unsigned int maxThreadsInThisService, unsigned int &nreservedThreads)
{
	unsigned int avthreads; //number of threads available
	
	nreservedThreads = 0;
	
	if( nservices >= maxServices )
		return DCT_RSC_MAX_SIM_SERVICES_REACHED;
	
	
	if( maxThreadsInThisService <= 0 || maxThreadsInThisService > this->maxThreadsByService)
		maxThreadsInThisService = this->maxThreadsByService;
	
	
	SEMAPH_nthreads.lock();
	{
		avthreads = maxThreads - nthreads;
		
		if( avthreads > 0 && nservices < maxServices )
		{
			nreservedThreads = DCT_min(maxThreadsInThisService, avthreads);
			nthreads += nreservedThreads;
			nservices++;
		}
	}
	SEMAPH_nthreads.unlock();
	
	#if DCT_SERVER_DEBUG_MODE
		assert(nthreads <= maxThreads);
	#endif
	
	if(avthreads <= 0)
		return DCT_RSC_NO_THREADS_AVAILABLE;
		
	
	return nreservedThreads;
}


void DCT_BBServer::freeThreads(unsigned int numberOfThreads)
{
	if(numberOfThreads == 0)
		return;
	
	#if DCT_SERVER_DEBUG_MODE
		assert(nthreads >= numberOfThreads);
	#endif
	
	SEMAPH_nthreads.lock(); //we use 2 here to always guarantee we will lock it correctly (or maybe do not lock)
		nthreads -= numberOfThreads;
		nservices--;
	SEMAPH_nthreads.unlock();
}






//function to be runned in a separated thread for each connection in a muriqui server
int DCT_connectionThreadHandler( DCT_BBServer *server, dctools::DCT_ServerConnection *connection, const long unsigned int serviceNumber, DCT_ServiceCoreGenerator *serviceCoreGenerator, std::unordered_map<long unsigned int, DCT_RunningThread> *openThreads)
{
	const int sizeBuffer = 1024;
	const int maxResponseSize = 5;
	
	int ioRetCode; //I/O return code
	unsigned int mynthreads = 0;
	size_t bytesTransfered;
	DCT_Int32 requestCode;
	DCT_UInt32 nRequestedNodes = 0;
	DCT_UInt64 sizeOfOpenNodesRep = 0;
	DCT_Byte buffer[sizeBuffer];
	
	DCT_Int32 response[maxResponseSize]; 
	//char *pbuffer;
	
	double lb = -INFINITY, ub = INFINITY;
	
	DCT_UInt64 nBasicInputParameters = 0;
	DCT_Byte *basicInputParameters = NULL;
	
	DCT_Byte *openNodesRep = NULL;
	
	//DCT_VarBounds varsBounds;
	DCT_FileNames inputFiles;
	DCT_AllGeneralParams allGeneralParams;
	DCT_ServiceCore *serviceCore = NULL;
	
	
	//DCT_BranchAndBound bb;
	//DCT_Algorithm *alg = &bb; //maybe, in the future, we allow any algorithm be used. By now, we only accept branch-and-bound
	
	
	//std::cerr << "oi 1" << std::endl;
	
	
	{
		int r = serviceCoreGenerator->generateServiceCore(serviceCore);
		if(r != 0)
		{
			std::cerr << DCT_PREPRINT "Error: callback user function to generate service core returned " << r << DCT_GETFILELINE;
			
			goto endloop;
		}
	}
	
	//we have a connection with a client. We must receive client request and act acoording them
	
	while(true)
	{
		//ioRetCode = connection->readData(buffer, sizeBuffer *sizeof(buffer[0]));
		
		//std::cerr << "oi 2\n" ;
		
		//first, we only read respecting to a DCT_ConectionCode
		ioRetCode = connection->readData(&requestCode, sizeof(requestCode), NULL);
		
		if(ioRetCode != 0)
		{
			if(ioRetCode == DCT_RC_CLOSED_CONNECTION)
			{
				//std::cout << "connection " << serviceNumber << " was closed by client\n";
				break;
			}
			else
			{
				//I think we never reach this point, because here we must be performing blocking reads
				break;
				//continue;
			}
		}
		
		{
			std::string strCode;
			
			DCT_enumToStr( (DCT_CONECTION_CODE) requestCode, strCode );
			
			std::cout << DCT_PREPRINT "service: " << serviceNumber << " requestCode: " << requestCode << " (" << strCode <<  ")\n" ;
		}
		
		//now, we must interpret codes from 
		
		if(requestCode == DCT_CC_SERVICE_REQUEST)
		{
			DCT_UInt32 requestedThreads;
			const int responseSize = 2;
			DCT_Int32 r;
			
			
			/*
			 * here, we read:
			 * 1 - the number of threads (uint32)
			 */
			
			/*
			 * here, we writte:
			 * 1 - the service response (int32)
			 * 2 - the number of threads or service refused code (int32)
			 */ 
			
			
			ioRetCode = connection->readData(&requestedThreads, sizeof(requestedThreads), NULL);
			
			if(ioRetCode != 0)
			{
				//connection was closed or we have other problem
				break;
			}
			
			
			if(mynthreads > 0) //mynthreads is the number of threads effectivelly reserved to this client. So it must be zero here, because we still go reserve the threads
			{
				//note, that is only possible if client already sent a DCT_SERVICE_REQUEST, anyway, we put it even so
				server->freeThreads(mynthreads);
				
				#if DCT_SERVER_DEBUG_MODE
					DCT_PRINTERRORMSG("mynthreads was already positive. Probably client send request service two or more times in the same service.\n");
				#endif
			}
			
			
			r = server->reserveThreads(requestedThreads, mynthreads);
			
			if( mynthreads > 0 )
			{
				//so, we can accept the service
				response[0] = DCT_CC_SERVICE_ACCEPTED;
				response[1] = mynthreads;
			}
			else
			{
				response[0] = DCT_CC_SERVICE_REFUSED;
				response[1] = r;
			}
			
			ioRetCode = connection->writeData(&response, responseSize *sizeof(response[0]), &bytesTransfered);
			
			if(ioRetCode != 0)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORNUMBER(ioRetCode);
				#endif
				break;
			}
			
			if(bytesTransfered == 0)
			{
				//connection was closed by client
				break;
			}
			
			
			if(mynthreads <= 0)
			{
				break;
			}
			
		}
		else if(requestCode == DCT_CC_VERSION_REQUEST)
		{
			const int responseSize = 4;
			
			/*
			 * here, we read:
			 * 
			 * nothing
			 */ 
			
			/*
			 * here we writte
			 * 
			 * 1 - version response code (int32)
			 * 2 - number of number being responsed (n) (int32), n >= 2
			 * 3 - server major version (int32)
			 * 4 - server minor version (int32)
			 * 5 - other remainder n-2 numbers
			 */ 
			
			
			response[0] = DCT_CC_VERSION_RESPONSE;
			response[1] = 2;
			response[2] = DCT_SERVER_MAJOR_VERSION;
			response[3] = DCT_SERVER_MINOR_VERSION;
			
			ioRetCode = connection->writeData(&response, responseSize *sizeof(response[0]), &bytesTransfered);
			
			if(ioRetCode != 0)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORNUMBER(ioRetCode);
				#endif
				break;
			}
			
			if(bytesTransfered == 0)
			{
				//connection was closed by client
				break;
			}
		}
		else if(requestCode == DCT_CC_SERVICE_STATUS_REQUEST)
		{
			const int responseSize = 2;
			
			/*
			 * here, we read:
			 * 
			 * nothing
			 */ 
			
			/*
			 * here we writte:
			 * 
			 * 1 - DCT_SERVICE_STATUS_RESPONSE connection code
			 * 2 - DCT_SERVICE_STATUS_RESPONSE code
			 *
			 */ 
			
			
			response[0] = DCT_CC_SERVICE_STATUS_RESPONSE;
			response[1] = DCT_SSC_SERVICE_NOT_STARTED;
			
			ioRetCode = connection->writeData(&response, responseSize *sizeof(response[0]), &bytesTransfered);
			
			if(ioRetCode != 0)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORNUMBER(ioRetCode);
				#endif
				break;
			}
			
			
			if(bytesTransfered <= 0)
			{
				//connection was closed by client
				break;
			}
		}
		else if( mynthreads == 0 )
		{
			const int responseSize = 2;
			
			//service was not accept. So, we do not response any aditional request to allocate resources
			std::cerr << DCT_PREPRINT "service was not accepeted and received the code: " << requestCode << "\n";
			std::cerr << DCT_PREPRINT "invalid code received by the server: " << requestCode << ". closing the connection.";
			
			
			response[0] = DCT_CC_SERVER_CLOSING_CONNECTION;
			response[1] = DCT_SCC_SERVICE_WAS_NOT_ACCEPTED;
			
			ioRetCode = connection->writeData(response, responseSize, &bytesTransfered);
			
			if(ioRetCode != 0)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORNUMBER(ioRetCode);
				#endif
				break;
			}
			
			if(bytesTransfered <= 0)
			{
				//connection was closed by client
				break;
			}
			
			
			break;
		}
		else if(requestCode == DCT_CC_START_SERVICE)
		{
			bool stopService = false;
			DCT_SERVER_CLOSING_CODE responseCodeToClientIfStop = DCT_SCC_REQUEST_BY_CLIENT;
			
			/*
			 * here we writte, beffore run:
			 * 
			 * nothing
			 * 
			 * after run, if stopService, is true, we writte:
			 * 
			 * 	1 - DCT_SERVER_CLOSING_CONNECTION (int32)
			 * 	2 - DCT_SERVER_CLOSING_CODES (int32)
			 */ 
			
			#if 0
			const int responseSize = 2;
			/*
			 * here we writte, beffore run:
			 * 
			 * 1 - SERVICE_STATUS connection code
			 * 2 - DCT_SERVICE_RUNNING
			 * 
			 */ 
			
			response[0] = DCT_SERVICE_STATUS;
			response[1] = DCT_SERVICE_RUNNING;
			
			ioRetCode = connection->writeData(response, responseSize*sizeof(response[0]), &bytesTransfered);
			
			if(ioRetCode != 0)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORNUMBER(ioRetCode);
				#endif
				break;
			}
			
			if(bytesTransfered <= 0)
			{
				//connection was closed by client
				break;
			}
			#endif
			
			connection->disableBlock();
			
			
			//int r = serviceCore->run(server, connection, serviceNumber, mynthreads, nBasicInputParameters, basicInputParameters, inputFiles, allGeneralParams, varsBounds, lb, ub, nRequestedNodes, stopService, responseCodeToClientIfStop);
			
			int r = serviceCore->run(server, connection, serviceNumber, mynthreads, nBasicInputParameters, basicInputParameters, inputFiles, allGeneralParams, sizeOfOpenNodesRep, openNodesRep, lb, ub, nRequestedNodes, stopService, responseCodeToClientIfStop);
			
			std::cout << DCT_PREPRINT " service core return code: " << r << "\n";
			
			nRequestedNodes = 0;
			
			if( stopService )
			{
				const int responseSize = 2;
			
				std::cerr << DCT_PREPRINT "Service core request stop the service. closing the connection.";
				
				
				response[0] = DCT_CC_SERVER_CLOSING_CONNECTION;
				response[1] = responseCodeToClientIfStop;
				
				ioRetCode = connection->writeData(response, responseSize, &bytesTransfered);
				
				if(ioRetCode != 0)
				{
					#if DCT_DEBUG_MODE
						DCT_PRINTERRORNUMBER(ioRetCode);
					#endif
				}
				
				break;
			}
			
			connection->enableBlock();
			
			
			DCT_secFree(openNodesRep);
			sizeOfOpenNodesRep = 0;
			
			
			lb = -INFINITY;
			ub = INFINITY;
			
			/*ioRetCode = connection->writeData(response, responseSize, bytesTransfered);
			
			if(ioRetCode != 0)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORNUMBER(ioRetCode);
				#endif
				break;
			}*/
			
		}
		else if(requestCode == DCT_CC_STOP_SERVICE)
		{
			/* here, we read: nothing
			 *
			 * here, we write:
			 * 
			 * 	1 - DCT_SERVER_CLOSING_CONNECTION (int32)
			 * 	2 - DCT_SSC_REQUEST_BY_CLIENT (int32)
			 */
			
			const int responseSize = 2;
			
			std::cerr << DCT_PREPRINT "Client request stop the service. closing the connection.";
			
			
			response[0] = DCT_CC_SERVER_CLOSING_CONNECTION;
			response[1] = DCT_SCC_REQUEST_BY_CLIENT;
			
			ioRetCode = connection->writeData(response, responseSize, &bytesTransfered);
			
			if(ioRetCode != 0)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORNUMBER(ioRetCode);
				#endif
			}
			
			break;
		}
		else if(requestCode == DCT_CC_BASIC_INPUT_PARAMETERS)
		{
			nBasicInputParameters = 0;
			DCT_secFree(basicInputParameters);
			
			/* here, we read:
			 * 1 - the number of bytes of basic input parameters (n) (uint64)
			 * 2 - and so, we read n bytes parameters */
			
			/*
			 * here, we write:
			 * 
			 * nothing
			 */ 
			
			ioRetCode = connection->readData(&nBasicInputParameters, sizeof(nBasicInputParameters), NULL);
			
			if(ioRetCode != 0)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORNUMBER(ioRetCode);
				#endif
				break;
			}
			
			
			if( nBasicInputParameters > 0 )
			{
				DCT_malloc(basicInputParameters, nBasicInputParameters);
				
				if( !basicInputParameters )
				{
					DCT_PRINTMEMERROR;
					break;
				}
				
				ioRetCode = connection->readData(basicInputParameters, nBasicInputParameters * sizeof(basicInputParameters[0]), NULL );
				
				if(ioRetCode != 0)
				{
					#if DCT_DEBUG_MODE
						DCT_PRINTERRORNUMBER(ioRetCode);
					#endif
					break;
				}
				
			}
			
			//here, we do not response anything
		}
		else if(requestCode == DCT_CC_INPUT_PROBLEM_FILE)
		{
			//Warning: here, the first parameter is an unsigned 64 bit integer!
			DCT_Int32 fileType;
			uint64_t fileSize = 0, totalBytesRead = 0;
			char inputFileName[L_tmpnam + 1];
			FILE *inputFile = NULL;
			
			
			/*
			 * here, we read:
			 *
			 * 1 - a 32 integer to specify the type o file (we use that as key in our maps)
			 * 2 - the file size (64 bits)
			 * 3 - the file.
			 */
			
			/*
			 * here, we write:
			 * 
			 * nothing
			 */
			
			
			
			
			ioRetCode = connection->readData( &fileType, sizeof(fileType), NULL );
			if(ioRetCode != 0)
			{
				DCT_PRINTERRORNUMBER(ioRetCode);
				break;
			}
			
			ioRetCode = connection->readData( &fileSize, sizeof(fileSize), NULL );
			if(ioRetCode != 0)
			{
				DCT_PRINTERRORNUMBER(ioRetCode);
				break;
			}
			
			tmpnam(inputFileName);
			
			
			//inputFile = tmpfile();
			inputFile = fopen(inputFileName, "wb");
			
			if(!inputFile)
			{
				DCT_PRINTERRORMSG("Failure to create temporary file.\n");
				break;
			}
			
			auto item = inputFiles.find(fileType);
			if(item != inputFiles.end())
			{
				remove(item->second.c_str());
			}
			
			try
			{
				inputFiles[fileType] = inputFileName; //inputFile;
			}
			catch(std::bad_alloc &ba)
			{
				DCT_PRINTMEMERROR;
				break;
			}
			
			
			while(totalBytesRead < fileSize)
			{
				size_t bytesToRead = DCT_min<decltype(fileSize)>( sizeBuffer, fileSize - totalBytesRead);
				size_t bytesRead;
				
				ioRetCode = connection->readData(buffer, bytesToRead, &bytesRead );
				
				if(ioRetCode != 0)
				{
					DCT_PRINTERRORNUMBER(ioRetCode);
					goto endloop;
				}
				
				totalBytesRead += bytesRead;
				
				size_t bytesWritten = fwrite(buffer, sizeof(buffer[0]), bytesRead, inputFile);
				
				if(bytesWritten != bytesRead*sizeof(buffer[0]) )
				{
					DCT_PRINTERRORMSG("Error to write file");
					goto endloop;
				}
			}
			
			#if DCT_DEBUG_MODE
				if( totalBytesRead != fileSize )
				{
					std::cerr << DCT_PREPRINT "inconsistency about file size detected. totalBytesRead: " << totalBytesRead << " fileSize: " << fileSize << DCT_GETFILELINE << "\n";
					
					break;
				}
			#endif
			
			fclose(inputFile);
		}
		else if(requestCode == DCT_CC_SOLVER_PARAMETERS)
		{
			DCT_Int32 key;
			DCT_GeneralParams *params;//new (std::nothrow) DCT_GeneralParams;
			
			/*here we read:
			 * 
			 * 1 - An DCT_Int32 specifying a key for the parameter group (client can use to specify the owner of its parameter. Can be solvers parameters or subsolvers parameters)
			 * 
			 * 
			 * 2 - An DCT_UInt32 specifying number of double paramenter (nd)
			 * 
			 * 3 - nd pairs of <size string (DCT_UInt32), string (char*)> and double (64 bits)
			 * 
			 * 
			 * 4 - An DCT_UInt32 specifying number of integer paramenter (ni)
			 *
			 * 5 - ni pairs of <size string (DCT_UInt32), string (char*)> and a int64_t (NOTE, 64 bits)
			 * 
			 * 
			 * 6 - An DCT_UInt32 specifying number of string paramenter (ns)
			 *
			 * 7 - ns pairs of <size string (DCT_UInt32), string (char*)> and <size string (DCT_UInt32), string (char*)>
			 * 
			 * 
			 * ATTENTION: size of strings here include "\0" characteres
			 */ 
			
			/*
			 * here, we write:
			 * 
			 * nothing
			 */ 
			
			
			
			ioRetCode = connection->readData(&key, sizeof(key), NULL);
			if(ioRetCode != 0)
			{
				DCT_PRINTERRORNUMBER(ioRetCode);
				break;
			}
			
			
			auto item = allGeneralParams.find(key);
			if(item != allGeneralParams.end() )
			{//so, there is already an element in this key. We delete the current element to put a new
				delete item->second;
				//allGeneralParams.erase(key);
			}
			
			params = new (std::nothrow) DCT_GeneralParams;
			if(!params)
			{
				DCT_PRINTMEMERROR;
				break;
			}
			
			
			try
			{
				allGeneralParams[key] = params; //we put it here to guaratee params will be 
			}
			catch( std::bad_alloc &ba )
			{
				delete params;
				DCT_PRINTMEMERROR;
				break;
			}
			
			
			{
				DCT_UInt32 nd;
				
				ioRetCode = connection->readData(&nd, sizeof(nd), NULL);
				if(ioRetCode != 0)
				{
					DCT_PRINTERRORNUMBER(ioRetCode);
					break;
				}
				
				for(DCT_UInt32 i = 0; i < nd; i++)
				{
					int r = DCT_readParameterToMap(*connection, params->dblParams, (char*) buffer, sizeBuffer);
					
					if(r != 0)
					{
						DCT_PRINTERRORNUMBER(r);
						goto endloop;
					}
				}
			}
			
			
			{
				DCT_UInt32 ni;
				
				ioRetCode = connection->readData(&ni, sizeof(ni), NULL);
				if(ioRetCode != 0)
				{
					DCT_PRINTERRORNUMBER(ioRetCode);
					break;
				}
				
				for(DCT_UInt32 i = 0; i < ni; i++)
				{
					int r = DCT_readParameterToMap(*connection, params->intParams, (char*) buffer, sizeBuffer);
					
					if(r != 0)
					{
						DCT_PRINTERRORNUMBER(r);
						goto endloop;
					}
				}
			}
			
			
			{
				DCT_UInt32 ns;
				
				ioRetCode = connection->readData(&ns, sizeof(ns), NULL);
				if(ioRetCode != 0)
				{
					DCT_PRINTERRORNUMBER(ioRetCode);
					break;
				}
				
				
				for(DCT_UInt32 i = 0; i < ns; i++)
				{
					int r = DCT_readParameterToMap(*connection, params->strParams, (char*) buffer, sizeBuffer);
					
					if(r != 0)
					{
						DCT_PRINTERRORNUMBER(r);
						goto endloop;
					}
				}
			}
			
			
		}
		#if 0
		else if(requestCode == DCT_CC_VAR_BOUNDS)
		{
			/* Here, we read:
			 * 
			 * 1 - The number of variable bounds to receive (n) (DCT_UInt32)
			 * 
			 * 2 - Read n times: 
			 * 		2.1 - An index (DCT_UInt32)
			 * 		2.2 - A lower bound (double)
			 * 		2.3 - A upper bound (double)
			 * 
			 * Here, we write:
			 * 		nothing
			 */
			
			assert(false);
			
			int r = DCT_readVarBounds(*connection, varsBounds); 
			
			if(r != 0)
			{
				DCT_PRINTERRORNUMBER(r);
				break;
			}
			
		}
		#endif
		else if(requestCode == DCT_CC_OPEN_NODE_RESPONSE)
		{
			/* Here, we read:
			 * 
			 * 1 - The number m of bytes representing the open nodes
			 * 2 - m bytes representing open nodes. We just pass this m bytes to user server and let it interpret this bytes
			 * 
			 * Here, we write:
			 * 		nothing
			 */
			
			auto oldSizeOfOpenNodesRep = sizeOfOpenNodesRep;
			
			ioRetCode = connection->readData(&sizeOfOpenNodesRep, sizeof(sizeOfOpenNodesRep));
			if(ioRetCode != 0)
			{
				DCT_PRINTERRORNUMBER(ioRetCode);
				break;
			}
			
			if( oldSizeOfOpenNodesRep != sizeOfOpenNodesRep )
			{
				DCT_realloc(openNodesRep, sizeOfOpenNodesRep);
				if(!openNodesRep)
				{
					DCT_PRINTMEMERROR;
					break;
				}
			}
			
			ioRetCode = connection->readData(openNodesRep, sizeOfOpenNodesRep);
			if(ioRetCode != 0)
			{
				DCT_PRINTERRORNUMBER(ioRetCode);
				break;
			}
			
			
			
			#if 0
			int r;
			ioRetCode = connection->readData(&lb, sizeof(lb));
			if(ioRetCode != 0)
			{
				DCT_PRINTERRORNUMBER(ioRetCode);
				break;
			}
			
			r = DCT_readVarBounds(*connection, varsBounds); 
			if(r != 0)
			{
				DCT_PRINTERRORNUMBER(r);
				break;
			}
			
			#endif
			
		}
		else if(requestCode == DCT_CC_NUMBER_OF_INITIAL_OPEN_NODES_REQUEST)
		{
			/*
			 * Here, we read:
			 * 		1 - The number of initial open nodes that should be separed to send to other servers (DCT_UInt32)
			 * 
			 * Here, we write:
			 * 		nothing
			 */
			
			ioRetCode = connection->readData( &nRequestedNodes, sizeof(nRequestedNodes));
			if(ioRetCode != 0)
			{
				DCT_PRINTERRORNUMBER(ioRetCode);
				break;
			}
			
		}
		else if(requestCode == DCT_CC_CLEAR_VAR_BOUNDS)
		{
			/*
			 * Here, we read:
			 * 		nothing
			 * 
			 * Here, we write:
			 * 		nothing
			 */
			
			DCT_secFree(openNodesRep);
			sizeOfOpenNodesRep = 0;
			
			//varsBounds.clear();
		}
		else if(requestCode == DCT_CC_OPEN_NODE_REQUEST)
		{
			/*
			 * Here, we read:
			 * 		nothing
			 * 
			 * Here, we write:
			 * 		1 - A response code DCT_NO_OPEN_NODE_SERVICE_NOT_STARTED (DCT_Int32)
			 */
			
			response[0] = DCT_CC_NO_OPEN_NODE_SERVICE_NOT_STARTED;
			
			ioRetCode = connection->writeData(&response, sizeof(response[0]), &bytesTransfered);
			
			if(ioRetCode != 0)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORNUMBER(ioRetCode);
				#endif
				break;
			}
			
		}
		else if(requestCode == DCT_CC_NEW_SOLUTION_FOUND)
		{
			/*
			 * We just put it here because server can receive a feasible solution before client knows it is finished. We just read the solution and do nothing
			 * 
			 * Here we read:
			 * 		1 - The number n of variables (uint32)
			 * 		2 - The objective value (double)
			 * 		3 - The values of n variables 
			 * 
			 */
			
			DCT_UInt32 n;
			double v;
			
			ioRetCode = connection->readData( &n, sizeof(n) );
			if(ioRetCode != 0)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORNUMBER(ioRetCode);
				#endif
				break;
			}
			
			for(unsigned int i = 0; i <= n; i++ ) //we need read n+1 numbers
			{
				ioRetCode = connection->readData( &v, sizeof(v));
				if(ioRetCode != 0)
				{
					#if DCT_DEBUG_MODE
						DCT_PRINTERRORNUMBER(ioRetCode);
					#endif
					goto endloop;
				}
			}
			
		}
		else if(requestCode == DCT_CC_LOWER_BOUND_RESPONSE)
		{
			/* Here, we read:
			 * 		1 - A lower bound for current partition (double)
			 * Here, we write: nothing
			 */ 
			
			ioRetCode = connection->readData( &lb, sizeof(double) );
			if(ioRetCode != 0)
			{
				DCT_PRINTERRORNUMBER(ioRetCode);
				break;
			}
		}
		else if(requestCode == DCT_CC_UPPER_BOUND_RESPONSE)
		{
			/* Here, we read:
			 * 		1 - An upper bound for current partition (double)
			 * Here, we write: nothing
			 */
			
			ioRetCode = connection->readData( &ub, sizeof(double) );
			if(ioRetCode != 0)
			{
				DCT_PRINTERRORNUMBER(ioRetCode);
				break;
			}
		}
		else if(requestCode == DCT_CC_NO_SOLUTION_AVAILABLE)
		{
			/* Here, we read: nothing
			 * Here, we write: nothing
			 * 
			 * we just put this if to avoid server abort because it does not know the code (ok, i just need put an else if, but i will let things in this way)
			 */ 
		}
		else
		{
			const int responseSize = 2;
			
			std::cerr << DCT_PREPRINT "invalid code received by the server: " << requestCode << ". closing the connection.";
			
			
			response[0] = DCT_CC_SERVER_CLOSING_CONNECTION;
			response[1] = DCT_SCC_INVALID_CONNECTION_CODE;
			
			ioRetCode = connection->writeData(response, responseSize, &bytesTransfered);
			
			if(ioRetCode != 0)
			{
				#if DCT_DEBUG_MODE
					DCT_PRINTERRORNUMBER(ioRetCode);
				#endif
				break;
			}
			
			/*if(bytesTransfered <= 0)
			{
				//connection was closed by client
				break;
			} */
			
			
			break;
		}
		
	}
	
	
endloop:
	
	//std::cout << "tchau 1\n" ;
	
	server->freeThreads(mynthreads);
	
	
	/*if(ioRetCode == 0)
	{
		std::cout << DCT_PREPRINT "connection " << serviceNumber << " was closed by client\n";
	} */
	
	
	
	std::cout << DCT_PREPRINT "closing connection " << serviceNumber << "\n";
	
	connection->closeSocket();
	
	
	
	
	
	for(auto pair : inputFiles)
	{
		const char *fileName = pair.second.c_str();
		remove(fileName);
		//if(file)
			//fclose(file);
	}
	
	for(auto pair : allGeneralParams )
	{
		delete pair.second;
	}
	
	
	if(basicInputParameters)  free(basicInputParameters);
	
	if(serviceCore)	delete serviceCore;
	
	if(openNodesRep)	free(openNodesRep);

	//deleting thread object to be deleted, so, we can free the memory
	openThreads->at(serviceNumber).running = false;
	
	return 0;
}










int DCT_BBServer::run(int portNumber)
{
	const int MAX_CLIENT_NAME = 200;
	
	int retCode, r;
	long unsigned int serviceNumber = 0;
	char clientName[MAX_CLIENT_NAME];
	
	DCT_ServerSocket socket;
	DCT_RunningThread runningThread;
	DCT_ServerConnection *connection = NULL;
	
	std::thread *mythread = NULL;
	
	std::unordered_map<long unsigned int, DCT_RunningThread> openThreads;
	
	
	
	
	r = socket.openSocket(portNumber);
	if(r != 0)
	{
		DCT_PRINTERRORMSGP("Error to open socket ", r);
		retCode = DCT_SRC_SOCKET_ERROR;
		
		goto termination;
	}
	
	
	std::cout << DCT_PREPRINT "Server is working on port " << portNumber << " with socket file descriptor " << socket.getFileDescriptor() << "\n";
	
	
	while(true)
	{
		connection = new (std::nothrow) DCT_ServerConnection;
		DCT_IFMEMERRORGOTOLABEL(!connection, retCode, termination );
		
		
		DCT_PRINTMSG("waiting for a conection\n");
		
		r = socket.acceptConnection(clientName, MAX_CLIENT_NAME, *connection);
		
		if( r != 0 )
		{
			DCT_PRINTERRORNUMBER(r);
			std::cout << DCT_PREPRINT "A connection gott failure to be built. Error code: " << r << "\n";
			
			continue;
		}
		
		
		std::cout << DCT_PREPRINT "Connection with " << clientName << " with file descriptor: " << connection->getFileDescriptor() << " service number: " << serviceNumber << "\n";
		
		mythread = new (std::nothrow) std::thread;
		DCT_IFMEMERRORGOTOLABEL(!mythread, retCode, termination);
		
		
		runningThread.running = true;
		runningThread.pthread = mythread;
		
		try
		{
			openThreads.insert( {serviceNumber, runningThread} );
		}
		catch(std::bad_alloc &ba)
		{
			if(mythread)	delete mythread;
			
			DCT_PRINTMEMERROR;
			retCode = DCT_RC_MEMORY_ERROR;
			
			goto termination;
		}
		
		
		
		
		*mythread = std::thread( DCT_connectionThreadHandler, this, connection, serviceNumber, serviceCoreGenerator, &openThreads );
		
		
		#if 1
		//removing ended threads from openThreads
		{
			std::unordered_set<long unsigned int> remove;
			
			
			for(auto &pairThreads : openThreads)
			{
				auto &serviceNumber = pairThreads.first;
				DCT_RunningThread &rthreads = pairThreads.second;
				
				if( rthreads.running == false )
				{
					rthreads.pthread->join();
					
					delete rthreads.pthread;
					remove.insert( serviceNumber );
				}
			}
			
			//we just remove here because I am fair to remove while I am iterating
			for(auto &snumber : remove)
				openThreads.erase(snumber);
		}
		#endif
		
		
		serviceNumber++;
	}
	
	//std::cout << "Sai do laco!\n";
	
	
	retCode = 0;
	
termination:
	
	
	//we must guarantee services still running could not be aborted. So, we run join method, and after, we delete the thread object
	for( auto &pairThreads : openThreads )
	{
		std::thread *p = pairThreads.second.pthread;
		
		p->join();
		
		//after join, the thread is finished
		delete p;
	}
	
	
	if(connection)	delete connection;
	
	return retCode;
}


