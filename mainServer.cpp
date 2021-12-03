// main to perform test on distributed computation


#include <cstdio>
#include <cstring>
#include <climits>
#include <csignal>


#include <iostream>
#include <new>
#include <exception>


#include "muriqui.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_advanced.hpp"
#include "DCT_bbserver.hpp"



#include "MRQ_server.hpp"


using namespace dctools;
using namespace branchAndBound;
using namespace muriqui;


/*
* this function is to threat a possible segmentation fault signal. We use it
* for the case where a service generates a segmentation fault. We want avoid
* all the server be aborted. So, we throw a exception to abort only the service
* responsable by the segmentation signal
*/
void MRQ_segmentationFaultHanlder(int signal)
{
	MRQ_PRINTERRORMSG("Segmentation Fault handler was called");
	throw MRQ_SegFaultException();
}


/*
* this function is to threat a possible abort signal. We use it
* for the case where a service generates an abort signal (for example, by using some function assert). We want avoid
* all the server be aborted. So, we throw a exception to abort only the service
* responsable by the abort signal
*/
void MRQ_aborttHanlder(int signal)
{
	MRQ_PRINTERRORMSG("Abort handler was called");
	throw MRQ_AbortException();
}






int main(int argc, char **argv)
{
	const int MAX_PORT_VALUE = 65535;
	int retCode;
	
	unsigned int portNumber = -1;
	unsigned int maxThreads = 0, maxThreadsService = 0;
	unsigned int maxServices = UINT_MAX;
	char *autoIpFile = NULL;
	
	MRQ_ServerServiceCoreGenerator serverCore;
	DCT_BBServer server(serverCore);
	
	
	
	MRQ_helloOnceTime();
	
	
	if( argc < 3 )
	{
		std::cout << "usage:\n\n"
		
		<< argv[0] << " -p <port number> [-t <max threads>] [-m <max threads by service>] [-s <max simultaneous services>] [-a <autorized ip file>]\n\n"
		
		"a negative value to -t <max threads> indicates a minimum number of cores that should not be used by the server. For example, if the hardware has 8 cores and user specifies -t -2, so, the server will use only 6 cores in the maximum, letting so 2 cores free.\n\n"
		
		"number of cores detected in the system: " << DCT_getNumCores() << "\n";
		
		
		retCode = MRQ_VALUE_ERROR;
		goto termination;
	}
	
	
	for(int i = 1; i < argc-1; i++)
	{
		if( strcmp("-p", argv[i]) == 0 )
		{
			int r = sscanf(argv[i+1], "%u", &portNumber);
			
			if(r == 0 || portNumber < 0 || portNumber > MAX_PORT_VALUE )
			{
				MRQ_PRINTERRORMSGP("Invalid port number: ", argv[i+1]);
				
				retCode = MRQ_VALUE_ERROR;
				goto termination;
			}
			
			//we increment i to skip this parameter 
			i++;
		}
		else if( strcmp("-t", argv[i]) == 0 )
		{
			sscanf(argv[i+1], "%u", &maxThreads);
			
			//we increment i to skip tis parameter 
			i++;
		}
		else if( strcmp("-m", argv[i]) == 0 )
		{
			sscanf(argv[i+1], "%u", &maxThreadsService);
			
			//we increment i to skip this parameter 
			i++;
		}
		else if( strcmp("-s", argv[i]) == 0 )
		{
			sscanf(argv[i+1], "%u", &maxServices);
			
			//we increment i to skip this parameter 
			i++;
		}
		else if( strcmp("-a", argv[i]) == 0 )
		{
			//we increment i to skip this parameter MAX_PORT_VALUE
			autoIpFile = argv[i+1];
			
			i++;
		}
		
	}
	
	
	if(portNumber > MAX_PORT_VALUE)
	{
		MRQ_PRINTMSG("No valid port was specified!");
		retCode = MRQ_VALUE_ERROR;
		goto termination;
	}
	
	
	
	server.setParams(maxThreads, maxThreadsService, maxServices);
	
	
	#if MRQ_HANDLE_TERMINATE_SIGNALS_ON_SERVER
	{
		#ifdef SIGSEGV
			//setting our function to handle the segmentation fault signal
			signal(SIGSEGV, MRQ_segmentationFaultHanlder);
		#endif
		
		#ifdef SIGABRT
			//setting our function to handle the abort signal (it can be happen by assert or some other called to abort function)
			signal(SIGABRT, MRQ_aborttHanlder);
		#endif
	}
	#endif
	
	server.run(portNumber);
	
	
	
	
	retCode = 0;
termination:
	
	
	return retCode;
}


