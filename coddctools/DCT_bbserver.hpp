

#ifndef DCT_BBSERVER_HPP
#define DCT_BBSERVER_HPP

#include <string>

#include "DCT_dctools.hpp"
#include "DCT_sockets.hpp"



namespace dctools
{
	
	#define DCT_SERVER_DEBUG_MODE 1
	
	//version of server is defined by MAJOR_VERSION.MINOR_VERSION ex: 0.3, 1.6, 2.1 etc 
	
	#define DCT_SERVER_MAJOR_VERSION DCT_MAJOR_VERSION
	#define DCT_SERVER_MINOR_VERSION DCT_MINOR_VERSION
	
	
	
	
	
	
	
	
	
	
	class DCT_BBServer;
	
	
	
	/*
	 * user should derive this class to implement the server
	 */
	class DCT_ServiceCore
	{
	protected:
		
		int sendResume(DCT_ServerConnection *connection, DCT_FinalResults &finalResults, DCT_UInt32 nvars, double *solution);
		
		
	public:
		
		virtual ~DCT_ServiceCore(){}
		
		//if stopService will set as true, DCT_BBServer will close the connection with the client sending DCT_SERVER_CLOSING_CONNECTION and send responseCodeToClientIfStop
		
		virtual int run(DCT_BBServer *server, DCT_ServerConnection *connection, const long unsigned int serviceNumber, const unsigned int nthreads, DCT_Int64 nBasicInputParameters, DCT_Byte *basicInputParameters, DCT_FileNames &inputFiles, DCT_AllGeneralParams &allGeneralParams, DCT_Int64 sizeOfOpenNodeRep, DCT_Byte *openNodeRep, double lowerBound, double upperBound, DCT_UInt32 nRequestedNodes, bool &stopService, DCT_SERVER_CLOSING_CODE &responseCodeToClientIfStop ) = 0; 
		
		
		//virtual int run(DCT_BBServer *server, DCT_ServerConnection *connection, const long unsigned int serviceNumber, const unsigned int nthreads, int32_t nBasicInputParameters, int32_t *basicInputParameters, DCT_FileNames &inputFiles, DCT_AllGeneralParams &allGeneralParams, DCT_VarBounds &varsBounds, double lowerBound, double upperBound, DCT_UInt32 nRequestedNodes, bool &stopService, DCT_SERVER_CLOSING_CODE &responseCodeToClientIfStop ) = 0;
		
	};
	
	
	//this class should generate an object ServiceCore which its method run will be called to threat the service. For each service, a new object must be generated. At the end of service, object will be deleted, i.e., dctools takes ownership about the object generated
	class DCT_ServiceCoreGenerator
	{
	public:
		
		virtual ~DCT_ServiceCoreGenerator(){}
		
		virtual int generateServiceCore(DCT_ServiceCore* &serviceCore) = 0;
	};
	
	
	
	class DCT_BBServer
	{
		unsigned int maxThreads;
		unsigned int nthreads;
		
		unsigned int maxThreadsByService;
		
		unsigned int maxServices; //maximum of simultaneous services
		unsigned int nservices;
		
		
		DCT_Mutex SEMAPH_nthreads;
		
	public:
		
		DCT_ServiceCoreGenerator *serviceCoreGenerator;
		
		
		
		DCT_BBServer(DCT_ServiceCoreGenerator &serviceCoreGenerator);
		
		~DCT_BBServer();
		
		void setParams(int maxThreads, unsigned int maxThreadsByService, unsigned int maxSimultaneuousServices);
		
		
		void freeThreads(unsigned int numberOfThreads);
		
		//user request a max number of threads, but server is not obligated to allocte all threads requested. Number of reserved thread is returned
		int reserveThreads( unsigned int maxThreads, unsigned int &nreservedThreads);
		
		int run(int portNumber);
		
	};
	
	
	
	
	//function to read a value from a array of bytes. ATTENTION: Value is read and pointer is shift to next data.
	template <class myClass>
	void DCT_readValueFromBuffer(DCT_Byte* &pbuffer, myClass &value )
	{
		value = *((myClass*) pbuffer) ;
		
		pbuffer += sizeof(value);
	}
	
	
	
	
	
	
	
	
	
}





#endif
