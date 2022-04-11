

#ifndef DCT_DCTOOLS_HPP
#define DCT_DCTOOLS_HPP

#include <cstdio>
#include <cstdint>

#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <iostream>

#include "DCT_config.hpp"





#if DCT_CPP_MULTITHREADING
    #include <mutex>
#endif




namespace dctools
{
	#define DCT_ATTSTR(att)   att, #att
	
	#define DCT_DEBUG_MODE 1
	
	
	//#define DCT_DEFAULT_PORT 22022
	#define DCT_DEFAULT_MAXTHREADS 0
	#define DCT_DEFAULT_NUMBER_OF_COMPONENTS_TO_BE_USED 1
	
	#define DCT_XML_COMPONENT_NAME "component"
	#define DCT_XML_NUMBEROF_NAME "numberof"
	#define DCT_XML_IPPADDRESS_NAME "ipaddress"
	#define DCT_XML_PORT_NAME "port"
	#define DCT_XML_MAXTHREADS_NAME "maxThreads"
	#define DCT_XML_NUMBERINNUMBEROF_NAME "numberof"
	
	
	#define DCT_MAJOR_VERSION 0
	#define DCT_MINOR_VERSION 1
	
	#define DCT_SAVE_CLIENT_SERVER_LOGS 1
	#define DCT_CLIENT_SERVER_LOGS_FILE_NAME_PREFIX "dct_output_bbclient_server_"
	
	#define DCT_RESEND_DATA_IF_EGAIN_IS_GOTTEN 1
	
	#define DCT_SECONDS_TO_SLEEP_WAITING_FOR_BOUNDS 30
	
	
	typedef int8_t 		DCT_Byte;
	typedef int32_t 	DCT_Int32;
	typedef int64_t 	DCT_Int64;
	typedef uint32_t	DCT_UInt32;
	typedef uint64_t 	DCT_UInt64;
	
	#define DCT_Dict std::unordered_map
	
	
	enum DCT_RETURN_CODE
	{
		DCT_RC_SUCCESS					= 0,
		DCT_RC_MEMORY_ERROR				= -1,
		DCT_RC_NAME_ERROR				= -2,
		DCT_RC_VALUE_ERROR				= -3,
		DCT_RC_IO_ERROR					= -4,
		DCT_RC_ADDRESS_ERROR			= -5,
		DCT_RC_NETWORK_ERROR			= -6,
		DCT_RC_BUFFER_ERROR				= -7,
		DCT_RC_SERVER_ERROR				= -8,
		DCT_RC_REFUSED_SEVICE			= -9,
		DCT_RC_NO_SERVER_AVAILABLE 		= -10,
		DCT_RC_CLOSED_CONNECTION		= -11,
		DCT_RC_NO_DATA_TO_READ			= -12,
		DCT_RC_UNDEFINED_ERROR 			= -13,
		DCT_RC_CALLBACK_FUNCTION_ERROR 	= -14,
	};
	
	
	//codes to optimization jobs performed by serves
	enum DCT_OPTIMIZATION_RETURN_CODE
	{
		DCT_ORC_OPTIMAL_SOLUTION		= 4000,
		DCT_ORC_INFEASIBLE_PROBLEM,
		DCT_ORC_UNBOUNDED_PROBLEM,
		DCT_ORC_CLIENT_ERROR,
		DCT_ORC_EVALUATION_ERROR,
		DCT_ORC_LIBRARY_NOT_AVAILABLE,
		DCT_ORC_MAX_ITERATIONS,
		DCT_ORC_MAX_TIME,
		DCT_ORC_MILP_SOLVER_ERROR,
		DCT_ORC_NLP_SOLVER_ERROR,
		DCT_ORC_MEMORY_ERROR,
		DCT_ORC_SERVER_CONNECTION_ERROR,
		DCT_ORC_SUBROUTINE_ERROR,
		DCT_ORC_SUBSOLVER_ERROR, //this a generic code compared to DCT_NLP_SOLVER_ERROR and DCT_MILP_SOLVER_ERROR
		DCT_ORC_UNDEFINED_ERROR
	};
	
	int DCT_enumToStr(const DCT_OPTIMIZATION_RETURN_CODE code, std::string &str);
	
	
	enum DCT_SERVER_RETURN_CODE
	{
		DCT_SRC_INVALID_PARAMETERS	= -1000,
		DCT_SRC_SOCKET_ERROR
	};
	
	
	//codes to be used in the connection between client and server
	enum DCT_CONECTION_CODE
	{
		DCT_CC_SERVICE_REQUEST		= 5000,
		DCT_CC_SERVICE_ACCEPTED,
		DCT_CC_SERVICE_REFUSED,
		
		DCT_CC_VERSION_REQUEST,
		DCT_CC_VERSION_RESPONSE,
		
		DCT_CC_START_SERVICE,
		DCT_CC_STOP_SERVICE,
		
		DCT_CC_SERVICE_STATUS_REQUEST, //running, finished, etc
		DCT_CC_SERVICE_STATUS_RESPONSE,
		
		DCT_CC_BASIC_INPUT_PARAMETERS,    //basic parameters like alg choica, solvers etc. actually, it is possible pass any parameter in a byte array
		DCT_CC_INPUT_PROBLEM_FILE,
		DCT_CC_OUTPUT_PARAMETERS,
		
		DCT_CC_SOLVER_PARAMETERS,  //we use also to send subsolver parameters
		
		DCT_CC_LOWER_BOUND_REQUEST,
		DCT_CC_LOWER_BOUND_RESPONSE,
		
		DCT_CC_UPPER_BOUND_RESPONSE,
		
		DCT_CC_OPEN_NODE_REQUEST,
		DCT_CC_OPEN_NODE_RESPONSE,
		DCT_CC_NO_OPEN_NODE_AVAILABLE,
		DCT_CC_NO_OPEN_NODE_SERVICE_NOT_STARTED, //ok, we could use DCT_NO_OPEN_NODE_AVAILABLE in all situations where we use DCT_NO_OPEN_NODE_SERVICE_NOT_STARTED, but, have a different code for the case where the service had no started gave me luck because helped-me to discover and fix bug!
		
		DCT_CC_VAR_BOUNDS, //deprecated connection code. Use DCT_OPEN_NODE_RESPONSE
		DCT_CC_CLEAR_VAR_BOUNDS,
		
		DCT_CC_BEST_SOLUTION_REQUEST, //if server request the best solution, we answer with the DCT_NEW_SOLUTION_FOUND code
		DCT_CC_NEW_SOLUTION_FOUND, //to send data to client, lower bound is passed also
		DCT_CC_NO_SOLUTION_AVAILABLE,
		
		DCT_CC_FINAL_RESULTS,   //use to send final results of optimization
		
		//DCT_SERVICE_NOT_STARTED,
		
		DCT_CC_SERVER_CLOSING_CONNECTION, //in this code, we pass after an DCT_SERVERCLOSINGCODES to explain why we are closing the connection
		
		DCT_CC_NUMBER_OF_INITIAL_OPEN_NODES_REQUEST, //to inform to the first server the number of nodes that should be separed in the beginning of exploration to be sent to other servers.
		
		DCT_CC_INITIAL_OPEN_NODE_RESPONSE, //to send an initial open node. An initial open node is node requested by means of DCT_CC_NUMBER_OF_INITIAL_OPEN_NODES_REQUEST. Note, in this code, we just send one open node by time. 
	};
	
	int DCT_enumToStr(const DCT_CONECTION_CODE code, std::string &str);
	
	enum DCT_REFUSED_SERVICE_CODE
	{
		DCT_RSC_NO_THREADS_AVAILABLE = 6000,
		DCT_RSC_MAX_SIM_SERVICES_REACHED
	};
	
	enum DCT_SERVICE_STATUS_CODE
	{
		DCT_SSC_SERVICE_NOT_STARTED = 7000,
		DCT_SSC_SERVICE_RUNNING,
		DCT_SSC_CORE_SERVICE_FINISHED
	};
	
	enum DCT_SERVER_CLOSING_CODE
	{
		DCT_SCC_INVALID_CONNECTION_CODE = 8000,
		DCT_SCC_SERVICE_WAS_NOT_ACCEPTED,
		DCT_SCC_REQUEST_BY_CLIENT,
		DCT_SCC_SERVER_ERROR,
	};
	
	
	
	
	//a component is a server to run threads to help us solve MINLP problems.
	class DCT_Component
	{
		
	public:
		
		unsigned int maxThreads;
		unsigned int port;
		char *ipaddress;
		
		
		
		DCT_Component(const char *ipaddress, const unsigned int maxThreads, const unsigned int port);
		
		~DCT_Component();
		
		void desallocate();
		
		int initialize(const char *ipaddress, const unsigned int maxThreads, const unsigned int port);
		
		void print(std::ostream &out = std::cout);
		
	};
	
	//DCT_Components is a vector of components
	class DCT_Components : public std::vector<DCT_Component*>
	{
	public:
		
		//std::vector<DCT_Component*> comps;
		
		virtual ~DCT_Components();
		
		int addComponent(const char *ipaddress, const unsigned int maxThreads, const unsigned int port );
		
		void desallocate();
		
		void print(std::ostream &out = std::cout) const;
		
	};
	
	
	//DCT_NumberOfComponents is also a vector. The diference, is DCT_NumberOfComponents is a group of components where we just use a number of componentes, not all of them.
	class DCT_NumberOfComponents : public DCT_Components
	{
	public:
		unsigned int number; //number of components to be used
		
		DCT_NumberOfComponents(unsigned int number = DCT_DEFAULT_NUMBER_OF_COMPONENTS_TO_BE_USED);
		
		virtual ~DCT_NumberOfComponents();
		
		
		void initialize(unsigned int number);
	};
	
	
	//here, we have a vector of number of components
	class DCT_SeveralNumberOfComponents : public std::vector<DCT_NumberOfComponents*>
	{
	public:
		//std::vector<DCT_NumberOfComponents*> comps;
		
		
		~DCT_SeveralNumberOfComponents();
		
		//object takes the ownership of the pointer
		void addNumberOfComponents(DCT_NumberOfComponents *nofcomponents);
		
		void desallocate();
		
		//return the maximum number of servers that can be used, i.e., the sum of min(number, i.size) for i in  DCT_NumberOfComponents
		unsigned int getMaxNumberOfServers();
		
		void print(std::ostream &out = std::cout);
	};
	
	
	
	class DCT_ComponentsFileReader
	{
		
	public:
		
		unsigned int defaultPort;
		unsigned int defaultMaxThreads;
		
		
		//Function to read a xml file decribing components to perform our distributed computation
		DCT_ComponentsFileReader(const unsigned int defaultPort, const unsigned int defaultMaxThreads = DCT_DEFAULT_MAXTHREADS);
		
		
		int read(const char* fileName, DCT_Components& components, DCT_SeveralNumberOfComponents& numberOfComponents);
	};
	
	
	
	//mutex to implement mutual exclusion on multithread proceduring. In this way, the rest of code do not need worry about haw library to multithreading is being used...
	class DCT_Mutex
	{
		
	public:
		
		#if DCT_CPP_MULTITHREADING
			std::mutex mymutex;
		#else
			#if DCT_OMP_MULTITHREADING
				omp_lock_t mutex;
			#endif
		#endif
		
		
		DCT_Mutex();
		
		void initialize();
		
		int lock( );
		
		int tryLock( );
		
		int unlock( );
		
		void destroy();
		
		~DCT_Mutex();
	};
	
	
	
	/*
	 * This class is inspired on optsolvers::OPT_GeneralSolverParams. We just copy its definition because we do not want include our external librares here
	 * 
	 */
	class DCT_GeneralParams
	{
	public:
		std::map<std::string, DCT_Int64> intParams;
		std::map<std::string, double> dblParams;
		std::map<std::string, std::string> strParams;
		
		
		int addIntegerParameter(const char *name, const int value);
		
		int addDoubleParameter(const char *name, const double value);
		
		int addStringParameter(const char *name, const char *value);
		
		void desallocate();
		
		void print(std::ostream &out = std::cout);
		
		~DCT_GeneralParams();
	};
	
	
	class DCT_Bounds
	{
	public:
		double lb, ub;
	};
	
	
	
	
	class DCT_FinalResults
	{
	public:
		
		DCT_Byte feasSolFound;	//true if a feasible solution was found
		
		DCT_Int32 algorithm;	//algorithm code
		DCT_Int32 optCode;   	//return code
		DCT_UInt32 nthreads;	//number of threads
		
		DCT_UInt64 nServerCalls; //number of server calls
		DCT_UInt64 nIters;		//number of iterations
		
		double lowerBound;		//lower bound
		double upperBound;		//upper bound
		double objBestSol;		//objective value on best solution found
		double objAtFirstRelax; //objective value at first relaxation
		double cpuTime;			//cpu time
		double clockTime;		//clock time
		
		
		DCT_FinalResults();
		
		void initialize();
	};
	
	
	
	//now, we are prefering using maps instead of unordered_map to save memory
	
	typedef std::map<uint32_t, DCT_Bounds> DCT_VarBounds; //we just use map instead of unordered_map because we need put indices in order in int MRQ_ServerBBCallbacks::BB_generateRootNode( MRQ_NewBBNode* &rootNode ). Someday maybe I will fix it sorting the indices by muy self and replacing this map by an unordered_map...
	typedef std::map<int32_t, DCT_GeneralParams*> DCT_AllGeneralParams;
	typedef std::map<int32_t, std::string> DCT_FileNames;
	
	
	long unsigned int DCT_getFileSize(FILE *file);
	
}



#endif

