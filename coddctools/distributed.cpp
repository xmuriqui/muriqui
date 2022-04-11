
#include <cstdio>
#include <climits>
#include <cmath>
#include <cassert>


#include <new>
#include <iostream>

#include "pugixml.hpp"

#include "DCT_dctools.hpp"
#include "DCT_tools.hpp"




using namespace dctools;








DCT_Component::DCT_Component(const char *ipaddress, const unsigned int maxThreads, const unsigned int port)
{
	this->ipaddress = NULL;
	initialize(ipaddress, maxThreads, port);
}


DCT_Component::~DCT_Component()
{
	desallocate();
}


void DCT_Component::desallocate()
{
	DCT_secFree(ipaddress);
}


int DCT_Component::initialize(const char *ipaddress, const unsigned int maxThreads, const unsigned int port)
{
	unsigned int sizeip = 1; //saving space for \0
	
	desallocate();
	
	this->maxThreads = maxThreads;
	this->port = port;
	
	if(ipaddress)
	{
		sizeip += DCT_trimSize(ipaddress);
	}
	
	DCT_malloc(this->ipaddress, sizeip);
	if(!ipaddress)
	{
		#if DCT_DEBUG_MODE
			DCT_PRINTMEMERROR;
		#endif
		
		return DCT_RC_MEMORY_ERROR;
	}
	
	
	DCT_copyTrimString(ipaddress, this->ipaddress);
	
	
	
	return 0;
}


void DCT_Component::print(std::ostream &out)
{
	out << "ipaddress: " << ipaddress << " port: " << port << " maxThreads: " << maxThreads;
}


DCT_Components::~DCT_Components()
{
	desallocate();
}


int DCT_Components::addComponent(const char *ipaddress, const unsigned int maxThreads, const unsigned int port)
{
	DCT_Component *component = new (std::nothrow) DCT_Component(ipaddress, maxThreads, port);
	
	if(!component)
	{
		#if DCT_DEBUG_MODE
			DCT_PRINTMEMERROR;
		#endif
		
		return DCT_RC_MEMORY_ERROR;
	}
	
	push_back(component);
	
	return 0;
}


void DCT_Components::desallocate()
{
	for(auto component : *this )
		delete component;
	
	clear();
}


void DCT_Components::print(std::ostream &out) const
{
	for(auto c: *this)
	{
		c->print(out);
		out << "\n";
	}
}


DCT_NumberOfComponents::DCT_NumberOfComponents(unsigned int number):DCT_Components()
{
	initialize(number);
}


DCT_NumberOfComponents::~DCT_NumberOfComponents()
{
}


void DCT_NumberOfComponents::initialize(unsigned int number)
{
	this->number = number;
}


DCT_SeveralNumberOfComponents::~DCT_SeveralNumberOfComponents()
{
	desallocate();
}


void DCT_SeveralNumberOfComponents::addNumberOfComponents( DCT_NumberOfComponents* nofcomponents)
{
	push_back(nofcomponents);
}




void DCT_SeveralNumberOfComponents::desallocate()
{
	for(auto ncomps : *this)
		delete ncomps;
	
	clear();
}



unsigned int DCT_SeveralNumberOfComponents:: getMaxNumberOfServers()
{
	unsigned int v = 0u;
	
	for(auto nofc : *this)
	{
		v += DCT_min<unsigned int>(nofc->number, nofc->size() );
	}
	
	return v;
}



void DCT_SeveralNumberOfComponents::print( std::ostream &out)
{
	for(auto nofc : *this)
	{
		out << "Number of components:\n";
		nofc->print(out);
	}
}



/* function to get components: this function shoulb be inside DCT_componentsFileReader, but I want hide pugi::xml_node from users 
 */
static int DCT_getComponents(pugi::xml_node ancestral, DCT_Components &components, const unsigned int defaultPort, const unsigned int defaultMaxThreads)
{
	for( pugi::xml_node node : ancestral.children(DCT_XML_COMPONENT_NAME)  )
	{
		unsigned int port = defaultPort;
		unsigned int maxThreads = defaultMaxThreads;
		
		auto ip = node.child_value(DCT_XML_IPPADDRESS_NAME);
		
		
		sscanf(node.child_value(DCT_XML_PORT_NAME), "%u", &port);
		
		sscanf(node.child_value(DCT_XML_MAXTHREADS_NAME), "%u", &maxThreads);
		
		
		int r = components.addComponent(ip, maxThreads, port);
		
		if(r != 0)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTMEMERROR;
			#endif
			return r;
		}
		
		
		std::cout << "Minhas vars: port: " << port << " maxThreads: " << maxThreads << "\n";
		
		
		std::cout << "name: " << node.name() << " ";
		
		std::cout << "ip value: " << node.child_value(DCT_XML_IPPADDRESS_NAME) << "\n";
		
		std::cout << "port: " << node.child_value(DCT_XML_PORT_NAME) << "\n";
		
		std::cout << "maxThreads: " << node.child_value(DCT_XML_MAXTHREADS_NAME) << "\n\n";
		
	}
	
	return 0;
}




DCT_ComponentsFileReader::DCT_ComponentsFileReader(const unsigned int defaultPort, const unsigned int defaultMaxThreads)
{
	this->defaultPort = defaultPort;
	this->defaultMaxThreads = defaultMaxThreads;
}



/*
 * Function to read a xml file decribing components to perform our distributed computation
 */
int DCT_ComponentsFileReader::read(const char *fileName, DCT_Components &components, DCT_SeveralNumberOfComponents &numberOfComponents)
{
	int r;
	pugi::xml_document doc;
	
	
	if(!doc.load_file(fileName))
	{
		DCT_PRINTERRORMSGP("Error to load xml component file ", fileName);
		return DCT_RC_NAME_ERROR;
	}
	
	
	r = DCT_getComponents(doc, components, defaultPort, defaultMaxThreads);
	if(r != 0)
	{
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORNUMBER(r);
		#endif
		return r;
	}
	
	
	for( pugi::xml_node node : doc.children(DCT_XML_NUMBEROF_NAME) )
	{
		int r;
		unsigned int number = DCT_DEFAULT_NUMBER_OF_COMPONENTS_TO_BE_USED;
		
		sscanf(node.child_value(DCT_XML_NUMBERINNUMBEROF_NAME), "%u", &number);
		
		DCT_NumberOfComponents *nofComps = new (std::nothrow) DCT_NumberOfComponents(number);
		
		if(!nofComps)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTMEMERROR;
			#endif
			return DCT_RC_MEMORY_ERROR;
		}
		
		numberOfComponents.addNumberOfComponents(nofComps);
		
		r = DCT_getComponents(node, *nofComps, defaultPort, defaultMaxThreads);
		if(r != 0)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORNUMBER(r);
			#endif
			return r;
		}
		
	}
	
	
	return 0;
}






DCT_Mutex::DCT_Mutex()
{
	initialize();
}


void DCT_Mutex::initialize()
{
	#if DCT_CPP_MULTITHREADING
		
	#else	
		#if DCT_OMP_MULTITHREADING
			omp_init_lock(&mutex);
		#endif
	#endif
}


int DCT_Mutex::lock( )
{
	
	#if DCT_CPP_MULTITHREADING
		//if( nthreads > 1 )
			mymutex.lock();
	#else
		#if DCT_OMP_MULTITHREADING
			//if( nthreads > 1 )
				omp_set_lock(&mutex);
		#endif
	#endif
	
	return 0;
}


int DCT_Mutex::tryLock(  )
{
	#if DCT_CPP_MULTITHREADING
		//if( nthreads > 1 )
			return mymutex.try_lock() == true ? 0 : DCT_RC_UNDEFINED_ERROR;
	#else
		#if DCT_OMP_MULTITHREADING
			//if( nthreads > 1 )
				return omp_test_lock(&mutex) == 1 ? 0 : DCT_RC_UNDEFINED_ERROR;
		#endif
	#endif
	
	return 0;
}


int DCT_Mutex::unlock( )
{
	#if DCT_CPP_MULTITHREADING
		//if( nthreads > 1 )
			mymutex.unlock();
	#else
		#if DCT_OMP_MULTITHREADING
			//if( nthreads > 1 )
				omp_unset_lock(&mutex);
		#endif
	#endif
	
	return 0;
}


void DCT_Mutex::destroy()
{
	#if DCT_OMP_MULTITHREADING
		omp_destroy_lock(&mutex);
	#endif
}


DCT_Mutex::~DCT_Mutex()
{
	destroy();
}





int DCT_GeneralParams::addIntegerParameter(const char *name, const int value)
{
	try
	{
		intParams[name] = value;
	}
	catch (std::bad_alloc &e)
	{
		return DCT_RC_MEMORY_ERROR;
	}
	
	return 0;
}


int DCT_GeneralParams::addDoubleParameter(const char* name, const double value)
{
	
	try
	{
		dblParams[name] = value;
	}
	catch (std::bad_alloc &e)
	{
		return DCT_RC_MEMORY_ERROR;
	}
	
	return 0;
}


int DCT_GeneralParams::addStringParameter(const char *name, const char *value)
{
	
	try
	{
		strParams[name] = value;
	}
	catch (std::bad_alloc &e)
	{
		return DCT_RC_MEMORY_ERROR;
	}
	
	return 0;
}


void DCT_GeneralParams::desallocate()
{
    intParams.clear();
    dblParams.clear();
    strParams.clear();
}


DCT_GeneralParams::~DCT_GeneralParams()
{
    desallocate();
}



void DCT_GeneralParams::print(std::ostream &out)
{
	for(auto &pair: intParams)
		out << pair.first << ": " << pair.second << "\n";
	
	for(auto &pair: dblParams)
		out << pair.first << ": " << pair.second << "\n";
	
	for(auto &pair: strParams)
		out << pair.first << ": " << pair.second << "\n";
}



DCT_FinalResults::DCT_FinalResults()
{
	initialize();
}



void DCT_FinalResults::initialize()
{
	feasSolFound = false;
	
	optCode = DCT_ORC_UNDEFINED_ERROR;
	algorithm = INT_MIN;
	nthreads = 0;
	nServerCalls = 0;
	nIters = 0;
	
	lowerBound = -INFINITY;
	upperBound = INFINITY;
	objBestSol = INFINITY;
	objAtFirstRelax = NAN;
	
	cpuTime = 0.0;
	clockTime = 0.0;
}



long unsigned int dctools::DCT_getFileSize(FILE *file)
{
	const long unsigned int prev = ftell(file);
	long unsigned int lSize;
	
	// obtain file size:
	
	fseek(file , 0L , SEEK_END);
	lSize = ftell(file);
	fseek(file, prev, SEEK_SET); //go back to where we were//rewind(file);
	
	return lSize;
}



