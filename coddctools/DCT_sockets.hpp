
#ifndef DCT_SOCKETS_HPP
#define DCT_SOCKETS_HPP

#include <cstdio>

#include "DCT_dctools.hpp"
#include "DCT_tools.hpp"


#if DCT_HAVE_POSIX
	#include <fcntl.h> 
	#include <errno.h>
	#include <signal.h>
	#include <unistd.h>

	#include <netinet/in.h>
	#include <sys/socket.h>

	/*for addr2name*/
	#include <netdb.h>
	#include <arpa/inet.h>
	/*end of for addr2name*/
#endif


namespace dctools 
{

class DCT_AddressAndNames
{
	
public:
	
	static void addr2name(struct in_addr addr, char *name, int namelen);
	
	static int name2addr(char *name, in_addr_t *addrp);
};



//int DCT_writeData(int fd, const void *buf, size_t size);
//int DCT_readData(int fd, void *buf, size_t size);


int DCT_ignoreSigpipe();



//base class for sockets. Users should not instantiate it directly, only instatiate the derived classes.
class DCT_BaseSocket
{
protected:
	int socketfd; //socket file descriptor
	
public:
	
	
	DCT_BaseSocket();
	
	virtual ~DCT_BaseSocket();
	
	
	void closeSocket();
	
	//turn operations over this file description nonblocking. That is useful if you wanna perform a read but you are not sure if there is data in the pipe and you do not want to be blocked due to it
	void disableBlock();
	
	//turn operations over this file description blocking (that is the default)
	void enableBlock();
	
	int getFileDescriptor();
	
	//return number of bytes read
	int readData(void* buffer, size_t size, size_t* bytesRead = 0);
	
	int writeFile(FILE *file);
	
	int writeFile(const char *fileName);
	
	void setFileDescriptor(int fd);
	
	//return 0 if get success
	int writeData(const void* buf, size_t size, size_t* byteswritten = NULL);
	
};



//conections accepeted by a server
class DCT_ServerConnection : public DCT_BaseSocket
{
	
public:
	
	DCT_ServerConnection();
	
	virtual ~DCT_ServerConnection();
};




class DCT_ServerSocket : public DCT_BaseSocket
{
	
public:
	
	DCT_ServerSocket();
	
	virtual ~DCT_ServerSocket();
	
	
	//Wait for a connection request from a host
	// hostn = a buffer that will hold the name of the remote host
	// hostnsize = size of hostn buffer
	int acceptConnection(char* hostn, int hostnsize, DCT_ServerConnection& connection);
	
	int openSocket(int port);
};





class DCT_ClientSocket : public DCT_BaseSocket
{
	
public:
	
	DCT_ClientSocket();
	
	virtual ~DCT_ClientSocket();
	
	
	int connectToAServer(int port, char *hostn);
	
};



template <class myClass>
inline int DCT_writeData(DCT_BaseSocket &socket, const myClass &data)
{
	return socket.writeData( &data, sizeof(myClass) );
}


//write the size of string, the string and, finally, '\0' character. Size includes '\0'
int DCT_writeStringWithSizeAnd0( dctools::DCT_BaseSocket& socket, const char* data, size_t *bytesWritten = NULL);



inline int DCT_errno2DCT_RETURN_CODE(int error)
{
	int code;
	
	
	switch(error)
	{
	case 0:
		code = DCT_RC_SUCCESS;
		break;
		
	case EADDRINUSE:
	case EADDRNOTAVAIL:
		code = DCT_RC_ADDRESS_ERROR;
		break;
		
	case ENETDOWN:
	case ENETUNREACH:
	case ENETRESET:
	case EISCONN:
	case ENOTCONN:
	case ESHUTDOWN:
		code = DCT_RC_NETWORK_ERROR;
		break;
		
	case ENOBUFS:
		code = DCT_RC_BUFFER_ERROR;
		break;
		
	default:
		code = DCT_RC_UNDEFINED_ERROR;
	}
	
	
	return code;
}


	
template <class myClass>
inline int DCT_readParameterToMap( dctools::DCT_BaseSocket &socket,  std::map<std::string, myClass> &map, char *auxBuffer, const unsigned int sizeAuxBuffer )
{
	int retCode;
	size_t bytesRead;
	int ioRetValue, totalBytesRead = 0;
	DCT_UInt32 sizeString;
	myClass value;
	char *myname = NULL;
	char *pname;
	
	
	ioRetValue = socket.readData(&sizeString, sizeof(sizeString), &bytesRead ); //size of string must include "\0"
	
	totalBytesRead += bytesRead;
	
	if( ioRetValue != 0 )
	{
		DCT_PRINTERRORNUMBER(ioRetValue);
		retCode = DCT_SRC_SOCKET_ERROR;
		goto termination;
	}
	
	
	
	if(sizeString > sizeAuxBuffer)
	{
		DCT_malloc( myname, sizeString );
		if(!myname)
		{
			DCT_PRINTMEMERROR;
			
			retCode = DCT_RC_MEMORY_ERROR;
			goto termination;
		}
		
		pname = myname;
	}
	else
	{
		pname =  auxBuffer;
	}
	
	
	ioRetValue = socket.readData(pname, sizeString *sizeof(pname[0]), &bytesRead );
	totalBytesRead += bytesRead;
	if( ioRetValue != 0 )
	{
		DCT_PRINTERRORNUMBER(ioRetValue);
		retCode = DCT_SRC_SOCKET_ERROR;
		goto termination;
	}
	
	
	ioRetValue = socket.readData(&value, sizeof(value), &bytesRead);
	totalBytesRead += bytesRead;
	if( ioRetValue != 0 )
	{
		DCT_PRINTERRORNUMBER(ioRetValue);
		retCode = DCT_SRC_SOCKET_ERROR;
		goto termination;
	}
	
	try
	{
		map[pname] = value;
	}
	catch (std::bad_alloc &e)
	{
		DCT_PRINTMEMERROR;
		
		retCode = DCT_RC_MEMORY_ERROR;;
		goto termination;
	}
	
	
	retCode = 0;
	
termination:
	
	if(myname)	free(myname);
	
	return retCode;
}




template <>
inline int DCT_readParameterToMap<std::string>( dctools::DCT_BaseSocket &socket,  std::map<std::string, std::string> &map, char *auxBuffer, const unsigned int sizeAuxBuffer )
{
	int retCode;
	int ioRetValue, totalBytesRead = 0;
	DCT_UInt32 sizeString;
	char *myname = NULL;
	char *pname;
	size_t bytesRead;
	const int nstrings = 2;
	std::string names[nstrings];
	
	
	for(int i = 0; i < nstrings; i++)
	{
		ioRetValue = socket.readData(&sizeString, sizeof(sizeString), &bytesRead ); //size of string must include "\0"
		
		totalBytesRead = bytesRead;
		
		if( ioRetValue != 0 )
		{
			DCT_PRINTERRORNUMBER(ioRetValue);
			retCode = DCT_SRC_SOCKET_ERROR;
			goto termination;
		}
		
		
		
		if(sizeString > sizeAuxBuffer)
		{
			DCT_malloc( myname, sizeString );
			if(!myname)
			{
				DCT_PRINTMEMERROR;
				
				retCode = DCT_RC_MEMORY_ERROR;
				goto termination;
			}
			
			pname = myname;
		}
		else
		{
			pname = auxBuffer;
		}
		
		ioRetValue = socket.readData(pname, sizeString *sizeof(pname[0]), &bytesRead );
		
		totalBytesRead += bytesRead;
		
		if( ioRetValue != 0 )
		{
			DCT_PRINTERRORNUMBER(ioRetValue);
			retCode = DCT_SRC_SOCKET_ERROR;
			goto termination;
		}
		
		
		try
		{
			//names[i].copy(pname, sizeString-1);
			names[i] = pname;
		}
		catch (std::bad_alloc &e)
		{
			DCT_PRINTMEMERROR;
			
			retCode = DCT_RC_MEMORY_ERROR;
			goto termination;
		}
		
	}
	
	
	try
	{
		map[names[0]] = names[1];
	}
	catch (std::bad_alloc &e)
	{
		DCT_PRINTMEMERROR;
		
		retCode = DCT_RC_MEMORY_ERROR;
		goto termination;
	}
	
	
	
	retCode = 0;
	
termination:
	
	if(myname)	free(myname);
	
	return retCode;
	
}


#if 0
class DCT_VarBoundsWriter
{
public:
	unsigned int sizeBuffer;
	DCT_Byte *buffer;
	
	DCT_VarBoundsWriter();
	
	~DCT_VarBoundsWriter();
	
	void desallocate();
	
	void initialize();
	
	int reallocateBuffer(unsigned int newSize);
	
	int write(DCT_BaseSocket &socket, const DCT_VarBounds &varsBounds, size_t *bytesWritten = NULL);
	
	//int writeOpenNode(DCT_BaseSocket &socket, const DCT_VarBounds &varsBounds, size_t *bytesWritten = NULL);
};
#endif


#if 0
inline int DCT_writeVarBounds(DCT_BaseSocket &socket, DCT_VarBounds &varsBounds, size_t *bytesWritten = NULL)
{
	/*
	 * Here, we write:
	 *  1 - The DCT_VAR_BOUNDS  code (int32)
	 *  2 - The number n of tuple of bounds tranfesred (uint32)
	 *  3 - n tuple of < index (uint32), lb (double), ub (double)>
	 *
	 */ 
	
	const unsigned int nrequest = 2;
	const unsigned int nbounds = 2; //lower and upper bound
	DCT_Int32 request[nrequest] = {DCT_VAR_BOUNDS, varsBounds.size()};
	
	int code, r;
	unsigned int k = 0;
	size_t bytesTransfered, totalbytesTransfered;
	
	
	r = socket.writeData(request, nrequest* sizeof(request[0]), &bytesTransfered );
	
	totalbytesTransfered = bytesTransfered;
	
	if(r != 0)
	{
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORNUMBER(r);
		#endif
		code = r;
		goto termination;
	}
	
	
	for(auto &pairVarBounds : varsBounds)
	{
		const DCT_UInt32 ind = pairVarBounds.first;
		const double bds[nbounds] = { pairVarBounds.second.lb, pairVarBounds.second.ub};
		
		r = socket.writeData(&ind, sizeof(ind), &bytesTransfered);
		totalbytesTransfered += bytesTransfered;
		if(r != 0)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORNUMBER(r);
			#endif
			code = r;
			goto termination;
		}
		
		r = socket.writeData(bds, nbounds*sizeof(bds[0]), &bytesTransfered);
		totalbytesTransfered += bytesTransfered;
		if(r != 0)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORNUMBER(r);
			#endif
			code = r;
			goto termination;
		}
		
		k++;
	}
	
	#if DCT_DEBUG_MODE
		assert((int) k == request[1]);
	#endif
	
	code = 0;
	
termination:
	
	if(bytesWritten)
		*bytesWritten = totalbytesTransfered;
	
	return code;
}
#endif


#if 0
inline int DCT_readVarBounds( DCT_BaseSocket &socket, DCT_VarBounds &varsBounds )
{
	DCT_UInt32 nbounds = 0;
	int ioRetCode;
	
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
	
	
			
	ioRetCode = socket.readData(&nbounds, sizeof(nbounds), NULL);
	
	if(ioRetCode != 0)
	{
		DCT_PRINTERRORNUMBER(ioRetCode);
		return ioRetCode;
	}
	
	
	for( DCT_UInt32 i = 0; i < nbounds; i++ )
	{
		DCT_UInt32 index;
		const int sizeBounds = 2;
		double values[sizeBounds];
		DCT_Bounds bounds;
		
		ioRetCode = socket.readData(&index, sizeof(index), NULL);
		if(ioRetCode != 0)
		{
			DCT_PRINTERRORNUMBER(ioRetCode);
			return ioRetCode;
		}
		
		
		ioRetCode = socket.readData(values, sizeBounds* sizeof(values[0]), NULL);
		if(ioRetCode != 0)
		{
			DCT_PRINTERRORNUMBER(ioRetCode);
			return ioRetCode;
		}
		
		
		bounds.lb = values[0];
		bounds.ub = values[1];
		
		try
		{
			varsBounds[index] = bounds;
		}
		catch (std::bad_alloc &e)
		{
			DCT_PRINTMEMERROR;
			
			return DCT_RC_MEMORY_ERROR;
		}
		
	}
	
	return 0;
}


#endif






















}


#endif

