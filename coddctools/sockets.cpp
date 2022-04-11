/*
 * Function implementation to sockets.
 * 
 * Reference: Kay RObbins, Steven Robbins, Unix System Programming. Communication, concurrency and tthreads. Prentice Hall PTR.
 * 
 */ 

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <iostream>
#include <new>

#include "DCT_sockets.hpp"



using namespace dctools;



/* Convert struct in_addr to a host name */
void DCT_AddressAndNames::addr2name(struct in_addr addr, char *name, int namelen)
{
	struct hostent *hostptr;
	
	hostptr = gethostbyaddr((char *)&addr, 4, AF_INET);
	
	if (hostptr == NULL)
		strncpy(name, inet_ntoa(addr), namelen-1);
	else
		strncpy(name, hostptr->h_name, namelen-1);
	
	name[namelen-1] = 0;
}



/* Return -1 on error, 0 on success */
int DCT_AddressAndNames::name2addr(char *name, in_addr_t *addrp)
{
	struct hostent *hp;
	
	if (isdigit((int)(*name)))
	{
		*addrp = inet_addr(name);
	}
	else
	{
		hp = gethostbyname(name);
		if (hp == NULL)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORMSG("Error to get host by name.");
			#endif
			return DCT_RC_VALUE_ERROR;
		}
		
		memcpy((char *)addrp, hp->h_addr_list[0], hp->h_length);
	}
	
	return 0;
}





DCT_BaseSocket::DCT_BaseSocket()
{
	socketfd = -1;
}


DCT_BaseSocket::~DCT_BaseSocket()
{
	closeSocket();
}


void DCT_BaseSocket::closeSocket()
{
	if(socketfd >= 0)
	{
		//close(socketfd);
		
		while ((close(socketfd) == -1) && (errno == EINTR));
		
		socketfd = -1;
	}
}


void DCT_BaseSocket::disableBlock()
{
	int options;
	
	options = fcntl(socketfd, F_GETFL, 0);
	fcntl(socketfd, F_SETFL, options | O_NONBLOCK);
}


void DCT_BaseSocket::enableBlock()
{
	int options;
	
	options=fcntl(socketfd, F_GETFL, 0);
	
	//options = options & (~O_NONBLOCK);
	
	fcntl(socketfd, F_SETFL, options & (~O_NONBLOCK));
}


int DCT_BaseSocket::getFileDescriptor()
{
	return socketfd;
}


int dctools::DCT_ignoreSigpipe()
{
	struct sigaction act;
	int r;
	
	r = sigaction(SIGPIPE, NULL, &act);
	
	if(r != 0)
	{
		int error = errno;
		
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORNUMBER(r);
			
			if(error != 0)
			{
				DCT_PRINTERRORNUMBER(error);
				DCT_PRINTERRORMSGP("errno msg: ", strerror(error));
			}
		#endif
		return DCT_errno2DCT_RETURN_CODE(error);
	}
	
	
	//SIG_DFL is the default action to be taken if we got a SIgnal. By default, it terminates our program. So, we repalce default action by SIG_IGN to ignore this signal (and do, do not terminate our program)
	if (act.sa_handler == SIG_DFL) 
	{
		act.sa_handler = SIG_IGN;
		
		int r = sigaction(SIGPIPE, &act, NULL);
		if(r != 0)
		{
			int error = errno;
			
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORNUMBER(r);
				if(error != 0)
				{
					DCT_PRINTERRORNUMBER(error);
					DCT_PRINTERRORMSGP("errno msg: ", strerror(error));
				}
			#endif
			return DCT_errno2DCT_RETURN_CODE(error);
		}
	}

	return 0;
}


//restart version of read
int DCT_BaseSocket::readData( void* buffer, size_t size, size_t* bytesRead)
{
	int code = 0;
	ssize_t retval;
	size_t mybytesRead = 0;
	size_t remainderSize;
	int myerrno;
	
	DCT_Byte *pbuffer = (DCT_Byte*) buffer;
	
	//if user writes bytes calling writeData several times, we can just receive the first bytes. So, beyond the restart, we peform a loop to read all requested bytes
	
	do
	{
		remainderSize = size - mybytesRead;
		
		retval = read(socketfd, &pbuffer[mybytesRead], remainderSize);
		myerrno = errno;
		
		if( retval > 0 )
			mybytesRead += retval;
		
	}while( (retval == -1 && myerrno == EINTR) || (retval > 0 && mybytesRead < size) ); //note, if retval is 0, connection can have been closed by other point
	
	
	
	
	//std::cout << "retval: " << retval << " errno: " << errno << "\n";
	
	if(retval == 0 && size > 0)
	{
		//we assume connection was closed by the other point
		code = DCT_RC_CLOSED_CONNECTION;
	}
	else if( retval == -1 && myerrno == EAGAIN)
	{
		code = DCT_RC_NO_DATA_TO_READ;
	}
	else if( retval == -1 && myerrno == ECONNRESET )
	{
		//we assume connection was closed by the other point
		code = DCT_RC_CLOSED_CONNECTION;
	}
	else if( mybytesRead < size )
	{
		//we could not transfer all bytes
		std::cout << DCT_PREPRINT " Error to read data. expected size: " << size << " return value: " << retval << " errno: " << myerrno << DCT_GETFILELINE << "\n";
		code = DCT_RC_UNDEFINED_ERROR;
	}
	else
	{
		#if DCT_DEBUG_MODE
			assert(mybytesRead == size);
		#endif
	}
	
	
	
	if(bytesRead)
		*bytesRead = mybytesRead;
	
	return code;
}





int DCT_BaseSocket::writeFile(FILE *file)
{
	const int sizeBuffer = 1048576; //1 MB
	size_t bytesRead, bytesWritten;
	
	int retCode;
	long unsigned int fileSize = DCT_getFileSize(file);
	
	char *buffer = NULL;
	
	
	DCT_malloc(buffer, (int) DCT_min<long unsigned int>(sizeBuffer, fileSize)) ;
	
	if(!buffer)
	{
		#if DCT_DEBUG_MODE
			DCT_PRINTMEMERROR;
		#endif
		
		retCode = DCT_RC_MEMORY_ERROR;
		goto termination;
	}
	
	//std::cout << "fileSize: " << fileSize << DCT_GETFILELINE << "\n";
	
	
	while(fileSize > 0)
	{
		long unsigned int count = DCT_min<long unsigned int>(fileSize, sizeBuffer);
		
		bytesRead = fread( buffer, sizeof(char), count, file );
		
		//std::cout << "bytesRead: " << bytesRead << " count: " << count << "\n";
		
		
		#if DCT_DEBUG_MODE
			assert( bytesRead == count );
		#endif
		
		if(bytesRead == 0)
		{ //something get wrong...
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORMSG("Error to read file");
			#endif
			
			retCode = DCT_RC_UNDEFINED_ERROR;
			goto termination;
		}
		
		writeData(buffer, bytesRead, &bytesWritten);
		
		if(bytesWritten != bytesRead)
		{
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORMSG("Error to write file");
			#endif
			
			return DCT_RC_UNDEFINED_ERROR;
			goto termination;
		}
		
		
		fileSize -= bytesRead;
	}
	
	
	retCode = 0;
	
termination:
	
	
	if(buffer)	free(buffer);
	
	return retCode;
}


int DCT_BaseSocket::writeFile(const char *fileName)
{
	int retCode;
	FILE *file = fopen(fileName, "rb");
	
	
	if(!file)
	{
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORMSGP("Failure to open file: ", fileName);
		#endif
		return DCT_RC_VALUE_ERROR;
	}
	
	
	retCode = writeFile(file);
	
	
	
	if(file) fclose(file);
	
	return retCode;
}


void DCT_BaseSocket::setFileDescriptor(int fd)
{
	socketfd = fd;
}



int DCT_BaseSocket::writeData(const void *buf, size_t size, size_t *byteswritten)
{
	const char *bufp = (const char*) buf;
	int code = 0;
	int myerrno;
	size_t bytestowrite;
	ssize_t nbyteswritten;
	size_t mybyteswritten;
	
	
	mybyteswritten = 0;
	
	
	for (bytestowrite = size; bytestowrite > 0; bufp += nbyteswritten, bytestowrite -= nbyteswritten)
	{
		nbyteswritten = write(socketfd, bufp, bytestowrite);
		myerrno = errno;
		
		if (nbyteswritten == -1)
		{
			if( (myerrno == EINTR) || (myerrno == EAGAIN) )
			{
				std::cerr << DCT_PREPRINT "Warning: error to write Bytes. nbyteswritten: " << nbyteswritten << " errno: " << myerrno << " . Trying again!\n";
			}
			else
			{
				#if DCT_DEBUG_MODE
					std::cerr << DCT_PREPRINT "Error to write Bytes. nbyteswritten: " << nbyteswritten << " errno: " << myerrno << DCT_GETFILELINE << "\n";
				#endif
				
				if( myerrno == EPIPE )
					code = DCT_RC_CLOSED_CONNECTION;
				else
					code = DCT_RC_IO_ERROR;
				goto termination;
			}
		}
		
		if (nbyteswritten < 0)
			nbyteswritten = 0;
		
		
		mybyteswritten += nbyteswritten;
	}
	
	
	if( mybyteswritten == 0 && size > 0)
	{
		//we assume connection was closed by the other point
		code = DCT_RC_CLOSED_CONNECTION;
	}
	else if( mybyteswritten < size )
	{
		//we could not transfer all bytes
		code = DCT_RC_UNDEFINED_ERROR;
	}
	else
	{
		#if DCT_DEBUG_MODE
			assert(mybyteswritten == size);
		#endif
	}
	
	
termination:
	
	if(byteswritten)
		*byteswritten = mybyteswritten;
	
	return code;
}




DCT_ServerConnection::DCT_ServerConnection():DCT_BaseSocket()
{
}


DCT_ServerConnection::~DCT_ServerConnection()
{
}





DCT_ServerSocket::DCT_ServerSocket():DCT_BaseSocket()
{
}



DCT_ServerSocket::~DCT_ServerSocket()
{
}



/*
 *                           acceptConnection
 * Wait for a connection request from a host on a specified port.
 *
 * parameters:
 *      socketNumber = file descriptor previously bound to listening port
 *      hostn = a buffer that will hold the name of the remote host
 *      hostnsize = size of hostn buffer
 * returns:  a communication file descriptor on success
 *              hostn is filled with the name of the remote host.
 *           -1 on error with errno set
 *
 * comments: This function is used by the server to wait for a
 * communication.  It blocks until a remote request is received
 * from the port bound to the given file descriptor.
 * hostn is filled with an ASCII string containing the remote
 * host name.  It must point to a buffer of size at least hostnsize.
 * If the name does not fit, as much of the name as is possible is put
 * into the buffer.
 * If hostn is NULL or hostnsize <= 0, no hostname is copied.
 */
int DCT_ServerSocket::acceptConnection(char *hostn, int hostnsize, DCT_ServerConnection &connection)
{
	int len = sizeof(struct sockaddr);
	int retval;
	struct sockaddr_in netclient;
	
	
	do
	{
		retval = accept(socketfd, (struct sockaddr *)(&netclient), (socklen_t *) &len);
		
	} while( (retval == -1) && (errno == EINTR) );
	
	
	//while ( ((retval = accept(socketNumber, (struct sockaddr *)(&netclient), (socklen_t *) &len)) == -1)  &&  (errno == EINTR) );


	if ((retval == -1) || (hostn == NULL) || (hostnsize <= 0))
	{
		int error = errno;
		
		if(error != 0)
		{
			DCT_PRINTERRORNUMBER(error);
			DCT_PRINTERRORMSGP("errno message: ", strerror(error));
		}
		
		return retval;
	}
	
	DCT_AddressAndNames::addr2name(netclient.sin_addr, hostn, hostnsize);
	
	connection.setFileDescriptor(retval);
	
	return 0;
}



/*
 * Open a new socket
 * 
 * port is the number of port to dind the socket
 * 
 * socketNumber is the file descriptor of new socket
 * 
 */
int DCT_ServerSocket::openSocket(const int port)
{
	struct sockaddr_in server; //we need a sockaddr to specify options to binding proccess. However, is it easir work on sockaddr_in and convert after to sockaddr when necessary.
	
	const int MAXBACKLOG = 50;
	int socketNumber;
	int r;
	
	
	//We have to disbale default action for sigpipe, since it will terminate the program if we got error for write in a closed pipe.
	r = DCT_ignoreSigpipe();
	if(r != 0)
	{
		int error = errno;
		
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORNUMBER(r);
			if(error != 0)
			{
				DCT_PRINTERRORNUMBER(error);
				DCT_PRINTERRORMSGP("errno msg: ", strerror(error));
			}
		#endif
		return DCT_errno2DCT_RETURN_CODE(error);
	}
	
	
	/*
	 * First parameter to socket funcion is the domain, indicate the protocol family to be used. We use AF_INET, indicating we are using IPv4.
	 *
	 * Second argument is type. We adopt SOCK_STREAM to use a sequenced, relaiable, connection-based byte streams
	 * 
	 * Third parameter (protocol) specifies the protocol. In general, there is only one protocol available for each case (e.g. TCP for SOCK_STREAM or UDP fpr SPCK_DGRAM). So, protocol is usually 0
	 * 
	 * Returns the number of descriptor of our socket.
	 * 
	*/
	socketNumber = socket(AF_INET, SOCK_STREAM, 0);
	if( socketNumber < 0 )
	{
		int error = errno;
		
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORMSG("Error to open socket.");
			if(error != 0)
			{
				DCT_PRINTERRORNUMBER(error);
				DCT_PRINTERRORMSGP("errno msg: ", strerror(error));
			}
		#endif
		return DCT_errno2DCT_RETURN_CODE(error);
	}
	
	
	/*
	 * Set socket option
	 * 
	 * First argument is the socket.
	 * 
	 * Second argument specifies the protocol level. SOL_SOCKET indicates we are seting a option in the socket itself
	 * 
	 * third argument, is the SO_REUSEADDR, indicating we are set/changing this option. SO_REUSEADDR permit (or no) the server to be restarted immediately, using the same port.
	 * 
	 * Fourth argument is new value for obpiton being changed, and Fifth argument is the size of new value for the optin being changed. We pass1, indicating true, i.e., we allow reuse the port.
	 * 
	 */ 
	{
		int one = 1; //this value will be used as true
		
		r = setsockopt(socketNumber, SOL_SOCKET, SO_REUSEADDR, &one, sizeof(one));
		
		if(r != 0)
		{
			int error = errno;
			
			#if DCT_DEBUG_MODE
				DCT_PRINTERRORNUMBER(r);
				if(error != 0)
				{
					DCT_PRINTERRORNUMBER(error);
					DCT_PRINTERRORMSGP("errno msg: ", strerror(error));
				}
			#endif
			return DCT_errno2DCT_RETURN_CODE(error);
		}
	}
	
	
	//options for binding proccess
	server.sin_family = AF_INET;  //for internet communication, we use AF_NET
	server.sin_addr.s_addr = htonl(INADDR_ANY); //sin_addr is IP address. INADDR_ANY is Address to accept any incoming messages. htonl  is a function to  convert address to network byte order.
	server.sin_port = htons((short)port); //port number. htons convert the port number to network byte order (We could have problem about the order of bytes in different architetures if we do not use htonl and htons)
	
	
	/*
	 * binding procces: associates the handle for a socket communication endpoint with a specific logical network connection.
	 * 
	 * First argument is the socket.
	 * 
	 * Second argument is a structure contating informations (parameters);
	 * 
	 * Thir argument is the number of bytes of second parameter
	 * 
	 */
	r = bind(socketNumber, (struct sockaddr *)&server, sizeof(server));
	if(r != 0)
	{
		int error = errno;
		
		
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORNUMBER(r);
			if(error != 0)
			{
				DCT_PRINTERRORNUMBER(error);
				DCT_PRINTERRORMSGP("errno msg: ", strerror(error));
			}
			if(error == EADDRINUSE)
			{
				DCT_PRINTERRORMSG("check if this port is not being used by some application!");
			}
		#endif
		
		
		while ((close(socketNumber) == -1) && (errno == EINTR));
		errno = error;
		
		return DCT_errno2DCT_RETURN_CODE(error);
	}
	
	
	/*
	 * listen proccess: causes the underlying system network infraestructure to allocate queues to hold pending requests.
	 * 
	 * Second argument is the maximum number of connections before a conection be rejected
	 * 
	 */ 
	r = listen(socketNumber, MAXBACKLOG);
	if(r != 0)
	{
		int error = errno;
		
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORNUMBER(r);
			if(error != 0)
			{
				DCT_PRINTERRORNUMBER(error);
				DCT_PRINTERRORMSGP("errno msg: ", strerror(error));
			}
		#endif
		
		while ((close(socketNumber) == -1) && (errno == EINTR));
		errno = error;
		
		return DCT_errno2DCT_RETURN_CODE(error);
	}
	
	
	socketfd = socketNumber;
	return 0;
}




DCT_ClientSocket::DCT_ClientSocket():DCT_BaseSocket()
{
}


DCT_ClientSocket::~DCT_ClientSocket()
{
}


int DCT_ClientSocket::connectToAServer(int port, char *hostn)
{
	int error;
	int retval;
	struct sockaddr_in server;
	int sock;
	fd_set sockset;
	
	
	//check if hostn is a valid address and put the binary address in server
	if( DCT_AddressAndNames::name2addr( hostn, &(server.sin_addr.s_addr)) == -1)
	{
		errno = EINVAL;
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORMSG("Error to resolve address name.");
		#endif
		return DCT_RC_UNDEFINED_ERROR;
	}
	
	server.sin_port = htons((short)port); //htons convert the port number to network byte order (We could have problem about the order of bytes in different architetures if we do not use htonl and htons)
	server.sin_family = AF_INET; //for internet communication, we use AF_NET
	
	
	
	//We have to disbale default action for sigpipe, since it will terminate the program if we got error for write in a closed pipe.
	retval = DCT_ignoreSigpipe();
	if( retval != 0 )
	{
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORNUMBER(retval);
		#endif
		return DCT_RC_UNDEFINED_ERROR;
	}
	
	
	/*
	 * First parameter to socket funcion is the domain, indicate the protocol family to be used. We use AF_INET, indicating we are using IPv4.
	 *
	 * Second argument is type. We adopt SOCK_STREAM to use a sequenced, relaiable, connection-based byte streams
	 * 
	 * Third parameter (protocol) specifies the protocol. In general, there is only one protocol available for each case (e.g. TCP for SOCK_STREAM or UDP fpr SPCK_DGRAM). So, protocol is usually 0
	 * 
	 * Returns the number of descriptor of our socket.
	 * 
	*/
	sock = socket(AF_INET, SOCK_STREAM, 0);
	
	if( sock  == -1 )
	{
		#if DCT_DEBUG_MODE
			DCT_PRINTERRORMSG("Error to open socket.");
		#endif
		
		return DCT_RC_UNDEFINED_ERROR;
	}
	
	
	retval = connect(sock, (struct sockaddr *) &server, sizeof(server));
	
	if ( (retval == -1) && ((errno == EINTR) || (errno == EALREADY)) )
	{
		//if conect was interrupted by a signal, we should not restart synce network subsystem has already initiated the TCP-3 way handshake. So, we have to use select to detect if 
		
		do
		{
			FD_ZERO(&sockset);
			FD_SET(sock, &sockset);
			
			retval = select(sock+1, NULL, &sockset, NULL, NULL);
			
		}while( (retval == -1) && (errno == EINTR) );
		
		/*FD_ZERO(&sockset);
		FD_SET(sock, &sockset);
		while ( ((retval = select(sock+1, NULL, &sockset, NULL, NULL)) == -1) &&
				(errno == EINTR) ) {
			FD_ZERO(&sockset);
			FD_SET(sock, &sockset);
		}*/
	}
	
	
	if(retval == -1)
	{
		error = errno;
		while ((close(sock) == -1) && (errno == EINTR));
		errno = error;
		
		#if DCT_DEBUG_MODE
			std::cout << DCT_PREPRINT " Error to stablish a connection. erro: " << error << DCT_GETFILELINE << "\n";
			//DCT_PRINTERRORMSG("Error to create a connection.");
		#endif
		
		return DCT_RC_UNDEFINED_ERROR;
	}
	
	
	
	socketfd = sock;
	return 0;
}




//write the size of string, the string and, finally, '\0' character. Size includes '\0'
int dctools:: DCT_writeStringWithSizeAnd0( DCT_BaseSocket &socket, const char *data, size_t *bytesWritten)
{
	const DCT_UInt32 size = strlen(data) + 1;
	int r;
	size_t myBytesWritten;
	
	
	if(bytesWritten)
		bytesWritten = 0;
	
	r = socket.writeData( &size, sizeof(size), &myBytesWritten );
	if(bytesWritten)
		bytesWritten += myBytesWritten;
	
	if(r != 0)
	{
		DCT_PRINTERRORNUMBER(r);
		return DCT_RC_SERVER_ERROR;
	}
	
	r = socket.writeData(data, size * sizeof(data[0]), &myBytesWritten );
	if(bytesWritten)
		bytesWritten += myBytesWritten;
	
	if(r != 0)
	{
		DCT_PRINTERRORNUMBER(r);
		return DCT_RC_SERVER_ERROR;
	}
	
	//just checking if last character is '\0'
	#if DCT_DEBUG_MODE
		assert( data[size-1] == '\0' );
	#endif
	
	return 0;
}



#if 0
DCT_VarBoundsWriter::DCT_VarBoundsWriter()
{
	initialize();
}

DCT_VarBoundsWriter::~DCT_VarBoundsWriter()
{
	desallocate();
}

void DCT_VarBoundsWriter::desallocate()
{
	DCT_secFree(buffer);
	sizeBuffer = 0;
}

void DCT_VarBoundsWriter::initialize()
{
	sizeBuffer = 0;
	buffer = NULL;
}

int DCT_VarBoundsWriter::reallocateBuffer(unsigned int newSize)
{
	if( newSize == 0 )
	{
		DCT_secFree(buffer);
	}
	else
	{
		int r = DCT_realloc(buffer, newSize);
		if(r != 0)
		{
			DCT_PRINTERRORNUMBER(r);
			return r;
		}
	}
	
	sizeBuffer = newSize;
	
	return 0;
}


int DCT_VarBoundsWriter::write(DCT_BaseSocket& socket, const DCT_VarBounds& varsBounds, size_t* bytesWritten)
{
	/*
	 * Here, we write:
	 *  1 - The DCT_VAR_BOUNDS  code (int32)
	 *  2 - The number n of tuple of bounds tranfesred (uint32)
	 *  3 - n tuple of < index (uint32), lb (double), ub (double)>
	 *
	 * 
	 * we try send all data in just one call to speed up node transfer between servers and clients
	 */ 
	
	const unsigned int varBoundsSize = varsBounds.size();
	const unsigned int totalBytes = sizeof(DCT_Int32) + sizeof(DCT_UInt32) + varBoundsSize * ( sizeof(DCT_UInt32) + 2*sizeof(double) ) ;
	DCT_UInt32 *beg;
	double *bounds;
	
	int r;
	DCT_Byte *pBuffer;
	
	
	if(bytesWritten)
		*bytesWritten = 0;
	
	if( totalBytes > sizeBuffer )
	{
		r = reallocateBuffer(2*totalBytes); //we take advatnage the realloc to allocate a little more memory
		if(r != 0)
		{
			DCT_PRINTERRORNUMBER(r);
			return r;
		}
	}
	
	
	beg = (DCT_UInt32*) buffer;
	beg[0] = DCT_CC_VAR_BOUNDS;
	beg[1] = varBoundsSize;
	
	pBuffer = (DCT_Byte*) &beg[2];
	
	for(const auto &pairVarBounds : varsBounds)
	{
		const DCT_UInt32 ind = pairVarBounds.first;
		const double lb = pairVarBounds.second.lb;
		const double ub = pairVarBounds.second.ub;
		
		beg = (DCT_UInt32*) pBuffer;
		beg[0] = ind;
		
		bounds = (double*) &beg[1];
		bounds[0] = lb;
		bounds[1] = ub;
		
		pBuffer = (DCT_Byte*) &bounds[2];
	}
	
	//now, pBuffer points to end of buffer
	#if DCT_DEBUG_MODE
		assert( &buffer[totalBytes-1] == pBuffer-1 );
	#endif
	
	r = socket.writeData(buffer, totalBytes, bytesWritten);
	if(r != 0)
	{
		DCT_PRINTERRORNUMBER(r);
		return r;
	}
	
	return 0;
}

#endif




