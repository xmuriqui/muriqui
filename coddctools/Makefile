

include ../make.inc


LIBNAME = libdctools.a
SHAREDLIBNAME = libdctools.so

OBJS = distributed.o sockets.o bbserver.o bbclient.o constants.o

INCL = $(PUGIXMLINC)




all: $(LIBNAME)

.SUFFIXES: .cpp .c .o .obj

$(LIBNAME): $(OBJS)
	$(MAKE_CMD) -C $(PUGIXMLDIR)
	$(AR_CMD)$(LIBNAME) $(OBJS)
	- $(COPY_CMD) $(PUGIXMLDIR)$(SYS_SEP)libpugixml.a  libpugixml.a;


clean:
	- $(RM_CMD) $(LIBNAME) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(ADDCXXFLAGS) $(INCL) -c $(COMPILER_OUTPUT_OBJ_FLAG)$@ $<


