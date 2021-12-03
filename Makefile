#Makefile for gnu make. Please run this using GNU Make. Iy you are in a Windows system, you can download gnu make for windows.

include make.inc


ADDCXXFLAGS =

MYSUBDIRLIBS = $(MINLPPROBDIR) $(BBLDIR) $(OPTSDIR)  $(DCTDIR)


# CHANGEME: This should be the name of your executable
EXE = muriqui
LIBNAME = libmuriqui.a
MYLIBSDIR = lib
MYINCSDIR = include

OBJS_LIB = ampl.o dataStructures.o algClasses.o tools.o extCutPlan.o outerApp.o feasibilityPump.o continuousRelax.o constants.o diving.o oaFeasibilityPump.o igma0.o solvers.o igma1.o advanced.o branchAndBound.o MRQ_bb.o extSupHipPlane.o bblpnlp_oa.o igma1bb.o igma2.o bonminHybrid.o milpCallbacks.o rens.o bblp_ecp.o bblp_esh.o gams.o testRun.o ssrounding.o

OBJS = $(OBJS_LIB) main.o


OBJS_SERVER = $(OBJS_LIB) MRQ_server.o mainServer.o
EXE_SERVER = $(EXE)dcserver

OBJS_CLIENT = mainClient.o
EXE_CLIENT = $(EXE)dcclient


#Additional libraries
ADDLIBS = 

#Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =


# Include directories

INCL = $(SPMINC) $(MINLPPROBINC) $(BBLINC) $(OPTSINC) $(DCTINC) $(NMCINC) $(ADDINCFLAGS)







# Linker flags


ifeq ("$(CXX)" , "cl")
	
	LIBS = $(BBLLIBNAME)  $(OPTSLIBNAME) $(SPMLIBNAME) $(NMCLIBNAME)  $(MINLPPROBLIBNAME)
	
	LIB_SSERVER_CLIENT = $(LIBS) $(DCTLIBNAME)

else
# g++ and icpc
	LIBS = $(BBLLIB) $(OPTSSUBLIB) $(SPMLIB)  $(NMCLIB) $(MINLPPROBLIB)     

	LIB_SSERVER_CLIENT = $(LIBS) $(DCTLIB)
endif






# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

#to tell $(MAKE_CMD)file examples is not-file-related targets
.PHONY: examples

all: lib cleanexe $(EXE)

alltools: cleanexe all examples server client 

.SUFFIXES: .cpp .c .o .obj
	

$(EXE): $(OBJS) sublibs
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@  $(OBJS)  $(ADDLIBS) $(LIBS)

clean: cleanexe cleanclient cleanserver
	- $(RM_CMD) $(OBJS) $(LIBNAME)

cleanexamples:
	$(MAKE_CMD) clean -C "examples/t1gross"
	$(MAKE_CMD) clean -C "examples/t2gross"

cleanexe:
	- $(RM_CMD) $(EXE) $(EXE).exe $(EXE_CLIENT) $(EXE_CLIENT).exe

cleanclient:
	- $(RM_CMD) $(EXE_CLIENT) $(EXE_CLIENT).exe $(OBJS_CLIENT)

cleanserver:
	- $(RM_CMD) $(EXE_SERVER) $(EXE_SERVER).exe $(OBJS_SERVER)

cleantools: cleanclient cleanserver

cleanall: clean cleansublibs cleantools cleanexamples


cleansublibs:
	$(MAKE_CMD) clean -C $(MINLPPROBDIR)
	
	$(MAKE_CMD) clean -C  $(BBLDIR)
	
	$(MAKE_CMD) clean -C  $(OPTSDIR) 
#dctools is not distributed by default. Se we can get an error and we should skip it in this case
	- $(MAKE_CMD) clean -C  $(DCTDIR) 
	
	$(MAKE_CMD) clean -C  $(NMCDIR) 


examples:
	$(MAKE_CMD) -C "examples/t1gross"
	$(MAKE_CMD) -C "examples/t2gross"

install: 
	
	- $(MKDIR_CMD) $(MYLIBSDIR)
	- $(MKDIR_CMD) $(MYINCSDIR)
	
	- $(COPY_CMD) $(MINLPPROBDIR)$(SYS_SEP)libminlpprob.a $(MYLIBSDIR)$(SYS_SEP)libminlpprob.a;
	- $(COPY_CMD) $(MINLPPROBDIR)$(SYS_SEP)*.h* $(MYINCSDIR)
	
	- $(COPY_CMD) $(BBLDIR)$(SYS_SEP)libbbl.a $(MYLIBSDIR)$(SYS_SEP)libbbl.a;
	- $(COPY_CMD) $(BBLDIR)$(SYS_SEP)*.h* $(MYINCSDIR)
	
	- $(COPY_CMD) $(OPTSDIR)$(SYS_SEP)liboptsolvers.a $(MYLIBSDIR)$(SYS_SEP)liboptsolvers.a;
	- $(COPY_CMD) $(OPTSDIR)$(SYS_SEP)*.h* $(MYINCSDIR)
	
	- $(COPY_CMD) $(NMCDIR)$(SYS_SEP)libnumcomp.a $(MYLIBSDIR)$(SYS_SEP)libnumcomp.a; 
	- $(COPY_CMD) $(NMCDIR)$(SYS_SEP)*.h* $(MYINCSDIR)
	
	- $(COPY_CMD) $(DCTDIR)$(SYS_SEP)libdctools.a $(MYLIBSDIR)$(SYS_SEP)libdctools.a; 
	- $(COPY_CMD) $(DCTDIR)$(SYS_SEP)*.h* $(MYINCSDIR)
	
	- $(COPY_CMD) $(SPMDIR)$(SYS_SEP)*.h* $(MYINCSDIR)
	
	- $(COPY_CMD) $(LIBNAME) $(MYLIBSDIR)$(SYS_SEP)$(LIBNAME)
	- $(COPY_CMD) *.h* $(MYINCSDIR)

lib: $(OBJS_LIB)
	$(AR_CMD)$(LIBNAME)  $(OBJS_LIB)


$(EXE_CLIENT): cleanclient lib sublibs dctools $(OBJS_CLIENT)
	
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $(OBJS_CLIENT) $(ADDLIBS) $(LIB_SSERVER_CLIENT)

client: dctools  $(EXE_CLIENT)


$(EXE_SERVER): lib sublibs dctools $(OBJS_SERVER)
	
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $(OBJS_SERVER) $(ADDLIBS) $(LIB_SSERVER_CLIENT)

server: dctools $(EXE_SERVER)


sublibs:
	$(MAKE_CMD) -C $(MINLPPROBDIR)
	$(MAKE_CMD) -C $(BBLDIR)
	$(MAKE_CMD) -C $(OPTSDIR) 
	$(MAKE_CMD) -C $(NMCDIR) 
	
	
dctools:
	$(MAKE_CMD) -C $(DCTDIR)
	

.cpp.o:
	$(CXX) $(CXXFLAGS) $(ADDCXXFLAGS) $(INCL) -c $(COMPILER_OUTPUT_OBJ_FLAG)$@ $<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(ADDCXXFLAGS) $(INCL) -c $(COMPILER_OUTPUT_OBJ_FLAG)$@ `$(CYGPATH_W) '$<'`
