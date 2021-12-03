include ../make.inc


LIBNAME = libminlpprob.a

OBJS = MIP_minlpProblem.o MIP_ampl.o MIP_preprocessing.o MIP_preprocessing_old.o MIP_tools.o MIP_gams.o


INCL = -I../$(SPMDIR) $(ASLINC) $(GAMSINC)




all: $(LIBNAME)

.SUFFIXES: .cpp .c .o .obj

$(LIBNAME): $(OBJS)
	$(AR_CMD)$(LIBNAME) $(OBJS)

clean:
	- $(RM_CMD) $(LIBNAME) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(ADDCXXFLAGS) $(INCL) -c $(COMPILER_OUTPUT_OBJ_FLAG)$@ $<


