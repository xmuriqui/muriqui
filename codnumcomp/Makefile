

include ../make.inc


LIBNAME = libnumcomp.a
SHAREDLIBNAME = libnumcomp.so

OBJS = numComp.o

INCL = 



all: $(LIBNAME)

.SUFFIXES: .cpp .c .o .obj

$(LIBNAME): $(OBJS)
	$(AR_CMD)$(LIBNAME) $(OBJS)


clean:
	- $(RM_CMD) $(LIBNAME) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(ADDCXXFLAGS) $(INCL) -c $(COMPILER_OUTPUT_OBJ_FLAG)$@ $<


