

include ../make.inc


LIBNAME = libbbl.a
SHAREDLIBNAME = libbbl.so

OBJS = BBL_branchAndBound.o bb.o node.o


ADDCXXFLAGS = -std=c++11



all: $(LIBNAME)

.SUFFIXES: .cpp .c .o .obj

$(LIBNAME): $(OBJS)
	$(AR_CMD)$(LIBNAME) $(OBJS)

clean:
	- $(RM_CMD) $(LIBNAME) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(ADDCXXFLAGS)  -c $(COMPILER_OUTPUT_OBJ_FLAG)$@ $<


