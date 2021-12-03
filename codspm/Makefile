
LIBNAME = libspm.a

OBJS = SPM_SparseMatrix.o

CXX = icpc
CXXFLAGS = -O0 -g -Wall
CXXFLAGS = -O3 -march=native -Wall
INCL = 




all: $(LIBNAME)

.SUFFIXES: .cpp .c .o .obj

$(LIBNAME): $(OBJS)
	ar rv $(LIBNAME) $(OBJS)

clean:
	rm -rf $(LIBNAME) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<


