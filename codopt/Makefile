include ../make.inc

LIBNAME = liboptsolvers.a

OBJS = solvers.o mosek.o cplex.o gurobi.o glpk.o xpress.o mynlp.o ipopt.o worhp.o iquad.o knitro.o algencan.o tools.o cbc.o osi.o optizelle.o



INCL =  $(OPTSINC)  -I../$(MINLPPROBDIR)  -I../$(SPMDIR)




all: $(LIBNAME)

.SUFFIXES: .cpp .c .o .obj

$(LIBNAME): $(OBJS)
	$(AR_CMD)$(LIBNAME) $(OBJS)

clean:
	- $(RM_CMD) $(LIBNAME) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(ADDCXXFLAGS) $(INCL) -c $(COMPILER_OUTPUT_OBJ_FLAG)$@ $<


