


# C++ Compiler command (lower case, please)
CXX = g++ 
#CXX = clang
#CXX = icpc

PARALEL = -fopenmp  #for g++
#PARALEL = -liomp5   #for icpc


# C++ Compiler options
#CXXFLAGS = -O0 -m64 -g -Wall     
CXXFLAGS = -O3 -m64  -march=native  -Wall 

#for microsoft vcms
#CXXFLAGS = -O2 -nologo -favor:INTEL64
#CXXFLAGS = -O2 -nologo -favor:AMD64


# additional C++ Compiler options for linking
CXXLINKFLAGS = 



#######################################################
########       Operating System setting         #######
#######################################################


ifeq ($(OS),Windows_NT)
	MAKE_CMD = make  #yes, use gnu-make
	#we use strip because gnu make interpret \ as a special caracter
	SYS_SEP =$(strip \)
	RM_CMD = del
	#no link command at ms-dos
	COPY_CMD = copy   
	MKDIR_CMD = mkdir
	AR_CMD = lib -out:
	COMPILER_OUTPUT_OBJ_FLAG = -Fo
else
	MAKE_CMD = make
	SYS_SEP =/
	RM_CMD = rm -rf
	COPY_CMD = cp
	MKDIR_CMD = mkdir -p
	#please, let a blank space at the end of AR_CMD for non windows systems
	AR_CMD = ar rv 
	COMPILER_OUTPUT_OBJ_FLAG = -o
endif





#######################################################
########       Modeling languages setting       #######
#######################################################

#AMPL Solver Library (ASL) variables to compilation. You can let them empty if you want compile whitout ASL
#ASLDIR = /opt/asl/
#ASLINC = -I$(ASLDIR)
#ASLLIB = -L$(ASLDIR)/ -lamplsolver -ldl -lgomp -fopenmp   -lirc -limf

#for windows
#ASLDIR = C:\Users\wendel\Documents\muriqui\solvers
#ASLINC = -I$(ASLDIR)
#ASLLIB = "$(ASLDIR)\amplsolv.lib"



#GAMS variables to compilation. You can let them empty if you want compile whitout GAMS
#GAMSDIR = /opt/gams/
#GAMSINC = -I$(GAMSDIR)/apifiles/C/api/
#GAMSLIB = $(GAMSDIR)/apifiles/C/api/gmomcc.c  $(GAMSDIR)/apifiles/C/api/gevmcc.c  -ldl






#######################################################
##########      Lapack and Blas setting      ##########
#######################################################


#BLAS variables to compilation
#BLASDIR = /opt/BLAS
#BLASLIB_MKL = -L/opt/intel/lib/intel64/ -L/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/  -lmkl_intel_lp64  -lmkl_sequential   -lmkl_core  -liomp5  -lpthread -lm
#BLASLIB = -L$(BLASDIR) -lblas


#Lapack variables to compilation
#LAPACKDIR = /opt/lapack
#LAPACKLIB_MKL = $(BLASLIB_MKL)
#LAPACKLIB = $(LAPACKLIB_MKL) # -L$(LAPACKDIR) -llapack -lblas -lrt -lgfortran 






#######################################################
############        Solvers setting       #############
#######################################################



#Algencan variables to compilation. You can let them empty if you want compile whitout Algencan
#ALGENCANDIR = /opt/algencan
#ALGENCANINC = 
#for gcc, you do not need use -lifport nor -lifcoremt
#ALGENCANLIB = -L$(ALGENCANDIR)/lib   -lalgencan  -lifport  -lifcoremt  -lgfortran -lm


#CBC variables to compilation. You can let them empty if you want compile whitout CBC
#CBCDIR = /opt/cbc
#CBCINC = -I$(CBCDIR)/include/coin/
#CBCLIB = -L$(CBCDIR)/lib/  -lCbcSolver -lCbc -lOsiCbc   -lCgl  -lOsiClp -lOsi -lCoinUtils  -lClpSolver -lClp


#CPLEX variables to compilation. You can let them empty if you want compile whitout CPLEX
#CPLEXDIR = /opt/cplex
#CPLEXINC = -I$(CPLEXDIR)/include/ilcplex/
#CPLEXLIB = -L$(CPLEXDIR)/lib/x86-64_linux/static_pic/ -lcplex -lpthread -pthread -lm -ldl

#for windows
#CPLEXDIR = "C:\Program Files\IBM\ILOG\CPLEX_Studio201\cplex"
#CPLEXINC = -I$(CPLEXDIR)\include\ilcplex
#CPLEXDIRLIB = "C:\Program Files\IBM\ILOG\CPLEX_Studio201\cplex\lib\x64_windows_msvc14\stat_mda"
#CPLEXLIB = $(CPLEXDIRLIB)\cplex2010.lib  
#$(CPLEXDIRLIB)\ilocplex.lib


#CSDP variables to compilation. You can let them empty if you want compile whitout CSDP
#CSDPDIR = /opt/Csdp-6.1.1
#CSDPINC = -I$(CSDPDIR)/include/ -DNOSHORTS
#CSDPLIB = -L$(CSDPDIR)/lib/ -lsdp


#GLPK variables to compilation. You can let them empty if you want compile whitout GLPK
#GLPKDIR = /opt/glpk
#GLPKINC = -I$(GLPKDIR)/include/
#GLPKLIB = -L$(GLPKDIR)/lib  -lglpk


#GUROBI variables to compilation. You can let them empty if you want compile whitout GUROBI
#GUROBIDIR = /opt/gurobi
#GUROBIINC = -I$(GUROBIDIR)/include/
#GUROBILIB = -L$(GUROBIDIR)/lib  -lgurobi -lgfortran $(PARALEL) -lpthread -lm


#IPOPT variables to compilation. You can let them empty if you want compile whitout IPOPT
#IPOPTDIR = /opt/Ipopt
#IPOPTINC = -I$(IPOPTDIR)/include/coin-or/
#IPOPTLIB = -L$(IPOPTDIR)/lib -lipopt -L$(IPOPTDIR)/ThirdParty-HSL/lib -lcoinhsl  -L/opt/intel/lib/intel64/ -lirc -limf -lirng -lintlc -lsvml  $(LAPACKLIB)  $(BLASLIB)


#KNITRO variables to compilation. You can let them empty if you want compile whitout KNITRO
#KNITRODIR = /opt/knitro
#KNITROINC = -I$(KNITRODIR)/include
#KNITROLIB = -L$(KNITRODIR)/lib  -lknitro
#KNITROLIB = -L/opt/ampl  -lknitro1210


#MOSEK variables to compilation. You can let them empty if you want compile whitout MOSEK
#MOSEKDIR = /opt/mosek/8/tools/platform/linux64x86
#MOSEKINC = -I$(MOSEKDIR)/h/
#MOSEKLIB = -L$(MOSEKDIR)/bin -lmosek64 -Wl,-rpath-link,$(MOSEKDIR)/bin  -pthread -lc -lm


#Optizelle variables to compilation. You can let them empty if you want compile whitout Optizelle
#OPTIZELLEDIR = /opt/optizelle
#OPTIZELLEINC = -I$(OPTIZELLEDIR)/include/    -I$(OPTIZELLEDIR)/thirdparty/include
#OPTIZELLELIB = -L$(OPTIZELLEDIR)/lib  -loptizelle  -L$(OPTIZELLEDIR)/thirdparty/lib -ljsoncpp  -lgfortran


#WORHP variables to compilation. You can let them empty if you want compile whitout WORHP
#WORHPDIR = /opt/worhp
#WORHPINC = -I$(WORHPDIR)/include/worhp/
#WORHPLIB = -L$(WORHPDIR)/lib  -lworhp


#XPRESS variables to compilation. You can let them empty if you want compile whitout XPRESS
#XPRESSDIR = /opt/xpressmp
#XPRESSINC = -I$(XPRESSDIR)/include/
#XPRESSLIB = -L$(XPRESSDIR)/lib/  -lxprs -lxprl -lm








#######################################################
############      XML Library settings    #############
#######################################################



#pugixml variables to compilation
PUGIXMLDIR = ./pugixml/
PUGIXMLINC = -I$(PUGIXMLDIR)
PUGIXMLLIB = -L$(PUGIXMLDIR) -lpugixml












#######################################################
############       Internal Packages      #############
#######################################################





#our internal third packages (only change if you know what you are doing)
SPMDIR = codspm
SPMINC = -I$(SPMDIR)
#SPMLIBNAME =
#SPMLIB = #-L$(SPMDIR) -lspm

MINLPPROBDIR = codminlpp
MINLPPROBINC = -I$(MINLPPROBDIR) $(SPMINC)  $(ASLINC)  $(GAMSINC)
MINLPPROBLIBNAME = $(MINLPPROBDIR)/libminlpprob.a  $(ASLLIB) $(GAMSLIB)
MINLPPROBSUBLIB = -L$(MINLPPROBDIR) -lminlpprob  $(SPMLIB) $(ASLLIB)
MINLPPROBLIB = $(MINLPPROBSUBLIB) $(GAMSLIB)		#that is just to include GAMSLIB

BBLDIR = codbb
BBLINC = -I$(BBLDIR)
BBLLIBNAME = $(BBLDIR)/libbbl.a
BBLLIB = -L$(BBLDIR) -lbbl  -lpthread


OPTSDIR = codopt
OPTSINC = -I$(OPTSDIR) $(SPMINC) $(MINLPPROBINC) $(CPLEXINC) $(MOSEKINC) $(GUROBIINC) $(CBCINC) $(XPRESSINC) $(GLPKINC) $(XPRESSINC)  $(WORHPINC) $(KNITROINC) $(IPOPTINC) $(ALGENCANINC) $(OPTIZELLEINC)
OPTSLIBNAME = $(OPTSDIR)/liboptsolvers.a $(CPLEXLIB)  $(GUROBILIB) $(GLPKLIB) $(CBCLIB) $(XPRESSLIB)  $(WORHPLIB)  $(ALGENCANLIB) $(MOSEKLIB) $(IPOPTLIB) $(KNITROLIB)  $(OPTIZELLELIB)

#ipopt should came before knitro due to mkl implementations
OPTSSUBLIB = -L$(OPTSDIR) -loptsolvers   $(CPLEXLIB)  $(GUROBILIB) $(GLPKLIB) $(CBCLIB) $(XPRESSLIB)  $(WORHPLIB)  $(ALGENCANLIB) $(MOSEKLIB) $(IPOPTLIB) $(KNITROLIB)  $(OPTIZELLELIB) $(MINLPPROBSUBLIB) $(SPMLIB) 
OPTSLIB = $(OPTSSUBLIB) $(MINLPPROBLIB)   #that is just to include GAMSLIB


DCTDIR = coddctools
DCTINC = -I$(DCTDIR)
DCTLIBNAME =  $(DCTDIR)/libdctools.a  $(DCTDIR)/libpugixml.a
DCTLIB = -L$(DCTDIR) -ldctools -lpugixml


NMCDIR = codnumcomp
NMCINC = -I$(NMCDIR)
NMCLIBNAME = $(NMCDIR)/libnumcomp.a
NMCLIB = -L$(NMCDIR) -lnumcomp


#IQUADDIR = ~/researches/codiquad/lastversion
#IQUADINC = -I$(IQUADDIR) $(BBLINC) $(OPTSINC) $(MINLPPROBINC) $(SPMINC)
#IQUADSUBLIB = -L$(IQUADDIR)/libs -liquad  $(CSDPLIB) $(BBLLIB) $(OPTSLIB) $(SPMLIB)  $(MINLPPROBSUBLIB)
#IQUADLIB = $(IQUADSUBLIB)  $(MINLPPROBLIB)	#that is just to include GAMSLIB


MURIQUIDIR = ~/researches/codminlp/lastversion
MURIQUIONLYINC = $(BBLINC) $(OPTSINC) $(MINLPPROBINC) $(SPMINC)
MURIQUIINC = -I$(MURIQUIDIR)/include $(MURIQUIONLYINC)
MURIQUIONLYLIB = -lmuriqui $(BBLLIB) $(OPTSSUBLIB)  $(NMCLIB) $(SPMLIB)  $(MINLPPROBSUBLIB)
MURIQUISUBLIB = -L$(MURIQUIDIR)/lib  $(MURIQUIONLYLIB)
MURIQUILIB = $(MURIQUISUBLIB)  $(MINLPPROBLIB)	#that is just to include GAMSLIB


MYSTEINERDIR = ~/researches/steiner/code/base_library/lastversion/
MYSTEINERINC = -I$(MYSTEINERDIR)
MYSTEINERLIB = -L$(MYSTEINERDIR) -lestp


SAMBADIR = ~/researches/steiner/code/samba/lastversion
SAMBAINC = -I$(MYSTEINERINC) -I$(SAMBADIR)
SAMBALIB = -L$(SAMBADIR) -lsamba $(MYSTEINERLIB) $(MURIQUILIB)


