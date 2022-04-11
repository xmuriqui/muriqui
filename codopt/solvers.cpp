


#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <ctime>

#include <new>
#include <iostream>
#include <exception>

#include "MIP_minlpProblem.hpp"

#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"



#define OPT_MAX_SIZE_OFSTRING_TO_READ_PARAMS 512




using namespace std;
using namespace optsolvers;
using namespace newspm;
using namespace minlpproblem;




bool optsolvers::OPT_isSolverAvailable(const int solverCode)
{
    //#define OPT_SOLVER_CODE_COUNTER 0
    
    switch(solverCode)
    {
        case OPT_GLPK:
            //#define OPT_SOLVER_CODE_COUNTER 1
        #if OPT_HAVE_GLPK
            return true;
        #else
            return false;
        #endif
        
        
        case OPT_CBC:
            //#define OPT_SOLVER_CODE_COUNTER 2
        #if OPT_HAVE_CBC
            return true;
        #else
            return false;
        #endif
        
        
        case OPT_CPLEX:
            //#define OPT_SOLVER_CODE_COUNTER 3
        #if OPT_HAVE_CPLEX
            return true;
        #else
            return false;
        #endif
        
        
        case OPT_GUROBI:
            //#define OPT_SOLVER_CODE_COUNTER 4
        #if OPT_HAVE_GUROBI
            return true;
        #else
            return false;
        #endif
        
        
        case OPT_XPRESS:
            //#define OPT_SOLVER_CODE_COUNTER 5
        #if OPT_HAVE_XPRESS
            return true;
        #else
            return false;
        #endif
        
        
        case OPT_MOSEK:
            //#define OPT_SOLVER_CODE_COUNTER 6
        #if OPT_HAVE_MOSEK
            return true;
        #else
            return false;
        #endif
        
        
        case OPT_KNITRO:
            //#define OPT_SOLVER_CODE_COUNTER 7
        #if OPT_HAVE_KNITRO
            return true;
        #else
            return false;
        #endif
        
        
        case OPT_IPOPT:
            //#define OPT_SOLVER_CODE_COUNTER 8
        #if OPT_HAVE_IPOPT
            return true;
        #else
            return false;
        #endif
        
        
        case OPT_ALGENCAN:
            //#define OPT_SOLVER_CODE_COUNTER 9
        #if OPT_HAVE_ALGENCAN
            return true;
        #else
            return false;
        #endif
        
        
        case OPT_WORHP:
            //#define OPT_SOLVER_CODE_COUNTER 10
        #if OPT_HAVE_WORHP
            return true;
        #else
            return false;
        #endif
        
            
        case OPT_OPTIZELLE:
            //#define OPT_SOLVER_CODE_COUNTER 11
        #if OPT_HAVE_OPTIZELLE
            return true;
        #else
            return false;
        #endif
        
        
        case OPT_IQUAD:
            #define OPT_SOLVER_CODE_COUNTER 12
        #if OPT_HAVE_IQUAD
            return true;
        #else
            return false;
        #endif
            
        case OPT_UNDEFINEDSOLVER:
            return false;
        
        default:
            OPT_PRINTERRORMSGP("Invalid solver code: ", solverCode);
            return false;
            
        #if OPT_SOLVER_CODE_COUNTER < OPT_NUMBER_OF_SOLVERS
            Compilation error! You forget to include some solver code here! (do not comment this)
        #endif
    }
}




int OPT_GeneralSolverParams::storeIntegerParameter(const char *name, const int64_t value)
{
    try
    {
        intParams[name] = value;
    }
    catch (std::bad_alloc &e)
    {
        return OPT_MEMORY_ERROR;
    }
    
    return 0;
}


int OPT_GeneralSolverParams::storeDoubleParameter(const char* name, const double value)
{
    
    try
    {
        dblParams[name] = value;
    }
    catch (std::bad_alloc &e)
    {
        return OPT_MEMORY_ERROR;
    }
    
    return 0;
}


int OPT_GeneralSolverParams::storeStringParameter(const char *name, const char *value)
{
    
    try
    {
        strParams[name] = value;
    }
    catch (std::bad_alloc &e)
    {
        return OPT_MEMORY_ERROR;
    }
    
    return 0;
}


void OPT_GeneralSolverParams::desallocate()
{
    intParams.clear();
    dblParams.clear();
    strParams.clear();
}


OPT_GeneralSolverParams::~OPT_GeneralSolverParams()
{
    desallocate();
}



void OPT_GeneralSolverParams::print( std::ostream &out)
{
    for(auto &pair: intParams)
        out << pair.first << ": " << pair.second << "\n";
    
    for(auto &pair: dblParams)
        out << pair.first << ": " << pair.second << "\n";
    
    for(auto &pair: strParams)
        out << pair.first << ": " << pair.second << "\n";
}



void OPT_GeneralSolverParams::removeIntegerParameter( const char *name)
{
    if( intParams.count(name) >= 1 )
        intParams.erase(name);
}


void OPT_GeneralSolverParams::removeDoubleParameter(const char *name)
{
    if( dblParams.count(name) >= 1 )
        dblParams.erase(name);
}


void OPT_GeneralSolverParams::removeStringParameter( const char *name)
{
    if( strParams.count(name) >= 1 )
        strParams.erase(name);
}


int OPT_GeneralSolverParams::storeParametersFromFile( const char *fileName, const bool printErrorMsgs, const bool printFileOpenError)
{
    OPT_GeneralSolverParamsParameterSetter setter(this);
    
    return OPT_readParametersWithTypeFromFile(fileName, printErrorMsgs, printFileOpenError, setter);
}




OPT_Solver::OPT_Solver()
{
    initialize();
}


OPT_Solver::~OPT_Solver()
{
    deallocateMemory();
}


void OPT_Solver::initialize()
{
    naux = 0;
    maux = 0;
    auxIndex = NULL;
    auxIndex2 = NULL;
    auxValues = NULL;
    auxValues2 = NULL;
    

    feasSol = false;
    retCode = OPT_UNDEFINED;
    objValue = dualObjValue = NAN;
    sol = NULL;
    constr = NULL;
    dualSolC = NULL;
    dualSolV = NULL;
    
    numberOfWarningsByIterLimit = 0;
    maxNumberOfWarningsByIterLimit = OPT_DEFAULT_MAX_NUMBER_OF_WARNINGS_BY_ITERATION_LIMIT;
}



int OPT_Solver::addConstraints(const int nm)
{
    int r;

    r = allocateConstrStructures( maux + nm );

    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    //maux += nm;
    
    r = __addConstraints(nm);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    
    return r;
}



int OPT_Solver::addVariables(const int nn, const bool initFree)
{
    int r;

    r = allocateVarStructures( naux + nn );
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return r;
    }

    //naux += nn;

    return __addVariables(nn, initFree);
}



int OPT_Solver:: allocateAuxStructures(const int size)
{
    int r;
    /*int *auxi;
    double *auxd;

    auxd = (double *) realloc( auxValues, 2 * size * sizeof(double) );
    if( !auxd )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    auxValues = auxd; */
    
    r = OPT_realloc(auxValues, 2*size);
    OPT_IFERRORRETURN(r, r);
    
    auxValues2 = &auxValues[size];
    
    
    /*auxi = (int *) realloc( auxIndex, 2 * size * sizeof(int) );
    if( !auxi )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    auxIndex = auxi; */
    
    r = OPT_realloc(auxIndex, 2*size);
    OPT_IFERRORRETURN(r, r);
    
    auxIndex2 = &auxIndex[size];
    
    
    return 0;
}


int OPT_Solver::allocateVarStructures(const int n)
{
    int r;
    //int *auxi;
    //double *auxd;


    /*auxd = (double *) realloc( sol, n * sizeof(double) );
    if( !auxd )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    sol = auxd;*/

    r = OPT_realloc(sol, n);
    OPT_IFERRORRETURN(r, r);
    

    /*auxd = (double *) realloc( dualSolV, 2*n * sizeof(double) );
    if( !auxd )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    dualSolV = auxd; */
    
    r = OPT_realloc(dualSolV, 2*n);
    OPT_IFERRORRETURN(r, r);
    

    if( n > naux )
    {
        const int s = n - naux;
        //naux is the old value
        OPT_setAllArray<double>(s, &sol[naux], NAN );
        OPT_setAllArray<double>(s, &dualSolV[naux], NAN);
    }


    if( n > maux )
    {
        int r = allocateAuxStructures(n);

        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
    }


    naux = n;



    return 0;
}


int OPT_Solver::checkConstraintIndex(const int index)
{
    int r, m;
    
    r = getNumberOfConstraints(m);
    OPT_IFERRORRETURN(r, r);
    
    if( index < 0 || index >= m )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORMSGP("Invalid index for constraint: ", index);
        #endif
        return OPT_BAD_INPUT;
    }
    
    return 0;
}


int OPT_Solver::checkConstraintIndices(const unsigned int size, const int *indices)
{
    int r, m;
    
    r = getNumberOfConstraints(m);
    OPT_IFERRORRETURN(r, r);
    
    for(unsigned int i = 0; i < size; i++)
    {
        if( indices[i] < 0 || indices[i] >= m )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORMSGP("Invalid index for variable: ", indices[i]);
            #endif
            return OPT_BAD_INPUT;
        }
    }
    
    return 0;
}


int OPT_Solver::checkVariableIndex(const int index)
{
    int r, n;
    
    r = getNumberOfVars(n);
    OPT_IFERRORRETURN(r, r);
    
    if( index < 0 || index >= n )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORMSGP("Invalid index for variable: ", index);
        #endif
        return OPT_BAD_INPUT;
    }
    
    return 0;
}


int OPT_Solver::checkVariableIndices(const unsigned int size, const int *indices)
{
    int r, n;
    
    r = getNumberOfVars(n);
    OPT_IFERRORRETURN(r, r);
    
    for(unsigned int i = 0; i < size; i++)
    {
        if( indices[i] < 0 || indices[i] >= n )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORMSGP("Invalid index for variable: ", indices[i]);
            #endif
            return OPT_BAD_INPUT;
        }
    }
    
    return 0;
}


int OPT_Solver::allocateConstrStructures(const int m)
{
   //double *auxd;

    int r;

    /*auxd = (double *) realloc( constr, m * sizeof(double) );
    if( !auxd )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    constr = auxd;*/
    
    r = OPT_realloc(constr, m);
    OPT_IFERRORRETURN(r, r);
    

    /*auxd = (double *) realloc( dualSolC, m * sizeof(double) );
    if( !auxd )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    dualSolC = auxd; */
    
    
    r = OPT_realloc(dualSolC, m);
    OPT_IFERRORRETURN(r, r);


    if( m > maux )
    {
        const int s = m - maux;
        //maux is the old value
        OPT_setAllArray<double>(s, &constr[maux], NAN );
        OPT_setAllArray<double>(s, &dualSolC[maux], NAN );
    }


    if( m > naux )
    {
        const int r = allocateAuxStructures(m);

        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
    }


    maux = m;

    return 0;
}




void OPT_Solver::deallocateMemory()
{
    OPT_secFree( auxIndex );
    OPT_secFree( auxValues );

    deallocateSol();
}


void OPT_Solver::deallocateSol()
{
    OPT_secFree( sol );
    OPT_secFree( constr );
    OPT_secFree( dualSolC );
    OPT_secFree( dualSolV );
}


bool OPT_Solver::getMinusLambdaOnLagran()
{
    return true;
}


/*void OPT_Solver::desallocateSolverEnv()
{
} */


void OPT_Solver::getDualSolution( double *dualConstrs, double *dualVarBounds, const bool correctSignal )
{
    const bool minusLambdaOnLagran = getMinusLambdaOnLagran();
    
    if(dualConstrs)
    {
        int m = 0;
        getNumberOfConstraints(m);
        
        if( minusLambdaOnLagran || !correctSignal )
        {
            OPT_copyArray(m, dualSolC, dualConstrs);
        }
        else
        {
            OPT_copyArrayTimesFactor(m, -1, dualSolC, dualConstrs);
        }
    }
    
    
    if(dualVarBounds)
    {
        int n = 0;
        getNumberOfVars(n);
        
        if( minusLambdaOnLagran || !correctSignal )
        {
            OPT_copyArray(2*n, dualSolV, dualVarBounds);
        }
        else
        {
            OPT_copyArrayTimesFactor(2*n, -1, dualSolV, dualVarBounds);
        }
    }
}


double OPT_Solver::getObjValue()
{
    return objValue;
}

double OPT_Solver::getDualObjValue()
{
    return dualObjValue;
}

void OPT_Solver::getSolution( double *solution, double *constraints )
{
    if(solution)
    {
        int n = 0;
        getNumberOfVars(n);
        
        OPT_copyArray(n, sol, solution);
    }
    
    
    if(constraints)
    {
        int m = 0;
        getNumberOfConstraints(m);
        
        OPT_copyArray(m, constr, constraints);
    }
}


std::string OPT_Solver::getSolverName()
{
    return OPT_getSolverName( getSolverCode() );
}


bool OPT_Solver::isMyNLPClass()
{
    return false;
}


int OPT_Solver::removeConstraintsByRange(const int begin, const int end)
{
    const int ninds = end - begin + 1;
    int *indices = auxIndex2; //we use auxIndex2 because some solvers like cplex and glpk already use auxIndex to remove constraints.
    
    OPT_setSequentialValuesArray(indices, begin, end);
    
    const int r = removeConstraints(ninds, indices);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    return 0;
}


void OPT_Solver::printDblParamErrorMsg(const int error, const char* param, const double value )
{
    std::cerr << OPT_PREPRINT "Error " << error << " at setting " << getSolverName()  << " double parameter: " << param << " at the value: " << value << std::endl;
}



void OPT_Solver::printIntParamErrorMsg(const int error, const char *param, const int value )
{
    std::cerr << OPT_PREPRINT "Error " << error << " at setting " << getSolverName() << " integer parameter: " << param << " at the value: " << value << std::endl;
}



void OPT_Solver::printStrParamErrorMsg(const int error, const char* param, const char* value )
{
    std::cerr << OPT_PREPRINT "Error " << error << " at setting " << getSolverName()  << " string parameter: " << param << " at the value: " << value << std::endl;
}




void OPT_Solver::resetSol()
{
    feasSol = false;
    retCode = OPT_UNDEFINED;
    objValue = dualObjValue = NAN;
    origSolverRetCode = INT_MAX;

    OPT_setAllArray<double>(naux, sol, NAN);
    OPT_setAllArray<double>(naux, dualSolV, NAN);
    OPT_setAllArray<double>(maux, constr, NAN);
    OPT_setAllArray<double>(maux, dualSolC, NAN);
}



//by default, we set max cpu time. Multhithread solvers classes should overwrite this method
int OPT_Solver::setMaxTime(const double time)
{
    return setMaxCPUTime(time);
}



int OPT_Solver:: setParameters( OPT_GeneralSolverParams &params )
{
    int r, code = 0;


    for(auto &pair : params.intParams)
    {
        r = setIntegerParameter(pair.first.c_str(),  pair.second);
        if ( r != 0 )
        {
            code = OPT_BAD_INPUT;
        }
    }

    for(auto &pair : params.dblParams)
    {
        r = setDoubleParameter(pair.first.c_str(),  pair.second);
        if ( r != 0 )
        {
            code = OPT_BAD_INPUT;
        }
    }

    for(auto &pair : params.strParams)
    {
        r = setStringParameter(pair.first.c_str(),  pair.second.c_str());

        if ( r != 0 )
        {
            code = OPT_BAD_INPUT;
        }
    }

    return code;
}








void OPT_Solver::setThreadNumber(const unsigned int threadNumber)
{
    //we just only use it to do evaluations in NLP solver...
}


int OPT_Solver::setVariableType( const int index, const OPT_VARTYPE varType )
{
    #if OPT_DEBUG_MODE
        OPT_PRINTNONSUPMSG;
    #endif
    return OPT_OPERATION_NOT_SUPPORTED;
}


int OPT_Solver::setObjCutLowerBound(const double objLBound)
{
    #if OPT_DEBUG_MODE
        //OPT_PRINTNONSUPMSG;
    #endif
    return OPT_OPERATION_NOT_SUPPORTED;
}


int OPT_Solver::setObjCutUpperBound(const double objUBound)
{
    #if OPT_DEBUG_MODE
        //OPT_PRINTNONSUPMSG;
    #endif
    return OPT_OPERATION_NOT_SUPPORTED;
}


int OPT_Solver::solveAndGetTime(double *cpuTime, double *clockTime, const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
{
    int r;
    
    clock_t initialCpuClockTime;
    double initialTime;
    
    
    if(cpuTime)
        initialCpuClockTime = clock();
    
    if(clockTime)
        initialTime = OPT_getTime();
    
    
    r = solve(resetSol, storeSol, storeConstrs, storeDualSol);
    
    
    if(cpuTime)
        *cpuTime = OPT_calcCPUTtime(initialCpuClockTime);
    
    if(clockTime)
        *clockTime = OPT_getTime() - initialTime;
    
    
    return r;
}


int OPT_Solver::warmUp()
{
    #if OPT_DEBUG_MODE
        OPT_PRINTNONSUPMSG;
    #endif
    return OPT_OPERATION_NOT_IMPLEMENTED;
}



OPT_LPSolver::OPT_LPSolver():OPT_Solver()
{
}


OPT_LPSolver::~OPT_LPSolver()
{
}



int OPT_LPSolver::addVariablesFromProb(const MIP_MINLPProb& prob, const bool setObj, const bool setVarBounds, const bool setVarType, const int naddvars)
{
    auto n = prob.getNumberOfVars();
    int r;
    
    r = addVariables(n + naddvars);
    OPT_IFERRORRETURN(r, r);
    
    
    if( setVarBounds )
    {
        const double *lx = prob.lx;
        const double *ux = prob.ux;
        
        setnVariablesBounds(n, lx, ux);
        OPT_IFERRORRETURN(r, r);
    }
    
    
    if(setVarType && prob.getNumberOfIntegerVars() > 0)
    {
        for(decltype(n) i = 0; i < n; i++)
        {
            MIP_VARTYPE vt;
            
            prob.getVariableType(i, vt);
            
            if( vt == MIP_VT_INTEGER )
            {
                r = setVariableType(i, OPT_VT_INTEGER );
                OPT_IFERRORRETURN(r, r);
            }
        }
    }
    
    if( setObj )
    {
        //if( !setStrictPartOnly || (prob.Q.getNumberOfElements() == 0 && !prob.hasObjNLTerm() ) )
        {
            const double of = prob.getObjFactor();
            const double aof = of; //OPT_abs(of);
            
            
            r = setObjSense(optsolvers::OPT_MINIMIZE); 
            OPT_IFERRORRETURN(r, r);
            
            
            setObjConstant( aof * prob.getObjConstant() );
            
            
            if( prob.hasLinCoefObj() && aof != 0.0 )
            {
                //int nzs = 0;
                const double *c = prob.c;
                
                //prob.getObjLinearCoefs(auxValues);
                
                /*for(int i = 0; i < n; i++)
                {
                    if( auxValues[i] != 0.0)
                    {
                        auxIndex2[nzs] = i;
                        auxValues2[nzs] = auxValues[i];
                        nzs++;
                    }
                }
                
                if( aof != 1.0 )
                {
                    for(int i = 0; i < nzs; i++)
                        auxValues2[i] *= aof;
                }
                
                
                r = setObjLinearCoefs( nzs, auxIndex2, auxValues2 );
                OPT_IFERRORRETURN(r, r); */
                
                
                for(int i = 0; i < n; i++)
                {
                    if( c[i] != 0.0 )
                    {
                        r = setObjLinearCoef(i, aof*c[i]);
                        OPT_IFERRORRETURN(r, r);
                    }
                }
            }
        }
    }
    
    return 0;
}





int OPT_LPSolver::copyLPPartFrom(OPT_LPSolver &other, const bool copyObj, const bool copyConstrs, const bool copyVarBounds, const bool copyVarTypes)
{
    int n, m;
    int r, retCode;
    OPT_OPTSENSE objSense;
    double objConstant;
    
    int *rows = NULL, *cols = NULL;
    double *values = NULL;
    
    
    r = other.getNumberOfVars(n);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    r = other.getNumberOfConstraints(m);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    deallocateMemory();
    
    r = initSolverEnv(m, n);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    r = addVariables(n);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    //std::cout << "Copying objective\n";
    
    if( copyObj )
    {
        r = other.getObjSense(objSense);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        r = other.getObjConstant(objConstant);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        r = setObjSense(objSense);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        setObjConstant(objConstant);
       
        
        //setting linear objective coeficients
        for(int i = 0; i < n; i++)
        {
            double v;
            
            r = other.getObjLinearCoef(i, v);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            r = setObjLinearCoef(i, v);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        }
    }
    
    //std::cout << "Copying var bounds\n";
    
    if( copyVarBounds )
    {
        //setting var bounds
        for(int i = 0; i < n; i++)
        {
            double lx, ux;
            r = other.getVariableBounds(i, lx, ux);
            OPT_IFERRORRETURN(r, r);
            
            r = setVariableBounds(i, lx, ux);
            OPT_IFERRORRETURN(r, r);
        }
    }
    
    //std::cout << "Copying var types\n";
    
    if( copyVarTypes )
    {
        //setting var types
        for(int i = 0; i < n; i++)
        {
            OPT_VARTYPE vt;
            
            r = other.getVariableType(i, vt);
            OPT_IFERRORRETURN(r, r);
            
            //if we got an error of operation not supported, we do not care
            r = setVariableType(i, vt);
            if(r != OPT_OPERATION_NOT_SUPPORTED)
                OPT_IFERRORRETURN(r, r);
        }
    }
    
    
    if( copyConstrs )
    {
        r = addConstraints(m);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        {
            //here, we try copy all coeficients togheter to optimize the copy.
            int nzA = 0;
            int nzA2 = 0;
            
            //std::cout << "Counting matrix A and copyng RHS\n";
            
            for(decltype(m) i = 0; i < m; i++)
            {
                double lc, uc;
                
                r = other.getConstraintBounds(i, lc, uc);
                OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
                
                r = setConstraintBounds(i, lc, uc);
                OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            }
            
            r = other.getNumberOfLinearCoefsInConstraints(nzA);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            OPT_malloc(rows, nzA);
            OPT_malloc(cols, nzA);
            OPT_malloc(values, nzA);
            
            OPT_IFMEMERRORGOTOLABEL(!rows || !cols || !values, retCode, termination);
            
            
            //std::cout << "Reading matriz A\n";
            
            r = other.getLinearCoefsInConstraints(nzA2, rows, cols, values);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            
            //std::cout << "Storing matriz A\n";
            
            //setting the entire matrix
            r = setConstraintsLinearCoefs(nzA2, rows, cols, values);
            OPT_IFERRORRETURN(r, r);
            
        }
        
    }
    
    retCode = 0;
    
termination:
            
    if(rows)    free(rows);
    if(cols)    free(cols);
    if(values)  free(values);
    
    return retCode;
}



int OPT_LPSolver::generateModelFile(const char* fileName)
{
    #if OPT_DEBUG_MODE
        OPT_PRINTNONSUPMSG;
    #endif
    return OPT_OPERATION_NOT_IMPLEMENTED;
}


int OPT_LPSolver::getFullConstraintLinearPart(const int constrIndex, double *values)
{
    int n, nzs;
    
    int r = getConstraintLinearPart( constrIndex, nzs, auxIndex2, auxValues2 );
    OPT_IFERRORRETURN(r, r);
    
    r = getNumberOfVars(n);
    OPT_IFERRORRETURN(r, r);
    
    
    if( n != nzs )
        OPT_setAllArray(n, values, 0.0);
    
    for(int i = 0; i < nzs; i++)
        values[ auxIndex2[i] ] = auxValues2[i];
    
    return 0;
}


int OPT_LPSolver::getLinearCoefsInConstraints(int &nzs, int *rows, int *cols, double *values)
{
    int r, m;
    
    nzs = 0;
    
    r = getNumberOfConstraints(m);
    OPT_IFERRORRETURN(r, r);
    
    for(int i = 0; i < m; i++)
    {
        int nz = 0;
        
        r = getConstraintLinearPart(i, nz, &cols[nzs], &values[nzs] );
        OPT_IFERRORRETURN(r, r);
        
        OPT_setAllArray(nz, &rows[nzs], i);
        
        nzs += nz;
    }
    
    return 0;
}



int OPT_LPSolver:: getNumberOfLinearCoefsInConstraints( int &nzs)
{
    int r, m;
    
    nzs = 0;
    
    r = getNumberOfConstraints(m);
    OPT_IFERRORRETURN(r, r);
    
    for(int i = 0; i < m; i++)
    {
        int nz = 0;
        r = getNumberOfConstraintLinearCoefs(i, nz);
        OPT_IFERRORRETURN(r, r);
        
        nzs += nz;
    }
        
    return 0;
}




OPT_SOLVERTYPE OPT_LPSolver::getSolverType()
{
    return OPT_LP;
}



/*int OPT_LPSolver::resetLinearConstraintPart(const int constrIndex, spm::SPM_SparseRow<double> &row )
{
    const int nz = row.getStructureAndValues( auxIndex2, auxValues2);
    
    return resetLinearConstraintPart(constrIndex, nz, auxIndex2, auxValues2 );
}*/



int OPT_LPSolver::setLinearColumn( const int varIndex, const int nzs, const int *rows, const double *values)
{
    int ret = 0;


    for(int i = 0; i < nzs; i++)
    {
        ret += setConstraintLinearCoef( rows[i], varIndex, values[i] );
    }

    if( ret != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(ret);
        #endif
        return OPT_UNDEFINED_ERROR;
    }

    return 0;
}



/*int OPT_LPSolver::setConstraintBounds( const int index, const double lb, const double ub )
{
    const int ret = setConstraintLowerBound( index, lb) + setConstraintUpperBound(index, ub);

    if( ret != 0 )
        return OPT_UNDEFINED_ERROR;

    return 0;
} */



int OPT_LPSolver::setConstraintsLinearCoefs( const int nz, const int* rows, const int* cols, const double* values )
{
    int ret = 0;

    for(int i = 0; i < nz; i++)
        ret += setConstraintLinearCoef( rows[i], cols[i], values[i] );


    if(ret != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(ret);
        #endif
        
        return OPT_UNDEFINED_ERROR;
    }

    return 0;
}



int OPT_LPSolver::setConstraintLinearCoefs( const int constrIndex, const int nz, const int *cols, const double *values)
{
    int ret = 0;

    for(int i = 0; i < nz; i++)
        ret += setConstraintLinearCoef( constrIndex, cols[i], values[i] );


    if( ret != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(ret);
        #endif
        
        return OPT_UNDEFINED_ERROR;
    }

    return 0;
}



int OPT_LPSolver::setLinearObjAndConstraintsFrom(const minlpproblem::MIP_MINLPProb &prob, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int naddvars , const int naddconstrs)
{
    auto n = prob.getNumberOfVars();
    auto m = prob.getNumberOfConstraints();
    
    const auto ml = prob.getNumberOfLinearConstraints();
    
    int r, code;
    int *rowsA = NULL, *colsA = NULL;
    double *valuesA = NULL;
    
    
    deallocateSolverEnv();
    deallocateMemory();
    
    r = initSolverEnv( 2*(n + naddvars), 2*(ml + naddconstrs) );
    OPT_IFERRORGOTOLABEL(r, code, r, termination);
    
    
    r = addVariablesFromProb(prob, setObj, setVarBounds, setVarType, naddvars);
    
    r = addConstraints(ml + naddconstrs);
    OPT_IFERRORGOTOLABEL(r, code, r, termination);
    
    
    if( setConstrs && ml > 0 )
    {
        const MIP_SparseMatrix &A = prob.A;
        const unsigned int maxnz = A.getNumberOfElements();
        unsigned int mynz = 0, myrow = 0;
        
        if( maxnz > 0 )
        {
            const double *lc = prob.lc;
            const double *uc = prob.uc;
            
            OPT_malloc(rowsA, maxnz);
            OPT_malloc(colsA, maxnz);
            OPT_malloc(valuesA, maxnz);
            OPT_IFMEMERRORGOTOLABEL( !rowsA || !colsA || !valuesA, code, termination );
            
            for( decltype(m) i = 0; i < m; i++)
            {
                {
                    bool linear;
                    prob.isLinearConstraint(i, linear);
                    if( !linear )
                        continue;
                }
                
                r = setConstraintBounds( myrow, lc[i], uc[i] );
                OPT_IFERRORGOTOLABEL(r, code, r, termination);
                
                const unsigned int nzrow = A.getNumberOfElementsAtRow(i);
                
                OPT_setAllArray<int>( nzrow, &rowsA[mynz], myrow );
                
                OPT_copyArray( nzrow, A[i], &colsA[mynz] );
                
                OPT_copyArray( nzrow, A(i), &valuesA[mynz] );
                
                
                mynz += nzrow;
                myrow++;
            }
            
            #if OPT_DEBUG_MODE
                assert( mynz <= maxnz );
                assert( myrow == (unsigned) ml );
            #endif
            
            r = setConstraintsLinearCoefs(mynz, rowsA, colsA, valuesA);
            OPT_IFERRORGOTOLABEL(r, code, r, termination);
            
        }
        
    }
    
    
    
    code = 0;
    
termination:
    
    if(rowsA)   free(rowsA);
    if(colsA)   free(colsA);
    if(valuesA) free(valuesA);
    
    return code;
}



int OPT_LPSolver::setObjLinearCoefs( const int nz, const int *cols, const double *values )
{
    int ret = 0;

    for(int i = 0; i < nz; i++)
        ret += setObjLinearCoef( cols[i], values[i] );

    if( ret != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(ret);
        #endif
            
        return OPT_UNDEFINED_ERROR;
    }

    return 0;
}



int OPT_LPSolver::setObjLinearPart( const int n, const double *values )
{
    int ret = 0;

    for(int i = 0; i < n; i++)
        ret += setObjLinearCoef(i, values[i]);

    if( ret != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(ret);
        #endif
            
        return OPT_UNDEFINED_ERROR;
    }

    return 0;
}


//we prefer do not declare this method like virtual because we assume if it is called from a LP pointer, user only want LP part of the problem...
int OPT_LPSolver::setProblemFrom( const MIP_MINLPProb& prob, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int naddvars, const int naddconstrs )
{
    const int n = prob.getNumberOfVars();
    const int m = prob.getNumberOfConstraints();
    
    //const int size = OPT_max(n, m);
    
    //we have to allocate arrays for indices and values. We cannot use auxIndex and auxValues from the class because they are already used by the methods being called...
    
    
    int r, code;
    
    
    deallocateSolverEnv();
    deallocateMemory();
    
    
    
    code = initSolverEnv( 2*(n + naddvars), 2*(m + naddconstrs) );
    if( code != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(code);
        #endif
        
        goto termination;
    }
    
    #if 0
    {
    code = addVariables(n + naddvars);
    if( code != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(code);
        #endif
        
        goto termination;
    }
    
    
    if(setVarBounds)
    {
        prob.getVariableLowerBounds(auxValues);
        prob.getVariableUpperBounds(auxValues2);
        
        for(int i = 0; i < n; i++)
        {
            double lb = auxValues[i] > -MIP_INFINITY ? auxValues[i] : -OPT_INFINITY;
            double ub = auxValues2[i] < MIP_INFINITY ? auxValues2[i] : OPT_INFINITY;
            
            code = setVariableBounds(i, lb, ub);
            
            if( code != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(code);
                #endif
                
                goto termination;
            }
        }
    }
    
    
    if(setVarType && prob.getNumberOfIntegerVars() > 0)
    {
        for(int i = 0; i < n; i++)
        {
            MIP_VARTYPE vt;
            
            prob.getVariableType(i, vt);
            
            
            if( vt == MIP_VT_INTEGER )
            {
                //code = setVariableType(i, vt == MIP_VT_INTEGER ? OPT_VT_INTEGER : OPT_VT_CONTINUOUS);
                code = setVariableType(i, OPT_VT_INTEGER );
                
                if( code != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(code);
                    #endif
                    
                    goto termination;
                }
            }
        }
    }
    
    
    if( setObj )
    {
        //if( !setStrictPartOnly || (prob.Q.getNumberOfElements() == 0 && !prob.hasObjNLTerm() ) )
        {
            const double of = prob.getObjFactor();
            const double aof = of; //OPT_abs(of);
            
            
            code = setObjSense(optsolvers::OPT_MINIMIZE); //setObjSense( of >= 0 ? OPT_MINIMIZE : OPT_MAXIMIZE );
            if( code != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(code);
                #endif
                
                goto termination;
            }
            
            
            setObjConstant( aof * prob.getObjConstant() );
            
            
            if( prob.hasLinCoefObj() )
            {
                int nzs = 0;
                
                prob.getObjLinearCoefs(auxValues);
                
                for(int i = 0; i < n; i++)
                {
                    if( auxValues[i] != 0.0)
                    {
                        auxIndex2[nzs] = i;
                        auxValues2[nzs] = auxValues[i];
                        nzs++;
                    }
                }
                
                if( aof != 1.0 )
                {
                    for(int i = 0; i < nzs; i++)
                        auxValues2[i] *= aof;
                }
                
                
                code = setObjLinearCoefs( nzs, auxIndex2, auxValues2 );
                if( code != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(code);
                    #endif
                    
                    goto termination;
                }
            }
        }
    }
    }
    #endif
    
    
    r = addVariablesFromProb(prob, setObj, setVarBounds, setVarType, naddvars);
    OPT_IFERRORGOTOLABEL(r, code, r, termination);
    
    
    if( setConstrs )
    {
        int k, rm = 0;
        const MIP_SparseMatrix &A = prob.A;
        
        
        /*if( setStrictPartOnly )
        {
            for(int i = 0; i < m; i++)
            //for( MIP_SparseMatrixRowIndexIterator Ait = A.beginRowIndex(); *Ait < m; ++Ait)
            {
                //const int i = *Ait; 
                if( prob.isLinearConstraint(i) )
                    rm++;
            }
        }
        else */
        {
            rm = m;
        }
        
        
        code = addConstraints(rm + naddconstrs);
        if( code != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(code);
            #endif
            goto termination;
        }
        
        
        {
            const int nzs = A.getNumberOfElements();
            if(nzs > 0)
            {
                int r;
                
                int *rows = (int*) malloc( nzs * sizeof(int) );
                if(!rows)
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTMEMERROR;
                    #endif
                    goto termination;
                }
                
                A.getStructure(rows, NULL);
                
                r = setConstraintsLinearCoefs(nzs, rows, A[0], A(0));
                
                free(rows);
                if(r != 0)
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(r);
                    #endif
                    goto termination;
                }
            }
        }
        
        
        
        
        k = 0;
        //for( MIP_SparseMatrixRowIndexIterator Ait = A.beginRowIndex(); *Ait < m; ++Ait) we cannot run on row index iterator because we still need set lb and ub for nonlinear constraint...
        for(int i = 0; i < m; i++)
        {
            //const int i = *Ait;
            
            /*if( setStrictPartOnly && !prob.isLinearConstraint(i) )
                continue;*/
            
            
            double lb, ub;
            
            /*const int nzs = A.getNumberOfElementsAtRow(i);
            
            if( nzs > 0 )
            {
                const int* const rcols = A[i];
                const double* const rvalues = A(i);
                
                code = resetLinearConstraintPart( k, nzs, rcols, rvalues );
                
                if( code != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(code);
                    #endif
                    goto termination;
                }
            } */
            
            prob.getConstraintBounds(i, lb, ub);
            
            if(lb <= -MIP_INFINITY)
                lb = -OPT_INFINITY;
            
            if(ub >= MIP_INFINITY)
                ub = OPT_INFINITY;
            
            code = setConstraintBounds( k, lb, ub);
            if( code != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(code);
                #endif
                
                goto termination;
            }
            
            
            k++;
        }
        
        
        #if OPT_DEBUG_MODE
            assert(k == rm);
        #endif
        
        
        /*prob.getConstraintLowerBounds(auxValues);
        prob.getConstraintUpperBounds(auxValues2);
        
        for(int i = 0; i < m; i++)
        {
            double lb = auxValues[i] > -MIP_INFINITY ? auxValues[i] : -OPT_INFINITY;
            double ub = auxValues2[i] < MIP_INFINITY ? auxValues2[i] : OPT_INFINITY;
            
            code = setConstraintBounds( i, lb, ub);
            if( code != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(code);
                #endif
                
                goto termination;
            }
        } */
    
    }
    
    
    
    
    code = 0;
    
termination:
    
    if(code != 0)
    {
        deallocateSolverEnv();
        deallocateMemory();
    }
    
    return code;
}



/* int OPT_LPSolver::setVariableBounds( const int index, const double lb, const double ub )
{
    const int ret = setVariableLowerBound(index, lb) + setVariableUpperBound(index, ub);

    return ret == 0 ? 0 : OPT_UNDEFINED_ERROR;
}
*/


int OPT_LPSolver::setnVariablesBounds( const int n, const double *lb, const double *ub )
{
    int r, code = 0;
    
    for(int i = 0; i < n; i++)
    {
        r = setVariableBounds( i, lb[i], ub[i] );
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                    
            code = r;
        }
    }
    
    return code;
}


int OPT_LPSolver::setVariablesBounds( const int ninds, const int *inds, const double *lb, const double *ub )
{
    int r, code = 0;
    
    for(int i = 0; i < ninds; i++)
    {
        r = setVariableBounds( inds[i], lb[i], ub[i] );
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                    
            code = r;
        }
    }
    
    return code;
}







OPT_QPSolver::OPT_QPSolver():OPT_LPSolver()
{
}


OPT_QPSolver::~OPT_QPSolver()
{
}


int OPT_QPSolver::copyQPPartFrom( OPT_QPSolver &other, const bool copyObj, const bool copyConstrs, const bool copyVarBounds, const bool copyVarTypes)
{
    int r, nzs, nzs2, retCode;
    int *rows = NULL, *cols = NULL;
    double *vals = NULL;
    
    r = copyLPPartFrom(other, copyObj, copyConstrs, copyVarBounds, copyVarTypes);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    
    if( copyObj )
    {
        r = other.getNumberOfQuadObjTerms(nzs);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        OPT_malloc(rows, nzs);
        OPT_malloc(cols, nzs);
        OPT_malloc(vals, nzs);
        
        OPT_IFMEMERRORGOTOLABEL(!rows || !cols || !vals , retCode, termination);
        
        r = other.getObjQuadPart( nzs2, rows, cols, vals);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        #if MRQ_DEBUG_MODE
            assert(nzs == nzs2);
        #endif
        
        r = setObjQuadMatrix(nzs, rows, cols, vals);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    }
    
    retCode = 0;
    
termination:
    
    if(rows)	free(rows);
    if(cols)	free(cols);
    if(vals)	free(vals);
    
    return retCode;
}



OPT_SOLVERTYPE OPT_QPSolver::getSolverType()
{
    return OPT_QP;
}


int OPT_QPSolver::setProblemFrom( const MIP_MINLPProb& prob, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int naddvars, const int naddconstrs )
{
    int code;
    double *values = NULL;
    
    code = OPT_LPSolver::setProblemFrom(prob, setObj, setConstrs, setVarBounds, setVarType, naddvars, naddconstrs);
    if( code != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(code);
        #endif
        
        goto termination;
    }
    
    
    if( setObj )
    {
        const MIP_SparseMatrix &OQ = prob.Q;
        const unsigned int nzs = OQ.getNumberOfElements(); //prob.getNumberOfObjQuadTerms( );
        
        #if 0
        if( nzs > 0 )
        {
            rows = (int *) malloc(nzs* sizeof(int));
            
            if( !rows )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                
                code = OPT_MEMORY_ERROR;
                goto termination;
            }
            
            //cols = &rows[nzs];
            
            
            Q.getStructure(rows, NULL);
            
            const int* const pcols = Q[0];
            const double* pvalues = Q(0);
            
            //prob.getObjQuadCoefsMatrix( &nzs, rows, cols, values );
            
            /*
            #if OPT_DEBUG_MODE
                if( r != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(r);
                    #endif
                        
                    code = OPT_SOLVER_ERROR;
                    goto termination;
                }
            #endif*/
            
            
            if( prob.objFactor != 1.0 )
            {
                const double of = prob.objFactor;
                
                values= (double *) malloc( nzs* sizeof(double) );
                if( !values )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTMEMERROR;
                    #endif
                    
                    code = OPT_MEMORY_ERROR;
                    goto termination;
                }
                
                #pragma ivdep
                #pragma GCC ivdep
                for(unsigned int i = 0; i < nzs; i++)
                    values[i] = pvalues[i]*of;
                
                pvalues = values;
            }
            
            code = setQuadObjMatrix( nzs, rows, pcols, pvalues );
            if( code != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(code);
                #endif
                
                goto termination;
            }
        }
        #endif
        
        
        if( nzs > 0 )
        {
            const double of = prob.objFactor;
            double *pvalues = OQ(0);
            
            if( of != 1.0 )
            {
                values = (double *) malloc( nzs* sizeof(double) );
                if( !values )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTMEMERROR;
                    #endif
                    code = OPT_MEMORY_ERROR;
                    goto termination;
                }
                
                #pragma ivdep
                #pragma GCC ivdep
                for(unsigned int i = 0; i < nzs; i++)
                    values[i] = of*pvalues[i];
                
                pvalues = values;
            }
            
            code = setObjQuadMatrix( (int*) OQ.offset, OQ[0], pvalues );
            if( code != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(code);
                #endif
                
                goto termination;
            }
        }
        
    }
    
    
    code = 0;
    
termination:
    
    if( code != 0 )
    {
        deallocateSolverEnv();
        deallocateMemory();
    }
    
    if( values )	free(values);
    
    
    return code;
}



int OPT_QPSolver::setObjQuadMatrix( const int nzs, const int *rows, const int *cols, const double *values )
{
    int ret = 0;

    for(int i = 0; i < nzs; i++)
        ret += setObjQuadCoef( rows[i], cols[i], values[i] );

    return ret == 0? 0 : OPT_UNDEFINED_ERROR;
}



int OPT_QPSolver::setObjQuadMatrix(const int *rowStart, const int *cols, const double *values)
{
    int k, code;
    int n;
    getNumberOfVars(n);
    
    const int nzs = rowStart[n]-rowStart[0];
    int *rows = NULL;
    
    
    rows = (int *) malloc( nzs * sizeof(int) );
    if( !rows )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    k = 0;
    for(int i = 0; i < n; i++)
    {
        const int rnzs = rowStart[i+1] - rowStart[i];
        
        OPT_setAllArray(rnzs, &rows[k], i);
        k += rnzs;
    }
    
    #if OPT_DEBUG_MODE
        assert(k == nzs);
    #endif
    
    code = setObjQuadMatrix(nzs, rows, &cols[rowStart[0]], &values[rowStart[0]] );
    if( code != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(code);
        #endif
        goto termination;
    }
    
    code = 0;
    
termination:
    
    if(rows)	free(rows);
    
    return code;
}







OPT_QCPSolver::OPT_QCPSolver():OPT_QPSolver()
{
}


OPT_QCPSolver::~OPT_QCPSolver()
{
}



OPT_SOLVERTYPE OPT_QCPSolver::getSolverType()
{
    return OPT_QCP;
}



int OPT_QCPSolver::copyQCPPartFrom( OPT_QCPSolver &other, const bool copyObj, const bool copyConstrs, const bool copyVarBounds, const bool copyVarTypes)
{
    int r, retCode;
    int *rows = NULL, *cols = NULL;
    double *vals  = NULL;
    
    r = copyQPPartFrom(other, copyObj, copyConstrs, copyVarBounds, copyVarTypes);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    if( copyConstrs )
    {
        int m;
        int nzs = 0;
        
        r = other.getNumberOfConstraints(m);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        for(int i = 0; i < m; i++)
        {
            int inzs, inzs2;
            
            r = other.getNumberOfConstraintQuadTerms(i, inzs);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            if( inzs <= 0 )
                continue;
            
            if( inzs > nzs )
            {
                OPT_secFree(rows);
                OPT_secFree(cols);
                OPT_secFree(vals);
                
                OPT_malloc( rows, inzs );
                OPT_malloc( cols, inzs );
                OPT_malloc( vals, inzs );
                
                OPT_IFMEMERRORGOTOLABEL( !rows || !cols || !vals, retCode, termination );
                
                nzs = inzs;
            }
            
            r = other.getConstraintQuadMatrix(i, inzs2, rows, cols, vals);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            #if MRQ_DEBUG_MODE
                assert(inzs == inzs2);
            #endif
            
            r = setConstraintQuadMatrix(i, inzs, rows, cols, vals);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        }
    }
    
    retCode = 0;
    
termination:
    
    if(rows)	free(rows);
    if(cols)	free(cols);
    if(vals)	free(vals);
    
    return retCode;
}



int OPT_QCPSolver::setProblemFrom( const minlpproblem:: MIP_MINLPProb &prob, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int naddvars, const int naddconstrs )
{
    const int m = prob.getNumberOfConstraints();
    int r, code;
    int maxnzs = 0;
    int *rows = NULL;//, *cols;
    //double *values = NULL;
    
    code = OPT_QPSolver::setProblemFrom(prob, setObj, setConstrs, setVarBounds, setVarType, naddvars, naddconstrs);
    if( code != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(code);
        #endif
        
        goto termination;
    }
    
    
    if( setConstrs )
    {
        for(int i = 0; i < m; i++)
        {
            int nzs;
            
            r = prob.getNumberOfConstraintQuadCoefMatrixTerms(i, nzs);
            
            #if OPT_DEBUG_MODE
                if( r != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(r);
                    #endif
                        
                    return OPT_SOLVER_ERROR;
                }
            #endif
            
            if( nzs > maxnzs )
                maxnzs = nzs;
        }
        
        if( maxnzs > 0) 
        {
            rows = (int *) malloc( maxnzs * sizeof(int) );
            //values = (double *) malloc( maxnzs * sizeof(double) );
            
            if( !rows )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                code = OPT_MEMORY_ERROR;
                goto termination;
            }
            
            
            for(int i = 0; i < m; i++)
            {
                int nzs;
                
                r = prob.getNumberOfConstraintQuadCoefMatrixTerms(i, nzs);
                #if OPT_DEBUG_MODE
                    if( r != 0 )
                    {
                        #if OPT_DEBUG_MODE
                            OPT_PRINTERRORNUMBER(r);
                        #endif
                            
                        return OPT_SOLVER_ERROR;
                    }
                #endif
                
                
                if( nzs > 0 )
                {
                    const MIP_SparseMatrix &QCi = prob.QC[i];
                    
                    QCi.getStructure(rows, NULL);
                    const int* const cols = QCi[0];
                    const double* const values = QCi(0);
                    
                    /*r = prob.getConstraintQuadCoefMatrix( i, &nzs, rows, cols, values );
                    
                    #if OPT_DEBUG_MODE
                        if( r != 0 )
                        {
                            #if OPT_DEBUG_MODE
                                OPT_PRINTERRORNUMBER(r);
                            #endif
                                
                            return OPT_SOLVER_ERROR;
                        }
                    #endif */
                    
                    
                    code = setConstraintQuadMatrix( i, nzs, rows, cols, values );
                    if( code != 0 )
                    {
                        #if OPT_DEBUG_MODE
                            OPT_PRINTERRORNUMBER(code);
                        #endif
                        goto termination;
                    }
                }
            }
            
            
        }
        
    }
    
    
    code = 0;
    
termination:
    
    if( code != 0 )
    {
        deallocateMemory();
        deallocateSolverEnv();
    }
    
    
    if( rows )		free(rows);
    //if( values )	free(values);
    
    return code;
}


int OPT_QCPSolver::setConstraintQuadMatrix( const int index, const int *qrowStart, const int *qcols, const double *qvalues )
{
    int n, k, code;
    getNumberOfVars(n);
    
    const int nzs = qrowStart[n];
    int *rows = NULL;
    
    rows = (int *) malloc( nzs * sizeof(int) );
    if( !rows )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    k = 0;
    for(int i = 0; i < n; i++)
    {
        const int rnzs = qrowStart[i+1] - qrowStart[i];
        
        OPT_setAllArray(rnzs, &rows[k], i);
        k += rnzs;
    }
    
    #if OPT_DEBUG_MODE
        assert(k == nzs);
    #endif
    
    code = setConstraintQuadMatrix(index, nzs, rows, qcols, qvalues);
    if( code != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(code);
        #endif
        goto termination;
    }
    
    
    code = 0;
    
termination:
        
    if(rows)	free(rows);
    
    return 0;
}


OPT_NLPSolver::OPT_NLPSolver():OPT_QCPSolver()
{
    initialize();
}



OPT_NLPSolver::~OPT_NLPSolver()
{
    deallocateMemory();
}


int OPT_NLPSolver::allocateConstrStructures(const int m)
{
    //bool *auxb;
    
    
    
    if( m > maux )
    {
        /*auxb = (bool *) realloc( auxCEval, m*sizeof(bool) );
        if( !auxb )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            return OPT_MEMORY_ERROR;
        }
        
        auxCEval = auxb;
        
        for(int i = maux; i < m; i++)
            auxCEval[i] = false; */
        
        int r = OPT_realloc( auxCEval, m );
        OPT_IFERRORRETURN(r, r);
        
        if( m > maux)
            OPT_setAllArray(m-maux, &auxCEval[maux], false);
    }
    
    if( useSPMtoJacAndHess )
    {
        const int mJ = J.getNumberOfRows();
        
        if( m > mJ )
        {
            const int r = J.addNewRows(m - mJ);
            if( r != 0  )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                return OPT_MEMORY_ERROR;
            }
        }
        else
        {
            for(int i = m; i < mJ; i++)
                J.deleteRowStructure(i);
        }
    }
    
    //we are calling allocateConstrStructures at end because we use variable maux 
    return OPT_QCPSolver::allocateConstrStructures(m);
}


int OPT_NLPSolver::allocateVarStructures(const int n)
{
    if( useSPMtoJacAndHess )
    {
        const int nH = lagH.getNumberOfRows();
        
        if( n > nH )
        {
            const int r = lagH.addNewRows( n - nH );
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                return OPT_MEMORY_ERROR;
            }
        }
        else
        {
            for(int i = n; i < nH; i++)
                lagH.deleteRowStructure(i);
        }
        
        
        J.setNumberOfColumns(n);
        lagH.setNumberOfColumns(n);
    }
    
    //we are calling allocateConstrStructures at end because we use variable maux 
    
    const int r = OPT_QCPSolver::allocateVarStructures(n);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
            OPT_PRINTERRORNUMBER(r);
    #endif
    
    
    return r;
}


/*int OPT_NLPSolver::allocateSPMatrices()
{
    J = new (nothrow) OPT_SparseMatrix(0, 0, false);
    H = new (nothrow) OPT_SparseMatrix(0, 0, true);
    
    
    if( !J || !H )
        return OPT_MEMORY_ERROR;
    
    return 0;
} */



int OPT_NLPSolver::copyNLPartFrom( OPT_NLPSolver &other, const bool copyObj, const bool copyConstrs, const bool copyVarBounds, const bool copyVarTypes, const bool shareEvaluationObject)
{
    int r, retCode;
    
    int nzs = 0;
    int *rows = NULL, *cols = NULL;
    
    
    r = copyQCPPartFrom(other, copyObj, copyConstrs, copyVarBounds, copyVarTypes);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    if( copyConstrs )
    {
        int m;
        
        r = other.getNumberOfConstraints(m);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        for(int i = 0; i < m; i++)
        {
            bool nlConstr;
            
            r = other.getConstrNLFlag(i, nlConstr);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            if( nlConstr )
            {
                r = setConstrNLFlag(i, nlConstr);
                OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            }
        }
        
        //copying jacobian structure
        if( useSPMtoJacAndHess || isMyNLPClass() )
        {
            int *rowStart;
            int *jacCols;
            OPT_SparseMatrix *pJ; 
            
            if( useSPMtoJacAndHess )
            {
                pJ = &other.J;
            }
            else
            {
                #if OPT_DEBUG_MODE
                    assert( isMyNLPClass() ); //if we do not pass in this assert, we need at least else if here to threat a new case
                #endif
                
                pJ = &( (OPT_MyNLPSolver*) &other)->prob.J;
            }
            
            rowStart = (int *) pJ->offset;
            jacCols = pJ->baseStructure.cols;
            
            
            r = setJacobianStructure(rowStart, jacCols);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        }
        else
        { //so, we use our methods to set...
            int nzJac;
            
            r = other.getNumberOfNonZerosInJacobian(nzJac);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            if( nzJac > 0 )
            {
                #if OPT_DEBUG_MODE
                    assert(nzs == 0);
                #endif
                
                OPT_malloc(rows, nzJac);
                OPT_malloc(cols, nzJac);
                
                OPT_IFMEMERRORGOTOLABEL(!rows || !cols , retCode, termination );
                
                
                nzs = nzJac;
            
                r = other.getJacobianStructure(nzJac, rows, cols);
                OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
                
                #if OPT_DEBUG_MODE
                    assert(nzs == nzJac);
                #endif
                
                r = setJacobianStructure(nzJac, rows, cols);
                OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            }
        }
        
    }
    
    
    if( copyObj )
    {
        bool objNLFlag;
        
        r = other.getObjNLFlag(objNLFlag);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        r = setObjNLFlag(objNLFlag);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    }
    
    if( copyObj || copyConstrs )
    { //so, we must set the lagranian hessian structure
        
        if( useSPMtoJacAndHess || isMyNLPClass() )
        {
            int *rowStart;
            int *hessCols;
            OPT_SparseMatrix *pH;
            
            if( useSPMtoJacAndHess )
            {
                pH = &other.lagH;
            }
            else
            {
                #if OPT_DEBUG_MODE
                    assert( isMyNLPClass() );
                #endif
                
                pH = &( (OPT_MyNLPSolver*) &other)->prob.lagH;
            }
            
            rowStart = (int *) pH->offset;
            hessCols = pH->baseStructure.cols;
            
            r = setLagrangianHessianStructure(rowStart, hessCols);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        }
        else
        {
            int nzHess;
            
            r = other.getNumberOfNonZerosInLagrangianHessian(nzHess);
            OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            
            if( nzHess > 0 )
            {
                if(nzHess > nzs)
                {
                    free(rows);
                    free(cols);
                    
                    OPT_malloc(rows, nzs);
                    OPT_malloc(cols, nzs);
                    
                    OPT_IFMEMERRORGOTOLABEL(!rows || !cols , retCode, termination );
                }
                
                r = other.getLagrangianHessianStructure(nzs, rows, cols);
                OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
                
                r = setLagrangianHessianStructure(nzs, rows, cols);
                OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
            }
        }
    }
    
    
    if( shareEvaluationObject )
    {
        OPT_NonLinearEval *nlEval;
        
        r = other.getNonLinearEvalObjectPointer( nlEval);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        r = setNonLinearEvalObject(nlEval);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    }
    
    
    retCode = 0;
    
termination:
    
    if(rows)	free(rows);
    if(cols)	free(cols);
    
    return retCode;
}



void OPT_NLPSolver::deallocateMemory()
{
    OPT_secFree(auxCEval);
    
    J.desallocateMemory();
    lagH.desallocateMemory();

    OPT_QCPSolver::deallocateMemory();
}



void OPT_NLPSolver::initialize()
{
    OPT_QCPSolver::initialize();
    
    in_nl_obj_factor = 1.0;
    
    in_absolute_feas_tol = 1e-6;
    in_relative_feas_tol = 1e-6;
    
    threadNumber = 0;
    auxCEval = NULL;
    nlEval = NULL;
    
    useSPMtoJacAndHess = true;
    
    lagH.setSymmetricFlag(true);
}


OPT_SOLVERTYPE OPT_NLPSolver::getSolverType()
{
    return OPT_NLP;
}



int OPT_NLPSolver::getLagrangianHessianStructure(int &nzs, int* rows, int* cols)
{
    if(!useSPMtoJacAndHess)
        return OPT_OPERATION_NOT_IMPLEMENTED; //it must be implemented in derived classes
    
    nzs = lagH.getStructure(rows, cols);
    return 0;
}


int OPT_NLPSolver::getJacobianStructure(int &nzs, int* rows, int* cols)
{
    if(!useSPMtoJacAndHess)
        return OPT_OPERATION_NOT_IMPLEMENTED; //it must be implemented in derived classes
    
    nzs = J.getStructure(rows, cols);
    return 0;
}


int OPT_NLPSolver::getNonLinearEvalObjectPointer( OPT_NonLinearEval* &nlEval )
{
    nlEval = this->nlEval;
    return 0;
}


int OPT_NLPSolver::getNumberOfNonZerosInLagrangianHessian(int &nzs)
{
    if(!useSPMtoJacAndHess)
        return OPT_OPERATION_NOT_IMPLEMENTED; //it must be implemented in derived classes
    
    nzs = lagH.getNumberOfElements();
    return 0;
}


int OPT_NLPSolver::getNumberOfNonZerosInJacobian(int &nzs)
{
    if(!useSPMtoJacAndHess)
        return OPT_OPERATION_NOT_IMPLEMENTED; //it must be implemented in derived classes
    
    nzs = J.getNumberOfElements();
    return 0;
}


int OPT_NLPSolver::removeConstraints(const int ninds, const int* indices )
{
    J.removeRows(ninds, indices);
    
    
    return __removeConstraints(ninds, indices);
}



int OPT_NLPSolver::setJacobianStructure(const int nzs, const int* rows, const int* cols)
{
    J.deleteStructure();
    
    const int r = J.setStructureAndValues( nzs, rows, cols );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return r == SPM_MEMORY_ERROR ? OPT_MEMORY_ERROR : OPT_UNDEFINED_ERROR;
    }
    
    return 0;
}


int OPT_NLPSolver::setJacobianStructure(const int* rowStart, const int* cols)
{
    J.deleteStructure();
    
    const int r = J.setStructureAndValues(rowStart, cols);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return r == SPM_MEMORY_ERROR ? OPT_MEMORY_ERROR : OPT_UNDEFINED_ERROR;
    }
    
    return 0;
}


//warning: that method store the jacobian of a  nonlinear constraint. If that line have already been defined, it will be overwritten
int OPT_NLPSolver::setJacobianRowStructure(const int row, const int nzs, const int* cols)
{
    int r = J.setRowStructure( row, nzs, cols );
    
    if( r != 0 )
    {
        //std::cout << "J: \n";
        //std::cout << "row: " << row << " nzs: " << nzs << "\n";
        
        
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
                
        return r == SPM_MEMORY_ERROR ? OPT_MEMORY_ERROR : OPT_UNDEFINED_ERROR;
    }
    
    return 0;
}



int OPT_NLPSolver::setLagrangianHessianStructure(const int nzs, const int* rows, const int* cols)
{
    lagH.deleteStructure();
    
    const int r = lagH.setStructureAndValues(nzs, rows, cols);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
                
        return r == SPM_MEMORY_ERROR ? OPT_MEMORY_ERROR : OPT_UNDEFINED_ERROR;
    }
    
    return 0;
}



int OPT_NLPSolver::setLagrangianHessianStructure(const int *rowStart, const int* cols)
{
    lagH.deleteStructure();
    
    const int r= lagH.setStructureAndValues(rowStart, cols);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
                
        return r == SPM_MEMORY_ERROR ? OPT_MEMORY_ERROR : OPT_UNDEFINED_ERROR;
    }
    
    return 0;
}



int OPT_NLPSolver::setLagrangianHessianRowStructure( const int row, const int nzs, const int* cols)
{
    int r = lagH.setRowStructure(row, nzs, cols);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
                
        return r == SPM_MEMORY_ERROR ? OPT_MEMORY_ERROR : OPT_UNDEFINED_ERROR;
    }
    
    return 0;
}



int OPT_NLPSolver::setNonLinearEvalObject( OPT_NonLinearEval* nlEval )
{
    this->nlEval = nlEval;
    return 0;
}



int OPT_NLPSolver::setProblemFrom( const MIP_MINLPProb& prob, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int naddvars, const int naddconstrs )
{
    const int m = prob.getNumberOfConstraints();
    const int n = prob.getNumberOfVars();
    int code;
    
    
    code = OPT_QCPSolver::setProblemFrom(prob, setObj, setConstrs, setVarBounds, setVarType, naddvars, naddconstrs );
    if( code != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(code);
        #endif
        goto termination;
    }
    
    
    if( setObj && prob.hasObjNLTerm() )
    {
        setObjNLFlag(true);
        in_nl_obj_factor = prob.getObjFactor(); //we already considered the objfactor to set Q and c
    }
    
    
    if( setConstrs )
    {
        if( prob.hasNLConstraints() )
        {
            for(int i = 0; i < m; i++)
            {
                if( prob.hasConstraintNLTerm(i) )
                    setConstrNLFlag(i, true);
            }
        }
        
        if( prob.getNumberOfJacobianNonZeros() > 0 )
        {
            int mym;
            const int r = getNumberOfConstraints(mym);
            
            #if OPT_DEBUG_MODE
                assert(r == 0 && mym >= m);
            #endif
            
            
            /*for(int i = 0; i < m; i++)
            {
                int nzs;
                
                r = prob.getJacobianRowStructure( i, nzs, auxIndex  );
                
                #if OPT_DEBUG_MODE
                    if( r != 0 )
                    {
                        #if OPT_DEBUG_MODE
                            OPT_PRINTERRORNUMBER(r);
                        #endif
                            
                        return OPT_SOLVER_ERROR;
                    }
                #endif
                
                code = setJacobianRowStructure(i, nzs, auxIndex);
                if( code != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(code);
                    #endif
                    
                    goto termination;
                }
            } */
            
            //TODO: resolve it someday: jac_rowStart should be an unsigned int because prob.lagH.offset is unsigned. However, setLagrangianHessianStructure rceieve an int array... 
            const int *jac_rowStart = (int*) prob.J.offset;
            const int *jac_cols = prob.J.getRowColsPointer(0);
            
            if(mym > m)
            { //in this case, the OPT_solver object has more constraints than original MIP_MINLPProb. We cannot pass the prob.J.offset directly because it will generate a segmentation fault when spm try copy the extra rows
                
                OPT_copyArray(m, prob.J.offset, auxIndex);
                OPT_setAllArray<int>(mym-m+1, &auxIndex[m], prob.J.offset[m]);
                jac_rowStart = auxIndex;
            }
            
            
            
            code = setJacobianStructure(jac_rowStart, jac_cols);
            if( code != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(code);
                #endif
                
                goto termination;
            }
            
        }
    }
    
    
    if( setObj || setConstrs )
    {
        setNonLinearEvalObject( prob.getNonLinearEvaluationObject() );
        
        
        if( (setObj == false && prob.hasNLConstraints() == false) || (setConstrs == false && prob.hasObjNLTerm() == false) ) 
        {
            //we do not have nlp lagrangian hessian...
        }
        else
        {
            if( prob.getNumberOfLagrangianHessianNonZeros() > 0 )
            {
                int myn;
                const int r = getNumberOfVars(myn);
                
                #if OPT_DEBUG_MODE
                    assert(r == 0 && myn >= n);
                #endif
                
                
                /*for(int i = 0; i < n; i++ )
                {
                    int nzs;
                    
                    r = prob.getLagrangianHessianRowStructure(i, &nzs, auxIndex);
                    #if OPT_DEBUG_MODE
                        if( r != 0 )
                        {
                            #if OPT_DEBUG_MODE
                                OPT_PRINTERRORNUMBER(r);
                            #endif
                                
                            return OPT_SOLVER_ERROR;
                        }
                    #endif
                    
                    code = setLagrangianHessianRowStructure( i, nzs, auxIndex );
                    if( code != 0 )
                    {
                        #if OPT_DEBUG_MODE
                            OPT_PRINTERRORNUMBER(code);
                        #endif
                        
                        goto termination;
                    }
                } */
                
                
                //const unsigned int lagh_nzs = prob.lagH.getNumberOfElements();
                
                //TODO: resolve it someday: lagh_rowStart should be an unsigned int because prob.lagH.offset is unsigned. However, setLagrangianHessianStructure rceieve an int array... 
                const int *lagh_rowStart = (int*) prob.lagH.offset;
                const int* cols = prob.lagH.getRowColsPointer(0);
                
                if( myn > n )
                {
                    OPT_copyArray(n, prob.lagH.offset, auxIndex);
                    OPT_setAllArray<int>(myn-n+1, &auxIndex[n], prob.lagH.offset[n]);
                    lagh_rowStart = auxIndex;
                }
                
                
                code = setLagrangianHessianStructure( lagh_rowStart, cols );
                if( code != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(code);
                    #endif
                    
                    goto termination;
                }
            }
        }
    
    }
    
    
    
    code = 0;
    
termination:
    
    
    if( code != 0 )
    {
        deallocateSolverEnv();
        deallocateMemory();
    }
    
    return code;
}



void OPT_NLPSolver::setThreadNumber(const unsigned int threadNumber)
{
    this->threadNumber = threadNumber;
}





























OPT_LPSolver* optsolvers::OPT_newLPSolver( const int solver )
{
    //#define OPT_LP_SOLVER_CODE_COUNTER 0
    
    if( solver == OPT_LP_GLPK )
        //#define OPT_LP_SOLVER_CODE_COUNTER 1
        return new (nothrow) OPT_Glpk;
    else if( solver == OPT_LP_CBC )
        //#define OPT_LP_SOLVER_CODE_COUNTER 2
        return new (nothrow) OPT_Cbc;
    else if( solver == OPT_LP_CPLEX )
        //#define OPT_LP_SOLVER_CODE_COUNTER 3
        return new (nothrow) OPT_Cplex;
    else if( solver == OPT_LP_GUROBI )
        //#define OPT_LP_SOLVER_CODE_COUNTER 4
        return new (nothrow) OPT_Gurobi;
    else if( solver == OPT_LP_XPRESS )
        //#define OPT_LP_SOLVER_CODE_COUNTER 5
        return new (nothrow) OPT_Xpress;
    else if( solver == OPT_LP_MOSEK )
        //#define OPT_LP_SOLVER_CODE_COUNTER 6
        return new (nothrow) OPT_Mosek;
    else if( solver == OPT_LP_KNITRO )
        //#define OPT_LP_SOLVER_CODE_COUNTER 7
        return new (nothrow) OPT_Knitro;
    else if( solver == OPT_LP_IPOPT )
        //#define OPT_LP_SOLVER_CODE_COUNTER 8
        return new (nothrow) OPT_Ipopt;
    else if( solver == OPT_LP_ALGENCAN )
        //#define OPT_LP_SOLVER_CODE_COUNTER 9
        return new (nothrow) OPT_Algencan;
    else if( solver == OPT_LP_WORHP )
        //#define OPT_LP_SOLVER_CODE_COUNTER 10
        return new (nothrow) OPT_Worhp;
    else if( solver == OPT_LP_OPTIZELLE )
        //#define OPT_LP_SOLVER_CODE_COUNTER 11
        return nullptr;// new (nothrow) OPT_Optizelle;
    else if( solver == OPT_LP_IQUAD )
        #define OPT_LP_SOLVER_CODE_COUNTER 12
        return new (nothrow) OPT_Iquad;
    
    #if OPT_LP_SOLVER_CODE_COUNTER <  OPT_NUMBER_OF_SOLVERS
        Compile error. You forgot consider some solver here! (do not comeemt this)
    #endif
    else
    {
        #if OPT_DEBUG_MODE
            std::cerr << OPT_PREPRINT << "Invalid solver code: " << solver << OPT_GETFILELINE ;
        #endif
    }

    return NULL;
}


OPT_QPSolver* optsolvers::OPT_newQPSolver( const int solver )
{
    if( solver == OPT_QP_CPLEX )
        //#define OPT_QP_SOLVER_CODE_COUNTER 3
        return new (nothrow) OPT_Cplex;
    else if( solver == OPT_QP_GUROBI )
        //#define OPT_QP_SOLVER_CODE_COUNTER 4
        return new (nothrow) OPT_Gurobi;
    else if( solver == OPT_QP_XPRESS )
        //#define OPT_QP_SOLVER_CODE_COUNTER 5
        return new (nothrow) OPT_Xpress;
    else if( solver == OPT_QP_MOSEK )
        //#define OPT_QP_SOLVER_CODE_COUNTER 6
        return new (nothrow) OPT_Mosek;
    else if( solver == OPT_QP_KNITRO )
        //#define OPT_QP_SOLVER_CODE_COUNTER 7
        return new (nothrow) OPT_Knitro;
    else if( solver == OPT_QP_IPOPT )
        //#define OPT_QP_SOLVER_CODE_COUNTER 8
        return new (nothrow) OPT_Ipopt;
    else if( solver == OPT_QP_ALGENCAN )
        //#define OPT_QP_SOLVER_CODE_COUNTER 9
        return new (nothrow) OPT_Algencan;
    else if( solver == OPT_QP_WORHP )
        //#define OPT_QP_SOLVER_CODE_COUNTER 10
        return new (nothrow) OPT_Worhp;
    else if( solver == OPT_QP_OPTIZELLE )
        //#define OPT_QP_SOLVER_CODE_COUNTER 11
        return nullptr; //new (nothrow) Optizelle;
    else if( solver == OPT_QP_IQUAD )
        #define OPT_QP_SOLVER_CODE_COUNTER 12
        return new (nothrow) OPT_Iquad;
    #if OPT_QP_SOLVER_CODE_COUNTER < OPT_NUMBER_OF_SOLVERS
        Compile error. You forgot consider some solver here! (do not comeemt this)
    #endif
    else
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORMSGP("Invalid solver code: ", solver);
        #endif
    }

    return NULL;
}


OPT_QCPSolver* optsolvers::OPT_newQCPSolver( const int solver )
{
    if( solver == OPT_QCP_XPRESS )
        //#define OPT_QCP_SOLVER_CODE_COUNTER 5
        return new (nothrow) OPT_Xpress;
    else if( solver == OPT_QCP_MOSEK )
        //#define OPT_QCP_SOLVER_CODE_COUNTER 6
        return new (nothrow) OPT_Mosek;
    else if( solver == OPT_QCP_KNITRO )
        //#define OPT_QCP_SOLVER_CODE_COUNTER 7
        return new (nothrow) OPT_Knitro;
    else if( solver == OPT_QCP_IPOPT )
        //#define OPT_QCP_SOLVER_CODE_COUNTER 8
        return new (nothrow) OPT_Ipopt;
    else if( solver == OPT_QCP_ALGENCAN )
        //#define OPT_QCP_SOLVER_CODE_COUNTER 9
        return new (nothrow) OPT_Algencan;
    else if( solver == OPT_QCP_WORHP )
        //#define OPT_QCP_SOLVER_CODE_COUNTER 10
        return new (nothrow) OPT_Worhp;
    else if( solver == OPT_QCP_OPTIZELLE )
        //#define OPT_QCP_SOLVER_CODE_COUNTER 11
        return nullptr; //new (nothrow) OPT_Optizelle;
    else if( solver == OPT_QCP_IQUAD )
        #define OPT_QCP_SOLVER_CODE_COUNTER 12
        return new (nothrow) OPT_Iquad;
    #if OPT_QCP_SOLVER_CODE_COUNTER < OPT_NUMBER_OF_SOLVERS
        Compile error. You forgot consider some solver here! (do not comeemt this)
    #endif
    else
    {
        #if OPT_DEBUG_MODE
            std::cerr << OPT_PREPRINT << "Invalid solver code: " << solver << OPT_GETFILELINE ;
        #endif
    }
    
    return NULL;
}


OPT_NLPSolver* optsolvers::OPT_newNLPSolver( const int solver )
{
    if( solver == OPT_NLP_MOSEK )
        //#define OPT_NLP_SOLVER_CODE_COUNTER 6
        return new (nothrow) OPT_Mosek;
    else if( solver == OPT_NLP_KNITRO )
        //#define OPT_NLP_SOLVER_CODE_COUNTER 7
        return new (nothrow) OPT_Knitro;
    else if( solver == OPT_NLP_IPOPT )
        //#define OPT_NLP_SOLVER_CODE_COUNTER 8
        return new (nothrow) OPT_Ipopt;
    else if( solver == OPT_NLP_ALGENCAN )
        //#define OPT_NLP_SOLVER_CODE_COUNTER 9
        return new (nothrow) OPT_Algencan;
    else if( solver == OPT_NLP_WORHP )
        //#define OPT_NLP_SOLVER_CODE_COUNTER 10
        return new (nothrow) OPT_Worhp;
    else if( solver == OPT_NLP_OPTIZELLE )
        //#define OPT_NLP_SOLVER_CODE_COUNTER 11
        return nullptr; //new (nothrow) OPT_Optizelle;
    else if( solver == OPT_NLP_IQUAD )
        #define OPT_NLP_SOLVER_CODE_COUNTER 12
        return new (nothrow) OPT_Iquad;
    #if OPT_NLP_SOLVER_CODE_COUNTER < OPT_NUMBER_OF_SOLVERS
        Compile error. You forgot consider some solver here! (do not comeemt this)
    #endif
    else
    {
        #if OPT_DEBUG_MODE
            std::cerr << OPT_PREPRINT << "Invalid solver code: " << solver << OPT_GETFILELINE ;
        #endif
    }
    
    return NULL;
}




int optsolvers::OPT_setMINLPProblem(const MIP_MINLPProb& prob, OPT_Solver* solver, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int naddvars, const int naddconstrs )
{
    int code; 
    const int ptype = prob.getProblemType();
    
    
    if( MIP_isLinearProblemType(ptype) ) 
    {
        code = ( (OPT_LPSolver *) solver )->setProblemFrom(prob, setObj, setConstrs, setVarBounds, setVarType, naddvars, naddconstrs); //problem is linear. There is no difference the value of parameter set setStrictPartOnly
    }
    else if( MIP_isQuadOnlyProblemType(ptype) )
    {
        code = ( (OPT_QPSolver *) solver )->setProblemFrom(prob, setObj, setConstrs, setVarBounds, setVarType, naddvars, naddconstrs);
    }
    else if( MIP_isQuadConstrProblemType(ptype) )
    {
        code = ( (OPT_QCPSolver *) solver )->setProblemFrom(prob, setObj, setConstrs, setVarBounds, setVarType, naddvars, naddconstrs);
    }
    else //if(MIP_isNonlinearProblemType(ptype))
    {
        code = ( (OPT_NLPSolver *) solver )->setProblemFrom(prob, setObj, setConstrs, setVarBounds, setVarType, naddvars, naddconstrs);
    }
    
    return code;
}



int optsolvers::OPT_setQCPProblemOnCplex(const minlpproblem::MIP_MINLPProb& prob, OPT_Cplex *optCplex, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int naddvars, const int naddconstrs )
#if OPT_HAVE_CPLEX
{
    int r;
    
    CPXENVptr  &cpxEnv = optCplex->env;
    CPXLPptr   &cpxProb = optCplex->prob;
    
    
    const int m = prob.getNumberOfConstraints();
    int code, nzq = 0;
    int *qck_rows = NULL; //*cols = NULL;
    double *qkc_values = NULL;
    const minlpproblem::MIP_SparseMatrix *QC = prob.QC;
    const minlpproblem::MIP_SparseMatrix &A = prob.A;
    
    
    
    
    r = optCplex->setProblemFrom( prob, setObj, false, setVarBounds, setVarType, naddvars, naddconstrs );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return r;
    }
    
    if( !setConstrs )
        return 0;
    
    
    
    
    
    
    
    for( int i = 0; i < m; i++)
    {
        const int nzQi = QC[i].getNumberOfElements();
        
        if( nzQi > nzq )
            nzq = nzQi;
        
        /*int nzar;
        prob.getNumberOfLinearCoefsInConstr( i, nzar );
        if( nzar > nza )
            nza = nzar;*/
    }
    
    
    if(nzq > 0)
    {
        qck_rows = (int *) malloc( nzq *sizeof(int) );
        //cols = (int *) malloc( (nzq + nza) * sizeof(int) );
        qkc_values = (double *) malloc( nzq * sizeof(double) );
        
        if( !qck_rows || !qkc_values )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            
            code = OPT_MEMORY_ERROR;
            goto termination;
        }
    }
    
    
    
    for(int k = 0; k < m; k++)
    {
        double lb, ub;
        prob.getConstraintBounds(k, lb, ub);
        
        
        //const int nzqk = QC[k].getStructureAndValues( rows, cols, vals );
        
        //we take advantage internal pointers to getcols
        const int* const qck_cols = QC[k].getRowColsPointer(0);
        const int nzqk = QC[k].getStructure(qck_rows, NULL); 
        const int nzqk2 = QC[k].getValues(qkc_values);
        
        #if OPT_DEBUG_MODE
            assert( nzqk == nzqk2 );
        #endif
        
        
        for( int i = 0; i < nzqk; i++ )
        {
            if( qck_rows[i] == qck_cols[i] )
                qkc_values[i] = qkc_values[i]*0.5;
        }
        
        //prob.getConstraintLinearPart( k, nzak, colsa, valsa );
        
        const int* const colsa = A.getRowColsPointer(k); //&(cols[nzqk]);
        const double* const valsa = A.getRowValuesPointer(k); //&(vals[nzqk]);
        const int nzak = A.getNumberOfElementsAtRow(k);
        
        
        
        if( nzqk > 0 )
        {
            #if 0
            //unfortunatelly, this code does not work... I hate cplex
            
            const char sense = optCplex->constraintSense(lb, ub);
            const double rhs = sense == 'L' ? ub : lb;
                
            
            r = CPXaddqconstr( cpxEnv, cpxProb, nzak, nzqk, rhs, sense, colsa, valsa, rows, cols, vals, NULL );
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                code = OPT_SOLVER_ERROR;
                goto termination;
            }
            
            
            if( sense == 'R' )
            {
                double rngval = ub -lb;// second explanations in CPXnewrows (cplex manual), ranged constraints are between rhs and rhs + rngval.
                
                puts("droga 1");
                
                r = CPXchgrngval( cpxEnv, cpxProb, 1, &k, &rngval );
                
                puts("droga 2");
                
                if( r != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(r);
                    #endif
                        
                    code = OPT_SOLVER_ERROR;
                    goto termination;
                }
            }
            
            #endif
        
        
        
            if( lb > -OPT_INFINITY )
            {
                r = CPXaddqconstr( cpxEnv, cpxProb, nzak, nzqk, lb, 'G', colsa, valsa, qck_rows, qck_cols, qkc_values, NULL );
            
                if( r != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(r);
                    #endif
                    
                    code = OPT_SOLVER_ERROR;
                    goto termination;
                }
            }
            
            
            //unfortunatelly, we have to add two constaints if the constraint is double bounded...
            
            if( ub < OPT_INFINITY )
            {
                r = CPXaddqconstr( cpxEnv, cpxProb, nzak, nzqk, ub, 'L', colsa, valsa, qck_rows, qck_cols, qkc_values, NULL );
            
                if( r != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(r);
                    #endif
                    
                    code = OPT_SOLVER_ERROR;
                    goto termination;
                }
            }
        }
        else 
        {
            const char sense = optCplex->constraintSense(lb, ub);
            const double rhs = sense == 'L' ? ub : lb;
            
            qck_rows[0] = 0;
            
            
            r = CPXaddrows( cpxEnv, cpxProb, 0, 1, nzak, &rhs, &sense, qck_rows, colsa, valsa, NULL, NULL );
            
            
            if( sense == 'R' )
            {
                double rngval = ub -lb;// second explanations in CPXnewrows (cplex manual), ranged constraints are between rhs and rhs + rngval.
                
                const int index = CPXgetnumrows( cpxEnv, cpxProb )-1;
                
                r = CPXchgrngval( cpxEnv, cpxProb, 1, &index, &rngval );
                
                if( r != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(r);
                    #endif
                        
                    code = OPT_SOLVER_ERROR;
                    goto termination;
                }
            }
            
        }
        
    }
    
    
    
    code = 0;
    
termination:
    
    if(qck_rows)		free(qck_rows);
    //if(cols)	free(cols);
    if(qkc_values)	free(qkc_values);
    
    return code;
}
#else
{
    #if OPT_DEBUG_MODE
        OPT_PRINTERRORNUMBER(OPT_LIBRARY_NOT_AVAILABLE);
    #endif
    return OPT_LIBRARY_NOT_AVAILABLE;
}
#endif






int OPT_ObjCutSetter::setObjCut( OPT_LPSolver *solver, const int constrIndex, const minlpproblem::MIP_MINLPProb &prob, const double zu, minlpproblem::MIP_NonLinearEval *eval ) const
{
    const bool linCoefObj = prob.hasLinCoefObj();
    const bool objNLTerm = prob.hasObjNLTerm();
    const int n = prob.n;
    const int nzq = prob.getNumberOfObjQuadTerms();
    const double of = prob.objFactor;
    
    int code, r, nz;
    int *cols = NULL, *cols2;
    double *vals = NULL, *p;
    
    
    
    nz = OPT_max( nzq, linCoefObj || objNLTerm ? n : 0 );
    
    if(nz > 0)
    {
        cols = (int *) malloc( 2*nz*sizeof(int) );
        vals = (double *) malloc( nz*sizeof(double) );
        
        if( !cols || !vals )
        {
            OPT_PRINTMEMERROR;
            code = OPT_MEMORY_ERROR;
            goto termination;
        }
        
        cols2 = &cols[nz];
    }
    
    if( linCoefObj )
    {
        int onz = 0;
        
        if( of == 1.0 )
        {
            p = prob.c;
        }
        else
        {
            p = vals;
            for(int i = 0; i < n; i++)
                p[i] = prob.c[i] * of;
        }
        
        
        for( int i = 0; i < n; i++ )
        {
            if( p[i] != 0.0 )
            {
                cols[onz] = i;
                vals[onz] = p[i]; //vals can be the same array, but there is no problem...
                onz++;
            }
        }
        
        
        //r = solver->setLinearObjFunction(n, p);
        r = solver->setConstraintLinearCoefs( constrIndex, onz, cols, vals  );
        if( r != 0 )
        {
            OPT_PRINTERRORNUMBER(r);
            code = OPT_SOLVER_ERROR;
            goto termination;
        }
    }
    
    
    if( nzq > 0 )
    {
        int nzquad;
        
        prob.getObjQuadCoefsMatrix(&nzquad, cols, cols2, vals);
        /*#if OPT_DEBUG_MODE
            if( r != 0 )
            {
                OPT_PRINTERRORNUMBER(r);
                code = OPT_SOLVER_ERROR;
                goto termination;
            }
        #endif */
        
        if( of != 1.0 )
        {
            for(int i = 0; i < nzq; i++)
                vals[i] *= of;
        }
        
        r = ( (OPT_QCPSolver *) solver)->setConstraintQuadMatrix( constrIndex, nzq, cols, cols2, vals );
        
        if( r != 0 )
        {
            OPT_PRINTERRORNUMBER(r);
            code = OPT_SOLVER_ERROR;
            goto termination;
        }
    }
    
    
    if( objNLTerm )
    {
        ( (OPT_NLPSolver *) solver)->setNonLinearEvalObject( eval );
        
        for(int i = 0; i < n; i++)
            cols[i] = i;
        
        r = ( (OPT_NLPSolver *) solver)->setJacobianRowStructure( constrIndex, n, cols );
        if( r != 0 )
        {
            OPT_PRINTERRORNUMBER(r);
            code = OPT_SOLVER_ERROR;
            goto termination;
        }
    }
    
    
    r = updateObjCut(solver, constrIndex, prob, zu);
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            OPT_PRINTERRORNUMBER(r);
            code = OPT_SOLVER_ERROR;
            goto termination;
        }
    #endif
    
    
    code = 0;
    
termination:
    
    if(cols)	free(cols);
    if(vals)	free(vals);
    
    return code;
}




int OPT_ObjCutSetter::updateObjCut( OPT_LPSolver* solver, const int constrIndex, const MIP_MINLPProb& prob, const double zu  ) const
{
    return solver->setConstraintBounds( constrIndex, -OPT_INFINITY, zu - prob.d* prob.objFactor );
    
}



int OPT_ConstrsBoundsCalculator::setProblemBase( MIP_MINLPProb& prob, const double* olc, const double* ouc, OPT_LPSolver* solver, OPT_GeneralSolverParams* params, const bool calcEvenOnOriginalValues, const bool* calcConstr) 
{
    const int m = prob.m;
    int r, nAuxVars = 0;
    
        
    for(int i = 0; i < m; i++)
    {
        if( calcConstr )
            if( !calcConstr[i] )
                continue;
                
        if( (calcEvenOnOriginalValues || olc[i] <= -MIP_INFINITY || ouc[i] >= MIP_INFINITY )  )
            nAuxVars++;
    }
    
    
    r = OPT_setMINLPProblem(prob, solver, false, true, true, false, nAuxVars);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return r;
    }
    
    if( params )
        solver->setParameters(*params);
    
    r = solver->setObjSense(optsolvers::OPT_MINIMIZE);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return r;
    }
    
    
    
    return 0;
}



OPT_NLPNonObjEval::OPT_NLPNonObjEval( const int noriginal, const int moriginal, const int nzjacoriginal, const int nzhessoriginal, minlpproblem::MIP_NonLinearEval *originalEval )
{
    initialize(noriginal, moriginal, nzjacoriginal, nzhessoriginal, originalEval);
}


OPT_NLPNonObjEval::~OPT_NLPNonObjEval()
{
}


void OPT_NLPNonObjEval::initialize(const int noriginal, const int moriginal, const int nzjacoriginal, const int nzhessoriginal, minlpproblem::MIP_NonLinearEval *originalEval)
{
    mynewx = true;
    norig = noriginal;
    morig = moriginal;
    nzjacorig = nzjacoriginal;
    nzhessorig = nzhessoriginal;
    
    oeval = originalEval;
}


int OPT_NLPNonObjEval::eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value)
{
    if(newx)
        mynewx = true;
    
    value = 0.0;
    return 0;
}


int OPT_NLPNonObjEval::eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *constrEval, const double *x, double *values)
{
    int r;
    
    if(newx)
        mynewx = true;
    
    r = oeval->eval_nl_constrs_part( threadnumber, norig, morig, mynewx, constrEval, x, values);
    
    #if OPT_DEBUG_MODE
        if(r != 0)
            OPT_PRINTCALLBACKERRORNUMBER(r);
    #endif
    
    mynewx = false;
    
    return r;
}


int OPT_NLPNonObjEval::eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *values)
{
    if(newx)
        mynewx = true;
    
    OPT_setAllArray(n, values, 0.0);
    
    return 0;
}


int OPT_NLPNonObjEval::eval_grad_nl_constrs_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *constrEval, const double *x, minlpproblem::MIP_SparseMatrix& jacobian)
{
    int r;
    
    if(newx)
        mynewx = true;
    
    r = oeval->eval_grad_nl_constrs_part( threadnumber, norig, morig, nzjacorig, mynewx, constrEval, x, jacobian );
    
    #if OPT_DEBUG_MODE
        if(r != 0)
            OPT_PRINTCALLBACKERRORNUMBER(r);
    #endif
    
    
    mynewx =false;
    
    return r;
}


int OPT_NLPNonObjEval::eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double objFactor, const double *lambda, minlpproblem::MIP_SparseMatrix& hessian)
{
    int r;
    
    if(newx)
        mynewx = true;
    
    r = oeval->eval_hessian_nl_lagran_part( threadnumber, norig, morig, nzhessorig, mynewx, x, 0.0, lambda, hessian );
    
    #if OPT_DEBUG_MODE
        if(r != 0)
            OPT_PRINTCALLBACKERRORNUMBER(r);
    #endif
    
    mynewx = false;
    
    return r;
}




int OPT_ConstrsBoundsCalculator::calculate( minlpproblem::MIP_MINLPProb &prob, const double *olc, const double *ouc, OPT_LPSolver *solver, OPT_GeneralSolverParams *params, const bool calcEvenOnOriginalValues, const bool *calcConstr, double *lc, double *uc )
{
    const int n = prob.n;
    const int m = prob.m;
    
    int r, auxInd, retCode;
    
    OPT_NLPNonObjEval eval(n, m, prob.getNumberOfJacobianNonZeros(), prob.getNumberOfLagrangianHessianNonZeros(), prob.getNonLinearEvaluationObject());
    
    
    if(!olc)
        olc = prob.lc;
    
    if(!ouc)
        ouc = prob.uc;
    
    r = setProblemBase(prob, olc, ouc, solver, params, calcEvenOnOriginalValues, calcConstr);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        retCode = r;
        goto termination;
    }
    
    if( prob.hasNlConstrs )
        ( (OPT_NLPSolver *) solver)->setNonLinearEvalObject(&eval);
    
    auxInd = n;
    for(int i = 0; i < m; i++)
    {
        lc[i] = olc[i];
        uc[i] = ouc[i];
        
        if( calcConstr )
            if( !calcConstr[i] )
                continue;
            
        bool calci = false;
        
        if( calcEvenOnOriginalValues || olc[i] <= -MIP_INFINITY )
        {
            //calculatin a lower bound to constraint i
            int r = solver->setConstraintLinearCoef(i, auxInd, -1.0 );
            r += solver->setObjLinearCoef(auxInd, 1.0);
            r += solver->setConstraintBounds(i, -MIP_INFINITY, 0.0 );
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                retCode = r;
                goto termination;
            }
            
            r = solver->solve(false, false, false, false);
            
            printf("solve para lc[%d]. status: %d obj: %lf\n", i, r, solver->getObjValue());
            
            if( r == OPT_OPTIMAL_SOLUTION )
            {
                lc[i] = OPT_max( solver->getObjValue(), olc[i] );
            }
            else
            {
                #if OPT_DEBUG_MODE
                    std::cerr << OPT_PREPRINT << "warning: error " << r << "to calculate lower bound to constraint " << i << OPT_GETFILELINE << "\n";
                #endif
                
                lc[i] = olc[i];
            }
            
            calci = true;
        }
        
        if( calcEvenOnOriginalValues || ouc[i] >= MIP_INFINITY )
        {
            int r = solver->setConstraintLinearCoef(i, auxInd, -1.0 );
            r += solver->setObjLinearCoef(auxInd, -1.0);
            r += solver->setConstraintBounds(i, 0.0, MIP_INFINITY );
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                retCode = r;
                goto termination;
            }
            
            r = solver->solve(false, false, false, false);
            
            printf("solve para uc[%d]. status: %d obj: %lf\n", i, r, solver->getObjValue());
            
            if( r == OPT_OPTIMAL_SOLUTION )
            {
                //objective has the maximum value of constraint times  -1
                uc[i] = OPT_min( -solver->getObjValue(), ouc[i] );
            }
            else
            {
                #if OPT_DEBUG_MODE
                    std::cerr << OPT_PREPRINT << "warning: error " << r << "to calculate upper bound to constraint " << i << OPT_GETFILELINE << "\n";
                #endif
                
                uc[i] = ouc[i];
            }
            
            calci = true;
        }
        
        
        if( calci )
        {
            //restoring constraint dualVarBound
            r = solver->setConstraintBounds(i, olc[i], ouc[i]);
            
            //removing the auxiliary variable...
            r+= solver->setConstraintLinearCoef(i, auxInd, 0.0);
            r+= solver->setObjLinearCoef(auxInd, 0.0);
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                retCode = r;
                goto termination;
            }
            
            auxInd++;
        }
    }
    
    
    retCode = 0;
termination:
    
    //evaluate object is local...
    if( prob.hasNlConstrs )
        ( (OPT_NLPSolver *) solver)->setNonLinearEvalObject(&eval);
    
    return retCode;
}



OPT_MINLPProbOrSolverUnifier::OPT_MINLPProbOrSolverUnifier()
{
    prob = NULL;
    solver = NULL;
}


void OPT_MINLPProbOrSolverUnifier:: setAsMINLPProblem(OPT_MINLPProb *prob)
{
    this->prob = prob;
    this->solver = NULL;
}


void OPT_MINLPProbOrSolverUnifier:: setAsSolver(OPT_LPSolver *solver)
{
    this->prob = NULL;
    this->solver = solver;
}


int OPT_MINLPProbOrSolverUnifier:: setVariableBounds(int varIndex, double lb, double ub)
{
    int r;
    
    if(prob)
    {
        int rlb, rub;
        
        rlb = prob->setVariableLowerBound(varIndex, lb);
        if(rlb != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(rlb);
            #endif
        }
        
        rub = prob->setVariableUpperBound(varIndex, ub);
        if(rub != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(rub);
            #endif
        }
        
        r = rlb != 0 ? rlb : rub; 
    }
    else
    {
        r = solver->setVariableBounds(varIndex, lb, ub);
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
    }
        
    return r;
}


int OPT_MINLPProbOrSolverUnifier:: setConstraintBounds(int constrIndex, double lb, double ub)
{
    int r;
    
    if(prob)
    {
        int rlb, rub;
        
        rlb = prob->setConstraintLowerBound(constrIndex, lb);
        if(rlb != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(rlb);
            #endif
        }
        
        rub = prob->setConstraintUpperBound(constrIndex, ub);
        if(rub != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(rub);
            #endif
        }
        
        r = rlb != 0 ? rlb : rub; 
    }
    else
    {
        r = solver->setConstraintBounds(constrIndex, lb, ub);
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
    }
    
    return r;
}


int OPT_MINLPProbOrSolverUnifier:: setConstraintQuadMatrix( int constrIndex, int nvars, const int *rows, const int *cols, const double *values)
{
    return prob ?  
    prob->setConstraintQuadCoefsMatrix( constrIndex, nvars, rows, cols, values) : 
    ((OPT_QCPSolver*) solver)->setConstraintQuadMatrix( constrIndex, nvars, rows, cols, values);
}


int OPT_MINLPProbOrSolverUnifier:: setConstraintLinearCoefs( int constrIndex, int nvars, const int *cols, const double *values)
{
    return prob ?  
    prob->setConstraintLinearPart(constrIndex, nvars, cols, values) :
    solver->setConstraintLinearCoefs( constrIndex, nvars, cols, values);
}




int OPT_SolutionDistanceConstraintSetter:: setDistConstraint( OPT_MINLPProb &prob, const int constrIndex, const int nvars, const int *varIndices, const double *sol, const double maxDistance, int *auxInds, double *auxVars)
{
    OPT_MINLPProbOrSolverUnifier unifier;
    
    unifier.setAsMINLPProblem(&prob);
    
    return setDistConstraint(unifier, constrIndex, nvars, varIndices, sol, maxDistance, auxInds, auxVars);
}



int OPT_SolutionDistanceConstraintSetter:: setDistConstraint( OPT_QCPSolver &solver, const int constrIndex, const int nvars, const int *varIndices, const double *sol, const double maxDistance, int *auxInds, double *auxVars)
{
    OPT_MINLPProbOrSolverUnifier unifier;
    
    unifier.setAsSolver(&solver);
    
    return setDistConstraint(unifier, constrIndex, nvars, varIndices, sol, maxDistance, auxInds, auxVars);
}



int OPT_SolutionDistanceConstraintSetter:: setDistConstraint( OPT_MINLPProbOrSolverUnifier &unifier, const int constrIndex, const int nvars, const int* varIndices, const double* sol, const double maxDistance, int *auxInds, double* auxVars )
{
    int r, code;
    int nLinCoefs = 0;
    int *pAuxInds, *myAuxInds = NULL;
    double *pAuxVars, *myAuxVars = NULL;
    double rhs = 0.0;
    
    
    if( auxInds )
    {
        pAuxInds = auxInds;
    }
    else
    {
        OPT_malloc(myAuxInds, nvars);
        OPT_IFMEMERRORGOTOLABEL(!myAuxInds, code, termination);
        
        pAuxInds = myAuxInds;
    }
    
    
    if(auxVars)
    {
        pAuxVars = auxVars;
    }
    else
    {
        OPT_malloc(myAuxVars, nvars);
        OPT_IFMEMERRORGOTOLABEL(!myAuxVars, code, termination);
        
        pAuxVars = myAuxVars;
    }
    
    //setting quadratic terms
    OPT_setAllArray(nvars, pAuxVars, 2.0);
    
    r = unifier.setConstraintQuadMatrix( constrIndex, nvars, varIndices, varIndices, pAuxVars);
    
    if(r != 0)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        code = r;
        goto termination;
    }
    
    
    //setting linear terms
    for(int i = 0; i < nvars; i++)
    {
        const int ind = varIndices[i];
        
        if( sol[ind] != 0.0 )
        {
            pAuxInds[nLinCoefs] =  ind;
            pAuxVars[nLinCoefs] = -2.0*sol[ind];
            nLinCoefs++;
        }
        
        //pAuxVars[i] = -2.0*sol[ind];
    }
    
    
    r = unifier.setConstraintLinearCoefs( constrIndex, nLinCoefs, pAuxInds, pAuxVars);
    OPT_IFERRORGOTOLABEL(r, code, r, termination);
    
    //setting the rhs term
    
    for(int i = 0; i < nvars; i++)
    {
        const int ind = varIndices[i];
        if( sol[ind] != 0.0 )
            rhs += sol[ind]*sol[ind];
    }
    
    //std::cout << "maxDistance: " << maxDistance << " rhs: " << rhs << "\n";
    
    rhs = maxDistance*maxDistance - rhs;
    
    
    r = unifier.setConstraintBounds( constrIndex, -OPT_INFINITY, rhs);
    OPT_IFERRORGOTOLABEL(r, code, r, termination);
    
    
    code = 0;
termination:
    
    if(myAuxInds)   free(myAuxInds);
    if(myAuxVars)	free(myAuxVars);
    
    return code;
}



int OPT_SolutionDistanceBoxConstraintsSetter:: setDistConstraint( OPT_MINLPProb &prob, const double *olx, const double *oux, const int nvars, const int *varIndices, const double *sol, const double percentageToNeighborhood )
{
    OPT_MINLPProbOrSolverUnifier unifier;
    
    unifier.setAsMINLPProblem(&prob);
    
    return setDistConstraint(unifier, olx, oux, nvars, varIndices, sol, percentageToNeighborhood);
}


int OPT_SolutionDistanceBoxConstraintsSetter:: setDistConstraint( OPT_LPSolver &solver, const double *olx, const double *oux, const int nvars, const int *varIndices, const double *sol, const double percentageToNeighborhood )
{
    OPT_MINLPProbOrSolverUnifier unifier;
    
    unifier.setAsSolver(&solver);
    
    return setDistConstraint(unifier, olx, oux, nvars, varIndices, sol, percentageToNeighborhood);
}


int OPT_SolutionDistanceBoxConstraintsSetter:: setDistConstraint( OPT_MINLPProbOrSolverUnifier &unifier, const double *olx, const double *oux, const int nvars, const int *varIndices, const double *sol, const double percentageToNeighborhood )
{
    const double absValueIfBoxIsInf = 1.0;
    int code = 0;
    
    
    for(int k = 0; k < nvars; k++)
    {
        const int ind = varIndices[k];
        int r;
        double d;
        
        
        if( olx[ind] > -OPT_INFINITY && oux[ind] < OPT_INFINITY )
            d = oux[ind] - olx[ind];
        else
            d = OPT_abs(sol[ind]) + absValueIfBoxIsInf ;
        
        
        const double lb = OPT_max( sol[ind] - percentageToNeighborhood*d, olx[ind] ); 
        const double ub = OPT_min( sol[ind] + percentageToNeighborhood*d, oux[ind] );
        
        
        r = unifier.setVariableBounds(ind, lb, ub);
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            code = r;
        }
        
    }
    
    return code;
}



int optsolvers::OPT_copyConstraintLinearParts(int beginOriginIndex, int endOriginIndex, int beginDestinIndex, OPT_LPSolver &origin, OPT_LPSolver &destin )
{
    int n;
    int r, retCode;
    int *auxCols = NULL;
    double *auxVals = NULL;
    
    r = origin.getNumberOfVars(n);
    OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    
    
    OPT_malloc(auxCols, n);
    OPT_malloc(auxVals, n);
    
    OPT_IFMEMERRORGOTOLABEL( !auxCols || !auxVals , retCode, termination);
    
    
    for( int i = beginOriginIndex; i <= endOriginIndex; i++ )
    {
        int nzs;
        int dind = beginDestinIndex + i  -beginDestinIndex ;
        double lc, uc;
        
        r = origin.getConstraintBounds(i, lc, uc);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        r = origin.getConstraintLinearPart(i, nzs, auxCols, auxVals);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        r = destin.setConstraintBounds(dind, lc, uc);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
        
        r = destin.setConstraintLinearCoefs(dind, nzs, auxCols, auxVals);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    }
    
    retCode = 0;
    
termination:
    
    if(auxCols)		free(auxCols);
    if(auxVals)		free(auxVals);
    
    return retCode;
}



int optsolvers::OPT_readParametersWithTypeFromFile(const char *fileName, const bool printErrorMsgs, const bool printFileOpenError, OPT_ParameterSetter &psetter)
{
    unsigned int lineCounter = 0;
    int r;
    FILE *file;
    char type[OPT_MAX_SIZE_OFSTRING_TO_READ_PARAMS];
    char param[OPT_MAX_SIZE_OFSTRING_TO_READ_PARAMS]; 
    char value[OPT_MAX_SIZE_OFSTRING_TO_READ_PARAMS];
    
    const int sline = 3*OPT_MAX_SIZE_OFSTRING_TO_READ_PARAMS + 10;
    char line[sline];
    
    file = fopen(fileName, "r");
    if( !file )
    {
        if(printFileOpenError)
            std::cerr << OPT_PREPRINT "Error at opening paremeter file " << fileName << OPT_GETFILELINE << "\n";
            
        return OPT_BAD_INPUT;
    }
    
    
    while( !feof(file) )
    {
        int retset;
        char *p;
        
        lineCounter++;
        
        p = fgets(line, sline, file);
        
        if(p == NULL) //so, the end of file was reached before read any charcter. That is so strange because was not detect at while test...
            break;
        
        
        if( line[0] == OPT_CHARAC_COMENT_ON_PARAMS_FILE ) //we assume there is only one peer to line
            continue;
        
        if( OPT_isEmptyString(line) ) //empty line
            continue;
        
        
        type[0] = '\0';
        param[0] = '\0';
        value[0] = '\0';
        
        
        r = sscanf(line, 
            "%" OPT_EXPSTR(OPT_MAX_SIZE_OFSTRING_TO_READ_PARAMS) "s"
            "%" OPT_EXPSTR(OPT_MAX_SIZE_OFSTRING_TO_READ_PARAMS) "s " "%" OPT_EXPSTR(OPT_MAX_SIZE_OFSTRING_TO_READ_PARAMS) "s", type, param, value );
        
        if( r < 3 )
        {
            if( printErrorMsgs )
                std::cerr << OPT_PREPRINT "Error to process line " << lineCounter << " in the file " << fileName << ": " << line; //do not \n here because line already have \n
            continue;
        }
        
        
        if( strcmp(type, "dbl") == 0 )
        {
            double dvalue;
            
            r = sscanf(value, "%lf", &dvalue);
            
            if( r == 0 )
            {
                if( printErrorMsgs )
                    std::cerr << OPT_PREPRINT "Error to process line " << lineCounter << " in the file " << fileName << ": " << line << "\t -  specified value " << value << " is not a real number\n";
                continue;
            }
            
            
            retset = psetter.setDoubleParameter( param, dvalue );
            
        }
        else if( strcmp(type, "int") == 0 )
        {
            long int ivalue;
            
            r = sscanf(value, "%ld", &ivalue);
            
            if( r == 0 )
            {
                if( printErrorMsgs )
                    std::cerr << OPT_PREPRINT "Error to process line " << lineCounter << " in the file " << fileName << ": " << line << "\t -  specified value " << value << " is not a integer number\n";
                continue;
            }
            
            retset = psetter.setIntegerParameter(param, ivalue);
        }
        else if( strcmp(type, "str") == 0 )
        {
            retset = psetter.setStringParameter(param, value);
        }
        else
        {
            std::cerr << OPT_PREPRINT "Error to process line " << lineCounter << " in the file " << fileName << "\t - invalid type " << type << "\n"; 
            continue;
        }
        
        
        if( printErrorMsgs )
        {
            if( retset == 0 )
            {
                std::cout << OPT_PREPRINT "Parameter " << param << " set to value " << value << "\n";
            }
            else
            {
                std::cerr << OPT_PREPRINT "Erro to set parameter " << param << " to value " << value << "\n";
            }
        }
        
    }
    
    std::cout << "\n";
    
    fclose(file);
    
    return 0;
}


