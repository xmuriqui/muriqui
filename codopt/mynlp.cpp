

#include <cstdlib>
#include <climits>
#include <iostream>
#include <map>

#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"


using namespace std;

using namespace minlpproblem;
using namespace optsolvers;




OPT_MyNLPSolver::OPT_MyNLPSolver():OPT_NLPSolver()
{
}


OPT_MyNLPSolver::~OPT_MyNLPSolver()
{
    deallocateMemory();
    deallocateSolverEnv();
}



// __ methods from Solver __


int OPT_MyNLPSolver::allocateConstrStructures(const int m)
{
    /*double *auxd;
    bool *auxb; */
    int r;
    
    /*int n;
    getNumberOfVars(n);*/
    
    
    /*auxb = (bool *) realloc( constrChg, m * sizeof(bool) );
    if( !auxb )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    constrChg = auxb;
    
    for(int i = maux; i < m; i++)
        constrChg[i] = false; */
    
    r = OPT_realloc(constrChg, m);
    if( m > maux)
        OPT_setAllArray( m - maux, &constrChg[maux], false );
    
    
    /*auxd = (double *) realloc( lambdaInit, m* sizeof(double) );
    if( !auxd )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    lambdaInit = auxd; */
    
    r = OPT_realloc(lambdaInit, m);
    OPT_IFERRORRETURN(r, r);
    
    if( m > 0 )
        lambdaInit[0] = NAN;
    
    
    return OPT_NLPSolver::allocateConstrStructures(m);
}



int OPT_MyNLPSolver::allocateVarStructures(const int n)
{
    int r;
    //double *auxd;
    //bool *auxb;
    //int m;
    
    //getNumberOfConstraints(m);
    
    /*auxb = (bool *) realloc( quadObjChg, n * sizeof(bool) );
    
    if( !auxb )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    quadObjChg = auxb;
    for(int i = naux; i < n; i++)
        quadObjChg[i] = false; */
    
    r = OPT_realloc( quadObjChg, n );
    OPT_IFERRORRETURN(r, r);
    
    if( n > naux )
        OPT_setAllArray(n-naux, &quadObjChg[naux], false);
    
    
    /*auxb = (bool *) realloc(hessChg, n* sizeof(bool));
    if( !auxb )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    hessChg = auxb;
    for(int i = naux; i < n; i++)
        hessChg[i] = false;*/
    
    
    r = OPT_realloc(hessChg, n);
    OPT_IFERRORRETURN(r, r);
    
    if(n > naux)
        OPT_setAllArray(n-naux, &hessChg[naux], false);
    
    
    /*for(int i = 0; i < m; i++)
    {
        auxb = (bool *) realloc( quadConstrChg[i], n* sizeof(bool) );
        
        if( !auxb )
            return OPT_MEMORY_ERROR;
        quadConstrChg[i] = auxb;
        
        for( int j = naux; j < n; j++)
            auxb[j] = false; //taking advantage of auxb
    } */
    
    
    /*auxb = (bool *) realloc( rowQuadConstrChg, n * sizeof(bool) );
    if( !auxb )
        return OPT_MEMORY_ERROR;
    rowQuadConstrChg = auxb;
    
    for(int i = naux; i < n; i++)
        rowQuadConstrChg[i] = false;*/
    
    
    r = OPT_realloc(rowQuadConstrChg, n);
    OPT_IFERRORRETURN(r, r);
    if(n > naux)
        OPT_setAllArray(n-naux, &rowQuadConstrChg[naux], false);
    
    
    /*auxd = (double *) realloc( xInit, n * sizeof(double) );
    if( !auxd )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    xInit = auxd; */
    
    r = OPT_realloc(xInit, n);
    OPT_IFERRORRETURN(r, r);
    
    
    /*auxd = (double *) realloc( zInit, 2* n * sizeof(double) );
    if( !auxd )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    zInit = auxd;*/
    
    r = OPT_realloc(zInit, 2* n);
    OPT_IFERRORRETURN(r, r);
    
    
    if(n > 0)
    {
        xInit[0] = NAN;
        zInit[0] = NAN;
    }
    
    
    
    return OPT_NLPSolver::allocateVarStructures(n);
}



void OPT_MyNLPSolver::deallocateMemory()
{
    OPT_secFree(quadObjChg);
    OPT_secFree(constrChg);
    
    /*if( quadConstrChg )
    {
        for(int i = 0; i < maux; i++)
        {
            if( quadConstrChg[i] )
                free(quadConstrChg[i]);
        }
        
        free(quadConstrChg);
        quadConstrChg = NULL;
    } */
    
    OPT_secFree( rowQuadConstrChg );
    OPT_secFree( hessChg );
    
    
    OPT_secFree( xInit );
    OPT_secFree( zInit );
    OPT_secFree( lambdaInit );
    
    deallocateAuxDerivativeIndexStructures();
    
    OPT_NLPSolver::deallocateMemory();
}



void OPT_MyNLPSolver::deallocateSolverEnv()
{
    prob.deallocateMatrices();
    deallocateMemory();
}


void OPT_MyNLPSolver::getDualSolution( double *dualConstrs, double *dualVarBounds, const bool correctSignal )
{
    const bool minusLambdaOnLagran = getMinusLambdaOnLagran();
    const bool minusFactor = prob.objFactor < 0;
    
    const bool reverseSignal = correctSignal && minusLambdaOnLagran != minusFactor ; //exclusive or (XOR) between minusLambdaOnLagran and minusFactor
    
    
    if(dualConstrs)
    {
        int m = 0;
        getNumberOfConstraints(m);
        
        if( !reverseSignal )
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
        
        if( !reverseSignal )
        {
            OPT_copyArray(2*n, dualSolV, dualVarBounds);
        }
        else
        {
            OPT_copyArrayTimesFactor(2*n, -1, dualSolV, dualVarBounds);
        }
    }
    
    
}


int OPT_MyNLPSolver::getVariableType( const int index, OPT_VARTYPE &varType )
{
    int r;
    MIP_VARTYPE vt;
    
    r = prob.getVariableType( index, vt );
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    varType = vt == MIP_VT_INTEGER ? OPT_VT_INTEGER : OPT_VT_CONTINUOUS;
    
    return 0;
}


void OPT_MyNLPSolver::initialize()
{
    OPT_NLPSolver::initialize();
    
    useSPMtoJacAndHess = false; //we do not use the Sparse Matrices stored and OPT_NLPSolver...
    
    nmChg = false;
    constrsBndsChg = false;
    genQuadConstrChg = false;
    genConstrChg = false;
    genHessChg = false;
    varTypeChg = false;
    
    quadObjChg = NULL;
    constrChg = NULL;
    //quadConstrChg = NULL;
    rowQuadConstrChg = NULL;
    hessChg = NULL;
    
    xInit = NULL;
    zInit = NULL;
    lambdaInit = NULL;
    
    
    jacCols = NULL;
    jacRowStartIndex = NULL;
    hessCols = NULL;
    hessRowStartIndex = NULL;
}



int OPT_MyNLPSolver::setVariableType( const int index, const OPT_VARTYPE varType )
{
    const int r = prob.setVariableType(index, varType == OPT_VT_INTEGER ? MIP_VT_INTEGER : MIP_VT_CONTINUOUS );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    varTypeChg = true;
    
    return 0;
}





// __ methods from LPSolver __




int OPT_MyNLPSolver::__addConstraints(const int nm)
{
    deallocateAuxDerivativeIndexStructures();
    
    const int r = prob.addConstraints(nm);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    nmChg = true;
    
    return 0;
}



int OPT_MyNLPSolver::__addVariables(const int nn, const bool initFree)
{
    deallocateAuxDerivativeIndexStructures();
    const int r = prob.addVariables(nn);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    nmChg = true;
    
    return 0;
}




int OPT_MyNLPSolver::generateModelFile(const char* fileName)
{
    const unsigned int size = strlen(fileName);
    const unsigned int send3 = OPT_max<unsigned>( size -3u, 0);
    const unsigned int send4 = OPT_max<unsigned>( size -4u, 0);
    
    int r;
    
    if( strcmp( &fileName[send3],  ".lp" )   ==  0 )
        r = prob.writeMIQCPModelInLPFile(fileName);
    else if(  strcmp( &fileName[send4],   ".gms" )  == 0   )
        r = prob.writeMIQCPModelInGAMSFile(fileName, "optsolvers", NULL);
    else //(  strcmp( &fileName[send4],   ".mod" )  == 0   )
        r = prob.writeMIQCPModelInAMPLFile(fileName);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}


int OPT_MyNLPSolver::getConstraintBounds( const int index, double &lb, double &ub )
{
    const int r = prob.getConstraintBounds( index, lb, ub );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    if( lb <= -MIP_INFINITY )
        lb = -OPT_INFINITY;
    
    if( ub >= MIP_INFINITY )
        ub = OPT_INFINITY;
    
    return 0;
}


int OPT_MyNLPSolver::getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value)
{
    const int r = prob.getConstraintLinearCoef( constrIndex, varIndex, value);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}



int OPT_MyNLPSolver::getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values)
{
    const int r = prob.getConstraintLinearPart( constrIndex, &nzs, cols, values);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}



int OPT_MyNLPSolver::getObjConstant(double &objConstant) 
{
    objConstant = prob.getObjConstant(); //i think we must not muliply to objFactor because objFactor is only set Maximization or minimziation
    return 0;
}


int OPT_MyNLPSolver::getObjLinearCoef( const int index, double &value )
{
    const int r = prob.getObjLinearCoef(index, value);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}



int OPT_MyNLPSolver::getNumberOfConstraints(int &m)
{
    m = prob.getNumberOfConstraints();
    return 0;
}



int OPT_MyNLPSolver::getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs)
{
    const int r = prob.getNumberOfLinearCoefsInConstr( constrIndex, nzs);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}



int OPT_MyNLPSolver::getNumberOfIntVars(int &nI)
{
    nI = prob.getNumberOfIntegerVars();
    return 0;
}



int OPT_MyNLPSolver::getNumberOfVars(int &n)
{
    n = prob.getNumberOfVars();
    return 0;
}



int OPT_MyNLPSolver::getObjSense(OPT_OPTSENSE &sense)
{
    sense = prob.getObjFactor() >= 0 ? OPT_MINIMIZE : OPT_MAXIMIZE ;
    
    return 0;
}



int OPT_MyNLPSolver::getVariableBounds(const int index, double &lb, double &ub)
{
    
    const int r = prob.getVariableBounds( index, lb, ub );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    if( lb <= -MIP_INFINITY )
        lb = -OPT_INFINITY;
    
    if( ub >= MIP_INFINITY )
        ub = OPT_INFINITY;
    
    return 0;
}



int OPT_MyNLPSolver::removeVars(const int ninds, const int *indices )
{
    
    deallocateAuxDerivativeIndexStructures();
    
    const int r = prob.removeVars( ninds, indices );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    nmChg = true;
    
    
    return 0;
}



//that method is a disaster for our Sparse Matrix implementations (oriented by rows...)
int OPT_MyNLPSolver::setLinearColumn( const int varIndex, const int nzs, const int* rows, const double* values)
{
    const int r = prob.setConstraintsLinearColumn( varIndex, nzs, rows, values );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    genConstrChg = true;
    for(int i = 0; i < nzs; i++)
        constrChg[ rows[i] ] = true;
    
    return 0;
}




int OPT_MyNLPSolver::resetConstraintLinearPart( const int constrIndex, const int nzs, const int *cols, const double *values )
{
    int r, code = 0;
    
    r = prob.setConstraintLinearPart( constrIndex, nzs, cols, values );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR;
    }	
    
    genConstrChg = true;
    constrChg[ constrIndex ] = true;
    
    return code;
}



int OPT_MyNLPSolver::setConstraintBounds( const int index, const double lb, const double ub )
{
    int code = 0, r;
    
    r = prob.setConstraintLowerBound( index, lb > -OPT_INFINITY ? lb : -MIP_INFINITY );
    
    constrsBndsChg = true;
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR;
    }
    
    
    r = prob.setConstraintUpperBound( index, ub < OPT_INFINITY ? ub : MIP_INFINITY );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR;
    }
    
    
    return code;
}





int OPT_MyNLPSolver::setConstraintsLinearCoefs( const int nzs, const int* rows, const int* cols, const double* values )
{
    int anzs, r, code = 0;
    //double *auxValues2 = NULL;
    unsigned int nanzs;
    unsigned int *nzrows = NULL;
    
    
    anzs = prob.getNumberOfLinearCoefs( );
    
    if( anzs > 0 )
    {
        //const int n = prob.getNumberOfVars();
        const int m = prob.getNumberOfConstraints();
        
        OPT_calloc(nzrows, m);
        OPT_IFMEMERRORGOTOLABEL(!nzrows, code, termination);
        
        
        for(int k = 0; k < nzs; k++)
            nzrows[ rows[k] ]++;
        
        for(int i = 0; i < m; i++)
        {
            #if 0
            OPT_setAllArray(n, auxValues2, 0.0);
            
            for(int j = 0; j < nzs; j++)
            {
                if( rows[j] == i )
                    auxValues2[cols[j]] = values[j];
            }
            
            
            r = prob.getConstraintLinearPart( i, &anzs, auxIndex, auxValues);
            
            #if OPT_DEBUG_MODE
                if( r != 0 )
                {
                    #if OPT_DEBUG_MODE
                        OPT_PRINTERRORNUMBER(r);
                    #endif
                    code = OPT_BAD_INPUT;
                    continue;
                }
            #endif
            
            //mixing old and new coefficients...
            for(int j = 0; j < anzs; j++)
                auxValues2[ auxIndex[j] ] = auxValues[j];
            
            
            nanzs = 0;
            for( int j = 0; j < n; j++ )
            {
                if( auxValues2[j] != 0.0 )
                {
                    auxIndex[ nanzs ] = j; 
                    auxValues[nanzs] = auxValues2[j];
                    nanzs++;
                }
            }
            #endif
            
            if( nzrows[i] == 0 )
                continue;
            
            
            unsigned int mynzrow = 0;
            
            for(int j = 0; j < nzs; j++)
            {
                if( rows[j] == i )
                {
                    auxIndex2[mynzrow] = cols[j];
                    auxValues2[mynzrow] = values[j];
                    mynzrow++;
                }
            }
            
            #if OPT_DEBUG_MODE
                assert(mynzrow == nzrows[i]);
            #endif
            
            r = prob.getConstraintLinearPart( i, &anzs, auxIndex, auxValues);
            OPT_IFERRORGOTOLABEL(r, code, OPT_UNDEFINED_ERROR, termination);
            
            r = OPT_unifyArraysOfIndicesAndValues( anzs, auxIndex, auxValues, mynzrow, auxIndex2, auxValues2, nanzs, auxIndex, auxValues );
            OPT_IFERRORGOTOLABEL(r, code, r, termination);
            
            r = prob.setConstraintLinearPart(i, nanzs, auxIndex, auxValues );
            OPT_IFERRORGOTOLABEL(r, code, OPT_UNDEFINED_ERROR, termination);
            
            constrChg[i] = true;
            
        }
        
    }
    else
    {
        r = prob.setConstraintsLinearPart( nzs, rows, cols, values );
        OPT_IFERRORRETURN(r, OPT_UNDEFINED_ERROR);
        
        for(int i = 0; i < nzs; i++)
            constrChg[ rows[i] ] = true;
    }
    
    
termination:
    
    genConstrChg = true;
    
    if(nzrows)	free(nzrows);

    return code;
}



int OPT_MyNLPSolver::setConstraintLinearCoefs( const int constrIndex, const int nzs, const int *cols, const double *values)
{
    //const int n = prob.getNumberOfVars();
    int r, anzs;
    
    r = prob.getConstraintLinearPart( constrIndex, &anzs, auxIndex, auxValues );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    if( anzs > 0 )
    {
        unsigned int mynzs = 0;
        
        #if 0
        {
            OPT_setAllArray(n, auxValues2, 0.0);
            
            for(int i = 0; i < anzs; i++)
                auxValues2[ auxIndex[i] ] = auxValues[i];
            
            for(int i = 0; i < nzs; i++)
                auxValues2[ cols[i] ] = values[i];
            
            
            for( int i = 0; i < n; i++ )
            {
                if( auxValues2[i] != 0.0 )
                {
                    auxIndex[ mynzs ] = i;
                    auxValues[ mynzs ] = auxValues2[i];
                    mynzs++;
                }
            }
        }
        #endif
        
        r = OPT_unifyArraysOfIndicesAndValues( anzs, auxIndex, auxValues, nzs, cols, values, mynzs, auxIndex2, auxValues2 );
        OPT_IFERRORRETURN(r, r);
        
        r = prob.setConstraintLinearPart( constrIndex, mynzs, auxIndex2, auxValues2 );
        OPT_IFERRORRETURN(r, OPT_SOLVER_ERROR);
    }
    else
    {
        r = prob.setConstraintLinearPart( constrIndex, nzs, cols, values );
        OPT_IFERRORRETURN(r, OPT_SOLVER_ERROR);
    }
    
    
    genConstrChg = true;
    constrChg[ constrIndex ] = true;
    
    
    return 0;
}



int OPT_MyNLPSolver::setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value)
{
    bool inStruct;
    int r;
    double v;
    
    r = prob.getConstraintLinearCoef( constrIndex, varIndex, v, &inStruct );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    if( inStruct )
    {
        r = prob.setConstraintLinearCoefInStructure( constrIndex, varIndex, value );
        
    }
    else
    {
        int nzs;
        
        r = prob.getConstraintLinearPart( constrIndex, &nzs, auxIndex, auxValues );
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
        
        auxIndex[nzs] = varIndex;
        auxValues[nzs] = value;
        
        r = prob.setConstraintLinearPart( constrIndex, nzs+1, auxIndex, auxValues );
    }
    
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    genConstrChg = true;
    constrChg[ constrIndex ] = true;
    
    return 0;
}




int OPT_MyNLPSolver::setObjLinearCoef( const int index, const double value )
{
    const int r = prob.setObjLinearCoefficient(index, value);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}



int OPT_MyNLPSolver::setObjLinearCoefs( const int nzs, const int* cols, const double* values )
{
    const int r = prob.setObjLinearCoefficients(nzs, cols, values);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}



int OPT_MyNLPSolver::setObjLinearPart( const int n, const double *values )
{
    int r;
    const int on = prob.getNumberOfVars();
    
    
    if( on < n )
        return OPT_BAD_INPUT;
    
    
    for(int i = 0; i < on; i++)
        auxIndex[i] = i;
    
    
    r = prob.setObjLinearCoefficients(n, auxIndex, values);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_BAD_INPUT;
    }
    
    
    if( n != on )
    {
        for(int i = n; i < on; i++)
            auxValues[i] = 0.0;
        
        r = prob.setObjLinearCoefficients(on-n, &auxIndex[n], &auxValues[n]);
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_BAD_INPUT;
        }
    }
    
    
    return 0;
}




void OPT_MyNLPSolver::setObjConstant(const double value)
{
    prob.setObjConstant(value);
}



int OPT_MyNLPSolver::setObjSense( const OPT_OPTSENSE sense )
{
    prob.setObjFactor( sense == OPT_MINIMIZE ? 1.0 : -1.0 );
    
    return 0;
}



int OPT_MyNLPSolver::setnVariablesBounds( const int n, const double *lb, const double *ub )
{
    int r = prob.setVariableLowerBounds(n, lb);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_BAD_INPUT;
    }
    
    
    r = prob.setVariableUpperBounds(n, ub);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}



int OPT_MyNLPSolver::setVariableBounds( const int index, const double lb, const double ub )
{
    int r;
    
    r = prob.setVariableLowerBound(index, lb) + prob.setVariableUpperBound(index, ub);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}


int OPT_MyNLPSolver::setVariablesBounds( const int ninds, const int *inds, const double *lb, const double *ub )
{
    int r = prob.setVariableLowerBounds(ninds, inds, lb);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_BAD_INPUT;
    }
    
    
    r = prob.setVariableUpperBounds(ninds, inds, ub);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}


//__ methods from QPSolver __



int OPT_MyNLPSolver::getNumberOfQuadObjTerms(int &nzs)
{
    nzs = prob.getNumberOfObjQuadTerms( );
    
    return 0;
}



int OPT_MyNLPSolver::getObjQuadTerm( const int row, const int col, double &value)
{
    const int r = prob.getObjQuadCoef(row, col, value );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}



int OPT_MyNLPSolver::getObjQuadPart( int &nzs, int *rows, int *cols, double *values )
{
    /*const int r = prob.getObjQuadCoefsMatrix( &nzs, rows, cols, values );
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_BAD_INPUT;
    }*/
    
    prob.getObjQuadCoefsMatrix( &nzs, rows, cols, values );
    
    return 0;
}


int OPT_MyNLPSolver::setObjQuadCoef( const int row, const int col, const double value )
{
    int r;
    int nzs = 0;
    
    int myrow = row, mycol = col;
    
    if( myrow < mycol )
        OPT_swap( myrow, mycol );
    
    
    r = prob.getObjQuadCoefsMatrixRow( myrow, &nzs, auxIndex, auxValues );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    r = 0;
    
    for( int i = 0; i < nzs; i++ )
    {
        if( auxIndex[i] == mycol )
        {
            auxValues[i] = value;
            r = 1;
            break;
        }
    }
    
    if( r == 0 )
    {
        auxIndex[nzs] = mycol;
        auxValues[nzs] = value;
        nzs++;
    }
    
    r = prob.setObjQuadCoefsMatrixRow( myrow, nzs, auxIndex, auxValues );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    genHessChg = true;
    quadObjChg[myrow] = true;
    
    return 0;
}



int OPT_MyNLPSolver::setObjQuadMatrix( const int nzs, const int* rows, const int* cols, const double* values )
{
    int r; 
    
    #if 0
    int anzs = prob.getNumberOfObjQuadTerms();
    
    if( anzs > 0 )
    {
        int n;
        
        getNumberOfVars(n);
        
        
        for(int i = 0; i < n; i++)
        {
            r = prob.getObjQuadCoefsMatrixRow(i, anzs, auxIndex, auxValues );
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
                #endif
                
                return OPT_SOLVER_ERROR;
            }
            
            if( anzs > 0 )
            {
                OPT_setAllArray(n, auxValues2, 0.0);
                
                for(int j = 0; j < anzs; j++)
                    auxValues2[ auxIndex[j] ] = auxValues[j];
                
                
                for(int j = 0; j < nzs; j++)
                {
                    if( rows[j] == i )
                        auxValues2[ cols[j] ] = values[j];
                }
                
                
                anzs = 0;
                for(int j = 0; j < n; j++)
                {
                    if( auxValues2[j] != 0.0 )
                    {
                        auxIndex[anzs] = j;
                        auxValues[anzs] = auxValues2[j];
                        anzs++;
                    }
                }
                
            }
            else
            {
                anzs = 0;
                for(int j = 0; j < nzs; j++)
                {
                    if( rows[j] == i )
                    {
                        auxIndex[anzs] = cols[j];
                        auxValues[anzs] = values[j];
                        anzs++;
                    }
                }
            }
            
            
            r = prob.setObjQuadCoefsMatrixRow(i, anzs, auxIndex, auxValues);
                
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                return OPT_SOLVER_ERROR;
            }
            
        }
        
        genHessChg = true;
        OPT_setAllArray(n, quadObjChg, true);
        
    }
    else
    #endif
    {
        OPT_setNonzeroRows(prob.Q, quadObjChg);
        
        r = prob.setObjQuadCoefsMatrix( nzs, rows, cols, values );
        
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            return OPT_SOLVER_ERROR;
        }
        
        genHessChg = true;
        
        OPT_setNonzeroRows(prob.Q, quadObjChg);
    }
    
    return 0;
}



int OPT_MyNLPSolver::setObjQuadMatrix(const int *rowStart, const int *cols, const double *values)
{
    //const int n = prob.n;
    
    
    OPT_setNonzeroRows(prob.Q, quadObjChg);
    
    
    const int r = prob.setObjQuadCoefsMatrix(rowStart, cols, values);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return OPT_SOLVER_ERROR;
    }
    
    
    genHessChg = true;
    OPT_setNonzeroRows(prob.Q, quadObjChg);
    
    
    return 0;
}




// __ methods from QCPSolver __


int OPT_MyNLPSolver::getNumberOfConstraintQuadTerms( const int index, int &nzs)
{
    const int r = prob.getNumberOfConstraintQuadCoefMatrixTerms( index, nzs);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}



int OPT_MyNLPSolver::getConstraintQuadMatrix( const int index, int &nzs, int *rows, int *cols, double *values )
{
    const int r = prob.getConstraintQuadCoefMatrix(index, &nzs, rows, cols, values);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}



int OPT_MyNLPSolver::setConstraintQuadMatrix( const int index, const int nzs, const int *qrows, const int *qcols, const double *qvalues )
{
    OPT_setNonzeroRows(prob.QC[index], rowQuadConstrChg);
    
    
    const int r = prob.setConstraintQuadCoefsMatrix( index, nzs, qrows, qcols, qvalues );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return OPT_SOLVER_ERROR;
    }
    
    genConstrChg = true;
    genQuadConstrChg = true;
    genHessChg = true;
    
    constrChg[index] = true;
    
    /*for(int i = 0; i < nzs; i++)
        rowQuadConstrChg[ qrows[i] ] = true; //quadConstrChg[index][ qrows[i] ] = true; */
    
    OPT_setNonzeroRows(prob.QC[index], rowQuadConstrChg);
    
    return 0;
}


int OPT_MyNLPSolver::setConstraintQuadMatrix( const int index, const int *qrowStart, const int *qcols, const double *qvalues )
{
    OPT_setNonzeroRows(prob.QC[index], rowQuadConstrChg);
    
    const int r = prob.setConstraintQuadCoefsMatrix(index, qrowStart, qcols, qvalues);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return OPT_SOLVER_ERROR;
    }
    
    genConstrChg = true;
    genQuadConstrChg = true;
    genHessChg = true;
    
    constrChg[index] = true;
    
    OPT_setNonzeroRows(prob.QC[index], rowQuadConstrChg);
    
    return 0;
}



// __ methods from OPT_NLPSolver __


int OPT_MyNLPSolver::__removeConstraints(const int ninds, const int* indices )
{
    deallocateAuxDerivativeIndexStructures();
    
    const int r = prob.removeConstraints( ninds, indices );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        return OPT_SOLVER_ERROR;
    }
    
    nmChg = true;
    
    return 0;
}


int OPT_MyNLPSolver::getLagrangianHessianStructure(int &nzs, int* rows, int* cols)
{
    int r;
    
    r = prob.getLagrangianHessianStructure(&nzs, rows, cols);
    OPT_IFERRORRETURN(r, OPT_UNDEFINED_ERROR);
    
    return 0;
}


int OPT_MyNLPSolver::getJacobianStructure(int &nzs, int* rows, int* cols)
{
    int r;
    
    r = prob.getJacobianStructure(&nzs, rows, cols);
    OPT_IFERRORRETURN(r, OPT_UNDEFINED_ERROR);
    
    return 0;
}


int OPT_MyNLPSolver::getConstrNLFlag(const int index, bool &flag)
{
    const int r = prob.getConstraintsNonLinearTermFlag(index, flag);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}


int OPT_MyNLPSolver:: getNonLinearEvalObjectPointer( OPT_NonLinearEval* &nlEval )
{
    nlEval = prob.getNonLinearEvaluationObject();
    return 0;
}


int OPT_MyNLPSolver::getNumberOfNLConstraints(int& mnl)
{
    mnl = prob.getNumberOfNLConstraints();
    
    return 0;
}


int OPT_MyNLPSolver::getNumberOfNonZerosInLagrangianHessian(int &nzs)
{
    nzs = prob.getNumberOfLagrangianHessianNonZeros();
    return 0;
}


int OPT_MyNLPSolver::getNumberOfNonZerosInJacobian(int &nzs)
{
    nzs = prob.getNumberOfJacobianNonZeros();
    return 0;
}


int OPT_MyNLPSolver::getObjNLFlag(bool &flag)
{
    flag = prob.hasObjNLTerm();
    return 0;
}


//if this method return true, we implement the solver class derived from OPT_MyNLPSolver, i.e.,  using MIP_MINLPProb to store coefficients
bool OPT_MyNLPSolver::isMyNLPClass()
{
    return true;
}



int OPT_MyNLPSolver::setConstrNLFlag(const int index, const bool flag)
{
    const int r = prob.setConstraintNonLinearTermFlag( index, flag);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}


int OPT_MyNLPSolver::setJacobianStructure(const int nzs, const int* rows, const int* cols)
{
    //old coefficients will be deleted. So, we mark their respective rows like changing rows.
    
    OPT_setNonzeroRows(prob.J, constrChg);
    
    
    const int r = prob.setJacobianStructure( nzs, rows, cols );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    genConstrChg = true;
    OPT_setNonzeroRows(prob.J, constrChg);
    //for(int i = 0; i < nzs; i++)
        //constrChg[ rows[i] ] = true;
    
    
    return 0;
}



int OPT_MyNLPSolver::setJacobianStructure(const int* rowStart, const int* cols)
{
    //old coefficients will be deleted. So, we mark their respective rows like changing rows.
    OPT_setNonzeroRows(prob.J, constrChg);
    
    
    const int r = prob.setJacobianStructure(rowStart, cols);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    genConstrChg = true;
    
    OPT_setNonzeroRows(prob.J, constrChg);
    
    return 0;
}




int OPT_MyNLPSolver::setInitialSolution(const double* x, const double* dualConstrs, const double* dualVars)
{
    int n = 0, m = 0;
    
    getNumberOfVars(n);
    getNumberOfConstraints(m);
    
    
    OPT_setInicialSolution(n, m, x, dualConstrs, dualVars, xInit, lambdaInit, zInit);
    
    /*if( x )
    {
        OPT_copyArray(n, x, xInit);
    }
    else
    {
        if( n >= 1 )
            xInit[0] = NAN;
    }
    
    
    if( dualConstrs )
    {
        OPT_copyArray(m, dualConstrs, lambdaInit);
    }
    else
    {
        if( m >= 1 )
            lambdaInit[0] = NAN;
    }
    
    
    if( dualVars )
    {
        OPT_copyArray(n, dualVars, zInit);
    }
    else
    {
        if( n >= 1 )
            zInit[0] = NAN;
    }*/
    
    return 0;
}





int OPT_MyNLPSolver::setJacobianRowStructure( const int row, const int nzs, const int* cols)
{
    const int r = prob.setJacobianRowStructure( row, nzs, cols);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    genConstrChg = true;
    constrChg[ row ] = true;
    
    return 0;
}



int OPT_MyNLPSolver::setLagrangianHessianStructure( const int nzs, const int* rows, const int* cols)
{
    //old coefficients will be deleted. So, we mark their respective rows like changing rows.
    OPT_setNonzeroRows(prob.lagH, hessChg);
    
    const int r = prob.setLagrangianHessianStructure( nzs, rows, cols);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    genHessChg = true;
    OPT_setNonzeroRows(prob.lagH, hessChg);
    //for(int i = 0; i < nzs; i++)
        //hessChg[ rows[i] ] = true;
    
    return 0;
}


int OPT_MyNLPSolver::setLagrangianHessianStructure( const int *rowStart, const int* cols)
{
    //old coefficients will be deleted. So, we mark their respective rows like changing rows.
    OPT_setNonzeroRows(prob.lagH, hessChg);
    
    const int r = prob.setLagrangianHessianStructure(rowStart, cols);
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    genHessChg = true;
    OPT_setNonzeroRows(prob.lagH, hessChg);
    //for(int i = 0; i < nzs; i++)
        //hessChg[ rows[i] ] = true;
    
    return 0;
}


int OPT_MyNLPSolver::setLagrangianHessianRowStructure( const int row, const int nzs, const int* cols)
{
    const int r = prob.setLagrangianHessianRowStructure( row, nzs, cols );
    
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    genHessChg = true;
    hessChg[row] = true;
    
    return 0;
}



int OPT_MyNLPSolver::setNonLinearEvalObject( OPT_NonLinearEval* nlEval )
{
    prob.setNonLinearEvaluationObject(nlEval);
    return 0;
}



int OPT_MyNLPSolver::setObjNLFlag(const bool flag)
{
    prob.setObjNonLinearTermFlag(flag);
    
    return 0;
}













int OPT_MyNLPSolver::allocateAuxDerivativeIndexStructures(  )
{
    int n, m;
    
    getNumberOfVars(n);
    getNumberOfConstraints(m);
    
    
    
    
    jacRowStartIndex = (unsigned int *) calloc( m+1 , sizeof(unsigned int) );
    jacCols = (int **) malloc( m * sizeof(int *) );
    
    hessRowStartIndex = (unsigned int *) calloc( n+1 , sizeof(unsigned int) );
    hessCols = (int **) malloc( n * sizeof(int *) );
    
    if( !jacRowStartIndex || !jacCols || !hessRowStartIndex || !hessCols )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        return OPT_MEMORY_ERROR;
    }
    
    OPT_setAllArray<int *>( m, jacCols, NULL);
    OPT_setAllArray<int *>(n, hessCols, NULL);
    
    return 0;
}



void OPT_MyNLPSolver::deallocateAuxDerivativeIndexStructures()
{
    OPT_secFree( jacRowStartIndex );
    if( jacCols )
    {
        int m;
        getNumberOfConstraints(m);
        
        for( int i = 0; i < m; i++ )
        {
            //std::cout << i << endl;
            OPT_secFree( jacCols[i] );
        }
        
        free( jacCols );
        jacCols = NULL;
    }
    
    
    OPT_secFree( hessRowStartIndex );
    if( hessCols )
    {
        int n;
        getNumberOfVars(n);
        
        for( int i = 0; i < n; i++ )
            OPT_secFree( hessCols[i] );
        
        free( hessCols );
        hessCols = NULL;
    }
}


#if 0

int OPT_MyNLPSolver::setJacIndexRow( const int rowIndex, int &nzs )
{
    int r, n, mynzs;
    unsigned int *row;
    MIP_SparseMatrix &QC = prob.QC[rowIndex]; //ok, that is not good, but I think access QC directly is the best way to do it...
    
    getNumberOfVars(n);
    
    
    row = &old_jacIndexSh[ rowIndex*n ];
    
    OPT_setAllArray(n, row, UINT_MAX);
    
    
    r = prob.getConstraintLinearPart( rowIndex, mynzs, auxIndex, auxValues );
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            return OPT_UNDEFINED_ERROR;
        }
    #endif
    
    
    for(int j = 0; j < mynzs; j++)
        row[ auxIndex[j] ] = 0; 
    
    
    r = prob.getJacobianRowStructure( rowIndex, mynzs, auxIndex );
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
            #endif
            
            return OPT_UNDEFINED_ERROR;
        }
    #endif
    
    
    for(int j = 0; j < mynzs; j++)
        row[ auxIndex[j] ] = 0;
    
    
    if( QC.getNumberOfElements() > 0 )
    {
        int *nzCol = auxIndex2;
        
        QC.countRowsEachColumn(nzCol);
        
        for(int j = 0; j < n; j++)
        {
            //checking if line in QC has some nonzero element. Note Sparse Matrix Only Store lower triangle, So we we have to test if line and column j has a nozero element...
            if( nzCol[j] > 0 ||  QC.getNumberOfElementsAtRow(j) > 0 )
            {
                row[j] = 0;
            }
        }
    }
    
    
    //setting now the indexes for this row. Note, this indices are relatives to get the correct indixes for ipopt triple sparse format, we have to sum those relatice indixes with jacIndexBase[  ]. Note, positions in jacobian are marked as 0 (remainder indexes are marked with UINT_MAX)
    
    nzs = 0;
    for( int j = 0; j < n; j++ )
    {
        if( row[j] == 0 )
        {
            row[j] = nzs;
            nzs++;
        }
    }
    
    
    return 0;
}


int OPT_MyNLPSolver::setHessIndexRow( const int rowIndex, int &nzs )
{
    int r, n, m, mynzs;
    unsigned int *row;
    
    
    getNumberOfVars(n);
    getNumberOfConstraints(m);
    
    
    row = &old_hessIndexSh[ (rowIndex*(rowIndex+1))/2 ]; //just lower triangle...
    
    OPT_setAllArray(rowIndex + 1, row, UINT_MAX);
    
    
    r = prob.getLagrangianHessianRowStructure( rowIndex, mynzs, auxIndex );
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
            #endif
            
            return OPT_UNDEFINED_ERROR;
        }
    #endif
    
    for(int j = 0; j < mynzs; j++)
        row[ auxIndex[j] ] = 0;
    
    
    r = prob.getObjQuadCoefsMatrixRow(rowIndex, mynzs, auxIndex, auxValues);
    
    #if OPT_DEBUG_MODE
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
            #endif
            
            return OPT_UNDEFINED_ERROR;
        }
    #endif
    
    for(int j = 0; j < mynzs; j++)
        row[ auxIndex[j] ] = 0;
    
    
    
    for(int i = 0; i < m; i++)
    {
        r = prob.getConstraintQuadCoefMatrixRow(i, rowIndex, mynzs, auxIndex, auxValues );
        
        #if OPT_DEBUG_MODE
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
                #endif
                
                return OPT_UNDEFINED_ERROR;
            }
        #endif
        
        for(int j = 0; j < mynzs; j++)
            row[ auxIndex[j] ] = 0;
    }
    
    
    //setting now the indexes for this row. Note, this indices are relatives to get the correct indixes for ipopt triple sparse format, we have to sum those relatice indixes with jacIndexBase[  ]. Note, positions in jacobian are marked as 0 (remainder indexes are marked with UINT_MAX)
    
    nzs = 0;
    for( int j = 0; j <= rowIndex; j++ ) //symmetric matrix
    {
        if( row[j] == 0 )
        {
            row[j] = nzs;
            nzs++;
        }
    }
    
    
    return 0;
}

#endif




int OPT_MyNLPSolver::setJacIndexRow( const int i, int* &jacColsi , unsigned int &nzl )
{
    int n, m;
    
    
    const bool nlConstr = prob.nlConstr[i];
    bool *vars = (bool *)auxIndex2 ;
    int *auxp;
    const MIP_SparseMatrix *QC = prob.QC;
    //MIP_SparseRow &arow = prob.A[i];
    //MIP_SparseRow &jrow = prob.J[i];
    
    const MIP_SparseMatrix &A = prob.A;
    const MIP_SparseMatrix &J = prob.J;
    
    const int nzAi = A.getNumberOfElementsAtRow(i);
    const int nzJi = J.getNumberOfElementsAtRow(i);
    
    getNumberOfVars(n);
    getNumberOfConstraints(m);
    
    
    OPT_secFree( jacColsi );
    
    if( nlConstr && nzAi == 0 && QC[i].getNumberOfElements() == 0 )
    {
        nzl = nzJi;
        
        if( nzl > 0 )
        {
            jacColsi = (int *) malloc( nzl * sizeof(int) );
        
            if( !jacColsi )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                return OPT_MEMORY_ERROR;;
            }
            
            //jrow.getStructure( jacColsi );
            J.getRowStructure(i, jacColsi, NULL );
        }
    }
    else if( (!nlConstr || nzJi == 0)  && QC[i].getNumberOfElements() == 0 )
    {
        nzl = nzAi;
        
        if( nzl > 0 )
        {
            jacColsi = (int *) malloc( nzl * sizeof(int) );
        
            if( !jacColsi )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                return OPT_MEMORY_ERROR;;
            }
            
            A.getRowStructure(i, jacColsi, NULL);
        }
    }
    else
    {
        nzl = 0;
        
        OPT_setAllArray(n, vars, false);
        
        if( QC[i].getNumberOfElements() > 0 )
        {
            const unsigned int nrows = QC[i].getNumberOfRows();
            
            /*for( int j = 0; j < nrows; j++ )
            {
                //unsigned int nel = QC[i][j].getStructure( vars ); //for each nondiagonal position in QC[i], we have two positions to set in the gradient
                
                unsigned int nel = QC[i].getRowStructure(j, vars, NULL, true ); //for each nondiagonal position in QC[i], we have two positions to set in the gradient
                
                if( nel > 0u )
                    vars[j] = true; //if there is some column in line j, both the column and j must be set as true 
            }*/
            
            
            for(MIP_SparseMatrixRowIndexIterator it = QC[i].beginRowIndex(); *it < nrows ; ++it )
            {
                const unsigned int j = *it;
                
                unsigned int nel = QC[i].getRowStructure(j, vars, NULL, true ); //for each nondiagonal position in QC[i], we have two positions to set in the gradient
                
                if( nel > 0u )
                    vars[j] = true; //if there is some column in line j, both the column and j must be set as true 
            }
            
        }
        
        //arow.getStructure( vars );
        A.getRowStructure(i, vars, NULL, true);
        
        //prob.J[i].getStructure( vars );
        J.getRowStructure(i, vars, NULL, true);
        
        
        for( int j = 0; j < n; j++ )
        {
            if( vars[j] )
                nzl++;
        }
        
        
        if( nzl > 0 )
        {
            jacColsi = (int *) malloc( nzl * sizeof(int) );
            
            if( !jacColsi )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                return OPT_MEMORY_ERROR;;
            }
            
            auxp = jacColsi;
            nzl = 0;
            for( int j = 0; j < n; j++ )
            {
                if( vars[j] )
                {
                    auxp[nzl] = j;
                    nzl++;
                }
            }
        }
    
    }
    
    return 0;
}
    
    



int OPT_MyNLPSolver::setFullJacIndex( unsigned int *jacRowStart, int** jacCols  )
{
    int m, r;
    unsigned int nz = 0, nzl;
    
    
    getNumberOfConstraints(m);
    
    for( int i = 0; i < m; i++ )
    {
        jacRowStart[i] = nz;
        
        r = setJacIndexRow(i, jacCols[i], nzl);
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            return OPT_MEMORY_ERROR;
        }
        
        nz += nzl;
    }
    
    
    jacRowStart[m] = nz;
    
    return 0;
}




unsigned int OPT_MyNLPSolver::setHessIndexRow( const int i, const int mquad, const int *quadIndex, int* &hessColsi, unsigned int &nzl ) 
{
    int nqconsi = 0; //number of matrices in QC having some element in row i
    int indqconsi;
    int n, m;
    
    //unsigned int *auxp;
    bool *vars = (bool *)auxIndex2 ;
    
    const MIP_SparseMatrix *QC = prob.QC;
    const MIP_SparseMatrix &Q = prob.Q;
    const MIP_SparseMatrix &lagH = prob.lagH;
    
    
    getNumberOfVars(n);
    getNumberOfConstraints(m);
    
    //unsigned int aux;
    
    
    for(int j = 0; j < mquad; j++)
    {
        //if( QC[ quadIndex[j] ][i].getNumberOfElements() > 0 )
        if( QC[ quadIndex[j] ].getNumberOfElementsAtRow(i) > 0 )
        {
            nqconsi++;
            if( nqconsi > 1 )
                break;
            
            indqconsi = quadIndex[j];
        }
    }
    
    
    if( Q.getNumberOfElementsAtRow(i) == 0 && nqconsi == 0 )
    {
        //nzl = lagH[i].getNumberOfElements();
        nzl = lagH.getNumberOfElementsAtRow(i);
        
        if( nzl > 0 )
        {
            hessColsi = (int *) malloc( nzl * sizeof(int) );
            
            if( !hessColsi )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                return OPT_MEMORY_ERROR;
            }
            
            lagH.getRowStructure(i, hessColsi, NULL);
        }
    }
    else if( lagH.getNumberOfElementsAtRow(i) == 0 && nqconsi == 0 )
    {
        //nzl = Q[i].getNumberOfElements();
        nzl = Q.getNumberOfElementsAtRow(i);
        
        if( nzl > 0 )
        {
            hessColsi = (int *) malloc( nzl * sizeof(int) );
            
            if( !hessColsi )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                return OPT_MEMORY_ERROR;
            }
            
            //Q[i].getStructure( hessColsi );
            Q.getRowStructure(i, hessColsi, NULL );
        }
    }
    else if( Q.getNumberOfElementsAtRow(i) == 0 && lagH.getNumberOfElementsAtRow(i) == 0 && nqconsi == 1 )
    {
        const MIP_SparseMatrix &QCc = QC[ indqconsi ];
        
        //nzl = QCc[i].getNumberOfElements();
        nzl = QCc.getNumberOfElementsAtRow(i);
        
        if( nzl > 0 )
        {
            hessColsi = (int *) malloc( nzl * sizeof(int) );
            
            if( !hessColsi )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                return OPT_MEMORY_ERROR;
            }
            
            //QCc[i].getStructure( hessColsi );
            QCc.getRowStructure(i, hessColsi, NULL);
        }
    }
    else
    {
        nzl = 0;
        int nzs;
        
        OPT_setAllArray(n, vars, false);
        
        if( Q.getNumberOfElements() > 0 )
        {
            //aux = Q[i].getStructure(vars);
            Q.getRowStructure(i, vars, &nzs, true);
            if( nzs > 0 )
                vars[i] = true;
        }
        
        
        for( int j = 0; j < mquad; j++ )
        {
            const int ind = quadIndex[j];
            
            //if(QC[ind].getNumberOfElements() > 0)
            {
                //aux = QC[ind][i].getStructure(vars);
                QC[ind].getRowStructure(i, vars, &nzs, true);
                if( nzs > 0 )
                    vars[i] = true;
            }
        }
        
        if( lagH.getNumberOfElements() > 0 )
        {
            //aux = lagH[i].getStructure(vars);
            lagH.getRowStructure(i, vars, &nzs, true);
            if( nzs > 0 )
                vars[i] = true;
        }
        
        
        for(int k = 0; k < n; k++ )
        {
            if( vars[k] )
                nzl++;
        }
        
        
        OPT_secFree( hessColsi );
        
        if( nzl > 0 )
        {
            hessColsi = (int *) malloc( nzl * sizeof(int) );
            
            if( !hessColsi )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                return OPT_MEMORY_ERROR;
            }
            
            int *auxp = hessColsi;
            
            nzl = 0;
            for( int k = 0; k < n; k++)
            {
                if( vars[k] )
                {
                    auxp[nzl] = k;
                    nzl++;
                }
            }
        }
        
    }
    
    return 0;
}
    
    


unsigned int OPT_MyNLPSolver::setFullHessIndex(const int mquad, const int *quadIndex, unsigned int* hessRowStart, int** hessCols)
{
    int n, r;
    unsigned int nz = 0, nzl;
    
    getNumberOfVars(n);
    
    for(int i = 0; i < n; i++)
    {
        hessRowStart[i] = nz;
        
        r = setHessIndexRow( i, mquad, quadIndex, hessCols[i], nzl );
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            return r;
        }
        
        nz += nzl;
    }
    
    
    hessRowStart[n] = nz;
    
    return 0;
}











