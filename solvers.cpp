
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cmath>
#include <new>


#include "OPT_solvers.hpp"
#include "MRQ_constants.hpp"
#include "MRQ_solvers.hpp"
#include "MRQ_tools.hpp"
#include "MRQ_bb.hpp"


using namespace std;
using namespace optsolvers;
using namespace muriqui;




int MRQ_Norm1ConstrHandler::addAuxVariables( const int n, const int naux, MRQ_LPSolver *master)
{
    int r;
    
    r = master->addVariables( naux, true );
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    for(int i = 0; i < naux; i++)
    {
        r += master->setObjLinearCoef(n+i, 1.0);
    }
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}


//we assume auxilary variables are added after the n original variables from minp problem.  Sol can be NULL. In this case, RHS will not be set
int MRQ_Norm1ConstrHandler::addNorm1constr(const int n, const int nIntVars, const int* intVars, const double *sol, MRQ_LPSolver *master, const bool roundSol, int &begNormConstrs)
{
    int code = 0, r, k;
    
    int cols[4];
    int rows[4];
    double coefs[4] = {1.0, -1.0, -1.0, -1.0};
    //double rhs[2];
    
    
    r = master->getNumberOfConstraints( begNormConstrs);
    #if MRQ_DEBUG_MODE
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    #endif
    
    r = master->addConstraints( 2*nIntVars );
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    for(k = 0; k < nIntVars; k++)
    {
        const int i = intVars[k];
        const int c = begNormConstrs + 2*k;
        
        cols[0] = cols[2] = i;
        cols[1] = cols[3] = n+k;
        
        rows[0] = rows[1] = c;
        rows[2] = rows[3] = c+1;
        
        
        r = master->setConstraintsLinearCoefs(4, rows, cols, coefs);
        
        if(r != 0)
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERROR;
            #endif
            code = MRQ_MILP_SOLVER_ERROR;
        }
    }
    
    //seting the rhs
    if(sol)
    {
        r = changeNorm1constr(begNormConstrs, nIntVars, intVars, sol, master, roundSol);
        if(r != 0)
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERROR;
            #endif
            code = MRQ_MILP_SOLVER_ERROR;
        }
    }
    
    
    return code;
}


int MRQ_Norm1ConstrHandler::changeNorm1constr( const int begNormConstrs,  const int nIntVars, const int* intVars, const double *sol, MRQ_LPSolver *master, const bool roundSol)
{
    int code = 0, k, r;
    //int indrows[2];
    //double rhs[2];
    
    
    #if MRQ_DEBUG_MODE
        if( begNormConstrs < 0 )
        {
            MRQ_PRINTERRORMSG("Error: No norm1constraint added. Call method addNorm1constr before that!");
            getchar();
            
            return MRQ_UNDEFINED_ERROR;
        }
    #endif
    
    for(k = 0; k < nIntVars; k++)
    {
        const int i = intVars[k];
        const double rsol = roundSol ? round(sol[i]) : sol[i];
        const double rhs0 = rsol;
        const double rhs1 = -rsol;
        const int c = begNormConstrs + 2*k;
        
        r = master->setConstraintBounds( c, -OPT_INFINITY, rhs0 );
        #if MRQ_DEBUG_MODE
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        #endif
        r = master->setConstraintBounds( c+1, -OPT_INFINITY, rhs1 );
        #if MRQ_DEBUG_MODE
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        #endif
    }
    
    
    return code;
}



MRQ_MasterMILPProb::MRQ_MasterMILPProb()
{
    initialize();
}


MRQ_MasterMILPProb::~MRQ_MasterMILPProb()
{
    desallocateMemory();
}


void MRQ_MasterMILPProb::desallocateMemory()
{
    MRQ_secFree(auxConstrEval);
    MRQ_secFree(auxCols);
    MRQ_secFree(auxVals);
    MRQ_secFree(auxConstr);
    
    MRQ_secDelete(master);
    MRQ_secDelete(gradEval);
}


void MRQ_MasterMILPProb::initialize()
{
    begNormConstrs = -1;
    
    auxConstrEval = NULL;
    auxCols = NULL;
    auxVals = NULL;
    auxConstr = NULL;
    master = NULL;
    gradEval = NULL;
}



int MRQ_MasterMILPProb::allocateMemory(const int sizeInds, const int sizeConstrs, const bool allocateGradEval, MRQ_MINLPProb *prob)
{
    if( sizeInds > 0 )
    {
        MRQ_malloc(auxCols, sizeInds); //auxCols = (int *) malloc( sizeInds * sizeof(int) );
        MRQ_malloc(auxVals, sizeInds); //auxVals = (double *) malloc( sizeInds * sizeof(double) );
        MRQ_IFMEMERRORRETURN( !auxCols || !auxVals);
    }
    
    if( sizeConstrs > 0 )
    {
        auxConstrEval = (bool *) malloc( sizeConstrs * sizeof(bool) );
        auxConstr = (double *) malloc( sizeConstrs * sizeof(double) );
        MRQ_IFMEMERRORRETURN( !auxConstrEval || !auxConstr);
    }
    
    if( allocateGradEval )
    {
        gradEval = new (std::nothrow) MRQ_GradientsEvaluation;
        MRQ_IFMEMERRORRETURN( !gradEval);
        
        const int r = gradEval->initialize(thnumber, prob);
        MRQ_IFERRORRETURN(r, r);
    }
    
    return 0;
}




//set The variables, linear constraints, variables bounds, variable types and linear part of objective function...
int MRQ_MasterMILPProb::setProblemBase(const int thnumber, MRQ_MINLPProb &prob, const int solver, const bool setLinearObj, const bool setQuad, const bool setVarTypes, const double *lx, const double *ux, const int nauxvars, MRQ_GeneralSolverParams *params, const int nthreads )
{
    const int n = prob.n;
    const int m = prob.m;
    const int nzq = prob.Q.getNumberOfElements();
    const double of = prob.objFactor;
    
    const bool *nlConstr = prob.nlConstr;
    const int *xtype = prob.xtype;
    const double *lc = prob.lc;
    const double *uc = prob.uc;
    
    
    
    const bool linearObj = !prob.hasNlObj && ( prob.Q.getNumberOfElements() == 0 || setQuad ) ;
    
    int k, ml, mq, r;
    int mmaster, code;
    
    //int *rows = NULL;
    int *cols = NULL;
    double *values = NULL;
    double* qvalues = NULL;
    
    minlpproblem::MIP_SparseMatrix *QC = prob.QC;
    minlpproblem::MIP_SparseMatrix &A = prob.A;
    
    
    this->thnumber = thnumber;
    this->prob = &prob;
    
    
    const int size = n+1;
    
    MRQ_malloc(cols, size); //cols = (int *) malloc( size* sizeof(int) );
    MRQ_malloc(values, size);// values = (double *) malloc( size* sizeof(double) ); 
    master = OPT_newLPSolver(solver);
    r = allocateMemory( n +nauxvars, m, true, &prob );
    MRQ_IFMEMERRORGOTOLABEL(!cols || !values || !master || r, code, termination);
    
    prob.getConstraintStatistcs(&ml, &mq, &r);
    
    
    mmaster = ml + (setQuad ? mq : 0) + (linearObj && setLinearObj ? 1 : 0);
    
    if(setQuad)
    {
        if( mq > 0 && master->getSolverType() < optsolvers::OPT_QCP )
        {
            MRQ_PRINTERRORMSG("User parameter enforce quadratic setting, but solver does not support quadratic constraints!");
            
            code = MRQ_MILP_SOLVER_ERROR;
            goto termination;
        }
    }
    
    
    
    r = master->initSolverEnv(0, n+nauxvars);
    MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    
    //master->setNumberOfThreads(1);
    
    if(params)
    {
        r = master->setParameters(*params);
    
        if( r != 0 )
        {
            //if(in_print_level > 0)
                MRQ_PRINTERRORMSG("Warning: Error on setting parameters list.");
        }
    }
    
    r = master->setNumberOfThreads(nthreads); //if we got an error, we do nothing
    if(r != 0)
        MRQ_PRINTERRORNUMBER(r);
    
    
    r = master->addVariables( n+nauxvars ); // if the problem has a nonlinear objective, we need an auxiliary variable...
    MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    
    
    r = master->addConstraints( mmaster );
    MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    
    
    for( int i = 0; i < n; i++ )
    {
        r = master->setVariableBounds( i, lx[i] <= -MIP_INFINITY ? -OPT_INFINITY : lx[i], ux[i] >= MIP_INFINITY ? OPT_INFINITY : ux[i] );
        MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    }
    
    for( int i = 0; i < nauxvars; i++)
    {
        //bounds to auxiliary variable
        r = master->setVariableBounds( n+i, -OPT_INFINITY, OPT_INFINITY );
        MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    }
    
    
    if( setVarTypes )
    {
        for( int i = 0; i < n; i++ )
        {
            if( minlpproblem::MIP_isIntegerType( xtype[i] ) )
            {
                r = master->setVariableType(i, optsolvers::OPT_VT_INTEGER);
                MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
            }
        }
    }
    
    
    
    k = 0;
    for(int i = 0; i < m; i++)
    {
        const int nz = A.getNumberOfElementsAtRow(i);
        
        const int nzqi = QC[i].getNumberOfElements();
        
        if( nlConstr[i] == false)
        {
            bool setcons;
            
            if( setQuad )
                setcons = nzqi > 0 || nz > 0;
            else
                setcons = nzqi == 0 && nz > 0;
            
            
            if( setcons )
            {
                //A[i].getStructureAndValues(cols, values);
                //r = master->resetLinearConstraintPart(k, nz, cols, values);
                
                r = master->resetConstraintLinearPart(k, nz, A[i], A(i));
                MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
                
                r = master->setConstraintBounds(k, lc[i] > -MIP_INFINITY ? lc[i] : -OPT_INFINITY,  uc[i] < MIP_INFINITY ? uc[i] : OPT_INFINITY );
                MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
                
                if( setQuad && nzqi > 0 )
                {
                    OPT_QCPSolver *qcpip = (OPT_QCPSolver*) master;
                    
                    //QC[i].getStructureAndValues( rows, cols, values );
                    //r = qcpip->setQuadConstraintMatrix( k, nzqi, rows, cols, values );
                    
                    r = qcpip->setConstraintQuadMatrix( k, (int*) QC[i].offset, QC[i][0],  QC[i](0) );
                    MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
                }
                
                k++;
            }
        }
        
    }
    
    
    //setting objective function. We use alpha variable because after is easier set objective cut just setting the upper bound of alpha. So, even if the objective is linear, we set alpha and add a constraint. We have little problems if we set quadratic in objective, but ok...
    
    master->setObjSense( optsolvers::OPT_MINIMIZE );
    
    
    for(int i = n; i < n+nauxvars; i++)
        r += master->setObjLinearCoef(i, 1.0);
    MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    
    
    if( linearObj && setLinearObj )
    {
        //if the problem has a linear objective function, we still use the alpha auxiliary variable...
        
        k = 1;
        cols[0] = n;
        values[0] = -1.0;
        double rhs = -prob.getObjConstant();
        
        
        if( prob.hasLinCoefObj() )
        {
            const double *c = prob.c;
            
            #if MRQ_DEBUG_MODE
                //we must have at least one auxiliary variable
                assert(nauxvars > 0);
            #endif
            
            
            
            for( int i = 0; i < n; i++ )
            {
                if( c[i] != 0.0 )
                {
                    cols[k] = i;
                    values[k] = c[i];
                    k++;
                }
            }
        }

            
        if( of != 1.0 )
        {
            rhs *= of;
            //for( int i = 1; i < k; i++ )
                //values[i] *= of;
            MRQ_multiplyAllArray(k-1, &values[1], of );
        }
        
        r = master->resetConstraintLinearPart( mmaster-1, k, cols, values );
        MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
        
        r = master->setConstraintBounds( mmaster-1, -OPT_INFINITY, rhs );
        MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    
        
    }
    
    
    if(setQuad && nzq > 0)
    {
        const minlpproblem::MIP_SparseMatrix &objQ = prob.Q;
        OPT_QCPSolver *qcpip = (OPT_QCPSolver*) master;
        
        if( master->getSolverType() < optsolvers::OPT_QP )
        {
            MRQ_PRINTERRORMSG("User parameter enforce quadratic setting, but solver does not support quadratic objective!");
            
            code = MRQ_MILP_SOLVER_ERROR;
            goto termination;
        }
        
        double* pvalues;
        
        //prob.Q.getStructureAndValues(rows, cols, values);
        
        if( of != 1.0 )
        {
            MRQ_malloc(qvalues, nzq); //qvalues = (double *) malloc( nzq* sizeof(double) );
            MRQ_IFMEMERRORGOTOLABEL(!qvalues, code, termination);
            
            objQ.getValues(qvalues);
            
            for(int i = 0; i < nzq; i++)
                qvalues[i] *= of;
            
            pvalues = qvalues;
        }
        else
        {
            pvalues = objQ(0);
        }
        
        r = qcpip->setObjQuadMatrix( (int *) objQ.offset,  objQ[0], pvalues );
        MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    }
    
    
    code = 0;
    
termination:
    
    //if(rows)	free(rows);
    if(cols)	free(cols);
    if(values)	free(values);
    if(qvalues) free(qvalues);
    
    return code;
}



/*
* Calculates binary cut to linearisation based algorithms. cols and values describes columns indices and coeficients of this constraints. nz has the number of nonzeros in cols and values. This cosntraints is less than rhs.
* 
* We assume all non-fixed variables are binaries. The number of elements set in values is nI.
*/
void muriqui::MRQ_calculateBinCut(const int nI, const int *intVars, const double *sol, double *values, double &rhs)
{
    int n1 = 0;
    
    for(int j = 0; j < nI; j++)
    {
        const int i = intVars[j];
        
        if(MRQ_abs(sol[i] - 1.0) <= MRQ_BINTOL)
        {
            //cols[j] = i;
            values[j] = 1.0;
            n1++;
        }
        else if(MRQ_abs(sol[i]) <= MRQ_BINTOL)
        {
            //cols[j] = i;
            values[j] = -1.0;
        }
        else
        {
            //we can have general integer variable fixed in the original problem...
            //cols[j] = i;
            values[j] = 0.0;
        }
        
    }
    
    #if MRQ_DEBUG_MODE
        assert(n1 <= nI);
    #endif
    
    rhs = n1-1;
}



int MRQ_MasterMILPProb::addBinaryCut(const int nI, const int *intVars, const double *x, double *auxVals)
{
    double rhs;
    int r, ind;
    
    MRQ_calculateBinCut(nI, intVars, x, auxVals, rhs);
    
    r = master->getNumberOfConstraints(ind);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    r = master->addConstraints(1);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    r = master->resetConstraintLinearPart(ind, nI, intVars, auxVals);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    r = master->setConstraintBounds(ind, -OPT_INFINITY, rhs);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}




int MRQ_MasterMILPProb::addNorm1constr(const int n, const int nIntVars, const int* intVars, const double *sol)
{
    MRQ_Norm1ConstrHandler norm1Handler;
    
    return norm1Handler.addNorm1constr(n, nIntVars, intVars, sol, master, false, begNormConstrs);
}




int MRQ_MasterMILPProb::changeNorm1constr( const int nIntVars, const int* intVars, const double *sol)
{
    MRQ_Norm1ConstrHandler norm1Handler;
    
    return norm1Handler.changeNorm1constr( begNormConstrs, nIntVars, intVars, sol, master, false);
}


//auxConstrEval2, auxCols, auxVals and constrValues are optional output parameters. If you do not pass, internal structures will be used.
int MRQ_MasterMILPProb::addConstraintLinearisationPointsToMILP( const double epsActiveConstrToLinearisation, unsigned int *nConstrLinearsSaved, const int nPoints, double **points, const bool incQuadsInMaster, const int constrLinStrat, bool *constrEval, bool *auxConstrEval2, int *auxCols, double *auxVals, double *constrValues)
{
    //const bool &hasNlObj = prob.hasNlObj;
    //const bool setObj = hasNlObj || ( prob.Q.getNumberOfElements() > 0 && !incQuadsInMaster );
    
    //const int &n = prob.n;
    const int m = prob->m;
    //bool newx;
    int code = 0;
    
    
    if( auxCols == NULL)
        auxCols = this->auxCols;
    
    if( auxVals == NULL )
        auxVals = this->auxVals;
    
    
    if( auxConstrEval2 == NULL )
        auxConstrEval2 = this->auxConstrEval;
    
    if( constrValues == NULL )
        constrValues = this->auxConstr;
    
    
    if( nPoints > 0 )
    {// in the first point, we linearize all constraints (we have to guaratee all constraints are being approximated by at least one linearization...)
        
        for(int j = 0; j < m; j++)
            constrEval[j] =  prob->nlConstr[j] || (prob->QC[j].getNumberOfElements() > 0 && !incQuadsInMaster);
        
        int r = addLinearizedNLConstraints( true, points[0], incQuadsInMaster, constrEval, auxCols, auxVals, constrValues, false);
        MRQ_IFERRORRETURN(r, r);
    }
    
    
    
    //first, we linearize the constraints.
    for(int k = 1; k < nPoints; k++)
    {
        int r = addLinearizedNLConstraintsByStrategy( epsActiveConstrToLinearisation, nConstrLinearsSaved, true, points[k], incQuadsInMaster, constrLinStrat, constrEval, auxConstrEval2, auxCols, auxVals, constrValues, false );
        
        if( r != 0 )
        {
            //if( in_print_level > 0 )
                MRQ_PRINTERRORNUMBER(r);
            
            code = r;
        }
    }
    
    
    return code;
}






int MRQ_MasterMILPProb::addObjLinearisationPointsToMILP( const int nPoints, double **points, const double zu, const bool incQuadsInMaster, const int objLinStrat, int *auxCols, double *auxVals, MRQ_LAAPointsStoring *laps)
{
    const bool &hasNlObj = prob->hasNlObj;
    const bool setObj = hasNlObj || ( prob->Q.getNumberOfElements() > 0 && !incQuadsInMaster );
    
    const int &n = prob->n;
    //const int &m = prob.m;
    //bool newx;
    int code = 0, r;
    int mmaster;
    
    
    if( auxCols == NULL )
        auxCols = this->auxCols;
    
    if( auxVals == NULL )
        auxVals = this->auxVals;
    
    
    
    if( setObj )
    {
        if( objLinStrat == MRQ_OLS_NON_OBJ_CUT_POINTS )
        {
            r = master->getNumberOfConstraints( mmaster);
            #if MRQ_DEBUG_MODE 
                if(r != 0)
                {
                    //if(in_print_level > 0)
                        MRQ_PRINTERRORNUMBER(r);
                    
                    code = MRQ_MILP_SOLVER_ERROR;
                }
            #endif
        }
        
        
        for(int i = 0; i < nPoints; i++)
        {
            //newx = true;
            r = addLinearizedObjFunction( true, points[i], incQuadsInMaster, auxCols, auxVals, NULL);
            MRQ_IFERRORSETVAR(r, code, r);
            
            //if( hasNlObj )
                //newx = false;
            
            if( objLinStrat == MRQ_OLS_NON_OBJ_CUT_POINTS && laps)
            {
                //auxVars has the gradient of obj linearization, i.e., the coefiicents of obj linearizations constraint and the rhs at the end...
                r = laps->updateObjLinearizationsByNonObjCutPointsByNewPoint( *master, zu, auxVals );
                MRQ_IFERRORSETVAR(r, code, r);
                
                r = laps->addPoint(n, points[i]);
                MRQ_IFERRORSETVAR(r, code, r);
                
                laps->indMaster[ laps->npoints-1 ] = mmaster;
                
                mmaster++;
            }
        }
    }
    
    return code;
}


/*Here, it is implicit colas and values has prob.n+1 elements*/
int muriqui::MRQ_calculateObjLinearizedConstraint(MRQ_MINLPProb &prob, const int thnumber, bool newx, const double* x, const bool incQuadsInMaster, const double* objValue, int* cols, double* values, double &rhs)
{
    const int n = prob.n;
    int r, code;
    double aux;
    
    rhs = 0.0;
    
    if( !prob.hasNlObj )
    {
        if( prob.Q.getNumberOfElements() == 0 || incQuadsInMaster )
        {
            code = 0;
            goto termination;
        }
    }
    
    
    if( incQuadsInMaster )
    {
        //we must not consider the quadratic part in the gradient...
        if( prob.hasNlObj )
        {
            r = prob.nlEval->eval_nl_obj_part(thnumber, n, newx, x, rhs);
            MRQ_IFCALLBACKERRORGOTOLABEL(r, code, termination);
            
            rhs = -rhs -prob.d; //considering objective constant...
            
            newx = false;
            
            r = prob.nlEval->eval_grad_nl_obj_part( thnumber, n, newx, x, values);
            MRQ_IFCALLBACKERRORGOTOLABEL(r, code, termination);
        }
        else
        {
            MRQ_setAllArray(n, values, 0.0);
        }
        
        
        /*if(!incQuadsInMaster && prob.Q.getNumberOfElements() > 0) 
        {
            // accumlating the gradient of 0.5*xQx in values
            
            rhs += prob.Q.quadraticEvaluationAndGradient(x, auxVals, true);
        } */
        
        if( prob.hasLinCoefObj() )
        {
            const double* c = prob.c;
            #pragma ivdep
            #pragma GCC ivdep
            for(int i = 0; i < n; i++)
                values[i] += c[i];
        }
        
        
        if( prob.objFactor != 1.0)
        {
            const double of = prob.objFactor;
            
            //for(int i = 0; i < n; i++)
                //values[i] *= of;
            
            MRQ_multiplyAllArray(n, values, of);
            
            rhs *= of;
        }
        
    }
    else
    {
        if( objValue )
        {
            rhs = -(*objValue);
        }
        else
        {
            r =prob.objEval( thnumber, newx, x, rhs);
            MRQ_IFCALLBACKERRORGOTOLABEL(r, code, termination);
            
            rhs = -rhs;
            
            if( prob.hasNlObj )
                newx = false;
        }
        
        r = prob.objGradEval(thnumber, newx, x, values);
        MRQ_IFCALLBACKERRORGOTOLABEL(r, code, termination);
        
    }
    
    
    aux = 0.0;
    //we do that account in a separted variable to try avoid numerical problems...
    for(int i = 0; i < n; i++)
        aux += values[i]*x[i];
    
    rhs += aux;
    
    
    if(cols)
    {
        for(int i = 0; i <= n; i++) //considering alpha variable also
            cols[i] = i;
    }
    
    values[n] = -1.0; //alpha variable
    
    
    code = 0;
    
termination:
    
    return code;
}


int MRQ_MasterMILPProb::addLinearizedObjFunction( bool newx, const double* x, const bool incQuadsInMaster, const int in_obj_linearisation_strategy, const double zu, MRQ_LAAPointsStoring *laps, int* auxCols, double* auxVals, double* objValue)
{
    const int n = prob->n;
    
    int r = addLinearizedObjFunction(newx, x, incQuadsInMaster, auxCols, auxVals, objValue);
    MRQ_IFERRORRETURN(r, r);
    
    
    if( in_obj_linearisation_strategy == MRQ_OLS_NON_OBJ_CUT_POINTS )
    {
        int mmaster;
        
        int r = laps->updateObjLinearizationsByNonObjCutPointsByNewPoint( *master, zu, auxVals );
        MRQ_IFERRORRETURN(r, r);
        
        r = laps->addPoint(n,  x);
        MRQ_IFERRORRETURN(r, r);
        
        r = master->getNumberOfConstraints( mmaster);
        MRQ_IFERRORRETURN(r, r);
        
        laps->indMaster[laps->npoints-1] = mmaster;
    }
    
    
    return 0;
}


int MRQ_MasterMILPProb::addLinearizedObjFunction( bool newx, const double* x, const bool incQuadsInMaster, int* auxCols, double* auxVals, double* objValue)
{
    const int n = prob->n;
    int r, mmas, code;
    double rhs;
    
    
    if( auxCols == NULL )
        auxCols = this->auxCols;
    
    if( auxVals == NULL )
        auxVals = this->auxVals;
    
    
    if( !prob->hasNlObj )
    {
        if( prob->Q.getNumberOfElements() == 0 || incQuadsInMaster )
        {
            code = 0;
            goto termination;
        }
    }
    
    
    r = MRQ_calculateObjLinearizedConstraint(*prob, thnumber, newx, x, incQuadsInMaster, objValue, auxCols, auxVals, rhs);
    MRQ_IFERRORGOTOLABEL(r, code, r, termination);
    
    r = master->getNumberOfConstraints(mmas);
    MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    
    r += master->addConstraints(1);
    MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    
    r = master->resetConstraintLinearPart( mmas, n+1, auxCols, auxVals );
    MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    
    r = master->setConstraintBounds( mmas, -OPT_INFINITY, rhs );
    MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    
    //we just put rhs in the last position of auxVals because we need store the linearization after in the algorithms
    auxVals[n] = rhs;
    
    
    //linear part of objective function is already considered in the objective function os master problem, including obj constant...
    
    
    code = 0;
    
termination:
    
    return code;
}


int MRQ_MasterMILPProb::addLinearizedNLConstraints( bool newx, const double* x, const bool incQuadsInMaster, const bool* constrEval, int* auxCols, double* auxVals, double* constrValues, bool constrCalculated)
{
    const int n = prob->n;
    const int m = prob->m;
    int k, r, code;
    int indstart, nnewc = 0;
    double rhs, nuc, nlc;
    
    const bool *nlConstr = prob->nlConstr;
    const double *lc = prob->lc, *uc = prob->uc;
    MRQ_SparseMatrix *QC = prob->QC;
    
    
    
    if(auxCols == NULL)
        auxCols = this->auxCols;
    
    if(auxVals == NULL)
        auxVals = this->auxVals;
    
    if(constrValues == NULL)
        constrValues = this->auxConstr;
    
    
    for(int i = 0; i < m; i++)
    {
        if( constrEval[i] &&  ( nlConstr[i] || (QC[i].getNumberOfElements() > 0 && !incQuadsInMaster ) )  )
            nnewc++;
    }
    
    if(nnewc == 0)
    {
        code = 0;
        goto termination;
    }
    
    r = master->getNumberOfConstraints(indstart);
    MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    
    r = master->addConstraints(nnewc);
    MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    
    for(int i = 0; i < n; i++)
        auxCols[i] = i;
    
    if( !constrCalculated )
    {
        r = prob->constraintsEval(thnumber, newx,  constrEval, x, constrValues);
        MRQ_IFCALLBACKERRORGOTOLABEL(r, code, termination);
        
        if( prob->hasNlConstrs )
            newx = false;
    }
    
    
    
    if( prob->hasNlConstrs )
    {
        r = gradEval->evaluateJacobian(newx, constrEval, x);
        MRQ_IFCALLBACKERRORGOTOLABEL(r, code, termination);
        
        newx = false;
    }
    
    
    
    k = 0;
    for(int i = 0; i < m; i++)
    {
        if( constrEval[i] == false )
        {
            //if( nlConstr[i] || (QC[i].getNumberOfElements() > 0u && !incQuadsInMaster) )
                //std::cout << "Descartando linearização sobre restrição " << i << std::endl;
            
            continue;
        }
        
        if( nlConstr[i] == false && ( QC[i].getNumberOfElements() == 0u || incQuadsInMaster ) )
            continue;
        
        //let's linearize
        
        gradEval->constraintCompleteGradient(i, x, auxVals);
        
        
        r = master->resetConstraintLinearPart( indstart + k, n, auxCols, auxVals );
        MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
        
        
        rhs = MRQ_calculateDeltaRHSToConstraintLinearisation(constrValues[i], n, auxVals, x);
        
        //cout << "constrValues[" << i << "]: " << constrValues[i] << " rhs: " << rhs << " lc[" << i << "]: " << lc[i] << " uc[" << i << "]: " <<  uc[i] <<  endl;
        
        nuc = uc[i] < MIP_INFINITY ? rhs + uc[i] : OPT_INFINITY;
        
        nlc = lc[i] > -MIP_INFINITY ? rhs + lc[i] : -OPT_INFINITY;
        
        
        //printf("i: %d cValues[%d]: %0.14f nlc: %0.14f nuc:%0.14f\n", i, i, constrValues[i], nlc, nuc);
        
        
        r = master->setConstraintBounds(indstart + k, nlc, nuc );
        MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
        
        k++;
    }
    
    #if MRQ_DEBUG_MODE
        assert(k == nnewc);
    #endif
    
    code = 0;
    
termination:
    
    return code;
}



void muriqui::MRQ_calculateConstraintsToBeLinearizedByStrategy( MRQ_MINLPProb &prob, const int constrLinStrat, const double epsActiveConstrToLinearisation, const bool* constrEval, const double* constrValues, const double *masterConstrValues, const double *newlc, const double *newuc, unsigned int *nConstrLinearsSaved, bool* outputConstrEval)
{
    const int m = prob.m;
    unsigned int nConsSaved = 0;
    
    
    if( constrLinStrat == MRQ_CLS_ONLY_INFEAS_AND_ACTIVE || constrLinStrat == MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_MASTER_SOL_ALSO )
    {
        const double *lc = prob.lc, *uc = prob.uc;
        int nlc = 0;
        
        for(int i = 0; i < m; i++)
        {
            outputConstrEval[i] = false;
            
            if( constrEval[i] )
            {
                //we run two times. In the first we check the constriant values. in the second, we check the master constraint Values
                for(signed char j = 0; j < 2; j++)
                {
                    const double* pconstrValues = j == 0 ? constrValues : masterConstrValues;
                    
                    if( constrLinStrat != MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_MASTER_SOL_ALSO )
                        j = 10; //to stop the for after this iteration
                    
                    if( pconstrValues[i] < lc[i] || pconstrValues[i] > uc[i] )
                    { //constraint i is infeasible in this solution
                        outputConstrEval[i] = true;
                    }
                    else
                    {
                        //constraint i is feasible in this solution...
                        if(lc[i] > -MIP_INFINITY)
                        {
                            const double den = MRQ_max(MRQ_abs(lc[i]), 1.0); // it rhs is less than 1, we do not divide...
                            
                            const double f = (pconstrValues[i] - lc[i])/den;
                            
                            if( f <= epsActiveConstrToLinearisation )
                                outputConstrEval[i] = true;
                        }
                        
                        if(uc[i] < MIP_INFINITY)
                        {
                            const double den = MRQ_max(MRQ_abs(uc[i]), 1.0); // it rhs is less than 1, we do not divide...
                            
                            const double f = (uc[i] - pconstrValues[i])/den;
                            
                            if( f <= epsActiveConstrToLinearisation )
                                outputConstrEval[i] = true;
                        }
                    }
                    
                    if( outputConstrEval[i] )
                    {
                        nlc++;
                        j = 10; //to stop the for after this iteration
                    }
                    
                }
                
                if( !outputConstrEval[i] )
                    nConsSaved++;
            }
            
            //if( constrEval[i] && masterConstrValues )
                //std::cout << "i: " << i << " cEval: " << constrEval[i] << " lc: " << lc[i] << " uc: " << uc[i] << " cVal:" << constrValues[i] << " mastercVal: " << masterConstrValues[i] << " auxCEval2: " << auxConstrEval2[i] << "\n" ;
        }
        
        
        //if( nlc == 0 )
            //goto termination;
    }
    else if( constrLinStrat == MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_BY_BOX_FILL )
    {
        for(int i = 0; i < m; i++)
        {
            outputConstrEval[i] = false;
            
            if( constrEval[i] )
            {
                if( newlc[i] <= -MIP_INFINITY || newuc[i] >= MIP_INFINITY )
                {
                    outputConstrEval[i] = true;
                }
                else
                {
                    const double box = newuc[i] - newlc[i];
                    const double epszero = 1e-3;
                    
                    if( box <= epszero )
                        outputConstrEval[i] = true;
                    else
                    {
                        if( (newuc[i] - constrValues[i])/box < epsActiveConstrToLinearisation )
                            outputConstrEval[i] = true;
                    }
                }
                
                if( !outputConstrEval[i] )
                    nConsSaved++;
            }
        }
    }
    else
    {
        MRQ_copyArray(m, constrEval, outputConstrEval);
    }
    
    
    
    if( nConstrLinearsSaved )
        *nConstrLinearsSaved += nConsSaved;
}



int MRQ_MasterMILPProb::addLinearizedNLConstraintsByStrategy( const double epsActiveConstrToLinearisation, unsigned int *nConstrLinearsSaved, const bool newx, const double* x, const bool incQuadsInMaster, int constrLinStrat, const bool* constrEval, bool* auxConstrEval2, int* auxCols, double* auxVals, double* constrValues, const bool constrCalculated, const double *masterConstrValues, const double *newlc, const double *newuc)
{
    
    
    //const double EPS_DEM = 1.0e-5;
    
    int r, code = 0;
    
    
    if( auxConstrEval2 == NULL )
        auxConstrEval2 = this->auxConstrEval;
    
    if( auxCols == NULL )
        auxCols = this->auxCols;
    
    if( auxVals == NULL )
        auxVals = this->auxVals;
    
    if( constrValues == NULL )
        constrValues = this->auxConstr;
    
    const bool *pceval;
    
    if( !masterConstrValues && constrLinStrat == MRQ_CLS_ONLY_INFEAS_AND_ACTIVE_MASTER_SOL_ALSO )
        constrLinStrat = MRQ_CLS_ONLY_INFEAS_AND_ACTIVE;
    
    
    if( !constrCalculated )
    {
        int r = prob->constraintsEval( thnumber, newx, constrEval, x, constrValues );
        MRQ_IFCALLBACKERRORSETVAR(r, code);
    }
    
    /*for(int i = 0; i < m; i++)
    {
        if( prob->nlConstr[i] || ( prob->QC[i].getNumberOfElements() > 0u && !incQuadsInMaster) )
        {
            std::cout << "cValues["<<i<<"]: " << constrValues[i] << "   lc["<<i<<"]: " << ( prob->lc[i] <= -MIP_INFINITY ? -INFINITY : prob->lc[i] ) 
            
            << "   uc["<<i<<"]: " << ( prob->uc[i] >= MIP_INFINITY ? INFINITY : prob->uc[i] );
            
            if(newlc)
                std::cout << "   newlc: " << newlc[i] << "   newuc: " << newuc[i];
            
            std::cout << "   cEval["<<i<<"]: " << auxConstrEval2[i] << "\n";
        }
    }
    
    MRQ_getchar(); */
    
    //std::cout << "constrLinStrat: " << constrLinStrat << "\n";
    
    if(constrLinStrat == MRQ_CLS_ALL_CONSTRS)
    {
        pceval = constrEval;
    }
    else
    {
        MRQ_calculateConstraintsToBeLinearizedByStrategy(*prob, constrLinStrat, epsActiveConstrToLinearisation, constrEval, constrValues, masterConstrValues, newlc, newuc, nConstrLinearsSaved, auxConstrEval2);
        
        pceval = auxConstrEval2;
    }
    
    
    r = addLinearizedNLConstraints( newx, x, incQuadsInMaster, pceval, auxCols, auxVals, constrValues, true );
    MRQ_IFERRORSETVAR(r, code, r);
    
//termination:
    
    //if( nConstrLinearsSaved )
        //*nConstrLinearsSaved += nConsSaved;
    
    return code;
}





int MRQ_MasterMILPProb::setVarBounds(const int n, const double *lb, const double *ub)
{
    return MRQ_setVarBounds(master, n, lb, ub);
}



bool MRQ_MasterMILPProb::minlpProbhasNLConstraints()
{
    return prob->hasNLConstraints();
}


int MRQ_NlpFeasBuilder::setAuxVariables( optsolvers::OPT_LPSolver *solver, int *nAuxVars) 
{
    int code, r;
    int msolver = 0;
    int nsolver, nextaux, nauxvars, nz;
    int inds[2];
    double lci, uci;
    double coefs[2];
    
    r = solver->getNumberOfConstraints( msolver );
    r+= solver->getNumberOfVars( nsolver );
    
    #if MRQ_DEBUG_MODE
        assert( r == 0 );
    #endif
    
    //counting the number of auxiliary variables
    nauxvars = 0;
    for(int i = 0; i < msolver; i++)
    {
        r = solver->getConstraintBounds(i, lci, uci);
        
        #if MRQ_DEBUG_MODE
            assert( r == 0 );
        #endif
        
        //we need two auxiliary variables to double bound (maybe equalities) constraints. Ok, just linear constraints can be double bounded, but I do it anyway...
        
        if( lci > -OPT_INFINITY )
            nauxvars++;
        
        if( uci < OPT_INFINITY )
            nauxvars++;
    }
    
    r = solver->addVariables( nauxvars, false );
    MRQ_IFERRORGOTOLABEL(r, code, MRQ_NLP_SOLVER_ERROR, termination);
    
    
    if(nAuxVars)
        *nAuxVars = nauxvars;
    
    nextaux = nsolver;
    for( int i = 0; i < msolver; i++)
    {
        r = solver->getConstraintBounds(i, lci, uci);
        
        #if MRQ_DEBUG_MODE
            assert( r == 0 );
        #endif
        
        nz = 0;
        
        if( lci > -OPT_INFINITY )
        {
            coefs[0] = 1.0;
            inds[0] = nextaux;
            
            nextaux++;
            nz = 1;
        }
        
        
        if( uci < OPT_INFINITY )
        {
            coefs[nz] = -1.0;
            inds[nz] = nextaux;
            
            nextaux++;
            nz++;
        }
        
        r = solver->setConstraintLinearCoefs( i, nz, inds, coefs );
        MRQ_IFERRORGOTOLABEL(r, code, MRQ_NLP_SOLVER_ERROR, termination);
    }
    
    #if MRQ_DEBUG_MODE
        assert( nextaux -nsolver == nauxvars );
    #endif
    
    
    for( int i = nsolver; i < nextaux; i++ )
    {
        r = solver->setObjLinearCoef( i, 1.0 );
        r += solver->setVariableBounds( i, 0.0, OPT_INFINITY);
        
        #if MRQ_DEBUG_MODE
            assert( r == 0 );
        #endif
    }
    
    
    code = 0;
    
termination:
    
    return code;
}





MRQ_NLPFeasProb::MRQ_NLPFeasProb()
{
    initialize();
}


MRQ_NLPFeasProb::~MRQ_NLPFeasProb()
{
    desallocate();
}


void MRQ_NLPFeasProb::desallocate()
{
    MRQ_secDelete(eval);
    MRQ_secDelete(solver);
}


void MRQ_NLPFeasProb::initialize()
{
    eval = NULL;
    solver = NULL;
}


int MRQ_NLPFeasProb::setProblem(const int solverCode, MRQ_MINLPProb &prob, const double *lx, const double *ux, MRQ_GeneralSolverParams *params, const int thnumber, const bool setSpecialParams, const int nthreads, const double maxCpuTime, const double maxTime)
{
    int r;
    MRQ_NlpFeasBuilder fbuilder;
    
    solver = OPT_newLPSolver(solverCode);
    MRQ_IFMEMERRORRETURN(!solver);
    
    r = MRQ_setNLPRelaxProb(prob, lx, ux, NULL, NULL, solver, false, true, true, false, thnumber, setSpecialParams, params, nthreads, maxCpuTime, maxTime, 0, 0);
    MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    
    //r = solver->setnVariablesBounds( prob.n, lx, ux );
    //MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    
    solver->setObjSense( optsolvers::OPT_MINIMIZE ); // obj function is nondetermined in optsolvers
    
    
    r = fbuilder.setAuxVariables(solver);
    MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);	
    
    if( prob.hasNLConstraints() )
    {
        eval = new (std::nothrow) MRQ_NLPNonObjEval( prob.n, prob.m, prob.getNumberOfJacobianNonZeros(), prob.getNumberOfLagrangianHessianNonZeros(), prob.getNonLinearEvaluationObject() );
        MRQ_IFMEMERRORRETURN(!eval);
    }
    
    
    if( prob.hasNLConstraints() || prob.hasObjNLTerm() )
    //we could do it inside if above, bt we let here to test if there is some bug in OPT_SOLVERS...
        ( (OPT_NLPSolver *) solver )->setNonLinearEvalObject( eval );
    
    return 0;
}




MRQ_NLPIPProblem::MRQ_NLPIPProblem()
{
    initialize();
}


MRQ_NLPIPProblem::~MRQ_NLPIPProblem()
{
    desallocate();
}


void MRQ_NLPIPProblem::desallocate()
{
    MRQ_secDelete(eval);
    MRQ_secDelete(solver);
}


void MRQ_NLPIPProblem::initialize()
{
    eval = NULL;
    solver = NULL;
}


int MRQ_NLPIPProblem::setProblem(const int solverCode, MRQ_MINLPProb &prob, const double *lx, const double *ux, MRQ_GeneralSolverParams *params, const int thnumber, const bool setSpecialParams, const bool considerQuadToInterior, const int nthreads, const double maxCpuTime, const double maxTime)
{
    const int n = prob.n;
    const int m = prob.m;
    const int auxInd = n;
    const bool *nlConstr = prob.nlConstr;
    const double *lc = prob.lc, *uc = prob.uc;
    minlpproblem::MIP_SparseMatrix *QC = prob.QC;
    
    int r;
    
    
    solver = OPT_newLPSolver(solverCode);
    MRQ_IFMEMERRORRETURN(!solver);
    
    
    r = MRQ_setNLPRelaxProb(prob, lx, ux, NULL, NULL, solver, false, true, true, false, thnumber, setSpecialParams, params, nthreads, maxCpuTime, maxTime, 1, 0);
    MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    
    //r = solver->setnVariablesBounds( prob.n, lx, ux );
    //MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    
    solver->setObjSense( optsolvers::OPT_MINIMIZE ); // obj function is nondetermined in optsolvers
    
    
    //seting auxilary variable part. we assume nonlinear and quadratic constraints are not double bounded...
    for(int i = 0; i < m; i++)
    {
        if( lc[i] <= -MIP_INFINITY && uc[i] >= MIP_INFINITY )
            continue;
        
        
        if( nlConstr[i] || (QC[i].getNumberOfElements() > 0 && considerQuadToInterior) )
        {
            int r = 0;
            double auxcoef;
            
            if( lc[i] > -MIP_INFINITY )
            {
                auxcoef = MRQ_abs(lc[i]) + 1.0 ; //we try consider a scale to constraints...
                
                #if MRQ_DEBUG_MODE
                    assert( uc[i] >= MIP_INFINITY ); //we assume nonlinea and quadratic constraints are not double bounded here... if some day you need, you have to change it... :(
                #endif
            }
            else if( uc[i] < MIP_INFINITY )
            {
                auxcoef = -MRQ_abs(uc[i]) - 1.0; //we try consider a scale to constraints...
            }
            else
            {
                //if constraint is free, we do not need add variables. anyway, we set it to compiler do not complaina bout uninitialized value. Actually, program should never pass here.
                auxcoef = 0.0;
            }
            
            
            r = solver->setConstraintLinearCoef(i, auxInd, auxcoef);
            MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        }
        
    }
    
    solver->setObjSense(optsolvers::OPT_MINIMIZE);
    r = solver->setObjLinearCoef(auxInd, 1.0);
    MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    
    
    // end of setting auxialiary variable part
    
    
    
    if( prob.hasNLConstraints() )
    {
        eval = new (std::nothrow) MRQ_NLPNonObjEval(n, m, prob.getNumberOfJacobianNonZeros(), prob.getNumberOfLagrangianHessianNonZeros(), prob.getNonLinearEvaluationObject() );
        MRQ_IFMEMERRORRETURN(!eval)
    }
    
    
    if( prob.hasNLConstraints() || prob.hasObjNLTerm() )
    //we could do it inside if above, bt we let here to test if there is some bug in OPT_SOLVERS...
        ( (OPT_NLPSolver *) solver )->setNonLinearEvalObject( eval );
    
    return 0;
}



MRQ_NLPFeasPumpProb::MRQ_NLPFeasPumpProb()
{
    initialize();
}


MRQ_NLPFeasPumpProb::~MRQ_NLPFeasPumpProb()
{
    desallocate();
}


int MRQ_NLPFeasPumpProb::allocaten1varIndex(const int size)
{
    n1varIndex = (int *) malloc(size *sizeof(int));
    MRQ_IFMEMERRORRETURN( !n1varIndex );
    
    return 0;
}


void MRQ_NLPFeasPumpProb::desallocate()
{
    MRQ_secFree(n1varIndex);
    nn1vars = 0;
    MRQ_secDelete(eval);
    MRQ_secDelete(solver);
}


void MRQ_NLPFeasPumpProb::initialize()
{
    nn1vars = 0;
    n1varIndex = NULL;
    eval = NULL;
    solver = NULL;
}



int MRQ_NLPFeasPumpProb::setProblemBase(const int solverCode, MRQ_MINLPProb &prob, const double *lx, const double *ux, const double *lc, const double *uc, MRQ_GeneralSolverParams *params, const int thnumber, const bool setSpecialParams, const bool setLinearObjTermOnBinVars, const bool useNorm1, const int nthreads, const double maxCpuTime, const double maxTime)
{
    const int n = prob.n;
    const int m = prob.m;
    const int *xtype = prob.xtype;
    int r;
    
    this->useNorm1 = useNorm1;
    this->setLinearObjTermOnBinVars = setLinearObjTermOnBinVars;
    
    
    solver = OPT_newLPSolver(solverCode);
    MRQ_IFMEMERRORRETURN( !solver );
    
    
    r = MRQ_setNLPRelaxProb( prob, lx, ux, lc, uc, solver, false, true, true, false, thnumber, setSpecialParams, params, nthreads, maxCpuTime, maxTime, 0, 0 );
    MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    
    /*r = solver->setnVariablesBounds(n, lx, ux);
    MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    
    if( lc != NULL)
    {
        r = 0;
        for(int i = 0; i < m; i++)
            r += solver->setConstraintBounds(i, lc[i], uc[i]);
        
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    }*/
    
    
    if( prob.hasNLConstraints() )
    {
        eval = new (std::nothrow) MRQ_NLPNonObjEval(n, m, prob.getNumberOfJacobianNonZeros(), prob.getNumberOfLagrangianHessianNonZeros(), prob.getNonLinearEvaluationObject() );
        MRQ_IFMEMERRORRETURN( !eval );
    }
    
    
    if( prob.hasNLConstraints() || prob.hasObjNLTerm() )
    //we could do it inside if above, bt we let here to test if there is some bug in OPT_SOLVERS...
        ( (OPT_NLPSolver *) solver )->setNonLinearEvalObject( eval );
    
    
    //seting objective function
    solver->setObjSense( optsolvers::OPT_MINIMIZE ); // obj function is nondetermined in optsolvers
    
    
    if( setLinearObjTermOnBinVars )
    {
        nn1vars = 0;
        for(int i = 0; i < n; i++)
        {
            if( minlpproblem::MIP_isIntegerType( xtype[i] ) )
            {
                if( lx[i] != ux[i] && !MRQ_isBinary(lx[i], ux[i]) )
                {
                    nn1vars++;
                }
            }
        }
        
        if( nn1vars )
        {
            r = allocaten1varIndex(nn1vars);
            MRQ_IFERRORRETURN(r, r);
            
            int k = 0;
            for(int i = 0; i < n; i++)
            {
                if( minlpproblem::MIP_isIntegerType( xtype[i] ) )
                {
                    //we ignore fixed variables beyond binaries
                    if( ux[i] != lx[i] && !MRQ_isBinary(lx[i], ux[i]) )
                    {
                        n1varIndex[k] = i;
                        k++;
                    }
                }
            }
            
            #if MRQ_DEBUG_MODE
                assert( k == nn1vars );
            #endif
        }
        
    }
    else
    {
        nn1vars = 0;
        for(int i = 0; i < n; i++)
        {
            if( minlpproblem::MIP_isIntegerType( xtype[i] ) && ux[i] != lx[i] )
                nn1vars++;
        }
        
        r = allocaten1varIndex( nn1vars );
        MRQ_IFERRORRETURN(r, r);
        
        int k = 0;
        for(int i = 0; i < n; i++)
        {
            if( minlpproblem::MIP_isIntegerType( xtype[i] ) && ux[i] != lx[i] )
            {
                n1varIndex[k] = i;
                k++;
            }
        }
        
        #if MRQ_DEBUG_MODE
            assert( nn1vars == k );
        #endif
    }
    
    
    if( nn1vars > 0 && useNorm1 )
    {
        int begAuxVars;
        MRQ_Norm1ConstrHandler n1hand;
        
        //we just get the number of vars because gurobi is a disaster and add auxiliary variables to range constraints...
        r = solver->getNumberOfVars( begAuxVars );
        #if MRQ_DEBUG_MODE
            MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        #endif
        
        r = n1hand.addAuxVariables( begAuxVars, nn1vars, solver );
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        r = n1hand.addNorm1constr(begAuxVars, nn1vars, n1varIndex, NULL, solver, true, begNorm1Constrs);
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    }
    
    
    return 0;
}




//if setLinearObjTermOnBinVars is true, we set the terms in objective function related to bin vars as linear terms. Otherwise, we use a quadratic terms when setLinearObjTermOnBinVars is false, or for non binary variables
int MRQ_NLPFeasPumpProb::setObjective( const int nI, const int *intVars, const double *lx, const double *ux, const double *sol )
{
    int r, code = 0;
    
    if( setLinearObjTermOnBinVars )
    {
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            if( lx[ind] == ux[ind] )
                continue;
            
            if( MRQ_isBinary(lx[ind], ux[ind])  )
            {
                r = solver->setObjLinearCoef(ind, (sol[ind] < 0.5 ? 1.0 : -1.0) );
                MRQ_IFERRORSETVAR(r, code, MRQ_NLP_SOLVER_ERROR);
            }
            /*else
            {//we use a quadratic term...
                const double solind = round(sol[ind]);
                
                //we do not set the constant term...
                r = solver->setLinearObjCoef(ind, -2*solind);
                
                r += ((OPT_QPSolver*) solver)->setQuadObjCoef(ind, ind, 2.0);
            } */
            
        }
    }
    
    if( nn1vars > 0 )
    {
        if( useNorm1 )
        {
            MRQ_Norm1ConstrHandler n1hand;
            
            r = n1hand.changeNorm1constr( begNorm1Constrs, nn1vars, n1varIndex, sol, solver, true );
            MRQ_IFERRORSETVAR(r, code, MRQ_NLP_SOLVER_ERROR);
        }
        else
        {//we add a quadratic term to simulate a norm2
            
            for(int i = 0; i < nn1vars; i++)
            {
                const int ind = n1varIndex[i];
                
                const double solind = round(sol[ind]);
                
                //we do not set the constant term...
                r = solver->setObjLinearCoef(ind, -2*solind);
                
                r += ((OPT_QPSolver*) solver)->setObjQuadCoef(ind, ind, 2.0);
                
                MRQ_IFERRORSETVAR(r, code, MRQ_NLP_SOLVER_ERROR);
            }
            
        }
    }
    
    
    return code;
}




MRQ_NoObjWithObjCutNLPEval::MRQ_NoObjWithObjCutNLPEval( const bool nlConstrs, minlpproblem::MIP_NonLinearEval *originalEval, int origm, double objFactor)
{
    initialize(nlConstrs, originalEval, origm, objFactor);
}


MRQ_NoObjWithObjCutNLPEval::~MRQ_NoObjWithObjCutNLPEval()
{
}


void MRQ_NoObjWithObjCutNLPEval::initialize( const bool nlConstrs, minlpproblem::MIP_NonLinearEval *originalEval, int origm, double objFactor)
{
    mynewx = true;
    this->objFactor = objFactor;
    //nlObjCut = false;
    this->origm = origm;
    this->nlConstrs = nlConstrs;
    oeval = originalEval;
    
    //auxVals = NULL;
}


int MRQ_NoObjWithObjCutNLPEval::eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value)
{
    #if MRQ_DEBUG_MODE
        assert(false);
    #else
        if(newx)
            mynewx = true;
        
        value = 0.0;
    #endif
    
    return 0;
}


int MRQ_NoObjWithObjCutNLPEval::eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *constrEval, const double *x, double *values)
{
    const int objCutIndex = origm;
    int r1 = 0, r2 = 0;
    
    if(newx)
        mynewx = true;
    
    if( nlConstrs )
    {
        r1 = oeval->eval_nl_constrs_part(threadnumber, n, origm, mynewx, constrEval, x, values);
        
        #if MRQ_DEBUG_MODE
            if(r1 != 0)
                MRQ_PRINTCALLBACKERRORNUMBER(r1);
        #endif
        
        mynewx = false;
    }
    
    
    if( constrEval[objCutIndex] ) //objective cut
    {
        double obj;
        r2 = oeval->eval_nl_obj_part(threadnumber, n, mynewx, x, obj );
        
        #if MRQ_DEBUG_MODE
            if(r2 != 0)
                MRQ_PRINTCALLBACKERRORNUMBER(r2);
        #endif
        
        mynewx = false;
        
        if( objFactor != 1.0 )
            obj *= objFactor;
        
        values[objCutIndex] = obj;
    }
    
    
    return r1 == 0 ? r2 : r1;
}


int MRQ_NoObjWithObjCutNLPEval::eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *values)
{
    
    #if MRQ_DEBUG_MODE
        assert(false);
    #else
        if(newx)
            mynewx = true;
        
        for(int i = 0; i < n; i++)
            values[i] = 0.0;
    #endif
    
    return 0;
}


int MRQ_NoObjWithObjCutNLPEval::eval_grad_nl_constrs_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *constrEval, const double *x, minlpproblem::MIP_SparseMatrix& jacobian)
{
    const int objCutIndex = origm;
    int r1 = 0, r2 = 0;
    
    if(newx)
        mynewx = true;
    
    if( nlConstrs )
    {
        r1 = oeval->eval_grad_nl_constrs_part( threadnumber, n, origm, nz, mynewx, constrEval, x, jacobian );
        
        #if MRQ_DEBUG_MODE
            if(r1 != 0)
                MRQ_PRINTCALLBACKERRORNUMBER(r1);
        #endif
        
        mynewx = false;
    }
    
    
    if( constrEval[objCutIndex] )
    {
        double *grad = jacobian(objCutIndex);
        
        //r2 = oeval->eval_grad_nl_obj_part( threadnumber, n, mynewx, x, auxVals );
        r2 = oeval->eval_grad_nl_obj_part( threadnumber, n, mynewx, x, grad );
        
        #if MRQ_DEBUG_MODE
            if(r2 != 0)
                MRQ_PRINTCALLBACKERRORNUMBER(r2);
        #endif
        
        mynewx = false;
        
        if(objFactor != 1.0)
        {
            MRQ_multiplyAllArray(n, grad, objFactor);
            //for(int i = 0; i < n; i++)
                //grad[i] *= objFactor;
        }
        //jacobian[objCutIndex].setElementsByOrder(n, auxVals);
    }
    
    return r1 == 0 ? r2 : r1;
}


int MRQ_NoObjWithObjCutNLPEval::eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double objFactor, const double *lambda, minlpproblem::MIP_SparseMatrix& hessian)
{
    int r;
    const int objCutIndex = origm;
    
    const double objf = this->objFactor * lambda[objCutIndex];
    
    
    if(newx)
        mynewx = true;
    
    r = oeval->eval_hessian_nl_lagran_part( threadnumber, n, origm, nz, mynewx, x, objf, lambda, hessian );
    
    #if MRQ_DEBUG_MODE
        if(r != 0)
            MRQ_PRINTCALLBACKERRORNUMBER(r);
    #endif
    
    mynewx = false;
    
    return r;
}







MRQ_GapMinProb::MRQ_GapMinProb()
{
    initialize();
}



MRQ_GapMinProb::~MRQ_GapMinProb()
{
    desallocate();
}


int MRQ_GapMinProb::allocateAuxStructures(const int n)
{
    MRQ_malloc(auxCols, n+1); //auxCols = (int *) malloc( (n+1) * sizeof(int) );
    MRQ_malloc(auxVals, n+1); //auxVals = (double *) malloc( (n+1) * sizeof(double) );
    MRQ_IFMEMERRORRETURN( !auxCols || !auxVals );
    
    return 0;
}


void MRQ_GapMinProb::desallocate()
{
    MRQ_secFree(auxCols);
    MRQ_secFree(auxVals);
    MRQ_secDelete(solver);
    MRQ_secDelete(eval);
}



void MRQ_GapMinProb::initialize()
{
    //objCut = false;
    auxObjVarIndex = -1;
    intGapConstrIndex = -1;
    objCutIndex = -1;
    auxCols = NULL;
    auxVals = NULL;
    solver = NULL;
    eval = NULL;
}


//if setAuxVarObj, we add an auxiliary variable z to objezctive function and perform
//min z
//s.t.: sum_{i} -x_{i}² + x_i <= z
int MRQ_GapMinProb::setProblem(const int solverCode, MRQ_MINLPProb &prob, const double *lx, const double *ux, MRQ_GeneralSolverParams *params, const int thnumber, const bool setSpecialParams, const bool setGapExpOnConstr, const bool setGapUpperBound, const int nthreads, const double maxCpuTime, const double maxTime, const int naddvars, const int naddconstrs)
{
    const bool hasNlConstrs = prob.hasNLConstraints();
    
    const int n = prob.n;
    const int m = prob.m;
    const int *xtype = prob.xtype;
    
    bool setGapOnConstr = setGapExpOnConstr || setGapUpperBound;
    
    
    int nI = 0, r;
    
    
    r = allocateAuxStructures(n);
    solver = OPT_newQPSolver( solverCode );
    MRQ_IFMEMERRORRETURN( r != 0 || !solver );
    
    
    r = MRQ_setNLPRelaxProb(prob, lx, ux, NULL, NULL, solver, false, true, true, false, thnumber, setSpecialParams, params, nthreads, maxCpuTime, maxTime, naddvars, naddconstrs);
    MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    
    r = solver->setnVariablesBounds( n, lx, ux );
    MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    
    if( hasNlConstrs )
    {
        //eval = new (nothrow) MRQ_NLPMinGapEval( prob.getNumberOfConstraints(), hasNlConstrs, prob.getNonLinearEvaluationObject() );
        
        eval = new (std::nothrow) MRQ_NLPNonObjEval(n, m, prob.getNumberOfJacobianNonZeros(), prob.getNumberOfLagrangianHessianNonZeros(), prob.getNonLinearEvaluationObject() );
        MRQ_IFMEMERRORRETURN( !eval );
        
        //eval->auxVals = auxVals;
        
        ( (OPT_NLPSolver *) solver )->setNonLinearEvalObject( eval );
    }
    
    
    //setting objective function
    
    
    for(int i = 0; i < n; i++)
    {
        //const int ind = xtype[i];
        
        //we disconsider binary variables already fixed...
        if( minlpproblem::MIP_isIntegerType( xtype[i] )  &&  (lx[i] > -2.0 && lx[i] < 1.0)  &&  (ux[i] > 0.0 && ux[i] < 2.0)  )
        {
            auxCols[nI] = i;
            nI++;
        }
    }
    
    
    
    solver->setObjSense( optsolvers::OPT_MINIMIZE ); // obj function is undetermined in optsolvers
    
    MRQ_setAllArray(nI, auxVals, 1.0);
    
    
    if( setGapOnConstr )
    {
        //seting problem in this way
        
        //min z
        //s.t.: sum_{i} -x_{i}² + x_i <= z
        // 0<= z <= zb
        
        const double zb = 0.75*(0.25*nI); //0.25 is the maximum gap in our gap function. We allow a measure of gap at most 10 percent of variables having the maximum gap. We hope prune some nodes by infeasibility 
        
        r = solver->addVariables(1);
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        r = solver->addConstraints(1);
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        auxObjVarIndex = prob.n;
        r = solver->setObjLinearCoef( auxObjVarIndex, 1.0 );
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        intGapConstrIndex = prob.m;
        
        r = solver->setConstraintBounds( intGapConstrIndex, -OPT_INFINITY, 0.0 );
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        r = solver->setVariableBounds( auxObjVarIndex, 0.0, setGapUpperBound ? zb : OPT_INFINITY );
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        
        auxCols[nI] = auxObjVarIndex;
        auxVals[nI] = -1.0;
        
        r = solver->setConstraintLinearCoefs( intGapConstrIndex, nI+1, auxCols, auxVals );
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    }
    else
    {
        r = solver->setObjLinearCoefs(nI, auxCols, auxVals);
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    }
    
    
    MRQ_setAllArray(nI, auxVals, -2.0);
    
    if( setGapOnConstr )
    {
        OPT_QCPSolver *qcp = (OPT_QCPSolver *) solver;
        
        r = qcp->setConstraintQuadMatrix( intGapConstrIndex, nI, auxCols, auxCols, auxVals );
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    }
    else
    {
        r = solver->setObjQuadMatrix(nI, auxCols, auxCols, auxVals);
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        #if OPT_HAVE_KNITRO
        if( solver->getSolverCode() == optsolvers::OPT_KNITRO )
        {
            ((OPT_Knitro *)solver)->objFnType = KTR_FNTYPE_NONCONVEX; //there is no difference set it for continuous problems, but we set it anyway...
        }
        #endif
    }
    
    
    //solver->generateModelFile( "igma2.mod" );
    //MRQ_getchar();
    
    
    return 0;
}




int MRQ_GapMinProb::setObjCutConstr(MRQ_MINLPProb& prob, const double zu)
{
    int r, objCutIndex;
    optsolvers::OPT_ObjCutSetter ocset;
    
    
    if( prob.hasObjNLTerm() )
    {
        if(eval)	delete eval;
        
        eval = new (std::nothrow) MRQ_NoObjWithObjCutNLPEval( prob.hasNLConstraints(),  prob.getNonLinearEvaluationObject(), prob.m, prob.objFactor );
        MRQ_IFMEMERRORRETURN(!eval);
        if( !eval )
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTMEMERROR;
            #endif
            return MRQ_MEMORY_ERROR;
        }
        
        //( (MRQ_NLPMinGapEval *) eval)->auxVals = auxVals;
    }
    
    r = solver->getNumberOfConstraints( objCutIndex );
    #if MRQ_DEBUG_MODE
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        /*if( objCutIndex != prob.m )
        {// by now, that assertion should be true. But if some day you add cuts or aditional constraints, that test can fail (remove this assert)
            
            MRQ_PRINTERRORMSG("warning: objective cut added in a index different index than prob.m.");
        } */
    #endif
    
    r = solver->addConstraints(1);
    MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    
    r = ocset.setObjCut( solver, objCutIndex, prob, zu, eval);
    MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
    
    
    this->objCutIndex= objCutIndex;
    
    //objCut = true;
    
    return 0;
}



int MRQ_GapMinProb::updateObjCutConstr(MRQ_MINLPProb& prob, const double zu)
{
    optsolvers::OPT_ObjCutSetter ocset;
    
    if( objCutIndex >= 0 )
    {
        int r = ocset.updateObjCut( solver, objCutIndex, prob, zu );
        
        #if MRQ_DEBUG_MODE
            MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        #endif
        
        return 0;
    }
    else
    {
        return setObjCutConstr(prob, zu);
    }
    
}



MRQ_Rounding::~MRQ_Rounding()
{
}



bool MRQ_Rounding::roundSolution( const int thnumber, MRQ_MINLPProb &prob, const int nI, const int *intVars, const double absFeasTol, const double relFeasTol, const double* x, const double zu, const bool localSearchIfFeas, const bool localSearchIfInfeas, const double* nlx, const double* nux, double *auxConstr, double *auxVars,  MRQ_NLPSolver& nlp, double* sol, double& fsol )
{
    bool feas = roundSolution(thnumber, prob, nI, intVars, absFeasTol, relFeasTol, x, nlx, nux, auxConstr, sol, fsol);
    
    
    
    //printf("Entrei no arredondamento. feas: %d localSearchIfFeas: %d localSearchIfInfeas: %d\n", (int) feas, (int) localSearchIfFeas, (int) localSearchIfInfeas  );
    
    const bool localSearch = feas ? (localSearchIfFeas && fsol < zu): localSearchIfInfeas;
    
    if(localSearch)
    {
        feas = localSearchOverNLPSOlver(prob, nI, intVars, nlx, nux, auxVars, nlp, sol, fsol);
    }
    
    return feas;
}


//return true if we found a feasible solution. Note sol and fsol are input and output arguments. Local search is made around sol fixing integer variables at integer values.
bool MRQ_Rounding::localSearchOverNLPSOlver(MRQ_MINLPProb &prob, const int nI, const int *intVars, const double* nlx, const double* nux, double *auxVars, MRQ_NLPSolver& nlp, double* sol, double& fsol )
{
    const int n = prob.n;
    //const int m = prob.m;
    
    //saving the current state of nlp
    const bool feasNlpOrig = nlp.feasSol;
    const int rcNlpOrig = nlp.retCode;
    const double objNlpOrig = nlp.objValue;
    const double dualNlpOrig = nlp.dualObjValue;
    
    bool feas = false;
    int r;
    
    
    
    MRQ_copyArray( n, nlp.sol, auxVars );
    //MRQ_copyArray( m, nlp.constr, auxConstr );
    
    //fixing integer variables
    for(int i = 0; i < nI; i++)
    {
        //sol was already rounded by roundSolution()
        const int k = intVars[i];
        
        r = nlp.setVariableBounds(k, sol[k], sol[k]);
        #if MRQ_DEBUG_MODE
            if(r != 0)
                MRQ_PRINTERRORNUMBER(r);
        #endif
    }
    
    nlp.solve(false, true, false, false); //in this way, we do not overwrite dual variables and constraint values... 
    
    
    if( nlp.feasSol )
    {
        feas = true;
        fsol = nlp.objValue;
        MRQ_copyArray(n, nlp.sol, sol);
    }
    
    
    //restoring initial state of nlp
    nlp.feasSol = feasNlpOrig;
    nlp.objValue = objNlpOrig;
    nlp.dualObjValue = dualNlpOrig;
    nlp.retCode = rcNlpOrig;
    MRQ_copyArray( n, auxVars, nlp.sol );
    //MRQ_copyArray( m, auxConstr, nlp.constr );
    
    for(int i = 0; i < nI; i++)
    {
        const int k = intVars[i];
        
        r = nlp.setVariableBounds(k, nlx[k], nux[k]);
        #if MRQ_DEBUG_MODE
            if(r != 0)
                MRQ_PRINTERRORNUMBER(r);
        #endif
    }
    
    
    return feas;
}




bool MRQ_Rounding::roundSolution( const int thnumber, MRQ_MINLPProb& prob, const int nI, const int* intVars, const double absFeasTol, const double relFeasTol, const double* x, const double* nlx, const double* nux, double *auxConstr, double* sol, double& fsol )
{
    const int n = prob.n;
    //const int m = prob.m;
    
    bool answer;
    int r;
    
    
    MRQ_copyArray(n, x, sol);
    for(int i = 0; i < nI; i++)
    {
        int k = intVars[i];
        sol[k] = round( sol[k] );
    }
    
    
    
    //we pass auxEval as NULL to evaluate all constraints.
    prob.isFeasibleToConstraints( thnumber, sol, true, NULL, absFeasTol, relFeasTol, answer, auxConstr );
    
    if( answer )
    {
        r = prob.objEval(thnumber, !prob.hasNlConstrs, sol, fsol);
        
        if(r != 0)
            answer = false;
    }
    
    
    return answer;
}


MRQ_RandomRounding::~MRQ_RandomRounding()
{
}


MRQ_RandomRounding::MRQ_RandomRounding(long int seed) : MRQ_Rounding()
{
    random.setSeed(&seed);
}


bool MRQ_RandomRounding::roundSolution( const int thnumber, muriqui::MRQ_MINLPProb& prob, const int nI, const int* intVars, const double absFeasTol, const double relFeasTol, const double* x, const double* nlx, const double* nux, double* auxConstr, double* sol, double& fsol )
{
    const int n = prob.n;
    //const int m = prob.m;
    
    bool answer;
    int r;
    
    
    
    MRQ_copyArray(n, x, sol);
    for(int i = 0; i < nI; i++)
    {
        const int k = intVars[i];
        
        const double prob = sol[k] - floor(sol[k]) ;
        
        //the gap of variable k is the probability to be rounded up
        if( random.random() < prob )
        {
            sol[k] = ceil( sol[k] );
            //we can have bizarous error here if the solution sol[k] pass nux[k] by some decimal places on nlp relaxation solution
            if( sol[k] > nux[k] )
                sol[k] = nux[k];
        }
        else
        {
            sol[k] = floor( sol[k] );
            //we can have bizarous error here if the solution sol[k] pass nlx[k] by some decimal places on nlp relaxation solution
            if( sol[k] < nlx[k] )
                sol[k] = nlx[k];
        }
        
    }
    
    
    //we pass auxCEval as NULL to evaluate all constarints
    r = prob.isFeasibleToConstraints( thnumber, sol, true, NULL, absFeasTol, relFeasTol, answer, auxConstr );
    if(r != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTCALLBACKERRORNUMBER(r);
        #endif
    }
    
    if( answer )
    {
        r = prob.objEval(thnumber, !prob.hasNlConstrs, sol, fsol);
        
        if(r != 0)
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTCALLBACKERRORNUMBER(r);
            #endif
            answer = false;
        }
    }
    
    
    return answer;
}


MRQ_LPboundsUpdater::MRQ_LPboundsUpdater()
{
    lp = NULL;
}



MRQ_LPboundsUpdater::~MRQ_LPboundsUpdater()
{
    deallocate();
}



int MRQ_LPboundsUpdater::buildProblem(const int solverCode, MRQ_MINLPProb &prob, const double zu)
{
    int r;
    
    lp = OPT_newLPSolver(solverCode);
    MRQ_IFMEMERRORRETURN( !lp );
    
    //set variables and constraints
    //r = lp->setProblemFrom(prob, false, true, false, false, true);
    r = lp->setLinearObjAndConstraintsFrom(prob, false, true, false, false);
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    //seting objective cut
    if( prob.hasObjNLTerm() == false && prob.getNumberOfObjQuadTerms() == 0 )
    {
        int m;
        OPT_ObjCutSetter objCutSetter;
        
        r = lp->getNumberOfConstraints(m);
        #if MRQ_DEBUG_MODE
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        #endif
        
        r = lp->addConstraints(1);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        r = objCutSetter.setObjCut( lp, m, prob, zu, NULL );
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    }
    
    return 0;
}



int MRQ_LPboundsUpdater::__optOnVariable(const int index, double &value)
{
    int r, code;
    
    
    r = lp->setObjLinearCoef(index, 1.0);
    #if MRQ_DEBUG_MODE
        MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    #endif
    
    r = lp->solve(false, false, false, false);
    
    if( r == OPT_OPTIMAL_SOLUTION )
    {
        value = lp->objValue;
        code = MRQ_OPTIMAL_SOLUTION;
    }
    else if( r == OPT_INFEASIBLE_PROBLEM )
    {
        code = MRQ_INFEASIBLE_PROBLEM;
    }
    else
    {
        code = MRQ_UNDEFINED_ERROR;
    }
    
    std::cout << "index: " << index << " lp.retCode: " << lp->retCode << " obj: " << lp->objValue << "\n";
    
    //lp->generateModelFile("lpmodelnode.lp");
    //MRQ_getchar();
    
    
termination:
    
    lp->setObjLinearCoef(index, 0.0);
    
    return code;
}


int MRQ_LPboundsUpdater::calcMinValueOnVariable( const int index, double &value)
{
    lp->setObjSense( optsolvers::OPT_MINIMIZE );
    
    return __optOnVariable(index, value);
}


int MRQ_LPboundsUpdater::calcMaxValueOnVariable( const int index, double &value)
{
    lp->setObjSense( optsolvers::OPT_MAXIMIZE );
    
    return __optOnVariable(index, value);
}


void MRQ_LPboundsUpdater::deallocate()
{
    MRQ_secDelete(lp);
}



int MRQ_LPboundsUpdater::setVariablesBounds( const int n, const double *lx, const double *ux )
{
    const int r = lp->setnVariablesBounds( n, lx, ux );
    #if MRQ_DEBUG_MODE
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    #endif
    
    return 0;
}


int MRQ_LPboundsUpdater::updateObjectiveCut( MRQ_MINLPProb &prob, const double zu)
{
    if( prob.hasObjNLTerm() == false && prob.getNumberOfObjQuadTerms() == 0 )
    {
        int m;
        OPT_ObjCutSetter ocs;
        
        lp->getNumberOfConstraints(m);
        
        //std::cout << "zu: " << zu << "\n";
        //MRQ_getchar();
        
        const int r = ocs.updateObjCut(lp, m-1, prob, zu);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    }
    
    return 0;
}





MRQ_BoundsUpdaterSolver::MRQ_BoundsUpdaterSolver()
{
    objCutConstrIndex = UINT_MAX;
    lp = nullptr;
    objCutEval = nullptr;
}



MRQ_BoundsUpdaterSolver::~MRQ_BoundsUpdaterSolver()
{
    deallocate();
}



int MRQ_BoundsUpdaterSolver::buildProblem(const int solverCode, const MRQ_MINLPProb &prob, const bool onlyLinearPart, const bool setObjCut, const double zu, const unsigned thnumber, const unsigned int nThreads, MRQ_GeneralSolverParams *solverParams, const double *lc, const double *uc, const bool setSpecialNlpSolverParams  )
{
    int r;
    
    
    lp = OPT_newLPSolver(solverCode);
    MRQ_IFMEMERRORRETURN( !lp );
    
    
    if( onlyLinearPart )
    {
        r = lp->setLinearObjAndConstraintsFrom(prob, false, true, false, false, 0, 1);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        r = lp->setNumberOfThreads(nThreads);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
        
        if(solverParams)
            lp->setParameters(*solverParams);
        
        
        //seting objective cut
        if( setObjCut )
        {
            if( prob.hasObjNLTerm() == false && prob.getNumberOfObjQuadTerms() == 0 )  //by now, objective cut only is implemented with linear objective function
            {
                int m;
                OPT_ObjCutSetter objCutSetter;
                
                r = lp->getNumberOfConstraints(m);
                MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
                
                objCutConstrIndex = m-1;
                
                r = objCutSetter.setObjCut( lp, objCutConstrIndex, prob, zu, nullptr );
                MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            }
            
        }
        
    }
    else
    {
        const auto origm = prob.m;
        
        r = MRQ_setNLPRelaxProb(prob, NULL, NULL, lc, uc, lp, false, true, true, false, thnumber, setSpecialNlpSolverParams, solverParams, nThreads, INFINITY, INFINITY, 0, 1);
        MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
        
        
        if( setObjCut )
        { //OBJECTIVE CUT WITH NONLINEAR PROBLEMAS IS NOT TESTED ENOUGH...
            bool hasNlConstrs = prob.hasNLConstraints();
            OPT_ObjCutSetter objCutSetter;
            
            objCutConstrIndex = origm;
            
            if( objCutEval )
                delete objCutEval;
            
            if( prob.hasNlObj || hasNlConstrs )
            {
                objCutEval = new (std::nothrow) MRQ_NoObjWithObjCutNLPEval( hasNlConstrs, prob.nlEval, origm, prob.objFactor );
                MRQ_IFMEMERRORRETURN( !objCutEval );
            }
            
            /*r = lp->getNumberOfConstraints(m);
            #if MRQ_DEBUG_MODE
                MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            #endif
            
            r = lp->addConstraints(1);
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR); */
            
            r = objCutSetter.setObjCut( lp, objCutConstrIndex, prob, zu, objCutEval );
            MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
            
            
            //MRQ_PRINTERRORMSG("We are sorry, but, by now, objective cut is only implemented to linear objetive!");
            //assert(false);
        }
        
    }
    
    
    return 0;
}



int MRQ_BoundsUpdaterSolver::__optOnVariable(const int index, double &value)
{
    int r, code;
    
    
    r = lp->setObjLinearCoef(index, 1.0);
    #if MRQ_DEBUG_MODE
        MRQ_IFERRORGOTOLABEL(r, code, MRQ_MILP_SOLVER_ERROR, termination);
    #endif
    
    r = lp->solve(false, false, false, false);
    
    if( r == OPT_OPTIMAL_SOLUTION )
    {
        value = lp->objValue;
        code = MRQ_OPTIMAL_SOLUTION;
    }
    else if( r == OPT_INFEASIBLE_PROBLEM )
    {
        code = MRQ_INFEASIBLE_PROBLEM;
    }
    else
    {
        code = MRQ_UNDEFINED_ERROR;
    }
    
    //std::cout << "index: " << index << " lp.retCode: " << lp->retCode << " obj: " << lp->objValue << "\n";
    
    //lp->generateModelFile("lpmodelnode.lp");
    //MRQ_getchar();
    
    
termination:
    
    lp->setObjLinearCoef(index, 0.0);
    
    return code;
}


int MRQ_BoundsUpdaterSolver::calcMinValueOnVariable( const int index, double &value)
{
    lp->setObjSense( optsolvers::OPT_MINIMIZE );
    
    return __optOnVariable(index, value);
}


int MRQ_BoundsUpdaterSolver::calcMaxValueOnVariable( const int index, double &value)
{
    lp->setObjSense( optsolvers::OPT_MAXIMIZE );
    
    return __optOnVariable(index, value);
}


void MRQ_BoundsUpdaterSolver::deallocate()
{
    MRQ_secDelete(lp);
    MRQ_secDelete(objCutEval);
}



int MRQ_BoundsUpdaterSolver::setVariablesBounds( const int n, const double *lx, const double *ux )
{
    const int r = lp->setnVariablesBounds( n, lx, ux );
    MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    
    return 0;
}


int MRQ_BoundsUpdaterSolver::updateObjectiveCut( const MRQ_MINLPProb &prob, const double zu)
{
    if( objCutConstrIndex != UINT_MAX )
    {
        OPT_ObjCutSetter ocs;
        
        //std::cout << "zu: " << zu << "\n";
        //MRQ_getchar();
        
        const int r = ocs.updateObjCut(lp, objCutConstrIndex, prob, zu);
        MRQ_IFERRORRETURN(r, MRQ_MILP_SOLVER_ERROR);
    }
    
    return 0;
}


int MRQ_BoundsUpdaterSolver::calculateNewVarBounds( const unsigned int n, const int nIndVars, const int *indVars, double *lx, double *ux, const bool allIntVars, bool &infeasible )
{
    int r, rlb, rub;
    double vlb, vub;
    const double tol = 1e-4;
    const double tolToUpdating = 1e-2;
    
    infeasible = false;
    
    r = setVariablesBounds(n, lx, ux);
    MRQ_IFERRORRETURN(r, r);
    
    for( int i = 0; i < nIndVars; i++ )
    {
        const int ind = indVars[i];
        
        if( lx[ind] == ux[ind] )
            continue;
        
        rlb = calcMinValueOnVariable(ind, vlb);
        
        //printf("ind: %d lx: %lf rlb: %d vlb: %lf", ind, lx[ind], rlb, vlb);
        
        if( rlb == MRQ_OPTIMAL_SOLUTION )
        {
            #if MRQ_DEBUG_MODE
                assert( vlb >= lx[ind] - tol );
            #endif
            
            if( vlb > lx[ind] + tolToUpdating )
            {
                if( allIntVars )
                    lx[ind] = ceil( vlb );
                else
                    lx[ind] = vlb;
                
                //printf(" ATUALIZEI");
                //MRQ_getchar();
                
                r = lp->setVariableBounds(ind, lx[ind], ux[ind]);
                MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
            }
        }
        else if( rlb == MRQ_INFEASIBLE_PROBLEM )
        {
            //printf("\ndetectei inviabilidade!\n");
            //MRQ_getchar();
            infeasible = true;
            break;
        }
        
        
        rub = calcMaxValueOnVariable(ind, vub);
        
        //printf("\t\t\t    ux: %lf rub: %d vub: %lf", ux[ind], rub, vub);
        
        if( rub == MRQ_OPTIMAL_SOLUTION )
        {
            #if MRQ_DEBUG_MODE
                assert( vub <= ux[ind] + tol );
                
                if( rlb == MRQ_OPTIMAL_SOLUTION )
                    assert( vlb <= vub + tol );
            #endif
            
            if( vub < ux[ind] - tolToUpdating )
            {
                if( allIntVars )
                    ux[ind] = floor( vub );
                else
                    ux[ind] = vub;
                
                r = lp->setVariableBounds(ind, lx[ind], ux[ind]);
                MRQ_IFERRORRETURN(r, MRQ_NLP_SOLVER_ERROR);
                
                //printf(" ATUALIZEI");
                //MRQ_getchar();
            }
            
            
            
        }
        else if( rub == MRQ_INFEASIBLE_PROBLEM )
        {
            /*printf("\ndetectei inviabilidade!\n");
            MRQ_getchar();*/
            infeasible = true;
            break;
        }
        
        
        if( MRQ_abs(ux[ind] - lx[ind]) <= tol )
        { //we assume variable is fixed
            double v = 0.5*(ux[ind] + lx[ind]);
            
            ux[ind] = v;
            lx[ind] = v;
        }
        else if( ux[ind] < lx[ind] - tol )
        {
            /*printf("\ndetectei inviabilidade!\n");
            MRQ_getchar();*/
            infeasible = true;
            break;
        }
        
        //printf("\n");
    }
    
    //MRQ_getchar();
    
    return 0;
}




/*lx and ux are input/output arguments
 * zu should be an objective value of a feasible solution to the minlp problem
 * 
 */
int MRQ_BinVarsOptNlpRelaxSolFixer:: fixBinVarsFromNlpRelaxSol( const int nI, const int *intVars, const double zu, const double objRelaxSol, const double *duallx, const double *dualux, int &nFixed, double *lx, double *ux, MRQ_NewBBNode *node)
{
    const double ZERO = 1e-4;
    
    branchAndBound::BBL_ClassUnionNodeBoundsSolPointer *myBounds = node ? &(node->myBounds) : NULL;
    
    nFixed = 0;
    
    
    for( int k = 0; k < nI; k++ )
    {
        const int ind = intVars[k];
        bool fixed = false;
        
        if( lx[ind] == ux[ind] )
            continue; //variable already fixed
        
        if( lx[ind] > -1.0 && ux[ind] < 2.0 ) //binary var
        {
            
            //std::cout << "ind: " << ind << " - lx: " <<lx[ind] << " ux: " << ux[ind] << " zu: " << zu << " obj relax: " << objRelaxSol << " duallx: " << duallx[ind] << " dualux: " << dualux[ind] << "\n";
            
            
            //testing if we should fix to zero.
            if( duallx[ind] > ZERO  && (zu - objRelaxSol)/duallx[ind] < 1.0 )
            {
                #if MRQ_DEBUG_MODE
                    assert(lx[ind] < 1.0);
                    //we test now 
                #endif
                
                ux[ind] = 0.0;
                fixed = true;
            }
            
            //do no replace this else by a else ifsince we would like to test if there is some case attending both if's
            if( dualux[ind] > ZERO && (1.0 - (zu - objRelaxSol)/dualux[ind]) > 0.0 )
            {
                //if(fixed)
                    //std::cout << " Double fixed!\n";
                
                /*#if MRQ_DEBUG_MODE
                    assert(ux[ind] > 0.0);
                #endif*/
                
                if( ux[ind] == 0.0 )
                    MRQ_PRINTERRORMSG("Warning: two bounds updated for the same variable in dual valuefixing! This should be checked!");
                
                lx[ind] = 1.0;
                fixed = true;
            }
            
            //I think it is not possible fixing at 0 and 1 at same time because only one of these constraints can be active (maibe neither of them). So, only the active bound could have nonzero dual value
            
            if(fixed)
            {
                nFixed++;
                
                if( myBounds )
                {
                    
                    int r = myBounds->allocateElements( node->nMyBounds + 1 );
                    MRQ_IFERRORRETURN(r, MRQ_MEMORY_ERROR);
                    
                    myBounds->addElementOnSortedArray(node->nMyBounds, ind, lx[ind], ux[ind], NAN );
                }
            }
        }
        
    }
    
    
    return 0;
}



