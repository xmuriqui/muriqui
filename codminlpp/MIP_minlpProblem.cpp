/*
* minlpProb.cpp
*
*  Created on: 27/08/2013
*      Author: Wendel Melo
*/


#include <cmath>

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <climits>

#include <new>
#include <iostream>
#include <iomanip>
#include <fstream> //For writting in files



#include "MIP_minlpProblem.hpp"
#include "MIP_tools.hpp"


#define MIP_N_POINTS_TO_APP_DERIVATIVE 5




using namespace std;
using namespace newspm;
using namespace minlpproblem;




//we need 4 points to do the approx, two before, and two after the point. So, we assume f has 5 indeces. The first 2 points are the points after, the third point is the point that we want approximate and the last 2 points are the two points after
inline double MIP_AppDerivativeByCenteredDif5(const double *f, const double h)
{
    //that is wrong! We need the fifth dervative to calculate that...
    return (f[0] + 8*(f[3] - f[1]) -f[4])/( 12*h ); //we put the parenthesis inside numerator to avoid numerical problems
}



inline double MIP_AppDerivativeByCenteredDif(const int npoints, const double *f, const double h)
{
    if( npoints == 5 )
        return MIP_AppDerivativeByCenteredDif5(f, h);
    else
        return (f[2] - f[0])/(2*h);
}


//to vectorial functions. (We use to approximate hessians having the gradients)
inline void MIP_AppDerivativeByCenteredDif(const int n, const int npoints, double **f, const double h, double *grad)
{
    int i;
    
    if( npoints == 5 )
    {
        //(f[0] + 8*(-f[1] + f[3]) -f[4])/( 12*h );
        
        for(i = 0; i < n; i++)
            grad[i] = (f[0][i] + 8*(f[3][i] - f[1][i]) - f[4][i])/(12*h) ;
    }
    else
    {
        for(i = 0; i < n; i++)
            grad[i] = (f[2][i] - f[0][i])/(2*h);
    }
    
}




MIP_MINLPProb::MIP_MINLPProb()
{
    initialize();
}


void MIP_MINLPProb::initialize()
{
    objLinearPart = false;
    hasNlObj = false;
    hasNlConstrs = false;
    nlConstr = NULL;
    n = nI = m = 0;
    xtype = NULL;
    xprior = NULL;
    
    lx = ux = NULL;
    c = NULL;
    d = 0.0;
    
    lc = uc = NULL;
    x = NULL;
    //xInit = NULL;
    objFactor = 1.0;
    objValue = NAN;
    
    QC = NULL;
    nlEval = NULL;
    
    Q.setSymmetricFlag(true);
    lagH.setSymmetricFlag(true);
}


MIP_MINLPProb::~MIP_MINLPProb()
{
    deallocateMatrices();
}


int MIP_MINLPProb::addConstraints(const int ncons)
{
    int i, ret, code;
    //bool *auxb;
    //double *auxd;
    //MIP_SparseMatrix *auxM;
    
    if( ncons <= 0 )
        return 0;
    
    
    ret = MIP_realloc(nlConstr, (ncons + m)); //auxb = (bool *) realloc(nlConstr, (ncons + m)*sizeof(bool) );
    if( ret != 0 )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    //nlConstr = auxb;
    
    for(i = m; i < m + ncons; i++)
        nlConstr[i] = false;
    
    
    ret = MIP_realloc(lc, (ncons + m)); //auxd = (double *) realloc(lc, (ncons + m)*sizeof(double) );
    if( ret != 0 )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    //lc = auxd;
    
    MIP_setAllArray(ncons, &lc[m], -MIP_INFINITY); //for(i = m; i < m + ncons; i++)	lc[i] = -MIP_INFINITY;
    
    
    ret = MIP_realloc(uc, (ncons + m)); //auxd = (double *) realloc(uc, (ncons + m)*sizeof(double) );
    if( ret != 0 )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    //uc = auxd;
    
    MIP_setAllArray(ncons, &uc[m], MIP_INFINITY); //for(i = m; i < m + ncons; i++) uc[i] =  MIP_INFINITY;
    
    
    ret = MIP_realloc(QC, (ncons + m)); //auxM = (MIP_SparseMatrix *) realloc(QC, (ncons + m)*sizeof(MIP_SparseMatrix));
    if( ret != 0 )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    //QC = auxM;
    
    
    for(i = m; i < m + ncons; i++)
        QC[i].initialize(0, n, true);
    
    
    
    ret = A.addNewRows(ncons);
    ret += J.addNewRows(ncons);
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    
    
    m += ncons;
    
    code = 0;
    
termination:
    
    
    if( code != 0 )
    {
        //that is not so efficient, but we only execute this if we had some error, and probably, we will have to abort the main procedure anyway...
            
        unsigned int end = m+ncons;
        
        if( A.getNumberOfRows() == end )
        {
            for(unsigned int i = m; i < end; i++)
                A.removeRows(1, &i); 
        }
        
        if( J.getNumberOfRows() == end )
        {
            for(unsigned int i = m; i < end; i++)
                J.removeRows(1, &i);
        }
        
    }
    
    
    #if MIP_DEBUG_MODE
        assert( A.getNumberOfRows() == (unsigned int) m ); //making sure the number will not be incorrect if we have some memory error
        assert( J.getNumberOfRows() == (unsigned int) m );
    #endif
    
    
    return code;
}



int MIP_MINLPProb::addVariables(const int nvars)
{
    const int ns = n + nvars;// new size
    int r, code;
    
    //bool *auxb;
    int *auxi;
    double *auxd;
    
    
    auxi = (int *) realloc( xprior, ns * sizeof(int) );
    if(!auxi)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    xprior = auxi;
    MIP_setAllArray(nvars, &xprior[n], 0);
    
    
    auxi = (int *) realloc( xtype, ns * sizeof(int) );
    if(!auxi)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    xtype = auxi;
    MIP_setAllArray(nvars, &xtype[n], (int) MIP_VT_CONTINUOUS);
    //for(int i = n; i < ns; i++)
        //xtype[i] = MIP_VT_CONTINUOUS;
    
    
    auxd = (double *) realloc( lx, ns * sizeof(double) );
    if( !auxd )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    lx = auxd;
    MIP_setAllArray(nvars, &lx[n], -MIP_INFINITY);
    //for(int i = n; i < ns; i++)
        //lx[i] = -MIP_INFINITY;
    
    
    auxd = (double *) realloc( ux, ns * sizeof(double) );
    if( !auxd )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    ux = auxd;
    MIP_setAllArray(nvars, &ux[n], MIP_INFINITY);
    //for(int i = n; i < ns; i++)
        //ux[i] = MIP_INFINITY;
    
    
    auxd = (double *) realloc( c, ns* sizeof(double) );
    if( !auxd )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    c = auxd;
    MIP_setAllArray(nvars, &c[n], 0.0);
    //for(int i = n; i < ns; i++)
        //c[i] = 0.0;
    
    
    auxd = (double *) realloc( x, ns * sizeof(double) );
    if( !auxd )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    x = auxd;
    MIP_setAllArray<double>(nvars, &x[n], NAN);
    //for(int i = n; i < ns; i++ )
        //x[i] = NAN;
    
    
    A.setNumberOfColumns(ns);
    J.setNumberOfColumns(ns);
    
    
    Q.setNumberOfColumns(ns);
    lagH.setNumberOfColumns(ns);
    
    r = Q.addNewRows( nvars );
    r += lagH.addNewRows( nvars );
    
    if( r != 0 )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    
    for(int i = 0; i < m; i++)
    {
        if(QC[i].getNumberOfRows() > 0)
        {
            #if MIP_DEBUG_MODE
                assert( (int) QC[i].getNumberOfRows() == n );
            #endif
            
            r = QC[i].addNewRows( nvars );
            
            if( r != 0 )
            {
                #if MIP_DEBUG_MODE
                    MIP_PRINTMEMERROR;
                #endif
                code = MIP_MEMORY_ERROR;
                goto termination;
            }
        }
        
        QC[i].setNumberOfColumns( ns );
    }
    
    
    
    n = ns;
    code = 0;
    
termination:
    
    if( code != 0 )
    {//we try retrieve sparse matrix states...
        
        cout << "code: " << code << endl;
        
        A.setNumberOfColumns(n);
        J.setNumberOfColumns(n);
        
        
        if( Q.getNumberOfRows() == (unsigned int) ns )
        {
            for(int i = n; i < ns; i++)
                Q.removeRows(1, &i);
            
            Q.setNumberOfColumns(n);
        }
        
        if( lagH.getNumberOfRows() == (unsigned int) ns )
        {
            for(int i = n; i < ns; i++)
                lagH.removeRows(1, &i);
            
            lagH.setNumberOfColumns(n);
        }
        
        
        for(int i = 0; i < m; i++)
        {
            if( QC[i].getNumberOfRows() == (unsigned int) ns )
            {
                for(int j = n; j< ns; j++)
                    QC[i].removeRows(1, &j);
                
                QC[i].setNumberOfColumns(n);
            }
        }
        
    }
    
    return code;
}


//steps[i] has the step value (h) for i-th cordinate. If steps == NULL, all coordinates will use step as step value (h)
int MIP_MINLPProb::checkFisrtDerivatives(bool checkObj, bool checkConstr, const double* x, const double step, const double* steps, const double tolerance, bool& answer) const
{
    const int npoints = MIP_N_POINTS_TO_APP_DERIVATIVE;
    
    bool newx = true;
    int i, j, k, ret, code;
    double h = step, appDer;
    double *fgrad = NULL, *xp = NULL;
    double *f = NULL, **g = NULL, **ggrad = NULL;
    MIP_SparseMatrix Jac;
    
    
    if( !hasNlObj )
        checkObj = false;
    
    if( !hasNlConstrs )
        checkConstr = false;
    
    
    answer = true;
    
    
    if( !checkObj && !checkConstr )
    {
        //we do not have nonlinear terms. So, nothing to check... 
        code = 0;
        goto termination;
    }
    
    
    //we have to call initialize because we will perform function evaluations...
    i = nlEval->initialize(1, n, m, J.getNumberOfElements(), lagH.getNumberOfElements() );
    MIP_IFCALLBACKERRORGOTOLABEL(i, code, termination);
    
    
    MIP_malloc(xp, n);
    MIP_malloc(f, npoints); 
    MIP_IFMEMERRORGOTOLABEL(!xp || !f, code, termination);
    
    
    if( checkObj )
    {
        MIP_malloc(fgrad, n);
        MIP_IFMEMERRORGOTOLABEL(!fgrad, code, termination);
        
        ret = nlEval->eval_grad_nl_obj_part(0, n, true, x, fgrad);
        MIP_IFCALLBACKERRORGOTOLABEL(ret, code, termination);
        
        newx = false;
    }
    
    
    if( checkConstr )
    {
        //g = (double **) malloc( npoints * sizeof(double *) );
        //ggrad = (double **) malloc( m * sizeof(double *) );
        ret = MIP_allocateMatrix(g, npoints, m);
        ret += MIP_allocateZeroMatrix(ggrad, m, n);
        ret += Jac.copyStructureFrom( J );
        
        MIP_IFMEMERRORGOTOLABEL(ret, code, termination);
        
        
        /*g[0] = (double *) malloc( npoints*m*sizeof(double) );
        ggrad[0] = (double *) calloc( m*n, sizeof(double) );
        
        if( !g[0] || !ggrad[0] )
        {
            code = MIP_MEMORY_ERROR;
            goto termination;
        }
        
        
        for(i = 1; i < npoints; i++)
            g[i] = &g[0][i*m];
        
        for(i = 1; i < m; i++)
            ggrad[i] = &ggrad[0][i*n]; */
        
        
        //maybe it is not necessary, but we avoid write directly in J
        ret = nlEval->eval_grad_nl_constrs_part(0, n, m, J.getNumberOfElements(), newx, nlConstr, x, Jac);
        MIP_IFCALLBACKERRORGOTOLABEL(ret, code, termination);
        
        
        Jac.copyMatrixTo(ggrad, false, false, 1.0); //Jac.copyMatrixTo(ggrad, false);
        
        
        /*for(int i = 0; i < m; i++)
        {
            for(int j = 0; j < n; j++)
                printf("%0.2f ", ggrad[i][j]);
            printf("\n");
        } */
        
    }
    
    
    
    MIP_copyArray(n, x, xp);
    
    
    for(i = 0; i < n; i++)
    {
        if(steps)
            h = steps[i];
        
        ret = 0;
        for(j = 1; j <= npoints/2; j++)
        {
            xp[i] = x[i] + j*h;
            ret += nlObjAndConstraintsEval( checkObj, checkConstr, 0, true, NULL, xp, f[npoints/2 + j], g[npoints/2 + j] );
            
            xp[i] = x[i] - j*h;
            ret += nlObjAndConstraintsEval( checkObj, checkConstr, 0, true, NULL, xp, f[npoints/2 - j], g[npoints/2 - j] );
        }
        
        
        
        /*xp[i] += h;
        ret = nlObjAndConstraintsEval( checkObj, checkConstr, 0, true, NULL, xp, f[3], g[3] );
        
        
        xp[i] += h;
        ret += nlObjAndConstraintsEval( checkObj, checkConstr, 0, true, NULL, xp, f[4], g[4] );
        
        
        xp[i] = x[i] -h;
        ret += nlObjAndConstraintsEval( checkObj, checkConstr, 0, true, NULL, xp, f[1], g[1] );
        
        
        xp[i] -= h;
        ret += nlObjAndConstraintsEval( checkObj, checkConstr, 0, true, NULL, xp, f[0], g[0] ); */
        
        /*xp[i] = x[i] + h;
        ret = nlObjAndConstraintsEval( checkObj, checkConstr, 0, true, NULL, xp, f[2], g[2] );
        
        xp[i] = x[i] - h;
        ret = nlObjAndConstraintsEval( checkObj, checkConstr, 0, true, NULL, xp, f[0], g[0] ); */
        
        if(ret != 0)
        {
            code = MIP_CALLBACK_FUNCTION_ERROR;
            goto termination;
        }
        
        /*for(int i = 0; i < npoints; i++)
        {
            if(i == npoints/2  )
                continue;
            
            printf("P%d =>", i);
            for(int j = 0; j < m; j++)
                    printf(" %f", g[i][j]);
            
            printf("\n");
        }
        getchar(); */
        
        
        
        //checking obj function
        if( checkObj )
        {
            appDer = MIP_AppDerivativeByCenteredDif(npoints, f, h);
            
            if( MIP_abs(fgrad[i] - appDer) > tolerance )
            {
                printf("coord: %d - exact grad obj: %f app grad obj: %f dif: %f\n", i, fgrad[i], appDer, MIP_abs(fgrad[i] - appDer));
                answer = false;
            }
        }
        
        
        
        
        
        if( checkConstr )
        {
            for(j = 0; j < m; j++)
            {
                if( nlConstr[j] )
                {
                    for(k = 0; k < npoints; k++)
                        f[k] = g[k][j];
                    
                    appDer = MIP_AppDerivativeByCenteredDif(npoints, f, h);
                    
                    
                    if( MIP_abs(ggrad[j][i] - appDer) > tolerance )
                    {
                        printf("coord: %d - constraint: %d exact grad: %f  app grad: %f dif: %f", i, j, ggrad[j][i], appDer, MIP_abs(ggrad[j][i] - appDer) );
                        
                        if( ggrad[j][i] == 0.0 )
                        {
                            if( !Jac.hasIndex(j, i) )
                                printf(" Index not described in jacobian structure");
                        }
                        
                        printf("\n");
                        
                        answer = false;
                    }
                }
            }
        }
        
        
        xp[i] = x[i];
        
        //getchar();
    }
    
    
    
    
    code = 0;
    
termination:
    
    if( checkObj || checkConstr )
    {
        nlEval->finalize(1, n, m, J.getNumberOfElements(), lagH.getNumberOfElements() );
    }
    
    if( fgrad )		free(fgrad);
    if( xp )		free(xp);
    if(f)			free(f);
    
    if(g)
    {
        if(g[0])	free(g[0]);
        free(g);
    }
    
    if(ggrad)
    {
        if(ggrad[0])	free(ggrad[0]);
        free(ggrad);
    }
    
    
    if(code != 0)
        answer = false;
    
    
    return code;
}



//check obj and Constr Derivatives. If lambda is NULL, we set lambda
int MIP_MINLPProb::checkSecondDerivatives(bool checkObj, bool checkConstr, const double* x, const double objFactor, const double* lambda, const double step, const double* steps, const double tolerance, bool& answer) const
{
    const int npoints = MIP_N_POINTS_TO_APP_DERIVATIVE;
    
    bool newx = true;
    int i, j, k, w, ret, code;
    double h = step, lamb;
    double *p, *xp = NULL, *myLambda = NULL, *Hline = NULL;
    double **fgrad = NULL, **fHess = NULL, ***gHess = NULL;
    
    
    MIP_SparseMatrix H, *Jacs = NULL;
    
    if( !hasNlObj || objFactor == 0.0 )
        checkObj = false;
    
    if( !hasNlConstrs )
        checkConstr = false;
    
    answer = true;
    
    if( !checkObj && !checkConstr )
    {
        code = 0;
        goto termination;
    }
    
    //we have to call initialize because we will perform function evaluations...
    i = nlEval->initialize(1, n, m, J.getNumberOfElements(), lagH.getNumberOfElements());
    if(i != 0)
    {
        code = MIP_CALLBACK_FUNCTION_ERROR;
        goto termination;
    }
    
    
    MIP_malloc(xp, n); //xp = (double *) malloc(n * sizeof(double));
    MIP_calloc(myLambda, m); //myLambda = (double *) calloc( m , sizeof(double) );
    MIP_malloc(Hline, n); //Hline = (double *) malloc( n * sizeof(double) );
    
    ret = MIP_allocateMatrix(fgrad, npoints, n);
    ret += MIP_allocateZeroMatrix(fHess, n, n);
    ret += H.copyStructureFrom( lagH );
    
    if( !xp || !myLambda || !Hline || ret != 0 )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    
    if( checkObj )
    {
        ret = nlEval->eval_hessian_nl_lagran_part(0, n, m, lagH.getNumberOfElements(), true, x, objFactor, myLambda, H);
        
        if( ret )
        {
            code = MIP_CALLBACK_FUNCTION_ERROR;
            goto termination;
        }
        
        H.copyMatrixTo(fHess, true, false, 1.0); //H.copyMatrixTo(fHess, false);
        
        newx = false;
    }
    
    
    if( checkConstr )
    {
        Jacs = new (nothrow) MIP_SparseMatrix[npoints];
        MIP_malloc(gHess, m); //gHess = (double ***) malloc( m * sizeof(double **) ); 
        if(!Jacs || !gHess)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTMEMERROR;
            #endif
            code = MIP_MEMORY_ERROR;
            goto termination;
        }
        
        
        for(i = 0; i < m; i++)
        {
            ret = MIP_allocateZeroMatrix(gHess[i], n, n);
            if(ret)
            {
                code = MIP_MEMORY_ERROR;
                goto termination;
            }
        }
        
        
        for(i = 0; i < m; i++)
        {
            if(lambda)
                myLambda[i] = lambda[i];
            else
                myLambda[i] = 1.0;
            
            ret = nlEval->eval_hessian_nl_lagran_part(0, n, m, lagH.getNumberOfElements(), newx, x, 0.0, myLambda, H);
            
            if(ret)
            {
                code = MIP_CALLBACK_FUNCTION_ERROR;
                goto termination;
            }
            
            H.copyMatrixTo( gHess[i], true, true, 1.0 );
            
            myLambda[i] = 0.0;
        }
        
        
        
        ret = 0;
        for(i = 0; i < npoints; i++)
            ret += Jacs[i].copyStructureFrom(J);
        
        if(ret)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTERRORNUMBER(ret);
            #endif
            code = MIP_MEMORY_ERROR;
            goto termination;
        }
        
    }
    
    
    MIP_copyArray(n, x, xp);
    
    for(i = 0; i < n; i++)
    {
        if(steps)
            h = steps[i];
        
        ret = 0;
        for(j = 1; j <= npoints/2; j++)
        {
            xp[i] = x[i] + j*h;
            
            if( checkObj )
            {
                ret += nlEval->eval_grad_nl_obj_part(0, n, true, xp, fgrad[npoints/2 + j] );
                
                newx = false;
            }
            else
                newx = true;
            
            if( checkConstr )
                ret += nlEval->eval_grad_nl_constrs_part(0, n, m, J.getNumberOfElements(), newx, nlConstr, xp, Jacs[npoints/2 + j] );
            
            
            xp[i] = x[i] - j*h;
            
            if( checkObj )
            {
                ret += nlEval->eval_grad_nl_obj_part(0, n, true, xp, fgrad[npoints/2 - j] );
                
                newx = false;
            }
            else
                newx = true;
            
            if( checkConstr )
                ret += nlEval->eval_grad_nl_constrs_part(0, n, m, J.getNumberOfElements(), newx, nlConstr, xp, Jacs[npoints/2 - j] );
        }
        
        
        if(ret != 0)
        {
            code = MIP_CALLBACK_FUNCTION_ERROR;
            goto termination;
        }
        
        
        if( checkObj )
        {
            MIP_AppDerivativeByCenteredDif(n, npoints, fgrad, h, Hline);
            
            //we only check lower triangle
            for(j = 0; j <= i; j++)
            {
                if( MIP_abs( objFactor*Hline[j] - fHess[i][j] ) > tolerance )
                {
                    printf("obj hessian (%d, %d) exact: %f app: %f dif: %f", i, j, fHess[i][j], Hline[j], MIP_abs( objFactor*Hline[j] - fHess[i][j] ));
                    
                    if( fHess[i][j] == 0 )
                    {
                        if( !H.hasIndex(i,j) )
                            printf(" Index not described in lagrangean hessian");
                    }
                    
                    printf("\n");
                    
                    answer = false;
                }
            }
            
        }
        
        
        if( checkConstr )
        {
            for(k = 0; k < m; k++)
            {
                
                for(j = 0; j < npoints; j++)
                {
                    p = fgrad[j];
                    
                    for(w = 0; w < n; w++)
                        p[w] = 0.0;
                }
                
                
                for(j = 0; j < npoints; j++)
                    Jacs[j].copyRowTo(k,  fgrad[j]); //Jacs[j][k].copyTo( fgrad[j] ); //accessing rows[k] from Jacs[j]
                
                
                lamb = lambda ? lambda[i] : 1.0;
                
                
                MIP_AppDerivativeByCenteredDif(n, npoints, fgrad, h, Hline);
                
                //we only check lower triangle
                for(j = 0; j <= i; j++)
                {
                    if( MIP_abs(lamb*Hline[j] - gHess[k][i][j]) > tolerance )
                    {
                        printf("constr %d hessian (%d, %d) exact: %f app: %f dif: %f", k, i, j, gHess[k][i][j], Hline[j], MIP_abs( lamb*Hline[j] - gHess[k][i][j] ) );
                        
                        if( gHess[k][i][j] == 0 )
                        {
                            if( !H.hasIndex(i, k) )
                                printf(" Index not described in lagrangean hessian");
                        }
                        
                        printf("\n");
                    
                        answer = false;
                    }
                }
            }
            
            
        }
        
        
        //getchar();
        
        
        xp[i] = x[i];
    }
    
    
    
    
    
    
    
    
    
    code = 0;
    
termination:
    
    if( checkObj || checkConstr )
    {
        nlEval->finalize(1, n, m, J.getNumberOfElements(), lagH.getNumberOfElements());
    }
    
    if(xp)			free(xp);
    if(myLambda)	free(myLambda);
    if(Hline)		free(Hline);
    
    MIP_freeMatrix( fgrad );
    MIP_freeMatrix( fHess );
    
    if(Jacs)	delete[] Jacs;
    
    if( gHess )
    {
        for(i = 0; i < m; i++)
        {
            if( gHess[i] )
                MIP_freeMatrix( gHess[i] );
            else
                break;
        }
        
        free( gHess );
    }
    
    
    
    if(code != 0)
        answer = false;
    
    return code;
}




int MIP_MINLPProb::copyProblemFrom(const MIP_MINLPProb& other)
{
    int r, code;
    deallocateMatrices();
    
    
    r = addVariables( other.n );
    MIP_IFERRORGOTOLABEL(r, code, r, termination);
    
    r = addConstraints( other.m );
    MIP_IFERRORGOTOLABEL(r, code, r, termination);
    
    
    objLinearPart = other.objLinearPart;
    nI = other.nI;
    hasNlObj = other.hasNlObj;
    hasNlConstrs = other.hasNlConstrs;
    
    
    d = other.d;
    objFactor = other.objFactor;
    objValue = other.objValue;
    nlEval = other.nlEval;
    
    
    MIP_copyArray(n, other.lx, lx);
    MIP_copyArray(n, other.ux, ux);
    
    MIP_copyArray(n, other.xprior, xprior);
    MIP_copyArray(n, other.xtype, xtype);
    MIP_copyArray(n, other.x, x);
    
    if( objLinearPart )
        MIP_copyArray(n, other.c, c);
    
    
    MIP_copyArray(m, other.nlConstr, nlConstr);
    MIP_copyArray(m, other.lc, lc);
    MIP_copyArray(m, other.uc, uc);
    
    
    for(int i = 0; i < m; i++)
    {
        if( other.QC[i].getNumberOfElements() > 0 )
        {
            r = QC[i].copyMatrixFrom( other.QC[i] );
            MIP_IFERRORGOTOLABEL(r, code, MIP_UNDEFINED_ERROR, termination);
            
        }
    }
    
    
    if( other.Q.getNumberOfElements() > 0 )
    {
        r = Q.copyMatrixFrom( other.Q );
        MIP_IFERRORGOTOLABEL(r, code, MIP_UNDEFINED_ERROR, termination);
    }
    
    if( other.A.getNumberOfElements() > 0 )
    {
        r = A.copyMatrixFrom( other.A );
        MIP_IFERRORGOTOLABEL(r, code, MIP_UNDEFINED_ERROR, termination);
    }
    
    if( other.J.getNumberOfElements() > 0 )
    {
        r = J.copyMatrixFrom( other.J );
        MIP_IFERRORGOTOLABEL(r, code, MIP_UNDEFINED_ERROR, termination);
    }
    
    
    if( other.lagH.getNumberOfElements() > 0 )
    {
        r = lagH.copyMatrixFrom( other.lagH );
        MIP_IFERRORGOTOLABEL(r, code, MIP_UNDEFINED_ERROR, termination);
    }
    
    
    code = 0;
    
termination:

    return code;
}




void MIP_MINLPProb::deallocateMatrices()
{
    int i;
    
    MIP_secFree(nlConstr);
    MIP_secFree(xprior);
    MIP_secFree(xtype);
    //MIP_secFree(ctype);
    MIP_secFree(lx);
    MIP_secFree(ux);
    MIP_secFree(c);
    //MIP_secFree(b);
    MIP_secFree(lc);
    MIP_secFree(uc);
    MIP_secFree(x);
    //MIP_secFree(xInit);
    
    if(QC)
    {
        for(i = 0; i < m; i++)
            QC[i].desallocateMemory();
            
        //delete[] QC;
        free(QC);
        QC = NULL;
    }
    
    Q.desallocateMemory();
    J.desallocateMemory();
    A.desallocateMemory();
}



/*template <class myclass>
inline void MIP_shiftUpArray(const int size, myclass *a)
{
    int i;
    
    for(i = 1; i < size; i++)
        a[i-1] = a[i];
} */



int MIP_MINLPProb::deleteConstraints(const int ncons, const int* cons)
{
    int i, j, code;
    bool *remove = NULL;
    //double *auxd;
    //MIP_SparseMatrix *auxM;
    
    if(ncons > m)
    {
        code = MIP_BAD_DEFINITIONS;
        goto desallocate_memory;
    }
    
    
    MIP_calloc(remove, m); //remove = (bool *) calloc(m, sizeof(bool));
    
    if(!remove)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    
    
    for(i = 0; i < ncons; i++)
    {
        if( cons[i] >= m || cons[i] < 0 )
        {
            code = MIP_BAD_DEFINITIONS;
            goto desallocate_memory;
        }
        
        remove[ cons[i] ] = true;
    }
    
    
    
    
    for(i = m-1; i >= 0; i--)
    {
        if( remove[i] )
        {
            QC[i].desallocateMemory();
            
            
            for(j = i + 1; j < m; j++)
            {
                lc[j-1] = lc[j];
                uc[j-1] = uc[j];
                nlConstr[j-1] = nlConstr[j];
                QC[j-1] = QC[j];
            }
            
            m--;
        }
    }
    
    
    /*auxd = (double *) realloc( lc, m*sizeof(double) );
    if(auxd)
        lc = auxd; //if we got an memory error, there is no diference... */
    MIP_realloc(lc, m); //if we got an memory error, there is no diference...
    
    /*auxd = (double *) realloc( uc, m*sizeof(double) );
    if(auxd)
        uc = auxd; //if we got an memory error, there is no diference...*/
    MIP_realloc(uc, m); //if we got an memory error, there is no diference...
    
    /*auxb = (bool *) realloc( nlConstr, m*sizeof(bool) );
    if(auxb)
        nlConstr = auxb; //if we got an memory error, there is no diference... */
    MIP_realloc( nlConstr, m ); //if we got an memory error, there is no diference...
    
    /*auxM = (MIP_SparseMatrix *) realloc(QC, m*sizeof(MIP_SparseMatrix) );
    if(auxM)
        QC = auxM; //if we got an memory error, there is no diference... */
    MIP_realloc(QC, m); //if we got an memory error, there is no diference...
    
    
    J.removeRows(remove);
    A.removeRows(remove);
    
    
    updateNlConstrsFlag();
    
    code = 0;
    
    
desallocate_memory:
    
    if(remove)	free(remove);
    
    return code;
}



void MIP_MINLPProb::deleteJacobianStructure()
{
    J.deleteStructure();
}


void MIP_MINLPProb::deleteJacobianStructureLine(const int line)
{
    J.deleteRowStructure(line);
}


void MIP_MINLPProb::deleteLagrangianStructure()
{
    lagH.deleteStructure();
}


void MIP_MINLPProb::deleteLagrangianStructureLine(const int line)
{
    lagH.deleteRowStructure(line);
}



int MIP_MINLPProb::getConstraintBounds(const int constrIndex, double& lb, double& ub) const
{
    if( constrIndex < 0 || constrIndex >= m )
        return MIP_INDEX_FAULT;
    
    lb = lc[constrIndex];
    ub = uc[constrIndex];
    
    return 0;
}



void MIP_MINLPProb::getConstraintLowerBounds(double* lbs) const
{
    MIP_copyArray(m, lc, lbs);
}


void MIP_MINLPProb::getConstraintUpperBounds(double* ubs) const
{
    MIP_copyArray(m, uc, ubs);
}



//function get the number of linear, quadratic and nonlinear constraints
void MIP_MINLPProb::getConstraintStatistcs(int* ml, int* mq, int* mnl) const
{
    decltype(m) myml = 0, mymq = 0, mymnl = 0;
    
    for(decltype(m) i = 0; i < m; i++)
    {
        if( nlConstr[i] )
            mymnl++;
        else if( QC[i].getNumberOfElements() > 0 )
            mymq++;
        else
            myml++;
    }
    
    if(ml)
        *ml = myml;
    if(mq)
        *mq = mymq;
    if(mnl)
        *mnl = mymnl;
}


int MIP_MINLPProb::getNumberOfLinearConstraints(void) const
{
    decltype(m) ml = 0;
    
    for(decltype(m) i = 0; i < m; i++)
    {
        if( !nlConstr[i] && QC[i].getNumberOfElements() == 0 )
            ml++;
    }
    
    return ml;
}



int MIP_MINLPProb::getNumberOfNLConstraints(void) const
{
    decltype(m) i, mnl = 0;
    
    for(i = 0; i < m; i++)
    {
        if( nlConstr[i] )
            mnl++;
    }
    
    return mnl;
}



int MIP_MINLPProb::getNumberOfNLEqualityConstraints(void) const
{
    decltype(m) i, r = 0;
    
    
    for(i = 0; i < m; i++)
    {
        if( nlConstr[i] && (lc[i] == uc[i]) )
            r++;
    }
    
    return r;
}


int MIP_MINLPProb:: getNumbersOfNLConstraints( int &nNLEqualityConstraints, int &nNLInequalityConstraints, int &nNLFreeConstraints) const
{
    nNLEqualityConstraints = 0;
    nNLInequalityConstraints = 0;
    nNLFreeConstraints = 0;
    
    for(decltype(m) i = 0; i < m; i++)
    {
        if( nlConstr[i] )
        {
            if( lc[i] == uc[i] )
                nNLEqualityConstraints++;
            else if( lc[i] <= -MIP_INFINITY && uc[i] >= MIP_INFINITY)
                nNLFreeConstraints++;
            else
                nNLInequalityConstraints++;
        }
    }
    
    return 0;
}


int MIP_MINLPProb::getNumberOfQuadMatricesInConstrs( void) const
{
    decltype(m) mq = 0;
    
    for(decltype(m) i = 0; i < m; i++)
    {
        if( QC[i].getNumberOfElements() > 0 )
            mq++;
    }
    
    return mq;
}


int MIP_MINLPProb::getNumberOfQuadConstraints() const
{
    decltype(m) mq = 0;
    
    for(decltype(m) i = 0; i < m; i++)
    {
        if( !nlConstr[i] && QC[i].getNumberOfElements() > 0)
            mq++;
    }
    
    return mq;
}


int MIP_MINLPProb::getNumberOfQuadEqualityConstraints() const
{
    decltype(m) mq = 0;
    
    for(decltype(m) i = 0; i < m; i++)
    {
        if( !nlConstr[i] && QC[i].getNumberOfElements() > 0 && lc[i] == uc[i])
            mq++;
    }
    
    return mq;
}


int MIP_MINLPProb::getNumberOfQuadCoefsInConstr(const int constrIndex, int& nzs) const
{
    if( constrIndex < 0 || constrIndex >= m )
        return MIP_INDEX_FAULT;
    
    nzs = QC[constrIndex].getNumberOfElements();
    
    return 0;
}




double MIP_MINLPProb::getObjFactor() const
{
    return objFactor;
}


double MIP_MINLPProb::getObjConstant() const
{
    return d;
}



MIP_PROBLEMTYPE MIP_MINLPProb::getProblemType(void) const
{
    bool aux;
    int i;
    
    if( hasNlObj || hasNlConstrs )
    {//nonlinear programming...
        return nI == 0 ? MIP_PT_NLP : MIP_PT_MINLP;
    }
    else
    {
        aux = false;
        for(i = 0; i < m; i++)
        {
            if (QC[i].getNumberOfElements() > 0  )
            {
                aux = true;
                break;
            }
        }
        
        if( aux ) //quadratically constrained programming
            return nI == 0 ? MIP_PT_QCP : MIP_PT_MIQCP;
        else
        {
            if( Q.getNumberOfElements() == 0 )
                return nI == 0 ? MIP_PT_LP : MIP_PT_MILP;
            else
                return nI == 0 ? MIP_PT_QP : MIP_PT_MIQP;
        }
    }
}



bool MIP_MINLPProb::hasConstraintNLTerm(int index) const
{
    if( index < 0 || index >= m )
        return false;
    
    return nlConstr[index];
}



bool MIP_MINLPProb::hasLinCoefObj(void) const
{
    return objLinearPart;
}


bool MIP_MINLPProb::hasNLConstraints(void) const
{
    return hasNlConstrs;
}


bool MIP_MINLPProb::hasNLTerm(void) const
{
    return hasNlConstrs || hasNlObj;
}


bool MIP_MINLPProb::hasObjNLTerm(void) const
{
    return hasNlObj;
}


bool MIP_MINLPProb:: hasQuadMatrixInSomeConstraint(void) const
{
    for(int i = 0; i < m; i++)
    {
        if( QC[i].getNumberOfElements() > 0 )
            return true;
    }
    
    return false;
}



bool MIP_MINLPProb::isBinaryProblem(void) const
{
    int i;
    
    for(i = 0; i < n; i++)
    {
        if( MIP_isIntegerType(xtype[i]) )
        {
            if( -1.0 >= lx[i] || lx[i] >= 2.0 || -1.0 >= ux[i] || ux[i] >= 2.0 )
                return false;
        }
    }
    
    return true;
}



bool MIP_MINLPProb::isConstrValuesFeasible( const double absTol, const double relTol, const double* constrValues ) const
{
    
    for(int i = 0; i < m; i++)
    {
        
        if( lc[i] > -MIP_INFINITY )
        {
            const double ltol = absTol + MIP_abs( relTol *lc[i]);
            
            //note, we use not operator (!) because constrValues can be NAN. In this case, test inside parentehsys would be evaluated like false... We use note to consider this case correctly like infeasible
            if( !(lc[i] <= constrValues[i] + ltol) )
            {
                //std::cout << "Nao e viavel para restricao " << i << " lc: " << lc[i] << " uc: " << uc[i] << " value: " << constrValues[i] << " ltol: " << ltol << std::endl;
                //std::cout << "absTol: " << absTol << " relTol: " << relTol << "\n";
                
                return false;
                break;
            }
        }
        
        
        if( uc[i] < MIP_INFINITY )
        {
            const double utol = absTol + MIP_abs( relTol *uc[i]);
            
            //note, we use not operator (!) because constrValues can be NAN. In this case, test inside parentehsys would be evaluated like false... We use note to consider this case correctly like infeasible
            if( !(constrValues[i] - utol <= uc[i]) )
            {
                //std::cout << "Nao e viavel para restricao " << i << " lc: " << lc[i] << " uc: " << uc[i] << " value: " << constrValues[i] << " utol: " << utol << std::endl;
                //std::cout << "absTol: " << absTol << " relTol: " << relTol << "\n";
                
                return false;
                break;
            }
        }
        
    }
    
    
    return true; //i == m; //if i is equal to m, we have did not break loop above and so, constrValues are feasible
}



int MIP_MINLPProb::isFeasibleToConstraints(const int thnumber, const double* x, const bool newx, const bool* constrEval, const double absTol, const double relTol, bool& answer, double* constrValues) const
{
    bool valuesAlloc = false;
    int r, code;
    
    answer = false;
    
    if(constrValues == NULL)
    {
        MIP_malloc(constrValues, m); //constrValues = (double *) malloc( m * sizeof(double) );
        if(!constrValues)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTMEMERROR;
            #endif
            
            code = MIP_MEMORY_ERROR;
            goto desallocate_memory;
        }
        valuesAlloc = true;
    }
    
    
    r = constraintsEval(thnumber, newx, constrEval, x, constrValues);
    if(r != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        
        code = r;
        goto desallocate_memory;
    }
    
    
    answer = true;
    for(int i = 0; i < m; i++)
    {
        if( constrEval == NULL || constrEval[i] )
        {
            
            if( lc[i] > -MIP_INFINITY )
            {
                const double ltol = absTol + MIP_abs( relTol *lc[i]);
                
                //note, we use not operator (!) because constrValues can be NAN. In this case, test inside parentehsys would be evaluated like false... We use not to consider this case correctly like infeasible
                if( !(lc[i] <= constrValues[i] + ltol) )
                {
                    //printf("Restricao %d nao e viavel. lc: %f ltol: %f value: %f\n", i, lc[i], ltol, constrValues[i] );
                    answer = false;
                    break;
                }
            }
            
            
            if( uc[i] < MIP_INFINITY )
            {
                const double utol = absTol + MIP_abs( relTol *uc[i]);
                
                //note, we use not operator (!) because constrValues can be NAN. In this case, test inside parentehsys would be evaluated like false... We use not to consider this case correctly like infeasible
                if( !(constrValues[i] - utol <= uc[i]) )
                {
                    //printf("Restricao %d nao e viavel. uc: %f utol: %f value: %f\n", i, uc[i], utol, constrValues[i] );
                    answer = false;
                    break;
                }
            }
            
        }
    }
    
    
    
    code = 0;
    
desallocate_memory:
    
    if(valuesAlloc)
        free(constrValues);
    
    return code;
}







int MIP_MINLPProb::getConstraintLinearCoef(const int row, const int col, double& value, bool* inStructure) const
{
    if( row < 0 || row >= m || col < 0 || col >= n )
        return MIP_INDEX_FAULT;
    
    int r;
    
    r = A.getElement(row, col, value);
    
    if( inStructure )
        *inStructure = r == 0;
        
    
    return 0;
}




int MIP_MINLPProb::getConstraintLinearPart(const int row, int *nzs, int* cols, double* values) const
{
    if( row < 0 || row >= m )
        return MIP_INDEX_FAULT;
    
    //A.getRowStructureAndValues(row, nzs, cols, values);
    
    int nzs1, nzs2;
    int r = A.getRowStructure(row, cols, &nzs1);
    r += A.getRowValues(row, values, &nzs2);
    
    if(r)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        return MIP_UNDEFINED_ERROR;
    }
    
    
    if(nzs)
        *nzs=nzs1;
    
    return 0;
}


int MIP_MINLPProb::getConstraintsNonLinearTermFlag( const int constrIndex, bool& flag) const
{
    if( constrIndex < 0 || constrIndex >= m )
        return MIP_INDEX_FAULT;
    
    flag = nlConstr[ constrIndex ];
    
    return 0;
}






int MIP_MINLPProb::getConstraintQuadCoef(const int constrIndex, const int row, const int col, double& value, bool* inStructure) const
{
    if(  row < 0 || row >= n || col < 0 || col >= n )
        return MIP_INDEX_FAULT;
    
    int r;
    
    value  = 0.0;
    if(inStructure)
        *inStructure = false;
    
    
    if( QC[constrIndex].getNumberOfElements() > 0 )
    {
        r = QC[constrIndex].getElement( row, col, value );
        
        if( inStructure )
            *inStructure = r == 0;
    }
    
    return 0;
}



int MIP_MINLPProb::getConstraintQuadCoefMatrix(const int constrIndex, int *nzs, int* rows, int* cols, double* values) const
{
    if( constrIndex < 0 || constrIndex >= m )
        return MIP_INDEX_FAULT;
    
    if( QC[constrIndex].getNumberOfElements() > 0 )
    {
        //nzs = QC[constrIndex].getStructureAndValues( rows, cols, values );
        const int nzs1 = QC[constrIndex].getStructure(rows, cols);
        const int nzs2 = QC[constrIndex].getValues(values);
        
        #if MIP_DEBUG_MODE
            assert(nzs1 == nzs2);
        #endif
        
        if(nzs)
            *nzs = nzs1;
    }
    else
        nzs = 0;
    
    return 0;
}


int MIP_MINLPProb::getConstraintQuadCoefMatrix( const int constrIndex, int *rowStart, int *cols, double *values ) const
{
    if( constrIndex < 0 || constrIndex >= m )
        return MIP_INDEX_FAULT;
    
    if( QC[constrIndex].getNumberOfRows() ==  0 )
    {
        MIP_setAllArray(n+1, rowStart, 0);
    }
    else
    {
        QC[constrIndex].getStructureByCompressedRowFormat(rowStart, cols, values);
    }
    
    return 0;
}


int MIP_MINLPProb::getConstraintQuadCoefMatrixRow(const int constrIndex, const int row, int *nzs, int* cols, double* values) const
{
    if( constrIndex < 0 || constrIndex >= m )
        return MIP_INDEX_FAULT;
    
    if( row < 0 || row >= n )
        return MIP_INDEX_FAULT;
    
    
    if( QC[constrIndex].getNumberOfElements() > 0 )
    {
        int nzs1, nzs2;
        
        //QC[constrIndex].getRowStructureAndValues(row, nzs, cols, values);
        QC[constrIndex].getRowStructure(row, cols, &nzs1);
        QC[constrIndex].getRowValues(row, values, &nzs2);
        
        #if MIP_DEBUG_MODE
            assert(nzs1 == nzs2);
        #endif
        
        if(nzs)
            *nzs = nzs1;
    }
    else
        nzs = 0;
    
    
    return 0;
}


void MIP_MINLPProb::getConstraintsLinearPart(int* nzs, int *rows, int *cols, double *values) const
{
    A.getStructure(rows, cols);
    if(values)
        A.getValues(values);
}


//get linear terms by compressed row
void MIP_MINLPProb::getConstraintsLinearPart(int *rowStart, int *cols, double *values) const
{
    A.getStructureByCompressedRowFormat(rowStart, cols, values);
}


int MIP_MINLPProb::getContinuousIndices( int *inds) const
{
    int nc = 0;
    
    for(int i = 0; i < n; i++)
    {
        if( !MIP_isIntegerType( xtype[i] ) )
        {
            inds[nc] = i;
            nc++;
        }
    }
    
    #if MIP_DEBUG_MODE
        assert( nc == n - nI );
    #endif
    
    return nc;
}


int MIP_MINLPProb::getIntegerIndices( int* inds ) const
{
    unsigned int i;
    int k;
    unsigned int n = this->n;
    
    for(i = k = 0; i < n; i++)
    {
        if( MIP_isIntegerType( xtype[i] ) )
        {
            inds[k] = i;
            k++;
        }
    }
    
    #if MIP_DEBUG_MODE
        assert( k == nI );
    #endif
    
    return k;
}




double MIP_MINLPProb::getIntGap(const double* sol) const
{
    int i;
    double gap = 0.0;
    
    for(i = 0; i < n; i++)
    {
        if( MIP_isIntegerType(xtype[i]) )
            gap +=  MIP_abs( sol[i] - round(sol[i]) );
    }
    
    return gap;
}




int MIP_MINLPProb::getJacobianRowStructure(const int row, int *nzs, int* cols) const
{
    int nzs1;
    
    if( row < 0 || row >= m )
        return MIP_INDEX_FAULT;
    
    const int r = J.getRowStructure( row, cols, &nzs1 );
    if(r != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        return MIP_UNDEFINED_ERROR;
    }
    
    if(nzs)
        *nzs = nzs1;
    
    return 0;
}


int MIP_MINLPProb::getJacobianStructure(int *nzs, int* rows, int* cols) const
{
    int nzs1 = J.getStructure( rows, cols );
    
    if(nzs)
        *nzs = nzs1;
    
    return 0;
}


int MIP_MINLPProb::getJacobianStructure(int* rowStart, int* cols) const
{
    J.getStructureByCompressedRowFormat(rowStart, cols, NULL);
    
    return 0;
}


int MIP_MINLPProb::getLagrangianHessianRowStructure(const int row, int *nzs, int* cols) const
{
    if( row < 0 || row >= n )
        return MIP_INDEX_FAULT;
    
    
    const int r = lagH.getRowStructure( row, cols, nzs );
    
    if(r)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        return MIP_UNDEFINED_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::getLagrangianHessianStructure(int* nzs, int* rows, int* cols) const
{
    int mynzs;
    
    mynzs = lagH.getStructure(rows, cols);
    
    if(nzs)
        *nzs = mynzs;
    
    return 0;
}


int MIP_MINLPProb::getLagrangianHessianStructure(int *rowStart, int *cols) const
{
    lagH.getStructureByCompressedRowFormat(rowStart, cols, NULL);
    
    return 0;
}


MIP_NonLinearEval* MIP_MINLPProb::getNonLinearEvaluationObject(void) const
{
    return nlEval;
}




int MIP_MINLPProb::getNumberOfBinaryVars(void) const
{
    int i, nb = 0;
    
    for(i = 0; i < n; i++)
    {
        if( MIP_isIntegerType( xtype[i] ) && -1.0 < lx[i] && lx[i] <= 1.0 && 0.0 <= ux[i] && ux[i] < 2.0 )
            nb++;
    }
    
    return nb;
}


int MIP_MINLPProb::getNumberOfConstraints(void) const
{
    return m;
}


int MIP_MINLPProb::getNumberOfConstraintQuadCoefMatrixTerms(const int constrIndex, int& nzs) const
{
    if( constrIndex < 0 || constrIndex >= m )
        return MIP_INDEX_FAULT;
    
    nzs = QC[constrIndex].getNumberOfElements();
    
    return 0;
}


int MIP_MINLPProb::getNumberOfIntegerVars(void) const
{
    return nI;
}



int MIP_MINLPProb::getNumberOfJacobianNonZeros(void) const
{
    return J.getNumberOfElements();
}



int MIP_MINLPProb::getNumberOfJacobianNonZerosAtRow( const int row, int& nzs) const
{
    if( row < 0 || row >= m )
        return MIP_INDEX_FAULT;
    
    nzs = J.getNumberOfElementsAtRow(row); //J[row].getNumberOfElements();
    
    return 0;
}


int MIP_MINLPProb:: getNumberOfLagrangianHessianNonZeros( void) const
{
    return lagH.getNumberOfElements();
}


int MIP_MINLPProb:: getNumberOfLagrangianHessianNonZerosAtRow(const int row, int& nzs) const
{
    if( row < 0 || row >= n )
        return MIP_INDEX_FAULT;
    
    nzs = lagH.getNumberOfElementsAtRow(row);  // lagH[row].getNumberOfElements();
    
    return 0;
}



int MIP_MINLPProb::getNumberOfLinearCoefs(void) const
{
    return A.getNumberOfElements();
}


int MIP_MINLPProb::getNumberOfLinearCoefsInConstr(const int constrIndex, int& nzs) const
{
    if( constrIndex < 0 || constrIndex >= m )
        return MIP_INDEX_FAULT;
    
    nzs = A.getNumberOfElementsAtRow(constrIndex);  // A[constrIndex].getNumberOfElements();
    
    return 0;
}



int MIP_MINLPProb::getNumberOfVars(void) const
{
    return n;
}


int MIP_MINLPProb::getNumberOfObjQuadTerms(void) const
{
    return Q.getNumberOfElements();
}


int MIP_MINLPProb::getNumberOfObjQuadCoefMatrixRowTerms(int row, int& nzs) const
{
    if( row < n || row >= n )
        return MIP_INDEX_FAULT;
    
    
    nzs = Q.getNumberOfElementsAtRow(row);  // Q[row].getNumberOfElements();
    
    return 0;
}




int MIP_MINLPProb::getObjQuadCoef(const int row, const int col, double& value, bool* inStructure) const
{
    if( row < 0 || row >= n || col < 0 || col >= n )
        return MIP_INDEX_FAULT;
    
    value = 0.0;
    
    
    
    if( Q.getNumberOfElements() > 0 )
    {
        int r = Q.getElement(row, col, value);
        
        if(inStructure)
            *inStructure = r == 0;
    }
    else
    {
        if(inStructure)
            *inStructure = false;
    }
    
    
    return 0;
}



int MIP_MINLPProb::getObjQuadCoefsMatrixRow(const int row, int *nzs, int* cols, double* values) const
{
    int nzs1, nzs2;
    
    if( row < 0 || row >= n )
        return MIP_INDEX_FAULT;
    
    
    //Q.getRowStructureAndValues( row, nzs, cols, values );
    
    int r = Q.getRowStructure(row, cols, &nzs1);
    r += Q.getRowValues(row, values, &nzs2 );
    
    #if MIP_DEBUG_MODE
        if(r)
        {
            MIP_PRINTERRORNUMBER(r);
            return MIP_UNDEFINED_ERROR;
        }
        
        assert(nzs1 == nzs2);
    #endif
    
    if(nzs)
        *nzs = nzs1;
    
    return 0;
}


void MIP_MINLPProb::getObjQuadCoefsMatrix(int *nzs, int* rows, int* cols, double* values) const
{
    //nzs = Q.getStructureAndValues(rows, cols, values);
    
    const int nzs1 = Q.getStructure(rows, cols);
    const int nzs2 = Q.getValues(values);
    
    #if MIP_DEBUG_MODE
        assert(nzs1 == nzs2);
    #endif
    
    if(nzs)
        *nzs = nzs1;
}


void MIP_MINLPProb::getObjQuadCoefsMatrix(int *rowStart, int *cols, double *values) const
{
    Q.getStructureByCompressedRowFormat(rowStart, cols, values);
}


int MIP_MINLPProb::getObjLinearCoef(const int varIndex, double& value) const
{
    if( varIndex < 0 || varIndex >= n )
        return MIP_INDEX_FAULT;
    
    value = c[varIndex];
    
    return 0;
}


void MIP_MINLPProb::getObjLinearCoefs(double* value) const
{
    MIP_copyArray( n, c, value );
}



int MIP_MINLPProb::getQuadMatrixConstraintInds(int* indices) const
{
    int mq = 0;
    
    
    for(int i = 0; i < m; i++)
    {
        if( QC[i].getNumberOfElements() > 0 )
        {
            indices[mq] = i;
            mq++;
        }
    }
    
    return mq;
}



int MIP_MINLPProb::getReverseIntegerIndices(int* inds ) const
{
    unsigned int i;
    int k;
    unsigned int n = this->n;
    
    for(i = k = 0; i < n; i++)
    {
        if( MIP_isIntegerType( xtype[i] ) )
        {
            inds[i]= k;
            k++;
        }
        else
        {
            inds[i] = -1;
        }
    }
    
    #if MIP_DEBUG_MODE
        assert( k == nI );
    #endif
    
    return k;
}


int MIP_MINLPProb::getVariableBounds(const int varIndex, double& lb, double& ub) const
{
    if( varIndex < 0 || varIndex >= n )
        return MIP_INDEX_FAULT;
    
    lb = lx[varIndex];
    ub = ux[varIndex];
    
    return 0;
}



void MIP_MINLPProb::getVariableLowerBounds(double* lbs) const
{
    MIP_copyArray(n, lx, lbs);
}


void MIP_MINLPProb::getVariableUpperBounds(double* ubs) const
{
    MIP_copyArray(n, ux, ubs);
}


int MIP_MINLPProb::getVariableType(const int index, MIP_VARTYPE& type) const
{
    if(index < 0 || index >= n)
        return MIP_INDEX_FAULT;
    
    type = xtype[index] == MIP_VT_INTEGER ? MIP_VT_INTEGER : MIP_VT_CONTINUOUS;
    
    return 0;
}


bool MIP_MINLPProb::isIntegerSolution(const double *sol, const double int_tol) const
{
    int i;
    
    for(i = 0; i < n; i++)
    {
        if( MIP_isIntegerType(xtype[i]) && MIP_abs( sol[i] - round(sol[i]) ) > int_tol )
        {
            return false;
        }
    }
    
    return true;
}


int MIP_MINLPProb::isLinearConstraint( const int index, bool &response) const
{
    if(index < 0 || index >= m)
        return MIP_INDEX_FAULT;
    
    response = QC[index].getNumberOfElements() == 0 && !nlConstr[index];
    
    return 0;
}


int MIP_MINLPProb::isQuadraticConstraint(const int index, bool &response) const
{
    if(index < 0 || index >= n)
        return MIP_INDEX_FAULT;
    
    response = QC[index].getNumberOfElements() > 0;
    
    return 0;
}



static inline void MIP_printMatrixWithVars( std::ostream &out, const MIP_SparseMatrix &M, const double factor = 1.0)
{
    /*const unsigned int nrows = QC[i].getNumberOfRows();
    
    for(unsigned int j = 0; j < nrows; j++ )
    {
        const unsigned int ncols = QC[i][j].getNumberOfElements();
        
        for(unsigned int k = 0; k < ncols; k++ )
        {
            col = QC[i][j][k].getColumn();
            value = QC[i][j][k].getValue();
            
            //if( value >= 0.0 )
                //cout << "+";
            
            if( j == col )
                out << std::showpos << 0.5*value << std::noshowpos << "x_" << col << "²  ";
            else// if( j > col )
                out << std::showpos << value << std::noshowpos << "x_" << col << "_" << j << "  ";
        }
    }*/
    
    MIP_SparseMatrixIterator it = M.begin(), itend = M.end();
    
    for( ; it != itend; ++it)
    {
        const unsigned int row = (*it).getRow();
        const unsigned int col = (*it).getColumn();
        const double value = factor * (*it).getValue();
        
        if( row == col )
            out << std::showpos << 0.5*value << std::noshowpos << "(x_" << col << ")²  ";
        else// if( j > col )
            out << std::showpos << value << std::noshowpos << "x_" << col << "*x_" << row << "  ";
    }
}


void MIP_MINLPProb::print( std::ostream &out, const bool printDerivativeStructures) const
{
    //int i, aux;
    //unsigned int line, col;
    //double value;
    //MIP_SparseMatrix::iterator itrows;
    //MIP_SparseRow::iterator itcols;
    
    out << "Minimize ";
    
    if(c)
    {
        //out << fixed << setprecision(6);
        
        for(int i = 0; i < n; i++)
        {
            if(c[i] == 0.0)
                continue;
            
            const double coef = objFactor*c[i];
            
            out << std::showpos << coef << std::noshowpos << "x_" << i << "  ";  //printf("%0.2fx_%d  ", c[i], i);
        }
    }
    
    
    /*for( itrows.initialize( &Q ) ; !itrows.isAtEnd() ; itrows++ )
    {
        line = itrows.getPosition(); //get the index of the iterator...
        
        for( itcols.initialize( &(*itrows) ) ; !itcols.isAtEnd() ; itcols++ )
        {
            col = itcols.getCurrentColumn();
            value = itcols.getCurrentValue();
            
            if( value >= 0.0 )
                out << "+";
            
            if( line == col )
                out << 0.5*value << "x_" << col << "²  ";
            else //if( line > col )
                out << value << "x_" << col << "_" << line << "  ";
            
        }
    } */
    
    
    /*for(int i = 0; i < Q.nrows; i++)
    {
        const MIP_SparseElement *colAux = Q.rows[i].columns;
        unsigned int aux = Q.getNumberOfElementsAtRow(i);  //Q.rows[i].nElements;
        
        for(unsigned int j = 0; j < aux; j++)
        {
            if( colAux[j].getValue() >= 0.0 )
                out << "+";
            
            if( i == colAux[j].getColumn() )
                out << 0.5*colAux[j].v << "x_" << i << "²  " ;
            else if( i > colAux[j].col) //in this way, we avoid the repetion of terms in symmetric positions. Note we do not multiply the value by 0.5
                out << colAux[j].v << "x_" << colAux[j].col << "x_" << i << "  ";
        }
    }*/
    MIP_printMatrixWithVars(out, Q, objFactor);
    
    
    
    if(hasNlObj)
        out << "  " << std::showpos << objFactor << "f(x)" << std::noshowpos;
    
    if(d != 0)
        out << "  " << std::showpos << objFactor*d << std::noshowpos;
    
    out << "\nsubject to:\n";
    
    MIP_SparseMatrixIterator itA = A.begin();
    const MIP_SparseMatrixIterator itAend = A.end();
    
    
    
    for(int i = 0; i < m; i++ )
    {
        out << "constr" << i << ":  ";
        
        if( lc[i] > -MIP_INFINITY )
            out << lc[i] << " <=  ";
        
        
        auto Aelement = (*itA);
        
        while( Aelement.getRow() == (unsigned int) i && itA != itAend )
        {
            out << std::showpos << Aelement.getValue() << std::noshowpos << "x_" << Aelement.getColumn() << "  ";
            
            ++itA;
            Aelement = (*itA);
        }
        
        
        if( QC )
        {
            MIP_printMatrixWithVars(out, QC[i], 1.0);
        }
        
        
        if(nlConstr[i])
            out << " + g(x)";
        
        if( uc[i] < MIP_INFINITY )
            out << " <=  " << uc[i];
        
        out << "\n";
        
    }
    
    
    /*for(i = 0; i < m; i++)
    {
        if( lc[i] > -MIP_INFINITY )
            printf("%0.2f <= ", lc[i]);
        
        //if(A)
        //{
            colAux = A.rows[i].columns;
            aux = A.rows[i].nElements;
            
            for(j = 0; j < aux; j++)
            {
                //if( colAux[j].v >= 0.0 )
                    //printf("+");
            
                printf("%+0.2fx_%d  ", colAux[j].v, colAux[j].col);
            }
        //}
    
        if(QC)
        {
            for(j = 0; j < QC[i].nrows; j++)
            {
                colAux = QC[i].rows[j].columns;
                aux = QC[i].rows[j].nElements;
                
                for(k = 0; k < aux; k++)
                {
                    //if(colAux[k].v >= 0.0)
                    //printf("+");
                    
                    if(j == colAux[k].col)
                        printf("%+0.2fx_%d²  ", 0.5*colAux[k].v, j);
                    else if(j > colAux[k].col) //in this way, we avoid the repetion of terms in symmetric positions. Note we do not multiply the value by 0.5
                        printf("%+0.2fx_%dx_%d  ",colAux[k].v, colAux[k].col, j);
                }
            }
        }
        
        if(nlConstr[i])
            printf(" + g(x)");
        
        if( uc[i] < MIP_INFINITY )
            printf(" <= %0.2f", uc[i]);
        
        printf("\n");
    } */
    
    
    
    
    for(int i = 0; i < n; i++)
    {
        //printf("%0.2f <= x_%d <= %0.2f,  ", lx[i], i, ux[i] );
        out << lx[i] << " <= x_" << i << " <= " << ux[i] << ",  "  ;
    }
    out << "\n"; //printf("\n");
    
    if(nI > 0)
    {
        out << "integer variables: ";
        for(int i = 0, aux = 0; i < n; i++)
        {
            if( MIP_isIntegerType(xtype[i]) )
            {
                aux++;
                
                out << "x_" << i;
                
                if(aux < nI)
                {
                    out << ", ";
                }
                else
                {
                    out << "\n";
                    break;
                }
            }
        }
    }
    
    
    {
        int nz = 0;
        
        for(int i = 0; i < n; i++)
        {
            if( xprior[i] )
                nz++;
        }
        
        if(nz > 0)
        {
            out << "variable priorities: ";
        
            for(int i = 0, k = 0; i < n; i++)
            {
                if( xprior[i] )
                {
                    k++;
                    out << "x_" << i << ": " << xprior[i];
                    
                    if(k < nz)
                    {
                        out << ", ";
                    }
                    else
                    {
                        out << "\n";
                        break;
                    }
                }
            }
        }
    }
    
    /*if(xI)
    {
        printf("Initial values:\n");
        
        for(i = 0; i < n; i++)
            printf("x_%d: %0.2f\n", i, xI[i]);
    } */
    
    if(printDerivativeStructures)
    {
        std::cout << "Jacobian matrix:\n";
        J.printSparseMatrix();
        
        std::cout << "Hessian matrix:\n";
        lagH.printSparseMatrix();
    }
    
    out << "\n";
}





int MIP_MINLPProb::readMIQCPModelInFile(const char* fileName)
{
    int i, j, code;
    int n, mqc;
    int nzQObj, nzA, *nzQCons = NULL;
    int *rows = NULL, *cols = NULL;
    double *vals = NULL;
    FILE *inFile = NULL;
    
    
    inFile = fopen(fileName, "r");
    if(!inFile)
    {
        code = MIP_UNDEFINED_ERROR;
        goto desallocate_memory;
    }
    
    fscanf(inFile, "%d%d", &n, &mqc);
    
    i = setParametersAndAllocate(n, mqc);
    if(i != 0)
    {
        code = i;
        goto desallocate_memory;
    }
    
    
    {
        double objConstant;
        
        fscanf(inFile, "%lf", &objConstant);
        setObjConstant(objConstant);
    }
    
    for(i = 0; i < n; i++)
    { 
        double coef;
        fscanf(inFile, "%lf", &coef);
        
        if( coef != 0.0 )
            setObjLinearCoefficient(i, coef);
    }
    
    MIP_malloc(nzQCons, mqc); //nzQCons = (int *) malloc( mqc * sizeof(int) );
    if(!nzQCons)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    fscanf(inFile, "%d", &nzQObj);
    
    j = MIP_max(mqc, 1);
    j = MIP_max(j, nzQObj);
    for(i = 0; i < mqc; i++)
    {
        fscanf(inFile, "%d", &nzQCons[i]);
        j = MIP_max(j, nzQCons[i]);
    }
    
    if( mqc == 0 )
    {
        nzA = 0;
    }
    else
    {
        fscanf(inFile, "%d", &nzA);
        j = MIP_max(j, nzA);
    }
    
    MIP_malloc(rows, j); //rows = (int *) malloc( j * sizeof(rows) );
    MIP_malloc(cols, j); //cols = (int *) malloc( j * sizeof(cols) );
    MIP_malloc(vals, j); //vals = (double *) malloc( j * sizeof(double) );
    
    if(!rows || !cols || !vals)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    //reading quadratic coefficients at obj function
    for(int i = 0; i < nzQObj; i++)
        fscanf( inFile, "%u%u%lf", &rows[i], &cols[i], &vals[i] );
    
    i = setObjQuadCoefsMatrix(nzQObj, rows, cols, vals);
    if(i != 0)
    {
        code = i;
        goto desallocate_memory;
    }
    
    //reading quadratic coefficients at constraints
    for(i = 0; i < mqc; i++)
    {
        for(j = 0; j < nzQCons[i]; j++)
            fscanf( inFile, "%u%u%lf", &rows[j], &cols[j], &vals[j] );
        
        j =  setConstraintQuadCoefsMatrix(i, nzQCons[i], rows, cols, vals);
        if(j != 0)
        {
            code = j;
            goto desallocate_memory;
        }
    }
    
    free(nzQCons);
    nzQCons = NULL;
    
    
    //reading coefficients of A
    for(i = 0; i < nzA; i++)
        fscanf( inFile, "%u%u%lf", &rows[i], &cols[i], &vals[i] );
    
    //i = setNewConstraintsLinearElements(nzA, rows, cols, vals);
    i = setConstraintsLinearPart(nzA, rows, cols, vals);
    if(i != 0)
    {
        code = i;
        goto desallocate_memory;
    }
    
    free(cols);
    cols = NULL;
    
    
    //reading constraints' RHS 
    for(i = 0; i < mqc; i++)
        fscanf(inFile, "%lf", &vals[i]);
    
    i = setConstraintLowerBounds(vals);
    if(i != 0)
    {
        code = i;
        goto desallocate_memory;
    }
    
    
    //reading constraints' RHS 
    for(i = 0; i < mqc; i++)
        fscanf(inFile, "%lf", &vals[i]);
    
    i = setConstraintUpperBounds(vals);
    if(i != 0)
    {
        code = i;
        goto desallocate_memory;
    }
    
    
    //reading lower bounds to variables
    fscanf(inFile, "%d", &j);
    for(i = 0; i < j; i++ )
    {
        fscanf(inFile, "%d%lf", &rows[0], &vals[0]);
        lx[rows[0]] = vals[0];
    }
    
    //reading upper bounds to variables
    fscanf(inFile, "%d", &j);
    for(i = 0; i < j; i++ )
    {
        fscanf(inFile, "%d%lf", &rows[0], &vals[0]);
        ux[rows[0]] = vals[0];
    }
    
    //reading integer variables
    fscanf(inFile, "%d", &j);
    for(i = 0; i < j; i++)
    {
        fscanf(inFile, "%d", &rows[0]);
        rows[0] = setVariableType(rows[0], MIP_VT_INTEGER);
        if(rows[0] != 0)
        {
            code = rows[0];
            goto desallocate_memory;
        }
    }
    
    code = 0;
    
    
desallocate_memory:
    
    if(inFile) fclose(inFile);
    if(rows) free(rows);
    if(cols) free(cols);
    if(vals) free(vals);
    if(nzQCons) free(nzQCons);
    
    return code;
}



int MIP_MINLPProb::removeConstraints(const int nconstrs, const int* indexes)
{
    bool *cons = NULL;
    int r, code, mrem = 0;
    
    
    MIP_calloc(cons, m); //cons = (bool *) calloc( m, sizeof(bool) );
    
    if( !cons )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    for(int i = 0; i < nconstrs; i++)
        cons[ indexes[i] ] = true;
    
    //we perform this check to avoi dproblems if user duplicate indexes...
    for(int i = 0; i < m; i++)
    {
        if( cons[i] )
            mrem++; 
    }
    
    
    #if MIP_DEBUG_MODE
        assert( mrem <= nconstrs );
    #endif
    
    
    //removing lines from sparse matrices...
    
    A.removeRows( cons );
    J.removeRows( cons );
    
    //shifting the variables arrays
    
    r = MIP_shitfValuesInArray(m, nlConstr, cons);
    
    #if MIP_DEBUG_MODE
        assert( r == mrem );
    #endif
    
    MIP_shitfValuesInArray(m, lc, cons);
    MIP_shitfValuesInArray(m, uc, cons);
    
    
    for(int i = 0; i < m; i++)
    {
        if( cons[i] )
            QC[i].desallocateMemory();
    }
    
    MIP_shitfValuesInArray(m, QC, cons);
    
    m -= mrem;
    
    updateNlConstrsFlag();
    
    
    
    code = 0;
    
termination:
    
    if( cons )	free(cons);
    
    return code;
}




int MIP_MINLPProb::removeVars(const int nvars, const int* indexes)
{
    bool *vars = NULL;
    int r, code, nrem = 0;
    
    
    MIP_calloc(vars, n); //vars = (bool *) calloc( n, sizeof(bool) );
    
    if( !vars )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    
    for(int i = 0; i < nvars; i++)
        vars[ indexes[i] ] = true;
    
    
    //we perform this check to avoi dproblems if user duplicate indexes...
    for(int i = 0; i < n; i++)
    {
        if( vars[i] )
            nrem++; 
    }
    
    #if MIP_DEBUG_MODE
        assert( nrem <= nvars );
    #endif
    
    
    //removing lines from sparse matrices...
    
    if( lagH.getNumberOfRows() > 0 )
        lagH.removeRows( vars );
    
    if( Q.getNumberOfRows() > 0 )
        Q.removeRows( vars );
    
    for( int i = 0; i < m; i++ )
    {
        if( QC[i].getNumberOfRows() > 0 )
            QC[i].removeRows( vars );
        
        QC[i].removeCols( vars );
    }
    
    //removing cols from sparse matricex
    Q.removeCols( vars );
    lagH.removeCols( vars );
    J.removeCols( vars );
    A.removeCols( vars );
        
    
    
    
    
    //shifting the variables arrays
    r = MIP_shitfValuesInArray( n, c, vars );
    
    #if MIP_DEBUG_MODE
        assert( r == nrem );
    #endif
    
    MIP_shitfValuesInArray( n, xprior, vars );
    MIP_shitfValuesInArray( n, xtype, vars);
    MIP_shitfValuesInArray( n, lx, vars );
    MIP_shitfValuesInArray( n, ux, vars );
    MIP_shitfValuesInArray( n, x, vars );
    
    
    
    n -= nrem;
    
    updateNumberOfIntegerVars();
    
    
    code = 0;
    
termination:
    
    if( vars )	free(vars);
    
    return code;
}



int MIP_MINLPProb::setConstraintLinearCoefInStructure( const int constrIndex, const int varIndex, const double value)
{
    
    if( constrIndex < 0 || constrIndex >= m || varIndex < 0 || varIndex >= m )
    {
        return MIP_INDEX_FAULT;
    }
    
    bool printError = false;
    
    #if MIP_DEBUG_MODE
        printError=true;
    #endif
    
    
    const int r = A.setElement( constrIndex, varIndex, value, printError);
    if(r != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        return MIP_BAD_DEFINITIONS;
    }
    
    
    return 0;
}




int MIP_MINLPProb::setConstraintsLinearColumn(const int col, const int nzs, const int* rows, const double* values)
{
    int code = 0;
    int *inds = NULL;
    double v, *vals = NULL;
    
    
    MIP_malloc(inds, n); //inds = (int *) malloc( n * sizeof(int) );
    MIP_malloc(vals, n); //vals = (double *) malloc( n * sizeof(double) );
    if( !inds || !vals )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    
    for( int i = 0; i < nzs; i++ )
    {
        
        int ret = A.getElement( rows[i], col, v);
        
        if( ret == 0 )
        {//this element is already in Sparse Matrix (Nice!)
            
            ret = A.setElement( rows[i], col, values[i], true);
            if( ret != 0 )
            {
                #if MIP_DEBUG_MODE
                    MIP_PRINTERRORNUMBER(ret);
                #endif
                assert(false);
            }
            
        }
        else
        {
            //this element is not in Sparse Matrix (crap!)
            
            int rnzs;
            
            //A.getRowStructureAndValues( rows[i], rnzs, inds, vals );
            
            ret = A.getRowStructure( rows[i], inds, &rnzs );
            ret += A.getRowValues(rows[i], vals, NULL);
            
            if( ret != 0 )
            {
                #if MIP_DEBUG_MODE
                    MIP_PRINTERRORNUMBER(ret/2);
                #endif
                code = MIP_MEMORY_ERROR;
            }
            
            
            inds[rnzs] = col;
            vals[rnzs] = values[i];
            
            ret = A.setRowStructure( rows[i], rnzs + 1, inds, vals );
            
            if( ret != 0 )
            {
                #if MIP_DEBUG_MODE
                    MIP_PRINTERRORNUMBER(ret);
                #endif
                code = MIP_MEMORY_ERROR;
            }
        }
    }
    
    
termination:
    
    if( inds )		free(inds);
    if( vals )		free(vals);
    
    return code;
}





int MIP_MINLPProb::setConstraintLinearPart(const int row, const double *a, const int ncols, const double zeroTol)
{
    if( row < 0 || row >= m )
        return MIP_INDEX_FAULT;
    
    const int ret = A.setRowStructureAndValues(row, a, ncols, zeroTol);
    
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}



int MIP_MINLPProb::setConstraintLinearPart(const int row, const int nz, const int* cols, const double* vals)
{
    if( row < 0 || row >= m )
        return MIP_INDEX_FAULT;
    
    //const int ret = A.setRowStructureAndValues(row, nz, cols, vals);
    const int ret = A.setRowStructure(row, nz, cols, vals);
    
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


/*int MIP_MINLPProb::setConstraintLinearPart(const int row, MIP_SparseRow &a)
{
    if( row < 0 || row >= m )
        return MIP_INDEX_FAULT;
    
    const int ret = A.setRowStructureAndValues(row, a);
    
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}*/



int MIP_MINLPProb::setConstraintsLinearPart(MIP_SparseMatrix& A)
{
    if( A.getNumberOfRows() != (unsigned int) m || A.getNumberOfColumns() != (unsigned int) n || A.getSymmetricFlag() )
        return MIP_BAD_DEFINITIONS;
    
    const int ret = this->A.copyMatrixFrom(A);
    
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}



int MIP_MINLPProb::setConstraintsLinearPart(const int nzs, const int* rows, const int* cols, const double* values)
{
    int r;
    
    A.deleteStructure();
    
    r = A.setStructureAndValues(nzs, rows, cols, values);
    
    if( r != 0 )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        return MIP_UNDEFINED_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::setConstraintsLinearPart(const int *rowStart, const int *cols, const double *values)
{
    A.deleteStructure();
    
    const int r = A.setStructureAndValues(rowStart, cols, values);
    
    if( r != 0 )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        return MIP_UNDEFINED_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::setConstraintNonLinearTermFlag(const int index, const bool answer)
{
    int i;
    if( index < 0 || index >= m )
        return MIP_INDEX_FAULT;
        
    nlConstr[index] = answer;
    
    if(answer)
    {
        hasNlConstrs = true;
    }
    else
    {
        hasNlConstrs = false;
        for(i = 0; i < m; i++)
        {
            if(nlConstr[i])
            {
                hasNlConstrs = true;
                break;
            }
        }
    }
        
    
    return 0;
}


int MIP_MINLPProb::setConstraintsNonLinearTermFlag(const bool *answer)
{
    int i;
    
    for(i = 0; i < m; i++)
    {
        nlConstr[i] = answer[i];
        if(nlConstr[i])
            hasNlConstrs = true;
    }
    
    return 0;
}


int MIP_MINLPProb::setConstraintLowerBound(const int index, const double lb)
{
    if( index < 0 || index >= m )
        return MIP_INDEX_FAULT;
    
    lc[index] = lb;
    
    return 0;
}


int MIP_MINLPProb::setConstraintLowerBounds(const double *lb)
{
    MIP_copyArray(m, lb, lc);
    
    return 0;
}


int MIP_MINLPProb::setConstraintUpperBound(const int index, const double ub)
{
    if( index < 0 || index >= m )
        return MIP_INDEX_FAULT;
    
    uc[index] = ub;
    
    return 0;
}


int MIP_MINLPProb::setConstraintUpperBounds(const double *ub)
{
    MIP_copyArray(m, ub, uc);
    
    return 0;
}



/*int MIP_MINLPProb::setInitialSolutions(const double *xI)
{
    MIP_copyVector(n, xI, xInit);
    return 0;
}


int MIP_MINLPProb::setInitialSolution(const int index, const double value)
{
    if( index < 0 || index >= m )
        return MIP_INDEX_FAULT;
    
    xInit[index] = value;
    return 0;
} */


int MIP_MINLPProb::setJacobianStructure(const int nzs, const int* rows, const int* cols)
{
    int code;
    
    J.deleteStructure();
    
    code = J.setStructureAndValues(nzs, rows, cols);
    if(code != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(code);
        #endif
        code = MIP_BAD_DEFINITIONS;
        goto termination;
    }
    
    
    code = 0; //we do not need it here, but ok...
termination:

    updateNlConstrsFlag();
    
    return code;
}


int MIP_MINLPProb::setJacobianStructure(const int* rowStart, const int* cols)
{
    int code;
    
    J.deleteStructure();
    
    const int r = J.setStructureAndValues(rowStart, cols, NULL);
    
    if(r != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        code = MIP_BAD_DEFINITIONS;
        goto termination;
    }
    
    code = 0;
    
termination:
    
    updateNlConstrsFlag();
    
    return code;
}


void MIP_MINLPProb::updateNlConstrsFlag()
{
    int i;
    
    for(i = 0; i < m; i++)
    {
        if( J.getNumberOfElementsAtRow(i) > 0)
        {
            hasNlConstrs = true;
            nlConstr[i] = true;
        }
    }
}


void MIP_MINLPProb::updateNumberOfIntegerVars()
{
    nI = 0;
    
    for(int i = 0; i < n; i++)
    {
        if( MIP_isIntegerType(xtype[i]) )
            nI++;
    }
}




int MIP_MINLPProb::setJacobianStructure(MIP_SparseMatrix& M)
{
    int code;
    if( M.getSymmetricFlag() || M.getNumberOfRows() != (unsigned int) m || M.getNumberOfColumns() != (unsigned int) n )
        return MIP_BAD_DEFINITIONS; //we do not need goto here...
    
    
    code = J.copyMatrixFrom(M);
    
    if(code != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(code);
        #endif
        
        code = MIP_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    
    code = 0; //we do not need it here, but ok...
    
desallocate_memory:

    updateNlConstrsFlag();
    return code;
}





int MIP_MINLPProb::setJacobianRowStructure(const int row, const int nzs, const int* cols)
{
    int i;
    
    i = J.setRowStructure(row, nzs, cols);
    if( i != 0 )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(i);
        #endif
        return MIP_INDEX_FAULT;
    }
    
    if(nzs > 0)
    {
        hasNlConstrs = true;
        nlConstr[row] = true;
    }
    else if( nzs == 0 )
        updateNlConstrsFlag();
    
    
    return 0;
}



/*int MIP_MINLPProb::setJacobianRowStructure( const int index, MIP_SparseRow &row )
{
    if( index < 0 || index >= m )
        return MIP_INDEX_FAULT;
    
    const int r = J.setRowStructure(index, row);
    
    if(r)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    if(row[index].getNumberOfElements() > 0)
    {
        hasNlConstrs = true;
        nlConstr[row] = true;
    }
    else
        updateNlConstrsFlag();
    
    
    return 0;
}*/


int MIP_MINLPProb::setLagrangianHessianStructure(const int nzs, const int* rows, const int* cols)
{
    lagH.deleteStructure();
    
    const int ret = lagH.setStructureAndValues(nzs, rows, cols);
    
    if(ret)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::setLagrangianHessianStructure(const int *rowStart, const int* cols)
{
    lagH.deleteStructure();
    
    const int ret = lagH.setStructureAndValues(rowStart, cols, NULL);
    
    if(ret)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::setLagrangianHessianStructure(MIP_SparseMatrix& M)
{
    if( !M.getSymmetricFlag() || M.getNumberOfRows() != (unsigned int) n || M.getNumberOfColumns() != (unsigned int) n )
        return MIP_BAD_DEFINITIONS;
    
    const int ret = lagH.copyMatrixFrom(M);
    
    if(ret)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}



/*int MIP_MINLPProb::setLagrangianHessianRowStructure( const int index, MIP_SparseRow &row )
{
    int r;
    
    if( index < 0 || index >= n )
        return MIP_INDEX_FAULT;
    
    r = lagH.setRowStructure(index, row);
    
    if(r)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}*/



int MIP_MINLPProb::setLagrangianHessianRowStructure(const int row, const int nzs, const int* cols)
{
    if( row < 0 || row >= n )
        return MIP_INDEX_FAULT;
    
    const int ret = lagH.setRowStructure(row, nzs, cols);
    
    if(ret)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


#if 0
int MIP_MINLPProb::setNewConstraintsLinearElements(const int nzs, const int* rows, const int* cols, const double* vals)
{
    const int ret = A.setStructureAndValues(nzs, rows, cols, vals);
    
    if(ret)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::setNewConstraintLinearElements(const int index, const double* a)
{
    const int ret = A.setRowStructureAndValues(index, a);
    
    if(ret)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}
#endif


void MIP_MINLPProb::setNonLinearEvaluationObject(MIP_NonLinearEval *nl)
{
    nlEval = nl;
}


void MIP_MINLPProb::setObjConstant(const double constant)
{
    d = constant;
}


void MIP_MINLPProb::setObjFactor(const double factor)
{
    objFactor = factor;
}


int MIP_MINLPProb::setObjLinearCoefficient(const int index, const double c)
{
    if( index < 0 || index >= n )
        return MIP_INDEX_FAULT;
    
    this->c[index] = c;
    
    objLinearPart = true;
    return 0;
}


int MIP_MINLPProb::setObjLinearCoefficients(const double *c)
{
    MIP_copyArray(n, c, this->c);
    
    objLinearPart = true;
    return 0;
}


int MIP_MINLPProb::setObjLinearCoefficients(const int nzs, const int *cols, const double *c)
{
    int i;
    
    for(i = 0; i < nzs; i++)
        this->c[ cols[i] ] = c[i];
    
    objLinearPart = true;
    return 0;
}


void MIP_MINLPProb::setObjNonLinearTermFlag(const bool answer)
{
    hasNlObj = answer;
}


int MIP_MINLPProb::setObjQuadCoefsMatrix(const int nzs, const int* rows, const int* cols, const double* vals)
{
    int ret;
    
    if( Q.getNumberOfRows() == 0 )
    {
        //ret = Q.allocateSparseRows(n);
        ret = Q.initialize(n, n, true);
        if(ret != 0)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTERRORNUMBER(ret);
            #endif
            return MIP_MEMORY_ERROR;
        }
    }
    else
        Q.deleteStructure();
    
    
    ret = Q.setStructureAndValues(nzs, rows, cols, vals);
    
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::setObjQuadCoefsMatrix(const int *rowStart, const int *cols, const double *vals)
{
    if( Q.getNumberOfRows() == 0 )
    {
        //const int ret = Q.allocateSparseRows(n);
        const int ret = Q.initialize(n, n, true);
        if(ret != 0)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTERRORNUMBER(ret);
            #endif
            return MIP_MEMORY_ERROR;
        }
    }
    else
        Q.deleteStructure();
    
    
    const int ret = Q.setStructureAndValues(rowStart, cols, vals);
    
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::setObjQuadCoefsMatrix( const double *M, const double zeroTol, const int nrows, const int ncols )
{
    const int ret = Q.setStructureAndValues(M, zeroTol, nrows, ncols);
    
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::setObjQuadCoefsMatrix(const MIP_SparseMatrix& M)
{
    if( !M.getSymmetricFlag() || M.getNumberOfRows() != (unsigned int) n || M.getNumberOfColumns() != (unsigned int) n )
        return MIP_BAD_DEFINITIONS;
    
    
    const int ret = Q.copyMatrixFrom(M);
    
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}



int MIP_MINLPProb::setObjQuadCoefsMatrixRow(const int row, const int nzs, const int* cols, const double* values)
{
    if( row < 0 || row >= n )
        return MIP_INDEX_FAULT;
    
    if( nzs < 0 || nzs > row + 1) //just lower triangle
        return MIP_BAD_DEFINITIONS;
    
    
    for(int i = 0; i < nzs; i++)
    {
        if( cols[i] < 0 || cols[i] > row )
            return MIP_INDEX_FAULT;
    }
    
    const int r = Q.setRowStructure( row, nzs, cols, values );
    if( r != 0 )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    return 0;
}




int MIP_MINLPProb::setParametersAndAllocate(const int n, const int m)
{
    initialize();
    int r;
    
    if( n > 0 )
    {
        r = addVariables(n);
        if( r != 0)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
    }
    
    if(m > 0)
    {
        r = addConstraints(m);
        if( r != 0)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
    }
    
    
    /*this->n = n;
    this->m = m;
    
    allocateMatrices(); */
    
    return 0;
}


int MIP_MINLPProb:: setConstraintQuadCoefsMatrix( const int index, const int nzs, const int* rows, const int* cols, const double* vals)
{
    int ret;
    
    if( index < 0 || index >= m )
        return MIP_INDEX_FAULT;
    
    
    if( QC[index].getNumberOfRows() == 0 )
    {
        ret = QC[index].initialize(n, n, true);
        
        if(ret != 0)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTERRORNUMBER(ret);
            #endif
            return MIP_MEMORY_ERROR;
        }
    }
    else
        QC[index].deleteStructure();
    
    
    ret = QC[index].setStructureAndValues(nzs, rows, cols, vals);
    
    
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb:: setConstraintQuadCoefsMatrix(const int index, const int* rowStart, const int* cols, const double* vals)
{
    if( index < 0 || index >= m )
        return MIP_INDEX_FAULT;
    
    
    if( QC[index].getNumberOfRows() == 0 )
    {
        //const int ret = QC[index].allocateSparseRows(n);
        const int ret = QC[index].initialize(n, n, true);
        
        if(ret != 0)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTERRORNUMBER(ret);
            #endif
            return MIP_MEMORY_ERROR;
        }
    }
    else
        QC[index].deleteStructure();
    
    
    const int ret = QC[index].setStructureAndValues( rowStart, cols, vals);
    
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::setConstraintQuadCoefsMatrix(const int index, const MIP_SparseMatrix& M)
{
    if( !M.getSymmetricFlag() || M.getNumberOfRows() != (unsigned int) n || M.getNumberOfColumns() != (unsigned int) n )
        return MIP_BAD_DEFINITIONS;
    
    if( index < 0 || index >= m )
        return MIP_INDEX_FAULT;
    
    
    const int ret = QC[index].copyMatrixFrom(M);
    
    if(ret != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(ret);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::setConstraintQuadCoefsMatrix(const int index, const double* M, const double zeroTol, const int nrows, const int ncols)
{
    if( index < 0 || index >= m )
        return MIP_INDEX_FAULT;
    
    int r;
    
    
    if( QC[index].getNumberOfRows() == 0 )
    {
        //r = QC[index].allocateSparseRows(n);
        r = QC[index].initialize(n, n, true);
        
        if(r != 0)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTMEMERROR;
            #endif
            return MIP_MEMORY_ERROR;
        }
    }
    
    
    r = QC[index].setStructureAndValues(M, zeroTol, nrows, ncols);
    
    if(r != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


int MIP_MINLPProb::setConstraintQuadCoefsMatrixRow( const int constrIndex, const int row, const int nzs, const int* cols, const double* values)
{
    if( constrIndex < 0 || constrIndex >= m || row < 0 || row >= n )
        return MIP_INDEX_FAULT;
    
    
    for(int i = 0; i < nzs; i++)
    {
        if( cols[i] < row )
            return MIP_INDEX_FAULT;
    }
    
    
    
    int r;
    
    if( QC[constrIndex].getNumberOfRows() == 0 )
    {
        //r = QC[constrIndex].allocateSparseRows(n);
        r = QC[constrIndex].initialize(n, n, true);
        
        if(r != 0)
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTERRORNUMBER(r);
            #endif
            return MIP_MEMORY_ERROR;
        }
    }
    
    
    r = QC[constrIndex].setRowStructure( row, nzs, cols, values );
    if(r != 0)
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTERRORNUMBER(r);
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    
    return 0;
}




int MIP_MINLPProb::setVariableLowerBounds(const int n, const double *lb)
{
    if( n < 0 || n > this->n )
        return MIP_INDEX_FAULT;
    
    MIP_copyArray(n, lb, lx);
    return 0;
}


int MIP_MINLPProb::setVariableLowerBound(const int index, const double lb)
{
    if( index < 0 || index >= n )
        return MIP_INDEX_FAULT;
    
    lx[index] = lb;
    
    return 0;
}



int MIP_MINLPProb::setVariableLowerBounds(const int ninds, const int *inds, const double *lb)
{
    for(int i = 0; i < ninds; i++)
    {
        if( inds[i] < 0 || inds[i] >= n )
            return MIP_INDEX_FAULT;
    }
    
    for(int i = 0; i < ninds; i++)
        lx[ inds[i] ] = lb[i];
    
    return 0;
}



int MIP_MINLPProb::setVariableUpperBounds(const int n, const double *ub)
{
    if( n < 0 || n > this->n )
        return MIP_INDEX_FAULT;
    
    MIP_copyArray(n, ub, ux);
    return 0;
}


int MIP_MINLPProb::setVariableUpperBound(const int index, const double ub)
{
    if( index < 0 || index >= n )
        return MIP_INDEX_FAULT;
    
    ux[index] = ub;
    
    return 0;
}


int MIP_MINLPProb::setVariableUpperBounds(const int ninds, const int *inds, const double *ub)
{
    for(int i = 0; i < ninds; i++)
    {
        if( inds[i] < 0 || inds[i] >= n )
            return MIP_INDEX_FAULT;
    }
    
    for(int i = 0; i < ninds; i++)
        ux[ inds[i] ] = ub[i];
    
    return 0;
}


int MIP_MINLPProb::setVariablePriority(const int index, const int priority)
{
    if( index < 0 || index >= n )
        return MIP_INDEX_FAULT;
    
    xprior[index] = n;
    
    return 0;
}


int MIP_MINLPProb::setVariableType(const int index, const MIP_VARTYPE type)
{
    if( index < 0 || index >= n )
        return MIP_INDEX_FAULT;
    
    if( MIP_isValidVarType(type) == false )
        return MIP_BAD_DEFINITIONS;
    
    
    if( MIP_isIntegerType( xtype[index] ) )
        nI--;
    
    xtype[index] = type;
    
    if( MIP_isIntegerType(type) )
        nI++;
    
    
    return 0;
}


int MIP_MINLPProb::setVariableTypes(const MIP_VARTYPE *types)
{
    int i;
    
    for(i = 0; i < n; i++)
    {
        if( MIP_isValidVarType(types[i]) == false )
            return MIP_BAD_DEFINITIONS;
    }
    
    nI = 0;
    
    for(i = 0; i < n; i++)
    {
        xtype[i] = types[i];
        
        if( MIP_isIntegerType(types[i]) )
            nI++;
    }
    
    return 0;
}


//if constrEval is NULL, we evaluate all constraints having nonlinear terms...
int MIP_MINLPProb::nlObjAndConstraintsEval(const bool evalObj, const bool evalConstrs, const int threadnumber, bool newx, const bool* constrEval, const double* x, double& objValue, double* constrValues) const
{
    const bool *cEval = constrEval ? constrEval : nlConstr;
    
    int r = 0;
    
    if( evalObj )
    {
        r = nlEval->eval_nl_obj_part( threadnumber, n, newx, x, objValue);
        if( r != 0 )
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTCALLBACKERRORNUMBER(r);
            #endif
        }
        newx = false;
    }
    
    
    if(evalConstrs)
    {
        const int ret = nlEval->eval_nl_constrs_part( threadnumber, n, m, newx, cEval, x, constrValues);
        if( ret != 0 )
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTCALLBACKERRORNUMBER(ret);
            #endif
        }
        
        r += ret;
    }
    
    
    return r == 0 ? 0 : MIP_CALLBACK_FUNCTION_ERROR;
}




int MIP_MINLPProb::objEval(const int threadnumber, const bool newx, const double* x, double& value,  double aditionalNlObjFactor) const
{
    double v;
    
    
    if( hasNlObj )
    {
        int r = nlEval->eval_nl_obj_part(threadnumber, n, newx, x, value);
        if(r != 0)
        {
            #if MIP_DEBUG_MODE && MIP_PRINT_NLP_CALLBACK_FUNCTION_ERROR
                MIP_PRINTCALLBACKERRORNUMBER(r);
            #endif
            return r;
        }
        
        if( aditionalNlObjFactor != 1.0 )
            value *= aditionalNlObjFactor;
        
        value += d;
    }
    else
    {
        value = d;
    }
    
    
    /*printf("Solucao para linearizacao: \n");
    for(int i = 0; i < n; i++)
        printf("x[%d]: %f\n", i, x[i]); */
    
    
    if( objLinearPart )
    {
        v = 0.0; //we do a separated sum to minimize numerical errors
        for(int i = 0; i < n; i++)
            v += x[i]*c[i];
        
        //printf("Obj Linear part: %f\n", v);
        value += v ;
    }
    
    
    if( Q.getNumberOfElements() > 0 )
    {
        value += Q.quadraticEvaluation(x);
        //printf("Obj Quadratic part: %f\n", Q.quadraticEvaluation(x) );
    }
    
    value = objFactor*value;
    
    return 0;
}


int MIP_MINLPProb::objGradEval(const int threadnumber, const bool newx, const double *x, double *values, double additionalNlObjFactor) const
{
    int r;
    
    
    if( hasNlObj )
    {
        r = nlEval->eval_grad_nl_obj_part(threadnumber, n, newx, x, values);
        
        if(r != 0)
        {
            #if MIP_DEBUG_MODE && MIP_PRINT_NLP_CALLBACK_FUNCTION_ERROR
                MIP_PRINTCALLBACKERRORNUMBER(r);
            #endif
            return r;
        }
        
        if( additionalNlObjFactor != 1.0 )
            MIP_multiplyAllArray(n, values, additionalNlObjFactor);
        
        if( objLinearPart )
        {
            MIP_accumulateArray(n, c, values);
        }
        
    }
    else
    {
        if( objLinearPart )
            MIP_copyArray(n, c, values);
        else
            MIP_setAllArray<double>(n, values, 0.0);
    }
    
    
    if( Q.getNumberOfElements() > 0 )
        //accumlating the gradient of 0.5*xQx in values
        Q.quadraticGradientEvaluation(x, values, true);
    
    
    if(objFactor != 1.0)
    {
        /*for(int i = 0; i < n; i++)
            values[i] = objFactor * values[i];*/
        
        MIP_multiplyAllArray(n, values, objFactor);
    }
    
    
    return 0;
}


//if constrEval is NULL, we evaluate all constraints
int MIP_MINLPProb::constraintsEval(const int threadnumber, const bool newx, const bool* constrEval, const double* x, double* values) const
{
    
    if( hasNlConstrs )
    {
        //if pCEval points to nlConstr, we evaluate all nonlinear constraints here
        const bool *pCEval = constrEval ? constrEval : nlConstr;
        
        int r = nlEval->eval_nl_constrs_part(threadnumber, n, m, newx, pCEval, x, values);
        MIP_IFCALLBACKERRORRETURN(r);
        
        if(r != 0)
        {
            #if MIP_DEBUG_MODE && MIP_PRINT_NLP_CALLBACK_FUNCTION_ERROR
                MIP_PRINTCALLBACKERRORNUMBER(r);
            #endif
            return MIP_CALLBACK_FUNCTION_ERROR;
        }
    }
    
    
    for(int i = 0; i < m; i++)
    {
        if(!constrEval || constrEval[i])
        {
            if( !nlConstr[i] )
                values[i] = 0.0; //initializing the value with zeros...
            
            if( QC[i].getNumberOfElements() > 0 )
                values[i] += QC[i].quadraticEvaluation(x);
            
            values[i] += A.rowEvaluation(i, x);
        }
        else
            values[i] = 0.0;
    }
    
    
    return 0;
}


int MIP_MINLPProb::constraintLinearPartEvaluation(const int nIndices, const int *constrIndices, const double *x, double *values )
{
    for(int i = 0; i < nIndices; i++)
    {
        const auto ind = constrIndices[i];        
        values[i] = A.rowEvaluation(ind, x);
    }
    
    return 0;
}





//that function just call the Nonlinear evaluation object to calculate the Jacobian. Maybe we do not need it, but in the filter will be easier do alterations... 
int MIP_MINLPProb::nlJacobianEval(const int threadnumber, const bool newx, const bool *constrEval, const double *x, MIP_SparseMatrix& jacobian) const
{
    if( J.getNumberOfElements() > 0 )
    {
        const int r = nlEval->eval_grad_nl_constrs_part(threadnumber, n, m, J.getNumberOfElements(), newx, constrEval, x, jacobian);
        
        if(r != 0)
        {
            #if MIP_DEBUG_MODE && MIP_PRINT_NLP_CALLBACK_FUNCTION_ERROR
                MIP_PRINTCALLBACKERRORNUMBER(r);
            #endif
            return r;
        }
    }
    
    return 0;
}








//that function just call the Nonlinear evaluation object to calculate the hessian of the lagrangian. Maybe we do not need it, but in the filter will be easier do alterations... 
int MIP_MINLPProb::nlpHessianEval(const int threadnumber, const bool newx, const double* x, const double objFactor, const double* lambda, MIP_SparseMatrix& hessian) const
{
    if( lagH.getNumberOfElements() > 0 && (hasNlObj || hasNlConstrs) )
    {
        const int r = nlEval->eval_hessian_nl_lagran_part(threadnumber, n, m, lagH.getNumberOfElements(), newx, x, this->objFactor* objFactor, lambda, hessian);
        
        if(r != 0)
        {
            #if MIP_DEBUG_MODE && MIP_PRINT_NLP_CALLBACK_FUNCTION_ERROR
                MIP_PRINTCALLBACKERRORNUMBER(r);
            #endif
            return r;
        }
    }
    
    return 0;
}


static inline void MIP_printMatrixToAMPLFile(ofstream &fout, const MIP_SparseMatrix &M, const double factor, int &nTermsInLine, const int maxTermsInLine )
{
    /*for(unsigned int i = 0; i < Q.getNumberOfRows(); i++ )
    {
        MIP_SparseRow &row = Q[i];
        unsigned int rne = row.getNumberOfElements();
        
        for(unsigned int j = 0; j < rne; j++)
        {
            int col = row[j].getColumn();
            double v = row[j].getValue();
            
            double f = ( (col == i) ? 0.5 : 1.0 ) * objFactor;
            
            fout << " " << showpos << f*v << noshowpos  << "*x" << i << "*x" << col;
            
            nTermsInLine++;
            if( nTermsInLine >= maxTermsInLine )
            {
                nTermsInLine = 0;
                fout << "\n";
            }
        }
    }*/
    
    
    for(MIP_SparseMatrixIterator it=M.begin(); it != M.end(); ++it)
    {
        const unsigned int row = (*it).getRow();
        const unsigned int col = (*it).getColumn();
        const double v = (col == row ? 0.5 : 1.0) * factor * (*it).getValue();
        
        fout << " " << showpos << v << noshowpos  << "*x" << row << "*x" << col;
        
        nTermsInLine++;
        if( nTermsInLine >= maxTermsInLine )
        {
            nTermsInLine = 0;
            fout << "\n";
        }
    }
}


int MIP_MINLPProb::writeMIQCPModelInAMPLFile(const char* outFileName, const char* solverOption)
{
    //unsigned int rne;
    int nTermsInLine = 0, code;
    ofstream fout;
    
    const int maxTermsInLine = 20;
    
    
    fout.open( outFileName, ios::out );
    
    if( !fout.is_open() )
    {
        #if MIP_DEBUG_MODE
            std::cerr << MIP_PREPRINT << "Failure to open file '" << outFileName << "'.\n";
        #endif
        
        code = MIP_UNDEFINED_ERROR;
        goto termination;
    }
    
    
    
    //fout << std::showpos; //to show + signal in positive numbers... 
    
    fout << setprecision(15);
    
    fout << "#AMPL model file generated by minlp problem library, developed by "  MIP_AUTHOR ", " MIP_AUTHOR_FILIATION "\n\n\n";
    
    
    for(int i = 0; i < n; i++)
    {
        fout << "var x" << i;
        
        if( xtype[i] == MIP_VT_INTEGER )
            fout << " integer";
        
        if( lx[i] > -MIP_INFINITY )
            fout << " >= " << lx[i];
        if( ux[i] < MIP_INFINITY )
            fout << " <= " << ux[i];
        
        fout << ";\n";
    }
    
    fout << endl;
    
    if( hasNlObj )
        fout << "#objective function still has a nonlinear part, but, unfortunatelly, this part couldnot be represented in this file.\n\n";
    
    fout << "minimize obj:";
    
    if( objLinearPart )
    {
        unsigned int nTerms = 0;
        nTermsInLine = 0;
        
        for(int i = 0; i < n; i++)
        {
            if( c[i] != 0.0 )
            {
                fout << " " << std::showpos << objFactor*c[i] << "*x" << std::noshowpos << i;
                
                nTermsInLine++;
                nTerms++;
                if( nTermsInLine >= maxTermsInLine )
                {
                    nTermsInLine = 0;
                    fout << "\n";
                }
            }
        }
        
        if(nTerms == 0)
            objLinearPart = false;
    }
    
    
    if( Q.getNumberOfElements() > 0 )
    {
        MIP_printMatrixToAMPLFile(fout, Q, objFactor, nTermsInLine, maxTermsInLine);
    }
    
    if( d != 0 || (!objLinearPart && Q.getNumberOfElements() == 0) )
        fout << " " << showpos << objFactor*d << noshowpos;
    
    fout << ";\n\n";
    
    if( hasNlConstrs )
    {
        int mnl = getNumberOfNLConstraints();
        
        fout << "#constraint" << (mnl > 1 ? "s" : "");
        
        for(int i = 0; i < m; i++)
        {
            if( nlConstr[i] )
                fout << " " << i;
        }
        
        fout << " " << (mnl > 1 ? "have" : "has") << " nonlinear components, but, unfortunatelly, those nonlinear components cannot be represented in this file.\n";
    }
    
    
    for(int i = 0; i < m; i++)
    {
        fout << "subject to cons" << i << ": ";
        
        if( lc[i] > -MIP_INFINITY )
            fout << lc[i] << " <= ";
        
        
        nTermsInLine = 0;
        
        /*MIP_SparseRow &arow = A[i];
        rne = arow.getNumberOfElements();
        
        for(unsigned int j = 0; j < rne; j++)
        {
            fout << " " << showpos << arow[j].getValue() << "*x" << noshowpos << arow[j].getColumn();
            
            nTermsInLine++;
            if( nTermsInLine >= maxTermsInLine )
            {
                nTermsInLine = 0;
                fout << "\n";
            }
        }*/
        
        for(MIP_SparseMatrixIterator it = A.begin(i); it != A.end(i); ++it)
        {
            fout << " " << showpos << (*it).getValue() << "*x" << noshowpos << (*it).getColumn();
            
            nTermsInLine++;
            if( nTermsInLine >= maxTermsInLine )
            {
                nTermsInLine = 0;
                fout << "\n";
            }
        }
        
        
        if( QC[i].getNumberOfElements() > 0 )
        {
            
            /*for( int k= 0; k < n; k++ )
            {
                MIP_SparseRow &qrow = QC[i][k];
                rne = qrow.getNumberOfElements();
                
                for(unsigned int j = 0; j < rne; j++)
                {
                    int col = qrow[j].getColumn();
                    double v = qrow[j].getValue();
                    
                    double f = (col == k) ? 0.5 : 1.0 ;
                    
                    fout << " " << showpos << f*v << noshowpos <<  "*x" << k << "*x" << col;
                    
                    nTermsInLine++;
                    if( nTermsInLine >= maxTermsInLine )
                    {
                        nTermsInLine = 0;
                        fout << "\n";
                    }
                }
            }*/
            
            MIP_printMatrixToAMPLFile(fout, QC[i], 1, nTermsInLine, maxTermsInLine);
        }
        
        
        if( A.getNumberOfElementsAtRow(i) == 0 && QC[i].getNumberOfElements() == 0 ) //if( A[i].getNumberOfElements() == 0 && QC[i].getNumberOfElements() == 0 )
            fout << " 0";
        
        
        if( uc[i] < MIP_INFINITY )
            fout << " <= " << uc[i];
        
        
        fout << ";\n";
    }
    
    
    if(solverOption)
        fout << "\n\noption solver " << solverOption << ";\n\nsolve;\n\n";
    
    
    code = 0;
    
termination:
    
    if( fout.is_open() ) 
        fout.close();
    
    return code;
}



static inline void MIP_printMatrixToGAMSFile(ofstream &fout, const MIP_SparseMatrix &M, const double factor, int &valuesPerRowCounter, const int maxValuesPerRow)
{
    //const unsigned int nzrows = M.getNumberOfNonZeroRows();
    
    for(MIP_SparseMatrixIterator it = M.begin(); it != M.end(); ++it)
    {
        const unsigned int row = (*it).getRow();
        const unsigned int col = (*it).getColumn();
        double v = (*it).getValue() * (col == row ? 0.5*factor : factor);
        
        fout <<  " " << showpos << v << noshowpos << "*x" << row << "*x" << col;
        
        if( ++valuesPerRowCounter >= maxValuesPerRow )
        {
            fout << "\n";
            valuesPerRowCounter = 0;
        }
    }
    
    /*for(unsigned int i = 0; i < Q.getNumberOfRows(); i++)
    {
        MIP_SparseRow &row = Q[i];
        unsigned int rne = row.getNumberOfElements(); 
        
        
        for(unsigned int j = 0; j < rne; j++)
        {
            int col = row[j].getColumn();
            double v = row[j].getValue() * (col == i ? 0.5*objFactor : objFactor ) ;
            
            fout <<  " " << showpos << v << noshowpos << "*x" << i << "*x" << col;
            
            if( ++k >= maxValuesPerRow )
            {
                fout << endl;
                k = 0;
            }
        }
    }*/
}



int MIP_MINLPProb::writeMIQCPModelInGAMSFile(const char* outFileName, const char* probName, const char* solverOption, const bool optfile, const double optca, const double optcr)
{
    char *problemName = NULL;
    //unsigned int rne;
    int k, code;
    ofstream fout;
    
    const int maxValuesPerRow = 10;
    
    
    MIP_malloc(problemName, (strlen(probName) + 1 )); //problemName = (char *) malloc( (strlen(probName) + 1 )*sizeof(char) );
    if(!problemName)
    {
        code = MIP_MEMORY_ERROR;
        goto termination;
    }
    
    
    k = strlen(probName);
    for(int i = 0; i < k; i++)
    {
        if( probName[i] == '.' || probName[i] == '-' )
            problemName[i] = '_';
        else
            problemName[i] = probName[i];
    }
    problemName[ k ] = '\0';
    
    
    
    
    
    fout.open( outFileName, ios::out );
    
    if( !fout.is_open() )
    {
        #if MIP_DEBUG_MODE
            std::cerr << MIP_PREPRINT  "Failure to open file '" << outFileName << "'" << MIP_GETFILELINE << ".\n";
        #endif
        code = MIP_UNDEFINED_ERROR;
        goto termination;
    }
    
    fout << setprecision(15);
    
    fout << "*GAMS model file generated by minlp problem library, developed by "  MIP_AUTHOR ", " MIP_AUTHOR_FILIATION "\n\n\n";
    
    
    fout << "variables ";
    
    for(int i = 0; i < n; i++)
    {
        if( !MIP_isIntegerType( xtype[i] ) )
            fout << "x" << i << ", ";
    }
    
    fout << "objF;\n\n";
    
    
    if( nI > 0 )
    {
        fout << "integer variables";
        
        for(int i = 0, v = 0; i < n; i++)
        {
            if( MIP_isIntegerType( xtype[i] ) )
            {
                v++;
                fout << " x" << i;
                
                if( v < nI )
                    fout << ",";
            }
        }
        
        fout << ";\n\n";
    }
    
    
    fout << "equations ";
    for( int i = 0; i < m; i++ )
        fout << "cons" << i << ", ";
    
    fout << "objective;\n";
    
    fout << "*objective function: \n";
    
    if( hasNlObj )
        fout << "*objective function still has a nonlinear part, but, unfortunatelly, this part cannot be represented in this file.\n";
    
    fout << "objective.. objF =e=";
    
    k = 0;
    if( objLinearPart )
    {
        
        for(int i = 0; i < n; i++)
        {
            if( c[i] == 0.0 )
                continue;
            
            fout << " " << showpos <<  objFactor*c[i] << noshowpos << "*x" << i;
            
            if( ++k >= maxValuesPerRow )
            {
                fout << "\n";
                k = 0;
            }
        }
    }
    
    if( Q.getNumberOfElements() > 0 )
    {
        MIP_printMatrixToGAMSFile(fout, Q, objFactor, k, maxValuesPerRow);
        /*for(unsigned int i = 0; i < Q.getNumberOfRows(); i++)
        {
            MIP_SparseRow &row = Q[i];
            unsigned int rne = row.getNumberOfElements(); 
            
            
            for(unsigned int j = 0; j < rne; j++)
            {
                int col = row[j].getColumn();
                double v = row[j].getValue() * (col == i ? 0.5*objFactor : objFactor ) ;
                
                fout <<  " " << showpos << v << noshowpos << "*x" << i << "*x" << col;
                
                if( ++k >= maxValuesPerRow )
                {
                    fout << endl;
                    k = 0;
                }
            }
        }*/
    }
    
    
    //we have to guarantee at least the constant will be put on objective function if we have no coefficients...
    if( d != 0 || (objLinearPart == false && Q.getNumberOfElements() == 0) )
        fout << " " << showpos << objFactor*d << noshowpos;
    
    fout << ";\n\n" ;
    
    
    fout << "*constraints:\n";
    if( hasNlConstrs )
    {
        int mnl = getNumberOfNLConstraints();
        
        fout << "*constraint" << (mnl > 1 ? "s" : "");
        
        for(int i = 0; i < m; i++)
        {
            if( nlConstr[i] )
                fout << " " << i;
        }
        
        fout << (mnl > 1 ? "have" : "has") << " nonlinear components, but those nonlinear components cannot be represented in this file.\n"; 
    }
    
    
    for(int i = 0;  i < m; i++)
    {
        fout << "cons" << i << ".. ";
        
        
        if( lc[i] > -MIP_INFINITY )
            fout << lc[i] << " =l= ";
        
        
        //MIP_SparseRow &arow = A[i];
        /*unsigned int rne; //arow.getNumberOfElements();
        int *rcols;
        double *rvalues;
        
        A.getRowPointers(i, rne, rcols, rvalues);
        
        k = 0;
        for(int j = 0; j < rne; j++)
        {
            fout << " " << showpos << rvalues[j] << noshowpos << "*x" << rcols[j];
            
            if(++k >= maxValuesPerRow)
            {
                fout << "\n";
                k = 0;
            }
        }*/
        
        for(MIP_SparseMatrixIterator it = A.begin(i); it != A.end(i); ++it)
        {
            fout << " " << showpos << (*it).getValue() << noshowpos << "*x" << (*it).getColumn();
            
            if(++k >= maxValuesPerRow)
            {
                fout << "\n";
                k = 0;
            }
        }
        
        
        if( QC[i].getNumberOfElements() > 0 )
        {
            MIP_printMatrixToGAMSFile(fout, QC[i], 1.0, k, maxValuesPerRow);
            
            /*for(int j = 0; j < n; j++)
            {
                MIP_SparseRow &qrow = QC[i][j];
                rne = qrow.getNumberOfElements();
                
                for(int p = 0; p < rne; p++)
                {
                    int col = qrow[p].getColumn();
                    double v = qrow[p].getValue() * (j == col ? 0.5 : 1.0) ;
                    
                    fout << " " << showpos << v << noshowpos << "*x" << j << "*x" << col;
                    
                    if(++k >= maxValuesPerRow)
                    {
                        fout << endl;
                        k = 0;
                    }
                    
                }
            }*/
        }
        
        
        if( A.getNumberOfElementsAtRow(i) == 0 && QC[i].getNumberOfElements() == 0 )
            fout << " 0";
        
        
        if( uc[i] < MIP_INFINITY || lc[i] <= -MIP_INFINITY )
            fout << " =l= " << uc[i];
        
        fout << ";" << endl;
    }
    
    fout << endl;
    
    fout << "*Variable bounds:\n";
    
    for(int i = 0; i < n; i++)
    {
        if(lx[i] > -MIP_INFINITY)
            fout << "x" << i << ".lo = " << lx[i] << ";\n";
        
        if(ux[i] < MIP_INFINITY)
            fout << "x" << i << ".up = " << ux[i] << ";\n";
    }
    
    fout << "\n";
    
    fout << "*Options:\n";
    
    if(solverOption)
        fout << "\noption MIQCP = " << solverOption << ";\n";
    
    
    //defining a GAMS program variable called problemName. It is not necessary to run GAMS, but it helps me a lot!
    
    fout << "\n$set problemName '" << problemName << "';\n\n";
    
    
    fout << "model " << problemName << " /all/ ;\n\n";
    
    
    fout << problemName << ".optcr = " << optcr << ";\n";
    fout << problemName << ".optca = " << optca << ";\n";
    
    
    if(optfile)
        fout << problemName << ".optfile=1;\n\n";
    
    
    fout << "solve " << problemName << " using MIQCP minimizing objF;\n\n";
    
    
    
    code = 0;
    
termination:
    
    if( fout.is_open() ) 
        fout.close();
    
    if( problemName )	free(problemName);
    
    
    return code;
}





static inline void MIP_printMatrixMIQCPModelFormat(FILE *outFile, const MIP_SparseMatrix &M, const double factor)
{
    /*for(unsigned int i = 0; i < nrows; i++)
    {
        unsigned int ncols = Q[i].getNumberOfElements();k
        
        for(unsigned int j = 0; j < ncols; j++)
            fprintf( outFile, "%u %u %0.14f \n", i, Q[i][j].getColumn(), Q[i][j].getValue() );
    }*/
    
    for(MIP_SparseMatrixIterator it=M.begin(); it != M.end(); ++it)
        fprintf(outFile, "%u %u %0.14f \n", (*it).getRow(), (*it).getColumn(), (*it).getValue()*factor );
    
    
}



/*
* 
* Min c'x + 0.5 x'Qx + d
* s. t.:
* 		l_c <= A x + 0.5 x'QQ x <= u_c
* 		
* 		l <= x <= u
* 
* format of output file:
* 
* n mqc
* d
* c
* number of nonzeros of Q number of nonzeros of Q_i (i = 0 ... mqc-1) number of nonzeros of A
* Nonzeros of lower triangle of Q (triple sparse)
* for i = 0 ... mqc-1
*     Nonzeros of lower triangle of Q_i (triple sparse) (Q_i is the matrix of constraint i)
* Nonzeros of A (triple sparse)
* signals of all constraints {<, =, >}
* b
* Number of variables having finite lower bound
* finite lower bounds ( <index> <lower bound> )
* Number of variables having finite upper bound
* finite upper bounds ( <index> <upper bound> )
* Number of integer variables
* indexes of integer variables
*  
*/

int MIP_MINLPProb::writeMIQCPModelInFile(const char* fileName)
{
    int aux, code;
    //unsigned int nrows;
    FILE *outFile = NULL;
    
    outFile = fopen(fileName, "w");
    
    if(!outFile)
    {
        code = MIP_UNDEFINED_ERROR;
        goto desallocate_memory;
    }
    
    
    fprintf(outFile, "%u %u\n\n", (unsigned int)n, (unsigned int)m);
    
    fprintf(outFile, "%f\n\n", objFactor*d);
    
    for(int i = 0; i < n; i++)
        fprintf(outFile, "%0.14f ", objFactor*c[i]);
    
    fprintf(outFile, "\n\n");
    
    fprintf(outFile, "%u" , Q.getNumberOfElements() );
    
    
    for(int k = 0; k < m; k++)	
        fprintf(outFile, " %u", QC[k].getNumberOfElements() );
    
    fprintf(outFile, " %u", A.getNumberOfElements() );
    
    fprintf(outFile, "\n\n");
    
    /*nrows = Q.getNumberOfRows();
    for(unsigned int i = 0; i < nrows; i++)
    {
        unsigned int ncols = Q[i].getNumberOfElements();
        
        for(unsigned int j = 0; j < ncols; j++)
            fprintf( outFile, "%u %u %0.14f \n", i, Q[i][j].getColumn(), Q[i][j].getValue() );
    } */
    MIP_printMatrixMIQCPModelFormat(outFile, Q, objFactor);
    
    
    //quadratic coeficients in the cosntraints
    for(int k = 0; k < m; k++)
    {
        //fprintf(outFile, "%d\n", qc->Q[k].nElements );
        
        if( QC[k].getNumberOfElements() > 0)
        {
            fprintf(outFile, "\n");
            MIP_printMatrixMIQCPModelFormat(outFile, QC[k], 1);
        }
        
        /*nrows = QC[k].getNumberOfRows();
        
        for(unsigned int i = 0; i < nrows; i++ )
        {
            unsigned int ncols = QC[k][i].getNumberOfElements();
            
            for(unsigned int j = 0; j < ncols; j++)
                fprintf(outFile, "%u %u %0.14f \n", i,  QC[k][i][j].getColumn(), QC[k][i][j].getValue() );
        }*/
    }
    
    
    
    //linear coeficients matrix in the constraints
    //fprintf(outFile, "%d\n", qc->A->nElements  );
    fprintf(outFile, "\n");
    
    /*for(int i = 0; i < m; i++)
    {
        unsigned int ncols = A[i].getNumberOfElements();
        
        for(unsigned int j = 0; j < ncols; j++)
            fprintf(outFile, "%u %u %0.14f \n", i, A[i][j].getColumn(), A[i][j].getValue() );
    }*/
    MIP_printMatrixMIQCPModelFormat(outFile, A, 1);
    
    fprintf(outFile, "\n");
    
    
    for(int i = 0; i < m; i++)
        fprintf(outFile, "%0.14f ", lc[i]);
    fprintf(outFile, "\n");
    
    for(int i = 0; i < m; i++)
        fprintf(outFile, "%0.14f ", uc[i]);
    
    fprintf(outFile, "\n");
    
    
    
    
    /*variable bounds*/
    //lower bounds...
    aux = 0;
    for(int i = 0; i < n; i++)
    {
        if( lx[i] > -MIP_INFINITY )
            aux++;
    }
    fprintf(outFile, "\n");
    fprintf(outFile, "%u\n", (unsigned int) aux);
    
    for(int i = 0; i < n; i++)
    {
        if( lx[i] > -MIP_INFINITY )
            fprintf(outFile, "%u %0.14f\n", i, lx[i]);
    }
    
    aux = 0;
    for(int i = 0; i < n; i++)
    {
        if( ux[i] < MIP_INFINITY )
            aux++;
    }
    
    fprintf(outFile, "\n");
    fprintf(outFile, "%u\n", (unsigned int) aux);
    
    for(int i = 0;i < n; i++)
    {
        if( ux[i] < MIP_INFINITY )
            fprintf(outFile, "%u %0.14f\n", i, ux[i]);
    }
    
    //number of integer variables
    fprintf(outFile, "\n");
    fprintf(outFile, "%u\n", (unsigned int) nI);
    for(int i = 0; i < n; i++)
    {
        if( MIP_isIntegerType(xtype[i]) )
            fprintf(outFile, "%u ", i);
    }
    
    
    code = 0;
    
desallocate_memory:

    if(outFile)    fclose(outFile);
    
    return code;
}




static inline void MIP_printMatrixToLPFormat( const MIP_SparseMatrix &M, MIP_StreamMBLine &out, const double factor ) 
{
    int size;
    char s[100];
    
    out.print("  + [");
            
    /*for( int k= 0; k < n; k++ )
    {
        MIP_SparseRow &qrow = QC[i][k];
        rne = qrow.getNumberOfElements();
        
        for(unsigned int j = 0; j < rne; j++)
        {
            int col = qrow[j].getColumn();
            double v = qrow[j].getValue();
            
            if(col != k)
                v *= 2;
            
            
            sprintf( s, " %+0.10f x%d", v, k);
            
            if( col == k )
                size = sprintf(s, "%s ^ 2", s);
            else
                size = sprintf(s, "%s * x%d", s, col );
            
            out.print(s, size);
        }
    }*/
    
    for( MIP_SparseMatrixIterator it = M.begin(); it != M.end() ; ++it )
    {
        const unsigned int k = (*it).getRow();
        const unsigned int col = (*it).getColumn();
        double v = (*it).getValue()*factor;
        
        if(col != k)
            v *= 2;
        
        sprintf( s, " %+0.10f x%d", v, k);
        
        if( col == k )
            size = sprintf(s, "%s ^ 2", s);
        else
            size = sprintf(s, "%s * x%d", s, col );
        
        out.print(s, size);
    }
    
    
    out.print(" ]/2 ");
}



int MIP_MINLPProb::writeMIQCPModelInLPFile(const char* fileName)
{
    const int maxbytesinline = 255;
    
    char s[100];
    int code, size;
    //unsigned int rne;
    ofstream fout;
    MIP_StreamMBLine out(&fout, maxbytesinline);
    
    
    fout.open( fileName, ios::out );
    
    if( !fout.is_open() )
    {
        code = MIP_UNDEFINED_ERROR;
        goto termination;
    }
    
    //fout << setprecision(15);
    
    out.print( "\\LP model file generated by minlp problem library, developed by " );
    out.print( MIP_AUTHOR ", " MIP_AUTHOR_FILIATION );
    out.endline();
    out.endline();
    
    if( hasNlObj )
    {
        out.print("\\objective function still has a nonlinear part, but, unfortunatelly, ");
        out.endline();
        out.print("\\this part couldnot be represented in this file.");
        out.endline();
        out.endline();
    }
    
    
    out.print("Minimize ");
    out.endline();
    
    if( d )
        out.print( " + OBJFIX" );
    
    
    //we have to set even 0.0 coeficients in objective function
    for(int i = 0; i < n; i++)
    {
        size = sprintf( s, " %+0.10f x%d", objFactor*c[i], i);
        out.print(s, size);
    }
    
    
    if( Q.getNumberOfElements() > 0 )
    {
        /*out.print(" + [");
        
        for(unsigned int i = 0; i < Q.getNumberOfRows(); i++)
        {
            MIP_SparseRow &row = Q[i];
            unsigned int rne = row.getNumberOfElements();
            
            for(unsigned int j = 0; j < rne; j++)
            {
                int col = row[j].getColumn();
                double v = row[j].getValue();
                
                if(col != i)
                    v *= 2;
                    
                    
                sprintf( s, " %+0.10f x%d", v, i );
                
                if( col == i )
                    size = sprintf(s, "%s ^ 2", s);
                else
                    size = sprintf(s, "%s * x%d", s, col );
                
                out.print(s, size);
            }
        }
        
        out.print(" ]/2"); */
        
        MIP_printMatrixToLPFormat(Q, out, objFactor);
    }
    
    out.endline();
    out.endline();
    
    out.print("Subject to");
    out.endline();
    
    
    for(int i = 0; i < m; i++)
    {
        size = sprintf( s, " cons%04d:", i);
        out.print(s, size);
        
        //MIP_SparseRow &arow = A[i];
        //rne = arow.getNumberOfElements();
        
        unsigned int rne;
        int *rcols;
        double *rvalues;
        
        A.getRowPointers(i, rne, rcols, rvalues);
        
        for(unsigned int j = 0; j < rne; j++)
        {
            size = sprintf(s, " %+0.10f x%d", rvalues[j], rcols[j]);
            
            out.print(s, size);
        }
        
        
        if( QC[i].getNumberOfElements() > 0 )
        {
            MIP_printMatrixToLPFormat(QC[i], out, 1.0);
        }
        
        
        if( lc[i] > -MIP_INFINITY && uc[i] < MIP_INFINITY )
        {
            if( lc[i] == uc[i] )
                size = sprintf(s, " = %0.10f", uc[i]);
            else
                size = sprintf(s, " - W%d = 0", i);
            
        }
        else
        {
            if( lc[i] > -MIP_INFINITY )
                size = sprintf( s, " >= %0.10f", lc[i] );
            else
                size = sprintf( s, " <= %0.10f", uc[i] );
        }
        
        out.print(s, size);
        out.endline();
    }
    
    out.endline();
    out.print("Bounds");
    out.endline();
    
    if(d)
    {
        size = sprintf(s, " %0.10f <= OBJFIX <= %0.10f", d, d);
        out.print(s, size);
        out.endline();
        out.endline();
    }
    
    
    for(int i = 0; i < n; i++)
    {
        if( lx[i] <= -MIP_INFINITY )
            sprintf(s, " -infinity <= x%d <=", i);
        else
            sprintf(s, " %0.10f <= x%d <=", lx[i], i);
        
        if( ux[i] >= MIP_INFINITY )
            size = sprintf(s, "%s +infinity", s);
        else
            size = sprintf(s, "%s %0.10f", s, ux[i]);
        
        out.print(s, size);
        out.endline();
    }
    
    out.endline();
    
    for(int i = 0; i < m; i++)
    {
        if( lc[i] > -MIP_INFINITY && uc[i] < MIP_INFINITY && lc[i] != uc[i] )
        {
            size = sprintf(s, " %0.10f <= W%d <= %0.10f", lc[i], i, uc[i] );
            
            out.print(s, size);
            out.endline();
        }
    }
    
    
    if( nI > 0 )
    {
        out.print("Generals");
        out.endline();
        
        
        for( int i = 0; i < n; i++ )
        {
            if( MIP_isIntegerType(xtype[i]) )
            {
                size = sprintf(s, " x%d", i);
                out.print(s, size);
            }
        }
        
        out.endline();
    }
    
    out.print("End");
    out.endline();
    
    
    
    code = 0;
    
termination:
    
    if( fout.is_open() )
        fout.close();
    
    return code;
}





MIP_IntegratedHessian::MIP_IntegratedHessian(MIP_MINLPProb* prob)
{
    initialize(prob);
}



MIP_IntegratedHessian::~MIP_IntegratedHessian()
{
    desallocate();
}



int MIP_IntegratedHessian::allocate( const unsigned int n )
{
    MIP_malloc(hessIndex, ( (n*(n+1))/2 )); //hessIndex = (unsigned int *) malloc( ( (n*(n+1))/2 ) * sizeof(unsigned int) );
    
    if( !hessIndex )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


int MIP_IntegratedHessian::allocateTripleSparseArrays(const unsigned int size)
{
    MIP_malloc(rows, 2*size); //rows = (int *) malloc( 2*size *sizeof(int) );
    MIP_malloc(values, size); //values = (double *) malloc( size*sizeof(double) );
    
    if( !rows || !values )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    cols = &rows[size];
    
    return 0;
}



static void inline  MIP_markSPMStructureInLowerTriangle( const MIP_SparseMatrix &M, unsigned int *a, const unsigned int value = 1)
{
    const unsigned int nzrows = M.getNumberOfNonZeroRows();
    const unsigned int* nzrowinds = M.nzRowIndex;
    
    /*for( unsigned int i = 0; i < m; i++ )
    {
        MIP_SparseRow &row = M[i];
        unsigned int nel = row.getNumberOfElements();
        
        unsigned int *rowa = &a[ (i*(i+1))/2 ];
        
        for( unsigned int j = 0; j < nel; j++ )
            rowa[ row[j].getColumn() ] = value;
    }*/
    
    
    for( unsigned int k = 0; k < nzrows; k++ )
    {
        const unsigned int i = nzrowinds[k];
        
        unsigned int *rowa = &a[ (i*(i+1))/2 ];
        
        unsigned int rnz = M.getNumberOfElementsAtRow(i);
        int *cols = M.getRowColsPointer(i);
        
        for(unsigned int j = 0; j < rnz; j++)
            rowa[ cols[j] ] = value;
    }
    
}



int MIP_IntegratedHessian::buildStructures(const bool setLagHessianMatrixCopy, const bool setTripleSparArrays )
{
    const unsigned int n = prob->n;
    const unsigned int m = prob->m;
    
    int r;
    MIP_SparseMatrix *QC = prob->QC;
    
    
    r = allocate( n );
    if(r != 0)
        return r;
    
    
    if( setLagHessianMatrixCopy )
    {
        r = lagH.copyStructureFrom( prob->lagH );
        if(r != 0)
        {
            #if MIP_DEBUG_MODE
                cerr << MIP_PREPRINT << "Error " << r << MIP_GETFILELINE << endl;
            #endif
            
            return MIP_MEMORY_ERROR;
        }
    }
    
    
    MIP_setAllArray( (n*(n+1u))/2u , hessIndex, UINT_MAX );
    
    MIP_markSPMStructureInLowerTriangle(lagH, hessIndex);
    
    MIP_markSPMStructureInLowerTriangle(prob->Q, hessIndex);
    
    
    for( unsigned int i = 0; i < m; i++ )
        MIP_markSPMStructureInLowerTriangle( QC[i], hessIndex );
    
    
    nzs = 0;
    for(unsigned int i = 0; i < n; i++)
    {
        unsigned int *hrow = &hessIndex[ (i*(i+1))/2 ];
        
        for(unsigned int j = 0; j <= i; j++)
        {
            if( hrow[j] != UINT_MAX )
            {
                hrow[j] = nzs;
                nzs++;
            }
        }
    }
    
    
    if( setTripleSparArrays )
    {
        r = allocateTripleSparseArrays(nzs);
        if(r != 0)
        {
            #if MIP_DEBUG_MODE
                cerr << MIP_PREPRINT << "Error " << r << MIP_GETFILELINE << endl;
            #endif
            
            return MIP_MEMORY_ERROR;
        }
        
        setTripleSparseArrays();
    }
    
    
    return 0;
}


void MIP_IntegratedHessian::desallocate()
{
    nzs = 0;
    MIP_secFree(hessIndex);
    
    desallocateTripleSparseArrays();
}


void MIP_IntegratedHessian::desallocateTripleSparseArrays()
{
    MIP_secFree(rows);
    //MIP_secFree(cols);
    MIP_secFree(values);
}




int MIP_IntegratedHessian::evalCompleteHessian( const int thnumber, const bool newx, const double *x, const double objFactor, const double *lambda, MIP_SparseMatrix *lagH, int *rows, int *cols, double *values )
{
    const int &n = prob->n;
    const int &m = prob->m;
    
    int r;
    //int *myrows = rows ? rows : this->rows;
    //int *mycols = cols ? cols : this->cols;
    double *myvalues = values ? values : this->values;
    
    MIP_SparseMatrix *myLagH = lagH ? lagH : &(this->lagH);
    
    MIP_SparseMatrix &Q = prob->Q;
    MIP_SparseMatrix *QC = prob->QC;
    
    
    MIP_setAllArray( nzs, myvalues, 0.0 );
    
    
    if( prob->hasNlObj || prob->hasNlConstrs )
    {
        r = prob->nlpHessianEval( thnumber, newx, x, objFactor, lambda, *myLagH );
        
        if( r != 0 )
        {
            #if MIP_DEBUG_MODE
                MIP_PRINTERRORNUMBER(r);
            #endif
            
            return r;
        }
        
        
        /*for( int i = 0; i < n; i++ )
        {
            MIP_SparseRow &row = (*myLagH)[i];
            
            int rne = row.getNumberOfElements();
            
            const unsigned int *hindrow = &hessIndex[ (i*(i+1))/2 ];
            
            //we perform in the reverse order because in this way, we can avoi some problems if user set the same sparse structures two times. In this way, the last value will be the first...
            
            for( int j = rne-1; j >= 0; j-- )
            {
                const int ind = hindrow[ row[j].getColumn() ];
                
                myvalues[ ind ] = row[j].getValue();
            }
        } */
        
        { // here we are worried about performance. So, we run in a nonelegant way (no iterator)
            const unsigned int nnzrows = myLagH->getNumberOfNonZeroRows();
            const unsigned int *nzrowinds = myLagH->nzRowIndex;
            
            for(unsigned int k = 0; k < nnzrows; k++)
            {
                const unsigned int i = nzrowinds[k];
                
                const unsigned int *hindrow = &hessIndex[ (i*(i+1))/2 ];
                
                unsigned int rnz;
                int *rcols;
                double *rvalues;
                
                myLagH->getRowPointers(i, rnz, rcols, rvalues);
                
                for(unsigned int j = rnz-1; ;j--)
                {
                    myvalues[  hindrow[rcols[j]]  ] = rvalues[j];
                    
                    if(j==0)
                        break;
                }
            }
        }
        
        
    }
    
    
    if( Q.getNumberOfElements() > 0 && objFactor != 0 )
    {
        const double f = prob->objFactor * objFactor;
        
        //cout << "Vou setar Q na hessiana" << endl;
        
        setSPMonHessianValues(n, Q, f, myvalues);
    }
    
    //cout << "Vou setar Q das restricoes" << endl;
    
    for( int k = 0; k < m; k++ )
    {
        //cerr << "k: " << k << " m: " << m << endl;
        //cerr << "lambda[" << k << "]: " << lambda[k] << endl;
        if( QC[k].getNumberOfElements() > 0 && lambda[k] != 0.0 )
            setSPMonHessianValues(n, QC[k], lambda[k], myvalues);
    }
    
    
    return 0;
}



void MIP_IntegratedHessian::initialize(MIP_MINLPProb *prob)
{
    nzs = 0;
    
    rows = NULL;
    cols = NULL;
    values = NULL;
    
    hessIndex = NULL;
    
    this->prob = prob;
}



void MIP_IntegratedHessian::setTripleSparseArrays( int *rows, int *cols )
{
    const unsigned int n = prob->n;
    unsigned int k = 0;
    
    int *myrows = rows ? rows : this->rows ;
    int *mycols = cols ? cols : this->cols ;
    
    
    for(unsigned int i = 0; i < n; i++)
    {
        unsigned int *hrow = &hessIndex[(i*(i+1))/2];
        
        for(unsigned int j = 0; j <= i; j++)
        {
            if( hrow[j] != UINT_MAX )
            {
                #if MIP_DEBUG_MODE
                    assert( k == hrow[j] );
                #endif
                
                myrows[k] = i;
                mycols[k] = j;
                k++;
            }
        }
    }
    
    #if MIP_DEBUG_MODE
        assert( k == nzs );
    #endif
    
}



void MIP_IntegratedHessian::setSPMonHessianValues(const int n, const MIP_SparseMatrix& M, const double factor, double* values)
{
    const unsigned int nzrows = M.getNumberOfNonZeroRows();
    const unsigned int *nzrowinds = M.nzRowIndex;
    
    for(unsigned int k = 0; k < nzrows; k++)
    {
        const unsigned int i = nzrowinds[k];
        const unsigned int *hindrow = &hessIndex[ (i*(i+1))/2 ];
        
        unsigned int rnz;
        int *rcols;
        double *rvalues;
        
        M.getRowPointers(i, rnz, rcols, rvalues);
        
        for(unsigned int j = rnz-1; ; j++)
        {
            values[ hindrow[rcols[j]] ] += factor*rvalues[j] ;
            if(j == 0)
                break;
        }
    }
    
    /*for(int i = 0; i < n; i++)
    {
        MIP_SparseRow &row = M[i];
        
        const int rne = row.getNumberOfElements();
        
        const unsigned int *hindrow = &hessIndex[ (i*(i+1))/2 ];
        
        //we perform in the reverse order because in this way, we can avoid some problems if user set the same sparse structures two times. In this way, the last value will be the first...
        
        for(int j = rne-1; j >= 0; j--)
        {
            const unsigned int ind = hindrow[ row[j].getColumn() ];
            
            //cerr << "j: " << j << " ind: " << ind;
            //cerr << " values: " << values[ind] << endl;
            
            values[ ind ] += factor*row[j].getValue();
        }
    } */
}








int minlpproblem::MIP_checkRowStructure(const unsigned int nrows, const unsigned int ncols, const bool symmetric, unsigned int row, unsigned int nzs, int* cols, double* values)
{
    int code;
    
    const int r = newspm::SPM_checkRowStructure( nrows, ncols, symmetric, row, nzs, cols, values );
    
    switch(r)
    {
        case 0:
            code = 0;
            break;
        
        case newspm::SPM_INDEX_FAULT:
            code = MIP_INDEX_FAULT;
            break;
        
        case newspm::SPM_BAD_VALUE:
            code = MIP_BAD_VALUE;
            break;
        
        case newspm::SPM_REPETEAD_ELEMENT:
            code = MIP_REPETEAD_INDEXES;
            break;
            
        case newspm::SPM_UPPER_TRIANGLE:
            code = MIP_UPPER_TRIANGLE_INDEX;
            break;
        
        default:
            code = MIP_UNDEFINED_ERROR;
    }
    
    return code;
}






int minlpproblem::MIP_checkTripleSparseStructure(const int nrows, const int ncols, const bool symmetric, const unsigned int nzs, const int* rows, const int* cols, const double* values)
{
    int code;
    
    const int r = newspm::SPM_checkTripleSparseStructure( nrows, ncols, symmetric, nzs, rows, cols, values );
    
    switch(r)
    {
        case 0:
            code = 0;
            break;
        
        case newspm::SPM_INDEX_FAULT:
            code = MIP_INDEX_FAULT;
            break;
        
        case newspm::SPM_BAD_VALUE:
            code = MIP_BAD_VALUE;
            break;
        
        case newspm::SPM_REPETEAD_ELEMENT:
            code = MIP_REPETEAD_INDEXES;
            break;
            
        case newspm::SPM_UPPER_TRIANGLE:
            code = MIP_UPPER_TRIANGLE_INDEX;
            break;
        
        default:
            code = MIP_UNDEFINED_ERROR;
    }
    
    return code;
}




void minlpproblem::MIP_constrCompleteGrad(const MIP_MINLPProb &prob, const MIP_SparseMatrix &Jac, const int constr, const double *x, double *grad, const bool initializeWithZeros)
{
    if(initializeWithZeros)
        MIP_setAllArray(prob.n, grad, 0.0);
    
    
    if( prob.nlConstr[constr] )
    {
        //Jac.getFullRow(constr, grad, false, false, 1.0);
        
        Jac.copyRowTo(constr, grad, 1.0);
    }
    /*else
    {
        //const int n = prob.n;
        //for( int i = 0; i < n; i++)
            //grad[i] = 0.0;
        MIP_setAllArray(prob.n, grad, 0.0);
    }*/
    
    if( prob.QC[constr].getNumberOfElements() > 0 )
        prob.QC[constr].quadraticGradientEvaluation(x, grad, true);
    
    //prob.A.getFullRowAccumulation(constr, grad);
    prob.A.getFullRow(constr, grad, true, false, 1.0);
}




void minlpproblem::MIP_completeLagHessianRow(const MIP_MINLPProb &prob, const int mquad, const int *quadIndex, const MIP_SparseMatrix &lagH, const double objFactor, const double *lambda, const int rowIndex, double *rowValues, const bool initializeWithZeros)
{
    const double realObjFactor = prob.objFactor*objFactor;
    
    if(initializeWithZeros)
        MIP_setAllArray(rowIndex + 1, rowValues, 0.0);
    
    if( prob.hasNlObj || prob.hasNlConstrs )
    {
        lagH.copyRowTo(rowIndex, rowValues, 1.0);
    }
    
    
    if(realObjFactor != 0.0)
    {
        int r = prob.Q.accumulateRowInArray(rowIndex, rowValues, realObjFactor);
        #if MIP_DEBUG_MODE
            assert(r == 0);
        #endif
    }
    
    for(int i = 0; i < mquad; i++)
    {
        const int ind = quadIndex[i];
        
        if( lambda[ind] != 0.0 )
        {
            int r = prob.QC[ind].accumulateRowInArray(rowIndex, rowValues, lambda[ind]);
            #if MIP_DEBUG_MODE
                assert(r == 0);
            #endif
        }
    }
    
}





MIP_BinSumConstrsIndsByClass::MIP_BinSumConstrsIndsByClass()
{
    initialize();
}


MIP_BinSumConstrsIndsByClass::~MIP_BinSumConstrsIndsByClass()
{
    deallocate();
}


void MIP_BinSumConstrsIndsByClass::deallocate()
{
    //pointer will be set to null in initialize method
    for(unsigned int i = 0; i < nBinSumConstrClasses; i++)
    {
        if(classes[i])
            free(classes[i]);
    }
    
    if(knapsackInds)
    {
        for(decltype(nI) i = 0; i < nI; i++)
        {
            if(knapsackInds[i])
                free(knapsackInds[i]);
        }
        free(knapsackInds);
    }
    
    if(nKnapsackInds)
        free(nKnapsackInds);
    
    initialize();
}


void MIP_BinSumConstrsIndsByClass::initialize()
{
    for(unsigned int i = 0; i < nBinSumConstrClasses; i++)
        nClasses[i] = 0;
    
    MIP_setAllArray<unsigned int>( nBinSumConstrClasses, nClasses, 0);
    
    MIP_setAllArray<unsigned int *>( nBinSumConstrClasses, classes, NULL);
    
    knapsackInds = NULL;
    nKnapsackInds = NULL;
    
    nI = 0;
}


int MIP_BinSumConstrsIndsByClass:: addIndexToClassArray( const unsigned int classNumber, const unsigned int index, const unsigned maxNumberOfIndices )
{
    if( !classes[classNumber] )
    {
        MIP_malloc(classes[classNumber], maxNumberOfIndices );
        MIP_IFMEMERRORRETURN( !classes[classNumber] );
    }
    
    classes[classNumber][ nClasses[classNumber] ] = index;
    nClasses[classNumber]++;
    
    return 0;
}


int  MIP_BinSumConstrsIndsByClass:: addIndexToKnapsackIndsArray( const unsigned int intVarIndex, const unsigned int knapsackConstrIndex, const unsigned maxNumberOfIndices)
{
    if( !knapsackInds[intVarIndex] )
    {
        MIP_malloc(knapsackInds[intVarIndex], maxNumberOfIndices);
        MIP_IFMEMERRORRETURN( !knapsackInds[intVarIndex] );
    }
    
    knapsackInds[intVarIndex][ nKnapsackInds[intVarIndex] ] = knapsackConstrIndex;
    nKnapsackInds[intVarIndex]++;
    
    return 0;
}


//here, we reallocate class arrays to have the exact necessary size to store index
void MIP_BinSumConstrsIndsByClass:: adjustClassesArrays()
{  
    for(unsigned int i = 0; i < nBinSumConstrClasses; i++)
        MIP_realloc(classes[i], nClasses[i]);
}


//here, we reallocate knapsack arrays to have the exact necessary size to store index
void MIP_BinSumConstrsIndsByClass:: adjustKnapsacksArrays()
{
    if( knapsackInds )
    {
        for(unsigned int i = 0; i < nI; i++)
            MIP_realloc( knapsackInds[i], nKnapsackInds[i] );
    }
}


int MIP_BinSumConstrsIndsByClass:: calculateIndices( const MIP_MINLPProb &prob, const double *lx, const double *ux, const double *lc, const double *uc,  int *reverseIntVars, bool calculateKnapsakIndices, bool sortClass0And3IndicesByb, bool substituteFixVars )
{
    const int m = prob.m;
    const unsigned int knapsackClassNumber = 7;
    
    const bool *nlConstr = prob.nlConstr;
    const int *xtype = prob.xtype;

    const MIP_SparseMatrix *QC = prob.QC;
    const MIP_SparseMatrix &A = prob.A;
    
    const double ONE = 1.0;
    
    int r, retCode;
    
    //int *myIntVars = NULL;
    int *myReverseIntVars = NULL;
    
    
    if( !lx )
        lx = prob.lx;
    if( !ux )
        ux = prob.ux;
    if( !lc )
        lc = prob.lc;
    if( !uc )
        uc = prob.uc;
    
    
    deallocate();
    
    


    //detecting constraints: sum_{i \in I} x_i = 1
    //for(int i = 0; i < m; i++)
    for( minlpproblem::MIP_SparseMatrixRowIndexIterator rowait = A.beginRowIndex() ;  *rowait < (unsigned int) m ; ++rowait )
    {
        const int i = *rowait;
        
        if( lc[i] <= -MIP_INFINITY && uc[i] >= MIP_INFINITY )  //free constraint. This constraint was discarded by user or by the preprocessor.
            continue;
        
        if(nlConstr[i] || QC[i].getNumberOfElements() > 0)
            continue;
        
        
        
        const int* acols = A[i];
        const double* avalues = A(i);
        const unsigned int nel = A.getNumberOfElementsAtRow(i); 
        
        #if MRQ_DEBUG_MODE
            assert( nel > 0 ); //there must be at least some element in this row of A (we are running rowIndexIterator)
        #endif
        
        unsigned int nCoefOneInt = 0; //number of integer variables having one coeficient (including binaries)
        unsigned int nCoefMinusOneInt = 0; //number of integer variables having minus one coeficient (including binaries)
        unsigned int nCoefPosInt = 0; //number of positive coeficients in integer variables (including binaries)
        unsigned int nCoefNegInt = 0; //number of negatives coeficients in integer variables (including binaries)
        unsigned int nCoefZero = 0; //number of zero coefficients
        unsigned int nBin = 0, nGenInt = 0, nCont = 0; //number of binary, general integer and continuos vars...
        
        double subVarValue = 0.0; //value from variable substituiton
        
        for(unsigned int j = 0; j < nel; j++)
        {
            const auto var = acols[j]; //row[j].getColumn();
            const auto coef = avalues[j];
            
            
            if( coef == 0.0 )
            {
                nCoefZero++;
                continue; //this variable is not in the constraint by this coefficient, actually
            }
            
            
            if( substituteFixVars )
            {
                if( lx[var] == ux[var] )
                {
                    subVarValue += coef*var;
                    continue; //here, we do not care if variable is continuous or integer, neither about its coefficients because we will treat it like a constant 
                }
            }
            
            
            if( MIP_isIntegerType(xtype[var] ) )
            {
                if( lx[var] > -1.0 && ux[var] < 2.0 )
                    nBin++;
                else
                    nGenInt++;
                
                if( coef > 0.0)
                {
                    nCoefPosInt++;
                    if( coef == ONE )
                        nCoefOneInt++;
                }
                else
                {
                    nCoefNegInt++;
                    if( coef == -ONE )
                        nCoefMinusOneInt++;
                }
                
            }
            else
            {
                nCont++;
            }
            
        }
        
        
        if( nCont == 0 )
        {
        
            double rhs = ( uc[i] < MIP_INFINITY ? uc[i] : lc[i] )- subVarValue  ; //we are sure that is not a free constraint
            
            if( !MIP_isIntegerDouble(rhs) )
                continue;
            
            
            if( rhs > 0.0) //we put 0.0 because class5 we could have an rhs = 0.5, for example, and that would be valid.  
            {
                
                if( nGenInt == 0 && nCoefPosInt == nCoefOneInt && nCoefOneInt > 0 && nCoefNegInt == 0 )
                {
                    
                    //equality constraint, rhs >= 1, only binary variables, all coeficients are 1.0 
                    if( lc[i] == uc[i] && rhs >= 1 )
                    {
                        const unsigned int classNumber = 0;
                        r = addIndexToClassArray(classNumber, i, m);
                        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                    }
                    //lower than equal constraint, rhs >= 1, only binary variables, all coeficients are 1.0 
                    else if( lc[i] <= -MIP_INFINITY && rhs >= 1 )
                    {
                        #if MIP_DEBUG_MODE
                            assert( uc[i] < MIP_INFINITY );
                        #endif
                        
                        const unsigned int classNumber = 3;
                        r = addIndexToClassArray(classNumber, i, m);
                        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                    }
                    //greater than equal constraint, rhs > 0
                    else if( uc[i] >= MIP_INFINITY )
                    {
                        #if MIP_DEBUG_MODE
                            assert( lc[i] > -MIP_INFINITY );
                        #endif
                        
                        const unsigned int classNumber = 4;
                        r = addIndexToClassArray(classNumber, i, m);
                        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                    }
                    
                }
                
                
                if( nGenInt == 0 && nCoefPosInt > nCoefOneInt  && nCoefNegInt == 0  )
                {
                    if( lc[i] <= -MIP_INFINITY  )
                    {
                        //less than equal constraint, only binary variables, at least one positive coefficent different than one, no negative indices
                        #if MIP_DEBUG_MODE
                            assert( uc[i] < MIP_INFINITY );
                        #endif
                        
                        //0-1 knapsack
                        
                        const unsigned int classNumber = 7;
                        r = addIndexToClassArray(classNumber, i, m);
                        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                    }
                    else if( uc[i] >= MIP_INFINITY )
                    {
                        //greater than equal constraint, only binary variables, at least one positive coefficent different than one, no negative indices
                        
                        #if MIP_DEBUG_MODE
                            assert( lc[i] > -MIP_INFINITY );
                        #endif
                        
                        const unsigned int classNumber = 8;
                        r = addIndexToClassArray(classNumber, i, m);
                        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                    }
                    
                }
                
            }
            
            //now we allow classes 1, 2, 5 and 6 have nonnull rhs (now rhs can be zero, positive or negative)
            if( nGenInt == 0 && nCoefPosInt == nCoefOneInt && nCoefOneInt > 0 && nCoefMinusOneInt >= 1 && nCoefNegInt == nCoefMinusOneInt )
            {
                
                if( nCoefMinusOneInt == 1 )
                {
                    //equality constraint,  only binary variables, all positive coefs are 1, all negative coefas are -1, only one negative coef
                    if( lc[i] == uc[i] )
                    {
                        const unsigned int classNumber = 1;
                        r = addIndexToClassArray(classNumber, i, m);
                        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                    }
                    //lower than equal constraint, only binary variables, all positive coefs are 1, all negative coefas are -1, only one negative coef
                    else if( lc[i] <= -MIP_INFINITY )
                    {
                        #if MIP_DEBUG_MODE
                            assert( uc[i] < MIP_INFINITY );
                        #endif
                        
                        const unsigned int classNumber = 2;
                        r = addIndexToClassArray( classNumber, i, m);
                        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                    }
                    
                }
                else if( nCoefMinusOneInt > 1 )
                {
                    //equality constraint,  only binary variables, all positive coefs are 1, all negative coefas are -1, more than one negative coef
                    if( lc[i] == uc[i] )
                    {
                        
                        const unsigned int classNumber = 5;
                        r = addIndexToClassArray( classNumber, i, m);
                        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                    }
                    //less than equal constraint, only binary variables, all positive coes are 1, all negative coefs are -1, more than one negative coef
                    else if( lc[i] <= -MIP_INFINITY )
                    {
                        const unsigned int classNumber = 6;
                        r = addIndexToClassArray( classNumber, i, m);
                        MIP_IFERRORGOTOLABEL(r, retCode, r, termination);
                    }
                    //TODO: threat here a possible new class to reverse the class 6: x_{p_1} + x_{p_2} + x_{p_3} + ... + x_{p_k} - x_{n_1} - x_{n_2} - x_{n_3} - ... - x_{n_q} >= b
                    
                }
                
            } //end of if( nGenInt == 0 && nCoefPosInt == nCoefOneInt && nCoefOneInt > 0 && nCoefNegInt == nCoefMinusOneInt )
            
        
            
            
        }
        
    }
    
    adjustClassesArrays();
    
    
    if( sortClass0And3IndicesByb )
    {
        const double *b = prob.uc;
        
        if( nClasses[0] > 0 )
        {
            const unsigned int classNumber = 0;
            
            //selection sort
            MIP_sortIndicesByWeight( nClasses[classNumber], classes[classNumber], b );
        }
        
        if( nClasses[3] > 0 )
        {
            const unsigned int classNumber = 3;
            
            //selection sort
            MIP_sortIndicesByWeight( nClasses[classNumber], classes[classNumber], b );
        }
    }
    
    this->nI = prob.getNumberOfIntegerVars();
    
    if( calculateKnapsakIndices && nClasses[knapsackClassNumber] > 0 )
    {
        //int *pIntVars;
        int *pReverseIntVars;
        int mynI  = prob.getNumberOfIntegerVars();
        
        
        /*if( intVars )
        {
            mynI = nI;
            pIntVars = intVars;
        }
        else
        {
            mynI = prob.getNumberOfIntegerVars();
            
            MIP_malloc(myIntVars, mynI);
            MIP_IFMEMERRORGOTOLABEL(!myIntVars, retCode,  termination);
            
            prob.getIntegerIndices(myIntVars);
            
            pIntVars = myIntVars;
        } */
        
        if(reverseIntVars)
        {
            pReverseIntVars = reverseIntVars;
        }
        else
        {
            const int n = prob.n;
            
            MIP_malloc(myReverseIntVars, n);
            MIP_IFMEMERRORGOTOLABEL(!myReverseIntVars, retCode, termination);
            
            prob.getReverseIntegerIndices(myReverseIntVars);
            pReverseIntVars = myReverseIntVars;
        }
        
        
        MIP_calloc(knapsackInds, mynI);
        MIP_calloc(nKnapsackInds, mynI);
        MIP_IFMEMERRORGOTOLABEL(!knapsackInds || !nKnapsackInds, retCode, termination);
        
        
        auto nknapsackIndices = nClasses[knapsackClassNumber];
        auto knapsackIndices = classes[knapsackClassNumber];
        
        for(unsigned int k = 0; k < nknapsackIndices; k++)
        {
            auto i = knapsackIndices[k];
            
            const int* acols = A[i];
            const double* avalues = A(i);
            const unsigned int nel = A.getNumberOfElementsAtRow(i);
            
            #if MRQ_DEBUG_MODE
                assert( nel > 0 ); 
            #endif
            
            for(unsigned int j = 0; j < nel; j++)
            {
                const auto var = acols[j];
                const auto coef = avalues[j];
                
                if( !MIP_isIntegerType(xtype[var]) )
                    continue;
                
                if( coef == 0.0 )
                    continue;
                
                r = addIndexToKnapsackIndsArray( pReverseIntVars[var], i, nknapsackIndices );
                MIP_IFERRORGOTOLABEL( r, retCode, r, termination );
            }
        }
        
        adjustKnapsacksArrays();
        
    }
    
    
    retCode = 0;


termination:
    
    //if(myIntVars)   free(myIntVars);
    if(myReverseIntVars)    free(myReverseIntVars);
    
    return retCode;
}





