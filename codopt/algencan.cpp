/*
* Implementations to OPT_Algencan NLP solver.
* It is a chalenge because Algencan has no manual, information or nothing more to help us unless 3 very basic examples.
* 
* http://www.ime.usp.br/~egbirgin/tango/codes.php
* 
* To make algencan does not print the preamble when it run, it is necessary put a file called ".silent" in the current directory.
* 
* By Wendel Melo
* 16 October 2016
* 
*/ 


#include <cmath>
#include <cassert>
#include <climits>


#include <iostream>
#include <new>
#include <unordered_map>
#include <thread>


#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"


using namespace minlpproblem;
using namespace optsolvers;


typedef bool ALG_bool;


#if OPT_HAVE_ALGENCAN

extern "C"
{
    void c_algencan(
        void (*myevalf)(int n, double *x, double *f, int *flag),
        void (*myevalg)(int n, double *x, double *g, int *flag),
        void (*myevalh)(int n, double *x, int *hrow, int *hcol, double *hval, int *hnnz, int lim, ALG_bool *lmem, int *flag),
        void (*myevalc)(int n, double *x, int ind, double *c, int *flag),
        void (*myevaljac)(int n, double *x, int ind, int *jcvar, double *jcval, int *jcnnz, int lim, ALG_bool *lmem, int *flag),
        void (*myevalhc)(int n, double *x, int ind, int *hcrow, int *hccol, double *hcval, int *hcnnz, int lim, ALG_bool *lmem, int *flag),
        void (*myevalfc)(int n, double *x, double *f, int m, double *c, int *flag),
                    
        void (*myevalgjac)(int n, double *x, double *g, int m, int *jcfun, int *jcvar, double *jcval, int *jcnnz, int lim, ALG_bool *lmem, int *flag),
        void (*myevalgjacp)(int n, double *x, double *g, int m, double *p, double *q, char work, ALG_bool *gotj, int *flag),
        void (*myevalhl)(int n, double *x, int m, double *lambda, double scalef, double *scalec, int *hlrow, int *hlcol, double *hlval, int *hlnnz, int lim, ALG_bool *lmem, int *flag),
        void (*myevalhlp)(int n, double *x, int m, double *lambda, double scalef, double *scalec, double *p, double *hp, ALG_bool *goth, int *flag),
        int jcnnzmax,
        int hnnzmax,
        
        double *epsfeas,
        double *epsopt,
        double *efstin,
        double *eostin,
        double *efacc,
        double *eoacc,
        char *outputfnm,
        char *specfnm,
        
        int nvparam,
        char **vparam,
        int n,
        double *x,
        double *l,
        double *u,
        int m,
        double *lambda,
        ALG_bool *equatn,
        ALG_bool *linear,
        ALG_bool *coded,
        ALG_bool  checkder,
        double *f,
        double *cnorm,
        double *snorm,
        double *nlpsupn,
        int *inform
    );
}

#endif




namespace optsolvers{



//NOTE THIS WILL NOT WORK WITH MULTIPLE THREADS BECAUSE WE ARE OBLIGATED TO WORK WITH A GLOBAL VARIABLE TO HAVE AUXILIARY DATA ON EVALUATIONS

static std::unordered_map<std::thread::id, OPT_Algencan::OPT_MyData*> threadAlgData;




//Evaluates objective and constraints (myevalfc)
static void OPT_algencanEvalObjAndConstrs(int n, double *x, double *f, int m, double *c, int *flag)
{
    //std::cerr << "Entrei em OPT_algencanEvalObjAndConstrs\n";
    
    OPT_Algencan::OPT_MyData *algData = threadAlgData[std::this_thread::get_id()];
    const unsigned int thnumber = algData->thnumber;
    const bool *constrEval = algData->constrEval;
    double *g = algData->auxConstr;
    OPT_Algencan *nlp = algData->nlp;
    OPT_MINLPProb *prob = &nlp->prob;
    
    
    int r;
    
    //for(int i = 0; i < n; i++)
        //std::cout << "x["<<i<<"]: " << x[i] << "\n";
    
    r = prob->objEval(thnumber, true, x, *f, nlp->in_nl_obj_factor);
    if(r != 0)
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            OPT_PRINTERRORNUMBER(r);
        #endif
        *flag = r;
        return;
    }
    
    r = prob->constraintsEval(thnumber, !prob->hasNlObj, constrEval, x, g);
    if(r != 0)
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            OPT_PRINTERRORNUMBER(r);
        #endif
        *flag = r;
        return;
    }
        
    
    //unfortunatelly, all constraints in algencan has rhs as 0
    {
        const int om = prob->m;
        const double *lc = prob->lc;
        const double *uc = prob->uc;
        
        
        int dbconstr = prob->n; //double bounded constraint counter to auxiliary variable. First auxiliary variable has index n.
        
        int k = 0;
        for(int i = 0; i < om; i++)
        {
            if(!constrEval[i])
                continue;
            
            c[k] = g[i];
            
            if(uc[i] < MIP_INFINITY)
            {
                if(lc[i] > -MIP_INFINITY)
                { 
                    if(lc[i] == uc[i])
                    { //equality constraint
                        c[k] -= lc[i];
                    }
                    else
                    {
                        //double bounded constraint (we have set as equality constraint using our auxiliary variable)
                        c[k] -= x[dbconstr];
                        dbconstr++;
                    }
                }
                else
                { //only upper bound constraint
                    c[k] -= uc[i];
                }
            }
            else
            {
                if(lc[i] > -MIP_INFINITY)
                { //only lower bound constraint. Unfortunatelly algencan only work with upper bound constraints. So, we have to multiply by -1
                    c[k] = -c[k] + lc[i];
                }
                else
                { //free constraint (we no waste time evaluating)
                    c[k] = 0.0;
                }
            }
            
            k++;
        }
        
        #if OPT_DEBUG_MODE
            assert(k == m);
            assert(n == dbconstr);
        #endif
    }
    
    
    /*for(int i = 0; i < m; i++)
    {
        if( std::isnan(c[i]) || std::isinf(c[i]) || c[i] > 1e10 || c[i] < -1e10 )
        {
            printf("Achei nan or inf c[%d]: %f\n", i, c[i]);
            OPT_getchar();
        }
    } */
    
    //OPT_getchar();
    
    
    *flag = 0;
    //std::cerr << "Sai de OPT_algencanEvalObjAndConstrs\n";
}





//evaluates objective gradient and jacobian (myevalgjac)
static void OPT_algencanEvalObjAndConstrGrads(int n, double *x, double *grad, int m, int *jcfun, int *jcvar, double *jcval, int *jcnnz, int lim, ALG_bool *lmem, int *flag)
{
    //std::cerr << "Entrei em OPT_algencanEvalObjAndConstrGrads\n";
    
    OPT_Algencan::OPT_MyData *algData = threadAlgData[std::this_thread::get_id()];
    const unsigned int thnumber = algData->thnumber;
    const bool *constrEval = algData->constrEval;
    const int *sizeColsNzRowJac = algData->sizeColsNzRowJac;
    int * const *colsNzRowJac = algData->colsNzRowJac;
    double *auxVars = algData->auxVars;
    
    OPT_Algencan *nlp = algData->nlp;
    OPT_MINLPProb *prob = &nlp->prob;
    //const int orign = prob->n;
    
    bool newx = true;
    int r;
    
    
    *flag = 0;
    *lmem = false;
    *jcnnz = 0; //we need initialize it for the case where we have no memeory enough
    
    
    //evaluating jacobian
    {
        const double *lc = prob->lc;
        const double *uc = prob->uc;
        const bool *nlConstr = prob->nlConstr;
        //double *pg = auxVars;
        
        
        const int om = prob->m;
        
        const MIP_SparseMatrix *QC = prob->QC;
        MIP_SparseMatrix &J = prob->J;
        const MIP_SparseMatrix &A = prob->A;
        
        int dbconstrAuxIndex = prob->n; //double bounded constraint counter to auxiliary variable. First auxiliary variable has index n.
        
        
        if( prob->hasNlConstrs )
        {
            int r = prob->nlJacobianEval(thnumber, newx, constrEval, x, J);
            if(r != 0)
            {
                #if OPT_PRINT_CALLBACK_ERROR_MSG
                    OPT_PRINTERRORNUMBER(r);
                #endif
                *flag = r;
                return;
            }
            
            newx = false;
        }
        
        
        
        int jacIndex = 0;
        
        int k = 0;
        for(int i = 0; i < om; i++)
        {
            if(!constrEval[i])
                continue; //free constraint
            
            int nzs;
            
            if( nlConstr[i] && QC[i].getNumberOfElements() == 0 && A.getNumberOfElementsAtRow(i) == 0 )
            {
                nzs = J.getNumberOfElementsAtRow(i);
                
                if(jacIndex + nzs > lim)
                {
                    *lmem = true;
                    return;
                }
                
                J.getRowStructure(i, &jcvar[jacIndex], NULL);
                J.getRowValues(i, &jcval[jacIndex], NULL);
                
            }
            else if( QC[i].getNumberOfElements() == 0 && (!nlConstr[i] || J.getNumberOfElementsAtRow(i) == 0) )
            {
                nzs = A.getNumberOfElementsAtRow(i);
                
                if(jacIndex + nzs > lim)
                {
                    *lmem = true;
                    return;
                }
                
                A.getRowStructure(i, &jcvar[jacIndex], NULL);
                A.getRowValues(i, &jcval[jacIndex], NULL);
                
            }
            else
            {
                const int nminds = sizeColsNzRowJac[i];
                const int *pcols = colsNzRowJac[i];
                
                nzs = 0;
                
                //initializng the nonzero radient positions for this constraint
                for(int j = 0; j < nminds; j++)
                    auxVars[ pcols[j] ] = 0.0;
                
                
                MIP_constrCompleteGrad(*prob, J, i, x, auxVars, false);
                
                
                for(int j = 0; j < nminds; j++)
                {
                    int index = jacIndex + nzs;
                    const int col = pcols[j];
                    
                    if(index >= lim)
                    {
                        *lmem = true;
                        return;
                    }
                    
                    jcvar[index] = col;
                    jcval[index] = auxVars[col];
                    
                    //pg[col] = 0.0;
                    
                    nzs++;
                }
                
                
                /*for(int j = 0; j < n; j++)
                {
                    if(g[j] != 0.0)
                    {
                        int index = jacIndex + nzs;
                        
                        if(index >= lim)
                        {
                            *lmem = true;
                            return;
                        }
                        
                        
                        //jcfun[index] = i;
                        jcvar[index] = j;
                        jcval[index] = g[j];
                        
                        nzs++;
                    }
                }*/
                
            }
            
            
            if(lc[i] > -MIP_INFINITY)
            {
                if(uc[i] >= MIP_INFINITY)
                { //we have an lower bound constraint. We have to multiply by -1
                    OPT_multiplyAllArray(nzs, -1.0, &jcval[jacIndex] );
                }
                else if(lc[i] != uc[i])
                {
                    //double bounded constraint. We have to add an auxiliary variable
                    
                    int index = jacIndex + nzs;
                    
                    if(jacIndex + nzs > lim)
                    {
                        *lmem = true;
                        return;
                    }
                    
                    jcvar[index] = dbconstrAuxIndex;
                    jcval[index] = -1.0;
                    
                    nzs++;
                    dbconstrAuxIndex++;
                }
                
            }
            
            
            //setting the row indices
            OPT_setAllArray(nzs, &jcfun[jacIndex], k);
            
            k++;
            jacIndex += nzs;
        }
        
        *jcnnz = jacIndex;
        
        #if OPT_DEBUG_MODE
            assert(k == m);
            assert(dbconstrAuxIndex == n);
        #endif
    }
    
    
    //do not eval bjective gradient first because we use grad as an auxiliary array
    r = prob->objGradEval(thnumber, newx, x, grad, nlp->in_nl_obj_factor);
    if(r != 0)
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            OPT_PRINTERRORNUMBER(r);
        #endif
        *flag = r;
        return;
    }
    
    //seting zero for auxiliary variables...
    if(n > prob->n)
        OPT_setAllArray(n - prob->n, &grad[prob->n], 0.0);
    
    
    //std::cerr << "Sai de OPT_algencanEvalObjAndConstrGrads\n";
}





//evaluates lagrangian hessian*p (myevalhlp)
static void OPT_algencanEvalHessVector(int n, double *x, int m, double *lambda, double scalef, double *scalec, double *p, double *hp, ALG_bool *goth, int *flag)
{
    //std::cerr << "Entrei em OPT_algencanEvalHessVector\n";
    
    OPT_Algencan::OPT_MyData *algData = threadAlgData[std::this_thread::get_id()];
    const unsigned int thnumber = algData->thnumber;
    //const bool *constrEval = algData->constrEval;
    OPT_Algencan *nlp = algData->nlp;
    OPT_MINLPProb *prob = &nlp->prob;
    
    
    const double realObjF = prob->objFactor * scalef;
    
    //const bool evalLag = (prob->hasNlObj && realObjF != 0.0) || prob->hasNlConstrs;
    
    
    bool someNzLambda = false;
    const bool *nlConstr = prob->nlConstr;
    const double *lc = prob->lc;
    const double *uc = prob->uc;
    double *mylambda = algData->auxConstr;
    const bool *constrEval = algData->constrEval;
    
    MIP_SparseMatrix &lagH = prob->lagH;
    const MIP_SparseMatrix &Q = prob->Q;
    const MIP_SparseMatrix *QC = prob->QC;
    
    const int om = prob->m;
    
    *flag = 0;
    
    
    OPT_setAllArray(n, hp, 0.0);
    
    
    //considering obj Q
    if(realObjF != 0.0 && Q.getNumberOfElements() > 0)
    {
        Q.evalTimesxt(p, hp, true, realObjF);
    }
    
    
    //considering QC for each constraint (we calculate lagrangian multiplier also)
    int k = 0;
    for(int i = 0; i < om; i++)
    {
        
        if(!constrEval[i])
        {
            mylambda[i] = 0.0;
            continue;
        }
        
        
        
        mylambda[i] = lambda[k]*scalec[k];
        
        if(lc[i] > -MIP_INFINITY)
        {
            if(uc[i] >= MIP_INFINITY)
                mylambda[i] = -mylambda[i]; //lower bound constraint. We multiply by minus one
        }
        else
        {
            if(uc[i] >= MIP_INFINITY)
                mylambda[i] = 0.0;//free constraint
        }
        
        
        if(QC[i].getNumberOfElements() == 0)
        {
            if(!nlConstr[i])
                mylambda[i] = 0.0; //linear constraint
        }
        else //if(QC[i].getNumberOfElements() > 0)
        {
            if(mylambda[i] != 0.0)
            //we have to multiply hessian by p and put the result in hp.
                QC[i].evalTimesxt(p, hp, true, mylambda[i]);
        }
        
        if(mylambda[i] != 0.0)
            someNzLambda = true;
        
        k++;
    }
    
    #if OPT_DEBUG_MODE
        assert(k == m);
    #endif
    
    
    if((prob->hasNlObj && realObjF != 0.0) || (prob->hasNlConstrs && someNzLambda))
    {
        //prob->objFactor is already considered in prob->nlpHessianEval...
        
        int r = prob->nlpHessianEval(thnumber, true, x, prob->hasNlObj ? scalef * nlp->in_nl_obj_factor : 0.0, mylambda, lagH);
        if(r != 0)
        {
            #if OPT_PRINT_CALLBACK_ERROR_MSG
                OPT_PRINTERRORNUMBER(r);
            #endif
            *flag = r;
            return;
        }
        
        //lagrangian multipliers were already considered at hessian evaluation
        lagH.evalTimesxt(p, hp, true);
    }
    
    
    //std::cerr << "Sai de OPT_algencanEvalHessVector\n";
}






//void myevalf(int n, double *x, double *f, int *flag)
static void OPT_algencanObjEval(int n, double *x, double *f, int *flag)
{
    OPT_Algencan::OPT_MyData *algData = threadAlgData[std::this_thread::get_id()];
    const unsigned int thnumber = algData->thnumber;
    //const bool *constrEval = algData->constrEval;
    OPT_Algencan * nlp = algData->nlp;
    OPT_MINLPProb *prob = &(nlp->prob);
    
    *flag = 0;
    
    int r;
    
    r = prob->objEval(thnumber, true, x, *f, nlp->in_nl_obj_factor);
    if(r != 0)
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            OPT_PRINTERRORNUMBER(r);
        #endif
        *flag = r;
        return;
    }
    
    
    //assert(false);
}


//void myevalg(int n, double *x, double *g, int *flag)
static void OPT_algencanObjGradEval(int n, double *x, double *g, int *flag)
{
    OPT_Algencan::OPT_MyData *algData = threadAlgData[std::this_thread::get_id()];
    const unsigned int thnumber = algData->thnumber;
    //const bool *constrEval = algData->constrEval;
    OPT_Algencan *nlp = algData->nlp;
    OPT_MINLPProb *prob = &nlp->prob;
    //const int orign = prob->n;
    
    
    *flag = 0;
    
    int r = prob->objGradEval(thnumber, true, x, g, nlp->in_nl_obj_factor);
    if(r != 0)
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            OPT_PRINTERRORNUMBER(r);
        #endif
        *flag = r;
        return;
    }
    
    //seting zero for auxiliary variables...
    if(n > prob->n)
        OPT_setAllArray(n - prob->n, &g[prob->n], 0.0);
    
    
    //assert(false);
}


// myevalh
static void OPT_algencanObjHessEval(int n, double *x, int *hrow, int *hcol, double *hval, int *hnnz, int lim, ALG_bool *lmem, int *flag)
{
    *flag = -1;
    assert(false);
}


// myevalc
static void OPT_algencanConstrEval(int n, double *x, int ind, double *c, int *flag) {
*flag = -1;
assert(false);
}

// myevaljac
static void OPT_algencanJacEval(int n, double *x, int ind, int *jcvar, double *jcval, int *jcnnz, int lim, ALG_bool *lmem, int *flag) {
*flag = - 1;
assert(false);
}

// myevalhc
static void OPT_algencanConstrHessEval(int n, double *x, int ind, int *hcrow, int *hccol, double *hcval, int *hcnnz, int lim, ALG_bool *lmem, int *flag) {
*flag = - 1;
assert(false);
}

// myevalgjacp
static void OPT_algencanEvalgjacp(int n, double *x, double *g, int m, double *p, double *q, char work, ALG_bool *gotj, int *flag) {
*flag = -1;
assert(false);
}



// myevalhl
static void OPT_algencanEvalhl(int n, double *x, int m, double *lambda, double scalef, double *scalec, int *hlrow, int *hlcol, double *hlval, int *hlnnz, int lim, ALG_bool *lmem, int *flag)
{
    OPT_Algencan::OPT_MyData *algData = threadAlgData[std::this_thread::get_id()];
    const unsigned int thnumber = algData->thnumber;
    //const bool *constrEval = algData->constrEval;
    OPT_Algencan *nlp = algData->nlp;
    OPT_MINLPProb *prob = &nlp->prob;
    double *auxVars = algData->auxVars;
    
    const bool *constrEval = algData->constrEval;
    const int nNzRowsLagH = algData->nNzRowsLagH;
    const int *nzRowsLagH = algData->nzRowsLagH;
    const int *sizeColsNzLagH = algData->sizeColsNzLagH;
    int* const *colsNzRowLagH = algData->colsNzRowLagH;
    
    const double realObjF = prob->objFactor * scalef;
    
    const bool evalLag = (prob->hasNlObj && realObjF != 0.0) || prob->hasNlConstrs;
    
    
    const bool *nlConstr = prob->nlConstr;
    const double *lc = prob->lc;
    const double *uc = prob->uc;
    double *mylambda = algData->auxConstr;
    
    
    
    MIP_SparseMatrix &lagH = prob->lagH;
    const MIP_SparseMatrix &Q = prob->Q;
    const MIP_SparseMatrix *QC = prob->QC;
    
    const MIP_SparseMatrix *singleM = NULL;
    
    bool someNzLambda = false;
    
    const int mquad = algData->mquad;
    const int *quadIndex = algData->quadIndex;
    int newmquad = 0;
    
    int indQuad, hessIndex;
    
    
    
    *flag = 0;
    *lmem = 0;
    *hlnnz = 0; //we need initialize it for the case where we have no memeory enough
    
    if( mquad > 0 || evalLag )
    {
        const int om = prob->m;
        
        int k = 0;
        for(int i = 0; i < om; i++)
        {
            if(!constrEval[i])
            {
                mylambda[i] = 0.0;
                continue;
            }
            
            
            mylambda[i] = lambda[k]*scalec[k];
            
            if(lc[i] > -MIP_INFINITY)
            {
                if(uc[i] >= MIP_INFINITY)
                    mylambda[i] = -mylambda[i]; //lower bound constraint. We multiply by minus one
            }
            else
            {
                if(uc[i] >= MIP_INFINITY)
                    mylambda[i] = 0.0;//free constraint
            }
            
            
            if(QC[i].getNumberOfElements() == 0)
            {
                if(!nlConstr[i])
                    mylambda[i] = 0.0; //linear constraint
            }
            else //if( QC[i].getNumberOfElements() > 0 )
            {
                if(mylambda[i] != 0.0)
                {
                    newmquad++;
                    indQuad = i;
                }
            }
            
            if(mylambda[i] != 0.0)
                someNzLambda = true;
            
            k++;
        }
        
        #if OPT_DEBUG_MODE
            assert(k == m);
        #endif
    }
    
    
    
    
    
    if( (prob->hasNlObj && realObjF != 0.0) || (prob->hasNlConstrs && someNzLambda) )
    {
        //prob->objFactor is already considered in prob->nlpHessianEval...
        
        int r = prob->nlpHessianEval(thnumber, true, x, prob->hasNlObj ? scalef * nlp->in_nl_obj_factor: 0.0, mylambda, lagH);
        if(r != 0)
        {
            #if OPT_PRINT_CALLBACK_ERROR_MSG
                OPT_PRINTERRORNUMBER(r);
            #endif
            *flag = r;
            return;
        }
        
        
        if( (Q.getNumberOfElements() == 0 || realObjF == 0.0) && newmquad == 0) //we only have the nonlinear part
        {
            singleM = &lagH;
        }
    }
    else
    {
        if(newmquad == 0)
        { //so we have only Q (maybe no Q)
            
            if(realObjF == 0.0)
            {
                *hlnnz = 0; //nothing to evaluate
                goto termination; 
            }
            
            singleM = &Q;
        }
        
        
        if( newmquad == 1 && (Q.getNumberOfElements() == 0 || realObjF == 0.0) )
        {//we have only one quadratic constraint (havin nonzero lambda) in the hessian
            singleM = &QC[indQuad]; //indquad has the index of this unique quadratic constraint
        }
    }
    
    
    if(singleM)
    {
        const int nzs = singleM->getNumberOfElements();
        
        if(nzs > lim)
        {
            *lmem = 1;
            goto termination;
        }
        
        *hlnnz = nzs;
        singleM->getStructure(hlrow, hlcol);
        singleM->getValues(hlval);
        
        goto termination;
    }
    
    
    //so, we have hessian composed by more than one matrix
    hessIndex = 0;
    
    OPT_setAllArray<double>(n, auxVars, 0.0); //auxVars is reinitialized inside the loop
    
    
    {
        //int k = 0;
        
        //nzRowsLagH is a unordered set. So, we will not run this in the order, but there is no problem about that... :)
        //for(auto it = nzRowsLagH->begin(); it != nzRowsLagH->end(); ++it)
        for(int k = 0; k < nNzRowsLagH; k++)
        {
            const int i = nzRowsLagH[k]; //*it; //row index
            
            const int ncols = sizeColsNzLagH[k];
            const int *pcols = colsNzRowLagH[k];
            
            
            //initializing nonzero positions in the current hess row
            for(int j = 0; j < ncols; j++)
                auxVars[ pcols[j] ] = 0.0;
            
            
            MIP_completeLagHessianRow(*prob, mquad, quadIndex, lagH, scalef, mylambda, i, auxVars, false);
            
            
            OPT_setAllArray(ncols, &hlrow[hessIndex], i);
            
            for(int j = 0; j < ncols; j++)
            {
                const int col = pcols[j];
                
                //hlrow[hessIndex] = i;
                hlcol[hessIndex] = col;
                hlval[hessIndex] = auxVars[col];
                
                //auxVars[col] = 0.0; //we have to initialize auxVars
                
                hessIndex++;
            }
            
            //k++;
        }
    }
    
    
    *hlnnz = hessIndex;
    
    
    /*for(int w = 0; w < hessIndex; w++)
    {
        std::cout << "hlrow["<<w<<"]: " << hlrow[w] << " hlcol["<<w<<"]: " << hlcol[w] << " hlval["<<w<<"]: " << hlval[w] << "\n";
    }
    
    OPT_getchar();*/
    
    
termination:
    
    //assert(false);
    return;
}




}







OPT_Algencan::OPT_Algencan():OPT_MyNLPSolver()
{
    initialize();
}


OPT_Algencan::~OPT_Algencan()
{
    deallocateSolverEnv();
}


void OPT_Algencan::deallocateSolverEnv() 
{
    threadAlgData.erase( std::this_thread::get_id() );
    OPT_MyNLPSolver::deallocateSolverEnv();
}


void OPT_Algencan::clearParameterList()
{
    params.clear();
}


int OPT_Algencan::getNumberOfIterations(long unsigned int &niter) 
#if OPT_HAVE_ALGENCAN
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



OPT_LISTSOLVERS OPT_Algencan::getSolverCode() 
{
    return optsolvers::OPT_ALGENCAN;
}



//that function can receive estimatives of maximum number of variables, constraints and nonzeros in quadratic terms matrices, but it is not mandatory...
int OPT_Algencan::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz) 
#if OPT_HAVE_ALGENCAN
{
    deallocateSolverEnv();
    clearParameterList();
    
    checkder = false;
    //move that part to initialize method
    epsfeas  = 1.0e-06;
    epsopt   = 1.0e-08;
    
    
    
    //params.push_back("ITERATIONS-OUTPUT-DETAIL 0"); //to turn off prints from algencan
    
    
    //std::cout << "initSolverEnv algencan threadNumber: " << threadNumber << OPT_GETFILELINE << "\n";
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Algencan::setMaxCPUTime(const double time) 
#if OPT_HAVE_ALGENCAN
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Algencan::setNumberOfThreads(const int nthreads) 
#if OPT_HAVE_ALGENCAN
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Algencan::setOutputLevel( const int level ) 
#if OPT_HAVE_ALGENCAN
{
    const char fileName[] = ".silent"; //there must be a file .silent in the current directory to algencan do not output
    
    
    if( level <= 0 )
    {
        FILE *file = fopen(fileName, "w"); 
        
        if( !file )
            return OPT_OPERATION_NOT_SUPPORTED;
        
        fclose(file);
        return 0;
    }
    else
    {
        remove(fileName); //we can get a error here if the file does not exits. So, we do not know if the error is because file does not exists or file could not be deleted
    }
    
    
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Algencan::setRelativeDualTol( const double tol ) 
#if OPT_HAVE_ALGENCAN
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Algencan::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_ALGENCAN
{
    epsopt = tol;
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Algencan::setRelativePrimalTol( const double tol ) 
#if OPT_HAVE_ALGENCAN
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Algencan::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_ALGENCAN
{
    char svalue[100];
    
    if(strlen(param) > 60)
        return OPT_BAD_INPUT;
    
    
    sprintf(svalue, "%s %0.15E", param, value);
    
    params.push_back( svalue );
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Algencan::setIntegerParameter(const char *param, const int value ) 
#if OPT_HAVE_ALGENCAN
{
    char svalue[100];
    
    if(strlen(param) > 60)
        return OPT_BAD_INPUT;
    
    
    sprintf(svalue, "%s %d", param, value);
    
    params.push_back( svalue );
    
    return 0;
    
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



/*int OPT_Algencan::setParameters( OPT_GeneralSolverParams &params )
#if OPT_HAVE_ALGENCAN
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif */



int OPT_Algencan::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_ALGENCAN
{
    char svalue[100];
    
    if(strlen(param) + strlen(value) > 90)
        return OPT_BAD_INPUT;
    
    
    sprintf(svalue, "%s %s", param, value);
    
    params.push_back( svalue );
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Algencan::setVariableType( const int index, const OPT_VARTYPE varType )
#if OPT_HAVE_ALGENCAN
{
    //anyway we set var type in MIP_MINLPProb...
    int r = OPT_MyNLPSolver::setVariableType(index, varType);
    
    
    //algencan is a continuous solver...
    if( varType == OPT_VT_INTEGER )
        return OPT_OPERATION_NOT_SUPPORTED;
    
    return r;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




#if 0
//some day we can need this method. Its set jacobian index in the unordered_map<unsigned int, unsigned int> 
unsigned int OPT_Algencan::buildJacIndex()
{
    const int n = prob.n;
    const int m = prob.m;
    
    unsigned int totalNz = 0;
    
    bool* auxIndex = (bool*) this->auxIndex;
    
    
    const minlpproblem::MIP_SparseMatrix &J = prob.J;
    const minlpproblem::MIP_SparseMatrix &A = prob.A;
    
    
    for(int i = 0; i < m; i++)
    {
        
        OPT_setAllArray(n, auxIndex, false);
        
        
        A.getRowStructure(i, auxIndex, NULL, true);
        J.getRowStructure(i, auxIndex, NULL, true);
        
        int nzQi;
        
        int r = prob.getNumberOfQuadCoefsInConstr(i, nzQi);
        #if OPT_DEBUG_MODE
            if(r != 0)
                assert(false);
        #endif
        
        
        if(nzQi > 0)
        {
            const minlpproblem::MIP_SparseMatrix &QCi = prob.QC[i];
            
            
            for( minlpproblem::MIP_SparseMatrix::RowIndexIterator it = QCi.beginRowIndex() ; it != QCi.endRowIndex() ; ++it )
            {
                unsigned int row = *it;
                
                auxIndex[row] = true;
                QCi.getRowStructure(row, auxIndex, NULL, true);
            }
        }
        
        for(int j = 0; j < n; j++)
        {
            if(auxIndex[j])
            {
                jacIndex[i*m + j] = totalNz;
                totalNz++;
            }
        }
    }
    
    return totalNz;
}
#endif






int OPT_Algencan::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_ALGENCAN
{
    //beggining of algencar vars
    /*coded[0]  fsub     */
    /*coded[1]  gsub     */
    /*coded[2]  hsub     */
    /*coded[3]  csub     */
    /*coded[4]  jacsub   */
    /*coded[5]  hcsub    */
    /*coded[6]  fcsub    */
    /*coded[7]  gjacsub  */
    /*coded[8]  gjacpsub */
    /*coded[9]  hlsub    */
    /*coded[10] hlpsub   */
    ALG_bool coded[11] = {1,1,0,0,0,0,1,1,0,1,0};
    //char outputfnm[] = "optsolvers_algencan.out";
    //char specfnm[] = "optsolvers_algencan.dat";
    
    
    char *outputfnm = OPT_getPointerFromConstPointer( this->outputfnm.c_str());
    
    char *specfnm = OPT_getPointerFromConstPointer( this->specfnm.c_str());
    
    
    int jcnnzmax = prob.J.getNumberOfElements() + prob.A.getNumberOfElements(); //estimative for maximum number of nz in jacobian
    int hnnzmax = OPT_max(prob.Q.getNumberOfElements(), prob.lagH.getNumberOfElements()); //estimative for maximum number of nz in hessian
    
    int nvparam = this->params.size(); //number of param set as string.
    char **vparam = NULL;
    
    ALG_bool *equatn = NULL, *linear = NULL;
    double *lambda = NULL, *x = NULL, *lx = NULL, *ux = NULL;
    
    //end of algencan vars
    
    
    bool hasFreeConstrs = false;
    int n, m = 0, ret;
    int totaln; //number of auxiliary variables
    int mquad = 0;
    int nNzRowsLagH;
    int *nzRowsLagH = NULL, *sizeColsNzLagH = NULL, **colsNzRowLagH = NULL;
    int **colsNzRowJac = NULL, *sizeColsNzRowJac = NULL;
    
    int *quadIndex = auxIndex;
    
    //std::unordered_set<int> setNzRowsLagH;
    //std::unordered_set<int> *colsNzRowHess = NULL;
    
    const int om = prob.m;
    const double *lc = prob.lc, *uc = prob.uc;
    const double *olx = prob.lx, *oux = prob.ux;
    const MIP_SparseMatrix &Q = prob.Q;
    const MIP_SparseMatrix *QC = prob.QC;
    //const MIP_SparseMatrix &lagH = prob.lagH;
    
    
    getNumberOfVars(n);
    //getNumberOfConstraints(m);
    
    
    if(resetSol)
    {
        //origSolverRetCode = INT_MAX;
        this->resetSol();
    }
    else
    {
        feasSol = false;
    }
    
    
    for(int i = 0; i < om; i++)
    {
        auxCEval[i] = lc[i] > -MIP_INFINITY || uc[i] < MIP_INFINITY;
        
        if(auxCEval[i])
        {
            if( QC[i].getNumberOfElements() > 0)
            {
                quadIndex[mquad] = i;
                mquad++;
            }
            m++;
        }
        else
        {
            hasFreeConstrs = true;
        }
    }
    
    
    
    vparam = (char **) malloc(nvparam *sizeof(char*));
    lambda = (double*) malloc(m *sizeof(double));
    equatn = (ALG_bool*) malloc(m *sizeof(ALG_bool));
    linear = (ALG_bool*) malloc(m *sizeof(ALG_bool));
    
    if(!vparam || !lambda || !equatn || !linear)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        
        retCode = OPT_MEMORY_ERROR;
        goto termination;
    }
    
    
    
    {
        const char *pchar; //to auxiliary to convert a const char in a char
        
        for(int i = 0; i < nvparam; i++)
        {
            pchar = params[i].c_str();
            vparam[i] = *( (char **) ( (void **) &pchar ) );
        }
    }
    
    
    {
        int k = 0;
        //calculating number of aux variables and set equatn
        totaln = n;
        for(int i = 0; i < om; i++)
        {
            if(auxCEval[i])
            {
                if(lc[i] > -MIP_INFINITY && uc[i] < MIP_INFINITY)
                {
                    equatn[k] = true;
                    if(lc[i] != uc[i])
                        totaln++;
                }
                else
                {
                    equatn[k] = false;
                }
                k++;
            }
        }
        
        #if OPT_DEBUG_MODE
            assert(k == m);
        #endif
    }
    
    
    //setting linear and calculating hnnzmax
    {
        MIP_SparseMatrix *QC = prob.QC;
        const bool *nlConstr = prob.nlConstr;
        
        int k = 0;
        
        for(int i = 0; i < om; i++)
        {
            if(auxCEval[i])
            {
                const int nzQ = QC[i].getNumberOfElements();
                
                hnnzmax += nzQ;
                jcnnzmax += OPT_min(n, 2*nzQ); //we consider quadratic part could set consider all variables in this constraints. (the maximum number could se the nozo zeros *2 since each coefficient involves at most two variables)
                
                linear[k] = nzQ == 0 && !nlConstr[i];
                k++;
            }
        }
        
        #if OPT_DEBUG_MODE
            assert(k == m);
        #endif
        
        
        //hnnzmax cannot be greater than number of elements in lower triangle
        
        //we can get an overflow about n*m or (n*n -n)/2 + n. So, we put the maximum between this value and zero. (If an overflow is gotten, this value will be negative)
        
        
        const int v1 = (n*n -n)/2 + n;
        const int v2 = n*m;
        
        if( v1 > 0 ) //check if an overflow did occurr
            hnnzmax = OPT_min(hnnzmax, v1);
        
        if( v2 > 0 ) //check if an overflow did occurr
            jcnnzmax = OPT_min(jcnnzmax, v2) + (totaln - n); //we have to consider the auxiliary variables also
            
        if(jcnnzmax < 0) //overflow
            jcnnzmax = INT_MAX;
        if(hnnzmax < 0) //overflow
            hnnzmax = INT_MAX;
    }
    
    
    if(m > 0)
    {
        if(std::isnan(lambdaInit[0]))
        {
            OPT_setAllArray(m, lambda, 0.0);
        }
        else
        {
            if(hasFreeConstrs)
            {
                int k = 0;
                for(int i = 0; i < om; i++)
                {
                    if(auxCEval[i])
                    {
                        lambda[k] = lambdaInit[i];
                        k++;
                    }
                }
                
                #if OPT_DEBUG_MODE
                    assert(k == m);
                #endif
            }
            else
            {
                OPT_copyArray(m, lambdaInit, lambda);
            }
        }
    }
    
    
    x = (double*) malloc(totaln *sizeof(double));
    lx = (double*) malloc(totaln *sizeof(double));
    ux = (double*) malloc(totaln *sizeof(double));
    
    if(!x || !lx || !ux)
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTMEMERROR;
        #endif
        
        retCode = OPT_MEMORY_ERROR;
        goto termination;
    }
    
    
    OPT_copyArray(n, prob.lx, lx);
    OPT_copyArray(n, prob.ux, ux);
    
    //lx and ux for auxiliary variables
    {
        int k = n;
        int w = 0;
        for(int i = 0; i < om; i++)
        {
            if(auxCEval[i])
            {
                if(lc[i] > -MIP_INFINITY && uc[i] < MIP_INFINITY && lc[i] != uc[i])
                {
                    lx[k] = lc[i];
                    ux[k] = uc[i];
                    k++;
                }
                w++;
            }
        }
        
        #if OPT_DEBUG_MODE
            assert(w == m);
            assert(k == totaln);
        #endif
    }
    
    
    
    //seting initial solution
    if(std::isnan(xInit[0]))
    {
        for(int i = 0; i < n; i++)
            x[i] = olx[i] > -MIP_INFINITY ? olx[i] : (oux[i] < MIP_INFINITY ? oux[i] : 0) ;
    }
    else
    {
        OPT_copyArray(n, xInit, x);
    }
    
    //setting initial values for auxiliary variables
    OPT_setAllArray(totaln -n, &x[n], 0.0);
    
    
    //seting auxCEval (we do not evaluate free constraints) and adding nonzero row indices to nzRowsLagH
    
    //OPT_addNzRowToContainer(prob.lagH, nzRowsLagH);
    //OPT_addNzRowToContainer(prob.Q, nzRowsLagH);
    
    
    
    
    //testing if we have only one matrix composing jacobian
    {
        const OPT_SparseMatrix *singleJ;
        
        int r = OPT_getSingleJacPointerOrCalcJacIndices( prob, mquad, hasFreeConstrs, auxCEval, singleJ, jcnnzmax, sizeColsNzRowJac, colsNzRowJac);
        
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            retCode = r;
            goto termination;
        }
    }
    
    
    
    //testing if we have a hessian composed by more than one sparse matriz
    {
        const int hasH = prob.hasNlConstrs || prob.hasNlConstrs ? 1 : 0;
        const int hasQ = Q.getNumberOfElements() > 0 ? 1 : 0;
        
        if(hasQ + mquad + hasH > 1)
        {
            std::unordered_set<int> setNzRowsLagH;
            
            OPT_fillNzRowsLagH(prob, auxCEval, mquad, quadIndex, setNzRowsLagH);
            
            nNzRowsLagH = setNzRowsLagH.size();
            
            nzRowsLagH = (int*) malloc( nNzRowsLagH * sizeof(int) );
            
            if(!nzRowsLagH)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                
                retCode = OPT_MEMORY_ERROR;
                goto termination;
            }
            
            {
                int i = 0;
                for(auto it: setNzRowsLagH)
                {
                    nzRowsLagH[i] = it;
                    i++;
                }
            }
            
            
            sizeColsNzLagH = (int*) malloc( nNzRowsLagH * sizeof(int) );
            colsNzRowLagH = (int**) calloc( nNzRowsLagH , sizeof(int*) );
            
            if(!sizeColsNzLagH || !colsNzRowLagH)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                
                retCode = OPT_MEMORY_ERROR;
                goto termination;
            }
            
            
            ret = OPT_fillColsNzRowHess(prob, auxCEval, mquad, quadIndex, nNzRowsLagH, nzRowsLagH, colsNzRowLagH, sizeColsNzLagH);
            
            if(ret != 0)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(ret);
                #endif
                
                retCode = OPT_MEMORY_ERROR;
                goto termination;
            }
            
            hnnzmax = 0;
            for(int k = 0; k < nNzRowsLagH; k++)
                hnnzmax += sizeColsNzLagH[k];
            
        }
        
    }
    
    
    
    
    threadAlgData[std::this_thread::get_id()] = &data; // we must put it here because we can have a thread having more to one algencan object to solve (and object can be realloc and so change address also).
    
    data.hasFreeConstrs = hasFreeConstrs;
    data.thnumber = threadNumber;
    data.mquad = mquad;
    data.constrEval = auxCEval;
    data.quadIndex = quadIndex;
    data.auxVars = auxValues;
    data.auxConstr = auxValues2;
    data.nlp = this; // we must put it here because we can have a thread having more to one algencan object to solve...
    
    data.sizeColsNzRowJac = sizeColsNzRowJac;
    data.colsNzRowJac = colsNzRowJac;
    
    data.nNzRowsLagH = nNzRowsLagH;
    data.nzRowsLagH = nzRowsLagH;
    data.sizeColsNzLagH = sizeColsNzLagH;
    data.colsNzRowLagH = colsNzRowLagH;
    
    
    
    
    efacc = sqrt(epsfeas);
    eoacc = sqrt(epsopt);
    
    efstin = sqrt( epsfeas );
    eostin = pow( epsopt, 1.5 );
    
    feasSol = false;
    
    
    
    c_algencan(
        //function pointers
        OPT_algencanObjEval, OPT_algencanObjGradEval, OPT_algencanObjHessEval, OPT_algencanConstrEval, OPT_algencanJacEval, OPT_algencanConstrHessEval, OPT_algencanEvalObjAndConstrs, OPT_algencanEvalObjAndConstrGrads, OPT_algencanEvalgjacp, OPT_algencanEvalhl, OPT_algencanEvalHessVector,
        
        jcnnzmax, hnnzmax, &epsfeas, &epsopt, &efstin, &eostin, &efacc, &eoacc, outputfnm, specfnm,
        
        nvparam, vparam, totaln, x, lx, ux, m, lambda, equatn, linear, coded, checkder,
        &objValue, &cnorm, &snorm, &nlpsupn, &inform
    );
    
    
    switch(inform)
    {
        case 0:
            retCode = OPT_OPTIMAL_SOLUTION;
            feasSol = true;
            break;
            
        case 4:
            retCode = OPT_INFEASIBLE_PROBLEM; //maybe it is not infeasible... :(
            break;
            
        case 1: //max iterations
        case 2: //max iner iterations
            
            #if OPT_PRINT_MAX_ITER_WARNING
                if( numberOfWarningsByIterLimit < maxNumberOfWarningsByIterLimit )
                {
                    std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Algencan solving!\n";
                    numberOfWarningsByIterLimit++;
                    
                    if( numberOfWarningsByIterLimit == maxNumberOfWarningsByIterLimit )
                        std::cerr << OPT_PREPRINT "Warning: Maximum number of warnings by maximum iteration achieved! Stopping these warnings.\n";
                }
            #endif
            
            retCode = OPT_MAX_ITERATIONS;
            break;
        
        case 3:
            retCode = OPT_MAX_EVALUATIONS;
            break;
        
        case -95:
            retCode = OPT_INFEASIBLE_PROBLEM;
            break;
            
        case -90:
        case -91:
        case -92:
        case -96:
            retCode = OPT_CALLBACK_FUNCTION_ERROR;
            break;
        
        case -93:
        case -94:
            retCode = OPT_MEMORY_ERROR;
            break;
        
        default:
            retCode = OPT_UNDEFINED_ERROR;
            break;
    }
    
    
    origSolverRetCode = inform;
    
    
    if(prob.objFactor < 0)
        objValue = -objValue;
    
    
    
    if(storeSol)
    {
        OPT_copyArray(n, x, sol);
    }
    
    
    //checking feasibility and storing constraints
    {
        double *pconstr = storeConstrs ? constr : auxValues2;
        
        
        //algencan sometimes respond code 0 for feasible problems :(!!
        ret = prob.isFeasibleToConstraints( threadNumber, x, true, nullptr, epsfeas*10, epsfeas*10, feasSol, pconstr);
        
        if( ret != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(ret);
            #endif
                
            if(storeConstrs)
                OPT_setAllArray<double>(om, constr, NAN);
        }
        
        if(!feasSol && inform == 0)
            retCode = OPT_UNDEFINED_ERROR; //maybe should be OPT_INFEASIBLE_PROBLEM
        
    }
    
    
    
    
    
    if(storeDualSol)
    {
        if(hasFreeConstrs)
        {
            int k = 0;
            for(int i = 0; i < om; i++)
            {
                if(auxCEval[i])
                {
                    dualSolC[i] = lambda[k];
                    k++;
                }
                else
                {
                    dualSolC[i] = 0.0;
                }
            }
            #if OPT_DEBUG_MODE
                assert(k == m);
            #endif
        }
        else
        {
            OPT_copyArray(om, lambda, dualSolC);
        }
    }
    
    
    
    
    
    if(!feasSol && retCode != OPT_INFEASIBLE_PROBLEM)
    {
        feasSol = prob.isConstrValuesFeasible(in_absolute_feas_tol, in_relative_feas_tol, constr);
    }
    
    
    
termination:
    
    if(vparam)		free(vparam);
    if(equatn)		free(equatn);
    if(linear)		free(linear);
    if(lambda)		free(lambda);
    if(x)			free(x);
    if(lx)			free(lx);
    if(ux)			free(ux);
    //if(colsNzRowHess)	delete[] colsNzRowHess;
    
    if(nzRowsLagH)	free(nzRowsLagH);
    if(sizeColsNzLagH)	free(sizeColsNzLagH);
    if(sizeColsNzRowJac)	free(sizeColsNzRowJac);
    
    
    if(colsNzRowJac)
    {
        for(int i = 0; i < om; i++)
        {
            if(colsNzRowJac[i])
                free(colsNzRowJac[i]);
        }
        
        free(colsNzRowJac);
    }
    
    
    if(colsNzRowLagH)
    {
        for(int i = 0; i < nNzRowsLagH; i++)
        {
            if(colsNzRowLagH[i])
                free(colsNzRowLagH[i]);
        }
        
        free(colsNzRowLagH);
    }
    
    
    return retCode;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif








