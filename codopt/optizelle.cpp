

/* OPTIZELLE INTERFACE IS TO POOR:
* 
* 1 - No way to input box constraints directly. You must formulate each box constraint as TWO linear inequalities, including providing derivatives for box constraints!
* 
* 2 - No way to put nonlinear inequalities. You must formulate it as equality adding a slack variable.
* 
* 3 - Constraints should have rhs 0.0. You must reformulate your constraint put your RHS in the left side
* 
* 4 - In the callback evaluation, no way to sinalize if you are evaluating the same point from last call
* 
* 5 - No way to set parameter before insert the initial solution. That os bad or our optsolvers design where you are able to set parameters before input the problem
* 
* 6 - A bug of bad_allocate in my first tests. Joe told how fo tix it, but, if they know there is a bug like that, why not they just do not corrrect it?
* 
* 7 - After solving the bug above, optizelle send me data to evaluate having infinity and nan, i.e, another bug...
* 
* 8 - After the optimization there is no status to sinalize optimal solution or infeasibility. You have to interpret a very poor set of optimization status. That is bad, because you have to check the optimality conditions by yourself, including feasibility?!
* 
* 9 - There are different optimizers to the problem, deppending on the constraint sets nature. That is bug. There should have a unique routine to optimize when you just set the parameter choosing optimziers. THERE IS NO a base class to derive the optimizers and that make things a little bit more difficult
* 
* 
* 
* For those rasons, I am giving up to incorporate optizelle in optsolvers.
*/ 





#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"


#include <cmath>
#include <cstdlib>
#include <climits>
#include <new>
#include <vector>
#include <iostream>
#include <iomanip>


//using namespace std;
using namespace optsolvers;
using namespace newspm;
using namespace minlpproblem;




#if OPT_HAVE_OPTIZELLE


//using namespace Optizelle;
using Optizelle::Natural;



inline void OPT_preMultiplyNonSimmetrycMatrixByVectorToOptizelle( const OPT_SparseMatrix &M, const unsigned int mLinEqNLIneq, const unsigned int *linEqNLIneqConstrs, const double *dy, double *r)
{
    for(unsigned int i = 0; i < mLinEqNLIneq; i++)
    {
        const auto row = linEqNLIneqConstrs[i];
        
        const auto rnElements = M.getNumberOfElementsAtRow(row);
        const auto *rcols = M.getRowColsPointer(row);
        const auto *rvalues = M.getRowValuesPointer(row);
        
        for(unsigned int j = 0; j < rnElements; j++)
        {
            const auto col = rcols[j];
            const auto value = rvalues[j];
            
            r[col] += dy[i] * value; //here we have to use i instead of row because dy does not necessarilly consider all rows in our problem
        }
        
    }
    
}




/* defining objective function object*/

struct optsolvers::OPT_ObjToOptizelle: public Optizelle::ScalarValuedFunction<double, Optizelle::Rm> {
    
    typedef Optizelle::Rm <double> X;
    
    OPT_Optizelle *nlp;
    OPT_MINLPProb *prob;
    
    
    OPT_ObjToOptizelle(OPT_Optizelle *nlp, OPT_MINLPProb *prob) : Optizelle::ScalarValuedFunction<double, Optizelle::Rm>()
    {
        this->nlp = nlp;
        this->prob = prob;
    }
    
    
    double eval(const X::Vector& x) const 
    {
        const bool newx = true;
        const unsigned int thnumber = nlp->threadNumber;
        double objvalue;
        
        //printf("entering at OPT_ObjToOptizelle.eval\n");
        //OPT_getchar();
        
        int r = prob->objEval(thnumber, newx, x.data(), objvalue, nlp->in_nl_obj_factor);
        
        if( r != 0 )
        {
            OPT_PRINTERRORMSG("Error " << r << " at objective evaluation. Unfortunatelly Optizelle has no tool to inform the error");
            
            objvalue = NAN;
        }
        
        return objvalue;
    }
    
    
    void grad(X::Vector const & x, X::Vector & grad) const 
    {
        const bool newx = true;
        const unsigned int thnumber = nlp->threadNumber;
        const unsigned int n = prob->n;
        const unsigned int noptizelle = grad.size();
        
        double *agrad = grad.data();
        const double *ax = x.data();
        
        //printf("entering at OPT_ObjToOptizelle.grad. n: %u noptizelle: %u\n", n, noptizelle);
        //OPT_getchar();
        
        
        int r = prob->objGradEval(thnumber, newx, ax, agrad, nlp->in_nl_obj_factor);
        if( r != 0 )
        {
            OPT_PRINTERRORMSG("Error " << r << " at objective's gradient evaluation. Unfortunatelly Optizelle has no tool to inform the error");
            
            OPT_setAllArray<double>(noptizelle, agrad, NAN);
            return;
        }
        
        if( noptizelle > n )
            OPT_setAllArray<double>( noptizelle - n, &agrad[n], 0.0 );
        
        printf("obj grad: ");
        for(unsigned int i = 0; i < noptizelle; i++)
            printf("[%u]: %0.2lf \t", i, grad[i]);
        printf("\n\n");
    }
    
    
    // Hessian-vector product
    void hessvec( X::Vector const & x, X::Vector const & dx,
        X::Vector & H_dx ) const 
    {
        const bool newx = true;
        const unsigned int thnumber = nlp->threadNumber;
        const unsigned int n = prob->n;
        const double objFactor = prob->getObjFactor();
        const double *ax = x.data();
        const double *adx = dx.data();
        double *aH_dx = H_dx.data();
        
        auto &lagH = prob->lagH;
        
        //printf("entering at OPT_ObjToOptizelle.hessvec\n");
        //OPT_getchar();
        
        OPT_setAllArray<double>(H_dx.size(), aH_dx, 0.0);
        
        if( objFactor )
        {
            int r = prob->nlpHessianEval(thnumber, newx, ax, objFactor*nlp->in_nl_obj_factor, nullptr, lagH );
            
            if(r != 0)
            {
                OPT_PRINTERRORMSG("Error " << r << " at objective's hessian evaluation. Unfortunatelly Optizelle has no tool to inform the error");
                
                OPT_setAllArray<double>(H_dx.size(), aH_dx, NAN);
                return;
            }
            
            lagH.evalTimesxt(adx, aH_dx, true);
            
            if( prob->Q.getNumberOfElements() > 0 )
            {
                prob->Q.evalTimesxt(adx, aH_dx, true, objFactor);
            }
        }
        
        printf("obj hess vec: \n");
        for(unsigned int i = 0; i < H_dx.size(); i++)
            printf("x[%u]: %0.3lf, H_dx[%u]: %0.2lf \t", i, x[i], i, H_dx[i]);
        printf("\n\n");
    }
    
};


/*Optizelle requires all inequality constraints be linear. So, we must separate linear inequalities from remainder constraints.

Linear equalities will be put together nonlinear constraints. Nonlinear inequalities should be reformulated as equality constraints. (AFF!)
*/
struct optsolvers::OPT_IneqToOptizelle
    :public Optizelle::VectorValuedFunction<double, Optizelle::Rm, Optizelle::Rm>
{
    typedef Optizelle::Rm <double> X;
    typedef Optizelle::Rm <double> Y;
    
    
    OPT_Optizelle *nlp;
    OPT_MINLPProb *prob;
    
    
    OPT_IneqToOptizelle(OPT_Optizelle *nlp, OPT_MINLPProb *prob) :Optizelle::VectorValuedFunction<double, Optizelle::Rm, Optizelle::Rm>()
    {
        this->nlp = nlp;
        this->prob = prob;
    }
    
    
    // y=h(x)
    void eval( X::Vector const & x, Y::Vector & y ) const
    {
        const unsigned int n = prob->n;
        const unsigned int mLinIneq = nlp->mLinIneq;
        
        const unsigned int *linIneqConstrs = nlp->linIneqConstrs;
        
        const double *lc = prob->lc, *uc = prob->uc;
        
        const double *ax =  x.data();
        double *ay = y.data();
        
        //printf("entering at OPT_IneqToOptizelle.eval\n");
        //OPT_getchar();
        
        #if 0
        
        int r = prob->constraintsEval( thnumber, newx, evalLinearIneqConstrs, x.data(), auxConstr );
        
        if( r )
        {
             OPT_PRINTERRORMSG("Error " << r << " at linear inequality evaluation. Unfortunatelly Optizelle has no tool to inform the error");
             
             OPT_setAllArray<double>(y.size(), ay, NAN);
        }
        
        //optizelle only accepts constraints in formart h(x) <= 0 (AFF!!!!)
        
        unsigned int k = mLinIneq;
        const double *lc = prob->lc, *uc = prob->uc;
        for( unsigned int i = 0; i < mLinIneq; i++ )
        {
            const auto cind = linIneqConstrs[i];
            
            if( lc[cind] <= -OPT_INFINITY )
            { //we do not have the lower bound
                ay[i] = auxConstr[cind] - uc[cind];
            }
            else
            {
                ay[i] = -auxConstr[cind] + lc[cind]; //constarint should be in format h(x) <= 0
                
                if( uc[cind] < OPT_INFINITY )
                {
                    ay[k] = auxConstr[cind] - uc[cind];
                    k++;
                }
            }
        }
        
        #endif
        
        prob->constraintLinearPartEvaluation(mLinIneq, (const int*) linIneqConstrs, ax, ay);
        
        //optizelle only accepts constraints in formart h(x) >= 0 (AFF!!!!)
        
        unsigned int k = mLinIneq;
        
        for(unsigned int i = 0; i < mLinIneq; i++)
        {
            const auto cind = linIneqConstrs[i];
            
            if( lc[cind] <= -OPT_INFINITY )
            { //we do not have the lower bound
                y[i] =  -(y[i] - uc[cind]); 
            }
            else
            {
                if( uc[cind] < OPT_INFINITY )
                { //double bounded constraint
                    y[k] = -(y[i] - uc[cind]); //upper bound: do not put it after the lower bound constraint calculation
                    k++;
                }
                
                y[i] = y[i] - lc[cind]; //lower bound constraint: constarint should be in format h(x) >= 0
            }
        }
        
        //setting variable bounds
        {
            const double *lx = prob->lx;
            const double *ux = prob->ux;
            
            for(unsigned int i = 0; i < n; i++)
            {
                if( lx[i] > -OPT_INFINITY )
                {
                    y[k] = x[i] - lx[i];  //lower bound constraint: constarint should be in format h(x) >= 0
                    k++;
                }
                
                if( ux[i] < OPT_INFINITY )
                {
                    y[k] = -(x[i] - ux[i]);
                    k++;
                }
            }
        }
        
        
        //setting auxiliar variable bounds to inequalities constraints transformed in equalities
        {
            const unsigned int mNLIneq = nlp->mNLIneq;
            const unsigned int *NLIneqConstrs = nlp->NLIneqConstrs;
            
            for(unsigned int i = 0; i < mNLIneq; i++)
            {
                const auto cind = NLIneqConstrs[i];
                
                if( lc[cind] > -OPT_INFINITY )
                {
                    y[k] = x[n+i] - lc[cind]; //lower bound constraint: constraint should be in format h(x) >= 0
                    k++;
                }
                
                if( uc[cind] < OPT_INFINITY )
                {
                    y[k] = -(x[n+i] - uc[cind]);
                    k++;
                }
                
            }
            
        }
        
        #if OPT_DEBUG_MODE
            assert(k == y.size());
        #endif
        
        
        for(unsigned int i = 0; i < y.size(); i++)
            printf("ey[%u]: %0.3lf \t", i, y[i]);
        printf("\n\n");
        
    }
    
    
    // y=h'(x)dx
    void p( X::Vector const &x, X::Vector const &dx, Y::Vector &y) const 
    { //here, we must multiply the jacobiam (mxn) by an array nX1
        
        const unsigned int mLinIneq = nlp->mLinIneq;
        const unsigned int n = prob->n;
        //const bool *evalLinearIneqConstrs = nlp->evalLinearIneqConstrs;
        
        const unsigned int *linIneqConstrs = nlp->linIneqConstrs;
        
        const MIP_SparseMatrix &A = prob->A;
        
        const double *lc = prob->lc, *uc = prob->uc;
        
        const double *adx = dx.data();
        
        double *ay = y.data();
        
        
        //printf("entering at OPT_IneqToOptizelle.p\n");
        //OPT_getchar();
        
        
        unsigned int k = mLinIneq;
        for(unsigned int i = 0; i < mLinIneq; i++)
        {
            const auto cind = linIneqConstrs[i];
            
            const double rowdx = A.evalRowTimesxt(cind, adx);
            
            if( lc[cind] <= -OPT_INFINITY )
            {
                y[i] = -rowdx; //we only have upper bound (free constraints are not considered here). format is h(x) >= 0
            }
            else
            {
                y[i] = rowdx; //format is h(x) >= 0
                
                if( uc[cind] < OPT_INFINITY )
                {
                    y[k] = -rowdx;
                    k++;
                }
            }
        }
        
        //setting variable bounds
        {
            const double *lx = prob->lx, *ux = prob->ux;
            
            for(unsigned int i = 0; i < n; i++)
            {
                if( lx[i] > -OPT_INFINITY )
                {
                    y[k] = dx[i]; //1*adx: we have to multiply derivative by adx
                    k++;
                }
                
                if( ux[i] < OPT_INFINITY )
                {
                    y[k] = -dx[i]; //-1*adx[i]: we have to multiply derivative by adx
                    k++;
                }
            }
        }
        
        //setting auxiliar variable bounds to inequalities constraints transformed in equalities
        {
            const unsigned int mNLIneq = nlp->mNLIneq;
            const unsigned int *NLIneqConstrs = nlp->NLIneqConstrs;
            
            for(unsigned int i = 0; i < mNLIneq; i++)
            {
                const auto cind = NLIneqConstrs[i];
                
                if( lc[cind] > -OPT_INFINITY )
                {
                    y[k] = dx[n+i]; //1*adx[n+i]: we have to multiply derivative by adx
                    k++;
                }
                
                if( uc[cind] < OPT_INFINITY )
                {
                    y[k] = -dx[n+i]; //-1*adx[n+1]: we have to multiply derivative by adx
                    k++;
                }
                
            }
        }
        
        #if OPT_DEBUG_MODE
            assert(k == y.size());
        #endif
        
        for(unsigned int i = 0; i < y.size(); i++)
            printf("pdy[%u]: %0.3lf \t", i, y[i]);
        printf("\n\n");
    }
    
    
    // z=h'(x)*dy
    void ps( X::Vector const &x, Y::Vector const &dy, X::Vector &z) const
    {
        /*now, we must multiply transposed jacobian J' (nXm) by an array z (mx1)
         So, we multiply z' (1xm) by J (mxn) */
        
        const unsigned int mLinIneq = nlp->mLinIneq;
        const unsigned int n = prob->n;
        const unsigned int *linIneqConstrs = nlp->linIneqConstrs;
        
        const MIP_SparseMatrix &A = prob->A;
        
        const double *lc = prob->lc, *uc = prob->uc;
        
        const double *ady = dy.data();
        
        unsigned int k = mLinIneq;
        
        double *az = z.data();
        
        
        //printf("entering at OPT_IneqToOptizelle.ps\n");
        //OPT_getchar();
        
        
        for(unsigned int i = 0; i <  dy.size(); i++)
            printf("entrada-dy[%u]: %0.3f \t", i, dy[i]);
        printf("\n");
        
        OPT_setAllArray(z.size(), az, 0.0);
        
        for(unsigned int i = 0; i < z.size(); i++)
            printf("1z[%u]: %0.3lf \t", i, z[i]);
        printf("\n");
        
        for(unsigned int i = 0; i < mLinIneq; i++)
        {
            const auto row = linIneqConstrs[i];
            
            const auto rnElements = A.getNumberOfElementsAtRow(row);
            
            const auto *rcols = A.getRowColsPointer(row);
            const auto *rvalues = A.getRowValuesPointer(row);
            
            char nBounds; //1 ifconstr has only one bound and two if constraints is double bounded.
            double factor;
            
            if( lc[row] > -OPT_INFINITY )
            { //lower bound constraint. We must multiply to -1 because optizelle only accepets constraints in format h(x) >= 0  (aff)
                
                factor = 1.0;
                if(uc[row] < OPT_INFINITY)
                    nBounds = 2;
                else
                    nBounds = 1;
            }
            else
            { //we have only upper bound constraint
                factor = -1.0;
                nBounds = 1;
            }
            
            
            //Jacobian is not symmetric. So, we do not need worry about simmetry here
            
            const double factorxrow1 = factor * ady[i]; //here we have to use i instead of row because dy does not necessarilly consider all rows in our problem
            
            const double factorxrow2 = nBounds > 1 ? -factor* ady[k] : 0.0; //we have a double bound constraint and we have to consider upper bound also
            
            const double totalFactor = factorxrow1 + factorxrow2;
            
            for( unsigned int j = 0; j < rnElements; j++ )
            {
                const auto col = rcols[j];
                const auto value = rvalues[j];
                
                z[col] += totalFactor * value;
            }
            
            
            if( nBounds > 1 )
                k++;
        }
        
    
        
        for(unsigned int i = 0; i < z.size(); i++)
            printf("2z[%u]: %0.3lf \t", i, z[i]);
        printf("\n");
        
        
        //setting variable bounds
        {
            const double *lx = prob->lx;
            const double *ux = prob->ux;
            
            for(unsigned int i = 0; i < n; i++)
            {
                if( lx[i] > -OPT_INFINITY )
                {
                    z[i] += dy[k]; //we have to multiply derivative by adx
                    k++;
                }
                
                if( ux[i] < OPT_INFINITY )
                {
                    z[i] += -dy[k]; //we have to multiply derivative by adx
                    k++;
                }
            }
        }
        
        //setting auxiliar variable bounds to inequalities constraints transformed in equalities
        {
            const unsigned int mNLIneq = nlp->mNLIneq;
            const unsigned int *NLIneqConstrs = nlp->NLIneqConstrs;
            
            for(unsigned int i = 0; i < mNLIneq; i++)
            {
                const auto cind = NLIneqConstrs[i];
                
                if( lc[cind] > -OPT_INFINITY )
                {
                    z[n+i] += dy[k];
                    k++;
                }
                
                if( uc[cind] < OPT_INFINITY )
                {
                    z[n+i] += -dy[k];
                    k++;
                }
            }
        }
        
        #if OPT_DEBUG_MODE
            assert(k == dy.size());
        #endif
        
        printf("pdz: \n");
        for(unsigned int i = 0; i < x.size(); i++)
            printf("x[%u]: %0.3f \t", i, x[i]);
        printf("\n\n");
        for(unsigned int i = 0; i < z.size(); i++)
            printf("pdz[%u]: %0.3lf \t", i, z[i]);
        printf("\n\n");
        
    }
    
    
    // z=(h''(x)dx)*dy
    void pps( X::Vector const & x, X::Vector const & dx, Y::Vector const & dy, X::Vector & z) const 
    { //since constraints are linear, we do not have hessian here...
        
        //printf("entering at OPT_IneqToOptizelle.pps\n");
        //OPT_getchar();
        
        OPT_setAllArray<double>( z.size(), z.data(), 0.0 );
    }
    
};


/*Optizelle requires allinequalities be linear. So, we must reformulate nonlinear inequalities as equalities*/
struct optsolvers::OPT_EqToOptizelle : public Optizelle::VectorValuedFunction <double, Optizelle::Rm, Optizelle::Rm>
{
    typedef Optizelle::Rm <double> X;
    typedef Optizelle::Rm <double> Y;
    
    OPT_Optizelle *nlp;
    OPT_MINLPProb *prob;
    
    
    OPT_EqToOptizelle(OPT_Optizelle *nlp, OPT_MINLPProb *prob) :Optizelle::VectorValuedFunction<double, Optizelle::Rm,Optizelle::Rm>()
    {
        this->nlp = nlp;
        this->prob = prob;
    }
    
    // y = g(x)
    void eval(X::Vector const & x, Y::Vector & y) const
    {
        const bool newx = true;
        const unsigned int thnumber = nlp->threadNumber;
        const unsigned int n = prob->n;
        const unsigned int m = prob->m;
        const unsigned int mLinEqNLIneq = nlp->mLinEqNL;
        
        
        const bool *evalEqAndNLIneqConstrs = nlp->evalEqAndNLIneqConstrs;
        const unsigned int *linEqNLIneqConstrs = nlp->linEqNLConstrs;
        const double *lc = prob->lc, *uc = prob->uc;
       
        const double *ax = x.data();
        
        double *auxConstr = nlp->auxValues;
        double *ay = y.data();
        
        
        //printf("entering at OPT_EqToOptizelle.eval\n");
        //OPT_getchar();
        
        
        int r = prob->constraintsEval( thnumber, newx, evalEqAndNLIneqConstrs, ax, auxConstr );
        
        if( r )
        {
             OPT_PRINTERRORMSG("Error " << r << " at linear inequality evaluation. Unfortunatelly Optizelle has no tool to inform the error");
             
             OPT_setAllArray<double>(y.size(), ay, NAN);
             return;
        }
        
        for( unsigned int i = 0, j = n; i < mLinEqNLIneq; i++ )
        {
            const auto cind = linEqNLIneqConstrs[i];
            
            if( lc[cind] == uc[cind] )
            {  //equality constarint
                y[i] = auxConstr[cind] - lc[cind];
            }
            else
            { //inequality constraint. We have an artificial variable. Constraint bounds is considered as variable bound of artifical variable
                y[i] = auxConstr[cind] - x[j];
                j++;
            }
            
        }
        
        
    }
    
    // y=g'(x)dx
    void p(X::Vector const & x, X::Vector const & dx, Y::Vector & y) const 
    {  //here, we must multiply the jacobiam (mxn) by an array nX1
        const bool newx = true;
        const unsigned int thnumber = nlp->threadNumber;
        const unsigned int n = prob->n;
        const unsigned int m = prob->m;
        const unsigned int mLinEqNLIneq = nlp->mLinEqNL;
        
        
        const bool *evalEqAndNLIneqConstrs = nlp->evalEqAndNLIneqConstrs;
        const unsigned int *linEqNLIneqConstrs = nlp->linEqNLConstrs;
        const double *lc = prob->lc, *uc = prob->uc;
        
        const double *ax = x.data();
        const double *adx = dx.data();
        
        double *auxValues = nlp->auxValues;
        double *ay = y.data();
        
        const OPT_SparseMatrix *QC = prob->QC;
        const bool *nlConstr = prob->nlConstr;
        const OPT_SparseMatrix &A = prob->A;
        OPT_SparseMatrix &J = prob->J;
        
        
        //printf("entering at OPT_EqToOptizelle.p mLinEqNLIneq: %u y.size: %u\n", mLinEqNLIneq, (unsigned) y.size() );
        //OPT_getchar();
        
        
        int r = prob->nlJacobianEval(thnumber, true, evalEqAndNLIneqConstrs, ax, J);
        
        if( r )
        {
             OPT_PRINTERRORMSG("Error " << r << " at linear inequality evaluation. Unfortunatelly Optizelle has no tool to inform the error");
             
             OPT_setAllArray<double>(y.size(), ay, NAN);
             return;
        }
        
        
        OPT_setAllArray<double>(y.size(), ay, 0.0);
        
        
        for(unsigned int i = 0, j = n; i < mLinEqNLIneq; i++)
        {
            const auto cind = linEqNLIneqConstrs[i];
            
            if( nlConstr[cind] )
                y[i] = J.evalRowTimesxt(cind, adx);
            
            if( QC[cind].getNumberOfElements() > 0 )
            {
                //first we have to calculate w = QC*x
                QC[cind].evalTimesxt(ax, auxValues, false);
                
                //now, we must evaluate w*adx
                y[i] += OPT_innerProduct<double>(n, auxValues, adx);
            }
            
            if( A.getNumberOfElementsAtRow(cind) > 0 )
                y[i] += A.evalRowTimesxt(cind, adx); 
            
            if( lc[cind] != uc[cind] )
            { //we must consider the artificial variable in the constraint to turn it in an equality. Artificial variable having -1 coefficient
                y[i] += (-1 * dx[j]); // j is the index of artifical variable
                j++;
            }
        }
        
        
    }
    
    // z=g'(x)*dy
    void ps(X::Vector const & x, Y::Vector const & dy, X::Vector & z) const
    { 
        /*now, we must multiply transposed jacobian J' (nXm) by an array dy (mx1)
         So, we multiply dy' (1xm) by J (mxn) */
        
        const bool newx = true;
        //const bool hasQuadMatrixInConstrs = nlp->hasQuadMatrixInConstrs;
        const unsigned int thnumber = nlp->threadNumber;
        const unsigned int n = prob->n;
        const unsigned int m = prob->m;
        const unsigned int mLinEqNLIneq = nlp->mLinEqNL;
        
        
        const bool *evalEqAndNLIneqConstrs = nlp->evalEqAndNLIneqConstrs;
        const unsigned int *linEqNLIneqConstrs = nlp->linEqNLConstrs;
        const double *lc = prob->lc, *uc = prob->uc;
        
        const double *ax = x.data();
        const double *ady = dy.data();
        
        double *auxValues = nlp->auxValues;
        double *az = z.data();
        
        //const bool *nlConstr = prob->nlConstr;
        const OPT_SparseMatrix &A = prob->A;
        OPT_SparseMatrix &J = prob->J;
        const OPT_SparseMatrix *QC = prob->QC;
        
        
        //printf("entering at OPT_EqToOptizelle.ps\n");
        //OPT_getchar();
        
        
        int r = prob->nlJacobianEval(thnumber, true, evalEqAndNLIneqConstrs, ax, J);
        
        if( r )
        {
             OPT_PRINTERRORMSG("Error " << r << " at linear inequality evaluation. Unfortunatelly Optizelle has no tool to inform the error");
             
             OPT_setAllArray<double>(z.size(), az, NAN);
             return;
        }
        
        
        OPT_setAllArray<double>(z.size(), az, 0.0);
        
        OPT_preMultiplyNonSimmetrycMatrixByVectorToOptizelle( J, mLinEqNLIneq, linEqNLIneqConstrs, ady, az );
        
        OPT_preMultiplyNonSimmetrycMatrixByVectorToOptizelle( A, mLinEqNLIneq, linEqNLIneqConstrs, ady, az );
        
        
        for(unsigned int i = 0, j = n; i < mLinEqNLIneq; i++)
        {
            const auto cind = linEqNLIneqConstrs[i];
            
            if( QC[cind].getNumberOfElements() > 0 )
            {
                //first we have to calculate w = QC*x
                QC[cind].evalTimesxt(ax, auxValues, false);
                
                //now, we consider in the jacobian pre multiplying
                #pragma ivdep
                #pragma GCC ivdep
                for(unsigned int k = 0; k < n; k++)
                    z[k] += auxValues[k] * dy[i];
            }
            
            if( lc[cind] != uc[cind] )
            { //we must consider the artificial variable in the constraint to turn it in an equality. Artificial variable having -1 coefficient
                z[j] += (-1.0 * dy[i]) ; // j is the index of artifical variable
                j++;
            }
        }
        
    }
    
    
    // z=(h''(x)dx)*dy
    void pps( X::Vector const & x, X::Vector const & dx, Y::Vector const & dy, X::Vector & z) const
    {
        const bool hasQuadMatrixInConstrs = nlp->hasQuadMatrixInConstrs;
        const bool hasNLConstraints = prob->hasNLConstraints();
        
        const unsigned int mLinEqNLIneq = nlp->mLinEqNL;
        const unsigned int *linEqNLIneqConstrs = nlp->linEqNLConstrs;
        
        const double *adx = dx.data();
        const double *ady = dy.data();
        
        double *az = z.data();
        
        
        //printf("entering at OPT_EqToOptizelle.pps\n");
        //OPT_getchar();
        
        OPT_setAllArray<double>(z.size(), az, 0.0);
        
        
        if( hasNLConstraints )
        {
            const unsigned int thnumber = nlp->threadNumber;
            const bool newx = true;
            const double objFactor = 0.0;
            const unsigned int m = prob->m;
            const double *ax = x.data();
            
            double *lambda = nlp->auxValues2;
            OPT_SparseMatrix &lagH = prob->lagH;
            
            
            OPT_setAllArray<double>(m, lambda, 0.0);
            
            for(unsigned int i = 0; i < mLinEqNLIneq; i++)
                lambda[ linEqNLIneqConstrs[i] ] = dy[i];
            
            int r = prob->nlpHessianEval(thnumber, newx, ax, objFactor, lambda, lagH);
            
            if(r != 0)
            {
                OPT_PRINTERRORMSG("Error " << r << " at objective's hessian evaluation. Unfortunatelly Optizelle has no tool to inform the error");
                
                OPT_setAllArray<double>(z.size(), az, NAN);
                return;
            }
            
            lagH.evalTimesxt(adx, az, true);
        }
        
        
        if( hasQuadMatrixInConstrs )
        {
            const OPT_SparseMatrix *QC = prob->QC;
            
            for(unsigned int i = 0; i < mLinEqNLIneq; i++)
            {
                const auto cind = linEqNLIneqConstrs[i];
                const double lambdai = dy[i];
                
                if( QC[cind].getNumberOfElements() > 0 && lambdai != 0.0)
                {
                    QC[cind].evalTimesxt(adx, az, true, lambdai);
                }
            }
        }
        
        
    }
    
    
    
};

#endif



OPT_Optizelle::OPT_Optizelle():OPT_MyNLPSolver()
{
    initialize();
}
        
OPT_Optizelle::~OPT_Optizelle()
{
    deallocateMemory();
    deallocateSolverEnv();
}

void OPT_Optizelle::deallocateMemory()
{
    OPT_secFree(evalEqAndNLIneqConstrs);
    OPT_secFree(linIneqConstrs);
    OPT_secFree(linEqNLConstrs);
    OPT_secFree(NLIneqConstrs);
    
    OPT_MyNLPSolver::deallocateMemory();
}

void OPT_Optizelle::deallocateSolverEnv()
{
    OPT_MyNLPSolver::deallocateSolverEnv();
}


int OPT_Optizelle::getNumberOfIterations(long unsigned int &niter)
{
    OPT_PRINTERRORMSG("Number of iterations cannot be get from optizelle\n");
    return OPT_OPERATION_NOT_IMPLEMENTED;
}


OPT_LISTSOLVERS OPT_Optizelle::getSolverCode()
{
    return OPT_OPTIZELLE;
}


int OPT_Optizelle::getVariableType( const int index, OPT_VARTYPE &varType)
{
    varType = optsolvers::OPT_VT_CONTINUOUS;
    
    return 0;
}


void OPT_Optizelle::initialize()
{
    OPT_MyNLPSolver::initialize();
    
    //evalLinearIneqConstrs = nullptr;
    evalEqAndNLIneqConstrs = nullptr;
    linEqNLConstrs = nullptr;
    linIneqConstrs = nullptr;
    NLIneqConstrs = nullptr;
    
    mLinIneq = 0;
    //mLinEq = 0;
    mLinEqNL = 0;
    mNLIneq = 0;
    //mNLEq = 0;
}


int OPT_Optizelle::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
{
    return 0;
}


int OPT_Optizelle::setInitialSolution(const double* x, const double* dualConstrs, const double* dualVars)
{
    int n = 0, m = 0;
    
    getNumberOfVars(n);
    getNumberOfConstraints(m);
    
    
    if( xInit_vec.size() < (unsigned int) n )
        xInit_vec.resize(n);
    
        
    if( lambdaInit_vec.size() < (unsigned int) m )
        lambdaInit_vec.resize(m);
    
    if( zInit_vec.size() < (unsigned int) n )
        zInit_vec.resize(n);
    
    
    OPT_setInicialSolution(n, m, x, dualConstrs, dualVars, xInit_vec.data(), lambdaInit_vec.data(), zInit_vec.data());
    
    return 0;
}


int OPT_Optizelle::setMaxCPUTime(const double time)
{
    OPT_PRINTERRORMSG("Max CPU time is not being set in optzelle");
    return OPT_OPERATION_NOT_IMPLEMENTED;
}


int OPT_Optizelle::setNumberOfThreads(const int nthreads)
{
    OPT_PRINTERRORMSG("Number of threads is not being set in optzelle");
    return 0;
}


int OPT_Optizelle::setOutputLevel(const int level)
{
    OPT_PRINTERRORMSG("Output level is not being set in optzelle");
    return OPT_OPERATION_NOT_IMPLEMENTED;
}


int OPT_Optizelle::setRelativeDualTol( const double tol )
{
    OPT_PRINTERRORMSG("Relative dual tolerance is not being set in optzelle");
    return OPT_OPERATION_NOT_IMPLEMENTED;
}


int OPT_Optizelle::setRelativeOptimalityTol(const double tol)
{
    OPT_PRINTERRORMSG("Relative optimality tol is not being set in optzelle");
    return OPT_OPERATION_NOT_IMPLEMENTED;
}


int OPT_Optizelle::setRelativePrimalTol(const double tol)
{
    OPT_PRINTERRORMSG("Relative primal tol is not being set in optzelle");
    return OPT_OPERATION_NOT_IMPLEMENTED;
}


int OPT_Optizelle::setDoubleParameter(const char *param, const double value)
{
    OPT_PRINTERRORMSG("Double parameters are not being set in optzelle");
    return OPT_OPERATION_NOT_IMPLEMENTED;
}


int OPT_Optizelle::setIntegerParameter(const char *param, const int value )
{
    OPT_PRINTERRORMSG("Integer parameters are not being set in optzelle");
    return OPT_OPERATION_NOT_IMPLEMENTED;
}


int OPT_Optizelle::setStringParameter(const char *param, const char *value)
{
    OPT_PRINTERRORMSG("String parameters are not being set in optzelle");
    return OPT_OPERATION_NOT_IMPLEMENTED;
}


int OPT_Optizelle::setVariableType( const int index, const OPT_VARTYPE varType )
{
     //anyway we set var type in MIP_MINLPProb...
    int r = OPT_MyNLPSolver::setVariableType(index, varType);
    
    //optizelle is a continuous solver...
    if( varType == OPT_VT_INTEGER )
        return OPT_OPERATION_NOT_SUPPORTED;
    
    return r;
}


int OPT_Optizelle::allocateConstrAuxStructures(unsigned int m)
{
    
    OPT_realloc(evalEqAndNLIneqConstrs, m);
    //OPT_realloc(evalLinearIneqConstrs, m);
    OPT_realloc(linEqNLConstrs, m);
    OPT_realloc(linIneqConstrs, m);
    OPT_realloc(NLIneqConstrs, m);
    
    OPT_IFMEMERRORRETURN( !evalEqAndNLIneqConstrs || !linEqNLConstrs || !linIneqConstrs || !NLIneqConstrs );
    
    return 0;
}


int OPT_Optizelle::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_OPTIZELLE
{
    bool processConstrs = nmChg;
    int n, m;
    unsigned int mFree = 0;
    
    Optizelle::Unconstrained <double, Optizelle::Rm>::State::t *ustate = nullptr;
    Optizelle::EqualityConstrained <double, Optizelle::Rm, Optizelle::Rm>::State::t *eqstate = nullptr;
    Optizelle::InequalityConstrained <double, Optizelle::Rm, Optizelle::Rm>::State::t *ineqstate = nullptr;
    Optizelle::Constrained <double, Optizelle::Rm, Optizelle::Rm, Optizelle::Rm>::State::t *cstate = nullptr;
    
    
    
    getNumberOfVars(n);
    getNumberOfConstraints(m);
    
    
    if( nmChg )
    {
        int r = allocateConstrAuxStructures(m);
        OPT_IFERRORGOTOLABEL(r, retCode, r, termination);
    }
    
    
    if( genConstrChg )
        processConstrs = true;
    
    
    
    if( processConstrs )
    {
        const bool *nlConstr = prob.nlConstr;
        const double *lc = prob.lc, *uc = prob.uc;
        const double *lx = prob.lx, *ux = prob.ux;
        const OPT_SparseMatrix *QC = prob.QC;
        
        
        //OPT_setAllArray(m, evalLinearIneqConstrs, false);
        OPT_setAllArray(m, evalEqAndNLIneqConstrs, false);
        
        hasQuadMatrixInConstrs = false;
        mLinIneq = 0;
        //mLinEq = 0;
        mLinEqNL = 0;
        //mNLEq = 0;
        mNLIneq = 0;
        mFree = 0;
        
        
        for( decltype(m) i = 0; i < m; i++ )
        {
            
            if( lc[i] <= -OPT_INFINITY && uc[i] >= OPT_INFINITY )
            {
                mFree++;
            }
            else 
            {
                
                if( lc[i] == uc[i] )
                {
                    /*if( nlConstr[i] || QC[i].getNumberOfElements() > 0 )
                        mNLEq++;
                    else
                        mLinEq++;*/
                    
                    
                    linEqNLConstrs[mLinEqNL] = i;
                    evalEqAndNLIneqConstrs[i] = true;
                    mLinEqNL++;
                }
                else //inequality constraint
                {
                    const auto nQC = QC[i].getNumberOfElements();
                    
                    if( nlConstr[i] || nQC > 0 )
                    {
                        NLIneqConstrs[mNLIneq] = i;
                        mNLIneq++;
                        
                        linEqNLConstrs[mLinEqNL] = i;
                        mLinEqNL++;
                        
                        evalEqAndNLIneqConstrs[i] = true;
                        
                        if( nQC > 0 )
                            hasQuadMatrixInConstrs = true;
                    }
                    else
                    {
                        linIneqConstrs[mLinIneq] = i;
                        //evalLinearIneqConstrs[i] = true;
                        mLinIneq++;
                    }                    
                }
            }
            
        }
        
        #if OPT_DEBUG_MODE
            assert( mLinIneq + mLinEqNL + mFree == m );
        #endif
        
        OPT_realloc(NLIneqConstrs, mNLIneq);
        OPT_realloc(linEqNLConstrs, mLinEqNL);
        OPT_realloc(linIneqConstrs, mLinIneq);
        
        
        
        totalLinIneqConstrs = mLinIneq;
        
        //extra linear constraints to double bounded inequalities
        for(unsigned int i = 0; i < mLinIneq; i++)
        {
            const auto cind = linIneqConstrs[i];
            if( lc[cind] > -OPT_INFINITY && uc[cind] < OPT_INFINITY )
                totalLinIneqConstrs++; //we have to put two inequalities in this cases
        }
        
        //linear constraints to variable bounds
        for(unsigned int i = 0; i < n; i++)
        {
            if( lx[i] > -OPT_INFINITY )
                totalLinIneqConstrs++;
            if( ux[i] < OPT_INFINITY )
                totalLinIneqConstrs++;
        }
        
        //linear constraints to artificial variable bounds for nonlinear inequalities
        for(unsigned int i = 0; i < mNLIneq; i++)
        {
            const auto cind = NLIneqConstrs[i];
            if( lc[cind] > -OPT_INFINITY )
                totalLinIneqConstrs++;
            if( uc[cind] < OPT_INFINITY )
                totalLinIneqConstrs++;
        }
        
        totalConstrs = totalLinIneqConstrs + mLinEqNL;
    }
    
    
    
    {        
        double *constrValues = auxValues;
        
        
        if( n > 0 )
        {
            const auto *lx = prob.lx;
            const auto *ux = prob.ux;
            
            {
                auto xInitSize = xInit_vec.size();
                
                if( xInitSize < n)
                {
                    xInit_vec.resize(n + mNLIneq);
                    xInit_vec[0] = NAN;
                }
                else if( xInitSize < n + mNLIneq )
                {
                    xInit_vec.resize(n + mNLIneq);
                }
            }
            
            
            if( isnan( xInit_vec[0] ) )
            {
                
                for( decltype(n) i = 0; i < n; i++ )
                {
                    xInit_vec[i] = lx[i] <= 0.0 && 0.0 <= ux[i] ? 0.0 : ( lx[i] > -OPT_INFINITY ? lx[i] : ux[i]  ) ;
                }
            }
            
            //putting initial values in artificial variables to NL inequality constraints
            if( mNLIneq > 0 )
            {
                char k = 0;
                int r;
                bool *constrEval = (bool*) auxIndex;
                double *axInit_vec = xInit_vec.data();
                
                OPT_setAllArray(m, constrEval, false);
                
                for(unsigned int i = 0; i < mNLIneq; i++)
                    constrEval[ NLIneqConstrs[i] ] = true;
                
                do
                {
                    k++;
                    
                    r = prob.constraintsEval(threadNumber, true, constrEval, axInit_vec, constrValues );
                    
                    if(r != 0)
                    {
                        //we got error in constraint evaluation. So, we perturb the initial solution
                        for(unsigned int i = 0; i < n; i++)
                        {
                            if( lx[i] != ux[i] )
                            {
                                if( axInit_vec[i] + 0.5 <= ux[i] )
                                    axInit_vec[i] += 0.5;
                                else if( axInit_vec[i] - 0.3 >= lx[i] )
                                    axInit_vec[i] -= 0.3;
                                else //so, even lx[i] and ux[i] are not OPT_INFINITY
                                {
                                    const double step = (lx[i] + ux[i])/10.0;
                                    if( axInit_vec[i] + step <= ux[i] )
                                        axInit_vec[i] += step;
                                    else //( axInit_vec[i] - step >= lx[i] )
                                        axInit_vec[i] -= 0.9*step;
                                }
                                    
                            }
                        }
                    }
                    
                }while(r != 0 && k <= 10);
                
                
                for(unsigned int i = 0; i < mNLIneq; i++)
                {
                    xInit_vec[n+i] = constrValues[NLIneqConstrs[i]];
                }
                
                
            }
            
            
            /*if( zInit_vec.size() < n )
            {
                zInit_vec.resize(n);
                zInit_vec[0] = NAN;
            }
            
            if( isnan( zInit_vec[0] ) )
            {
                OPT_setAllArray( n, zInit_vec.data(), 0.0 );
            }*/
            
        }
        
        
        
        if( m > 0 )
        {
            if( lambdaInit_vec.size() < m )
            {
                lambdaInit_vec.resize(m);
                lambdaInit_vec[0] = NAN;
            }
            
            if( isnan( lambdaInit_vec[0] ) )
                OPT_setAllArray(m, lambdaInit_vec.data(), 0.0);
        }
        
        
        if( totalConstrs == 0 )
        {
            // Create a bundle of functions
            Optizelle::Unconstrained <double, Optizelle::Rm>::Functions::t fns;
            OPT_ObjToOptizelle *pMyObj;
            
            ustate = new (std::nothrow) Optizelle::Unconstrained <double, Optizelle::Rm>::State::t(xInit_vec);
            
            OPT_IFMEMERRORGOTOLABEL(!ustate, retCode, termination);
            
            
            pMyObj = new (std::nothrow) OPT_ObjToOptizelle(this, &prob);
            
            OPT_IFMEMERRORGOTOLABEL(!pMyObj, retCode, termination); //do not put ustate here, because pMyOBj cannot be free by us (optizelle takes owner ship n it)
            
            fns.f.reset(pMyObj);
            
            
            // Solve the optimization problem
            Optizelle::Unconstrained <double, Optizelle::Rm>::Algorithms::getMin (Optizelle::Messaging::stdout, fns, *ustate);
        }
        else
        {
            //so, we have to separate linear inequalities
                
            const double *lc = prob.lc, *uc = prob.uc;
            
            std::vector<double> lambdaLinIneq;
            std::vector<double> lambdaRemainder;
            
            lambdaLinIneq.resize(totalLinIneqConstrs);
            lambdaRemainder.resize(mLinEqNL);
            
            unsigned int k = mLinIneq;
            
            for(decltype(mLinIneq) i = 0; i < mLinIneq; i++)
            {
                const auto cind = linIneqConstrs[i];
                
                #if OPT_DEBUG_MODE
                    assert( !evalEqAndNLIneqConstrs[cind] );
                #endif
                
                if( lc[cind] <= -OPT_INFINITY  )
                {
                    lambdaLinIneq[i] = lambdaInit_vec[cind];
                }
                else
                {
                    lambdaLinIneq[i] = -lambdaInit_vec[cind];
                    if( uc[cind] < OPT_INFINITY )
                    { //double bounded constraint
                        lambdaLinIneq[k] = lambdaInit_vec[cind] ;
                        k++;
                    }
                }
            }
            
            
            
            for( decltype(mLinEqNL) i = 0; i < mLinEqNL; i++ )
            {
                const auto cind = linEqNLConstrs[i];
                lambdaRemainder[i] = lambdaInit_vec[cind];
                
                #if OPT_DEBUG_MODE
                    assert( evalEqAndNLIneqConstrs[cind] );
                #endif
            }
            
            
            if( mLinEqNL == 0 )
            {
                // Create a bundle of functions
                Optizelle::InequalityConstrained <double, Optizelle::Rm, Optizelle::Rm>::Functions::t fns;
                OPT_ObjToOptizelle *pMyObj;
                OPT_IneqToOptizelle *pMyIneq;
                
                
                ineqstate = new (std::nothrow) Optizelle::InequalityConstrained <double, Optizelle::Rm, Optizelle::Rm>::State::t(xInit_vec, lambdaLinIneq);
                OPT_IFMEMERRORGOTOLABEL(!ineqstate, retCode, termination);
                
                
                pMyObj = new (std::nothrow) OPT_ObjToOptizelle(this, &prob);
                pMyIneq = new (std::nothrow) OPT_IneqToOptizelle(this, &prob);
                
                //optizelle takes ownership about callback object pointers GRRRR! So, we cannot desalocate in our termination block
                
                if( !pMyObj || !pMyIneq )
                {
                    OPT_PRINTMEMERROR;
                    
                    if(pMyObj) delete pMyObj;
                    if(pMyIneq)delete pMyIneq;
                    
                    retCode = OPT_MEMORY_ERROR;
                    goto termination;
                }
                
                
                fns.f.reset(pMyObj);
                fns.h.reset(pMyIneq);
                
                // Solve the optimization problem
                Optizelle::InequalityConstrained <double, Optizelle::Rm, Optizelle::Rm>::Algorithms::getMin (Optizelle::Messaging::stdout, fns, *ineqstate);
                
            }
            else if( totalLinIneqConstrs == 0 )
            {
                // Create a bundle of functions
                Optizelle::EqualityConstrained <double, Optizelle::Rm, Optizelle::Rm>::Functions::t fns;
                OPT_ObjToOptizelle *pMyObj;
                OPT_EqToOptizelle *pMyEq;
                
                eqstate = new (std::nothrow) Optizelle::EqualityConstrained <double, Optizelle::Rm, Optizelle::Rm>::State::t(xInit_vec, lambdaRemainder);
                OPT_IFMEMERRORGOTOLABEL(!eqstate, retCode, termination);
                
                
                pMyObj = new (std::nothrow) OPT_ObjToOptizelle(this, &prob);
                pMyEq = new (std::nothrow) OPT_EqToOptizelle(this, &prob);
                
                //optizelle takes ownership about callback object pointers GRRRR! So, we cannot desalocate in our termination block
                
                if( !pMyObj || !pMyEq )
                {
                    OPT_PRINTMEMERROR;
                    
                    if(pMyObj) delete pMyObj;
                    if(pMyEq)  delete pMyEq;
                    
                    retCode = OPT_MEMORY_ERROR;
                    goto termination;
                }
                
                
                fns.f.reset(pMyObj);
                fns.g.reset(pMyEq);
                
                
                Optizelle::EqualityConstrained <double, Optizelle::Rm, Optizelle::Rm>::Algorithms::getMin(Optizelle::Messaging::stdout, fns, *eqstate);
            }
            else
            {
                // Create a bundle of functions
                Optizelle::Constrained <double, Optizelle::Rm, Optizelle::Rm, Optizelle::Rm>::Functions::t fns;
                
                OPT_ObjToOptizelle *pMyObj;
                OPT_EqToOptizelle *pMyEq;
                OPT_IneqToOptizelle *pMyIneq;
                
                
                //TODO: remove these asserts
                for(unsigned int i = 0; i < xInit_vec.size(); i++)
                    assert( !isnan(xInit_vec[i]) );
                
                //TODO: remove these asserts
                for(unsigned i = 0 ; i < lambdaRemainder.size(); i++)
                    assert( !isnan(lambdaRemainder[i]) );
                
                //TODO: remove these asserts
                for(unsigned i = 0 ; i < lambdaLinIneq.size(); i++)
                    assert( !isnan(lambdaLinIneq[i]) );
                
                cstate = new (std::nothrow) Optizelle::Constrained <double, Optizelle::Rm, Optizelle::Rm, Optizelle::Rm>::State::t (xInit_vec, lambdaRemainder, lambdaLinIneq);
                
                
                pMyObj = new (std::nothrow) OPT_ObjToOptizelle(this, &prob);
                pMyEq = new (std::nothrow) OPT_EqToOptizelle(this, &prob);
                pMyIneq = new (std::nothrow) OPT_IneqToOptizelle(this, &prob);
                
                //optizelle takes ownership about callback object pointers GRRRR! So, we cannot desalocate in our termination block
                
                if( !pMyObj || !pMyEq || !pMyIneq )
                {
                    OPT_PRINTMEMERROR;
                    
                    if(pMyObj) delete pMyObj;
                    if(pMyEq)  delete pMyEq;
                    if(pMyIneq)delete pMyIneq;
                    
                    retCode = OPT_MEMORY_ERROR;
                    goto termination;
                }
                
                
                
                //optizelle takes ownership about pointers: ARF!
                fns.f.reset(pMyObj );
                fns.g.reset(pMyEq );
                fns.h.reset(pMyIneq );
                
                Optizelle::Constrained <double, Optizelle::Rm, Optizelle::Rm, Optizelle::Rm>::Algorithms::getMin(Optizelle::Messaging::stdout, fns, *cstate);
                
                std::cout << "opt_stop: " << cstate->opt_stop << ": " << Optizelle::OptimizationStop::to_string(cstate->opt_stop) << "\n";
                
                
            }
        
        }
    }
    
    
    
termination:
    
    
    //liberar o state aqui (talvez)
    if(ustate)  delete ustate;
    if(ineqstate) delete ineqstate;
    if(eqstate) delete eqstate;
    if(cstate)  delete cstate;
    
    
    return retCode;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif
