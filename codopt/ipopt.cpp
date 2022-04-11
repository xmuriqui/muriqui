
#include <cmath>
#include <cstdlib>
#include <climits>
#include <iostream>
#include <new>


#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"

#if OPT_HAVE_IPOPT
    #include "IpIpoptData.hpp"

    #include "IpIpoptCalculatedQuantities.hpp"
    #include "IpTNLPAdapter.hpp"
    #include "IpOrigIpoptNLP.hpp"
#endif


#define OPT_INITIALIZE_IPOPT_INPUT_DATA 0 //initialize data to help debug

#define OPT_PRINT_IPOPT_INPUT_DATA 0 //printing the arrays to check with valgrind if there is something uninitialised




//using namespace std;
using namespace optsolvers;
using namespace newspm;
using namespace minlpproblem;




#if OPT_HAVE_IPOPT

using namespace Ipopt;

class optsolvers::OPT_MyProblemToIpopt : public TNLP
{
protected:
    
    bool newx;
    bool hasFreeConstrs;
    bool hasFreeNLConstrs;
    
public:
    
    bool storeSol, storeConstrs, storeDualSol;
    int m, mquad;
    int *quadIndex;
    
    int *sizeColsNzRowJac;
    int **colsNzRowJac;
    
    int *sizeColsNzLagH;
    int **colsNzRowLagH;
    
    int nNzRowsLagH;
    int *nzRowsLagH;
    
    OPT_Ipopt *nlp;
    MIP_MINLPProb *prob;
    
    
    
    
    
    OPT_MyProblemToIpopt();
    
    void initialize();
    
    //int allocate(const int m);
    
    void deallocate();
    
    
    virtual ~OPT_MyProblemToIpopt();
    
    
    /** Method to return some info about the nlp */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style);
    
    
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u);
    
    
    virtual bool get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda);
    
    
    /** Method to return the objective value */
    virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);
    
    
    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);
    
    
    /** Method to return the constraint residuals */
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);
    
    
    /** Method to return:
*   1) The structure of the jacobian (if "values" is NULL)
*   2) The values of the jacobian (if "values" is not NULL)
*/
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values);
    
    
    /** Method to return:
*   1) The structure of the hessian of the lagrangian (if "values" is NULL)
*   2) The values of the hessian of the lagrangian (if "values" is not NULL)
*/
    virtual bool eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values);
    
    
    virtual void finalize_solution(SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq);
    
    
    
    /** overload this method to return the constraint linearity.
    *  array has been allocated with length at least n. (default implementation
    *  just return false and does not fill the array).*/
    virtual bool get_constraints_linearity(Index m, LinearityType* const_types);
    
    
    /** Intermediate Callback method for the user.  Providing dummy
    *  default implementation.  For details see IntermediateCallBack
    *  in IpNLP.hpp. */
    virtual bool intermediate_callback(AlgorithmMode mode,
                                    Index iter, Number obj_value,
                                    Number inf_pr, Number inf_du,
                                    Number mu, Number d_norm,
                                    Number regularization_size,
                                    Number alpha_du, Number alpha_pr,
                                    Index ls_trials,
                                    const IpoptData* ip_data,
                                    IpoptCalculatedQuantities* ip_cq);
    
    
};




OPT_MyProblemToIpopt::OPT_MyProblemToIpopt()
{
    initialize();
}



OPT_MyProblemToIpopt::~OPT_MyProblemToIpopt()
{
    deallocate();
}


/*int OPT_MyProblemToIpopt::allocate(const int m)
{
    return 0;
}*/


void OPT_MyProblemToIpopt::deallocate()
{
    OPT_secFree(quadIndex);
    
    
    if(colsNzRowJac)
    {
        for(int i = 0; i < m; i++)
        {
            if(colsNzRowJac[i])
                free(colsNzRowJac[i]);
        }
        
        free(colsNzRowJac);
        colsNzRowJac = NULL;
    }
    
    OPT_secFree(sizeColsNzRowJac);
    
    if(colsNzRowLagH)
    {
        for(int i = 0; i < nNzRowsLagH; i++)
        {
            if(colsNzRowLagH[i])
                free(colsNzRowLagH[i]);
        }
        free(colsNzRowLagH);
        colsNzRowLagH = NULL;
    }
    
    OPT_secFree(sizeColsNzLagH);
    OPT_secFree(nzRowsLagH);
    
    mquad = 0;
    nNzRowsLagH = 0;
}


void OPT_MyProblemToIpopt::initialize()
{
    hasFreeConstrs = false;
    hasFreeNLConstrs = false;
    storeSol = true;
    storeConstrs = true;
    storeDualSol = true;
    
    
    mquad = 0;
    m = 0;
    quadIndex = NULL;
    
    sizeColsNzRowJac = NULL;
    colsNzRowJac = NULL;
    
    sizeColsNzLagH = NULL;
    colsNzRowLagH = NULL;
    
    nzRowsLagH = NULL;
}




/** Method to return some info about the nlp */
bool OPT_MyProblemToIpopt::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    //std::cout << "Entrei em get_nlp_info" << std::endl;
    const MIP_SparseMatrix *QC = prob->QC;

    bool *auxCEval = nlp->auxCEval;
    const bool *nlConstr = prob->nlConstr;
    const double *lc = prob->lc, *uc = prob->uc;
    //bool *flagCols = (bool*) nlp->auxIndex;
    int newm = 0, om = prob->m;

    deallocate(); //we put it here because derivative test does not free the memory


    newx = true;

    n = prob->n;
    m = om;
    this->m = prob->m;


    index_style = C_STYLE;

    hasFreeConstrs = false;
    hasFreeNLConstrs = false;

    for(int i = 0; i < om; i++)
    {
        if(lc[i] > -MIP_INFINITY || uc[i] < MIP_INFINITY)
        { 
            auxCEval[i] = true;
            newm++;
        }
        else //free constraint
        {
            //std::cout << "restricao " << i << " livre!\n";
            //OPT_getchar();
            hasFreeConstrs = true;
            auxCEval[i] = false;
            
            if( nlConstr[i] )
                hasFreeNLConstrs = true;
            
        }
    }

    if(hasFreeConstrs)
        m = newm;


    //calculating quadIndex
    {
        int k = 0;
        mquad = 0;
        
        
        for(int i = 0; i < om; i++)
        {
            if(auxCEval[i] && QC[i].getNumberOfElements() > 0)
                mquad++;
        }
        
        OPT_malloc(quadIndex, mquad); //quadIndex = (int*) malloc(mquad * sizeof(int));
        if(!quadIndex)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTMEMERROR;
            #endif
            
            return false;
        }
        
        
        for(int i = 0; i < om; i++)
        {
            if(auxCEval[i] && QC[i].getNumberOfElements() > 0)
            {
                quadIndex[k] = i;
                k++;
            }
        }
        
        #if OPT_DEBUG_MODE
            assert(k == mquad);
        #endif
    }

    //calculating number of nonzeros in jacobian
    {
        const OPT_SparseMatrix *singleJ;
        
        int r = OPT_getSingleJacPointerOrCalcJacIndices(*prob, mquad, hasFreeConstrs, auxCEval, singleJ, nnz_jac_g, sizeColsNzRowJac, colsNzRowJac);
        
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            return false;
        }
    }


    #if 0

    //calculating number of nonzeros in hessian
    {
        //testing if we have a hessian composed by more than one sparse matriz
        
        const OPT_SparseMatrix *singleH = OPT_getSingleLagHPointer(*prob, mquad, quadIndex);
        
        
        if(singleH)
        {
            nnz_h_lag = singleH->getNumberOfElements();
        }
        else
        {
            int ret, k = 0;
            std::unordered_set<int> setnzRowsLagH;
            
            
            
            OPT_fillNzRowsLagH(*prob, auxCEval, mquad, quadIndex, setnzRowsLagH);
            
            
            nNzRowsLagH = setnzRowsLagH.size();
            
            nzRowsLagH = (int *) malloc( nNzRowsLagH * sizeof(int) );
            if(!nzRowsLagH)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTMEMERROR;
                #endif
                
                return false;
            }
            
            
            {
                int i = 0;
                for(auto it: setnzRowsLagH)
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
                
                return false;
            }
            
            
            
            ret = OPT_fillColsNzRowHess(*prob, auxCEval, mquad, quadIndex, nNzRowsLagH , nzRowsLagH, colsNzRowLagH, sizeColsNzLagH);
            
            if(ret != 0)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(ret);
                #endif
                
                return false;
            }
            
            nnz_h_lag = 0;
            
            for(int k = 0; k < nNzRowsLagH; k++)
                nnz_h_lag += sizeColsNzLagH[k];
            
        }
        
        
    }

    #endif

    //calculating number of nonzeros in hessian
    {
        const OPT_SparseMatrix *singleH;
        
        int r = OPT_getSingleLagHPointerOrCalcLagHIndices(*prob, mquad, quadIndex, auxCEval, singleH, nnz_h_lag, nNzRowsLagH, nzRowsLagH, sizeColsNzLagH, colsNzRowLagH);
        
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            return false;
        }
        
        /*if(singleH)
        {
            std::cout << "singleH e nulo!\n";
        }
        else
        {
            std::cout << "singleH nao e nulo!\n";
            std::cout << "nNzRowsLagH: " << nNzRowsLagH << "\n";
            
            for(int i = 0; i < nNzRowsLagH; i++)
            {
                std::cout << "nzRowsLagH["<<i<<"]: " << nzRowsLagH[i] << "\n";
            }
            
        }
        
        OPT_getchar();*/
    }


    #if OPT_PRINT_IPOPT_INPUT_DATA //printing the arrays to check with valgrind if there is something uninitialised
    {
        std::cout << "n: " << n << "\t";
        std::cout << "m: " << m << "\t";
        std::cout << "nnz_jac_g: " << nnz_jac_g << "\t";
        std::cout << "nnz_h_lag: " << nnz_h_lag << "\t";
        std::cout << "index_style: " << index_style << "\t";
    }
    #endif

    //std::cout << "Sai de get_nlp_info" << std::endl;
    return true;
}


bool OPT_MyProblemToIpopt::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
    //std::cout << "Entrei em get_bounds_info" << std::endl;
    
    #if OPT_INITIALIZE_IPOPT_INPUT_DATA
    {
        OPT_setAllArray<double>(n, x_l, NAN);
        OPT_setAllArray<double>(n, x_u, NAN);
        
        OPT_setAllArray<double>(m, g_l, NAN);
        OPT_setAllArray<double>(m, g_u, NAN);
    }
    #endif
    

    prob->getVariableLowerBounds(x_l);
    prob->getVariableUpperBounds(x_u);
    

    if( hasFreeConstrs )
    {
        const int om = prob->m;
        const bool *auxCEval = nlp->auxCEval;
        const double *lc = prob->lc, *uc = prob->uc;
        
        int k = 0;
        
        for(int i = 0; i < om; i++)
        {
            if(auxCEval[i])
            {
                g_l[k] = lc[i];
                g_u[k] = uc[i];
                k++;
            }
            
            /*else
            {
                //actually, we do not evaluate this constraint in the pratice. (we always eval like 0)
                g_l[i] = -OPT_INFINITY;
                g_u[i] = OPT_INFINITY;
            }*/
        }
        
        #if OPT_DEBUG_MODE
            assert(k == m);
        #endif
    }
    else
    {
        //no redundant constraints...
        
        #if OPT_DEBUG_MODE
            assert( m == prob->m);
        #endif
        
        prob->getConstraintLowerBounds(g_l);
        prob->getConstraintUpperBounds(g_u);
    }


    #if OPT_PRINT_IPOPT_INPUT_DATA //printing the arrays to check with valgrind if there is something uninitialised
    {
        for(decltype(n) i = 0; i < n; i++)
        {
            std::cout << "x_l["<<i<<"]: " << x_l[i] << " ";
            std::cout << "x_u["<<i<<"]: " << x_u[i] << " \t";
            assert( !std::isnan(x_l[i]) );
            assert( !std::isinf(x_l[i]) );
            assert( !std::isnan(x_u[i]) );
            assert( !std::isinf(x_u[i]) );
        }
        
        for(decltype(m) j = 0; j < m; j++)
        {
            std::cout << "g_l["<<j<<"]: " << g_l[j] << " ";
            std::cout << "g_u["<<j<<"]: " << g_u[j] << " \t";
            assert( !std::isnan(g_l[j] ) );
            assert( !std::isinf(g_l[j] ) );
            assert( !std::isnan(g_u[j] ) );
            assert( !std::isinf(g_u[j] ) );
        }
    }
    #endif


    //std::cout << "Sai de get_bounds_info" << std::endl;
    return true;
}


bool OPT_MyProblemToIpopt::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda)
{
    //std::cout << "Entrei em get_starting_point" << std::endl;
    const double *xInit = nlp->xInit;
    const double *zInit = nlp->zInit;
    const double *lambdaInit = nlp->lambdaInit;
    
    #if OPT_INITIALIZE_IPOPT_INPUT_DATA
    {
        if( init_x )
        {
            OPT_setAllArray<double>(n, x, NAN);
        }
        
        if( init_z )
        {
            OPT_setAllArray<double>(n, z_L, NAN);
            OPT_setAllArray<double>(n, z_U, NAN);
        }
        
        if( init_lambda )
        {
            OPT_setAllArray<double>(m, lambda, NAN);
        }
    }
    #endif
    
    
    if( init_x )
    {
        const double *lx = prob->lx;
        const double *ux = prob->ux;
        
        if( std::isnan( xInit[0] ) )
        {
            //double l, u;
            
            #pragma ivdep
            #pragma GCC ivdep
            for(int i = 0; i < n; i++)
            {
                const double l = lx[i];
                const double u = ux[i];
                
                x[i] = l > -MIP_INFINITY ? l :  (u < MIP_INFINITY ? u : 0) ;
            }
            
        }
        else
        {
            OPT_copyArray(n, xInit, x);
            
            #if 0
            {
                //just by safe, we check if we have some value out of bounds or nan
                for(int i = 0; i < n; i++)
                {
                    const double l = lx[i];
                    const double u = ux[i];
                    
                    if( std::isnan( x[i] ) || x[i] < l || x[i] > u )
                    {
                        x[i] = l > -MIP_INFINITY ? l : (u < MIP_INFINITY ? u : 0);
                    }
                }
            }
            #endif
        }
        
        /*for(int i = 0; i < n; i++)
            printf("xInit[%d]: %0.12f\n", i, x[i]);
        
        OPT_getchar();*/
    }
    
    if( init_z )
    {
        if( std::isnan(zInit[0]) )
        {
            OPT_setAllArray( n, z_L, 0.0 );
            OPT_setAllArray( n, z_U, 0.0);
        }
        else
        {
            OPT_copyArray( n, zInit, z_L );
            OPT_copyArray( n, &zInit[n], z_U );
        }
        
    }
    
    
    if(init_lambda)
    {
        if( std::isnan(lambdaInit[0]) )
        {
            OPT_setAllArray( m, lambda, 0.0);
        }
        else
        {
            if(hasFreeConstrs)
            {
                int k = 0;
                const int om = prob->m;
                const bool *auxCEval = nlp->auxCEval;
                
                
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
            
            //OPT_copyArray( m, lambdaInit, lambda );
        }
    }
    
    
    #if OPT_PRINT_IPOPT_INPUT_DATA //printing the arrays to check with valgrind if there is something uninitialised
    {
        if( init_x )
        {
            for(decltype(n) i = 0; i < n; i++)
            {
                std::cout << "x["<<i<<"]: " << x[i] << " \t";
                assert( !std::isnan(x[i]) );
                assert( !std::isinf(x[i]) );
            }
        }
        
        if( init_z )
        {
            for(decltype(n) i = 0; i < n; i++)
            {
                std::cout << "z_L["<<i<<"]: " << z_L[i] << " ";
                std::cout << "z_U["<<i<<"]: " << z_U[i] << " \t";
                assert( !std::isnan(z_L[i]) );
                assert( !std::isinf(z_L[i]) );
                assert( !std::isnan(z_U[i]) );
                assert( !std::isinf(z_U[i]) );
            }
        }
        
        if( init_lambda )
        {
            for(decltype(m) j = 0; j < m; j++)
            {
                std::cout << "lambda["<<j<<"]: " << lambda[j] << " \t";
                assert( !std::isnan(lambda[j]) );
                assert( !std::isinf(lambda[j]) );
            }
        }
    }
    #endif
    
    
    //std::cout << "Sai de get_starting_point" << std::endl;
    return true;
}


/** Method to return the objective value */
bool OPT_MyProblemToIpopt::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    //std::cout << "Entrei em eval_f" << std::endl;
    const unsigned int thnumber = nlp->threadNumber;
    int r;
    
    #if OPT_INITIALIZE_IPOPT_INPUT_DATA
    {
        obj_value = NAN;
    }
    #endif
    
    
    if( new_x )
        newx = true;
    
    r = prob->objEval(thnumber, newx, x, obj_value, nlp->in_nl_obj_factor);
    
    //std:: cout << "obj_value: " << obj_value << std::endl;
    
    if( prob->hasNlObj )
        newx = false;
    
    if( r != 0 )
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return false;
    }
    
    
    #if OPT_PRINT_IPOPT_INPUT_DATA //printing the arrays to check with valgrind if there is something uninitialised
    {
        std::cout << "obj_value: " << obj_value << " \t";
        assert( !std::isnan(obj_value) );
        assert( !std::isinf(obj_value) );
    }
    #endif
    
    //std::cout << "Sai de eval_f" << std::endl;
    return true;
}



/** Method to return the gradient of the objective */
bool OPT_MyProblemToIpopt::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    //std::cout << "Entrei em eval_grad_f" << std::endl;
    const unsigned int thnumber = nlp->threadNumber;
    int r;
    
    #if OPT_INITIALIZE_IPOPT_INPUT_DATA
    {
        OPT_setAllArray<double>(n, grad_f, NAN);
    }
    #endif
    
    
    if( new_x )
        newx = true;
    
    r = prob->objGradEval(thnumber, newx, x, grad_f, nlp->in_nl_obj_factor);
    
    if( prob->hasNlObj )
        newx = false;
    
    //for(int i = 0; i < n; i++)
        //std::cout << "grad_f["<<i<<"]: " << grad_f[i] << std::endl; 
    
    if( r != 0 )
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return false;
    }
    
    
    #if OPT_PRINT_IPOPT_INPUT_DATA //printing the arrays to check with valgrind if there is something uninitialised
    {
        for(decltype(n) i = 0; i < n; i++)
        {
            std::cout << "grad_f["<<i<<"]: " << grad_f[i] << " \t";
            assert( !std::isnan(grad_f[i]) );
            assert( !std::isinf(grad_f[i]) );
        }
    }
    #endif
    
    
    //std::cout << "Sai de eval_grad_f" << std::endl;
    return true;
}



/** Method to return the constraint residuals */
bool OPT_MyProblemToIpopt::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    //std::cout << "Entrei em eval_g" << std::endl;
    const bool *auxCEval = nlp->auxCEval;
    const int thnumber = nlp->threadNumber;
    double *myg = nlp->auxValues2;
    double *pg  = hasFreeConstrs ? myg : g;
    int r;
    
    #if OPT_INITIALIZE_IPOPT_INPUT_DATA
    {
        OPT_setAllArray<double>(m, g, NAN);
    }
    #endif
    
    
    if( new_x )
        newx = true;
    
    
    r = prob->constraintsEval( thnumber, newx, auxCEval, x, pg );
    
    if( prob->hasNlConstrs )
        newx = false;
    
    if( hasFreeConstrs )
    {
        int k = 0;
        const int om = prob->m;
        
        for(int i = 0; i < om; i++)
        {
            if(auxCEval[i])
            {
                g[k] = pg[i];
                k++;
            }
        }
        #if OPT_DEBUG_MODE
            assert(k == m);
            assert(m < om);
        #endif
    }
    
    
    /*if( m < prob->m )
    {
        for(int i = 0; i < m; i++)
            g[i] = pg[ auxNewConstrIndex[i] ]; //we only can do it because auxIndex[i] is always greater or equal to i
    }*/
    
    /*for(int i = 0; i < OPT_max(prob->n, prob->m) ; i++)
    {
        if(i < prob->n)
            std::cout << "x["<<i<<"]: " << x[i] << " ";
        if(i < prob->m )
            std::cout << "\tg["<<i<<"]: " << g[i] << " auxCEval: " << auxCEval[i];
        std::cout << "\n";
    }*/
    
    if( r != 0 )
    {
        #if OPT_PRINT_CALLBACK_ERROR_MSG
            std::cerr << OPT_PREPRINT << "Callback function error " << r << std::endl;
        #endif
        return false;
    }
    
    
    #if OPT_PRINT_IPOPT_INPUT_DATA //printing the arrays to check with valgrind if there is something uninitialised
    {
        for(decltype(m) i = 0; i < m; i++)
        {
            std::cout << "g["<<i<<"]: " << g[i] << " \t";
            assert( !std::isnan(g[i]) );
            assert( !std::isinf(g[i]) );
        }
    }
    #endif
    
    
    //std::cout << "Sai de eval_g" << std::endl;
    return true;
}


/** Method to return:
*   1) The structure of the jacobian (if "values" is NULL)
*   2) The values of the jacobian (if "values" is not NULL)
*/
bool OPT_MyProblemToIpopt::eval_jac_g(Index n, const Number* x, bool new_x, Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values)
{
    //std::cout << "Entrei em eval_jac_g" << std::endl;
    
    const int om = prob->m;
    //const OPT_SparseMatrix &A = prob->A;
    //const OPT_SparseMatrix *QC = prob->QC;
    const bool *auxCEval = nlp->auxCEval;
    
    const OPT_SparseMatrix *singleJ = OPT_getSingleJacobianPointer(*prob, mquad);
    
    
    #if OPT_INITIALIZE_IPOPT_INPUT_DATA
    {
        if(values)
            OPT_setAllArray<double>(nele_jac, values, NAN);
        if(iRow)
            OPT_setAllArray<Index>(nele_jac, iRow, m+1);
        if(jCol)
            OPT_setAllArray<Index>(nele_jac, jCol, n+1);
    }
    #endif
    
    
    
    if( new_x )
        newx = true;
    
    newx = true;
    
    if( values )
    {
        const bool *auxCEval = nlp->auxCEval;
        //const bool *nlConstr = prob->nlConstr;
        
        int nzs;
        OPT_SparseMatrix &J = prob->J;
        
        
        if( prob->hasNlConstrs )
        {
            const int thnumber = nlp->threadNumber;
            
            const int r = prob->nlJacobianEval( thnumber, newx, auxCEval, x, J );
            
            newx = false;
            
            if( r != 0 )
            {
                #if OPT_PRINT_CALLBACK_ERROR_MSG && OPT_PRINT_CALLBACK_ERROR_MSG
                    std::cerr << OPT_PREPRINT "Callback function error " << r << "\n";
                #endif
                return false;
            }
        }
        
        
        if(singleJ && !hasFreeConstrs)
        {
            nzs = singleJ->getValues(values);
        }
        else
        {
            nzs = OPT_getValuesFromCompleteJacobian( *prob, J, sizeColsNzRowJac, colsNzRowJac, auxCEval, m, x, nlp->auxValues, values);
        }
        
        
        #if OPT_DEBUG_MODE
            assert(nzs == nele_jac);
        #endif
        
        
        #if 0
        {
            //checking derivatives:
            const int NP = 4, N = 2, S = NP-2;
            const double bigM = sqrt(2), LAMBDA = 0.0001;
            double A[NP][N] = { {1, 0}, {0, 0}, {1, 1}, {0, 1} };
            
            int mynzj = 0;
            
            //checking contsraint 1
            for(int i = 0; i < NP; i++)
            {
                for(int j = 0; j < S; j++)
                {
                    double vsquare = LAMBDA;
                    
                    for(int k = 0; k < N; k++)
                        vsquare += pow(A[i][k] - x[j*N + k], 2);
                    
                    double v1 = 0.5*pow(vsquare, -0.5);
                    
                    for(int k = 0; k < N; k++)
                    {
                        double dvk = v1 * 2 * (A[i][k] - x[j*N+k]);
                        
                        //printf("i: %d j: %d k: %d dvk: %f values[%d]: %f\n", i, j, k, dvk, mynzj, values[mynzj] );
                        assert( abs(dvk - values[mynzj]) < 1e-8  );
                        
                        mynzj++;
                    }
                    
                    //printf("dvk: %f values[%d]: %f\n", 1.0, mynzj, values[mynzj]);
                    assert( values[mynzj] == 1.0 ); // d variable
                    mynzj++;
                    
                    //printf("dvk: %f values[%d]: %f\n", -bigM, mynzj, values[mynzj]);
                    assert( abs(values[mynzj] - (-bigM)) < 1e-6 ); //integer variable multiplying big-M 
                    mynzj++;
                }
                
            }
            
            //checking constraint 2
            for(int i = 0; i < S; i++)
            {
                for(int j = i+1; j < S; j++)
                {
                    double vsquare = LAMBDA;
                    
                    for(int k = 0; k < N; k++)
                        vsquare += pow(x[i*N + k] - x[j*N + k], 2);
                    
                    double v1 = 0.5*pow(vsquare, -0.5);
                    
                    for(int k = 0; k < N; k++)
                    {
                        double dvk = v1 * 2 * (x[i*N + k] - x[j*N + k]) * (-1);
                        
                        assert( abs(dvk - values[mynzj]) < 1e-8  );
                        mynzj++;
                    }
                    
                    for(int k = 0; k < N; k++)
                    {
                        double dvk = v1 * 2 * (x[i*N + k] - x[j*N + k]);
                        
                        assert( abs(dvk - values[mynzj]) < 1e-8  );
                        mynzj++;
                    }
                    
                    assert( values[mynzj] == 1.0 ); // d variable
                    mynzj++;
                    
                    assert( abs(values[mynzj] - (-bigM)) < 1e-6 ); //integer variable multiplying big-M 
                    mynzj++;
                }
            }
            
        }
        #endif
        
        
    }
    else
    {
        //const OPT_SparseMatrix &J = prob->J;
        const OPT_SparseMatrix *singleJ = OPT_getSingleJacobianPointer(*prob, mquad);
        
        int nzs;
        
        
        if(singleJ && !hasFreeConstrs)
        {
            nzs = singleJ->getStructure(iRow, jCol);
        }
        else
        {
            nzs = 0;
            int k = 0;
            
            for(int i = 0; i < om; i++)
            {
                if( !auxCEval[i] )
                    continue;
                
                const int ncols = sizeColsNzRowJac[i];
                
                OPT_setAllArray(ncols, &iRow[nzs], k);
                OPT_copyArray(ncols, colsNzRowJac[i], &jCol[nzs]);
                
                nzs += ncols;
                k++;
            }
            
            #if OPT_DEBUG_MODE
                assert(k == m);
            #endif
        }
        
        
        /*for(int w = 0; w < nzs; w++)
            std::cout << "iRow["<<w<<"]: " << iRow[w] << " jCol["<<w<<"]: " << jCol[w] << "\n";*/
        
        
        #if OPT_DEBUG_MODE
            assert(nzs == nele_jac);
        #endif
        
        
        /*for(int i = 0; i < nele_jac; i++)
        {
            std::cout << "iRow["<<i<<"]: " << iRow[i];
            std::cout << " jCol["<<i<<"]: " << jCol[i] << std::endl;
        }*/
        
    }
    
    
    #if OPT_PRINT_IPOPT_INPUT_DATA //printing the arrays to check with valgrind if there is something uninitialised
    {
        if( values )
        {
            for(decltype(nele_jac) i = 0; i < nele_jac; i++)
            {
                std::cout << "values["<<i<<"]: " << values[i] << " \t";
                assert( !std::isnan(values[i]) );
                assert( !std::isinf(values[i]) );
            }
        }
        else
        {
            for(decltype(nele_jac) i = 0; i < nele_jac; i++)
            {
                std::cout << "iRow["<<i<<"]: " << iRow[i] << " ";
                std::cout << "jCol["<<i<<"]: " << jCol[i] << " \t";
                assert(iRow[i] >= 0 && iRow[i] < m);
                assert(jCol[i] >= 0 && jCol[i] < n);
            }
        }
    }
    #endif
    
    //std::cout << "Sai de eval_jac_g" << std::endl;
    return true;
}


/** Method to return:
*   1) The structure of the hessian of the lagrangian (if "values" is NULL)
*   2) The values of the hessian of the lagrangian (if "values" is not NULL)
*/
bool OPT_MyProblemToIpopt::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda, bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values)
{
    //std::cout << "Entrei em eval_h" << std::endl;
    
    #if OPT_INITIALIZE_IPOPT_INPUT_DATA
    {
        if(values)
            OPT_setAllArray<double>(nele_hess, values, NAN);
        if(iRow)
            OPT_setAllArray<Index>(nele_hess, iRow, n+1);
        if(jCol)
            OPT_setAllArray<Index>(nele_hess, jCol, n+1);
    }
    #endif
    
    
    if( new_x )
        newx = true;
    
    if( values )
    {
        const int thnumber = nlp->threadNumber;
        
        OPT_SparseMatrix &lagH = prob->lagH;
        const bool *auxCEval = nlp->auxCEval;
        int nzs;
        double *auxValues = nlp->auxValues;
        double *auxValues2 = nlp->auxValues2;
        
        
        int r = OPT_evalCompleteLagrangianHessian( thnumber, newx, x, *prob, lagH, mquad, quadIndex, obj_factor, nlp->in_nl_obj_factor, lambda, auxCEval, m, nNzRowsLagH, nzRowsLagH, sizeColsNzLagH, colsNzRowLagH, auxValues, auxValues2, nzs, values);
        
        if(r != 0)
        {
            #if OPT_DEBUG_MODE && OPT_PRINT_CALLBACK_ERROR_MSG
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return false;
        }
        
        #if OPT_DEBUG_MODE
            assert(nzs == nele_hess);
        #endif
        
        
        /*for(int i = 0; i < nele_hess; i++)
            std::cout << "values["<<i<<"]: " << values[i] << std::endl;*/
        
    }
    else
    {
        const OPT_SparseMatrix *singleH = OPT_getSingleLagHPointer(*prob, mquad, quadIndex);
        
        int nzs;
        
        if(singleH)
        {
            nzs = singleH->getStructure(iRow, jCol);
        }
        else
        {
            nzs = 0;
            
            for(int i = 0; i < nNzRowsLagH; i++)
            {
                const int ncols = sizeColsNzLagH[i];
                
                OPT_setAllArray(ncols, &iRow[nzs], nzRowsLagH[i]);
                
                OPT_copyArray(ncols, colsNzRowLagH[i], &jCol[nzs] );
                
                nzs += ncols;
            }
        }
        
        #if OPT_DEBUG_MODE
            assert(nzs == nele_hess);
        #endif
        
        
        
        /*for(int i = 0; i < nele_hess; i++)
        {
            std::cout << "iRow["<<i<<"]: " << iRow[i];
            std::cout << " jCol["<<i<<"]: " << jCol[i] << std::endl;
        }*/
        
    }
    
    
    #if OPT_PRINT_IPOPT_INPUT_DATA //printing the arrays to check with valgrind if there is something uninitialised
    {
        if( values )
        {
            for(decltype(nele_hess) i = 0; i < nele_hess; i++)
            {
                std::cout << "values["<<i<<"]: " << values[i] << " \t";
                assert( !std::isnan(values[i]) );
                assert( !std::isinf(values[i]) );
            }
        }
        else
        {
            for(decltype(nele_hess) i = 0; i < nele_hess; i++)
            {
                std::cout << "iRow["<<i<<"]: " << iRow[i] << " ";
                std::cout << "jCol["<<i<<"]: " << jCol[i] << " \t";
                assert(iRow[i] >= 0 && iRow[i] < n);
                assert(jCol[i] >= 0 && jCol[i] < n);
            }
        }
    }
    #endif
    
    //std::cout << "Sai de eval_h" << std::endl;
    return true;
}



void OPT_MyProblemToIpopt::finalize_solution( SolverReturn status, Index n, const Number* x, const Number* z_L, const Number* z_U, Index m, const Number* g, const Number* lambda, Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
{
    //std::cout << "Entrei em finalize_solution" << std::endl;
    
    
    nlp->objValue = obj_value;
    
    if(ip_data)
        nlp->numberOfIterations = ip_data->iter_count();
    
    
    if( prob->objFactor < 0 )
    {
        nlp->objValue = -nlp->objValue;
    }
    
    
    if(storeSol)
        OPT_copyArray(n, x, nlp->sol);
    
    if(storeDualSol)
    {
        OPT_copyArray(n, z_L, nlp->dualSolV);
        OPT_copyArray(n, z_U, &nlp->dualSolV[n]);
        
        if(hasFreeConstrs)
        {
            const int om = prob->m;
            const bool *auxCEval = nlp->auxCEval;
            double *dualSolC = nlp->dualSolC;
            
            int k = 0;
            
            for(int i = 0; i < om; i++)
            {
                if(auxCEval[i])
                {
                    dualSolC[i] = lambda[k];
                    k++;
                }
                else
                    dualSolC[i] = 0.0;
            }
            
            #if OPT_DEBUG_MODE
                assert(k == m);
            #endif
        }
        else
        {
            OPT_copyArray(m, lambda, nlp->dualSolC);
        }
    }
    
    
    if( storeConstrs )
    {
        if(hasFreeConstrs)
        {
            const int om = prob->m;
            const bool *auxCEval = nlp->auxCEval;
            double *constr = nlp->constr;
            
            int k = 0;
            
            for(int i = 0; i < om; i++)
            {
                if(auxCEval[i])
                {
                    constr[i] = g[k];
                    k++;
                }
                else
                    constr[i] = 0.0;
            }
            
            
            #if OPT_DEBUG_MODE
                assert(k == m);
            #endif
            
            //if( hasFreeConstrs && storeConstrs )
            {
                bool *ceval = (bool *) nlp->auxIndex;
                double *auxConstr = nlp->auxValues;
                
                for(int i = 0; i < om; i++)
                    ceval[i] = !auxCEval[i];
                
                int r = prob->constraintsEval( nlp->threadNumber, true, ceval, x, auxConstr );
                
                if(r != 0)
                    OPT_PRINTERRORMSG("warning: error evaluation on redundant constraints");
                
                for(int i = 0; i < om; i++)
                {
                    if( ceval[i] )
                        constr[i] = auxConstr[i];
                }
                
            }
        }
        else
        {
            OPT_copyArray(m, g, nlp->constr);
        }
    }
    
    
    
    deallocate();
    
    #if 0
    //checking if solution is feasible
    {
        const double absFeasTol = 1.0e-5;
        const double relFeasTol = 1.0e-5;
        
        bool &feasSol = nlp->feasSol;
        
        const double *lc = prob->lc;
        const double *uc = prob->uc;
        
        
        feasSol = true;
        
        for(int i = 0; i < mo; i++)
        {
            //note, we use not operator (!) here because constr[i] can be nan. In this case, the test inside parenthesis will fail and we consider that is not a feasible sol...
            
            
            if(lc[i] > -OPT_INFINITY)
            {
                const double ltol = absFeasTol + OPT_abs( lc[i]*relFeasTol );
                
                //note, we use not operator (!) here because constr[i] can be nan. In this case, the test inside parenthesis will fail and we consider that is not a feasible sol...
                if( !(lc[i] <= constr[i] + ltol) )
                {
                    std::cout << "Restricao " << i << "nao e viavel lc["<<i<<"]: " << lc[i] << " ltol: " << ltol << "g["<<i<<"]: " << constr[i] << "\n";
                    
                    feasSol = false;
                    break;
                }
            }
            
            if(uc[i] < OPT_INFINITY)
            {
                const double utol = absFeasTol + OPT_abs( uc[i]*relFeasTol );
                
                //note, we use not operator (!) here because constr[i] can be nan. In this case, the test inside parenthesis will fail and we consider that is not a feasible sol...
                if( !(constr[i] - utol <= uc[i]) )
                {
                    std::cout << "Restricao " << i << " nao e viavel g["<<i<<"]: " << constr[i] << " uc["<<i<<"]: " << uc[i] << " utol: " << utol << "\n";
                    
                    feasSol = false;
                    break;
                }
            }
            
            
            
            /*const double feasTol = OPT_abs( constr[i]*relFeasTol )  + absFeasTol;
            
            //note, we use not operator (!) here because constr[i] can be nan. In thsi case, the test inside parenthesis will fail and we consider that is not a feasible sol...
            if ( !( lc[i] <= constr[i] + feasTol  &&  constr[i] - feasTol <= uc[i] ) )
            {
                feasSol = false;
                break;
            } */
            
        }
        
        std::cout << "feasSol: " << feasSol << "\n";
    }
    #endif
    
    
    
    //std::cout << "OPT_MyProblemToIpopt::Sai de finalize_solution" << std::endl;
    //OPT_getchar();
}



bool OPT_MyProblemToIpopt::get_constraints_linearity(Index m, LinearityType* const_types)
{
    const int om = prob->m;
    const bool *nlConstr = prob->nlConstr;
    const double *lc = prob->lx, *uc = prob->uc;
    const MIP_SparseMatrix *QC = prob->QC;
    
    int k = 0;
    
    for(int i = 0; i < om; i++)
    {
        if(lc[i] <= -MIP_INFINITY && uc[i] >= MIP_INFINITY)
        { //free constraint
            continue; //const_types[i] = TNLP::LINEAR;
        }
        
        
        const_types[k] = nlConstr[i] || QC[i].getNumberOfElements() > 0 ? TNLP::NON_LINEAR : TNLP::LINEAR;
        k++;
    }
    
    #if OPT_DEBUG_MODE
        assert(k == m);
    #endif
    
    
    #if OPT_PRINT_IPOPT_INPUT_DATA
    {
        for(decltype(m) i = 0; i < m; i++)
        {
            std::cout << "i: " << i << " const_types: " << const_types[i] << "\n";
        }
    }
    #endif
        
    
    return true;
}



bool OPT_MyProblemToIpopt::intermediate_callback(AlgorithmMode mode, Index iter, Number obj_value, Number inf_pr, Number inf_du, Number mu, Number d_norm, Number regularization_size, Number alpha_du, Number alpha_pr, Index ls_trials, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
{
    
    
    if( nlp->ipoptInterCallback )
    {
        return nlp->ipoptInterCallback->intermediate_callback(mode, iter, obj_value, inf_pr, inf_du, mu, d_norm, regularization_size, alpha_du, alpha_pr, ls_trials, ip_data, ip_cq);
    }
    
    return true;
    
    #if 0
    //getting the current solution
    {
        const int n = prob->getNumberOfVars();
        
        
        Ipopt::TNLPAdapter* tnlp_adapter = NULL;
        if( ip_cq != NULL )
        { 
            Ipopt::OrigIpoptNLP* orignlp = NULL;
            orignlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
            
            
            if( orignlp != NULL )
                tnlp_adapter = dynamic_cast<TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
        }
        
        //If ipopt is in the restoration phase, tnlp_adapter will be NULL
        
        if( tnlp_adapter )
        {
            double x[n];
            
            tnlp_adapter->ResortX(*ip_data->curr()->x(), x);
            
            for(int i = 0; i < n; i++)
                std::cout << "x" << i << ": " << x[i] << "  ";
            std::cout << "\n";
        }
        
    }
    
    
    std::cout << "ipopt iter: " << iter << " mode: " << mode << " obj_value: " << obj_value << "  ";
    
    if(iter %2 == 0)
        std::cout << "\n";
    
    return true;
    #endif
}




#endif










OPT_Ipopt::OPT_Ipopt():OPT_MyNLPSolver()
{
    initialize();
}




OPT_Ipopt::~OPT_Ipopt()
{
    //desallocateMemory();
    deallocateSolverEnv();
}



// __methods from Solver __


void OPT_Ipopt::deallocateSolverEnv()
{
#if OPT_HAVE_IPOPT
    app = NULL; //app is a smart pointer. So, I hope that is enough to free the object pointed by it
    ipoptInterCallback = NULL;
#endif
    
    OPT_MyNLPSolver::deallocateSolverEnv();
}



bool OPT_Ipopt::getMinusLambdaOnLagran()
{
    return false;
}



/*int OPT_Ipopt::allocateVarStructures(const int n)
{
    double *auxd;
    
    auxd = (double *) realloc( xInit, n * sizeof(double) );
    if( !auxd )
        return OPT_MEMORY_ERROR;
    
    xInit = auxd;
    
    
    auxd = (double *) realloc( zInit, 2* n * sizeof(double) );
    if( !auxd )
        return OPT_MEMORY_ERROR;
    
    zInit = auxd;
    
    if(n > 0)
    {
        xInit[0] = NAN;
        zInit[0] = NAN;
    }
    
    
    return OPT_MyNLPSolver::allocateVarStructures(n);
}



int OPT_Ipopt::allocateConstrStructures(const int m)
{
    double *auxd;
    
    auxd = (double *) realloc( lambdaInit, m* sizeof(double) );
    if( !auxd )
        return OPT_MEMORY_ERROR;
    
    lambdaInit = auxd;
    
    if( m > 0 )
        lambdaInit[0] = NAN;
    
    
    return OPT_MyNLPSolver::allocateConstrStructures( m);
} */


void OPT_Ipopt::deallocateMemory()
{
    
    OPT_MyNLPSolver::deallocateMemory();
    //OPT_secFree(quadIndex);
    //OPT_secFree(nqcons);
}



void OPT_Ipopt::initialize()
{
    OPT_MyNLPSolver::initialize();
    
    #if OPT_HAVE_IPOPT
        origSolverRetCode = INT_MAX;
    #endif
    
    jac_c_const = false;
    jac_d_const = false;
    hessian_const = false;
    
    mquad = -1;
    numberOfIterations = 0;
    
    //quadIndex = NULL;
    //nqcons = NULL;
    
    enforceReoptimization = false;
}


int OPT_Ipopt::getNumberOfIterations(long unsigned int &niter)
{
    niter = numberOfIterations;
    return 0;
}


OPT_LISTSOLVERS OPT_Ipopt::getSolverCode()
{
    return OPT_IPOPT;
}



int OPT_Ipopt::getVariableType( const int index, OPT_VARTYPE &varType )
{
    varType = optsolvers::OPT_VT_CONTINUOUS;
    
    return 0;
}


int OPT_Ipopt::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
#if OPT_HAVE_IPOPT
{
    //int r;
    
    app = IpoptApplicationFactory();
    
    mynlp = NULL;
    
    const int r = app->Initialize();
    if( r != Solve_Succeeded )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    //app->Options()->SetStringValue( "honor_original_bounds", "yes");
    
    setOutputLevel(1);
    
    ipoptInterCallback = NULL;
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::setMaxCPUTime(const double time)
#if OPT_HAVE_IPOPT
{
    app->Options()->SetNumericValue( "max_cpu_time" , time );
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::setNumberOfThreads(const int nthreads)
#if OPT_HAVE_IPOPT
{
    return 0;
    //return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::setOutputLevel(const int level)
#if OPT_HAVE_IPOPT
{
    app->Options()->SetIntegerValue( "print_level" , level );
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::setRelativeDualTol( const double tol )
#if OPT_HAVE_IPOPT
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_IPOPT
{
    app->Options()->SetNumericValue( "tol" , tol );
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::setRelativePrimalTol( const double tol )
#if OPT_HAVE_IPOPT
{
    return OPT_OPERATION_NOT_SUPPORTED;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_IPOPT
{
    const bool r = app->Options()->SetNumericValue( param, value );
    
    if( !r )
    {
        //#if OPT_DEBUG_MODE
            //OPT_PRINTERRORNUMBER(r);
        //#endif
        
        printDblParamErrorMsg(!r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_IPOPT
{
    const bool r = app->Options()->SetIntegerValue( param, value );
    
    if( !r )
    {
        //#if OPT_DEBUG_MODE
            //OPT_PRINTERRORNUMBER(r);
        //#endif
        
        printIntParamErrorMsg(!r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_IPOPT
{
    const bool r = app->Options()->SetStringValue( param, value );
    
    if( !r )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        printStrParamErrorMsg(r, param, value);
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::setVariableType( const int index, const OPT_VARTYPE varType )
#if OPT_HAVE_IPOPT
{
    //anyway we set var type in MIP_MINLPProb...
    int r = OPT_MyNLPSolver::setVariableType(index, varType);
    
    
    //ipopt is a continuous solver...
    if( varType == OPT_VT_INTEGER )
        return OPT_OPERATION_NOT_SUPPORTED;
    
    return r;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol)
#if OPT_HAVE_IPOPT
{
    //const bool calcNqconsIndqcons = genQuadConstrChg || nmChg;
    bool reoptimize = false;
    int n, m;
    int r;
    ApplicationReturnStatus status;
    
    
    
    getNumberOfVars(n);
    getNumberOfConstraints(m);
    
    
    if(resetSol)
    {
        //origSolverRetCode = INT_MAX;
        this->resetSol();
        numberOfIterations = 0;
    }
    else
    {
        feasSol = false;
    }
    
    
    {
        const bool hasNlExp = (prob.hasNlObj || prob.hasNlConstrs);
        bool flagEqConstrsLinear = true, flagIneqConstrsLinear = true;
        bool change = false;
        
        
        
        for(int i = 0; i < m; i++)
        {
            if( prob.lc[i] == prob.uc[i] )
            {
                if( prob.nlConstr[i] || prob.QC[i].getNumberOfElements() > 0 )
                {
                    flagEqConstrsLinear = false;
                    if( flagIneqConstrsLinear == false )
                        break;
                }
            }
            else
            {
                if( prob.nlConstr[i] || prob.QC[i].getNumberOfElements() > 0 )
                {
                    flagIneqConstrsLinear = false;
                    if( flagEqConstrsLinear == false )
                        break;
                }
            }
        }
        
        if( flagEqConstrsLinear != jac_c_const )
        {
            app->Options()->SetStringValue("jac_c_constant", flagEqConstrsLinear ? "yes" : "no" ); //if flagEqConstrsLinear is true, all equalities constraints are linear...
            
            jac_c_const = flagEqConstrsLinear;
            change = true;
        }
        
        
        if( flagIneqConstrsLinear != jac_d_const )
        {
            app->Options()->SetStringValue("jac_d_constant", flagIneqConstrsLinear ? "yes" : "no" ); //if flagIneqConstrsLinear is true, all inequalities constraints are linear...
            
            jac_d_const = flagIneqConstrsLinear;
            change = true;
        }
        
        
        if( hessian_const != hasNlExp )
        {
            app->Options()->SetStringValue( "hessian_constant", hasNlExp ? "no" : "yes" );
            
            
            hessian_const = (prob.hasNlObj || prob.hasNlConstrs);
            change = true;
        }
        
        
        
        
        
        
        if( change )
        {
            //I believe app->Initialize process the options, so we need initialize again. However, maybe newer version of ipopt do not require this reinitialization (I do not know... in ipopt manual there is no information...)
            r = app->Initialize();
            if( r != Solve_Succeeded )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                return OPT_SOLVER_ERROR;
            } 
            
            //OPT_getchar();
        }
    }
    
    
    
    
    
    if( IsNull(mynlp) )
    {
        OPT_MyProblemToIpopt *p = new (std::nothrow) OPT_MyProblemToIpopt();
        OPT_IFMEMERRORGOTOLABEL( !p, retCode, termination );
        
        mynlp = p;
        mynlp->prob = &prob;
        mynlp->nlp = this;
    }
    
    mynlp->storeSol = storeSol;
    mynlp->storeConstrs = storeConstrs;
    mynlp->storeDualSol = storeDualSol;
    
    
    //prob.print();
    
    
    if(reoptimize || enforceReoptimization)
        status = app->ReOptimizeTNLP(mynlp);
    else
        status = app->OptimizeTNLP(mynlp);
    
    origSolverRetCode = status;
    
    //std::cout << OPT_PREPRINT << "ipopt status: " << status << std::endl;
    
    switch (status)
    {
        case Solve_Succeeded:
        case Solved_To_Acceptable_Level:
    
            retCode = OPT_OPTIMAL_SOLUTION;
            feasSol = true;
            break;
            
        case Restoration_Failed: 
            
            retCode = OPT_NO_FEASIBLE_SOLUTION_FOUND;
            break;
        
        case Infeasible_Problem_Detected:
            
            retCode = OPT_INFEASIBLE_PROBLEM;
            feasSol = false;
            break;
            
        case Feasible_Point_Found:
            
            retCode = OPT_FEASIBLE_SOLUTION;
            feasSol = true;
            break;
            
        case Maximum_Iterations_Exceeded :
            #if OPT_PRINT_MAX_ITER_WARNING
                if( numberOfWarningsByIterLimit < maxNumberOfWarningsByIterLimit )
                {
                    std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Ipopt solving!\n";
                    //std::cerr << OPT_PREPRINT "feasSol: " << feasSol << "\n";
                    numberOfWarningsByIterLimit++;
                    
                    if( numberOfWarningsByIterLimit == maxNumberOfWarningsByIterLimit )
                        std::cerr << OPT_PREPRINT "Warning: Maximum number of warnings by maximum iteration achieved! Stopping these warnings.\n";
                }
            #endif
            retCode = OPT_MAX_ITERATIONS;
            break;
            
        case Maximum_CpuTime_Exceeded:
            
            //std::cerr << OPT_PREPRINT "feasSol: " << feasSol << "\n";
            retCode = OPT_MAX_TIME;
            break;
            
        case Insufficient_Memory :
            
            retCode = OPT_MEMORY_ERROR;
            break;
        
        case Invalid_Number_Detected:
            
            retCode = OPT_CALLBACK_FUNCTION_ERROR;
            break;
            
        case User_Requested_Stop:
            retCode = OPT_STOP_REQUIRED_BY_USER;
            break;
            
        case Invalid_Option :
            
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORMSG("Invalid Option as Ipopt return code");
            #endif
            
            retCode = OPT_BAD_INPUT; 
            break;
            
        //case Invalid_Problem_Definition:
        default:
            
            retCode = OPT_UNDEFINED_ERROR;
            break;
    
    }
    
    
    if( retCode != OPT_OPTIMAL_SOLUTION && retCode != OPT_FEASIBLE_SOLUTION && retCode != OPT_INFEASIBLE_PROBLEM && retCode != OPT_NO_FEASIBLE_SOLUTION_FOUND )// && storeConstrs)
    {
        //sometimes, ipopt does not return the correct value of objective, so i am not sure if the same occours with the constraints... so, to be sure, I will evaluate constraints
        
        double *pcontrs = storeConstrs ? constr : auxValues;
        
        
        int r = prob.isFeasibleToConstraints(threadNumber, sol, true, NULL, in_absolute_feas_tol, in_relative_feas_tol, feasSol, pcontrs);
        
        if(r != 0)
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
        
        
        //feasSol = prob.isConstrValuesFeasible(in_absolute_feas_tol, in_relative_feas_tol, constr);
        
        
        
        if( feasSol )
        {
            int r = prob.objEval(threadNumber, !prob.hasNlConstrs, sol, objValue, in_nl_obj_factor );
            
            if(r != 0)
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                feasSol = false;
            }
        }
        
        
        /*for(int i = 0; i < prob.n; i++)
            printf("w[%d]: %0.12f\n", i, sol[i]);
        
        for(int i = 0; i < prob.m; i++)
            printf("g[%d]: %0.12f\n", i, constr[i]);*/
        
        //std::cout << "Testei solucao viavel! feasSol: " << feasSol << "\n";
        
        //OPT_getchar();
    }
    
    
    
termination:
    
    enforceReoptimization = false;
    
    return retCode;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Ipopt::warmUp()
#if OPT_HAVE_IPOPT
{
    enforceReoptimization = true;
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



#if OPT_HAVE_IPOPT
void OPT_Ipopt::setIntermediateCallbackPointer( OPT_IpoptIntermediateCallback *intermediateCallback)
{
    //OPT_getchar();
    
    ipoptInterCallback = intermediateCallback;
}
#endif

// __ methods from NLPSolver __

















