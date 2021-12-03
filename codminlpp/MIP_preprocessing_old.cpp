/*
* Implementation of class to do some preprocessings to reduce variable bounds on a MIP_minlpProblem
* 
* References: 
* 		Laurence Wolsey, Integer Programming, Lohn Wiley & Sons, 1998 
* 		My own head
* 
* 
* Author: Wendel Melo
* Date: 15-June-2015
* 
*/ 

#include <cmath>
#include <climits>
#include <iostream>
#include <iomanip>
#include <new>
#include "MIP_minlpProblem.hpp"
#include "MIP_tools.hpp"



using namespace std;

using namespace minlpproblem;






MIP_Preprocessing_old::MIP_Preprocessing_old(const MIP_MINLPProb *prob)
{
    initialize(prob);
}


MIP_Preprocessing_old::~MIP_Preprocessing_old()
{
    desallocateMemory();
}


int MIP_Preprocessing_old::allocateMemory( const int nvars, const int nconstrs )
{
    MIP_malloc(auxFlags, nconstrs); //auxFlags = (char *) malloc( nconstrs * sizeof(char) );
    MIP_malloc(auxValues, nvars); //auxValues = (double *) malloc( nvars * sizeof(double) );
    
    if( !auxFlags || !auxValues )
    {
        #if MIP_DEBUG_MODE
            MIP_PRINTMEMERROR;
        #endif
        return MIP_MEMORY_ERROR;
    }
    
    return 0;
}


void MIP_Preprocessing_old::desallocateMemory()
{
    MIP_secFree( auxFlags );
    MIP_secFree( auxValues );
}


void MIP_Preprocessing_old::initialize(const MIP_MINLPProb *prob)
{
    auxFlags = NULL;
    auxValues = NULL;
    this->prob = prob;
    
    resetParameters();
}


void MIP_Preprocessing_old::resetParameters()
{
    in_abs_feas_tol = 1.0e-6;
    in_rel_feas_tol = 1.0e-8;
    
    in_abs_bound_tol = 1.0e-4; //if we reduce more, we can have numerical problems
    in_rel_bound_tol = 1.0e-6; //if we reduce more, we can have numerical problems 
}


//we try update variable bound using rhsMinusSum/coef. Dependig of coef's signal, we try update lower or upper bound
inline bool MIP_Preprocessing_old::__tryUpdateBound( const double rhsMinusSum, const double coef, const int varType,  double &l, double &u)
{
    const double boundTol = in_abs_bound_tol;
    const double boundFactorTol = in_rel_bound_tol;
    
    bool updtBound = false;
    double newBound;
    
    
    if( coef > 0.0 )
    {
        newBound = rhsMinusSum/coef; //we have a new upper bound
        
        if( MIP_isIntegerType(varType) )
            newBound = floor(newBound + MIP_abs(newBound*boundFactorTol) + boundTol);
        
        if( newBound < u )
        {
            u = newBound;
            updtBound = true;
        }
    }
    else if( coef < 0.0 )
    {
        newBound = rhsMinusSum/coef; //we have a new lower bound
        
        if( MIP_isIntegerType(varType) )
            newBound = ceil( newBound - MIP_abs(newBound*boundFactorTol) - boundTol);
        
        if( newBound > l )
        {
            l = newBound;
            updtBound = true;
        }
    }
    
    return updtBound;
}


//if leqConstr is true, we assume constraint is <= b. Otherwise, we asusme is >= b. That method return true if some bound changes 
bool MIP_Preprocessing_old::__preprocLinConstr( const unsigned int nzsRow, const int *rowCols, const double *rowValues, double b, const bool leqConstr, const int *varTypes, double *lx, double *ux, bool &infeasible, bool &redundant )
{
    const double feasTol = in_abs_feas_tol;
    const double feasFactorTol = in_rel_feas_tol;
    
    //const MIP_SparseElement *colAux = row.columns;
    
    unsigned int ninf = 0;
    unsigned int indinf;
    
    bool updtBounds = false;
    double coef;
    double sum, sump = 0.0, summ = 0.0; //we separet sum of positive and negative terms, We just do it because some day in the past, a professor said that reduces numerical errors.
    
    infeasible = false;
    redundant = false;
    
    
    if( !leqConstr )
        b = -b;
    
    
    for(unsigned int j = 0; j < nzsRow; j++)
    {
        const unsigned int col = rowCols[j]; //colAux[j].getColumn();
        
        coef = rowValues[j]; //colAux[j].getValue();
        if( !leqConstr )
            coef = -coef;
        
        
        //printf("col: %u coef: %lf lx: %lf ux: %lf b: %lf\n", col, coef, lx[col], ux[col], b);
        
        
        if( coef > 0.0 )
        {
            if( lx[col] > -MIP_INFINITY  )
            {
                sump += coef*lx[col] ;//ok, there is no guarantee that term is really positive, but we sum here anyway... we just want reduce numerical errors, not be perfect ;) (I do not want perform a n if here to test is this term is positive or not...)
            }
            else
            {
                if(ninf > 0)
                    goto termination; //that is the second infinity. We cannot use that constraint to preprocess..
                    
                ninf = 1;
                indinf = j;
            }
        }
        else if( coef < 0.0 )
        {
            if( ux[col] < MIP_INFINITY )
            {
                summ += coef*ux[col]; //ok, there is no guaratee that term is really negative, but we sum here anyway... we just want reduce numerical errors, not be perfect ;) (I do not want perform a n if here to test is this term is positive or not...)
            }
            else
            {
                if(ninf > 0)
                    goto termination; //that is the second infinity. We cannot use that constraint to preprocess..
                    
                ninf = 1;
                indinf = j;
            }
        }
        
        
    }
    
    sum = sump + summ;
    
    
    if( ninf == 0 && sum > b + MIP_abs(feasFactorTol *b) + feasTol )
    {
        //that constraint is infeasible
        //printf( "sum: %0.16f b: %0.16f sump: %0.10f summ: %0.10f\n", sum, b, sump, summ );
        
        infeasible = true;
        goto termination;
    }
    
    
    if( ninf == 1 ) //ninf == 1
    {
        //we have only 1 variable making sum infinity.
        const unsigned int col = rowCols[indinf]; //colAux[indinf].getColumn();
        coef = rowValues[indinf]; //colAux[indinf].getValue();
        if( !leqConstr )
            coef = -coef;
        
        
        const bool aux = __tryUpdateBound(b - sum, coef, varTypes[col], lx[col], ux[col]);
        
        if( aux )
            updtBounds = true;
        
    }
    else if( ninf == 0)
    {
        //running the variables updating the bounds
        double constrp = 0.0, constrm = 0.0;
        
        
        for(unsigned int j = 0; j < nzsRow; j++)
        {
            const unsigned int col = rowCols[j]; //colAux[j].getColumn();
            coef = rowValues[j]; //colAux[j].getValue();
            
            if( !leqConstr )
                coef = -coef;
            
            if( coef > 0.0 )
            {
                //desconsiderating the variable in the sum
                const double sumw = sum -coef*lx[col];
                
                const bool aux = __tryUpdateBound(b-sumw, coef, varTypes[col], lx[col], ux[col]);
                
                if(aux)
                    updtBounds = true;
                
                
                constrp += coef * ux[col]; 
            }
            else if( coef < 0.0 )
            {
                //desconsiderating the variable in the sum
                const double sumw = sum - coef*ux[col];
                
                const bool aux = __tryUpdateBound(b-sumw, coef, varTypes[col], lx[col], ux[col] );
                
                if(aux)
                    updtBounds = true;
                
                constrm += coef*lx[col];
            }
        }
        
        if( constrp + constrm <= b )
        {
            redundant = true;
        }
        
    }
    
    
    
termination:
    
    return updtBounds;
}


static unsigned int oldGIter = 0;

int MIP_Preprocessing_old::preprocess(const bool preprocQuadConstrs, const bool preprocObj, const double zu, double* lx, double* ux, bool& updtVarBounds, bool &updtConstrBounds, const double* inlc, const double* inuc, double* outlc, double* outuc )
{
    const unsigned int MAX_ITERS = 10;
    
    const int n = prob->n;
    const int m = prob->m;
    
    const double boundTol = in_abs_bound_tol;
    const double boundFactorTol = in_rel_bound_tol;
    
    
    const bool *nlContsr = prob->nlConstr;
    const int *xtype = prob->xtype;
    const double *lc;
    const double *uc;
    bool updt, updtIter, infeas, redun;
    const MIP_SparseMatrix &A = prob->A;
    const MIP_SparseMatrix *QC = prob->QC;
    
    unsigned int iter = 0;
    
    updtVarBounds = false;
    updtConstrBounds = false;
    
    
    oldGIter++;
    
    
    if( outlc )
    {
        lc = outlc;
        MIP_copyArray(m, inlc ? inlc : prob->lc,  outlc);
    }
    else
    {
        lc = inlc ? inlc : prob->lc;
    }
    
    
    if( outuc )
    {
        uc = outuc;
        MIP_copyArray(m, inuc ? inuc : prob->uc, outuc);
    }
    else
    {
        uc = inuc ? inuc : prob->uc;
    }
    
    
    /*for(int i = 0; i < n; i++)
    {
        if( MIP_isIntegerType(xtype[i]) )
        {
            lx[i] = ceil( lx[i] );
            ux[i] = floor( ux[i] );
        }
    }*/
    
    
    MIP_setAllArray<char>(m, auxFlags, 3); //3 is the vale 011 in binary. We are using the last bit to flag if we must preprocess constraint having lb > -MIP_INFINITY and the bif before to flag if we must preprocess constraint having ub < MIP_INFINITY
    
    
    do
    {
        iter++;
        updtIter = false;
        
        for(int i = 0; i < m; i++)
        {
            //cout << "i: " << i << endl;
            
            if( nlContsr[i] )
                continue;
            
            const bool linConstr = QC[i].getNumberOfElements() == 0;
            
            if( !linConstr && !preprocQuadConstrs )
                continue;
            
            
            
            
            const unsigned int ainzs = A.getNumberOfElementsAtRow(i);
            int* const aiCols = A.getRowColsPointer(i);
            double* const aiValues = A.getRowValuesPointer(i);
            
            
            //second part of if test if last bit is one
            if( lc[i] > -MIP_INFINITY && (auxFlags[i] & (char)1)  )
            {
                if( linConstr )
                {
                    updt =  __preprocLinConstr( ainzs, aiCols, aiValues, lc[i], false, xtype, lx, ux, infeas, redun ); // __preprocLinConstr( A[i], lc[i], false, xtype, lx, ux, infeas, redun );
                }
                else
                {
                    int r = __preprocQuadConstr( QC[i], ainzs, aiCols, aiValues, lc[i], false, xtype, lx, ux, updt, redun );
                    
                    infeas = r == MIP_INFEASIBILITY;
                }
                
                /*if( oldGIter > 1300 )
                    printf("\tOLD global iter: %u  iter: %d lc updt: %d infeas: %d redun: %d\n", oldGIter, iter, (int) updt, (int) infeas, (int) redun); */
                
                if( updt )
                    updtVarBounds = updtIter = true;
                
                if( infeas )
                {
                    #if MIP_PRINT_PREPROC_CUT_MSG
                        std::cout << MIP_PREPRINT "old preproc " << iter << ": Infeasibility detected on constraint " << i << "! :D\n";
                    #endif
                    //MIP_getchar(); 
                    return MIP_INFEASIBILITY;
                }
                
                if( redun )
                {
                    auxFlags[i] -= 1;
                    #if MIP_DEBUG_MODE
                        assert( auxFlags[i] >= 0 );
                    #endif
                    updtConstrBounds = true;
                    
                    if( outlc )
                    {
                        outlc[i] = -MIP_INFINITY; //note in this case lc point to outlc. So, this constraint will not be preprocessed again
                    }
                    
                    /*std::cout << "preproc " << iter << ": Detectei redundancia na restricao " << i << "! :D" << std::endl;
                    MIP_getchar();*/
                }
            }
            
            //last part of this if test if second bit is one
            
            //std::cout << "i: " << i;
            //std::cout << " uc["<< i << "]: " << uc[i];
            //std::cout << " auxFlags[" << i << "]: " << (int) auxFlags[i] << std::endl;
            
            if( uc[i] < MIP_INFINITY && (auxFlags[i] & (char)2 ) )
            {
                
                if( linConstr )
                {
                    updt = __preprocLinConstr( ainzs, aiCols, aiValues, uc[i], true, xtype, lx, ux, infeas, redun );
                }
                else
                {
                    int r = __preprocQuadConstr( QC[i], ainzs, aiCols, aiValues, uc[i], true, xtype, lx, ux, updt, redun);
                    
                    infeas = r == MIP_INFEASIBILITY;
                }
                
                /*if( oldGIter > 1300 )
                    printf("\tOLD global iter: %u  iter: %d uc updt: %d infeas: %d redun: %d\n", oldGIter, iter, (int) updt, (int) infeas, (int) redun); */
                
                
                if( updt )
                    updtVarBounds = updtIter = true;
                
                if( infeas )
                {
                    #if MIP_PRINT_PREPROC_CUT_MSG
                        std::cout << "old preproc " << iter << ": Detectei inviabilidade na restricao " << i << "! :D\n";
                    #endif
                    //MIP_getchar();
                    return MIP_INFEASIBILITY;
                }
                
                if( redun )
                {
                    auxFlags[i] -= 2;
                    
                    #if MIP_DEBUG_MODE
                        assert( auxFlags[i] >= 0 );
                    #endif
                    updtConstrBounds = true;
                    
                    if( outuc )
                    {
                        outuc[i] = MIP_INFINITY; //note in this case lc point to outlc. So, this constraint will not be preprocessed again
                    }
                    
                    //std::cout << "preproc " << iter << ": Detectei redundancia na restricao " << i << "! :D" << std::endl;
                    //MIP_getchar();
                }
            }
            
        }
        
        
        if( preprocObj && !prob->hasNlObj && zu < MIP_INFINITY )
        {
            const int r = preprocessObjF( zu, lx, ux, updt );
            
            if( updt )
                updtVarBounds = updtIter = true;
            
            if( r == MIP_INFEASIBILITY )
            {
                #if MIP_PRINT_PREPROC_CUT_MSG
                    std::cout << "old preproc " << iter << ": Detectei inviabilidade na fc obj! :D\n";
                #endif
                //MIP_getchar();
                return MIP_INFEASIBILITY;
            }
        }
        
        
        if( updtIter )
        {
            
            /*std::cout << std::fixed << std::setprecision(6) << "prepeoc " << iter << ": Atualizei limites na iteracao " << iter << std::endl ;
            
            for(int j = 0; j < prob->n; j++)
                std::cout << "lx[" << j << "]: " << lx[j] << " ux[" <<  j << "]: " << ux[j] << std::endl;
            
            MIP_getchar(); */
            
            
            for(int k = 0; k < n; k++)
            {
                const double lxk = lx[k];
                const double uxk = ux[k];
                
                if( lxk > uxk )
                {
                    
                    if( lxk - uxk > MIP_abs(uxk*boundFactorTol) + boundTol )
                    {
                        #if MIP_PRINT_PREPROC_CUT_MSG
                        std::cout << "old prepoc " << iter << ": Detectei inviabilidade por limite de variaveis na iteracao " << iter << "\n";
                        #endif
                        //MIP_getchar();
                        return MIP_INFEASIBILITY;
                    }
                    else
                    {
                        //we consider this bounds like equal...
                        
                        const double avg = (lxk + uxk)*0.5; //we do not know what is the most correct, lx or ux. So, we calculate an average between both (althoug it can introduce more numerical errors... :( )
                        
                        lx[k] = avg;
                        ux[k] = avg;
                    }
                    
                }
            }
        }
        
        
        
    }while( updtIter && iter <= MAX_ITERS );
    
    printf("old iter: %d\n", iter);
    
    return 0;
}



//if Q is diagonal, calculate the minimum value that xQx can assume
inline int MIP_Preprocessing_old::__calcQminValue(const MIP_SparseMatrix &M, const double factor, const double *lx, const double *ux, double *coefs, double &constValue)
{
    const int n = M.getNumberOfColumns();
    const int m = M.getNumberOfRows(); //M should be square, but ok...
    
    double sump = 0.0, summ = 0.0;
    
    
    constValue = 0.0;
    MIP_setAllArray(n, coefs, 0.0);
    
    //if( M.getNumberOfElements() > n )
        //return MIP_NOT_APPLICABLE; // by now, we do not preprocess general quadratics, and of course, at least one element is not in the diagonal...
    
    
    for( int i = 0; i < m; i++ )
    {
        const bool ifix = lx[i] == ux[i];
        //MIP_SparseRow &row = M[i];
        //const MIP_SparseElement * colAux = M.rows[i].columns;
        
        const unsigned int nel = M.getNumberOfElementsAtRow(i);  //row.getNumberOfElements();
        int* const rcols = M.getRowColsPointer(i);
        double* const rvalues = M.getRowValuesPointer(i);
        
        for( unsigned int j = 0; j < nel; j++ )
        {
            const int col = rcols[j]; // colAux[j].getColumn(); //row[j].getColumn();
            double v = rvalues[j]*factor;  // colAux[j].getValue()*factor; //row[j].getValue()*factor;
            
            
            if( col == i )
            { //diagonal poisition
                
                if( v > 0.0 )
                {//we have a term in form v/2 x_i ^ 2
                    
                    if( lx[i] <= 0.0 && 0.0 <= ux[i] )
                    { //minimum value is zero
                    }
                    else
                    {
                        if( MIP_abs(lx[i]) < MIP_abs(ux[i]) )
                            sump += 0.5*v*(lx[i]*lx[i]);
                        else
                            sump += 0.5*v*(ux[i]*ux[i]);
                    }
                }
                else
                {
                    //that is a very bizarre case, because that case is nonconvex... anyway, we treat it. The minimum value is in the extremity
                    
                    if( MIP_abs(lx[i]) < MIP_abs(ux[i]) )
                        summ += 0.5*v*(ux[i]*ux[i]);
                    else
                        summ += 0.5*v*(lx[i]*lx[i]);
                }
                
            }
            else
            { 
                
                if(  lx[col] != ux[col] && !ifix )
                { //neither i nor col are fixed
                    
                    double vmin;
                    
                    if( v > 0 )
                    {
                        int p, q;
                        
                        if( lx[i] <= lx[col] )
                            p = i, q = col;
                        else
                            p = col, q = i;
                        
                        const double lp = lx[p], up = ux[p];
                        const double lq = lx[q], uq = ux[q];
                        //lp <= lq
                        
                        if( lp >= 0 )
                        {
                            #if MIP_DEBUG_MODE
                            //that is has no sense... lp is greater than zero and the lower between lp and lq. So, both are greater than zero
                            if( lp <= -MIP_INFINITY || lq <= -MIP_INFINITY )
                            {
                                assert(false);
                                return MIP_NOT_APPLICABLE;
                            }
                            #endif
                            
                            
                            vmin = lp*lq;
                        }
                        else //lp < 0
                        {
                            if( uq >= 0 )
                            {
                                if( lp <= -MIP_INFINITY || uq >= MIP_INFINITY )
                                    return MIP_NOT_APPLICABLE;
                                
                                vmin = lp*uq;
                            }
                            else //uq < 0
                            {
                                //since uq is lower than zero, lq is lower than zero also
                                if( up >= 0 ) 
                                {
                                    if( up >= MIP_INFINITY || lq <= -MIP_INFINITY )
                                        return MIP_NOT_APPLICABLE;
                                    
                                    vmin = up*lq;
                                }
                                else //up < 0
                                {
                                    if( up >= MIP_INFINITY || uq >= MIP_INFINITY )
                                        return MIP_NOT_APPLICABLE;
                                    
                                    vmin = up*uq;
                                }
                            }
                        }
                        
                        sump += v*vmin ;
                    }
                    else // v < 0
                    {
                        int p, q;
                        
                        if( ux[i] <= ux[col] )
                            p = i, q = col;
                        else
                            p = col, q = i;
                        
                        
                        const double lp = lx[p], up = ux[p];
                        const double lq = lx[q], uq = ux[q];
                        //up <= uq
                        
                        
                        if( up >= 0 )
                        {
                            if( up >= MIP_INFINITY || uq >= MIP_INFINITY )
                                return MIP_NOT_APPLICABLE;
                            
                            vmin = up*uq;
                        }
                        else //up < 0
                        {
                            if( lq <= 0 )
                            {
                                if(lp <= -MIP_INFINITY || lq <= -MIP_INFINITY)
                                    return MIP_NOT_APPLICABLE;
                                    
                                vmin = lp*lq;
                            }
                            else //lq > 0
                            {
                                if(up >= MIP_INFINITY || lq <= -MIP_INFINITY)
                                    return MIP_NOT_APPLICABLE;
                                
                                vmin = up*lq;
                            }
                        }
                        
                        summ += v*vmin;
                    }
                    
                }
                else
                {//or i, or col is fixed (maybe both)
                    unsigned int other;
                    
                    if( ifix )
                    {
                        v = v*lx[i];  //i is fixed
                        other = col;
                    }
                    else
                    {
                        v = v*lx[col]; // so, col is fixed
                        other = i;
                    }
                    
                    
                    /*if( lx[other] == ux[other] )
                    {
                        //both variables are fixed...
                        
                        if( v >= 0.0 )
                            sump += v*lx[other]; //note, we do not dive by 2
                        else
                            summ += v*ux[other]; //note, we do not dive by 2
                    }
                    else */
                    {
                        coefs[ other ] = v;
                    }
                    
                }
            }
            
        }
        
        
        
        #if 0
        if( nel > 1 )
        {
            return MIP_NOT_APPLICABLE;
        }
        else if( nel == 1 )
        {
            if( (int) row[0].getColumn() != i )
                return MIP_NOT_APPLICABLE;
            
            //we have a diagonal element...
            const double v = row[0].getValue() * factor;
            
            if( v > 0.0 )
            {
                //we have a term in form v/2 x_i ^ 2
                
                if( lx[i] <= 0.0 && 0.0 <= ux[i] )
                { //minimum value is zero
                }
                else
                {
                    if( MIP_abs(lx[i]) < MIP_abs(ux[i]) )
                        value += v*lx[i]*lx[i];
                    else
                        value += v*ux[i]*ux[i];
                }
            }
            else
            {
                //that is very bizarr case, because that case is nonconvex... anyway, we treat it. The minimum value is in the extremity
                
                if( MIP_abs(lx[i]) < MIP_abs(ux[i]) )
                    value += v*ux[i]*ux[i];
                else
                    value += v*lx[i]*lx[i];
            }
            
        }
        #endif
        
    }
    
    constValue = 0.5*constValue; //some day you remove this line
    
    constValue += sump + summ;
    
    
    return 0;
}




int MIP_Preprocessing_old::preprocessObjF(const double zu, double *lx, double *ux, bool &updtBounds )
{
    const bool quadMatrix = prob->Q.getNumberOfElements() > 0;
    const int n = prob->n;
    const double objFactor = prob->objFactor;
    const int *xtype = prob->xtype;
    const double *c = prob->c;
    const MIP_SparseMatrix &Q = prob->Q;
    
    bool infeas, redundant;
    int code;
    //int ninf = 0;
    //int indinf;
    double quad, b;
    //double sum, sump = 0.0, summ = 0.0; //we separate sum of positive and negative terms, We just do it because some day in the past, a professor said that reduces numerical errors.
    
    
    updtBounds = false;
    
    if( prob->hasNlObj )
    {
        code = 0;
        goto termination;
    }
    
    
    if( !quadMatrix )
    {
        MIP_setAllArray(n, auxValues, 0.0);
        quad = 0.0;
    }
    else
    {
        //we try do a preprocess only if we have terms in diagonal
        
        const int r = __calcQminValue(Q, objFactor, lx, ux, auxValues, quad );
        
        if(r != 0)
        {
            code = 0;
            goto termination;
        }
    }
    
    b = zu - prob->d*objFactor - quad;
    
    if( prob->hasLinCoefObj() )
    {
        for(int i = 0; i < n; i++)
            auxValues[i] += objFactor*c[i];
        
    }
    
    updtBounds = __preprocFullLinConstr( n, auxValues, b, xtype, lx, ux, infeas, redundant );
    
    if( infeas )
    {
        code = MIP_INFEASIBILITY;
        goto termination;
    }
    
    
    #if 0
    
    if( prob->hasLinCoefObj() || quadMatrix )
    {
        
        for(int i = 0; i < n; i++)
        {
            const double coef = objFactor*c[i] + auxValues[i];
            
            if( coef > 0.0 )
            {
                if( lx[i] > -MIP_INFINITY )
                {
                    sump += coef*lx[i]; //ok, there is no guaratee that term is really positive, but we sum here anyway... we just want reduce numerical errors, not be perfect ;) (I do not want perform an if here to test if this term is positive or not...)
                }
                else
                {
                    if(ninf > 0)
                    {
                        code = 0;
                        goto termination; //that is the second infinity. We cannot use that constraint to preprocess..
                    }
                    
                    ninf = 1;
                    indinf = i;
                }
            }
            else
            {
                if( ux[i] < MIP_INFINITY )
                {
                    summ += coef*ux[i]; //ok, there is no guaratee that term is really negative, but we sum here anyway... we just want reduce numerical errors, not be perfect ;) (I do not want perform a n if here to test is this term is positive or not...)
                }
                else
                {
                    if(ninf > 0)
                    {
                        code = 0;
                        goto termination; //that is the second infinity. We cannot use that constraint to preprocess..
                    }
                    
                    ninf = 1;
                    indinf = i;
                }
            }
        }
        
        
        sum = sump + summ;
        
        b = zu - prob->d*objFactor - quad;
        
        if( sum >= b )
        {
            code = MIP_INFEASIBILITY;
            goto termination;
        }
        
        if( ninf == 1 )
        {
            //we have only 1 variable making sum infinity.
            const double coef = c[indinf]*objFactor;
            
            const bool aux = __tryUpdateBound(b-sum, coef, xtype[indinf], lx[indinf], ux[indinf] );
            
            if( aux )
                updtBounds = true;
        }
        else if(ninf == 0)
        {
            //running the variables updating the bounds
            
            
            for(int i = 0; i < n; i++)
            {
                const double coef = c[indinf]*objFactor;
                
                
                if( coef > 0.0 )
                {
                    const double sumw = sum - coef*lx[i];
                    
                    const bool aux = __tryUpdateBound(b-sumw, coef, xtype[i], lx[i], ux[i]);
                    
                    if(aux)
                        updtBounds = true;
                }
                else if( coef < 0.0 )
                {
                    const double sumw = sum - coef*ux[i];
                    
                    const bool aux = __tryUpdateBound(b-sumw, coef, xtype[i], lx[i], ux[i]);
                    
                    if(aux)
                        updtBounds = true;
                }
            }
        }
    }
    
    #endif
    
    code = 0;
    
    
termination:
    
    return code;
}



int MIP_Preprocessing_old::__preprocQuadConstr( const MIP_SparseMatrix& Q, const unsigned int nzsRow, const int *rowCols, const double *rowValues, double b, const bool leqConstr, const int* varTypes, double* lx, double* ux, bool &updtBounds, bool& redundant )
{
    const int n = Q.getNumberOfColumns();
    //const unsigned int nzsRow = nzsRow; //a.getNumberOfElements();
    const double qfactor = leqConstr ? 1.0 : -1.0;
    //const MIP_SparseElement *colAuxa = a.columns;
    
    bool infeasible;
    int r;
    double quad;
    
    
    updtBounds = false;
    redundant = false;
    
    
    r = __calcQminValue(Q, qfactor, lx, ux, auxValues, quad);
    
    if(r != 0)
    {
        #if MIP_DEBUG_MODE
            //MIP_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    //cout << "qfactor: " << qfactor << " quad: " << quad << endl;
    //for(int i = 0; i < n; i++)
        //cout << "auxValues[" << i << "]: " << auxValues[i] << "  ";
    //cout << endl;
    
    
    b = qfactor*b - quad;
    
    //cout << "b: " << b << endl;
    
    
    //auxValues has deduced coeficient from quadratic part, and quad has a constant value deduced from  quadratic part.
    for( unsigned int j = 0; j < nzsRow; j++ )
        auxValues[ rowCols[j] ] += qfactor * rowValues[j]; //auxValues[ colAuxa[j].getColumn() ] += qfactor *colAuxa[j].getValue();
    
    updtBounds = __preprocFullLinConstr(n, auxValues, b, varTypes, lx, ux, infeasible, redundant);
    
    
    redundant = false;
    
    if( infeasible )
    {
        return MIP_INFEASIBILITY;
    }
    
    
    return 0;
}



bool MIP_Preprocessing_old::__preprocFullLinConstr( const int n, const double *row, double b, const int* varTypes, double* lx, double* ux, bool& infeasible, bool& redundant)
{
    const double feasTol = in_abs_feas_tol;
    const double feasFactorTol = in_rel_feas_tol;
    
    bool updtBounds = false;
    int ninf = 0, indinf;
    double sum, sump = 0.0, summ = 0.0;
    
    #if MIP_DEBUG_MODE
        indinf = INT_MAX;
    #endif
    
    //for(int i = 0; i < n; i++)
        //cout << "row["<<i<<"]: " << row[i] << "  ";
    //cout << endl;
    
    infeasible = false;
    redundant = false;
    
    for(int i = 0; i < n; i++)
    {
        const double coef = row[i];
        
        if( coef > 0.0 )
        {
            if( lx[i] > -MIP_INFINITY  )
            {
                sump += coef*lx[i]; //ok, there is no guaratee that term is really positive, but we sum here anyway... we just want reduce numerical errors, not be perfect ;) (I do not want perform an if here to test if this term is positive or not...)
            }
            else
            {
                if(ninf > 0)
                    goto termination; //that is the second infinity. We cannot use that constraint to preprocess..
                
                
                ninf = 1;
                indinf = i;
            }
        }
        else if( coef < 0.0 )
        {
            if( ux[i] < MIP_INFINITY )
            {
                summ += coef*ux[i]; //ok, there is no guaratee that term is really negative, but we sum here anyway... we just want reduce numerical errors, not be perfect ;) (I do not want perform a n if here to test is this term is positive or not...)
            }
            else
            {
                if(ninf > 0)
                    goto termination; //that is the second infinity. We cannot use that constraint to preprocess..
                
                ninf = 1;
                indinf = i;
            }
        }
        
    }
    
    
    sum = sump + summ;
    
    //cout << "ninf: " << ninf << " indinf: " << indinf << " sum: " << sum << " sump: " << sump << " summ: " << summ << " b: " << b << endl;
    
    if( ninf == 0 && sum > b  + MIP_abs(feasFactorTol *b) + feasTol )
    {
        //that constraint is infeasible
        infeasible = true;
        goto termination;
    }
    
    
    
    if( ninf == 1 )
    {
        //we have only 1 variable making sum infinity.
        
        const bool aux = __tryUpdateBound( b-sum, row[indinf], varTypes[indinf], lx[indinf], ux[indinf] );
        
        if( aux )
            updtBounds = true;
    }
    else if( ninf == 0 )
    {
        //running the variables updating the bounds
        double constrp = 0.0, constrm = 0.0;
        
        for(int i = 0; i < n; i++)
        {
            const double coef = row[i];
            
            if( coef > 0.0 )
            {
                const double sumw = sum - coef*lx[i];
                
                const bool aux = __tryUpdateBound( b-sumw, coef, varTypes[i], lx[i], ux[i] );
                
                if(aux)
                    updtBounds = true;
                
                constrp += coef * ux[i];
            }
            else if( coef < 0.0 )
            {
                const double sumw = sum - coef*ux[i];
                
                const bool aux = __tryUpdateBound( b-sumw, coef, varTypes[i], lx[i], ux[i]);
                
                if(aux)
                    updtBounds = true;
                
                constrm += coef * lx[i];
            }
        }
        
        
        if( constrp + constrm <= b )
            redundant = true;
    }
    
    
    
termination:
    
    return updtBounds;
}



