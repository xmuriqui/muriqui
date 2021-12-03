
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <iostream>

#include <new>
#include "SPM_SparseMatrix.hpp"





#ifdef __FILE__
    #ifdef __LINE__
        #define SPM_DEF_GETFILELINE 1
    #endif
#endif


#ifdef SPM_DEF_GETFILELINE

    #define SPM_GETFILELINE  \
        " on file: " << __FILE__ << " line: " << __LINE__
#else
    #define SPM_GETFILELINE ""
#endif


#define SPM_PREPRINT "sparse matrix: "



#define SPM_PRINTMEMERROR std::cerr << SPM_PREPRINT << "Memory error" << SPM_GETFILELINE << std::endl


#define SPM_PRINTERROR std::cerr << SPM_PREPRINT << "Error" << SPM_GETFILELINE << std::endl


#define SPM_PRINTERRORNUMBER(number) std::cerr << SPM_PREPRINT << "Error " << number << SPM_GETFILELINE << std::endl


#define SPM_PRINTERRORMSG(msg) std::cerr << SPM_PREPRINT << msg << SPM_GETFILELINE << std::endl


#define SPM_getchar()  SPM_PRINTERRORMSG("stopped in a getchar"), getchar()












template <class matrixType>
spm::SPM_SparseRow<matrixType>::SPM_SparseRow()
{
    initialize();
}


template <class matrixType>
spm::SPM_SparseRow<matrixType>::~SPM_SparseRow()
{
    desallocateColumns();
} 


template <class matrixType>
int spm::SPM_SparseRow<matrixType>::allocateColumns(const unsigned int ncols)
{
    spm::SPM_SparseElement<matrixType> *p;
    
    
    if(ncols == 0)
    {
        desallocateColumns(); //realloc has an undefined behaviour when size is zero in new c++ versions...
    }
    else
    {
        //if columns is NULL realloc acts as malloc...
        p = (spm::SPM_SparseElement<matrixType> *) realloc(columns, ncols * sizeof(spm::SPM_SparseElement<matrixType>) );
        
        if(!p)
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTMEMERROR;
            #endif
            return SPM_MEMORY_ERROR;
        }
        
        columns = p;
        
        this->nElements = ncols;
    }
    
    
    return 0;
}





template <class matrixType>
int spm::SPM_SparseRow<matrixType>::copyStructureFrom(spm::SPM_SparseRow<matrixType> &other )
{
    int r;
    
    
    if( nElements != other.nElements )
    {
        r = allocateColumns( other.nElements );
        if(r != 0)
            return r;
    }
    
    for(unsigned int i = 0; i < nElements; i++)
        columns[i].setColumn( other.columns[i].getColumn() );
    
    return 0;
}


template <class matrixType>
int spm::SPM_SparseRow<matrixType>::copyFrom(spm::SPM_SparseRow<matrixType> &other )
{
    int r;
    
    
    if( nElements != other.nElements )
    {
        r = allocateColumns( other.nElements );
        if(r != 0)
            return r;
    }
    
    for(unsigned int i = 0; i < nElements; i++)
        columns[i] = other.columns[i];
    
    return 0;
}


template <class matrixType>
void spm::SPM_SparseRow<matrixType>::copyValuesInOrderFrom( spm::SPM_SparseRow<matrixType> &other  )
{
    for(unsigned int i = 0; i < nElements; i++)
        columns[i].setValue( other.columns[i].getValue() );
}




template <class matrixType>
void spm::SPM_SparseRow<matrixType>::copyTo(double* row, const double factor) const
{
    unsigned int i;
    
    for(i = 0; i < nElements; i++)
        row[ columns[i].getColumn() ] = columns[i].getValue();
    
    
    if( factor != 1.0 )
    {
        for(i = 0; i < nElements; i++)
            row[ columns[i].getColumn() ] *= factor;
    }
    
}



template <class matrixType>
void spm::SPM_SparseRow<matrixType>::desallocateColumns()
{
    if( columns )
    {
        free(columns);
        columns = NULL;
        nElements = 0;
    }
}


template <class matrixType>
double spm::SPM_SparseRow<matrixType>::evalTimesxt(const double* x) const
{
    unsigned int k;
    double v = 0.0;
    
    for(k = 0; k < nElements; k++)
        v += columns[k].v * x[ columns[k].col ];
    
    return v;
}


template <class matrixType>
void spm::SPM_SparseRow<matrixType>::initialize()
{
    nElements = 0;
    columns = NULL;
}



template <class matrixType>
int  spm::SPM_SparseRow<matrixType>::getElement(const unsigned int col, matrixType &value) const
{
    for(unsigned int i = 0; i < nElements; i++)
    {
        if( columns[i].getColumn() == col )
        {
            value = columns[i].getValue();
            return 0;
        }
    }
    
    value = 0.0;
    
    //the element is not in the matrix
    return SPM_BAD_DEFINITIONS;
}



template <class matrixType>
template <class intClass>
unsigned int spm::SPM_SparseRow<matrixType>:: __getStructureAndValues( intClass* cols, matrixType* values) const
{
    unsigned int k;
    
    for(k = 0; k < nElements; k++)
    {
        cols[k] = columns[k].col;
        values[k] = columns[k].v;
    }
    
    return k;
}



template <class matrixType>
unsigned int spm::SPM_SparseRow<matrixType>::getStructureAndValues( int* cols, matrixType* values) const
{
    return __getStructureAndValues(cols, values);
}


template <class matrixType>
unsigned int spm::SPM_SparseRow<matrixType>::getStructureAndValues( unsigned int* cols, matrixType* values) const
{
    return __getStructureAndValues(cols, values);
}





template <class matrixType>
template <class intClass>
unsigned int spm::SPM_SparseRow<matrixType>::__getStructure( intClass* cols) const
{
    unsigned int k;
    
    for(k = 0; k < nElements; k++)
        cols[k] = columns[k].col;
    
    return k;
}


template <class matrixType>
unsigned int spm::SPM_SparseRow<matrixType>::getStructure( int* cols) const
{
    return __getStructure(cols);
}



template <class matrixType>
unsigned int spm::SPM_SparseRow<matrixType>::getStructure( unsigned int* cols) const
{
    return __getStructure(cols);
}




template <class matrixType>
unsigned int spm::SPM_SparseRow<matrixType>::getStructure( bool *cols) const
{
    unsigned int k;
    
    
    for(k = 0; k < nElements; k++)
        cols[ columns[k].col ] = true;
    
    return k;
}



template <class matrixType>
unsigned int spm::SPM_SparseRow<matrixType>::getValues(matrixType *values) const
{
    unsigned int k;
    
    for(k = 0; k < nElements; k++)
        values[k] = columns[k].v;
    
    return k;
}



template <class matrixType>
bool spm::SPM_SparseRow<matrixType>::hasColumn( const unsigned int col, unsigned int *index) const
{
    for(unsigned int i = 0; i < nElements; i++)
    {
        if( columns[i].getColumn() == col )
        {
            if(index)
                *index = i;
            
            return true;
        }
    }
    
    return false;
}


template <class matrixType>
void spm::SPM_SparseRow<matrixType>::print() const
{
    for(unsigned int i = 0; i < nElements; i++)
    {
        std::cout << columns[i].getColumn() << ": " << columns[i].getValue() << " ";
    }
    std::cout << "\n";
}



template <class matrixType>
void spm::SPM_SparseRow<matrixType>:: multiplyAllElements( const matrixType value )
{
    for(unsigned int i = 0; i < nElements; i++)
        columns[i].setValue( value * columns[i].getValue() );
}



template <class matrixType>
int spm::SPM_SparseRow<matrixType>::reallocateColumns( const unsigned int ncols )
{
    spm::SPM_SparseElement<matrixType> *aux;
    
    aux = ( spm::SPM_SparseElement<matrixType> *) realloc( columns, ncols * sizeof(spm::SPM_SparseElement<matrixType>) );
    
    if( !aux )
        return SPM_MEMORY_ERROR;
    
    columns = aux;
    nElements = ncols;
    
    return 0;
}



template <class matrixType>
void spm::SPM_SparseRow<matrixType>:: setAllElements(const matrixType value)
{
    for( unsigned int k = 0; k < nElements; k++ )
        columns[k].setValue( value ); //we could use (*this)[k] = value;
}


template <class matrixType>
int spm::SPM_SparseRow<matrixType>:: setElement(const unsigned int col, const matrixType value)
{
    for( unsigned int k = 0; k < nElements; k++ )
    {
        if( columns[k].getColumn() == col )
        {
            columns[k].setValue( value );
            return 0;
        }
    }
    
    //the element is not in the row
    return SPM_BAD_DEFINITIONS;
}



template <class matrixType>
int spm::SPM_SparseRow<matrixType>:: setElementsByOrder(const unsigned int nel, const matrixType *values)
{
    for(unsigned int k = 0; k < nel; k++)
        columns[k].setValue( values[k] );
    
    
    return 0;
}



template <class matrixType>
spm::SPM_SparseMatrix<matrixType>::SPM_SparseMatrix()
{
    nrows = ncols = nElements = 0;
    symmetric = false;
    rows = NULL;
    
    initialize(0, 0, false);
}


template <class matrixType>
spm::SPM_SparseMatrix<matrixType>::SPM_SparseMatrix(const unsigned int nrows, const unsigned int ncols, const bool symmetric)
{
    initialize(nrows, ncols, symmetric);
}



template <class matrixType>
spm::SPM_SparseMatrix<matrixType>::~SPM_SparseMatrix()
{
    /*
    printf("destruindo sparseMatrix\n");
    SPM_getchar();
    SPM_desalocateSparseMatrix(*this); */
    
    desallocateMemory();
}




template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::accumulateLineInVector(const unsigned int line, const double factor, double* v) const
{
    spm::SPM_SparseElement<matrixType> *colAux = rows[line].columns;
    const int nz = rows[line].nElements;
    int i;
    
    for(i = 0; i < nz; i++ )
        v[ colAux[i].col ] += factor * colAux[i].v;
    
}


template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::accumulatedSum(const double factor, spm::SPM_SparseMatrix<matrixType> &A)
{
    int i, j, p, aux, aux2;
    spm::SPM_SparseElement<matrixType> *colAux, *colAux2;
    
    
    for(i = 0; i < A.nrows; i++)
    {
        aux = A.rows[i].nElements;
        colAux = A.rows[i].columns;
        
        aux2 = rows[i].nElements;
        colAux2 = rows[i].columns;
        
        for(j = 0; j < aux; j++)
        {
            for(p = 0; p < aux2; p++)
            {
                if( colAux2[p].col == colAux[j].col )
                {
                    colAux2[p].v += factor * colAux[j].v;
                    break;
                }
            }
        }
    }
}




//that method accumulates the sum between this matrix and factor*A. All positions in A should be also in this matrix. 
template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>:: accumulatedSumToMatrix( const double factor, matrixType *M, const bool considerSymmetry) const
{
    
    for(unsigned int i = 0; i < nrows; i++)
    {
        matrixType *rrow = &M[i*ncols];
        
        SPM_SparseRow<matrixType> &srow = rows[i];
        
        unsigned int nel = srow.getNumberOfElements(); 
        
        for(unsigned int j = 0; j < nel; j++)
        {
            rrow[ srow[j].getColumn() ] += factor * srow[j].getValue();
        }
        
        
        if( considerSymmetry )
        {
            for(unsigned int j = 0; j < nel; j++)
            {
                unsigned int row = srow[j].getColumn();
                
                if( row != i )
                    M[ row*ncols + i ] += factor * srow[j].getValue();
            }
        }
    }
    
}




//This functions assumes that the position is in the matrix. The current value in the position will be added to value parameter.
template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::addToElement( unsigned int row, unsigned int col, const matrixType value)
{
    unsigned int i;
    spm::SPM_SparseElement<matrixType> *colAux;

    if(symmetric )
    {
        //we work only on the lower triangle
        if( col > row )
        {
            i = col;
            col = row;
            row = i;
        }
    }

    colAux = rows[row].columns;

    for(i = 0; i < rows[row].nElements; i++)
    {
        if( colAux[i].col == col )
        {
            colAux[i].v += value;
            return 0;
        }
    }

    #if SPM_DEBUG_MODE
        SPM_PRINTERRORMSG("Element was not found in the sparse matrix!");
        SPM_getchar();
    #endif
    
    //the element is not in the matrix
    return SPM_BAD_DEFINITIONS;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::addNewRows(const unsigned int nrows)
{
    //int i;
    spm::SPM_SparseRow<matrixType> *p;
    
    p = (spm::SPM_SparseRow<matrixType> *) realloc( rows, (this->nrows + nrows) * sizeof(spm::SPM_SparseRow<matrixType>) );
    
    if(!p)
        return SPM_MEMORY_ERROR;
    
    rows = p;
    
    for(unsigned int i = this->nrows; i < this->nrows + nrows; i++)
        rows[i].initialize();
    
    
    this->nrows += nrows;
    
    return 0;
}



template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: checkRowStructure( unsigned int row, unsigned int nzs, int *cols, double *values ) const
{
    return SPM_checkRowStructure( nrows, ncols, symmetric, row, nzs, cols, values);
}



template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: checkRowStructure( unsigned int row, unsigned int nzs, unsigned int *cols, double *values ) const
{
    return SPM_checkRowStructure( nrows, ncols, symmetric, row, nzs, cols, values);
}





template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: checkTripleSparseStructure( unsigned int nzs, const int *rows, const int *cols, const double *values ) const
{
    return SPM_checkTripleSparseStructure( nrows, ncols, symmetric, nzs, rows, cols, values );
}



template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: checkTripleSparseStructure( unsigned int nzs, const unsigned int *rows, const unsigned int *cols, const double *values ) const
{
    return SPM_checkTripleSparseStructure( nrows, ncols, symmetric, nzs, rows, cols, values );
}





template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::allocateSparseRows(const unsigned int m)
{
    //do not use new because we can need realloc
    
    if(rows)
        free(rows);
    
    rows = (spm::SPM_SparseRow<matrixType>*) malloc( m * sizeof(spm::SPM_SparseRow<matrixType>) );
    if( !(rows) )
        return SPM_MEMORY_ERROR;

    nrows = m;
    for(unsigned int i = 0; i < m; i++)
        rows[i].initialize();
    

    return 0;
}


template <class matrixType>
bool spm::SPM_SparseMatrix<matrixType>:: compareSparseMatrixStructure(const spm::SPM_SparseMatrix<matrixType> &other) const
{
    bool aux;
    int i, j, k;
    spm::SPM_SparseElement<matrixType> *colAux, *colAux2;
    
    //we do not compair constant_calculated
    //if( nrows != other.nrows || ncols != other.ncols || nElements != other.nElements || symmetric != other.symmetric || constant != other.constant )
    if( nrows != other.nrows || ncols != other.ncols || nElements != other.nElements || symmetric != other.symmetric )
        return false;
    
    if((rows == NULL && other.rows != NULL) || (rows != NULL && other.rows == NULL))
        return false;
    
    for(i = 0; i < nrows; i++)
    {
        if( rows[i].nElements != other.rows[i].nElements )
            return false;
        
        colAux = rows[i].columns;
        colAux2 = other.rows[i].columns;
        
        
            for(j = 0; j < rows[i].nElements; j++)
            {
                aux = false;
                //the order of the columns can be diferent...
                for(k = 0; k < rows[i].nElements; k++)
                {
                    if( colAux[j].col == colAux2[k].col )
                    {
                        aux = true;
                        break;
                    }
                }
                
                if(!aux)
                    return false;
            }
    }
    
    return true;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::copyLine(const unsigned int sourceLine, const unsigned int destLine, SPM_SparseMatrix< matrixType >& M)
{
    int j, aux, code;
    spm::SPM_SparseElement<matrixType> *colAux, *colAux2;
    
    
    aux = M.rows[sourceLine].nElements;
    
    if( rows[destLine].nElements != aux )
    {
        if( rows[destLine].nElements > 0 )
        {
            nElements -= rows[destLine].nElements;
            free( rows[destLine].columns );
        }
        
        if( aux == 0 )
        {
            rows[destLine].nElements = 0;
            rows[destLine].columns = NULL;
        }
        else
        {
            rows[destLine].columns = (spm::SPM_SparseElement<matrixType> *) malloc( aux * sizeof(spm::SPM_SparseElement<matrixType>) );
            
            if( !rows[destLine].columns )
            {
                rows[destLine].nElements = 0;
                code = SPM_MEMORY_ERROR;
                goto desallocate_memory;
            }
            
            rows[destLine].nElements = aux;
            nElements += aux;
        }
    }
    
    
    colAux = M.rows[sourceLine].columns;
    colAux2 = rows[destLine].columns;
    
    for(j = 0; j < aux; j++)
        colAux2[j] = colAux[j];
    
    
    if(M.ncols > ncols)
        ncols = M.ncols; //we copy a line from M, so if M has more columns that our matrix, we need update ncols
    
    
    code = 0;
    
desallocate_memory:
    
    return code;
}



template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::copyParametersFrom(const spm::SPM_SparseMatrix<matrixType> & other)
{
    symmetric = other.symmetric;
    //constant = other.constant;
    //constant_calculated = other.constant_calculated;
    nrows = other.nrows;
    ncols = other.ncols;
    nElements = other.nElements;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::copyStructureFrom(const SPM_SparseMatrix& M)
{
    //int i, j;
    //spm::SPM_SparseElement<matrixType> *colAux, *colAuxOther;
    
    
    if( nrows != M.nrows )
    {
        desallocateMemory();
        
        const int r = allocateSparseRows(M.nrows);
        if(r != 0)
            return r;
    }
    
    copyParametersFrom(M);
    
    nElements = 0;
    
    for(unsigned int i = 0; i < nrows; i++)
    {
        const int r = rows[i].copyStructureFrom( M.rows[i] );
        if( r != 0 )
            return SPM_MEMORY_ERROR;
        
        nElements += rows[i].getNumberOfElements();
    }
    
    #if SPM_DEBUG_MODE
        assert( nElements == M.nElements );
    #endif
    
    return 0;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::copyMatrixFrom(SPM_SparseMatrix& M)
{
    int r, code;
    //spm::SPM_SparseElement<matrixType> *colAux, *colAuxOther;
    
    
    if( nrows != M.nrows )
    {
        desallocateMemory();
        
        r = allocateSparseRows(M.nrows);
        if(r != 0)
        {
            code = r;
            goto termination;
        }
    }
    
    copyParametersFrom(M);
    
    nElements = 0;
    
    for(unsigned int i = 0; i < nrows; i++)
    {
        r = rows[i].copyFrom( M.rows[i] );
        if( r != 0 )
        {
            code = SPM_MEMORY_ERROR;
            goto termination;
        }
        
        nElements += rows[i].getNumberOfElements();
    }
    
    #if SPM_DEBUG_MODE
        assert( nElements == M.nElements );
    #endif
    
    
    code = 0;
    
termination:
    
    return code;
}



template <class matrixType>
matrixType * spm::SPM_SparseMatrix<matrixType>:: copyMatrixTo(matrixType *Matrix, bool initializeWithZero, const double factor ) const
{
    unsigned int i, j;
    matrixType *M;
    
    if(Matrix)
    {
        M = Matrix;
    }
    else
    {
        M = (matrixType *) calloc( nrows * ncols , sizeof(matrixType) );
        if( !M )
        {
            return NULL;
        }
        
        initializeWithZero = false;
    }
    
    
    
    
    if( initializeWithZero )
    {
        unsigned int aux = nrows*ncols;
        for(i = 0; i < aux; i++)
            M[i] = 0.0;
    }
    
    
    for(i = 0; i < nrows; i++)
        rows[i].copyTo( &M[i*ncols], factor );
    
    
    if( symmetric )
    {
        for(i = 1; i < nrows; i++)
        {
            matrixType *row = &M[i*ncols];
            
            for(j = 0; j < i; j++)
                M[ j*ncols + i ] = row[j];
        }
    }
    
    
    /*for(i = 0; i < nrows; i++ )
    {
        colAux = rows[i].columns;
        aux = rows[i].nElements;
        
        k = i*nrows;
        for(j = 0; j < aux; j++)
            M[ k + colAux[j].col ] = colAux[j].v;
    
        if(symmetric)
        {
            for(j = 0; j < aux; j++)
                M[ colAux[j].col * nrows + i ] = colAux[j].v;
        }
    } */
    
    return M;
}







template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>:: copyMatrixTo(matrixType** M, const bool initializeWithZero, const double factor) const
{
    unsigned int i, j;
    
    
    if( initializeWithZero )
    {
        for(i = 0; i < nrows; i++)
            for(j = 0; j < ncols; j++)
                M[i][j] = 0.0;
    }
    
    
    for(i = 0; i < nrows; i++)
        rows[i].copyTo( M[i], factor );
    
    
    if( symmetric )
    {
        for(i = 1; i < nrows; i++)
        {
            matrixType *row = M[i];
            for( j = 0; j < i; j++ )
                M[j][i] = row[j];
        }
    }
    
}



template <class matrixType>
template <class intClass>
void
spm::SPM_SparseMatrix<matrixType>::__countRowsEachColumn(intClass *counts, const bool accumulate) const
{
    //spm::SPM_SparseElement<matrixType> *colAux;
    unsigned int i, j, aux;
    
    
    if(!accumulate)
    {
        for(i = 0; i < ncols; i++)
            counts[i] = 0;
    }
    
    for(i = 0; i < nrows; i++ )
    {
        spm::SPM_SparseRow<matrixType> &row = rows[i];
        
        aux =  row.getNumberOfElements(); //  rows[i].nElements;	colAux = row.columns;
        
        for(j = 0; j < aux; j++)
            counts[ row[j].getColumn() ]++; //counts[ colAux[j].col ]++;
    }
    
}



template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::countRowsEachColumn( int *counts, const bool accumulate) const
{
    return __countRowsEachColumn(counts, accumulate);
}



template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::countRowsEachColumn(unsigned int *counts, const bool accumulate) const
{
    return __countRowsEachColumn(counts, accumulate);
}


template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>:: deleteStructure()
{
    for(unsigned int i = 0; i < nrows; i++)
    {
        rows[i].nElements = 0;
        if( rows[i].columns != NULL )
        {
            free( rows[i].columns );
            rows[i].columns = NULL;
        }
    }
    nElements = 0;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: deleteRowStructure(const unsigned int row)
{
    if(row < 0 || row >= nrows)
        return SPM_BAD_DEFINITIONS;

    nElements -= rows[row].nElements;
    
    rows[row].desallocateColumns();
    return 0;
}


template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::desallocateMemory(void)
{
    unsigned int i;

    if(rows)
    {
        for(i = 0; i < nrows; i++)
            rows[i].desallocateColumns();
        
        free( rows );
        rows = NULL;
        nrows = 0;
        
        nElements = 0;
    }
}



template <class matrixType>
int  spm::SPM_SparseMatrix<matrixType>::getElement(unsigned int row, unsigned int col, matrixType& value) const
{
    unsigned int i, aux;
    //spm::SPM_SparseElement<matrixType> *colAux;
    

    if( symmetric )
    {
        //we work only on the lower triangle
        if( col > row )
        {
            i = col;
            col = row;
            row = i;
        }
    }

    //colAux = rows[row].columns;
    
    
    spm::SPM_SparseRow<matrixType> &myrow = rows[row];
    
    aux = myrow.getNumberOfElements();
    
    for(i = 0; i < aux; i++)
    {
        if( myrow[i].getColumn() == col ) //( colAux[i].col == col )
        {
            value = myrow[i].getValue(); //colAux[i].v;
            return 0;
        }
    }
    
    value = 0.0;
    
    //the element is not in the matrix
    return SPM_BAD_DEFINITIONS;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getFullRow(const unsigned int index, matrixType* values, const bool considerSymmetry) const
{
    spm::SPM_SparseRow<matrixType> &myrow = rows[index];
    const unsigned int nel = myrow.getNumberOfElements();
    
    
    SPM_setAllArray( ncols, values, 0.0 );
    
    
    for( unsigned int j = 0; j < nel; j++)
        values[ myrow[j].getColumn() ] = myrow[j].getValue();
    
    
    if( considerSymmetry && symmetric )
    {
        
        //since we only store lower triangle, we need find the positions would be in the upper one. Note those positions are above row index
        
        for( unsigned int j = index + 1; j < nrows; j++ )
        {
            spm::SPM_SparseRow<matrixType> &myrow = rows[j];
            const unsigned int nel = myrow.getNumberOfElements();
            
            for(unsigned int k = 0; k < nel; k++)
            {
                if( myrow[k].getColumn() == index )
                {
                    values[j] = myrow[k].getValue();
                    break; //we assume indices are not duplicated
                }
            }
            
        }
        
    }
    
    
    return 0;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getFullRowAccumulation( const unsigned int row, matrixType* values) const
{
    unsigned int j, aux;
    spm::SPM_SparseRow<matrixType> &myrow = rows[row];
    
    
    aux = myrow.getNumberOfElements(); 
    
    for(j = 0; j < aux; j++)
        values[ myrow[j].getColumn() ] += myrow[j].getValue(); 
    
    return 0;
}





template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowStructureAndValues(const unsigned int row, int &nzs, int* jCol, matrixType* vals) const
{
    nzs = rows[row].getStructureAndValues(jCol, vals);
    
    return nzs;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowStructureAndValues(const unsigned int row, unsigned int &nzs, unsigned int* jCol, matrixType* vals) const
{
    nzs = rows[row].getStructureAndValues(jCol, vals);
    
    return nzs;
}



template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowStructure(const unsigned int line, unsigned int* jCol) const
{
    return rows[line].getStructure(jCol);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowStructure(const unsigned int line, int* jCol)  const
{
    return rows[line].getStructure(jCol);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowStructure(const unsigned int line, bool* cols, bool accumulate)  const
{
    if(!accumulate)
    {
        //for(unsigned int i = 0; i < ncols; i++)
            //cols[i] = false;
        SPM_setAllArray( ncols, cols, false );
    }
    
    return rows[line].getStructure(cols);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowValues(const unsigned int line, matrixType *vals) const
{
    return rows[line].getValues(vals);
}



/*template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: getNumberOfElements()
{
    return nElements;
}*/


template <class matrixType>
template <class intType>
int spm::SPM_SparseMatrix<matrixType>:: __getStructure(intType *iRow, intType *jCol)  const
{
    unsigned int i, j, k = 0, aux;
    //spm::SPM_SparseElement<matrixType> *colAux;
    
    
    for(i = 0; i < nrows; i++)
    {
        spm::SPM_SparseRow<matrixType> &row = rows[i];
        
        aux = row.getNumberOfElements(); //rows[i].nElements; colAux = rows[i].columns;
        
        for(j = 0; j < aux; j++)
        {
            iRow[k] = i;
            jCol[k] = row[j].getColumn(); //colAux[j].col;
            k++;
        }
    }
    
    return k;
}



template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: getStructure(int *iRow, int *jCol) const
{
    return __getStructure( iRow, jCol );
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: getStructure(unsigned int *iRow, unsigned int *jCol)  const
{
    return __getStructure( iRow, jCol );
}



template <class matrixType>
template <class intType>
int spm::SPM_SparseMatrix<matrixType>:: __getStructureAndValues(intType* iRow, intType* jCol, matrixType* vals) const
{
    unsigned int i, j, k = 0, aux;
    //spm::SPM_SparseElement<matrixType> *colAux;
    

    for(i = 0; i < nrows; i++)
    {
        spm::SPM_SparseRow<matrixType> &row = rows[i];
        
        
        aux = row.getNumberOfElements(); //rows[i].nElements; //colAux = rows[i].columns;
        
        for(j = 0; j < aux; j++)
        {
            iRow[k] = i;
            jCol[k] = row[j].getColumn(); //colAux[j].col;
            vals[k] = row[j].getValue(); //colAux[j].v;
            
            k++;
        }
    }
    
    return k;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: getStructureAndValues(int* iRow, int* jCol, matrixType* vals)  const
{
    return __getStructureAndValues(iRow, jCol, vals);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: getStructureAndValues(unsigned int* iRow, unsigned int* jCol, matrixType* vals)  const
{
    return __getStructureAndValues(iRow, jCol, vals);
}


/*
* If any of the vectors row, cols or values is NULL, we allocate memory for the one. Otherwise, we ASSUME that a enough quantity of memory is already allocated.
*/
template <class matrixType>
template <class intType>
int spm::SPM_SparseMatrix<matrixType>::__getTripleSparseFormat(intType** rows, intType** cols, matrixType** values) const
{
    unsigned int i, j, k = 0, aux;
    int code;
    //spm::SPM_SparseElement<matrixType> *colAux;
    
    if( *rows == NULL )
    {
        *rows = new (std::nothrow) intType[nElements];
        
        if(!*rows)
        {
            code = SPM_MEMORY_ERROR;
            goto termination;
        }
    }
    
    if( *cols == NULL )
    {
        *cols = new (std::nothrow) intType[nElements];
        
        if(!*cols)
        {
            code = SPM_MEMORY_ERROR;
            goto termination;
        }
    }
    
    if( *values == NULL )
    {
        *values = new (std::nothrow) double[nElements];
        
        if(!*values)
        {
            code = SPM_MEMORY_ERROR;
            goto termination;
        }
    }
    
    
    for(i = 0; i < nrows; i++)
    {
        spm::SPM_SparseRow<matrixType> &row = this->rows[i];
        
        aux = row.getNumberOfElements(); //this->rows[i].nElements; colAux = this->rows[i].columns;
        
        for(j = 0; j < aux; j++, k++)
        {
            (*rows)[k] = i;
            (*cols)[k] = row[j].getColumn(); //colAux[j].col;
            (*values)[k] = row[j].getValue(); //colAux[j].v;	
        }
    }
    
    code = k;
    
termination:
    
    return code;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getTripleSparseFormat(int** rows, int** cols, matrixType** values) const
{
    return __getTripleSparseFormat( rows, cols, values );
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getTripleSparseFormat(unsigned int** rows, unsigned int** cols, matrixType** values) const
{
    return __getTripleSparseFormat( rows, cols, values );
}






template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getValues(matrixType *values) const
{
    unsigned int i, j, k = 0, aux;
    //spm::SPM_SparseElement<matrixType> *colAux;

    for(i = 0; i < nrows; i++)
    {
        spm::SPM_SparseRow<matrixType> &row = rows[i];
        
        aux = row.getNumberOfElements(); //rows[i].nElements; colAux = rows[i].columns;
        for(j = 0; j < aux; j++)
        {
            values[k] = row[j].getValue(); //colAux[j].v;
            k++;
        }
    }
    
    return k;
}



template <class matrixType>
bool spm::SPM_SparseMatrix<matrixType>::hasIndex(const unsigned int row, const unsigned int col) const
{
    if( row >= nrows )
        return false;
    
    return rows[row].hasColumn(col);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::initialize(const unsigned int nrows, const unsigned int ncols, const bool symmetric)
{
    rows = NULL;
    this->symmetric = symmetric;
    //constant = constant_calculated = false;
    nElements = 0;
    
    this->ncols = ncols;
    
    
    if(nrows > 0)
        return allocateSparseRows(nrows);
    else
    {
        this->nrows = 0;
        return 0;
    }
    
}


template <class matrixType>
double spm::SPM_SparseMatrix<matrixType>::rowEvaluation(const unsigned int row, const matrixType* x) const
{
    return rows[row].evalTimesxt(x);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::mergeStructures(SPM_SparseMatrix &M, const bool storeOrigCoefs)
{
    int *auxCols = NULL;
    int i, j, k, aux, aux2, aux3, cols, nz, code;
    spm::SPM_SparseRow<matrixType> *p;
    spm::SPM_SparseElement<matrixType> *colAux, *colAux2, *colAux3;
    const int norows = nrows;
    const int newncols = SPM_max( M.ncols, ncols);
    
    if( M.nrows > nrows )
    {
        p = (spm::SPM_SparseRow<matrixType> *) realloc(rows, M.nrows * sizeof(spm::SPM_SparseRow<matrixType>) );
        if(!p)
        {
            code = SPM_MEMORY_ERROR;
            goto desallocate_memory;
        }
        rows = p;
        nrows = M.nrows;
    }
    
    
    auxCols = (int *) malloc( newncols * sizeof(int));
    if(!auxCols)
    {
        code = SPM_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    
    
    aux3 = SPM_min(norows, M.nrows);
    for(i = 0; i < aux3; i++)
    {
        aux = M.rows[i].nElements;
        if( aux > 0 )
        {
            colAux = M.rows[i].columns;
            
            if( rows[i].nElements == 0 )
            {
                copyLine(i, i, M);
            }
            else
            {
                cols = symmetric && M.symmetric ? i+1 : newncols;
                
                
                for(j = 0; j < cols; j++)
                    auxCols[j] = 0;
                
                for(j = 0; j < aux; j++ )
                    auxCols[ colAux[j].col ] = 1;
                
                aux2 = rows[i].nElements;
                colAux2 = rows[i].columns;
                
                for(j = 0; j < aux2; j++)
                    auxCols[ colAux2[j].col ] = 1;
                
                nz = 0;
                for(j = 0; j < cols; j++)
                    nz += auxCols[j];
                
                
                if( nz > aux2 )
                {
                    colAux3 = (spm::SPM_SparseElement<matrixType> *) malloc( nz * sizeof(spm::SPM_SparseElement<matrixType>) );
                    
                    if(!colAux3)
                    {
                        code = SPM_MEMORY_ERROR;
                        goto desallocate_memory;
                    }
                    
                    for(j = k = 0; j < cols; j++)
                    {
                        if( auxCols[j] )
                        {
                            colAux3[k].col = j;
                            colAux3[k].v = 0.0;
                            k++;
                        }
                    }
                    
                    if(storeOrigCoefs)
                    {
                        for(j = 0; j < aux; j++)
                        {
                            for(k = 0; k < nz; k++)
                            {
                                if( colAux3[k].col == colAux[j].col )
                                {
                                    colAux3[k].v = colAux[j].v;
                                    break;
                                }
                                else if( colAux3[k].col < colAux[j].col )
                                {
                                    break; //the positions are ordered...
                                }
                            }
                        }
                        
                        for(j = 0; j < aux2; j++)
                        {
                            for(k = 0; k < nz; k++)
                            {
                                if( colAux3[k].col == colAux2[j].col )
                                {
                                    colAux3[k].v += colAux2[j].v;
                                    break;
                                }
                                else if( colAux3[k].col > colAux2[j].col )
                                {
                                    break; //the positions are ordered...
                                }
                            }
                        }
                    }
                    
                    free(rows[i].columns);
                    rows[i].columns = colAux3;
                    nElements = nElements - rows[i].nElements + nz;
                    rows[i].nElements = nz;
                }
                else if(storeOrigCoefs)
                {//the columns in that line are the same for both matrices...
                    for(j = 0; j < aux; j++)
                        addToElement(i, colAux[j].col, colAux[j].v);
                }
            }
            
            
        }
    }
    
    if( norows < M.nrows )
    {
        for(i = norows; i < M.nrows; i++)
        {
            rows[i].columns = (spm::SPM_SparseElement<matrixType> *) malloc( M.rows[i].nElements * sizeof(spm::SPM_SparseElement<matrixType>) );
            rows[i].nElements = M.rows[i].nElements;
            nElements += rows[i].nElements;
            
            copyLine(i, i, M);
        }
    }
    
    
    if( M.ncols > ncols )
        ncols = M.ncols;
    
    code = 0;
    
desallocate_memory:
    
    if(auxCols) free(auxCols);
    
    return code;
}



template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::printAllMatrix(void) const
{
    unsigned int i, j;
    double c;
    
    for(i = 0; i < nrows; i++)
    {
        for(j = 0; j < ncols; j++)
        {
            getElement(i, j, c);
            printf("%0.6f, ", c);
        }
        printf(";\n");
    }
}


template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::printSparseMatrix(const bool showEmptyLines) const
{
    unsigned int i, j;
    
    std::cout << "SPM_SparseMatrix:: nrows: " << nrows << " ncols: " << ncols << " nElements: " << nElements << " symmetric: " << symmetric << "\n";
    
    for(i = 0; i < nrows; i++)
    {
        if( rows[i].getNumberOfElements() > 0 || showEmptyLines)
        {
            //printf("%d=> ", i);
            std::cout << i << "=> ";
            
            for(j = 0; j < rows[i].getNumberOfElements(); j++)
            {
                //printf("%d: %f   ", rows[i].columns[j].col, rows[i].columns[j].v);
                
                std::cout << rows[i][j].getColumn() << ": " << rows[i][j].getValue() << "   ";
            }
            std::cout << "\n";
        }
        else if( showEmptyLines )
        {
            //printf("%d=> \n", i);
            std::cout << i << "=> " << std::endl;
        }
    }
}



template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::multiplyAllElements(const matrixType value)
{
    unsigned int i, j, aux;
    spm::SPM_SparseElement<matrixType> *colAux;
    
    for(i = 0; i < nrows; i++)
    {
        aux = rows[i].nElements;
        colAux = rows[i].columns;
        
        for(j = 0; j < aux; j++)
            colAux[j].v *= value;
    }
}



template <class matrixType>
double spm::SPM_SparseMatrix<matrixType>::quadraticEvaluation(const matrixType *x) const
{
    unsigned int i, j, k, aux;
    spm::SPM_SparseElement<matrixType> *colAux;
    double v = 0.0;
    
    #if SPM_DEBUG_MODE
        assert(symmetric);
    #endif
    
    for(i = 0; i < nrows; i++)
    {
        aux = rows[i].nElements;
        colAux = rows[i].columns;
        
        for(j = 0; j < aux; j++)
        {
            k = colAux[j].col;
            if(i == k)
                v += 0.5* x[k] * x[k] * colAux[j].v;
            else
                v += x[i] * x[k] * colAux[j].v;
        }
    }
    
    return v;
}


template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::quadraticGradientEvaluation( const matrixType* x, matrixType* grad, const bool accumulate) const
{
    unsigned int i, j, k, aux;
    spm::SPM_SparseElement<matrixType> *colAux;
    
    
    if( !accumulate )
    {
        for(i = 0; i < nrows; i++)
            grad[i] = 0.0;
    }
    
    
    for(i = 0; i < nrows; i++)
    {
        aux = rows[i].nElements;
        colAux = rows[i].columns;
        
        for(j = 0; j < aux; j++)
        {
            k = colAux[j].col;
            
            grad[i] += x[k] * colAux[j].v;
            
            if(i != k) //nondiagonal positions
                grad[k] += x[i] * colAux[j].v;
        }
    }
}


template <class matrixType>
double spm::SPM_SparseMatrix<matrixType>:: quadraticEvaluationAndGradient( const matrixType* x, matrixType* grad, const bool accumulate) const
{
    double value = 0.0;
    spm::SPM_SparseElement<matrixType> *colAux;
    
    
    #if SPM_DEBUG_MODE
        assert(symmetric);
    #endif
    
    
    if( !accumulate )
    {
        for(unsigned int i = 0; i < nrows; i++)
            grad[i] = 0.0;
    }
    
    
    
    for(unsigned int i = 0; i < nrows; i++)
    {
        unsigned int aux = rows[i].nElements;
        colAux = rows[i].columns;
        
        for(unsigned int j = 0; j < aux; j++)
        {
            const unsigned int &k = colAux[j].col;
            
            const double v1 = x[k] *colAux[j].v;
            grad[i] += v1;
            
            if(i == k) 
            {
                value += 0.5 * x[k] * v1;
            }
            else //nondiagonal positions
            {
                value += x[i] * v1;
                grad[k] += x[i] * colAux[j].v;
            }
        }
    }
    
    return value;
}



template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::removeCols( const unsigned int ncols, const unsigned int *cols )
{
    return __removeCols( ncols, cols );
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::removeCols( const unsigned int ncols, const int *cols )
{
    return __removeCols( ncols, cols );
}


template <class matrixType>
template <class intType>
int spm::SPM_SparseMatrix<matrixType>::__removeCols( const unsigned int ncols, const intType *cols )
{
    bool *remove = NULL;
    int code;
    
    
    remove = (bool *) calloc( this->ncols, sizeof(bool) ) ;
    if(!remove)
    {
        code = SPM_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    
    for(unsigned int i = 0; i < ncols; i++)
        remove[ cols[i] ] = true;
    
    
    removeCols(remove);
    
    
    code = 0;
    
desallocate_memory:
    
    if(remove)	free(remove);
    
    return code;
}



template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::removeCols( const bool *cols )
{
    for( unsigned int i = 0; i < ncols; i++ )
        removeCol( i );
}




template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::removeCol( const unsigned int col )
{
    bool resize;
    unsigned int rstart = symmetric ? col : 0;
    
    
    
    for( unsigned int i = rstart; i < nrows; i++)
    {
        SPM_SparseRow<matrixType> &row = rows[i];	
        unsigned int aux = row.getNumberOfElements();
        
        
        for( unsigned int j = 0; j < aux; j++ )
        {
            const unsigned int c = row[j].getColumn();
            
            
            if( c > col )
            {
                row[j].setColumn( c - 1 ); //we have to decrease indexes of all columns greater than c...
            }
            if( c == col )
            {
                //shifting values...
                for(unsigned int k = j+1; k < aux; k++)
                    row[k-1] = row[k];
                
                aux--;
            }
        }
        
        if( aux < row.getNumberOfElements() )
        {
            row.reallocateColumns( aux ); //we do not care, if reallocation was not a success
        }
        
    }
    
    //return 0;
}







template <class matrixType>
template <class intType>
int spm::SPM_SparseMatrix<matrixType>::__removeRows(const unsigned int nlines, const intType* lines)
{
    bool *remove = NULL;
    int code;
    
    
    remove = (bool *) calloc( nrows, sizeof(bool) ) ;
    if(!remove)
    {
        code = SPM_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    
    for(unsigned int i = 0; i < nlines; i++)
        remove[ lines[i] ] = true;
    
    
    removeRows(remove);
    
    
    code = 0;
    
desallocate_memory:
    
    if(remove)	free(remove);
    
    return code;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::removeRows(const unsigned int nlines, const int* lines)
{
    return __removeRows(nlines, lines);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::removeRows(const unsigned int nlines, const unsigned int* lines)
{
    return __removeRows(nlines, lines);
}




template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::removeRows( const bool* lines )
{
    unsigned int i, k;
    spm::SPM_SparseRow<matrixType> *aux;
    
    for(i = nrows; i >= 1; i--)
    {
        unsigned int ind = i -1;
        
        if( lines[ind] )
        {
            nElements -= rows[ind].nElements;
            
            rows[ind].desallocateColumns();
            
            for(k = ind + 1; k < nrows; k++)
                rows[k-1] = rows[k];
            
            nrows--;
        }
    }
    
    
    
    if( nrows > 0 )
    {
        aux = (spm::SPM_SparseRow<matrixType> *) realloc( rows, nrows * sizeof(spm::SPM_SparseRow<matrixType>) );
        
        if(aux)
            rows = aux;
        //if we get fail to decrease the size, there is no problem...
    }
    else
    {
        free(rows);
        rows = NULL;
    }
    
    
    //symmetric = false; //we let user decides if change the flag...
}



template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::setAllElementsInARow(const unsigned int line, const matrixType value)
{
    return rows[line].setAllElements(value);
    
    
    /*int j, aux;
    spm::SPM_SparseElement<matrixType> *colAux;
    
    colAux = rows[line].columns;
    aux    = rows[line].nElements;
    
    for(j = 0; j < aux; j++)
        colAux[j].v = value; */
}



template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>:: setAllSparseMatrix(const matrixType value)
{
    for(unsigned int i = 0; i < nrows; i++)
        rows[i].setAllElements( value );
    
}



template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: setElement( unsigned int row, unsigned int col, const matrixType value, const bool printErrMsg)
{

    if( symmetric )
    {
        //we work only on the lower triangle
        if( col > row )
            SPM_swap(col, row);
    }
    
    return setElementNoCheckSymmetry(row, col, value, printErrMsg);
}


//in this method, we do not check if col <= row in symmetric matrices...
template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: setElementNoCheckSymmetry(unsigned int row, unsigned int col, const matrixType value, const bool printErrMsg)
{
    const int r = rows[row].setElement(col, value);
    
    if(r != 0 && printErrMsg)
    {
        std::cerr << "Warning: element not found on spm::SPM_SparseMatrix<matrixType>::NoCheckSymmetry. row:" << row << " col: " << col << std::endl;
        
        //this->printSparseMatrix();
        //SPM_getchar();
    }
    
    return r;
}



template <class matrixType>
template <class intType>
//warning: that method store a single line. If that line  already exists, it will be overwritten
int spm::SPM_SparseMatrix<matrixType>:: __setRowStructure(  const unsigned int row, const unsigned int nzs, const intType* cols)
{
    unsigned int i;
    int code;
    
    
    if( row >= nrows )
        return SPM_BAD_DEFINITIONS; //do not use goto here
    
    spm::SPM_SparseRow<matrixType> &myrow = rows[row];
    
    nElements -= myrow.nElements;
    
    
    if( myrow.nElements != nzs )
    {
        code = myrow.allocateColumns( nzs );
        if( code != 0 )
        {
            goto termination;
        }
    }
    
    for(i = 0; i < nzs; i++)
    {
        myrow[i].setColumn( cols[i] );
        myrow[i].setValue( 0.0 );
    }
    
    
    code = 0;
    
termination:

    nElements += rows[row].nElements;

    return code;
}



template <class matrixType>
//warning: that method store a single line. If that line  already exists, it will be overwritten
int spm::SPM_SparseMatrix<matrixType>::setRowStructure(const unsigned int row, const unsigned int sizeCols, const bool *cols)
{
    int r;
    unsigned int k = 0, nz = 0;
    SPM_SparseRow<matrixType> &sprow = rows[row];
    
    
    nElements -= sprow.nElements;
    
    sprow.desallocateColumns(); //deleteRowStructure(row);
    
    
    for(unsigned int j = 0; j < sizeCols; j++)
    {
        if( cols[j] )
            nz++;
    }
    
    r = rows[row].allocateColumns(nz);
    if( r != 0 )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    nElements += nz;
    
    for(unsigned int j = 0; j < sizeCols; j++)
    {
        if( cols[j] )
        {
            sprow[k].setColumn(j);
            sprow[k].setValue(0);
            k++;
        }
    }
    
    #if SPM_DEBUG_MODE
        assert( k == nz );
    #endif
    
    return 0;
}



template <class matrixType>
//warning: that method store a single line. If that line  already exists, it will be overwritten
int spm::SPM_SparseMatrix<matrixType>::setRowStructure( const unsigned int row, const unsigned int numberOfNZElements, const unsigned int* cols)
{
    return __setRowStructure( row, numberOfNZElements, cols );
}



template <class matrixType>
//warning: that method store a single line. If that line  already exists, it will be overwritten
int spm::SPM_SparseMatrix<matrixType>::setRowStructure(const unsigned int row, const unsigned int numberOfNZElements, const int* cols)
{
    return __setRowStructure( row, numberOfNZElements, cols );
}





template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::setRowStructure(const unsigned int rowIndex, SPM_SparseRow< matrixType >& row)
{
    int code;
    
    if(rowIndex >= nrows)
        return SPM_BAD_DEFINITIONS; //do not use goto here
    
    nElements -= rows[rowIndex].nElements;
    
    code = rows[rowIndex].copyStructureFrom( row );
    
    nElements += rows[rowIndex].nElements;
    
    return code;
}



template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::setRowStructureAndValues( const unsigned int rowIndex, SPM_SparseRow< matrixType >& row)
{
    int code;
    unsigned int nel = row.nElements;
    
    if(rowIndex >= nrows)
        return SPM_BAD_DEFINITIONS; //do not use goto here
    
    SPM_SparseRow<matrixType> &myrow = rows[rowIndex];
    
    nElements -= myrow.nElements;
    
    code = myrow.allocateColumns( nel );
    
    nElements += myrow.nElements;
    
    if( code != 0 )
        return code;
    
    
    for( unsigned int i = 0; i < nel; i++)
        myrow[i] = row[i];
    
    
    return 0;
}




template <class matrixType>
template <class intType>
//warning: that method store a single linear constraint. If that constraints have already exists, it will be overwritten
int spm::SPM_SparseMatrix<matrixType>:: __setRowStructureAndValues(const unsigned int row, const unsigned int nzs, const intType* cols, const matrixType* vals)
{
    int code;

    if(row >= nrows)
        return SPM_BAD_DEFINITIONS; //do not use goto here
    
    spm::SPM_SparseRow<matrixType> &myrow = rows[row];
    
    nElements -= myrow.getNumberOfElements();
    
    
    if( nzs != myrow.getNumberOfElements() )
    {
        code = myrow.allocateColumns( nzs );
        
        if( code != 0 )
        {
            goto termination;
        }
    }
    
    
    for(unsigned int i = 0; i < nzs; i++)
    {
        myrow[i].setColumn( cols[i] );
        myrow[i].setValue( vals[i] );
    }
    
    
    code = 0;
    
    
termination:

    nElements += rows[row].nElements;

    return code;
}



template <class matrixType>
//warning: that method store a single linear constraint. If that constraints have already exists, it will be overwritten
int spm::SPM_SparseMatrix<matrixType>:: setRowStructureAndValues(const unsigned int row, const unsigned int numberOfNZElements, const int* cols, const matrixType* vals)
{
    return __setRowStructureAndValues( row, numberOfNZElements, cols, vals );
}



template <class matrixType>
//warning: that method store a single linear constraint. If that constraints have already exists, it will be overwritten
int spm::SPM_SparseMatrix<matrixType>::setRowStructureAndValues(const unsigned int row, const unsigned int numberOfNZElements, const unsigned int* cols, const matrixType* vals)
{
    return __setRowStructureAndValues( row, numberOfNZElements, cols, vals );
}




template <class matrixType>
//warning: that method store a single line. If the line have already exists, it will be overwritten. If the matrix is symmetric, just the elements in lower triangle elements are setted
int spm::SPM_SparseMatrix<matrixType>:: setRowStructureAndValues(const unsigned int row, const matrixType* a, unsigned int ncols, const double zeroTol)
{
    int code;
    unsigned int aux = 0, realncols = 0;
    
    if(ncols == 0)
        ncols = this->ncols;
    
    
    ncols = symmetric ? row + 1 : ncols;
    
    
    if(row >= nrows)
        return SPM_BAD_DEFINITIONS; //do not use goto here
    
    spm::SPM_SparseRow<matrixType> &myrow = rows[row];
    
    
    nElements -= myrow.nElements;
    
    
    for( unsigned int i = 0; i < ncols; i++ )
    {
        if( SPM_abs( a[i] ) > zeroTol )
            realncols++;
    }
    
    
    if( myrow.getNumberOfElements() != realncols )
    {
        code = myrow.allocateColumns(realncols);

        if( code != 0 )
        {
            goto desallocateMemory;
        }
    }

    for(unsigned int i = 0; i < ncols; i++)
    {
        if( SPM_abs( a[i] ) > zeroTol )
        {
            myrow[aux].setColumn( i );
            myrow[aux].setValue( a[i] );
            aux++;
        }
    }
    
    
    
    #if SPM_DEBUG_MODE
        if( aux != rows[row].getNumberOfElements() )
        {
            std::cerr << "Para tudo! Inconsistencia na contagem de elementos na definicao de uma unica restricao linear" << std::endl;
            std::cerr << "Please report this mensage and this code for " << SPM_EMAIL << std::endl;
        }
    #endif
    
    
    code = 0;

desallocateMemory:

    nElements += rows[row].nElements;

    return code;
}



template <class matrixType>
template <class intType>
int spm::SPM_SparseMatrix<matrixType>::__setStructureAndValues(const unsigned int nzs, const intType* iRow, const intType* jCol, const matrixType* vals, const bool reset)
{
    bool hadels;
    
    int code;
    unsigned int row, col;
    unsigned int *aux, aux3;
    
    unsigned int pos;
    double v = 0.0;
    
    
    
    
    aux = (unsigned int*) calloc( nrows, sizeof(unsigned int) );
    if( !aux )
    {
        code = SPM_MEMORY_ERROR;
        goto desallocateMemory;
    }
    
    
    if(reset)
        deleteStructure();
    
    
    hadels = nElements > 0;
    
    
    
    for(unsigned i = 0; i < nzs; i++)
    {
        if( symmetric && jCol[i] > iRow[i] )
        {
            //we store only the lower triangle...
            row = jCol[i];
            col = iRow[i];
        }
        else
        {
            row = iRow[i];
            col = jCol[i];
        }
        
        
        if( !rows[row].hasColumn(col) )
        {
            aux[row]++;
        }
    }
    
    
    
    //now, we allocating space for the new elements
    for(unsigned int i = 0; i < nrows; i++)
    {
        if( aux[i] == 0 )
            continue;
        
        aux3 = rows[i].getNumberOfElements();
        code = rows[i].allocateColumns(aux[i] +aux3);
        
        if( code != 0 )
            goto desallocateMemory;
        
        aux[i] = aux3;//now, aux has the first empty index to use in each row
        
        rows[i].setNumberOfElements( aux3 );
        
        /*if( rows[i].nElements == 0 )
        {
            rows[i].columns = (spm::SPM_SparseElement<matrixType> *) malloc( aux[i]*sizeof(spm::SPM_SparseElement<matrixType>) );

            if( !rows[i].columns )
            {
                code = SPM_MEMORY_ERROR;
                goto desallocateMemory;
            }

            rows[i].nElements = aux[i];
            aux[i] = 0; //aux has the first empty index to use in each row.
        }
        else
        {
            aux3 = rows[i].nElements + aux[i];
            colAux = (spm::SPM_SparseElement<matrixType> *) realloc( rows[i].columns, aux3 * sizeof(spm::SPM_SparseElement<matrixType>) );

            if(!colAux)
            {
                code = SPM_MEMORY_ERROR;
                goto desallocateMemory;
            }

            rows[i].columns = colAux;

            aux[i] = rows[i].nElements; //aux has the first empty index to use in each row.
            rows[i].nElements = aux3;
        } */
    }
    
    
    //storing the new elements. Here, aux[i] has the first empty index to use in each row i.

    
    
    for(unsigned int i = 0; i < nzs; i++)
    {
        if( symmetric && jCol[i] > iRow[i] )
        {
            //we store only the lower triangle...
            row = jCol[i];
            col = iRow[i];
        }
        else
        {
            row = iRow[i];
            col = jCol[i];
        }
        
        if( vals )
            v = vals[i];
        
        spm::SPM_SparseRow<matrixType> &myrow = rows[row];
        
        
        if( hadels && myrow.hasColumn(col, &pos) )
        {
            myrow[pos].setValue( v );
        }
        else
        {
            myrow[ aux[row] ].setColumn( col );
            myrow[ aux[row] ].setValue( v );
            
            aux[ row ]++;
            
            myrow.nElements++;
            nElements++;
        }
    }
    
    
    
    
    #if SPM_DEBUG_MODE
        if( !hadels )
            assert( nElements == nzs );
    #endif
    
    
    
    #if SPM_DEBUG_MODE
        
        aux3 = 0;
        
        for(unsigned int i = 0; i < nrows; i++)
            aux3 += rows[i].getNumberOfElements();
        
        assert( aux3 == nElements );
    #endif
    
    
    code = 0;

desallocateMemory:

    if(aux)		free(aux);

    return code;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::setStructureAndValues(const unsigned int nzs, const int* iRow, const int* jCol, const matrixType* vals, const bool reset)
{
    return __setStructureAndValues(nzs, iRow, jCol, vals, reset);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::setStructureAndValues(const unsigned int nzs, const unsigned int* iRow, const unsigned int* jCol, const matrixType* vals, const bool reset)
{
    return __setStructureAndValues(nzs, iRow, jCol, vals, reset);
}





template <class matrixType>
//warning: that method overwrittes possible lines already alllocated
int spm::SPM_SparseMatrix<matrixType>::setStructureAndValues(matrixType** A, const double zero_tol, unsigned int nrows, unsigned int ncols)
{
    int r;
    
    if( nrows == 0 )
        nrows = this->nrows;
    
    if( ncols == 0 )
        ncols = this->ncols;
    
    
    nElements = 0;

    //counting the nonzeros in the matriz
    for(unsigned int i = 0; i < nrows; i++)
    {
        SPM_SparseRow<matrixType> &row = rows[i];
        row.setNumberOfElements(0);
        
        const double *ai = A[i];
        
        unsigned int maxj = symmetric ? i + 1: ncols;

        for(unsigned int j = 0; j < maxj; j++ )
        {
            if( SPM_abs( ai[j] ) > zero_tol )
                row.nElements++;
        }
    }
    
    r = 0;
    for(unsigned int i = 0; i < nrows; i++)
    {
        SPM_SparseRow<matrixType> &row = rows[i];
        
        int ri = row.allocateColumns( row.nElements );
        
        if( ri != 0 )
            row.nElements = 0;
        
        r += ri;
        
        nElements += row.nElements;
    }
    
    if( r != 0 )
        return SPM_MEMORY_ERROR;
    
    
    for(unsigned int i = 0; i < nrows; i++)
    {		
        unsigned int maxj = symmetric ? i + 1: ncols;
        
        SPM_SparseRow<matrixType> &row = rows[i];
        
        const double *ai = A[i];
        
        unsigned int aux = 0; //save the index for the columns
        for(unsigned int j = 0; j < maxj; j++)
        {
            //if( ai[j] != 0.0 )
            if( SPM_abs( ai[j] ) > zero_tol )
            {
                row[aux].setColumn( j );
                row[aux].setValue( ai[j] );
                aux++;
            }
        }
    }
    
    return 0;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: setStructureAndValues(const matrixType *A, const double zero_tol, unsigned int nrows, unsigned int ncols)
{
    int r;
    
    
    if( nrows == 0 )
        nrows = this->nrows;
    
    if( ncols == 0 )
        ncols = this->ncols;
    
    
    if( nrows > this->nrows )
    {
        r = allocateSparseRows( nrows );
        if( r != 0 )
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
    }
    
    if( ncols > this->ncols )
    {
        this->ncols = ncols;
    }
    
    
    nElements = 0;
    
    //counting the nonzeros in the matriz
    for(unsigned int i = 0; i < nrows; i++)
    {
        rows[i].setNumberOfElements( 0 );
        
        unsigned int maxj = symmetric ? i + 1: ncols;
        
        const double *ai = &A[i*ncols];

        for(unsigned int j = 0; j < maxj; j++ )
        {
            if( SPM_abs( ai[j] ) > zero_tol )
                rows[i].nElements++;
        }
    }
    
    r = 0;
    for(unsigned int i = 0; i < nrows; i++)
    {
        SPM_SparseRow<matrixType> &row = rows[i];
        
        int ri = row.allocateColumns( row.nElements );
        
        if( ri != 0 )
            row.nElements = 0;
        
        r += ri;
        
        nElements += row.nElements;
    }
    
    if( r != 0 )
        return SPM_MEMORY_ERROR;
    
    
    
    for(unsigned int i = 0; i < nrows; i++)
    {		
        unsigned int maxj = symmetric ? i + 1: ncols;
        
        
        const double *ai = &A[i*ncols];
        
        unsigned int aux = 0; //save the index for the columns
        for(unsigned int j = 0; j < maxj; j++)
        {
            if( SPM_abs( ai[j] ) > zero_tol )
            {
                rows[i].columns[aux].col = j;
                rows[i].columns[aux].v = ai[j];
                aux++;
            }
        }
    }
    
    return 0;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::setSymmetricStructureAndValuesByLowerTriangle(const matrixType* lowerTA, const double zeroTolerance)
{
    unsigned int i, j, aux;
    spm::SPM_SparseElement<matrixType> *colAux;
    
    nElements = 0;
    
    for(i = 0; i < nrows; i++)
    {
        for(j = 0, aux = 0; j <= i; j++)
        {
            if( SPM_abs( lowerTA[ ((i + 1)*i)/2 + j] ) > zeroTolerance ) //nonzero element
            {
                aux++;
            }
        }
    
        nElements += aux;
        
        rows[i].nElements = aux;
        
        if(rows[i].columns)
        {
            free(rows[i].columns);
            rows[i].columns = NULL;
        }
    
        if( aux > 0 )
        {
            rows[i].columns = (spm::SPM_SparseElement<matrixType> *)malloc(aux *sizeof(spm::SPM_SparseElement<matrixType>));
            
            if( !rows[i].columns )
                return SPM_MEMORY_ERROR;
            
            
            colAux = rows[i].columns;
            for(j = 0, aux = 0 ; j <= i ; j++)
            {
                if( SPM_abs( lowerTA[ ((i + 1)*i)/2 + j] ) > zeroTolerance ) //nonzero element
                {
                    colAux[aux].col = j;
                    colAux[aux].v = lowerTA[ ((i + 1)*i)/2 + j];
                    aux++;
                }
            }
        }
    
    }
    
    symmetric = true;
    
    return 0;
}




template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::sumAllLines(const double* factors, matrixType* v) const
{
    spm::SPM_SparseElement<matrixType> *colAux;
    int aux, i, j;
    double f;
    
    for(i = 0; i < ncols; i++)
        v[i] = 0.0;
    
    
    for(i = 0; i < nrows; i++)
    {
        f = factors ? factors[i] : 1.0 ;
        
        aux = rows[i].nElements;
        colAux = rows[i].columns;
        
        for(j = 0; j < aux; j++)
            v[ colAux[j].col ] += f*colAux[j].v;
    }
    
}









template <class matrixType>
int spm::SPM_mergeSparseMatricesStructures(const spm::SPM_SparseMatrix<matrixType> &M1, const spm::SPM_SparseMatrix<matrixType> &M2, spm::SPM_SparseMatrix<matrixType> &Res )
{
    unsigned int k;
    int code;
    bool *cols = NULL;
    spm::SPM_SparseElement<matrixType> *colAux, *colAux2;


    Res.desallocateMemory(); //SPM_desalocateSparseMatrix(Res);

    code = Res.allocateSparseRows(SPM_max( M1.nrows, M2.nrows ) );
    if(code != 0)
    {
        goto desallocateMemory;
    }

    Res.ncols = SPM_max( M1.ncols, M2.ncols );

    cols = (bool*) malloc( Res.ncols * sizeof(bool) );
    if(!cols)
    {
        code = SPM_MEMORY_ERROR;
        goto desallocateMemory;
    }
    Res.nElements = 0;

    for(unsigned int i = 0; i < Res.nrows; i++ )
    {
        if( i < M1.nrows && M1.rows[i].nElements > 0 && i < M2.nrows && M2.rows[i].nElements > 0)
        {
            //M1 and M2 have elements in this line

            for(unsigned int j = 0; j < Res.ncols; j++)
                cols[j] = false;
            

            colAux = M1.rows[i].columns;
            for(unsigned int j = 0; j < M1.rows[i].nElements; j++)
                cols[ colAux[j].col ] = true;
            

            colAux = M2.rows[i].columns;
            for(unsigned int j = 0; j < M2.rows[i].nElements; j++)
                cols[ colAux[j].col ] = true;
            

            k = 0;
            for(unsigned int j = 0; j < Res.ncols; j++)
            {
                if( cols[j] )
                    k++;
            }

            if(k > 0)
            {
                unsigned int j;
                
                Res.rows[i].columns = (spm::SPM_SparseElement<matrixType> *) malloc( k * sizeof(spm::SPM_SparseElement<matrixType>) ) ;
                if( !Res.rows[i].columns )
                {
                    code = SPM_MEMORY_ERROR;
                    goto desallocateMemory;
                }
                Res.rows[i].nElements = k;

                colAux = Res.rows[i].columns;
                for(j = k = 0; j < Res.ncols; j++)
                {
                    if( cols[j] )
                    {
                        colAux[k].col = j;
                        colAux[k].v = 0.0;
                        k++;
                    }
                }

                Res.nElements += Res.rows[i].nElements;

                /*
                if( k == Res.rows[i].nElements )
                {
                    puts("Numero de elementos bateu!");
                }
                else
                {
                    puts("Numero de elementos nao bateu!");
                    getchar();
                }
                */
            }
        }
        else if( (i >= M1.nrows || M1.rows[i].nElements == 0) && (i < M2.nrows && M2.rows[i].nElements > 0) )
        {
            //Only M2 has elements in this line
            Res.rows[i].columns = (spm::SPM_SparseElement<matrixType> *) malloc( M2.rows[i].nElements * sizeof( spm::SPM_SparseElement<matrixType> ) );
            if( !Res.rows[i].columns )
            {
                code = SPM_MEMORY_ERROR;
                goto desallocateMemory;
            }
            colAux = M2.rows[i].columns;
            colAux2 = Res.rows[i].columns;

            Res.rows[i].nElements = M2.rows[i].nElements;
            Res.nElements += Res.rows[i].nElements;

            for(unsigned int j = 0; j < Res.rows[i].nElements; j++)
                colAux2[j].col = colAux[j].col;
            
        }
        else if( ( i < M1.nrows && M1.rows[i].nElements > 0 ) && ( i >= M2.nrows || M2.rows[i].nElements == 0) )
        {
            //only M1 has elements in this line
            Res.rows[i].columns = (spm::SPM_SparseElement<matrixType> *) malloc( M1.rows[i].nElements * sizeof(spm::SPM_SparseElement<matrixType>) );

            if( !Res.rows[i].columns )
            {
                code = SPM_MEMORY_ERROR;
                goto desallocateMemory;
            }

            colAux = M1.rows[i].columns;
            colAux2 = Res.rows[i].columns;

            Res.rows[i].nElements = M1.rows[i].nElements;
            Res.nElements += Res.rows[i].nElements;

            for(unsigned int j = 0; j < Res.rows[i].nElements; j++)
                colAux2[j].col = colAux[j].col;
            
        }
    }
    
    
    
    
    
    code = 0;

desallocateMemory:

    if( cols )		free(cols);
    if( code != 0 )
    {
        Res.desallocateMemory();
    }

    return code;
}















