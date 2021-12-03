
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <new>
#include "SPM_SparseMatrix.hpp"






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
	
	//if columns is NULL realloc acts as malloc...
	p = (spm::SPM_SparseElement<matrixType> *) realloc(columns, ncols * sizeof(spm::SPM_SparseElement<matrixType>) );
	
	if(!p)
	{
		return SPM_MEMORY_ERROR;
	}
	
	columns = p;
	
	this->nElements = ncols;
	return 0;
}


template <class matrixType>
int spm::SPM_SparseRow<matrixType>::copyStructureFrom(spm::SPM_SparseRow<matrixType> &other )
{
	int i;
	
	
	if( nElements != other.nElements )
	{
		desallocateColumns();
		
		if( other.nElements > 0 )
		{
			i = allocateColumns( other.nElements );
			if(i != 0)
				return i;
		}
	}
	
	for(i = 0; i < nElements; i++)
		columns[i].col = other.columns[i].col;
	
	return 0;
}


template <class matrixType>
int spm::SPM_SparseRow<matrixType>::copyFrom(spm::SPM_SparseRow<matrixType> &other )
{
	int i;
	
	
	if( nElements != other.nElements )
	{
		if( other.nElements == 0 )
		{
			desallocateColumns();
		}
		else
		{
			i = allocateColumns( other.nElements );
			if(i != 0)
				return i;
		}
	}
	
	for(i = 0; i < nElements; i++)
		columns[i] = other.columns[i];
	
	return 0;
}



template <class matrixType>
void spm::SPM_SparseRow<matrixType>::copyTo(double* line)
{
	unsigned int i;
	
	for(i = 0; i < (unsigned int) nElements; i++)
		line[ columns[i].col ] = columns[i].v;
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
double spm::SPM_SparseRow<matrixType>::evalTimesxt(const double* x)
{
	int k;
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
int  spm::SPM_SparseRow<matrixType>::getElement(unsigned int col, matrixType &value)
{
    int i;
    spm::SPM_SparseElement<matrixType> *colAux;


    colAux = columns;
    for(i = 0; i < nElements; i++)
    {
		if( colAux[i].col == col )
		{
		    value = colAux[i].v;
		    return 0;
		}
    }
    
    value = 0.0;
    
    //the element is not in the matrix
    return SPM_BAD_DEFINITIONS;
}



template <class matrixType>
template <class intClass>
unsigned int spm::SPM_SparseRow<matrixType>:: __getStructureAndValues( intClass* cols, matrixType* values)
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
unsigned int spm::SPM_SparseRow<matrixType>::getStructureAndValues( int* cols, matrixType* values)
{
	return __getStructureAndValues(cols, values);
}


template <class matrixType>
unsigned int spm::SPM_SparseRow<matrixType>::getStructureAndValues( unsigned int* cols, matrixType* values)
{
	return __getStructureAndValues(cols, values);
}





template <class matrixType>
template <class intClass>
unsigned int spm::SPM_SparseRow<matrixType>::__getStructure( intClass* cols)
{
	unsigned int k;
	
	for(k = 0; k < nElements; k++)
		cols[k] = columns[k].col;
	
	return k;
}


template <class matrixType>
unsigned int spm::SPM_SparseRow<matrixType>::getStructure( int* cols)
{
	return __getStructure(cols);
}



template <class matrixType>
unsigned int spm::SPM_SparseRow<matrixType>::getStructure( unsigned int* cols)
{
	return __getStructure(cols);
}




template <class matrixType>
unsigned int spm::SPM_SparseRow<matrixType>::getStructure( bool *cols)
{
	unsigned int k;
	
	
	for(k = 0; k < nElements; k++)
		cols[ columns[k].col ] = true;
	
	return k;
}



template <class matrixType>
unsigned int spm::SPM_SparseRow<matrixType>::getValues(matrixType *values)
{
	unsigned int k;
	
	for(k = 0; k < nElements; k++)
		values[k] = columns[k].v;
	
	return k;
}


template <class matrixType>
void spm::SPM_SparseRow<matrixType>:: setAllElements(matrixType value)
{
	for( unsigned int k = 0; k < nElements; k++ )
		columns[k].setValue( value ); //we could use (*this)[k] = value;
}


template <class matrixType>
int spm::SPM_SparseRow<matrixType>:: setElement(unsigned int col, matrixType value)
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
    getchar();
    SPM_desalocateSparseMatrix(*this); */
    
    desallocateMemory();
}




template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::accumulateLineInVector(const unsigned int line, const double factor, double* v)
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


//This functions assumes that the position is in the matrix. The current value in the position will be added to value parameter.
template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::addToElement(unsigned int row, unsigned int col, const matrixType value)
{
    int i;
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
		printf("Element was not found in the sparse matrix!\n");
		getchar();
    #endif
    
    //the element is not in the matrix
    return SPM_BAD_DEFINITIONS;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::addNewRows(const unsigned int nrows)
{
	int i;
	spm::SPM_SparseRow<matrixType> *p;
	
	p = (spm::SPM_SparseRow<matrixType> *) realloc( rows, (this->nrows + nrows) * sizeof(spm::SPM_SparseRow<matrixType>) );
	
	if(!p)
		return SPM_MEMORY_ERROR;
	
	rows = p;
	
	for(i = this->nrows; i < this->nrows + nrows; i++)
		rows[i].initialize();
	
	
	this->nrows += nrows;
	
	return 0;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::allocateSparseRows(const unsigned int m)
{
    int i;
	
	//do not use new because we can need realloc
	
	if(rows)
		free(rows);
	
    rows = (spm::SPM_SparseRow<matrixType>*) malloc( m*sizeof(spm::SPM_SparseRow<matrixType>) );
    if( !(rows) )
		return SPM_MEMORY_ERROR;

    nrows = m;
    for(i = 0; i < m; i++)
		rows[i].initialize();
	

    return 0;
}


template <class matrixType>
bool spm::SPM_SparseMatrix<matrixType>::compareSparseMatrixStructure(const spm::SPM_SparseMatrix<matrixType> &other)
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
			    //the order os the columns can be diferent...
			    for(k = 0; k < rows[i].nElements; k++)
			    {
					if( colAux[j].col == colAux2[j].col )
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
int spm::SPM_SparseMatrix<matrixType>::copyStructureFrom(SPM_SparseMatrix& M)
{
	int i, j;
    //spm::SPM_SparseElement<matrixType> *colAux, *colAuxOther;
	
	
	if( nrows != M.nrows )
	{
		desallocateMemory();
		
		i = allocateSparseRows(M.nrows);
		if(i != 0)
			return i;
	}
	
	copyParametersFrom(M);
	
	
	for(i = 0; i < nrows; i++)
    {
		j = rows[i].copyStructureFrom( M.rows[i] );
		if( j != 0 )
			return SPM_MEMORY_ERROR;
    }
    
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

    for(unsigned int i = 0; i < nrows; i++)
    {
		r = rows[i].copyFrom( M.rows[i] );
		if( r != 0 )
		{
			code = SPM_MEMORY_ERROR;
			goto termination;
		}
    }
    
    
    
    code = 0;
	
termination:
    
    return code;
}



template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>:: copyMatrixTo(matrixType** M, const bool initializeWithZero)
{
	unsigned int i, j;
	matrixType* row;
	
	if( initializeWithZero )
	{
		for(i = 0; i < nrows; i++)
		{
			row = M[i];
			for(j = 0; j < ncols; j++)
				row[j] = 0.0;
		}
	}
	
	
	for(i = 0; i < nrows; i++)
		rows[i].copyTo( M[i] );
	
	
	if( symmetric )
	{
		for(i = 1; i < nrows; i++)
		{
			row = M[i];
			for( j = 0; j < i; j++ )
			{
				M[j][i] = row[j];
			}
		}
	}
	
}



template <class matrixType>
template <class intClass>
void
spm::SPM_SparseMatrix<matrixType>::__countRowsEachColumn(intClass *counts, const bool accumulate)
{
	//spm::SPM_SparseElement<matrixType> *colAux;
	unsigned int i, j, aux;
	
	
	if(!accumulate)
	{
		for(i = 0; i < ncols; i++)
			counts[i] = 0.0;
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
void spm::SPM_SparseMatrix<matrixType>::countRowsEachColumn( int *counts, const bool accumulate)
{
	return __countRowsEachColumn(counts, accumulate);
}



template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::countRowsEachColumn(unsigned int *counts, const bool accumulate)
{
	return __countRowsEachColumn(counts, accumulate);
}


template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>:: deleteStructure()
{
    int i;
    for(i = 0; i < nrows; i++)
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
    int i;

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
int  spm::SPM_SparseMatrix<matrixType>::getElement(unsigned int row, unsigned int col, matrixType& value)
{
    unsigned int i, aux;
    //spm::SPM_SparseElement<matrixType> *colAux;
	spm::SPM_SparseRow<matrixType> &myrow = rows[row];

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
int spm::SPM_SparseMatrix<matrixType>::getFullRow(const unsigned int row, matrixType* values)
{
	unsigned int j, aux;
	//spm::SPM_SparseElement<matrixType> *colAux;
	spm::SPM_SparseRow<matrixType> &myrow = rows[row];
	
	for(j = 0; j < ncols; j++)
		values[j] = 0.0;
	
	
	aux =  myrow.getNumberOfElements(); //rows[row].nElements;
	//colAux = rows[row].columns;
	
	for(j = 0; j < aux; j++)
		values[ myrow[j].getColumn() ] = myrow[j].getValue();
		//values[ colAux[j].col ] = colAux[j].v;
	
	return 0;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getFullRowAccumulation( const unsigned int row, matrixType* values)
{
	unsigned int j, aux;
	spm::SPM_SparseRow<matrixType> &myrow = rows[row];
	//spm::SPM_SparseElement<matrixType> *colAux;
	
	aux = myrow.getNumberOfElements(); // rows[line].nElements; colAux = rows[line].columns;
	
	for(j = 0; j < aux; j++)
		values[ myrow[j].getColumn() ] += myrow[j].getValue(); //colAux[j].v;
		//values[ colAux[j].col ] += colAux[j].v;
	
	return 0;
}


template <class matrixType>
matrixType * spm::SPM_SparseMatrix<matrixType>:: copyMatrixTo(matrixType *Matrix )
{
    int i, j, k, aux;
    matrixType *M;
    spm::SPM_SparseElement<matrixType> *colAux;
    
	if(Matrix)
	{
		M = Matrix;
	}
	else
	{
	    M = (matrixType *) malloc( nrows * ncols * sizeof(matrixType) );
	    if( !M )
	    {
			return NULL;
	    }
	}
	
    
    aux = 0;
    for(i = 0; i < nrows; i++)
		for(j = 0; j < ncols; j++, aux++)
		    M[ aux ] = 0.0;
    
    for(i = 0; i < nrows; i++ )
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
    }
    
    return M;
}






template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowStructureAndValues(const unsigned int line, int* jCol, matrixType* vals)
{
	return rows[line].getStructureAndValues(jCol, vals);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowStructureAndValues(const unsigned int line, unsigned int* jCol, matrixType* vals)
{
	return rows[line].getStructureAndValues(jCol, vals);
}



template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowStructure(const unsigned int line, unsigned int* jCol)
{
	return rows[line].getStructure(jCol);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowStructure(const unsigned int line, int* jCol)
{
	return rows[line].getStructure(jCol);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowStructure(const unsigned int line, bool* cols, bool accumulate)
{
	if(!accumulate)
	{
		for(int i = 0; i < ncols; i++)
			cols[i] = false;
	}
	
	return rows[line].getStructure(cols);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getRowValues(const unsigned int line, double* vals)
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
int spm::SPM_SparseMatrix<matrixType>:: __getStructure(intType *iRow, intType *jCol)
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
int spm::SPM_SparseMatrix<matrixType>:: getStructure(int *iRow, int *jCol)
{
	return __getStructure( iRow, jCol );
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: getStructure(unsigned int *iRow, unsigned int *jCol)
{
	return __getStructure( iRow, jCol );
}



template <class matrixType>
template <class intType>
int spm::SPM_SparseMatrix<matrixType>:: __getStructureAndValues(intType* iRow, intType* jCol, matrixType* vals)
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
int spm::SPM_SparseMatrix<matrixType>:: getStructureAndValues(int* iRow, int* jCol, matrixType* vals)
{
	return __getStructureAndValues(iRow, jCol, vals);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>:: getStructureAndValues(unsigned int* iRow, unsigned int* jCol, matrixType* vals)
{
	return __getStructureAndValues(iRow, jCol, vals);
}


/*
 * If any of the vectors row, cols or values is NULL, we allocate memory for the one. Otherwise, we ASSUME that a enough quantity of memory is already allocated.
 */
template <class matrixType>
template <class intType>
int spm::SPM_SparseMatrix<matrixType>::__getTripleSparseFormat(intType** rows, intType** cols, matrixType** values)
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
int spm::SPM_SparseMatrix<matrixType>::getTripleSparseFormat(int** rows, int** cols, matrixType** values)
{
	return __getTripleSparseFormat( rows, cols, values );
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getTripleSparseFormat(unsigned int** rows, unsigned int** cols, matrixType** values)
{
	return __getTripleSparseFormat( rows, cols, values );
}






template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::getValues(matrixType *values)
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
bool spm::SPM_SparseMatrix<matrixType>::hasIndex(const unsigned int line, const unsigned int col)
{
	const unsigned int nEls = rows[line].nElements;
	unsigned int i;
	const spm::SPM_SparseElement<matrixType> *cols = rows[line].columns;
	
	
	for(i = 0; i < nEls; i++)
	{
		if( cols[i].col == col )
			return true;
	}
	
	return false;
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
double spm::SPM_SparseMatrix<matrixType>::lineEvaluation(const unsigned int line, const matrixType* x)
{
	return rows[line].evalTimesxt(x);
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
void spm::SPM_SparseMatrix<matrixType>::printAllMatrix(void)
{
    int i, j;
    double c;
    
    for(i = 0; i < nrows; i++)
    {
		for(j = 0; j < ncols; j++)
		{
		    getElement(i, j, c);
		    printf("%0.6f, ", c);
		}
		printf("; ");
    }
}


template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::printSparseMatrix(const bool showEmptyLines)
{
    int i, j;

    for(i = 0; i < nrows; i++)
    {
		if( rows[i].nElements > 0 || showEmptyLines)
		{
			printf("%d=> ", i);
		
			
			for(j = 0; j < rows[i].nElements; j++)
			{
			    printf("%d: %f   ", rows[i].columns[j].col, rows[i].columns[j].v);
			}
			printf("\n");
		}
		else if( showEmptyLines )
		{
			printf("%d=> \n", i);
		}
    }
}



template <class matrixType>
void spm::SPM_SparseMatrix<matrixType>::multiplyAllElements(const matrixType value)
{
	int i, j, aux;
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
double spm::SPM_SparseMatrix<matrixType>::quadraticEvaluation(const matrixType *x)
{
    int i, j, k, aux;
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
void spm::SPM_SparseMatrix<matrixType>::quadraticGradientEvaluation( const matrixType* x, matrixType* grad, const bool accumulate)
{
	int i, j, k, aux;
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
template <class intType>
int spm::SPM_SparseMatrix<matrixType>::__removeRows(const unsigned int nlines, const intType* lines)
{
	bool *remove = NULL;
	int i, code;
	
	
	remove = (bool *) calloc( nrows, sizeof(bool) ) ;
	if(!remove)
	{
		code = SPM_MEMORY_ERROR;
		goto desallocate_memory;
	}
	
	
	for(i = 0; i < nlines; i++)
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
	int i, k, code;
	spm::SPM_SparseRow<matrixType> *aux;
	
	for(i = nrows-1; i >= 0; i--)
	{
		if( lines[i] )
		{
			for(k = i + 1; k < nrows; k++)
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
	
	
	symmetric = false;
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
int spm::SPM_SparseMatrix<matrixType>::setElement(unsigned int row, unsigned int col, const matrixType value)
{
    int r;

    if( symmetric )
    {
		//we work only on the lower triangle
		if( col > row )
		    SPM_swap(col, row);  //i = col;   col = row;    row = i;
    }
    
    
    r = rows[row].setElement(col, value);
    
    /*colAux = rows[row].columns;
    for(i = 0; i < rows[row].nElements; i++)
    {
		if( colAux[i].col == col )
		{
		    colAux[i].v = value;
		    return 0;
		}
    } */
    
	#if SPM_DEBUG_MODE
	if(r != 0)
		printf("Warning: element not found on spm::SPM_SparseMatrix<matrixType>::setElement. row: %d col: %d\n", row, col);
	#endif
	
	
    //the element is not in the matrix
    return r;
}



template <class matrixType>
template <class intType>
//warning: that method store a single line. If that line  already exists, it will be overwritten
int spm::SPM_SparseMatrix<matrixType>::__setRowStructure(  const unsigned int row, const unsigned int numberOfNZElements, const intType* cols)
{
    unsigned int i = 0;
	int ret, code;
	
	
    if( row >= nrows )
		return SPM_BAD_DEFINITIONS; //do not use goto here
	
	spm::SPM_SparseRow<matrixType> &myrow = rows[row];
	
	nElements -= myrow.nElements;
	
	if( myrow.nElements != numberOfNZElements  )
		myrow.desallocateColumns();
	
	
    if( numberOfNZElements > 0 )
    {
		/*rows[constraint].columns = (spm::SPM_SparseElement<matrixType> *) malloc( rows[constraint].nElements * sizeof(spm::SPM_SparseElement<matrixType>) ) ;

		if( !rows[constraint].columns )
		{
		    rows[constraint].nElements = 0;
		    code = SPM_MEMORY_ERROR;
		    goto desallocateMemory;
		} */
		
		
		ret = myrow.allocateColumns( numberOfNZElements );
		if( ret != 0 )
		{
			code = ret;
			goto termination;
		}
		
		//colAux = rows[constraint].columns;
		
		
		for(i = 0; i < numberOfNZElements; i++)
		{
			myrow[i].setColumn( cols[i] );
			myrow[i].setValue( 0.0 );
		}
    }
    else
    {
		myrow.nElements = 0;
    }
    
    
    
    
    #if SPM_DEBUG_MODE
		if( i != rows[row].nElements )
		{
		    puts("Para tudo! Inconsistencia na contagem de elementos na definicao de uma unica restricao linear");
		    printf("Please report this mensage and this code for %s\n", SPM_EMAIL);
		    getchar();
		}
    #endif

    code = 0;

termination:

    nElements += rows[row].nElements;

    return code;
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
int spm::SPM_SparseMatrix<matrixType>::setRowStructure(const unsigned int line, SPM_SparseRow< matrixType >& row)
{
	int i, code;
	//spm::SPM_SparseElement<matrixType> *colAux, *colAux2;
	
	if(line < 0 || line >= nrows)
		return SPM_BAD_DEFINITIONS; //do not use goto here
	
	
	nElements -= rows[line].nElements;
	
	i = rows[line].copyStructureFrom( row );
	
	if( i != 0 )
	{
		code = i;
		goto termination;
	}
	
	
	/* rows[line].columns = (spm::SPM_SparseElement<matrixType> *) malloc(  nEl * sizeof(spm::SPM_SparseElement<matrixType>) ) ;
	
	if( !rows[line].columns )
	{
		return SPM_MEMORY_ERROR;
	}
	
	
	colAux = rows[line].columns;
	colAux2 = row.columns;
	
	for(i = 0; i < nEl; i++)
		colAux[i] = colAux2[i];
	
	rows[line].nElements = nEl; */
	
	
	
	code = 0;
	
termination:
	
	nElements += rows[line].nElements;
	
	return code;
}




template <class matrixType>
template <class intType>
//warning: that method store a single linear constraint. If that constraints have already exists, it will be overwritten
int spm::SPM_SparseMatrix<matrixType>::__setRowStructureAndValues(const unsigned int row, const unsigned int numberOfNZElements, const intType* cols, const matrixType* vals)
{
    int code;
	unsigned int aux = 0;

    if(row < 0 || row >= nrows)
		return SPM_BAD_DEFINITIONS; //do not use goto here
    
    nElements -= rows[row].nElements;
	
    rows[row].desallocateColumns();

    if( numberOfNZElements > 0 )
	{
		spm::SPM_SparseRow<matrixType> &myrow = rows[row];
		
		code = myrow.allocateColumns( numberOfNZElements );
		
		if( code != 0 )
		{
		    code = SPM_MEMORY_ERROR;
		    goto desallocateMemory;
		}
		

		for(unsigned int i = 0; i < numberOfNZElements; i++)
		{
			myrow[i].setColumn( cols[i] );
			myrow[i].setValue( vals[i] );
			
		    //rows[row].columns[aux].col = cols[i];
		    //rows[row].columns[aux].v = vals[i];
		    aux++;
		}
    }
    
    
    

    #if SPM_DEBUG_MODE
		if( aux != rows[row].nElements )
		{
		    puts("Para tudo! Inconsistencia na contagem de elementos na definicao de uma unica restricao linear");
		    printf("Please report this mensage and this code for %s\n", SPM_EMAIL);
		    getchar();
		}
    #endif

    code = 0;

desallocateMemory:

    nElements += rows[row].nElements;

    return code;
}



template <class matrixType>
//warning: that method store a single linear constraint. If that constraints have already exists, it will be overwritten
int spm::SPM_SparseMatrix<matrixType>::setRowStructureAndValues(const unsigned int row, const unsigned int numberOfNZElements, const int* cols, const matrixType* vals)
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
//warning: that method store a single linea. If the line have already exists, it will be overwritten. If the matrix is symmetric, just the elements in lower triangle elements are setted
int spm::SPM_SparseMatrix<matrixType>::setRowStructureAndValues(const unsigned int row, const matrixType* a)
{
    int code;
	unsigned int ncols = symmetric ? row + 1 : this->ncols;
	unsigned int aux = 0, realncols = 0;
	

    if(row >= nrows || row < 0)
		return SPM_BAD_DEFINITIONS; //do not use goto here
	
	spm::SPM_SparseRow<matrixType> &myrow = rows[row];
	
	
	nElements -= myrow.nElements;
	myrow.desallocateColumns();
	
	
    for( unsigned int i = 0; i < ncols; i++ )
    {
		if( a[i] != 0.0 )
		    realncols++;
    }

    if( realncols > 0 )
    {
		code = myrow.allocateColumns( realncols );

		if( code )
		{
		    code = SPM_MEMORY_ERROR;
		    goto desallocateMemory;
		}

		for(unsigned int i = 0; i < ncols; i++)
		{
		    if( a[i] != 0.0 )
		    {
				myrow[aux].setColumn( i );
				myrow[aux].setValue( a[i] );
				//rows[row].columns[aux].col = i;
				//rows[row].columns[aux].v = a[i];
				aux++;
		    }
		}
    }
    
    
    #if SPM_DEBUG_MODE
		if( aux != rows[row].nElements )
		{
		    puts("Para tudo! Inconsistencia na contagem de elementos na definicao de uma unica restricao linear");
		    printf("\nPlease report this mensage and this code for %s\n", SPM_EMAIL);
		    getchar();
		}
    #endif

    code = 0;

desallocateMemory:

    nElements += rows[row].nElements;

    return code;
}



template <class matrixType>
template <class intType>
int spm::SPM_SparseMatrix<matrixType>::__setStructureAndValues(const unsigned int numberOfNZElements, const intType* iRow, const intType* jCol, const matrixType* vals)
{
    //bool found;
    int i, j, code, row, col;
    int *aux, aux3;
    spm::SPM_SparseElement<matrixType> *colAux;

    //printf("nrows: %d\n", nrows);
    aux = (int*) calloc( (nrows), sizeof(int) );
    if( !aux )
    {
		code = SPM_MEMORY_ERROR;
		goto desallocateMemory;
    }

    //for(i = 0; i < nrows; i++)
	//aux[i] = 0;

    for(i = 0; i < numberOfNZElements; i++)
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
		    
		    //aux[ iRow[i] ]++;
		}
	
		//checking if there is already that element
		//found = false;
		
		colAux = rows[ row ].columns;
		aux3 = rows[ row ].nElements;
	
	
		for(j = 0; j < aux3; j++)
		{
		    if(colAux[j].col == col)
		    {
				//we have already this element in sparse matrix
				code = SPM_BAD_DEFINITIONS;
				goto desallocateMemory;
		    }
		}
	
		aux[row]++;
    }
    
    
    //now, we allocating space for the new elements
    for(i = 0; i < nrows; i++)
    {
		if( aux[i] == 0 )
		    continue;

		if( rows[i].nElements == 0 )
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
		}
    }


    //storing the new elements. Here, aux has the first empty index to use in each row.
    for(i = 0; i < numberOfNZElements; i++)
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
		    //rows[ iRow[i] ].columns[ aux[ iRow[i] ] ].col = jCol[i];
		    //rows[ iRow[i] ].columns[ aux[ iRow[i] ] ].v = 0.0;
		    //aux[ iRow[i] ]++;
		}
	
		rows[ row ].columns[ aux[ row ] ].col = col;
		rows[ row ].columns[ aux[ row ] ].v =  vals ? vals[i] : 0.0;
		aux[ row ]++;
    }

    nElements += numberOfNZElements;
    //nzEleConstrFGrad = M.nElements;

    
    #if SPM_DEBUG_MODE
		for(i = 0; i < nrows; i++)
		{
		    if( rows[i].nElements != aux[i] && aux[i] > 0 )
		    {
				puts("Para tudo! Inconsistencia no numero de elementos setados para a matriz esparsa das restricoes lineares.");
				
				printf("i: %d rows[i].nElements: %d aux[i]: %d\n", i, rows[i].nElements, aux[i]);
				
				printf("Please report this code and this message for %s\n", SPM_EMAIL);
				getchar();
		    }
		}
    #endif

    code = 0;

desallocateMemory:

    if(aux)		free(aux);

    return code;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::setStructureAndValues(const unsigned int numberOfNZElements, const int* iRow, const int* jCol, const matrixType* vals)
{
	return __setStructureAndValues(numberOfNZElements, iRow, jCol, vals);
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::setStructureAndValues(const unsigned int numberOfNZElements, const unsigned int* iRow, const unsigned int* jCol, const matrixType* vals)
{
	return __setStructureAndValues(numberOfNZElements, iRow, jCol, vals);
}





template <class matrixType>
//warning: that method overwrittes possible lines already alllocated
int spm::SPM_SparseMatrix<matrixType>::setStructureAndValues(matrixType** A)
{
    int i, j, maxj;
    int aux;
    
    
    nElements = 0;

    //counting the nonzeros in the matriz
    for(i = 0; i < nrows; i++)
    {
		rows[i].nElements = 0;
		
		maxj = symmetric ? i + 1: ncols;

		for(j = 0; j < maxj; j++ )
		{
		    if( A[i][j] != 0.0 )
		    {
				rows[i].nElements++;
		    }
		}
	}
		
		
	for(i = 0; i < nrows; i++)
    {	
		nElements += rows[i].nElements;

		if( rows[i].columns != NULL )
		{
		    free(rows[i].columns);
		    rows[i].columns = NULL;
		}

		if( rows[i].nElements > 0)
		{
		    rows[i].columns = (spm::SPM_SparseElement<matrixType> *) malloc( rows[i].nElements * sizeof(spm::SPM_SparseElement<matrixType>) );

		    if( !rows[i].columns )
				return SPM_MEMORY_ERROR;
		}
	}
	
	for(i = 0; i < nrows; i++)
    {		
		maxj = symmetric ? i + 1: ncols;
		
	    aux = 0; //save the index for the columns
	    for(j = 0; j < maxj; j++)
	    {
			if( A[i][j] != 0.0 )
			{
			    rows[i].columns[aux].col = j;
			    rows[i].columns[aux].v = A[i][j];
			    aux++;
			}
	    }
    }
    
    return 0;
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::setStructureAndValues(const matrixType *A)
{
	int i, j, maxj;
	int aux;
	
	nElements = 0;
	
	//counting the nonzeros in the matriz
    for(i = 0; i < nrows; i++)
    {
		rows[i].nElements = 0;
		
		maxj = symmetric ? i + 1: ncols;

		for(j = 0; j < maxj; j++ )
		{
		    if( A[i*ncols + j] != 0.0 )
		    {
				rows[i].nElements++;
		    }
		}
	}
	
	
	for(i = 0; i < nrows; i++)
    {	
		nElements += rows[i].nElements;

		if( rows[i].columns != NULL )
		{
		    free(rows[i].columns);
		    rows[i].columns = NULL;
		}

		if( rows[i].nElements > 0)
		{
		    rows[i].columns = (spm::SPM_SparseElement<matrixType> *) malloc( rows[i].nElements * sizeof(spm::SPM_SparseElement<matrixType>) );

		    if( !rows[i].columns )
				return SPM_MEMORY_ERROR;
		}
	}
	
	
	for(i = 0; i < nrows; i++)
    {		
		maxj = symmetric ? i + 1: ncols;
		
	    aux = 0; //save the index for the columns
	    for(j = 0; j < maxj; j++)
	    {
			if( A[i*ncols + j] != 0.0 )
			{
			    rows[i].columns[aux].col = j;
			    rows[i].columns[aux].v = A[i*ncols + j];
			    aux++;
			}
	    }
    }
    
    return 0;
	
}


template <class matrixType>
int spm::SPM_SparseMatrix<matrixType>::setSymmetricStructureAndValuesByLowerTriangle(const matrixType* lowerTA, const double zeroTolerance)
{
    int i, j, aux;
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
void spm::SPM_SparseMatrix<matrixType>::sumAllLines(const double* factors, matrixType* v)
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
    int i, j, k, code;
    bool *cols = NULL;
    spm::SPM_SparseElement<matrixType> *colAux, *colAux2;


    Res.desallocateMemory(); //SPM_desalocateSparseMatrix(Res);

    j = Res.allocateSparseRows(SPM_max( M1.nrows, M2.nrows ) );
    if(j != 0)
    {
		code = j;
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

    for( i = 0; i < Res.nrows; i++ )
    {
		if( i < M1.nrows && M1.rows[i].nElements > 0 && i < M2.nrows && M2.rows[i].nElements > 0)
		{
		    //M1 and M2 have elements in this line

		    for(j = 0; j < Res.ncols; j++)
		    {
				cols[j] = false;
		    }

		    colAux = M1.rows[i].columns;
		    for(j = 0; j < M1.rows[i].nElements; j++)
		    {
				cols[ colAux[j].col ] = true;
		    }

		    colAux = M2.rows[i].columns;
		    for(j = 0; j < M2.rows[i].nElements; j++)
		    {
				cols[ colAux[j].col ] = true;
		    }

		    k = 0;
		    for(j = 0; j < Res.ncols; j++)
		    {
				if( cols[j] )
					k++;
		    }

		    if(k > 0)
		    {
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

		    for(j = 0; j < Res.rows[i].nElements; j++)
		    {
				colAux2[j].col = colAux[j].col;
		    }
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

			for(j = 0; j < Res.rows[i].nElements; j++)
			{
				colAux2[j].col = colAux[j].col;
			}
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















