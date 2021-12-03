
/***********************************************************
 ***********************************************************
 ****                                                   ****
 ****          my sparse matrix implementation          ****
 ****                                                   ****
 ***********************************************************
 ***********************************************************/

#ifndef SPM_SPARSE_MATRIX_HPP
#define SPM_SPARSE_MATRIX_HPP

#include <cstdlib>
#include <math.h>

namespace spm
{

#define SPM_DEBUG_MODE 1


#define SPM_MEMORY_ERROR -101
#define SPM_BAD_DEFINITIONS -102
#define SPM_REPETEAD_ELEMENT -103
#define SPM_BAD_VALUE -104
#define SPM_UPPER_TRIANGLE -105
#define SPM_INDEX_FAULT -106

#define SPM_EMAIL "wendelmelo@cos.ufrj.br"






template <class myClass>
inline myClass SPM_min( const myClass a, const myClass b )
{
    return (a <= b ? a : b);
}


template <class myClass>
inline myClass SPM_max( const myClass a, const myClass b )
{
    return (a >= b ? a : b);
}


template <class myClass>
inline myClass SPM_abs( const myClass num )
{
    return ( num >= 0 ? num : -num );
}


template <class myClass1, class myClass2>
inline void SPM_swap( myClass1 &a, myClass2 &b )
{
	myClass1 ca = a;
	
	a = b;
	b = ca;
}


template <class myClass>
inline void SPM_setAllArray( const unsigned int size, myClass *a, myClass value )
{
	#pragma ivdep
	#pragma GCC ivdep
	for( unsigned int i = 0; i < size; i++ )
		a[i] = value;
}



template <class baseClass>
class SPM_Iterator
{
public:
	
	unsigned int pos;
	baseClass *obj;
	
	inline unsigned int getPosition()
	{
		return pos;
	}
	
	
	inline void setPosition(unsigned int pos)
	{
		this->pos = pos;
	}
	
	
	inline void initialize( baseClass *basisObject, unsigned int initialPosition = 0)
	{
		setPosition( initialPosition );
		this->obj = basisObject;
	}
	
	
	
	SPM_Iterator( baseClass *basisObject = NULL, unsigned int initialPosition = 0)
	{
		initialize(basisObject, initialPosition);
	}
	
	
	inline bool operator== (const SPM_Iterator<baseClass> &iter)
	{
		return iter.obj == obj && iter.pos == pos;
	}
	
	
	inline SPM_Iterator<baseClass> begin(  )
	{
		return SPM_Iterator<baseClass>(obj);
	}
	
	
	inline SPM_Iterator<baseClass> end(  )
	{
		return SPM_Iterator<baseClass>(obj, obj->iteratorSize() );
	}
	
	
	inline bool isAtBegin()
	{
		return pos == 0;
	}
	
	
	inline bool isAtEnd()
	{
		return pos >= obj->iteratorSize();
	}
	
	
	inline void goToBegin()
	{
		setPosition( 0 );
	}
	
	inline void goToEnd()
	{
		setPosition( obj->iteratorSize() );
	}
	
	
	inline void operator++(int arg) //prefix version has no input parameter. posfix version receives 
	{
		pos++;
	}
	
	inline void operator--(int arg)
	{
		pos--;
	}
	
	
	
	inline void operator+=(const unsigned int n)
	{
		pos += n;
	}
	
	inline void operator-=(const unsigned int n)
	{
		pos -= n;
	}
	
};





/* Small class to represent a value in a sparse matrix. Each line of sparse matrix will be a collection (maybe an array) of objects from this class*/

template <class matrixType>
class SPM_SparseElement
{
public:

    unsigned int col;	//column of the element
    matrixType v;	//value of element

    //SPM_SparseElement(){} //void constructor
    
    inline unsigned int getColumn() const
	{
		return col;
	}
	
	
	inline matrixType getValue() const
	{
		return v;
	}
	
	
	inline void setColumn(const unsigned int column)
	{
		col = column;
	}
	
	
	inline void setValue(const matrixType value)
	{
		v = value;
	}
	
	
	//operator "=" just change the value
	inline SPM_SparseElement& operator =( const matrixType newValue )
	{
		v = newValue;
		return *this;
	}
	
};


template <class matrixType>
class SPM_SparseMatrix;


template <class matrixType>
class SPM_SparseRow
{
protected:
	
	//those methods are protected because if user call them directly, the counting of element in sparse matrix will be wrong...
	
	
	int allocateColumns(const unsigned int ncols);
	
	void desallocateColumns();
	
	int reallocateColumns( const unsigned int ncols );
	
	
	
	
	
public:
	
    SPM_SparseElement<matrixType>  *columns; //pointer to elements of the row
    unsigned int nElements;
	
	
	//we do not have methods to set the strucutre because SPM_SparseMatrix must have the control of total number of elements...
	
    SPM_SparseRow();
	
	~SPM_SparseRow();
	
	
	
	
	int copyStructureFrom(SPM_SparseRow<matrixType> &other );
	
	int copyFrom(SPM_SparseRow<matrixType> &other );
	
	void copyValuesInOrderFrom( SPM_SparseRow<matrixType> &other  );
	
	void copyTo(double* line, const double factor = 1.0) const;
	
	
	
    double evalTimesxt(const double *x) const;
    
    void initialize();
	
	inline unsigned int getNumberOfElements() const
	{
		return nElements;
	}
	
	int getElement(const unsigned int col, matrixType &value) const;
    
	
	//that method return the number of elements in the row. You must initialize by yourself all positions in cols as false 
    unsigned int getStructure(bool *cols) const;
    
	//those methods return the number of elements in the row
    unsigned int getStructure(int *cols) const;
	unsigned int getStructure(unsigned int *cols) const;
	
	//those methods return the number of elements in the row
    unsigned int getStructureAndValues(int *cols, matrixType *values) const;
	unsigned int getStructureAndValues(unsigned int *cols, matrixType *values) const;
	
	unsigned int getValues(matrixType *values) const;
	
	//this fucntion returns on index, the index position in columns where the col is found...
	bool hasColumn(const unsigned int col, unsigned int *index = NULL) const;
	
	void multiplyAllElements( const matrixType value );
	
	void print() const;
	
	void setAllElements(const matrixType value);
	
	int setElement(const unsigned int col, const matrixType value);
	
	
	//set first 'nel' elemets in this line using the first 'nel' values in array values. Note, this method does not look to columns indexes...
	int setElementsByOrder(const unsigned int nel, const matrixType *values);
	
	
	inline void setNumberOfElements(unsigned int nEls)
	{
		nElements = nEls;
	}
	
	
	inline SPM_SparseElement<matrixType>& operator[] (unsigned int index)
	{
		return columns[index];
	}
	
	
	class iterator: public SPM_Iterator< SPM_SparseRow<matrixType> >
	{
		
	public:
		
		iterator( SPM_SparseRow<matrixType> *row = NULL, unsigned int initialPosition = 0 ) : SPM_Iterator< SPM_SparseRow<matrixType> >( row, initialPosition )
		{
		}
		
		
		inline unsigned int getCurrentColumn()
		{
			const unsigned int &pos = SPM_Iterator< SPM_SparseRow<matrixType> >::pos;
			
			return SPM_Iterator< SPM_SparseRow<matrixType> >::obj->columns[pos].col;
		}
		
		inline matrixType getCurrentValue()
		{
			const unsigned int pos = SPM_Iterator< SPM_SparseRow<matrixType> >::pos;
			
			return SPM_Iterator< SPM_SparseRow<matrixType> >::obj->columns[pos].v;
		}
		
		
		inline SPM_SparseElement<matrixType>& operator *()
		{
			const unsigned int pos = SPM_Iterator< SPM_SparseRow<matrixType> >::pos;
			
			return SPM_Iterator< SPM_SparseRow<matrixType> >::obj->columns[pos];
		}
		
	};
	
	
	
protected:
	
	inline unsigned int iteratorSize() const
	{
		return nElements;
	}
	
	
	template <class intClass>
	unsigned int __getStructure(intClass *cols) const;
	
	
	template <class intClass>
	unsigned int __getStructureAndValues(intClass *cols, matrixType *values) const;
	
	
	friend class SPM_Iterator< SPM_SparseRow<matrixType> >;
	
	
	friend class SPM_SparseMatrix< matrixType >;
	
};







template <class matrixType>
int SPM_mergeSparseMatricesStructures(const SPM_SparseMatrix<matrixType> &M1, const SPM_SparseMatrix<matrixType> &M2, SPM_SparseMatrix<matrixType> &Res );





template <class matrixType>
class SPM_SparseMatrix
{
	void copyParametersFrom(const SPM_SparseMatrix& other);
	
	
public:
    bool symmetric;
    //bool constant;
    //bool constant_calculated;	//only used by constant matriz...
    unsigned int nrows, ncols;		//Number of rows and columns
    unsigned int nElements;		//Number of nElements
    SPM_SparseRow<matrixType> *rows;
	
	
public:
    
    SPM_SparseMatrix();
    

    SPM_SparseMatrix(const unsigned int nrows, const unsigned int ncols, const bool symmetric = false);
    
	
	//that method accumulates the sum between this matrix and factor*A. All positions in A should be also in this matrix. 
    void accumulatedSum(const double factor, SPM_SparseMatrix& A);
	
	
	//this method accumulates in M the current sparse matrix times factor ( M = M + facotr*this). Note, the result is stored in M and the sparse matrix is NOT changed...
	void accumulatedSumToMatrix( const double factor, matrixType *M, const bool considerSymmetry = true) const;
	
	//that method sum the coefficients in a line times a factor in a vector. The vector is not initialized
	void accumulateLineInVector(const unsigned int line, const double factor, double *v) const;
    
	//This functions assumes that the position is in the matrix. The current value in the position will be added to value parameter.
	int addToElement(unsigned int row, unsigned int col, const matrixType value);
	
	int addNewRows(const unsigned int nrows);
	
    int allocateSparseRows(const unsigned int m);
	
	
	int checkRowStructure( unsigned int row, unsigned int nzs, int *cols, double *values = NULL) const;
	
	int checkRowStructure( unsigned int row, unsigned int nzs, unsigned int *cols, double *values = NULL) const;
	
	int checkTripleSparseStructure( unsigned int nzs, const int *rows, const int *cols, const double *values = NULL ) const;
	
	int checkTripleSparseStructure( unsigned int nzs, const unsigned int *rows, const unsigned int *cols, const double *values = NULL ) const;
	
	
    
    bool compareSparseMatrixStructure(const SPM_SparseMatrix<matrixType>& other) const;
	
	int copyStructureFrom(const SPM_SparseMatrix& M);
	
	int copyMatrixFrom(SPM_SparseMatrix<matrixType> &M);
	
	
	void copyMatrixTo(matrixType **M, const bool initializeWithZero = true, const double factor = 1.0) const;
	
	
	//if Matrix == NULL, the method allocate an array  to store the coefficients and return. Otherwise, the method puts the coefficient in Matrix and return it
	matrixType* copyMatrixTo(matrixType* Matrix = NULL, bool initializeWithZero = true, const double factor = 1.0) const;
	
	
	//that method count, for each column, how many rows exists where the respective column appears in sparse matrix
	void countRowsEachColumn(int *counts, const bool accumulate = false) const;
	void countRowsEachColumn(unsigned int *counts, const bool accumulate = false) const;
	
	int deleteRowStructure(const unsigned int row);
	
    void deleteStructure();
	
	void desallocateMemory(void);
	
    int getElement(unsigned int row, unsigned int col, matrixType &value) const;
	
	//that method fill all ncols in vector values
	int getFullRow(const unsigned int index, matrixType *values, const bool considerSymmetry = false ) const;
	
	
	//taht method accumulate in value the coefficients. Observe that vector values is not initialized
	int getFullRowAccumulation(const unsigned int row, matrixType *values) const;
	
	
	inline unsigned int getNumberOfColumns() const
	{
		return ncols;
	}
	
	
	inline unsigned int getNumberOfElements() const
	{
		return nElements;
	}
	
	//do not forget: symmetric matrix just store lower triangle...
	inline unsigned int getNumberOfElementsAtRow(const unsigned int row) const
	{
		return rows[row].getNumberOfElements();
	}
	
	inline unsigned int getNumberOfRows() const
	{
		return nrows;
	}
	
	//return the number of elements in Row...
	int getRowStructure(const unsigned int line, unsigned int* jCol) const;
	
	//return the number of elements in Row...
	int getRowStructure(const unsigned int line, bool* cols, bool accumulate) const;
	
	
	
	//that method put true in cols[i] if i is a column in line. If accumulate is true, method does not initialize cols with false.
	int getRowStructure(const unsigned int line, int* jCol) const;
	
	
	
	int getRowStructureAndValues(const unsigned int row, int &nzs, int *jCol, matrixType *vals) const;
	int getRowStructureAndValues(const unsigned int row, unsigned int &nzs, unsigned int *jCol, matrixType *vals) const;
	
	
	
	int getRowValues(const unsigned int line, matrixType *vals) const;
	
	
	int getStructure(int* iRow, int* jCol) const;
	int getStructure(unsigned int* iRow, unsigned int* jCol) const;
	
	
	int getValues(matrixType *values) const;
	
	
	//return the size of Arrays...
    int getStructureAndValues(int* iRow, int* jCol, matrixType* vals) const;
	int getStructureAndValues(unsigned int* iRow, unsigned int* jCol, matrixType* vals) const;
	
	
	
	//If any of the vectors row, cols or values is NULL, we allocate memory for the one. Otherwise, we ASSUME that a enough quantity of memory is already allocated. Allocations for rows, cols and values are performed using new operator. It is your responsability free the memory after the usage.
	int getTripleSparseFormat(int** rows, int** cols, matrixType** values) const;
	int getTripleSparseFormat(unsigned int** rows, unsigned int** cols, matrixType** values) const;
	
	
	bool hasIndex(const unsigned int row, const unsigned int col) const;
	
	int initialize(const unsigned int nrows = 0, const unsigned int ncols= 0, const bool symmetric = false);
	
	
	int mergeStructures(spm::SPM_SparseMatrix<matrixType>& M, const bool storeOrigCoefs = false);
	
	
	void multiplyAllElements( const matrixType value );
	
	void printSparseMatrix(const bool showEmptyLines = false) const;
    
    void printAllMatrix(void) const;
	
	//That function calculate 0.5* x'M x to a symmetric matrix M
	double quadraticEvaluation(const matrixType* x) const;
	
	
	/* Function to calculate the gradient of 0.5* x'Mx in grad vector. We suppose M is symmetric
	 * 
	 *Warning: If accumulate == true, That function only accumulates the gradient in grad, i.e., you should initialize grad vector by yourself. */
	void quadraticGradientEvaluation(const matrixType* x, matrixType* grad, const bool accumulate = false) const;
	
	
	/* Function to calculate 0.5* x'M x and the gradient of 0.5* x'Mx in grad vector. We suppose M is symmetric
	 * 
	 *Warning: If accumulate == true, That function only accumulates the gradient in grad, i.e., you should initialize grad vector by yourself. */
	double quadraticEvaluationAndGradient( const matrixType* x, matrixType* grad, const bool accumulate) const;
	
	
	int removeCols( const unsigned int ncols, const int *cols );
	
	int removeCols( const unsigned int ncols, const unsigned int *cols );
	
	//we remove the columns marked as true
	void removeCols( const bool *cols);
	
	void removeCol( const unsigned int col );
	
	int removeRows( const unsigned int nrows, const int *rows );
	int removeRows( const unsigned int nrows, const unsigned int* rows );
	
	//we remove the lines marked as true
	void removeRows( const bool *rows);
	
	//that method calculates M[row] * x 
	double rowEvaluation(const unsigned int row, const matrixType* x) const;
	
	void setAllElementsInARow(const unsigned int line, const matrixType value);
	
    void setAllSparseMatrix(const matrixType value);
    
    //This functions assumes that the position is in the matrix
    int setElement(unsigned int row, unsigned int col, const matrixType value, const bool printErrMsg = true);
	
	//this fucntion no check if col < row to symmetric matrices...
	int setElementNoCheckSymmetry(unsigned int row, unsigned int col, const matrixType value, const bool printErrMsg = true);
	
	
	inline void setNumberOfColumns(unsigned int ncols)
	{
		this->ncols = ncols;
	}
	
	
	//set in the row, col makerd as true
	int setRowStructure(const unsigned int row, const unsigned int sizeCols, const bool *cols);
	
	
	int setRowStructure(const unsigned int row, const unsigned int numberOfNZElements, const unsigned int* cols);
	int setRowStructure(const unsigned int row, const unsigned int numberOfNZElements, const int* cols);
	
	
	int setRowStructure(const unsigned int rowIndex, SPM_SparseRow< matrixType >& row);
	
	
	int setRowStructureAndValues(const unsigned int rowIndex, SPM_SparseRow< matrixType >& row);
	
	int setRowStructureAndValues(const unsigned int row, const unsigned int numberOfNZElements, const int *cols, const matrixType *vals);
	int setRowStructureAndValues(const unsigned int row, const unsigned int numberOfNZElements, const unsigned int *cols, const matrixType *vals);
	
	
	//if ncols = 0, spm will use the internal value of ncols...
	int setRowStructureAndValues(const unsigned int row, const matrixType *a, unsigned int ncols = 0, const double zeroTol = 0.0);
    
    int setStructureAndValues(const unsigned int nzs, const int* iRow, const int* jCol, const matrixType* vals = 0, const bool reset = false);
	
	int setStructureAndValues(const unsigned int nzs, const unsigned int* iRow, const unsigned int* jCol, const matrixType* vals = 0, const bool reset = false);
	
	
    
	//the matrix is stored in a vector ordered by row. if nrows or ncols is zero, spm will use the values stored inside object
	int setStructureAndValues(const matrixType* A, const double zero_tol, unsigned int nrows = 0, unsigned int ncols = 0);
	
    int setStructureAndValues(matrixType **A, const double zero_tol, unsigned int nrows = 0, unsigned int ncols = 0);
	
	
    int setSymmetricStructureAndValuesByLowerTriangle(const matrixType* lowerTA, const double zeroTolerance = 0.0);
	
	
	inline bool getSymmetricFlag() const
	{
		return symmetric;
	}
	
	inline void setSymmetricFlag(const bool flag)
	{
		symmetric = flag;
	}
	
	
	//that method can be used to premultiply the matrix by a vector: factor'M
	
	//that method sum all lines times a respective factor in vector. That vector is initialized with zeros. factors can be a null pointer. In this case, we consider all factors like 1.0
	void sumAllLines(const double *factors, matrixType *v) const;
	
	
	//overwrite a line copying from other sparse matrix. Symmetry is ignored... That function copy sourceline from M in destline 
	int copyLine(const unsigned int sourceLine, const unsigned int destLine, SPM_SparseMatrix<matrixType> &M);
	
	
	
	inline SPM_SparseRow<matrixType>& operator[] (unsigned int index)
	{
		return rows[index];
	}
	
	
    //equality operator
    //SPM_SparseMatrix & operator = ( SPM_SparseMatrix& other);
    
    ~SPM_SparseMatrix();
    
    friend int SPM_mergeSparseMatricesStructures <matrixType> (const SPM_SparseMatrix<matrixType> &M1, const SPM_SparseMatrix<matrixType> &M2, SPM_SparseMatrix<matrixType> &Res );
	
	
	
	//iterators for lines...
	class iterator: public SPM_Iterator< SPM_SparseMatrix<matrixType> >
	{
		
	public:
		
		iterator( SPM_SparseMatrix<matrixType> *M = NULL, unsigned int initialPosition = 0 ) : SPM_Iterator< SPM_SparseMatrix<matrixType> >( M, initialPosition )
		{
		}
		
		inline SPM_SparseRow<matrixType>& operator*()
		{
			const unsigned int pos = SPM_Iterator< SPM_SparseMatrix<matrixType> >::pos;
			
			return SPM_Iterator< SPM_SparseMatrix<matrixType> >::obj->rows[ pos ];
		}
		
	};
	
	
	
protected:
	
	//we define the size of the matrix as the number of rows to keep compatibility with 
	inline unsigned int iteratorSize() const
	{
		return nrows;
	}
	
	
	
	//we need this template because we wanna provide a interface for cols as *int or *unsigned int
	template <class intType>
	void __countRowsEachColumn(intType* counts, const bool accumulate) const;
	
	template <class intType>
	int __getStructure(intType* iRow, intType* jCol) const;
	
	template <class intType>
	int __getStructureAndValues(intType *iRow, intType *jCol, matrixType *vals) const;
	
	template <class intType>
	int __getTripleSparseFormat( intType** rows, intType** cols, matrixType** values) const;
	
	template <class intType>
	int __removeCols( const unsigned int ncols, const intType *cols );
	
	template <class intType>
	int __removeRows( const unsigned int nlines, const intType *lines );
	
	template <class intType>
	int __setRowStructure(  const unsigned int row, const unsigned int nzs, const intType* cols );
	
	template <class intType>
	int __setRowStructureAndValues(const unsigned int line, const unsigned int numberOfNZElements, const intType *cols, const matrixType *vals);
	
	template <class intType>
	int __setStructureAndValues(const unsigned int nzs, const intType* iRow, const intType* jCol, const matrixType* vals, const bool reset);
	
	
	friend class SPM_Iterator< SPM_SparseMatrix<matrixType> >;
	
};











template <class intType>
inline int SPM_checkRowStructure(const unsigned int nrows, const unsigned int ncols, const bool symmetric, unsigned int row, unsigned int nzs, intType *cols, double *values = NULL)
{
	if( row < 0 || row >= nrows )
		return SPM_INDEX_FAULT;
	
	
	for(unsigned int i = 0; i < nzs; i++)
	{
		if( cols[i] < 0 || cols[i] >= nrows )
			return SPM_INDEX_FAULT;
		
		for(unsigned int j = i+1; j < nzs; j++)
		{
			if( cols[i] == cols[j] )
			{
				return SPM_REPETEAD_ELEMENT;
			}
		}
	}
	
	if( values )
	{
		for(unsigned int i = 0; i < nzs; i++)
		{
			if( isinf(values[i]) || isnan(values[i]) )
				return SPM_BAD_VALUE;
		}
	}
	
	if(symmetric)
	{
		for(unsigned int i = 0; i < nzs; i++)
		{
			if( cols[i] > row )
				return SPM_UPPER_TRIANGLE;
		}
	}
	
	return 0;
}




template <class intType>
inline int SPM_checkTripleSparseStructure(const unsigned int nrows, const unsigned int ncols, const bool symmetric,  unsigned int nzs, const intType *rows, const intType *cols, const double *values = NULL )
{
	//bool repeat = false, badvalue = false, //upperTriangle = false;
	
	
	for(unsigned int i = 0; i < nzs; i++)
	{
		for(unsigned int j = i+1; j < nzs; j++ )
		{
			if( rows[i] == rows[j] && cols[i] == cols[j] )
			{
				//repeat = true;
				return SPM_REPETEAD_ELEMENT;
			}
		}
	}
	
	
	for(unsigned int i = 0; i < nzs; i++)
	{
		if( rows[i] < 0 || rows[i] >= nrows )
			return SPM_INDEX_FAULT;
	}
	
	for(unsigned int i = 0; i < nzs; i++)
	{
		if( cols[i] < 0 || cols[i] >= ncols )
			return SPM_INDEX_FAULT;
	}
	
	
	
	if( values )
	{
		for(unsigned int i = 0; i < nzs; i++)
		{
			if( isinf(values[i]) || isnan(values[i]) )
			{
				//badvalue = true;
				return SPM_BAD_VALUE;
			}
		}
	}
	
	
	if( symmetric )
	{
		for( unsigned int i = 0; i < nzs; i++ )
		{
			if( rows[i] < cols[i] )
			{
				//upperTriangle = true;
				return SPM_UPPER_TRIANGLE;
			}
		}
	}
	
	
	return 0;
	
	
	
}





/*inline double SPM_min2( const double a, const double b );
inline int SPM_min2( const int a, const int b );
inline double SPM_max2( const double a, const double b );
inline int SPM_max2( const int a, const int b );
inline double SPM_abs2( const double num );
inline int SPM_abs2( const int num ); */





}


#include "SPM_SparseMatrixImp.hpp"


#endif
