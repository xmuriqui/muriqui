

/***********************************************************
***********************************************************
****                                                   ****
****          new sparse matrix implementation          ****
****                                                   ****
***********************************************************
***********************************************************/

#ifndef SPM_NEW_SPARSE_MATRIX_HPP
#define SPM_NEW_SPARSE_MATRIX_HPP

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include "../WAXM_config.h"


#include "SPM_tools.hpp"



namespace newspm
{





template <class indexType, class matrixType>
class SPM_ColsValues
{
public:
    
    indexType* cols;
    matrixType* values;
    
    SPM_ColsValues()
    {
        initialize();
    }
    
    // We prefer do not build a destructor. User has to call desallocateMemory if he wanna free the memory (we want change cols and values pointer to eprform some computations)
    
    void initialize()
    {
        cols = NULL;
        values = NULL;
    }
    
    int reallocateArrays(const unsigned int newSize)
    {
        if( newSize == 0 )
        {
            desallocateMemory();
        }
        else
        {
            indexType* pind;
            matrixType* pval;
            
            pind = (indexType*) realloc( cols, newSize * sizeof(indexType) );
            
            if( !pind )
            {
                #if SPM_DEBUG_MODE
                    SPM_PRINTMEMERROR;
                #endif
                return SPM_MEMORY_ERROR;
            }
            
            cols = pind;
            
            pval = (matrixType*) realloc( values, newSize * sizeof(matrixType) );
            
            if( !pval )
            {
                #if SPM_DEBUG_MODE
                    SPM_PRINTMEMERROR;
                #endif
                return SPM_MEMORY_ERROR;
            }
            
            values = pval;
        }
        
        return 0;
    }
    
    void desallocateMemory()
    {
        SPM_secFree(cols);
        SPM_secFree(values);
    }
    
};



template <class indexType, class matrixType>
class SPM_NewSparseRow
{
    /*
    * for improve eficiency with solvers, we represent in this way:
    * 
    * Index columuns array for a row is given by (baseStructure.cols)[__offset]. 
    * Values array for a row is given by (baseStructure.values)[__offset].
    */
    
    
    
    
public:
    
    indexType nElements;
    unsigned int __offset;
    SPM_ColsValues<indexType, matrixType> *baseStructure;
    
    SPM_NewSparseRow(SPM_ColsValues<indexType, matrixType> *baseStructure = NULL, const unsigned int offset = 0);
    
    ~SPM_NewSparseRow();
    
    inline indexType* getColsPointer(){ return  &(baseStructure->cols[__offset]);
    }
    
    inline matrixType* getValuesPointer(){
        return &(baseStructure->values[__offset]);
    }
    
    
    //this method sum the coefficient in this line times a factor in a array. The vector is not iniatialized (note there is no symmetry concecpt for a line operation).
    void accumulateInArray(matrixType *v, const matrixType factor = 1) const;
    
    int addToElement(const indexType col, const matrixType value);
    
    void addToAllElements(const matrixType value);
    
    void copyTo(matrixType* row, const matrixType factor = 1);
    
    void desallocate();
    
    double evalTimesxt(const matrixType* x) const;
    
    int getElement(const indexType col, matrixType &value) const;
    
    inline indexType getNumberOfElements() const
    {
        return nElements;
    }
    
    //that method return the number of elements in the row. You must initialize by yourself all positions in cols as false 
    indexType getStructure(bool* columns) const;
    
    //those methods return the number of elements in the row
    indexType getStructure(int* columns) const;
    indexType getStructure(unsigned int *cols) const;
    
    unsigned int getValues(matrixType *values) const;
    
    bool hasColumn(const indexType col, indexType* index = NULL, matrixType* value = NULL) const;
    
    void initialize(SPM_ColsValues<indexType, matrixType> *baseStructure = NULL, const unsigned int offset = 0);
    
    void multiplyAllElements(const matrixType value);
    
    void print( std::ostream &out = std::cout );
    
    void setAllElements(const matrixType value);
    
    int setElement(const indexType col, const matrixType value);
    
    //set first 'nel' elemets in this line using the first 'nel' values in array values. Note, this method does not look to columns indexes...
    int setElementsByOrder(indexType nel, const matrixType* values);
    
    
    
    //return the element of row in position pos
    matrixType& operator()(const unsigned int pos);
    
    //return the column of row in position pos
    indexType& operator[](const unsigned int pos);
    
    
    
private:
    //we do not let atributions. So, we declare the equality operator as private
    
    SPM_NewSparseRow& operator=(SPM_NewSparseRow &other)
    {
        assert(false);
    }
    
    void setOffset( unsigned int deslocamento )
    {
        deslocamento = deslocamento;
    }
    
    //friend class SPM_NewSparseMatrix<indexType, matrixType>;
};







template <class indexType, class matrixType>
class SPM_NewSparseMatrix
{
    
public:
    bool symmetric;
    unsigned int nrows, ncols;
    unsigned int nElements;
    
    //SPM_NewSparseRow<indexType, matrixType> *rows;
    
    unsigned int *offset; //offsets to rows
    
    
    unsigned int *nzRowIndex; //indices of rows having nonzero values.
    unsigned int nNzRowIndices; //number of rows having nonzero values
    
    SPM_ColsValues<indexType, matrixType> baseStructure;
    
    
    
    SPM_NewSparseMatrix(const unsigned int nrows = 0, const unsigned int ncols = 0, const bool symmetric = false);
    
    ~SPM_NewSparseMatrix();
    
    inline indexType* getRowColsPointer(unsigned int rowIndex) const
    {
        if( baseStructure.cols )
            return &(baseStructure.cols[ offset[rowIndex] ]);
        else
            return NULL;
    }
    
    inline matrixType* getRowValuesPointer(unsigned int rowIndex) const
    {
        if( baseStructure.values )
            return &(baseStructure.values[ offset[rowIndex] ]);
        else
            return NULL;
    }
    
    inline void getRowPointers( unsigned int rowIndex, unsigned int &nz, indexType* &cols, matrixType* &values ) const
    {
        nz = getNumberOfElementsAtRow(rowIndex);
        cols = getRowColsPointer(rowIndex);
        values = getRowValuesPointer(rowIndex);
    }
    
    //that method accumulates the sum between this matrix and factor*A. All positions in A should be also in this matrix. 
    //void accumulatedSum(const double factor, SPM_NewSparseMatrix<indexType, matrixType>& A);
    
    //this method accumulates in M the current sparse matrix times factor ( M = M + facotr*this). Note, the result is stored in M and the sparse matrix is NOT changed...
    void accumulatedSumToMatrix( matrixType *M, const matrixType factor = 1, const bool considerSymmetry = true) const;
    
    //that method sum the coefficients in a line times a factor in a vector. The vector is not initialized
    int accumulateRowInArray(const unsigned int rowIndex, matrixType* v, const matrixType factor) const;
    
    
    int addNewRows(const unsigned int nrows);
    
    
    void addToAllElements(const matrixType value);
    
    int addToAllElementsAtRow(const unsigned int rowIndex, const matrixType value);
    
    //This functions assumes that the position is in the matrix. The current value in the position will be added to value parameter.
    int addToElement(unsigned int row, indexType col, const matrixType value);
    
    
    bool dep_compareSparseMatrixStructure(const SPM_NewSparseMatrix<indexType, matrixType>& other) const;
    
    
    int copyMatrixFrom(const SPM_NewSparseMatrix<indexType, matrixType> &M);
    
    void copyMatrixTo(matrixType **M, const bool considerSymmetry, const bool initializeWithZero, const double factor) const;
    
    //if Matrix == NULL, the method allocate an array  to store the coefficients and return. Otherwise, the method puts the coefficient in Matrix and return it
    matrixType* copyMatrixTo(matrixType* Matrix, const bool considerSymmetry, const bool initializeWithZero, const double factor) const;
    
    void copyParametersFrom(const SPM_NewSparseMatrix<indexType, matrixType>& M);
    
    //this method does not consider symmetry and does not itialize the row array. If you need this resources, use getFullRow method.
    void copyRowTo(const unsigned int rowIndex, matrixType* row, const matrixType factor = 1) const;
    
    int copyStructureFrom(const SPM_NewSparseMatrix<indexType, matrixType>& M);
    
    //that method count, for each column, how many rows exists where the respective column appears in sparse matrix
    void countRowsEachColumn(int *counts, const bool accumulate = false) const;
    void countRowsEachColumn(unsigned int *counts, const bool accumulate = false) const;
    
    int deleteRowStructure(const unsigned int rowIndex);
    
    void deleteStructure();
    
    void desallocateMemory(void);
    
    //eval Matrix[rowIndex] *x'. Warning: this method will ignore simetry
    double evalRowTimesxt(const indexType rowIndex, const double *x) const;
    
    //eval factor *Matrix *x'
    void evalTimesxt(const double *x, double *v, const bool accumulate = false, double factor = 1.0) const;
    
    
    //eval factor x*Matrix
    void evalPretimesx(const double *x, double *v, const bool accumulate = false, double factor = 1.0) const;
    
    //eval factor x*Matrix only in some rows.
    void evalPretimesx(const unsigned int nRowInds, const unsigned *rowInds, const double *x, double *v, const bool accumulate = false, double factor = 1.0) const;
    
    
    int getElement(unsigned int rowIndex, indexType colIndex, matrixType& value) const;
    
    //that method fill all ncols in vector values
    int getFullRow(const unsigned int rowIndex, matrixType* values, const bool accumulate, const bool considerSymmetry, const double factor) const;
    
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
        //return rows[index].getNumberOfElements();
        return offset[row+1] - offset[row];
    }
    
    //return INDEX error if index is invalid
    inline int getNumberOfElementsAtRow(const unsigned int row, unsigned int &nzs) const
    {
        if(row >= nrows)
            return SPM_INDEX_FAULT;
            
        nzs = getNumberOfElementsAtRow(row);
        return 0;
    }
    
    inline unsigned int getNumberOfNonZeroRows() const
    {
        return nNzRowIndices;
    }
    
    inline unsigned int getNumberOfRows() const
    {
        return nrows;
    }
    
    
    int getRowStructure(const unsigned int rowIndex, indexType* columns, indexType* nzs) const;
    
    //that method put true in cols[i] if i is a column in line. If accumulate is true, method does not initialize cols with false.
    int getRowStructure(const unsigned int rowIndex, bool* columns, indexType* nzs, bool accumulate) const;
    
    //int getRowStructureAndValues(const unsigned int index, indexType *nzs, indexType *jCol, matrixType *vals) const;
    
    int getRowValues(const unsigned int rowIndex, matrixType* vals, indexType* nzs) const;
    
    int getStructure(int* iRow, indexType* cols) const;
    int getStructure(unsigned int* iRow, indexType* cols) const;
    
    void getStructureByCompressedRowFormat(int* rowStart, indexType* cols, matrixType* values = NULL) const;
    
    void getStructureByCompressedRowFormat( unsigned int *rowStart, indexType *cols, matrixType* values = NULL) const;
    
    int getValues(matrixType *values) const;
    
    /*//return the size of Arrays...
    int getStructureAndValues(int* iRow, indexType* jCol, matrixType* vals) const;
    int getStructureAndValues(unsigned int* iRow, indexType* jCol, matrixType* vals) const;*/
    
    bool hasIndex(const unsigned int rowIndex, const indexType col, matrixType* value = NULL, indexType* indexAtRowStructure = NULL) const;
    
    int initialize(const unsigned int nrows = 0, const indexType ncols= 0, const bool symmetric = false, const bool initialize = true );
    
    int mergeStructures(SPM_NewSparseMatrix<indexType, matrixType>& M, const bool storeOrigCoefs = false);
    
    void multiplyAllElements( const matrixType value );
    
    int multiplyAllElementsAtRow( const unsigned int rowIndex, const matrixType value );
    
    void printSparseMatrix(const bool showEmptyLines = false) const;
    
    void printAllMatrix(void) const;
    
    //That function calculate 0.5* x'M x to a symmetric matrix M
    double quadraticEvaluation(const double* x) const;
    
    /* Function to calculate the gradient of 0.5* x'Mx in grad vector. We suppose M is symmetric
    * 
    *Warning: If accumulate == true, That function only accumulates the gradient in grad, i.e., you should initialize grad vector by yourself. */
    void quadraticGradientEvaluation(const double* x, double* grad, const bool accumulate) const;
    
    /* Function to calculate 0.5* x'M x and the gradient of 0.5* x'Mx in grad vector. We suppose M is symmetric
    * 
    *Warning: If accumulate == true, That function only accumulates the gradient in grad, i.e., you should initialize grad vector by yourself. */
    double quadraticEvaluationAndGradient( const double* x, double* grad, const bool accumulate) const;
    
    //we remove the columns marked as true
    int removeCols(const bool* cols);
    
    int removeCols(const unsigned int ncols, const indexType* cols);
    
    //we remove the columns marked as true. Even if matrix is simmetric, we do not remove lines. We do not change ncols value. In the practice, is like you remove all coeficients in the marked columns
    void removeColsCoefficients( const bool *cols);
    
    int removeColsCoefficients( const unsigned int ncols, const indexType *cols );
    
    int removeRows( const unsigned int nrows, const int *rows );
    int removeRows( const unsigned int nrows, const unsigned int* rows );
    
    //we remove the lines marked as true
    int removeRows( const bool* rowsFlag);
    
    //that method calculates M[row] * x 
    double rowEvaluation(const unsigned int row, const matrixType* x) const;
    
    void setAllElementsInARow(const unsigned int row, const matrixType value);
    
    void setAllSparseMatrix(const matrixType value);
    
    //This functions assumes that the position is in the matrix
    int setElement(unsigned int rowIndex, indexType colIndex, const matrixType value, const bool printErrMsg = false);
    
    //this fucntion no check if col < row to symmetric matrices...
    int setElementNoCheckSymmetry(unsigned int row, indexType col, const matrixType value, const bool printErrMsg = true);
    
    inline void setNumberOfColumns(unsigned int ncols)
    {
        this->ncols = ncols;
    }
    
    //set in the row, col makerd as true
    int setRowStructure(const unsigned int rowIndex, const unsigned int sizeCols, const bool *cols);
    
    int setRowStructure(const unsigned int rowIndex, const indexType numberOfNZElements, const indexType* cols, const matrixType *values = NULL);
    
    int setRowStructure(const unsigned int rowIndex, const unsigned int sourceRowIndex, const SPM_NewSparseMatrix< indexType, matrixType >& M, const bool copyValues = true);
    
    //int setRowStructureAndValues(const unsigned int rowIndex, const unsigned int sourceRowIndex,  const SPM_NewSparseMatrix<indexType, matrixType>& M);
    
    //int setRowStructureAndValues(const unsigned int rowIndex, const indexType numberOfNZElements, const indexType* cols, const matrixType* values);
    
    //if ncols = 0, spm will use the internal value of ncols...
    int setRowStructureAndValues(const unsigned int rowIndex, const matrixType *a, unsigned int ncols = 0, const double zeroTol = 0.0);
    
    //set row values by order in the row. Values should have getNumberOfElementsAtRow(rowIndex) elements at least
    void setRowValues(const unsigned int rowIndex, const matrixType *values);
    
    //set row values taking a full row. Just columns in the row will be set 
    void setFullRowValues(const unsigned int rowIndex, const matrixType *values);
    
    
    int setStructureAndValues(const unsigned int nzs, const int* iRow, const indexType* jCol, const matrixType* vals = 0, const bool reset = false);
    
    int setStructureAndValues(const unsigned int nzs, const unsigned int* rows, const indexType* cols, const matrixType* vals = 0, const bool reset = false);
    
    //the matrix is stored in a vector ordered by row. if nrows or ncols is zero, spm will use the values stored inside object
    int setStructureAndValues(const matrixType* A, const double zeroTolerance, unsigned int nrows, unsigned int ncols);
    
    int setStructureAndValues(matrixType** A, const double zeroTolerance, unsigned int nrows, unsigned int ncols);
    
    //set by compressed row format
    int setStructureAndValues(const unsigned int *rowStart, const indexType *cols, const matrixType* vals = NULL);
    
    //set by compressed row format
    int setStructureAndValues(const int *rowStart, const indexType *cols, const matrixType* vals = NULL);
    
    int setSymmetricStructureAndValuesByLowerTriangle(const matrixType* lowerTA, const double zeroTolerance = 0.0);
    
    inline bool getSymmetricFlag() const
    {
        return symmetric;
    }
    
    inline void setSymmetricFlag(const bool flag)
    {
        symmetric = flag;
    }
    
    //that method sum all lines times a respective factor in vector. That vector is initialized with zeros. factors can be a null pointer. In this case, we consider all factors like 1.0
    void sumAllLines(const double *factors, matrixType *v) const;
    
    
    inline indexType* operator[] (unsigned int index) const
    {
        return getRowColsPointer(index);
    }
    
    inline matrixType* operator() (unsigned int index) const
    {
        return getRowValuesPointer(index);
    }
    
    //Note: this class iterates on nonzero elements in the  matrix
    class Iterator
    {
        class SPM_Element
        {
            unsigned int row;
            
        public:
            
            indexType *col;
            matrixType *value;
            
            SPM_Element()
            {
                row = -1; //maximum value
                col = NULL;
                value = NULL;
            }
            
            unsigned int getRow() const {return row;}
            
            indexType getColumn() const {return *col;}
            
            matrixType getValue() const {return *value;}
            
            void setColumn(indexType column) const { *col = column;}
            
            void setValue(matrixType value) const { *(this->value) = value;}
            
            
            friend class newspm::SPM_NewSparseMatrix<indexType, matrixType>::Iterator;
        };
        
        const SPM_NewSparseMatrix<indexType, matrixType> *matrix;
        
        unsigned int rownz;
        unsigned int pos;
        SPM_Element element;
        
    public:
        
        Iterator( const SPM_NewSparseMatrix<indexType, matrixType> *matrix, const unsigned int startRow = 0 )
        {
            initialize(matrix, startRow);
        }
        
        void initialize( const SPM_NewSparseMatrix<indexType, matrixType> *matrix, const unsigned int startRow = 0 )
        {
            const unsigned int nrows = matrix->nrows;
            const unsigned int nNzRowIndices = matrix->nNzRowIndices;
            const unsigned int *nzRowIndex = matrix->nzRowIndex;
            
            this->matrix = matrix;
            pos = 0;
            rownz = nNzRowIndices; //we initialize like end...
            
            if( startRow < nrows )
            {
                for( unsigned int i = 0; i < nNzRowIndices; i++ )
                {
                    if( nzRowIndex[i] >= startRow )
                    {
                        rownz = i;
                        break;
                    }
                }
            }
        }
        
        unsigned int getCurrentRow()
        {
            return matrix->nzRowIndex[rownz];
        }
        
        indexType getCurrentCol()
        {
            return (*matrix)[getCurrentRow()][pos];
        }
        
        matrixType getCurrentValue()
        {
            return (*matrix)(getCurrentRow())[pos];
        }
        
        
        void setCurrentCol(const indexType col)
        {
            (*matrix)[getCurrentRow()][pos] = col;
        }
        
        void setCurrentValue(const matrixType value)
        {
            (*matrix)(getCurrentRow())[pos] = value;
        }
        
        
        bool operator==( const Iterator &iter ) const
        {
            return matrix == iter.matrix && rownz == iter.rownz && pos == iter.pos;
        }
        
        bool operator!=( const Iterator &iter )const
        {
            return !(*this == iter);
        }
        
        SPM_Element& operator*()
        {
            element.row = getCurrentRow();//matrix->nzRowIndex[rownz];
            element.col = &matrix->getRowColsPointer(element.row)[pos];
            element.value = &matrix->getRowValuesPointer(element.row)[pos];
            
            return element;
        }
        
        
        Iterator& operator++()
        {
            if( rownz < matrix->nNzRowIndices )
            {
                const unsigned int nzcrow = matrix->getNumberOfElementsAtRow(matrix->nzRowIndex[rownz]);
                pos++;
                
                #if SPM_DEBUG_MODE
                    assert( pos <= nzcrow );
                #endif
                
                if( pos == nzcrow )
                {
                    pos = 0;
                    rownz++;
                }
            }
            
            return *this;
        }
        
        
        Iterator& operator--()
        {
            if(rownz > 0 || pos > 0)
            {
                if( pos == 0 )
                {
                    rownz--;
                    pos = matrix->getNumberOfElementsAtRow( matrix->nzRowIndex[rownz] );
                    
                    #if SPM_DEBUG_MODE
                        assert(pos > 0); //so we have a row having no elements in nzRowIndex
                    #endif
                }
                
                pos--;
            }
            
            return *this;
        }
        
    };
    
    //Note: This class Iterates only on row indices having some element.
    class RowIndexIterator
    {
        const SPM_NewSparseMatrix<indexType, matrixType> *matrix;
        
        unsigned int rownz;
        
    public:
        
        RowIndexIterator( const SPM_NewSparseMatrix<indexType, matrixType> *matrix, const unsigned int startRow = 0 )
        {
            initialize(matrix, startRow);
        }
        
        void initialize( const SPM_NewSparseMatrix<indexType, matrixType> *matrix, const unsigned int startRow = 0 )
        {
            const unsigned int nrows = matrix->nrows;
            const unsigned int nNzRowIndices = matrix->nNzRowIndices;
            const unsigned int *nzRowIndex = matrix->nzRowIndex;
            
            this->matrix = matrix;
            rownz = nNzRowIndices; //we initialize like end...
            
            if( startRow < nrows )
            {
                for( unsigned int i = 0; i < nNzRowIndices; i++ )
                {
                    if( nzRowIndex[i] >= startRow )
                    {
                        rownz = i;
                        break;
                    }
                }
            }
        }
        
        bool operator==(const RowIndexIterator &iter) const
        {
            return matrix == iter.matrix && rownz == iter.rownz;
        }
        
        bool operator!=(const RowIndexIterator &iter) const
        {
            return !(*this == iter);
        }
        
        unsigned int operator*() const
        {
            if( matrix->nNzRowIndices > 0 )
                return matrix->nzRowIndex[rownz];
            else //nzRowIndex has no indices
                return matrix->getNumberOfRows();
        }
        
        RowIndexIterator& operator++()
        {
            if( rownz < matrix->nNzRowIndices )
                rownz++;
            
            return *this;
        }
        
        RowIndexIterator& operator--()
        {
            if( rownz > 0 )
                rownz--;
            
            return *this;
        }
    };
    
    
    Iterator begin( const unsigned int indStartRow = 0 ) const { return Iterator(this, indStartRow);}
    
    //if you do not pass the end row, we consider the last row
    Iterator end( const unsigned int indEndRow = -1 ) const 
    { 
        if( indEndRow >= nrows )
            return Iterator(this, nrows);
        else
            return Iterator(this, indEndRow+1);
    }
    
    RowIndexIterator beginRowIndex(const unsigned int indStartRow = 0) const
    {
        return RowIndexIterator(this, indStartRow);
    }
    
    RowIndexIterator endRowIndex(const unsigned int indEndRow = -1 ) const 
    {
        if( indEndRow >= nrows )
            return RowIndexIterator(this, nrows);
        else
            return RowIndexIterator(this, indEndRow+1);
    }
    
    
    int checkConsistency(std::ostream& out = std::cerr) const;
    
private:
    
    int allocateSparseRows(const unsigned int nrows);
    
    int checkNElements() const;
    
    void updateNzRowIndex();
    
    int reallocateColsAndValues(const unsigned int newSize);
    
    
    int reallocateRowStructures(const unsigned int newNumberOfRows);
    
    void leftShiftColsValuesAndOffset(const unsigned int startRow, const unsigned int shift);
    
    void rightShiftColsValuesAndOffset(const unsigned int startRow, const unsigned int shift);
    
    //note leftShiftColsValuesAndOffset already call this function
    inline void shiftPlusRowOffset( const unsigned int initialRowIndex, const unsigned int shiftValue )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = initialRowIndex; i <= nrows; i++)
            offset[i] += shiftValue; //rows[i].__offset += shiftValue;
    }
    
    //note rightShiftColsValuesAndOffset already call this function
    inline void shiftMinusRowOffset( const unsigned int initialRowIndex, const unsigned int shiftValue )
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = initialRowIndex; i <= nrows; i++)
            offset[i] -= shiftValue; //rows[i].__offset -= shiftValue;
    }
    
    
    int reallocateRowSpace(const unsigned int rowIndex, const unsigned int newsize);
    
    
    inline const matrixType* getRowMatrixPointer(const matrixType* M, const unsigned int ncols, const unsigned int rowIndex, const bool lowerTriangle = false) const
    {
        return lowerTriangle ? &M[(rowIndex*(rowIndex+1))/2] : &M[rowIndex*ncols];
    }

    inline const matrixType* getRowMatrixPointer(const matrixType** M, const unsigned int ncols, const unsigned int rowIndex, const bool lowerTriangle = false) const
    {
        return M[rowIndex];
    }
    
    
    //we need this template because we wanna provide a interface for cols as *int or *unsigned int
    template <class intType>
    void __countRowsEachColumn(intType* counts, const bool accumulate) const;
    
    template <class intType>
    indexType __getStructure(intType* rows, indexType* cols) const;
    
    template <class intType>
    void __getStructureByCompressedRowFormat( intType* rowStart, indexType* cols, matrixType* vals) const;
    
    template <class intType>
    int __removeRows( const unsigned int nrows, const intType* rows );
    
    template <class intType>
    int __setStructureAndValues(const unsigned int nzs, const intType* rows, const indexType* cols, const matrixType* vals, const bool reset);
    
    template <class matrixTypePointer>
    int __setStructureAndValues(matrixTypePointer A, const double zero_tol, unsigned int nrows, unsigned int ncols, const bool lowerTriangle);
    
    template <class intType>
    int __setStructureAndValues(const intType *rowStart, const indexType *cols, const matrixType* vals);
};











template <class indexType, class matrixType>
SPM_NewSparseMatrix<indexType, matrixType>::SPM_NewSparseMatrix(const unsigned int nrows, const unsigned int ncols, const bool symmetric)
{
    initialize(nrows, ncols, symmetric);
}

template <class indexType, class matrixType>
SPM_NewSparseMatrix<indexType, matrixType>::~SPM_NewSparseMatrix()
{
    desallocateMemory();
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::updateNzRowIndex()
{
    nNzRowIndices = 0;
    
    for(unsigned int i = 0; i < nrows; i++)
    {
        if( getNumberOfElementsAtRow(i) > 0 ) //( rows[i].nElements > 0 )
        {
            nzRowIndex[ nNzRowIndices ] = i;
            nNzRowIndices++;
        }
    }
    
    nzRowIndex[nNzRowIndices] = nrows; //we sinalize the end of nzRowIndex with nrows value. It can be useful to iterators
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::reallocateColsAndValues(const unsigned int newSize)
{
    return baseStructure.reallocateArrays(newSize);
}


//this method accumulates in M the current sparse matrix times factor ( M = M + facotr*this). Note, the result is stored in M and the sparse matrix is NOT changed...
template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::accumulatedSumToMatrix(matrixType *M, const matrixType factor, const bool considerSymmetry) const
{
    if(factor == 0.0)
        return;
    
    for(decltype(nNzRowIndices) k = 0; k < nNzRowIndices; k++)
    {
        const auto rind = nzRowIndex[k];
        
        matrixType *rrow = &M[rind*ncols];
        
        accumulateRowInArray(rind, rrow, factor); //rows[rind ].accumulateInArray(rrow, factor);
    }
    
    
    if( symmetric && considerSymmetry )
    {
        for(decltype(nNzRowIndices) k = 0; k < nNzRowIndices; k++)
        {
            const auto rind = nzRowIndex[k];
            
            indexType* const cols = getRowColsPointer(rind);  //rows[rowIndex].getColsPointer();
        
            matrixType* const values = getRowValuesPointer(rind); // rows[rowIndex].getValuesPointer();
            
            const auto nel = getNumberOfElementsAtRow(rind); //rows[rowIndex].getNumberOfElements();
            
            
            for(unsigned int j = 0; j < nel; j++)
            {
                const decltype(rind) row = cols[j];
                
                if( row != rind )
                    M[ row*ncols + rind ] += factor * values[j];
            }
        }
    }
    
}

//that method sum the coefficients in a line times a factor in a vector. The vector is not initialized and symmetry is ignored
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::accumulateRowInArray(const unsigned int rowIndex, matrixType *v, const matrixType factor) const
{
    const unsigned int rnels = getNumberOfElementsAtRow(rowIndex);
    
    matrixType* const values = getRowValuesPointer(rowIndex);
    indexType* const cols = getRowColsPointer(rowIndex);
    
    
    if( factor == 1.0 ) 
    {
        for( unsigned int i = 0; i < rnels; i++)
            v[ cols[i] ] += values[i];
    }
    else
    {
        for( unsigned int i = 0; i < rnels; i++)
            v[ cols[i] ] += factor*values[i];
    }
    
    return 0;
}



template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::addNewRows(const unsigned int nrows)
{
    const unsigned int totalRows = this->nrows + nrows;
    /*unsigned int *uintp;
    newspm::SPM_NewSparseRow<indexType, matrixType> *p;
    
    p = (newspm::SPM_NewSparseRow<indexType, matrixType>) realloc( rows, totalRows * sizeof(SPM_NewSparseRow<indexType, matrixType>) );
    uintp = (unsigned int *) realloc( nzRowIndex, totalRows * sizeof(unsigned int) );
    
    if(!p || !uintp)
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTMEMERROR;
        #endif
        return SPM_MEMORY_ERROR;
    }
    
    rows = p;
    nzRowIndex = uintp;
    
    for(unsigned int i = this->nrows; i < totalRows; i++)
        rows[i].initialize(&baseStructure, nElements); */
    
    
    auto r = allocateSparseRows(totalRows);
    if(r != 0)
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    
    this->nrows = totalRows;
    
    
    return 0;
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::addToAllElements(const matrixType value)
{
    SPM_addToAllArray(nElements, baseStructure.values, value);
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::addToAllElementsAtRow(const unsigned int rowIndex, const matrixType value)
{
    if( rowIndex >= nrows )
        return SPM_INDEX_FAULT;
    
    SPM_addToAllArray(getNumberOfElementsAtRow(rowIndex), getRowValuesPointer(rowIndex), value);
    
    return 0;
}


//This functions assumes that the position is in the matrix. The current value in the position will be added to value parameter.
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::addToElement(unsigned int row, indexType col, const matrixType value)
{
    if(symmetric && (unsigned int) col > row)
        SPM_swap(col, row); //we work only in the lower triangle
    
    const auto rnels = getNumberOfElementsAtRow(row);
    indexType* const cols = getRowColsPointer(row);
    
    for(unsigned int i = 0; i < rnels; i++)
    {
        if( cols[i] == col )
        {
            getRowValuesPointer(row)[i] += value;
            return 0;
        }
    }
    
    return SPM_ELEMENT_NOT_PRESENT;
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::allocateSparseRows(const unsigned int nrows)
{
    /*//do not use new because we can need realloc
    
    rows = (SPM_NewSparseRow<indexType, matrixType> *) malloc( nrows * sizeof(SPM_NewSparseRow<indexType, matrixType>) );
    nzRowIndex = (unsigned int *) malloc( nrows * sizeof(unsigned int) );
    if( !rows || !nzRowIndex )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTMEMERROR;
        #endif
        return SPM_MEMORY_ERROR;
    }
    
    this->nrows = nrows;
    for(unsigned int i = 0; i < nrows; i++)
        rows[i].initialize(&baseStructure, nElements); */
    
    
    auto r = reallocateRowStructures(nrows );
    if(r != 0)
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    //initializing the new offsets
    
    /*#pragma ivdep
    #pragma GCC ivdep
    for(unsigned int i = this->nrows+1; i <= nrows; i++ )
        offset[i] = nElements;*/
    
    if( nrows > this->nrows )
    {
        SPM_setAllArray(nrows-this->nrows +1, &offset[this->nrows], nElements ); //offset now has size nrows+1. So, the last position is nrows. We start from nrows because if nrows is zero, the position nrows will be not initialized if we start from nrows+1
    }
    
    
    nzRowIndex[nNzRowIndices] = nrows;
    
    #if SPM_DEBUG_MODE
        if(nrows > nNzRowIndices)
            SPM_setAllArray(nrows-nNzRowIndices, &nzRowIndex[this->nNzRowIndices+1], (unsigned int)-1); //we set all index with the maximum index. So, if we use incorrectly, probably we will got an error. Note, nzRowIndex has nrows+1 positions and position nrows can be set as nrows value.
    #endif
    
    
    return 0;
}


template <class indexType, class matrixType>
bool SPM_NewSparseMatrix<indexType, matrixType>::dep_compareSparseMatrixStructure(const SPM_NewSparseMatrix<indexType, matrixType>& other) const
{
    if( nrows != other.nrows || ncols != other.ncols || nElements != other.nElements || symmetric != other.symmetric || nNzRowIndices != other.nNzRowIndices )
        return false;
    
    if((offset == NULL && other.offset != NULL) || (offset != NULL && other.offset == NULL))
        return false;
    
    //first, we check if getNumberOfElementsAtRow is the same for all rows. 
    for(unsigned int i = 0; i < nNzRowIndices; i++)
    {
        //we assume nzRowIndex is always in order...
        if( nzRowIndex[i] != other.nzRowIndex[i] )
            return false;
        
        const unsigned int rind = nzRowIndex[i];
        
        if( getNumberOfElementsAtRow(rind) != other.getNumberOfElementsAtRow(rind) ) //( rows[rowIndex].nElements != other.rows[rowIndex].nElements )
            return false;
    }
    
    
    for(unsigned int i = 0; i < nNzRowIndices; i++)
    {
        const auto rind = nzRowIndex[i];
        const unsigned int rnel = getNumberOfElementsAtRow(rind); // rows[rowIndex].nElements;
        
        indexType* const rcols = getRowColsPointer(rind); 
        indexType* const otherrcols = other.getRowColsPointer(rind);
        
        
        for(unsigned int j = 0; j < rnel; j++)
        {
            //if the elements are in the same order, this test will be true...
            if( rcols[j] == otherrcols[j] )
                continue;
            
            //the order of the columns can be diferent...
            if( ! other.hasIndex(rind, rcols[j]) )
                return false;
        }
    }
    
    return true;
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::checkConsistency(std::ostream &out) const
{
    int code = 0;
    bool *auxCols = NULL;
    unsigned int nEls = 0;
    unsigned int nnzrows = 0;
    
    if( symmetric )
    {
        if( nrows != ncols )
        {
            out << SPM_PREPRINT "Error! symmetric matrix having nrows different to ncols. nrows: " << nrows << " ncols: " << ncols << SPM_GETFILELINE << std::endl; 
            code = SPM_BAD_DEFINITIONS;
        }
    }
    
    
    auxCols = (bool *) malloc( ncols * sizeof(bool) );
    if( !auxCols )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTMEMERROR;
        #endif
        
        code = SPM_MEMORY_ERROR;
        goto termination;
    }
    
    
    if( nElements > 0 && (baseStructure.cols == NULL || baseStructure.values == NULL) )
    {
        out << SPM_PREPRINT "Error! Row:s an incorrect baseStructure pointer." << SPM_GETFILELINE << std::endl;
        code = SPM_BAD_DEFINITIONS;
    }
    
    
    if( offset )
    {
        if( nrows > 0 )
        {
            if( offset[0] != 0 )
            {
                out << SPM_PREPRINT "Error! first row offset is not zero. offset: " << offset[0] << SPM_GETFILELINE << std::endl;
                code = SPM_BAD_DEFINITIONS;
            }
        }
        
        for( decltype(nrows) i = 0; i < nrows; i++ )
        {
            
            if( offset[i+1] - offset[i] != getNumberOfElementsAtRow(i) )  //( rows[i].__offset - rows[i-1].__offset != rows[i].nElements )
            {
                out << SPM_PREPRINT " on row " << i << " inconsistency about offset and number of elements. nElements: " << getNumberOfElementsAtRow(i) << " __offset: " << offset[i] << " previous rows __offset: " << offset[i-1] << SPM_GETFILELINE << std::endl;
                
                code = SPM_BAD_DEFINITIONS;
            }
            
            
            SPM_setAllArray( (symmetric ? i+1 : ncols), auxCols, false );
            
            
            const auto rnEls = getNumberOfElementsAtRow(i); //rows[i].nElements;
            
            if( rnEls > ncols )
            {
                out << SPM_PREPRINT "Error! Number of elements in row " << i << " has a number of elements greater than number of columns. number of elements: " << rnEls << " number of columns: " << ncols << SPM_GETFILELINE << std::endl;
                
                code = SPM_BAD_DEFINITIONS;
            }
            
            
            //matrixType* const rvalues = getRowValuesPointer(i); //rows[i].getValuesPointer();
            indexType* const rcols = getRowColsPointer(i); //rows[i].getColsPointer();
            
            for(unsigned int j = 0; j < rnEls; j++)
            {
                const unsigned int col = rcols[j];
                
                if( col >= ncols || col < 0 )
                {
                    out << SPM_PREPRINT " Error! row: " << i << " has an invalid col. col: " << col << " ncols: " << ncols << " on index: " << j << SPM_GETFILELINE << std::endl;
                    
                    code = SPM_BAD_DEFINITIONS;
                }
                else if( symmetric && col > i )
                {
                    out << SPM_PREPRINT "Error! symmetric matrix having a non inf triangular nonzero. row: " << i << " col: "<< col << SPM_GETFILELINE << std::endl;
                    
                    code = SPM_BAD_DEFINITIONS;
                }
                else
                {
                    if( auxCols[col] )
                    {
                        out << SPM_PREPRINT "Error! row: " << i << " specify column " << col << "two or more times\n";
                        //rows[i].print(out);
                        out << SPM_GETFILELINE << std::endl;
                        
                        SPM_getchar();
                        
                        code = SPM_BAD_DEFINITIONS;
                    }
                    
                    auxCols[col] = true;
                }
                
            }
            
            if( rnEls > 0 )
            {
                nEls += rnEls;
                nnzrows++;
            }
        }
        
        if( offset[nrows] != nElements )
        {
            out << SPM_PREPRINT "Error! Inconsistency about offset[nrows]: " << offset[nrows] << " nElements: " << nElements << SPM_GETFILELINE << std::endl;
            
            code = SPM_BAD_DEFINITIONS;
        }
    }
    
    
    if( nnzrows != this->nNzRowIndices )
    {
        out << SPM_PREPRINT "Error! Inconsistency about nNzRowIndices and real number of rows having nonzeros. nNzRowIndices: " << nNzRowIndices << " real number of rows having nonzeros: " << nnzrows << SPM_GETFILELINE << std::endl;
        
        code = SPM_BAD_DEFINITIONS;
    }
    
    if( nzRowIndex )
    {
        if( nzRowIndex[nNzRowIndices] < nrows )
        {
            out << SPM_PREPRINT "Error! End of nNzRowIndices is not set as nrows or greater value!. auxNzRowIndex[nNzRowIndices]: " << nzRowIndex[nNzRowIndices] << " nrows: " << nrows << SPM_GETFILELINE << std::endl;
            code = SPM_BAD_DEFINITIONS;
        }
    }
    
    
    for( unsigned int i = 0; i < nNzRowIndices; i++ )
    {
        const auto rind = nzRowIndex[i];
        
        if( getNumberOfElementsAtRow(rind) == 0 ) //( rows[rind].nElements == 0 )
        {
            out << SPM_PREPRINT "Error! Row " << rind << " is counted in nzRowIndex, but has no elemtens." << SPM_GETFILELINE << std::endl;
            
            code = SPM_BAD_DEFINITIONS;
        }
        
        if( i > 0 )
        {
            if( nzRowIndex[i] <= nzRowIndex[i-1] )
            {
                out << SPM_PREPRINT "Error! nzRowIndex is not ordered! nzRowIndex["<<i<<"]: " << nzRowIndex[i] << " nzRowIndex["<<i-1<<"]: " << nzRowIndex[i-1] << SPM_GETFILELINE << std::endl;
                code = SPM_BAD_DEFINITIONS;
            }
        }
    }
    
    
    if( nEls != this->nElements )
    {
        printSparseMatrix();
        
        out << SPM_PREPRINT "Error! Inconsistency about nELements and number of elements in rows!. nElements: " << this->nElements << " number of elements in rows: " << nEls << SPM_GETFILELINE << std::endl;
        
        code = SPM_BAD_DEFINITIONS;
    }
    
    
    
termination:
    
    if(auxCols)	free(auxCols);
    
    
    return code;
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>:: checkNElements() const
{
    int code = 0;
    unsigned int nEl = 0;
    unsigned int nElnz = 0;
    
    
    for(unsigned int i = 0; i < nrows; i++)
        nEl += getNumberOfElementsAtRow(i);
    
    
    if( nEl != nElements )
    {
        #if SPM_DEBUG_MODE
            std::cout << SPM_PREPRINT << "Inconsistency error at number of elements. nElements: " << nElements << " real number of elements: " << nEl << SPM_GETFILELINE << "\n";
        #endif
        
        code = SPM_INTERNAL_ERROR;
    }
    
    //now, we check nzRowIndex
    for(unsigned int i = 0; i < nNzRowIndices; i++)
        nElnz += getNumberOfElementsAtRow(nzRowIndex[i]);  //rows[nzRowIndex[i]].nElements;
    
    if( nElnz != nEl )
    {
        #if SPM_DEBUG_MODE
            std::cout << SPM_PREPRINT << "Inconsistency error at number of elements at nzRowIndex. nElements: " << nElnz << " real number of elements: " << nEl << SPM_GETFILELINE << "\n";
        #endif
        
        code = SPM_INTERNAL_ERROR;
    }
    
    return code;
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::copyParametersFrom(const SPM_NewSparseMatrix<indexType, matrixType>& M)
{
    symmetric = M.symmetric;
    nrows = M.nrows;
    ncols = M.ncols;
    nElements = M.nElements;
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::copyRowTo(const unsigned int rowIndex, matrixType* row, const matrixType factor) const
{
    const auto rnEls = getNumberOfElementsAtRow(rowIndex); 
    
    matrixType* const values = getRowValuesPointer(rowIndex);
    indexType* const cols = getRowColsPointer(rowIndex);
    
    
    if( factor == 1 )
    {
        for(unsigned int i = 0; i < rnEls; i++)
            row[ cols[i] ] = values[i];
    }
    else
    {
        for(unsigned int i = 0; i < rnEls; i++)
            row[ cols[i] ] = factor*values[i];
    }
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::copyStructureFrom(const SPM_NewSparseMatrix<indexType, matrixType>& M)
{
    int r, code;
    
    if( nrows != M.nrows || M.nrows == 0 )
    {
        desallocateMemory();
        
        const int r = allocateSparseRows(M.nrows);
        if(r != 0)
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTERRORNUMBER(r);
            #endif
            code = r;
            goto termination;
        }
    }
    
    r = reallocateColsAndValues(M.nElements);
    if( r != 0 )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        code = r;
        goto termination;
    }
    
    copyParametersFrom(M);
    
    /*nElements = 0;
    for(decltype(nrows) i = 0; i < nrows; i++)
    {
        offset[i] = nElements;
        nElements += M.getNumberOfElementsAtRow(i);
    }
    //setting the end of offset
    offset[nrows] = nElements; 
    #if SPM_DEBUG_MODE
        assert(nElements == M.nElements);
    #endif
    */
    
    if (M.nrows == 0)
    {
        offset[0] = 0;
    }
    else
    {
        SPM_copyArray(nrows+1, M.offset, offset);
    }
    
    nElements = M.nElements;
    
    //for(int i = 0; i <= nrows; i++)
        //printf("offset[%d]: %d\n", i, offset[i]);
    
    
    SPM_copyArray(nElements, M.baseStructure.cols, baseStructure.cols);
    
    nNzRowIndices = M.nNzRowIndices;
    if( M.nzRowIndex )
    {
        SPM_copyArray(nNzRowIndices+1, M.nzRowIndex, nzRowIndex);
    }
    
    
    #if SPM_DEBUG_MODE
        r = checkNElements();
        if(r != 0)
        {
            SPM_PRINTERRORMSG(r);
            code = r;
            goto termination;
        }
    #endif
    
    code = 0;
    
termination:
    
    if(code != 0)
        desallocateMemory();
    
    return code;
}

template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::copyMatrixFrom(const SPM_NewSparseMatrix<indexType, matrixType> &M)
{
    const int r = copyStructureFrom(M);
    if( r != 0 )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    SPM_copyArray(nElements, M.baseStructure.values, baseStructure.values);
    
    //nzRowIndex is already copied in copyStructureFrom
    
    return 0;
}

template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::copyMatrixTo(matrixType **M, const bool considerSymmetry, const bool initializeWithZero, const double factor ) const
{
    if( initializeWithZero )
    {
        for(unsigned int i = 0; i < nrows; i++)
            SPM_setAllArray<matrixType>(ncols, M[i], 0);
    }
    
    if(factor == 0.0)
        return;
    
    for(unsigned int i = 0; i < nNzRowIndices; i++)
    {
        const unsigned int rind = nzRowIndex[i]; 
        copyRowTo( rind, M[rind], factor ); //rows[rowIndex].copyTo( M[rowIndex], factor );
    }
    
    if( considerSymmetry && symmetric )
    {
        for(unsigned int i = 0; i < nNzRowIndices; i++)
        {
            const unsigned int rind = nzRowIndex[i];
            
            const auto rnels = getNumberOfElementsAtRow(rind); //rows[rowIndex].nElements;
            
            indexType* const rcols = getRowColsPointer(rind); //rows[rowIndex].getColsPointer();
            matrixType* const rvalues = getRowValuesPointer(rind); //rows[rowIndex].getValuesPointer();
            
            for(unsigned int j = 0; j < rnels; j++)
                M[ rcols[j] ][rind] = rvalues[j];
            
            if( factor != 1 )
            {
                for(unsigned int j = 0; j < rnels; j++)
                    M[ rcols[j] ][rind] *= factor;
            }
        }
        
        #if SPM_DEBUG_MODE
        for(decltype(nrows) i = 1; i < nrows; i++)
        {
            for(unsigned int j = 0; j < i; j++)
                assert( M[i][j] == M[j][i] );
        }
        #endif
    }
    
}

//if Matrix == NULL, the method allocate an array  to store the coefficients and return. Otherwise, the method puts the coefficient in Matrix and return it
template <class indexType, class matrixType>
matrixType* SPM_NewSparseMatrix<indexType, matrixType>:: copyMatrixTo(matrixType* Matrix, const bool considerSymmetry, const bool initializeWithZero, const double factor) const
{
    matrixType *M;
    
    if(Matrix)
    {
        M = Matrix;
        
        if( initializeWithZero )
            SPM_setAllArray<matrixType>(nrows*ncols, M, 0);
    }
    else
    {
        M = (matrixType *) calloc( nrows*ncols, sizeof(matrixType) );
        if( !M )
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTMEMERROR;
            #endif
            
            return NULL;
        }
        
        //we already initialize with zeros
    }
    
    for(decltype(nNzRowIndices) i = 0; i < nNzRowIndices; i++)
    {
        const unsigned int rind = nzRowIndex[i];
        copyRowTo( rind, &M[rind*ncols], factor );
        //rows[rowIndex].copyTo( &M[rowIndex*ncols], factor );
    }
    
    if( considerSymmetry && symmetric )
    {
        for(decltype(nNzRowIndices) i = 0; i < nNzRowIndices; i++)
        {
            const unsigned int rind = nzRowIndex[i];
            
            const unsigned int rnels = getNumberOfElementsAtRow(rind); //rows[rowIndex].nElements;
            
            indexType* const rcols = getRowColsPointer(rind); //rows[rowIndex].getColsPointer();
            matrixType* const rvalues = getRowValuesPointer(rind); //rows[rowIndex].getValuesPointer();
            
            for(unsigned int j = 0; j < rnels; j++)
                M[rcols[j]*ncols + rind] = rvalues[j];
            
            if( factor != 1 )
            {
                for(unsigned int j = 0; j < rnels; j++)
                    M[rcols[j]*ncols + rind] *= factor;
            }
        }
        
        #if SPM_DEBUG_MODE
            for(decltype(nrows) i = 1; i < nrows; i++)
            {
                for(unsigned int j = 0; j < i; j++)
                    assert( M[i*ncols + j] == M[j*ncols + i] );
            }
        #endif
    }
    
    return M;
}




template <class indexType, class matrixType>
template <class intClass>
void SPM_NewSparseMatrix<indexType, matrixType>::__countRowsEachColumn(intClass *counts, const bool accumulate) const
{
    if( !accumulate )
        SPM_setAllArray<intClass>(ncols, counts, 0);
    
    for(decltype(nNzRowIndices) i = 0; i < nNzRowIndices; i++)
    {
        const auto rind = nzRowIndex[i];
        
        const auto rnel = getNumberOfElementsAtRow(rind); // rows[rowIndex].nElements;
        
        indexType* const rcols = getRowColsPointer(rind);
        
        for(unsigned int j = 0; j < rnel; j++)
            counts[ rcols[j] ]++;
    }
}


//that method count, for each column, how many rows exists where the respective column appears in sparse matrix
template <class indexType, class matrixType>
void newspm:: SPM_NewSparseMatrix<indexType, matrixType>::countRowsEachColumn(int *counts, const bool accumulate) const
{
    return __countRowsEachColumn(counts, accumulate);
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::countRowsEachColumn(unsigned int *counts, const bool accumulate) const
{
    return __countRowsEachColumn(counts, accumulate);
}

template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::deleteRowStructure(const unsigned int rowIndex)
{
    if(rowIndex >= nrows)
        return SPM_INDEX_FAULT;
    
    const auto rnel = getNumberOfElementsAtRow(rowIndex); //rows[row].nElements;
    
    if( rnel > 0 )
    {
        #if SPM_DEBUG_MODE
            const auto nnzr = nNzRowIndices;
        #endif
        
        //rows[row].desallocate();
        
        //shiftMinusRowOffset( rowIndex + 1, rnel);
        leftShiftColsValuesAndOffset( rowIndex+1, rnel );
        
        nElements -= rnel;
        
        updateNzRowIndex();
        
        #if SPM_DEBUG_MODE
            assert(nnzr == nNzRowIndices-1);
        #endif
    }
    
    return 0;
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::deleteStructure()
{
    if( offset )
        SPM_setAllArray<decltype(nElements)>( nrows+1, offset, 0 );
    
    baseStructure.desallocateMemory();
    
    nElements = 0;
    nNzRowIndices = 0;
    if(nzRowIndex) //not necessarilly nzRowIndex is allocated
        nzRowIndex[0] = nrows;
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::desallocateMemory(void)
{
    deleteStructure();
    SPM_secFree(offset);
    SPM_secFree(nzRowIndex);
    nrows = 0;
    initialize(0, 0, symmetric); //we keep the current value of symmetric flag passing it to the initialize value...
}



template <class indexType, class matrixType>
double SPM_NewSparseMatrix<indexType, matrixType>::evalRowTimesxt(const indexType rowIndex, const double *x) const
{
    const auto rnels = getNumberOfElementsAtRow(rowIndex);
    matrixType* const values = getRowValuesPointer(rowIndex);
    indexType* const cols = getRowColsPointer(rowIndex);
    
    double v = 0.0;
    
    for(unsigned int i = 0; i < rnels; i++)
        v += values[i]*x[cols[i]];
    
    return v;
}



template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::evalTimesxt(const double *x, double *v, const bool accumulate, const double factor) const
{
    if(!accumulate)
        SPM_setAllArray<double>(nrows, v, 0.0);
    
    
    if( factor == 0.0 )
        return;
    
    
    if(symmetric)
    {
        //we have to run elements 
        
        for(decltype(nNzRowIndices) i = 0; i < nNzRowIndices; i++ )
        {
            const auto row = nzRowIndex[i];
            const auto rnElements = getNumberOfElementsAtRow(row);
            
            const indexType *rcols = getRowColsPointer(row);
            const matrixType *rvalues = getRowValuesPointer(row);
            
            double vrow = 0.0;
            
            const double factorxrow = factor*x[row];
            
            for( unsigned int j = 0; j < rnElements; j++ )
            {
                const auto col = rcols[j];
                const auto value = rvalues[j];
                
                 vrow += value*x[col]; //we consider factor after the sum. Note, v can be being accumulated with previous values. So, do not save the sum in v[row] directly here.
                
                if(row != col)
                { //nondiagonal position. Since this matrix is symetric, we must consider the upper triangle also.
                    
                    v[col] += factorxrow*value; //we cannot remove factor multiplying from here because v can be being accumulated with previous values. 
                }
            }
            
            v[row] += factor*vrow;
        }
        
    }
    else
    {
        for(decltype(nNzRowIndices) i = 0; i < nNzRowIndices; i++)
        {
            const auto rowIndex = nzRowIndex[i];
            
            v[rowIndex] += factor* evalRowTimesxt(rowIndex, x);
        }
    }
}


//eval factor x*Matrix
template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::evalPretimesx(const double *x, double *v, const bool accumulate, const double factor) const
{
    return evalPretimesx(nNzRowIndices, nzRowIndex, x, v, accumulate, factor);
}


//eval factor x*Matrix only in some rows
template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::evalPretimesx(const unsigned int nRowInds, const unsigned *rowInds, const double *x, double *v, const bool accumulate, const double factor) const
{
    if(!accumulate)
        SPM_setAllArray<double>(ncols, v, 0.0);
    
    if( factor == 0.0 )
        return;
    
    
    for( unsigned int i = 0; i < nRowInds; i++ )
    {
        const auto row = rowInds[i];
        const auto rnElements = getNumberOfElementsAtRow(row);
        
        const indexType *rcols = getRowColsPointer(row);
        const matrixType *rvalues = getRowValuesPointer(row);
        
        
        const double factorxrow = factor * x[row];
        
        
        for( unsigned int j = 0; j < rnElements; j++ )
        {
            const auto col = rcols[j];
            const auto value = rvalues[j];
            
            
            v[col] += factorxrow * value ; //here, we should run matrix by columns to remove factor multiplying by here, but it is not possible
        }
        
        
        if( symmetric )
        {
            double vrow = 0.0;
            
            for( unsigned int j = 0; j < rnElements; j++ )
            {
                const auto col = rcols[j];
                
                if( col != row )
                { //nondiagonal position. Since this matrix is symetric, we must consider the upper triangle also.
                    const auto value = rvalues[j];
                    
                    vrow += x[col]*value;
                }
            }
            
            v[row] += factor*vrow;
        }
    }
    
    

    
    
}



template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::getElement(unsigned int rowIndex, indexType colIndex, matrixType &value) const
{
    if(symmetric && (unsigned int) colIndex > rowIndex)
        SPM_swap(rowIndex, colIndex);
    
    
    const auto rnels = getNumberOfElementsAtRow(rowIndex);
    matrixType* const values = getRowValuesPointer(rowIndex);
    indexType* const cols = getRowColsPointer(rowIndex);
    
    
    for(unsigned int i = 0; i < rnels; i++)
    {
        if( cols[i] == colIndex )
        {
            value = values[i];
            return 0;
        }
    }
    
    value = 0.0;
    return SPM_ELEMENT_NOT_PRESENT;
    //return rows[row].getElement(colIndex, value);
}


//that method fill all ncols in vector values
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::getFullRow(const unsigned int rowIndex, matrixType *values, const bool accumulate, const bool considerSymmetry, const double factor ) const
{
    if( rowIndex >= nrows )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTINDEXERROR;
        #endif
        return SPM_INDEX_FAULT;
    }
    
    
    if( accumulate )
    {
        accumulateRowInArray(rowIndex, values, factor); //rows[index].copyTo(values, factor);
    }
    else
    {
        SPM_setAllArray<matrixType>( ncols, values, 0 );
        copyRowTo(rowIndex, values, factor);
    }
    
    
    if( considerSymmetry && symmetric && nNzRowIndices > 0 )
    {
        //since we only store lower triangle, we need find the positions would be in the upper one. Note those positions are below row index
        
        //this work because we are sure nNzRowIndices is greater than  0
        for( decltype(nNzRowIndices) i = nNzRowIndices-1; ; i--)
        {
            const unsigned int rind = nzRowIndex[i];
            
            if( rind <= rowIndex )
                break; //we assume indices are in ascendet order in nzRowIndex, and we are running from the end
            
            //if index is the row, we set values on rowIndex colunm
            
            matrixType v;
            
            if( hasIndex(rind, rowIndex, &v, NULL) )
                values[rind] += factor*v;
            
            #if SPM_DEBUG_MODE
                assert(i != 0); // we cannot reach this point if i is equal to zero. In this, we would run all nzRowIndex and would have found the row rowIndex.
            //if(i == 0)
                //break; //we reach the end
            #endif
        }
    }
    
    return 0;
}




template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::getRowStructure(const unsigned int rowIndex, indexType* columns, indexType *nzs) const
{
    if( rowIndex >= nrows )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTINDEXERROR;
        #endif
        return SPM_INDEX_FAULT;
    }
    
    const auto mynzs = getNumberOfElementsAtRow(rowIndex); //rows[index].getStructure(cols);
    
    
    SPM_copyArray(mynzs, getRowColsPointer(rowIndex), columns);
    
    if(nzs)
        *nzs = mynzs;
    
    return 0;
}


//that method put true in cols[i] if i is a column in line. If accumulate is true, method does not initialize cols with false.
//return the number of elements in Row...
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::getRowStructure(const unsigned int rowIndex, bool* columns, indexType *nzs, bool accumulate) const
{
    if( rowIndex >= nrows )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTINDEXERROR;
        #endif
        return SPM_INDEX_FAULT;
    }
    
    if(!accumulate)
        SPM_setAllArray(ncols, columns, false);
    
    const auto mynzs = getNumberOfElementsAtRow(rowIndex); //rows[index].getStructure(cols);
    
    indexType* const cols = getRowColsPointer(rowIndex);
    
    #pragma ivdep
    #pragma GCC ivdep
    for(unsigned int i = 0; i < mynzs; i++)
        columns[ cols[i] ] = true;
    
    if(nzs)
        *nzs = mynzs;
    
    return 0;
}


/*template <class indexType, class matrixType>
int newspm:: SPM_NewSparseMatrix<indexType, matrixType>::getRowStructureAndValues(const unsigned int index, indexType *nzs, indexType *cols, matrixType *vals) const
{
    if( index >= nrows )
        return SPM_INDEX_FAULT;
    
    indexType mynzs = rows[index].getStructureAndValues(cols, vals);
    
    if(nzs)
        nzs = mynzs;
    
    return 0;
} */


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::getRowValues(const unsigned int rowIndex, matrixType *vals, indexType *nzs) const
{
    if( rowIndex >= nrows )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTINDEXERROR;
        #endif
        return SPM_INDEX_FAULT;
    }
    
    const auto mynzs = getNumberOfElementsAtRow(rowIndex); //rows[index].getValues(vals);
    
    SPM_copyArray( mynzs, getRowValuesPointer(rowIndex), vals );
    
    if(nzs)
        *nzs = mynzs;
    
    return 0;
}


template <class indexType, class matrixType>
template <class intType>
indexType SPM_NewSparseMatrix<indexType, matrixType>::__getStructure(intType* rows, indexType* cols) const
{
    if(rows)
    {
        decltype(nElements) k = 0;
        
        for(decltype(nNzRowIndices) i = 0; i < nNzRowIndices; i++)
        {
            const indexType rind = nzRowIndex[i];
            
            const indexType rnz =  getNumberOfElementsAtRow(rind);  //rows[rind].nElements;
            
            SPM_setAllArray<intType>( rnz, &rows[k], rind );
            
            k += rnz;
        }
        
        #if SPM_DEBUG_MODE
            assert(k == nElements);
        #endif
    }
    
    if(cols)
        SPM_copyArray( nElements, baseStructure.cols, cols );
    
    
    return nElements;
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::getStructure(int* rows, indexType* cols) const
{
    return __getStructure(rows, cols);
}

template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::getStructure(unsigned int* rows, indexType* cols) const
{
    return __getStructure(rows, cols);
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::getValues(matrixType *values) const
{
    SPM_copyArray(nElements, baseStructure.values, values);
    
    return nElements;
}


/*template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::getStructureAndValues(int* iRow, indexType* jCol, matrixType* vals) const
{
}
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::getStructureAndValues(unsigned int* iRow, indexType* jCol, matrixType* vals) const
{
} */



template <class indexType, class matrixType>
template <class intType>
void SPM_NewSparseMatrix<indexType, matrixType>::__getStructureByCompressedRowFormat( intType* rowStart, indexType* cols, matrixType* values) const
{
    if(rowStart)
        SPM_copyArray(nrows+1, offset, rowStart);
    
    if(cols)
        SPM_copyArray(nElements, baseStructure.cols, cols);
    
    if(values)
        SPM_copyArray(nElements, baseStructure.values, values);
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::getStructureByCompressedRowFormat(int *rowStart, indexType *cols, matrixType* values) const
{
    return __getStructureByCompressedRowFormat(rowStart, cols, values);
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::getStructureByCompressedRowFormat(unsigned int *rowStart, indexType *cols, matrixType* values) const
{
    return __getStructureByCompressedRowFormat(rowStart, cols, values);
}


template <class indexType, class matrixType>
bool SPM_NewSparseMatrix<indexType, matrixType>::hasIndex(const unsigned int rowIndex, const indexType col, matrixType *value, indexType *indexAtRowStructure) const
{
    if( rowIndex >= nrows )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTINDEXERROR;
        #endif
        return false;
    }
    
    const auto rnels = getNumberOfElementsAtRow(rowIndex);
    indexType* const rcols = getRowColsPointer(rowIndex);
    
    
    for(unsigned int i = 0; i < rnels; i++)
    {
        if( rcols[i] == col )
        {
            if(indexAtRowStructure)
                *indexAtRowStructure = i;
            if(value)
                *value = getRowValuesPointer(rowIndex)[i];
            return true;
        }
    }
    
    return false;
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::initialize(const unsigned int nrows, const indexType ncols, const bool symmetric, const bool initialize)
{
    offset = NULL;
    nElements = 0;
    this->nrows = 0; //do not remove it. We need initialize this->nrows as 0 and after put as nrows if we alloc the memory to set the offsets correctly
    
    nzRowIndex = NULL;
    nNzRowIndices = 0;
    
    
    this->symmetric = symmetric;
    this->ncols =ncols;
    
    
    baseStructure.initialize();
    
    if(initialize && nrows > 0)
    {
        const int r = allocateSparseRows(nrows);
        if( r != 0 )
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
    
        this->nrows = nrows;
    }
    
    
    
    return 0;
}

template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::mergeStructures(SPM_NewSparseMatrix<indexType, matrixType>& M, const bool storeOrigCoefs)
{
    //const unsigned int norows = nrows;
    //const unsigned int minrows = SPM_min( norows, M.nrows );
    const unsigned int maxcols = SPM_max(ncols, M.ncols);
    
    
    int code;
    matrixType *auxValues = NULL;
    
    if( M.nrows > nrows )
    {
        const int r = addNewRows(M.nrows - nrows);
        
        if(r != 0)
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
        
        #if SPM_DEBUG_MODE
            assert(M.nrows == nrows);
        #endif
    }
    
    
    auxValues = (matrixType *) malloc( maxcols * sizeof(matrixType) );
    if( !auxValues )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTMEMERROR;
        #endif
        code = SPM_MEMORY_ERROR;
        goto termination;
    }
    
    
    
    for(decltype(nNzRowIndices) i = 0; i < M.nNzRowIndices; i++)
    {
        const auto rowIndex = M.nzRowIndex[i];
        
        
        if( getNumberOfElementsAtRow(rowIndex) == 0) //( rows[rowIndex].nElements == 0 )
        {
            const int r = setRowStructureAndValues(rowIndex, rowIndex, M);
            
            if(r != 0)
            {
                #if SPM_DEBUG_MODE
                    SPM_PRINTERRORNUMBER(r);
                #endif
                code = r;
                goto termination;
            }
        }
        else
        {
            const decltype(ncols) cols = symmetric && M.symmetric ? rowIndex+1 : maxcols;
            //SPM_NewSparseRow<indexType, matrixType> *firstrow, *secondrow;
            
            
            SPM_setAllArray<matrixType>(cols, auxValues, 0);
            
            
            if(storeOrigCoefs)
            {
                //firstrow = &(M[rowIndex]);
                //secondrow = &(rows[rowIndex]);
                
                M.copyRowTo(rowIndex, auxValues, 1);
                copyRowTo(rowIndex, auxValues, 1);
            }
            else
            {
                //firstrow = &(rows[rowIndex]);
                //secondrow = &(M[rowIndex]);
                
                copyRowTo(rowIndex, auxValues, 1);
                M.copyRowTo(rowIndex, auxValues, 1);
            }
            
            
            //firstrow->copyTo(auxValues, 1);
            //secondrow->copyTo(auxValues, 1);
            
            
            //M[rowIndex].getStructure(auxCols);
            //rows[rowIndex].getStructure(auxCols);
            
            #if SPM_DEBUG_MODE
                unsigned int nz = 0;
                
                for(decltype(cols) j = 0; j < cols; j++)
                {
                    if( auxValues[j] != 0.0)
                        nz++;
                }
                
                assert(nz >= getNumberOfElementsAtRow(rowIndex)); //assert(nz >= rows[rowIndex].nElements);
                assert(nz >= M.getNumberOfElementsAtRow(rowIndex)); //assert(nz >= M.rows[rowIndex].nElements);
            #endif
            
            
            const int r = setRowStructureAndValues(rowIndex, auxValues, cols, 0.0);
            
            if( r != 0 )
            {
                #if SPM_DEBUG_MODE
                    SPM_PRINTERRORMSG(r);
                #endif
                code = r;
                goto termination;
            }
        }
        
    }
    
    
    if( M.ncols > ncols )
        ncols = M.ncols;
    
    if( !M.symmetric )
        symmetric = false;
    
    
    code = 0;
    
termination:
    
    updateNzRowIndex();
    
    if(auxValues)	free(auxValues);
    
    return code;
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::multiplyAllElements( const matrixType value )
{
    SPM_multiplyAllArray(nElements,  baseStructure.values, value);
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::multiplyAllElementsAtRow(const unsigned int rowIndex, const matrixType value)
{
    if( rowIndex >= nrows )
        return SPM_INDEX_FAULT;
    SPM_multiplyAllArray(getNumberOfElementsAtRow(rowIndex), getRowValuesPointer(rowIndex), value);
    
    return 0;
}


template <class indexType, class matrixType>
void  SPM_NewSparseMatrix<indexType, matrixType>::printSparseMatrix(const bool showEmptyLines) const
{
    std::cout << "SPM_SparseMatrix:: nrows: " << nrows << " ncols: " << ncols << " nElements: " << nElements << " symmetric: " << symmetric << " nonzero rows: " << nNzRowIndices << "\n";
    
    
    const unsigned int maxi = showEmptyLines ? nrows : nNzRowIndices ;
    
    for(unsigned int i = 0; i < maxi; i++)
    {
        const decltype(nrows) ind = showEmptyLines ? i : nzRowIndex[i];
        
        indexType* const rcols = getRowColsPointer(ind);
        matrixType* const rvalues = getRowValuesPointer(ind);
        const auto nel = getNumberOfElementsAtRow(ind); //rows[ind].getNumberOfElements();
        
        std::cout << ind << "=> ";
        
        for(unsigned int j = 0; j < nel; j++)
            std::cout << rcols[j] << ": " << rvalues[j] << "  "; //std::cout << rows[ind][j] << ": " << rows[ind](j) << " ";
        
        std::cout << "\n";
    }
    
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::printAllMatrix(void) const
{
    matrixType c;
    
    for(unsigned int i = 0; i < nrows; i++)
    {
        for(unsigned int j = 0; j < ncols; j++)
        {
            getElement(i, j, c);
            std::cout << c << ", ";
        }
        std::cout << ";\n";
    }
}


//That function calculate 0.5* x'M x to a symmetric matrix M
template <class indexType, class matrixType>
double SPM_NewSparseMatrix<indexType, matrixType>::quadraticEvaluation(const double* x) const
{
    double eval = 0.0;
    
    #if SPM_DEBUG_MODE
        assert(symmetric);
        assert(nrows == ncols);
    #endif
    
    for(unsigned int i = 0; i < nNzRowIndices; i++)
    {
        double reval = 0.0; //we accumulate the sum in a separated variable for each row to try decrease numerical error...
        
        const auto rind = nzRowIndex[i];
        
        const double xrind = x[rind];
        
        
        if( SPM_abs(xrind) <= SPM_MY_ZERO_TO_EVAL)
            continue; //this row will be turned zero by xrind 
            
        
        const auto nel = getNumberOfElementsAtRow(rind); //rows[rind].nElements;
        
        indexType* const rcols = getRowColsPointer(rind); //rows[rind].getColsPointer();
        matrixType* const rvalues = getRowValuesPointer(rind); // rows[rind].getValuesPointer();
        
        for(unsigned int j = 0; j < nel; j++)
        {
            const unsigned k = rcols[j];
            
            const double v = x[k] * rvalues[j]; //we multiply by reval at the end of evaluation
            
            if(rind == k)
                reval += 0.5*v;
            else
                reval += v;
        }
        
        eval += xrind *reval;
    }
    
    return eval;
}


/* Function to calculate the gradient of 0.5* x'Mx in grad vector. We suppose M is symmetric
    * 
    *Warning: If accumulate == true, That function only accumulates the gradient in grad, i.e., you should initialize grad vector by yourself. */
template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::quadraticGradientEvaluation(const double* x, double* grad, const bool accumulate) const
{
    #if SPM_DEBUG_MODE
        assert(symmetric);
        assert(nrows == ncols);
    #endif
    
    if( !accumulate )
        SPM_setAllArray<double>(nrows, grad, 0.0);
    
    for(decltype(nNzRowIndices) i = 0; i < nNzRowIndices; i++)
    {
        const auto rind = nzRowIndex[i];
        
        indexType* const rcols = getRowColsPointer(rind); //rows[rind].getColsPointer();
        
        matrixType* const rvalues = getRowValuesPointer(rind); //rows[rind].getValuesPointer();
        
        const auto rnel = getNumberOfElementsAtRow(rind); //rows[rind].nElements;
        
        for(unsigned int j = 0; j < rnel; j++)
        {
            const decltype(rind) k = rcols[j];
            
            grad[rind] += x[k] * rvalues[j];
            
            if(rind != k) //nondiagonal positions
                grad[k] += x[rind] * rvalues[j];
        }
    }
}


/* Function to calculate 0.5* x'M x and the gradient of 0.5* x'Mx in grad vector. We suppose M is symmetric
    * 
    *Warning: If accumulate == true, That function only accumulates the gradient in grad, i.e., you should initialize grad vector by yourself. */
template <class indexType, class matrixType>
double  SPM_NewSparseMatrix<indexType, matrixType>::quadraticEvaluationAndGradient( const double* x, double* grad, const bool accumulate) const
{
    double eval = 0.0;
    
    #if SPM_DEBUG_MODE
        assert(symmetric);
        assert(nrows == ncols);
    #endif
    
    if( !accumulate )
        SPM_setAllArray<double>(nrows, grad, 0);
    
    
    for(unsigned int i = 0; i < nNzRowIndices; i++)
    {
        double reval = 0.0; //we accumulate the sum in a separated variable for each row to try decrease numerical error...
        
        const unsigned int rind = nzRowIndex[i];
        
        const auto rnel = getNumberOfElementsAtRow(rind); //rows[rind].nElements;
        indexType* const rcols = getRowColsPointer(rind); //rows[rind].getColsPointer();
        matrixType* const rvalues = getRowValuesPointer(rind); //rows[rind].getValuesPointer();
        
        for(unsigned int j = 0; j < rnel; j++)
        {
            const unsigned int k = rcols[j];
            
            const double v1 = x[k] *rvalues[j];
            grad[rind] += v1;
            
            if(rind == k)
            {
                reval += 0.5 * v1;
            }
            else //nondiagonal positions
            {
                reval += v1;
                grad[k] += x[rind] * rvalues[j];
            }
        }
        
        eval += x[rind]*reval;
    }
    
    return eval;
}


template <class indexType, class matrixType>
int  SPM_NewSparseMatrix<indexType, matrixType>::reallocateRowSpace(const unsigned int rowIndex, const unsigned int newsize)
{
    const unsigned int oldsize = getNumberOfElementsAtRow(rowIndex);
    
    /*std::cout << "offset[0]: " << offset[0] << std::endl;
    std::cout << "offset[1]: " << offset[1] << std::endl;
    
    std::cout << "rowIndex: " << rowIndex << std::endl;
    std::cout << "newsize: " << newsize << std::endl;
    std::cout << "oldsize: " << oldsize << std::endl;
    
    SPM_getchar(); */
    
    if( newsize > oldsize )
    {
        const unsigned int shift = newsize - oldsize;
        const unsigned int newnnels = nElements + shift;
        
        const int r = reallocateColsAndValues( newnnels );
        if( r != 0 )
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
        
        rightShiftColsValuesAndOffset(rowIndex+1, shift);
        
        nElements = newnnels;
    }
    else if( oldsize > newsize )
    {
        const unsigned int shift = oldsize - newsize;
        const unsigned int newnnels = nElements - shift;
        
        leftShiftColsValuesAndOffset(rowIndex+1, shift);
        
        const int r = reallocateColsAndValues( newnnels );
        if( r != 0 )
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
        
        nElements = newnnels;
    }
    
    
    return 0;
}


template <class indexType, class matrixType>
int  SPM_NewSparseMatrix<indexType, matrixType>::reallocateRowStructures(const unsigned int newNumberOfRows)
{
    int code = 0;
    
    //if( newNumberOfRows > 0 )
    {
        decltype(offset) auxOffset;
        decltype(nzRowIndex) auxNzRowIndex;
        
        auxOffset = (decltype(auxOffset)) realloc( offset, (newNumberOfRows+1) * sizeof(offset[0]) );
        
        if( auxOffset )
            offset = auxOffset;
        else
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTMEMERROR;
            #endif
            code = SPM_MEMORY_ERROR;
        }
        
        auxNzRowIndex = (decltype(auxNzRowIndex)) realloc( nzRowIndex, (newNumberOfRows+1) * sizeof(nzRowIndex[0]) ); //we put one position more to sinalize last position with nrows value. It can be useful for iterators.
        
        if( auxNzRowIndex )
            nzRowIndex = auxNzRowIndex;
        else
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTMEMERROR;
            #endif
            code = SPM_MEMORY_ERROR;
        }
    }
    /*else
    {
        SPM_secFree(offset);
        SPM_secFree(nzRowIndex);
    }*/
    
    return code;
}
    



template <class indexType, class matrixType>
int  SPM_NewSparseMatrix<indexType, matrixType>::removeColsCoefficients( const unsigned int ncols, const indexType *cols )
{
    bool *fcols = NULL;
    
    fcols = (bool*) calloc( this->ncols, sizeof(bool) );
    if( !fcols )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTMEMERROR;
        #endif
        return SPM_MEMORY_ERROR;
    }
    
    #pragma ivdep
    #pragma GCC ivdep
    for(unsigned int i = 0; i < ncols; i++)
        fcols[ cols[i] ] = true;
    
    
    removeColsCoefficients(fcols);
    
    
    free(fcols);
    
    return 0;
}


//we remove the columns marked as true. Even if matrix is simmetric, we do not remove lines. We do not change ncols value. In the practice, is like you remove all coeficients in the marked columns
template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::removeColsCoefficients( const bool *cols)
{
    unsigned int firstCol = ncols;
    unsigned int totalRemoved = 0;
    
    if(nNzRowIndices == 0)
        return;
    
    for(unsigned int i = 0; i < ncols; i++)
    {
        if( cols[i] )
        {
            firstCol = i;
            break;
        }
    }
    
    //if firstCol == ncols, is because all cols array is false
    if(firstCol == ncols)
        return;
    
    
    for(auto i = nNzRowIndices-1; ; i--)
    {
        const auto rind = nzRowIndex[i];
        
        if( symmetric && firstCol > rind )
            continue; //we just store lower triangle
        
        unsigned int nremoved = 0;
        
        const auto rnel = getNumberOfElementsAtRow(rind); // rows[rind].nElements;
        
        
        if( rnel == 0 )
        {
            #if SPM_DEBUG_MODE
                assert(false); //this column should have some element...
            #endif
            continue;
        }
        
        indexType* const rcols = getRowColsPointer(rind); // rows[rind].getColsPointer();
        matrixType* const rvalues = getRowValuesPointer(rind); // rows[rind].getValuesPointer();
        
        for(auto j = rnel-1; ; j--)
        {
            if( cols[rcols[j]] )
            {
                const unsigned int size = rnel-j-1 -nremoved;
                
                #if SPM_DEBUG_MODE
                    assert(size < rnel); //if size get negative, it will be a larger number because it is a unsigned int
                #endif
                
                SPM_shiftLeftArray(size, &rcols[j+1]);
                SPM_shiftLeftArray(size, &rvalues[j+1]);
                
                nremoved++;
            }
            
            if(j == 0)
                break;
        }
        
        #if SPM_DEBUG_MODE
            assert(nremoved <= rnel);
        #endif
        
        //rows[rind].nElements -= nremoved;
        
        if( nremoved > 0 )
        {
            leftShiftColsValuesAndOffset(rind+1, nremoved);
        }
        
        #if SPM_DEBUG_MODE
            assert( getNumberOfElementsAtRow(rind) == rnel - nremoved );
        #endif
        
        totalRemoved += nremoved;
        
        if(i == 0)
            break;
    }
    
    if( totalRemoved > 0 )
    {
        #if SPM_DEBUG_MODE
            assert(totalRemoved <= nElements);
        #endif
        
        nElements -= totalRemoved;
        reallocateColsAndValues(nElements); //we do not care if realloction could not be done
    }
    
    updateNzRowIndex();
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::removeCols(const bool *cols)
{
    indexType* shiftCols = NULL; 
    
    
    shiftCols = (indexType *) calloc( ncols, sizeof(indexType) );
    if( !shiftCols )
    {
        #if SPM_MEMORY_ERROR
            SPM_PRINTMEMERROR;
        #endif
        return SPM_MEMORY_ERROR;
    }
    
    
    //first we remove the coefficients in the respective rows
    removeColsCoefficients(cols);
    
    //now, we fix the columns index
    for(unsigned int i = 0; i < ncols; i++)
    {
        if( cols[i] )
        {
            #pragma ivdep
            #pragma GCC ivdep
            for(unsigned int j = i+1; j < ncols; j++)
                shiftCols[j]++;
        }
    }
    
    indexType* const allCols = baseStructure.cols; //do not put it in the beginingbecause realloc in removeColsCoefficients can chang the pointer.
    
    #pragma ivdep
    #pragma GCC ivdep
    for(unsigned int i = 0; i < nElements; i++)
    {
        allCols[i] -= shiftCols[ allCols[i] ];
    }
    
    
    free(shiftCols);
    return 0;
}


template <class indexType, class matrixType>
int  SPM_NewSparseMatrix<indexType, matrixType>::removeCols(const unsigned int ncols, const indexType *cols)
{
    bool *fcols = NULL;
    
    
    fcols = (bool *) calloc(this->ncols, sizeof(bool));
    if( !fcols )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTMEMERROR;
        #endif
        return SPM_MEMORY_ERROR;
    }
    
    for(unsigned int i = 0; i < ncols; i++)
        fcols[ cols[i] ] = true;
    
    
    removeCols(fcols);
    
    
    free(fcols);
    
    #if SPM_DEBUG_MODE
        assert(checkConsistency() == 0);
    #endif
    
    return 0;
}


template <class indexType, class matrixType>
template <class intType>
int SPM_NewSparseMatrix<indexType, matrixType>::__removeRows( const unsigned int nrows, const intType *rows )
{
    if(nrows == 0)
        return 0;
    
    //there is no guarantee index are ordered. Anyway, we run from the end
    for( unsigned int i = nrows-1; ; i-- )
    {
        const unsigned int rind = rows[i];
        
        const auto rnel = getNumberOfElementsAtRow(rind);
        
        if( rnel > 0 )
        {
            leftShiftColsValuesAndOffset(rind+1, rnel);
            
            nElements -= rnel;
        }
        
        SPM_shiftLeftArray(this->nrows -rind-1, &offset[rind+1], 1); //offset has size this->nrows + 1
        
        this->nrows--;
        
        if(i == 0)
            break;
    }
    
    reallocateRowStructures(this->nrows);
    reallocateColsAndValues(nElements); //we do not care if fail
    updateNzRowIndex();
    
    #if SPM_DEBUG_MODE
        assert(checkConsistency() == 0);
    #endif
    
    return 0;
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::removeRows( const unsigned int nrows, const int *rows )
{
    return __removeRows(nrows, rows);
}

template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::removeRows( const unsigned int nrows, const unsigned int* rows )
{
    return __removeRows(nrows, rows);
}


//we remove the lines marked as true
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::removeRows( const bool *rowsFlag)
{
    unsigned int nrowsrem = 0;
    unsigned int *rowsInd = NULL;
    
    
    for(unsigned int i = 0; i < nrows; i++)
    {
        if( rowsFlag[i] )
            nrowsrem++;
    }
    
    rowsInd = (unsigned int *) malloc( nrowsrem * sizeof(unsigned int) );
    if( !rowsInd )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTMEMERROR;
        #endif
        return SPM_MEMORY_ERROR;
    }
    
    nrowsrem = 0;
    
    for(unsigned int i = 0; i < nrows; i++)
    {
        if( rowsFlag[i] )
        {
            rowsInd[nrowsrem] = i;
            nrowsrem++;
        }
    }
    
    removeRows(nrowsrem, rowsInd);
    
    return 0;
    
    /*SPM_NewSparseRow<indexType, matrixType>* paux;
    
    for(auto i = nrows-1;  ; i--) //nrows can be unsigned int
    {
        if( rowsFlag[i] )
        {
            const auto rnel = getNumberOfElementsAtRow(i); // rows[rind].nElements;
            
            if( rnel > 0 )
            {
                if( i < nrows-1 )
                    leftShiftColsValuesAndOffset(i+1, rnel);
                
                nElements -= rnel;
            }
            
            nrows--;
        }
        
        if(i == 0)
            break;
    }
    
    
    reallocateRowStructures(nrows);
    
    reallocateColsAndValues(nElements); //we do not care if fail
    updateNzRowIndex(); 
    
    #if SPM_DEBUG_MODE
        assert(checkConsistence() == 0);
    #endif */
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::leftShiftColsValuesAndOffset(const unsigned int startRow, const unsigned int shift)
{
    const unsigned int startIndex = offset[startRow]; // rows[startRow].__offset;
    
    #if SPM_DEBUG_MODE
        assert(startIndex >= shift); //we cannot shift from the first row and from a index lower than shift.
    #endif
    
    const unsigned int size = nElements-startIndex;
    
    if( size > 0 )
    {
        SPM_shiftLeftArray(size, &(baseStructure.cols[startIndex]), shift);
        SPM_shiftLeftArray(size, &(baseStructure.values[startIndex]), shift);
    }
    
    //for(decltype(nrows) i = startRow; i < nrows; i++)
        //offset[i] -= shift; //rows[i].__offset -= shift;
    
    shiftMinusRowOffset(startRow, shift);
}



template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::rightShiftColsValuesAndOffset(const unsigned int startRow, const unsigned int shift)
{
    const unsigned int startindex = offset[startRow]; //rows[startrow].__offset;
    
    const unsigned int size = nElements-startindex;
    
    if( size > 0 )
    {
        SPM_shiftRightArray(size, &(baseStructure.cols[startindex]), shift);
        SPM_shiftRightArray(size, &(baseStructure.values[startindex]), shift);
    }
    
    //for(decltype(nrows) i = startrow; i < nrows; i++)
        //offset[i] += shift; //rows[i].__offset += shift;
    
    shiftPlusRowOffset(startRow, shift);
}
    



//that method calculates M[row] * x
template <class indexType, class matrixType>
double SPM_NewSparseMatrix<indexType, matrixType>::rowEvaluation(const unsigned int rowIndex, const matrixType* x) const
{
    double v = 0.0;
    
    unsigned int rnels;
    indexType *rcols;
    matrixType *rvalues;
    
    getRowPointers(rowIndex, rnels, rcols, rvalues);
    
    for(unsigned int i = 0; i < rnels; i++)
        v += rvalues[i] *  x[ rcols[i] ];
    
    return v;
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::setAllElementsInARow(const unsigned int row, const matrixType value)
{
    SPM_setAllArray( getNumberOfElementsAtRow(row), getRowValuesPointer(row), value  );
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::setAllSparseMatrix(const matrixType value)
{
    SPM_setAllArray( nElements, baseStructure.values, value );
}


//This functions assumes that the position is in the matrix
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setElement(unsigned int rowIndex, indexType colIndex, const matrixType value, const bool printErrMsg)
{
    if( symmetric )
    {
        //we work only on the lower triangle
        if( (unsigned int) colIndex > rowIndex )
            SPM_swap(colIndex, rowIndex);
    }
    
    return setElementNoCheckSymmetry(rowIndex, colIndex, value, printErrMsg);
}


//this fucntion no check if col < row to symmetric matrices...
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setElementNoCheckSymmetry(unsigned int rowIndex, indexType colIndex, const matrixType value, const bool printErrMsg)
{
    unsigned int rnels;
    indexType *rcols;
    matrixType *rvalues;
    
    getRowPointers(rowIndex, rnels, rcols, rvalues);
    
    for( decltype(rnels) i = 0; i < rnels; i++ )
    {
        if( rcols[i] == colIndex )
        {
            rvalues[i] = value;
            return 0;
        }
    }
    
    if( printErrMsg )
        std::cerr << SPM_PREPRINT "Warning: element not found in sparse matrix. row:" << rowIndex << " col: " << colIndex << std::endl;
    
    //the element is not in the row
    return SPM_BAD_DEFINITIONS;
}


//set in the row, col makerd as true
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setRowStructure(const unsigned int rowIndex, const unsigned int sizeCols, const bool *cols)
{
    if(rowIndex >= nrows)
        return SPM_INDEX_FAULT;
    
    
    const unsigned int myncols = SPM_min(ncols, sizeCols);
    const unsigned int oldrnz = getNumberOfElementsAtRow(rowIndex);
    
    
    unsigned int newrnz = 0;
    
    
    for(unsigned int i = 0; i < myncols; i++)
    {
        if( cols[i] )
            newrnz++;
    }
    
    const int r = reallocateRowSpace(rowIndex, newrnz);
    if( r != 0 )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    
    if( newrnz > 0 )
    {
        indexType* const rcols = getRowColsPointer(rowIndex);
        
        newrnz = 0;
        for(unsigned int i = 0; i < myncols; i++)
        {
            if( cols[i] )
            {
                rcols[newrnz] = i;
                newrnz++;
            }
        }
        
        
        SPM_setAllArray<matrixType>(newrnz, getRowValuesPointer(rowIndex), 0);
    }
    
    
    #if SPM_DEBUG_MODE
        assert( getRowColsPointer(rowIndex+1) == &(getRowColsPointer(rowIndex)[newrnz]) );
    #endif
    
    if( bool(oldrnz) != bool(newrnz) ) //if( (oldrnz == 0 && newrnz > 0) || (oldrnz > 0 && newrnz == 0) )
        updateNzRowIndex();
    
    return 0;
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setRowStructure(const unsigned int rowIndex, const indexType numberOfNZElements, const indexType* cols, const matrixType *values)
{
    if(rowIndex >= nrows)
        return SPM_INDEX_FAULT;
    
    
    const unsigned int oldrnz = getNumberOfElementsAtRow(rowIndex);
    const int r = reallocateRowSpace(rowIndex, numberOfNZElements);
    if( r != 0 )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    
    SPM_copyArray( numberOfNZElements, cols, getRowColsPointer(rowIndex) );
    
    if( values )
        SPM_copyArray(numberOfNZElements, values, getRowValuesPointer(rowIndex));
    else
        SPM_setAllArray<matrixType>(numberOfNZElements, getRowValuesPointer(rowIndex), 0);
    
    
    if( bool(oldrnz) != bool(numberOfNZElements) ) //if( (oldrnz == 0 && numberOfNZElements > 0) || (oldrnz > 0 && numberOfNZElements == 0) )
        updateNzRowIndex();
    
    return 0;
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setRowStructure(const unsigned int rowIndex, const unsigned int sourceRowIndex,  const SPM_NewSparseMatrix<indexType, matrixType>& M, const bool copyValues)
{
    if(rowIndex >= nrows)
        return SPM_INDEX_FAULT;
    
    
    const unsigned int oldrnz = getNumberOfElementsAtRow(rowIndex);
    const unsigned int newrnz = M.getNumberOfElementsAtRow(sourceRowIndex);
    
    const int r = reallocateRowSpace(rowIndex, newrnz);
    if( r != 0 )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    SPM_copyArray( newrnz, M.getRowColsPointer(sourceRowIndex), getRowColsPointer(rowIndex) );
    
    if( copyValues )
        SPM_copyArray(newrnz, M.getRowValuesPointer(sourceRowIndex), getRowValuesPointer(rowIndex) );
    else
        SPM_setAllArray<matrixType>(newrnz, getRowValuesPointer(rowIndex), 0);
    
    
    
    if( bool(oldrnz) != bool(newrnz) ) //if( (oldrnz == 0 && newrnz > 0) || (oldrnz > 0 && newrnz == 0) )
        updateNzRowIndex();
    
    return 0;
}

#if 0
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setRowStructureAndValues(const unsigned int rowIndex, const unsigned int sourceRowIndex,  const SPM_NewSparseMatrix<indexType, matrixType>& M)
{
    if(rowIndex >= nrows)
        return SPM_INDEX_FAULT;
    
    
    const unsigned int newrnz = M.getNumberOfElementsAtRow(sourceRowIndex);
    
    const int r = setRowStructure(rowIndex, sourceRowIndex, M);
    if( r != 0 )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    SPM_copyArray(newrnz, M.getRowValuesPointer(sourceRowIndex), getRowValuesPointer(rowIndex) );
    
    //updateNzRowIndex is already called by setRowStructure
    
    return 0;
}
#endif


#if 0
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setRowStructureAndValues(const unsigned int rowIndex, const indexType numberOfNZElements, const indexType *cols, const matrixType *values)
{
    if(rowIndex >= nrows)
        return SPM_INDEX_FAULT;
    
    const int r = setRowStructure(rowIndex, numberOfNZElements, cols);
    if( r != 0 )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    SPM_copyArray(numberOfNZElements, values, getRowValuesPointer(rowIndex));
    
    //updateNzRowIndex is already called by setRowStructure
    
    return 0;
}
#endif


//if ncols = 0, spm will use the internal value of ncols...
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setRowStructureAndValues(const unsigned int rowIndex, const matrixType *a, unsigned int ncols, const double zeroTol)
{
    if(rowIndex >= nrows)
        return SPM_INDEX_FAULT;
    
    
    const unsigned int myncols = SPM_min(ncols, this->ncols);
    const unsigned int oldrnz = getNumberOfElementsAtRow(rowIndex);
    
    
    unsigned int newrnz = 0;
    
    for(unsigned int i = 0; i < myncols; i++)
    {
        if( SPM_abs(a[i]) > zeroTol )
            newrnz++;
    }
    
    const int r = reallocateRowSpace(rowIndex, newrnz);
    if( r != 0 )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    
    if( newrnz > 0 )
    {
        indexType* const rcols = getRowColsPointer(rowIndex);
        matrixType* const rvalues = getRowValuesPointer(rowIndex);
        
        newrnz = 0;
        for(unsigned int i = 0; i < myncols; i++)
        {
            if( SPM_abs(a[i]) > zeroTol )
            {
                rcols[newrnz] = i;
                rvalues[newrnz] = a[i];
                newrnz++;
            }
        }
    }
    
    #if SPM_DEBUG_MODE
        assert( getRowColsPointer(rowIndex+1) == &(getRowColsPointer(rowIndex)[newrnz]) );
    #endif
    
    
    if( bool(oldrnz) != bool(newrnz) ) //if( (oldrnz == 0 && newrnz > 0) || (oldrnz > 0 && newrnz == 0) )
        updateNzRowIndex();
    
    
    return 0;
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::setRowValues(const unsigned int rowIndex, const matrixType *values)
{
    SPM_copyArray<matrixType>( getNumberOfElementsAtRow(rowIndex), values, getRowValuesPointer(rowIndex) );
}


template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::setFullRowValues(const unsigned int rowIndex, const matrixType *values)
{
    const unsigned int rnz = getNumberOfElementsAtRow(rowIndex);
    const indexType *rcols = getRowColsPointer(rowIndex);
    matrixType *rvals = getRowValuesPointer(rowIndex);
    
    for(unsigned int k = 0; k < rnz; k++)
        rvals[k] = values[ rcols[k] ];
}


template <class indexType, class matrixType>
template <class intType>
int SPM_NewSparseMatrix<indexType, matrixType>::__setStructureAndValues(const unsigned int nzs, const intType* rows, const indexType* cols, const matrixType* vals, const bool reset)
{
    const bool hasvals = vals;
    int retCode;
    unsigned int *auxRows = NULL;
    unsigned int mynzs = 0;
    
    int r = reallocateColsAndValues(nzs);
    
    if( r != 0 )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    auxRows = (unsigned int *) calloc( nrows, sizeof(unsigned int) );
    if( !auxRows )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTMEMERROR;
        #endif
        return SPM_MEMORY_ERROR;
    }
    
    
    for(unsigned int i = 0; i < nzs; i++)
    {
        const unsigned int row = (symmetric && cols[i] > rows[i]) ? cols[i] : rows[i];
        
        if(row >= nrows)
        {
            #if SPM_DEBUG_MODE
                std::cerr << SPM_PREPRINT "Invalid index. rows["<<i<<"]: " << rows[i] << " cols["<<i<<"]: " << cols[i] << " nrows: " << nrows << "\n";
            #endif
            
            retCode = SPM_INDEX_FAULT;
            goto termination;
        }
        
        //std::cout << "iaia rows["<<i<<"]: " << rows[i] << " cols["<<i<<"]: " << cols[i] << "\n";
        
        auxRows[row]++;
    }
    
    for(unsigned int i = 0; i < nrows; i++)
    {
        offset[i] = mynzs;
        mynzs += auxRows[i];
    }
    
    
    offset[nrows] = nzs;
    
    #if SPM_DEBUG_MODE
        assert( mynzs == nzs );
    #endif
    
    SPM_setAllArray<unsigned int>(nrows, auxRows, 0);
    
    for(unsigned int i = 0; i < nzs; i++)
    {
        unsigned int row = rows[i];
        unsigned int col = cols[i];
        
        if( symmetric && col > row )
            SPM_swap(row, col);
        
        indexType* const rcols = getRowColsPointer(row);
        rcols[ auxRows[row] ] = col;
        
        if(hasvals)
        {
            matrixType* const rvalues = getRowValuesPointer(row);
            rvalues[ auxRows[row] ] = vals[i];
        }
        
        auxRows[row]++;
    }
    
    #if SPM_DEBUG_MODE
        for(unsigned int i = 0; i < nrows; i++)
        {
            //std::cout << "auxRows["<<i<<"]: " << auxRows[i] << " offset["<<i+1<<"]: " << offset[i+1] << "\n";
            assert( auxRows[i] == offset[i+1] -offset[i] );
        }
    #endif
    
    //SPM_getchar();
    
    nElements = nzs;
    
    updateNzRowIndex();
    
    
    retCode = 0;
    
termination:
    
    if(auxRows) free(auxRows);
    
    return retCode;
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setStructureAndValues(const unsigned int nzs, const int* iRow, const indexType* jCol, const matrixType* vals, const bool reset)
{
    return __setStructureAndValues(nzs, iRow, jCol, vals, reset);
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>:: setStructureAndValues(const unsigned int nzs, const unsigned int* rows, const indexType* cols, const matrixType* vals, const bool reset)
{
    return __setStructureAndValues(nzs, rows, cols, vals, reset);
}


//set matrix by compressed row format. Note, this method does not check abouy symmetry...
template <class indexType, class matrixType>
template <class intType>
int SPM_NewSparseMatrix<indexType, matrixType>::__setStructureAndValues( const intType *rowStart, const indexType *cols, const matrixType* vals)
{
    const unsigned int nzs = rowStart[nrows];
    
    if( rowStart[0] != 0 )
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORMSG("rowStart[0] is not zero!");
        #endif
        return SPM_BAD_DEFINITIONS;
    }
    
    
    /*if( nrows != this->nrows )
    {
        const int r = reallocateRowStructures(nrows);
        if(r != 0)
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
        
        this->nrows = nrows;
    }*/
    
    
    int r = reallocateColsAndValues(nzs);
    if(r != 0)
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    nElements = nzs;
    
    SPM_copyArray(nrows+1, rowStart, offset);
    
    if( cols )
        SPM_copyArray(nzs, cols, baseStructure.cols);
    
    if( vals )
        SPM_copyArray(nzs, vals, baseStructure.values);
    else
        SPM_setAllArray<matrixType>(nzs, baseStructure.values, 0);
    
    updateNzRowIndex();
    
    return 0;
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setStructureAndValues(const unsigned int *rowStart, const indexType *cols, const matrixType* vals)
{
    return __setStructureAndValues(rowStart, cols, vals);
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setStructureAndValues(const int *rowStart, const indexType *cols, const matrixType* vals)
{
    return __setStructureAndValues(rowStart, cols, vals);
}


template <class indexType, class matrixType>
template <class matrixTypePointer>
int SPM_NewSparseMatrix<indexType, matrixType>::__setStructureAndValues(matrixTypePointer A, const double zero_tol, unsigned int nrows, unsigned int ncols, const bool lowerTriangle)
{
    unsigned int nElements = 0;
    
    if(nrows == 0)
        nrows = this->nrows;
    
    if(ncols == 0)
        ncols = this->ncols;
    
    if( nrows > this->nrows )
    {
        const int r = allocateSparseRows( nrows );
        if( r != 0 )
        {
            #if SPM_DEBUG_MODE
                SPM_PRINTERRORNUMBER(r);
            #endif
            return r;
        }
    }
    
    if( ncols > this->ncols )
        this->ncols = ncols;
    
    
    //counting the nonzeros in the matriz
    for(unsigned int i = 0; i < nrows; i++)
    {
        const unsigned int maxcols = symmetric ? i+1: ncols;;
        const matrixType *Ai = getRowMatrixPointer(A, ncols, i, lowerTriangle); //&A[i*ncols];
        
        for(unsigned int j = 0; j < maxcols; j++)
        {
            if( SPM_abs(Ai[j]) > zero_tol )
                nElements++;
        }
    }
    
    const int r = reallocateColsAndValues(nElements);
    if(r != 0)
    {
        #if SPM_DEBUG_MODE
            SPM_PRINTERRORNUMBER(r);
        #endif
        return r;
    }
    
    this->nElements = nElements;
    
    nElements = 0;
    
    for(unsigned int i = 0; i < nrows; i++)
    {
        const matrixType *Ai = getRowMatrixPointer(A, ncols, i, lowerTriangle);
        const unsigned int maxcol = symmetric ? i+1 : ncols;
        
        unsigned int rnz = 0;
        offset[i] = nElements;
        
        indexType* const rcols = getRowColsPointer(i);
        matrixType* const rvalues = getRowValuesPointer(i);
        
        for(unsigned int j = 0; j < maxcol; j++)
        {
            if( SPM_abs(Ai[j]) > zero_tol )
            {
                rcols[rnz] = j;
                rvalues[rnz] = Ai[j];
                rnz++;
            }
        }
        
        nElements +=rnz;
    }
    
    offset[nrows] = nElements;
    
    #if SPM_DEBUG_MODE
        assert( nElements == this->nElements );
    #endif
    
    updateNzRowIndex();
    
    return 0;
}



//the matrix is stored in a vector ordered by row. if nrows or ncols is zero, spm will use the values stored inside object
template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setStructureAndValues(const matrixType* A, const double zeroTolerance, unsigned int nrows, unsigned int ncols)
{
    return __setStructureAndValues(A, zeroTolerance, nrows, ncols, false);
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setStructureAndValues(matrixType **A, const double zeroTolerance, unsigned int nrows, unsigned int ncols)
{
    return __setStructureAndValues(A, zeroTolerance, nrows, ncols, false);
}


template <class indexType, class matrixType>
int SPM_NewSparseMatrix<indexType, matrixType>::setSymmetricStructureAndValuesByLowerTriangle(const matrixType* lowerTA, const double zeroTolerance)
{
    return __setStructureAndValues(lowerTA, zeroTolerance, nrows, ncols, true);
}


//that method sum all lines times a respective factor in vector. That vector is initialized with zeros. factors can be a null pointer. In this case, we consider all factors like 1.0
template <class indexType, class matrixType>
void SPM_NewSparseMatrix<indexType, matrixType>::sumAllLines(const double *factors, matrixType *v) const
{
    SPM_setAllArray<matrixType>(ncols, v, 0);
    
    for(unsigned int i = 0; i < nNzRowIndices; i++)
    {
        const unsigned int rind = nzRowIndex[i];
        
        const unsigned int rnzs = getNumberOfElementsAtRow(rind);
        indexType* const rcols = getRowColsPointer(rind);
        matrixType* const rvalues = getRowValuesPointer(rind);
        
        if( factors )
        {
            const double f = factors[rind];
            
            for(unsigned int j = 0; j < rnzs; j++)
                v[ rcols[j] ] += f*rvalues[j];
        }
        else
        {
            for(unsigned int j = 0; j < rnzs; j++)
                v[ rcols[j] ] += rvalues[j];
        }
    }
}



/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/




template <class intType>
inline int SPM_checkRowStructure(const unsigned int nrows, const unsigned int ncols, const bool symmetric, unsigned int row, unsigned int nzs, intType *cols, double *values = NULL)
{
    if( row < 0 || row >= nrows )
        return SPM_INDEX_FAULT;
    
    
    for(unsigned int i = 0; i < nzs; i++)
    {
        if( cols[i] < 0 || (unsigned int) cols[i] >= nrows )
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
            if( std::isinf(values[i]) || std::isnan(values[i]) )
                return SPM_BAD_VALUE;
        }
    }
    
    if(symmetric)
    {
        for(unsigned int i = 0; i < nzs; i++)
        {
            if( (unsigned int) cols[i] > row )
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
        if( rows[i] < 0 || (unsigned int) rows[i] >= nrows )
            return SPM_INDEX_FAULT;
    }
    
    for(unsigned int i = 0; i < nzs; i++)
    {
        if( cols[i] < 0 || (unsigned int) cols[i] >= ncols )
            return SPM_INDEX_FAULT;
    }
    
    
    
    if( values )
    {
        for(unsigned int i = 0; i < nzs; i++)
        {
            if( std::isinf(values[i]) || std::isnan(values[i]) )
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








/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/




















template <class indexType, class matrixType>
newspm:: SPM_NewSparseRow<indexType, matrixType>:: SPM_NewSparseRow(SPM_ColsValues<indexType, matrixType> *baseStructure, const unsigned int offset)
{
    initialize(baseStructure, offset);
}


template <class indexType, class matrixType>
newspm:: SPM_NewSparseRow<indexType, matrixType>::~SPM_NewSparseRow()
{
    desallocate();
}


//this method sum the coefficient in this line times a factor in a array. The vector is not iniatialized (note there is no symmetry concecpt for a line operation).
template <class indexType, class matrixType>
void SPM_NewSparseRow<indexType, matrixType>::accumulateInArray(matrixType *v, const matrixType factor) const
{
    matrixType* const values = getValuesPointer();
    indexType* const cols = getColsPointer();
    
    for(indexType i = 0; i < nElements; i++)
        v[ cols[i] ] += factor*values[i];
}


template <class indexType, class matrixType>
int SPM_NewSparseRow<indexType, matrixType>::addToElement(const indexType col, const matrixType value)
{
    matrixType* const values = getValuesPointer();
    indexType* const cols = getColsPointer();
    
    for(indexType i = 0; i < nElements; i++)
    {
        if( cols[i] == col )
        {
            values[i] += value;
            return 0;
        }
    }
    
    return SPM_ELEMENT_NOT_PRESENT;
}


template <class indexType, class matrixType>
void SPM_NewSparseRow<indexType, matrixType>::addToAllElements(const matrixType value)
{
    SPM_addToAllArray(nElements, getValuesPointer(), value);
}


template <class indexType, class matrixType>
void SPM_NewSparseRow<indexType, matrixType>::copyTo(matrixType* row, const matrixType factor)
{
    matrixType* const values = getValuesPointer();
    indexType* const cols = getColsPointer();
    
    
    for(unsigned int i = 0; i < nElements; i++)
        row[ cols[i] ] = values[i];
    
    if( factor != 1 )
    {
        for( unsigned int i = 0; i < nElements; i++ )
            row[ cols[i] ] *= factor;
    }
}


template <class indexType, class matrixType>
void SPM_NewSparseRow<indexType, matrixType>::desallocate()
{
    //initialize();
    nElements = 0;
}


template <class indexType, class matrixType>
double SPM_NewSparseRow<indexType, matrixType>::evalTimesxt(const matrixType *x) const
{
    matrixType* const values = getValuesPointer();
    indexType* const cols = getColsPointer();
    
    matrixType v = 0;
    
    for(indexType i = 0; i < nElements; i++)
        v += values[i] * x[cols[i]];
    
    return v;
}


template <class indexType, class matrixType>
void SPM_NewSparseRow<indexType, matrixType>::initialize(SPM_ColsValues<indexType, matrixType> *baseStructure, const unsigned int offset )
{
    this->baseStructure = baseStructure;
    this->__offset = offset;
    nElements = 0;
}


template <class indexType, class matrixType>
int SPM_NewSparseRow<indexType, matrixType>::getElement(const indexType col, matrixType &value) const
{
    matrixType* const values = getValuesPointer();
    indexType* const cols = getColsPointer();
    
    for(indexType i = 0; i < nElements; i++)
    {
        if( cols[i] == col )
        {
            value = values[i];
            return 0;
        }
    }
    
    value = 0.0;
    return SPM_ELEMENT_NOT_PRESENT;
}


//that method return the number of elements in the row. You must initialize by yourself all positions in cols as false
template <class indexType, class matrixType>
indexType SPM_NewSparseRow<indexType, matrixType>::getStructure(bool *columns) const
{
    indexType* const cols = getColsPointer();
    
    #pragma ivdep
    #pragma GCC ivdep
    for(indexType i = 0; i < nElements; i++)
        columns[ cols[i] ] = true;
    
    return nElements;
}


//those methods return the number of elements in the row
template <class indexType, class matrixType>
indexType SPM_NewSparseRow<indexType, matrixType>::getStructure(int *columns) const
{
    SPM_copyArray(nElements, getColsPointer(), columns);
    return nElements;
}


template <class indexType, class matrixType>
indexType SPM_NewSparseRow<indexType, matrixType>::getStructure(unsigned int *columns) const
{
    SPM_copyArray(nElements, getColsPointer(), columns);
    return nElements;
}


template <class indexType, class matrixType>
unsigned int SPM_NewSparseRow<indexType, matrixType>::getValues(matrixType *values) const
{
    SPM_copyArray(nElements, getValuesPointer(), values);
    return nElements;
}


template <class indexType, class matrixType>
bool SPM_NewSparseRow<indexType, matrixType>::hasColumn(const indexType col, indexType* index, matrixType *value) const
{
    indexType* const cols = getColsPointer();
    
    for(indexType i = 0; i < nElements; i++)
    {
        if( cols[i] == col )
        {
            if(index)
                *index = i;
            if(value)
                *value = getValuesPointer()[i];
            return true;
        }
    }
    
    return false;
}


template <class indexType, class matrixType>
void SPM_NewSparseRow<indexType, matrixType>::multiplyAllElements(const matrixType value)
{
    SPM_multiplyAllArray(nElements, getValuesPointer(), value);
}






template <class indexType, class matrixType>
void SPM_NewSparseRow<indexType, matrixType>::print( std::ostream &out )
{
    matrixType* const values = getValuesPointer();
    indexType* const cols = getColsPointer();
    
    for(indexType i = 0; i < nElements; i++)
        out << cols[i] << ": " << values[i] << " ";
}


template <class indexType, class matrixType>
void SPM_NewSparseRow<indexType, matrixType>::setAllElements(const matrixType value)
{
    SPM_setAllArray(nElements, getValuesPointer(), value);
}


template <class indexType, class matrixType>
int SPM_NewSparseRow<indexType, matrixType>::setElement(const indexType col, const matrixType value)
{
    matrixType* const values = getValuesPointer();
    indexType* const cols = getColsPointer();
    
    for(indexType i = 0; i < nElements; i++)
    {
        if( cols[i] == col )
        {
            values[i] = value;
            return 0;
        }
    }
    
    #if SPM_DEBUG_MODE
        std::cerr << SPM_PREPRINT << "Element " << col << " not found in sparse row" << SPM_GETFILELINE;
    #endif
    
    return SPM_ELEMENT_NOT_PRESENT;
}


template <class indexType, class matrixType>
//set first 'nel' elemets in this line using the first 'nel' values in array values. Note, this method does not look to columns indexes...
int SPM_NewSparseRow<indexType, matrixType>::setElementsByOrder(indexType nel, const matrixType *values)
{
    SPM_copyArray(nel, values, getValuesPointer());
    return 0;
}


template <class indexType, class matrixType>
//return the element of row in position pos
matrixType& SPM_NewSparseRow<indexType, matrixType>::operator()(const unsigned int pos)
{
    return getValuesPointer()[pos];
}


template <class indexType, class matrixType>
//return the column of row in position pos
indexType& SPM_NewSparseRow<indexType, matrixType>::operator[](const unsigned int pos)
{
    return getColsPointer()[pos];
}




}










#endif
