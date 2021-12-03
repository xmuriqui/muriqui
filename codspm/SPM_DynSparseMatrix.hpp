/***********************************************************
***********************************************************
****                                                   ****
****      implementation of dynamic sparse matris      ****
**** the goal here, is implement a class to manage a	**** 
**** sparse matrix that can change frequenly using		****
**** C++ containers.									****
****                                                   ****
***********************************************************
***********************************************************/

#ifndef SPM_DYN_SPARSE_MATRIX_HPP
#define SPM_DYN_SPARSE_MATRIX_HPP

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include <set>
#include <map>

#include "../WAXM_config.h"



#include "SPM_tools.hpp"









namespace newspm{


typedef unsigned int SPM_dynIndex;
typedef double       SPM_dynValue;


//we are working with ordered map because they are implemented using binary trees while unordered are hash tables and get more memory space.
typedef std::map<SPM_dynIndex, SPM_dynValue> SPM_rowMap;


class SPM_DynSparseRow
{
public:
    
    
    SPM_rowMap element;
    
    
    virtual ~SPM_DynSparseRow()
    {
    }
    
    
    bool hasElement(const SPM_dynIndex &col)
    {
        return element.count(col) > 0;
    }
    
    
    int eraseCoefficient(const SPM_dynIndex &col)
    {
        if( hasElement(col) )
        {
            element.erase(col);
            return 0;
        }
        else
        {
            return SPM_INDEX_FAULT;
        }
    }
    
    
    //this method makes a shift in greater colum indices than col
    int removeColumns(const SPM_dynIndex ncols, const SPM_dynIndex *cols)
    {
        SPM_rowMap ecopy(element); //make a copy of our map
        std::set<SPM_dynIndex> setCols; //to avoid problems if some columns appears twice in cols
        
        //first, we remove the columns
        for(SPM_dynIndex i = 0; i < ncols; i++)
        {
            if( ecopy.count(cols[i]) > 0 )
                ecopy.erase( cols[i] );
            
            setCols.insert(cols[i]);
        }
        
        element.clear();
        
        for( auto it = ecopy.begin(); it != ecopy.end(); ++it )
        {
            SPM_dynIndex mycol = it->first;
            
            for(auto it2 = setCols.begin(); it2 != setCols.end(); ++it2)
            {
                if( mycol > *it2 )
                    mycol--;
            }
            
            #if SPM_DEBUG_MODE
                assert( element.count(mycol) == 0  );
            #endif
            
            element[mycol] = it->second;
        }
        
        return 0;
    }
    
    
    int setCoefficient(const SPM_dynIndex &col, const SPM_dynValue &value )
    {
        element[col] = value;
        return 0;
    }
    
    
    int getCoefficient(const SPM_dynIndex &col, SPM_dynValue &value)
    {
        if( hasElement(col) )
        {
            value = element[col];
            return 0;
        }
        else
        {
            value = 0;
            return SPM_INDEX_FAULT;
        }
    }
    
    
    //erase all elements
    int clear( )
    {
        element.clear();
        return 0;
    }
    
    
    int setCoefficients(const SPM_dynIndex &nElements, const SPM_dynIndex *cols, const SPM_dynValue *values = NULL)
    {
        
        for(SPM_dynIndex i = 0; i < nElements; i++)
        {
            int r = setCoefficient( cols[i],  values ? values[i] : 0.0 );
            SPM_IFERRORRETURN(r, r);
        }
        
        return 0;
    }
    
    
    //returns the number of elements
    int getCoefficients(SPM_dynIndex *cols, SPM_dynValue *values = NULL)
    {
        SPM_dynIndex k = 0 ;
        
        for( auto it = element.begin(); it != element.end(); ++it )
        {
            if(cols)
                cols[k] = it->first;
            
            if(values)
                values[k] = it->second;
            k++;
        }
        
        return k;
    }
    
    //return the full row (not sparse)
    int getFullRow(const SPM_dynIndex &ncols, SPM_dynValue* values, const bool accumulate, const double factor = 1.0)
    {
        if( !accumulate )
            SPM_setAllArray<SPM_dynValue>(ncols, values, 0.0);
        
        for( auto it = element.begin(); it != element.end(); ++it )
            values[it->first] += it->second;
        
        if(factor != 1.0)
            SPM_multiplyAllArray<SPM_dynValue>(ncols, values, factor);
        
        return 0;
    }
    
    
    SPM_dynIndex getNumberOfElements(  ) const
    {
        return element.size();
    }
    
    
};



typedef std::map<SPM_dynIndex, SPM_DynSparseRow*> SPM_matrixMap;

class SPM_DynSparseMatrix
{
public:
    
    SPM_dynIndex nrows_, ncols_;
    SPM_matrixMap rows_;
    
    
    SPM_DynSparseMatrix(const SPM_dynIndex &nrows = 0, const SPM_dynIndex &ncols = 0)
    {
        nrows_ = nrows;
        ncols_ = ncols;
    }
    
    
    virtual ~SPM_DynSparseMatrix()
    {	
        desallocateMemory();
    }
    
    
    void desallocateMemory()
    {
        for( auto it = rows_.begin(); it != rows_.end(); ++it )
        {
            delete it->second;
        }
        
        rows_.clear();
        nrows_ = 0;
    }
    
    
    int addRows(const SPM_dynIndex &numberOfNewRows)
    {
        nrows_ += numberOfNewRows;
        return 0;
    }
    
    
    int addCols(const SPM_dynIndex &numberOfNewCols)
    {
        ncols_ += numberOfNewCols;
        return 0;
    }
    
    
    SPM_dynIndex getNumberOfRows()
    {
        return nrows_;
    }
    
    
    int removeRows(const SPM_dynIndex &numberOfRows, const SPM_dynIndex *rows)
    {
        std::set<SPM_dynIndex> remove;
        SPM_dynIndex shift = 0;
        
        int r = checkRowIndices(numberOfRows, rows);
        SPM_IFERRORRETURN(r, r);
        
        
        for( SPM_dynIndex i = 0; i < numberOfRows; i++ )
            remove.insert( rows[i] );
        
        
        for( SPM_dynIndex i = 0; i < nrows_; i++ )
        {
            if( remove.count(i) > 0 && hasRow(i) )
            {
                delete rows_[i];
                rows_[i] = NULL;
                shift++;
            }
            else if(shift > 0)
            { //so, we have to shift the row above
                
                #if SPM_DEBUG_MODE
                    assert(rows_[i - shift] == NULL);
                #endif
                
                rows_[i - shift] = rows_[i];
                rows_[i] = NULL;
            }
            
        }
        
        
        for( SPM_dynIndex i = nrows_ - numberOfRows; i < nrows_ ; i++ )
        {
            if( hasRow(i) )
            {
                #if SPM_DEBUG_MODE
                    assert( rows_[i] == NULL );
                #endif
                rows_.erase(i);
            }
        }
        
        nrows_ -= remove.size(); //remove is a set, i.e., no repetition of elements
        
        
        #if SPM_DEBUG_MODE
            assert( rows_.size() <= nrows_ );
        #endif
        
        return 0;
    }
    
    
    int removeColumns(const SPM_dynIndex &numberOfCols, const SPM_dynIndex *cols, const bool checkRepetitions = false)
    {
        int r;
        
        r = checkColIndices(numberOfCols, cols);
        
        for( auto it = rows_.begin(); it != rows_.end(); ++it )
        {
            r = (it->second)->removeColumns(numberOfCols, cols); //SPM_DynSparseRow::removeColumns already check repetitions
            SPM_IFERRORRETURN(r, r);
        }
        
        if( checkRepetitions )
        {
            std::set<SPM_dynIndex> setCols;
            
            for(SPM_dynIndex k = 0; k < numberOfCols; k++)
                setCols.insert( cols[k] );
            
            ncols_ -= setCols.size();
        }
        else
        {
            ncols_ -= numberOfCols;
        }
        
        return 0;
    }
    
    
    bool hasRow(const SPM_dynIndex &row)
    {
        return rows_.count(row) > 0;
    }
    
    
    int generateRow(const SPM_dynIndex &row)
    {
        int r = checkRowIndices(1, &row);
        SPM_IFERRORRETURN(r, r);
        
        if( !hasRow(row) )
        {
            SPM_DynSparseRow *r = new (std::nothrow) SPM_DynSparseRow;
            SPM_IFMEMERRORRETURN(!r);
            
            rows_[row] = r;
        }
        
        return 0;
    }
    
    
    int getRowPointer(const SPM_dynIndex &rowIndex, SPM_DynSparseRow* &row  )
    {
        int r = checkRowIndices(1, &rowIndex);
        SPM_IFERRORRETURN(r, r);
        
        row = hasRow(rowIndex) ? rows_[rowIndex] : NULL;
        return 0;
    }
    
    
    int getCoefficient(const SPM_dynIndex &row, const SPM_dynIndex &col, SPM_dynValue &value)
    {
        int r;
        
        r = checkIndices(1, &row, &col);
        SPM_IFERRORRETURN(r, r);
        
        if( !hasRow(row) )
        {
            value = 0;
            return SPM_INDEX_FAULT;
        }
        
        r = rows_[row]->getCoefficient(col, value);
        SPM_IFERRORRETURN(r, r);
        
        return 0;
    }
    
    
    int setCoefficient(const SPM_dynIndex &row, const SPM_dynIndex &col, const SPM_dynValue &value)
    {
        int r;
        
        r = checkIndices(1, &row, &col);
        SPM_IFERRORRETURN(r, r);
        
        r = generateRow(row);  //if row already exists, does nothing
        SPM_IFERRORRETURN(r, r);
        
        r = rows_[row]->setCoefficient(col, value);
        SPM_IFERRORRETURN(r, r);
        
        return 0;
    }
    
    
    int getCoefficientsInARow(const SPM_dynIndex &row, SPM_dynIndex &nElements, SPM_dynIndex *cols, SPM_dynValue *values)
    {
        int r = checkRowIndices(1, &row);
        SPM_IFERRORRETURN(r, r);
        
        
        if( hasRow(row) )
            nElements = rows_[row]->getCoefficients(cols, values);
        else
            nElements = 0;
        
        return 0;
    }
    
    
    int setCoefficientsInARow(const SPM_dynIndex &row, const SPM_dynIndex &nElements, const SPM_dynIndex *cols, const SPM_dynValue *values = NULL)
    {
        int r;
        
        r = checkRowIndices(1, &row);
        SPM_IFERRORRETURN(r, r);
        
        r = generateRow(row);  //if row already exists, does nothing
        SPM_IFERRORRETURN(r, r);
        
        
        r = rows_[row]->setCoefficients(nElements, cols, values);
        SPM_IFERRORRETURN(r, r);
        
        return 0;
    }
    
    
    int getFullRow(const SPM_dynIndex &row, SPM_dynValue* values, const bool accumulate, const double factor = 1.0)
    {
        int r = checkRowIndices(1, &row);
        SPM_IFERRORRETURN(r, r);
        
        if( hasRow(row) )
        {
            r = rows_[row]->getFullRow(ncols_, values, accumulate, factor);
            SPM_IFERRORRETURN(r, r);
        }
        else
        {
            if(!accumulate)
                SPM_setAllArray<SPM_dynValue>(ncols_, values, 0.0);
        }
        
        return 0;
    }
    
    
    int checkRowIndices(const SPM_dynIndex &nElements, const SPM_dynIndex *rows)
    {
        for(SPM_dynIndex i = 0; i < nElements; i++)
        {
            if( rows[i] < 0 || rows[i] >= nrows_ )
            {
                SPM_PRINTINDEXERRORP(rows[i]);
                return SPM_INDEX_FAULT;
            }
        }
        
        return 0;
    }
    
    
    int checkColIndices(const SPM_dynIndex &nElements, const SPM_dynIndex *cols)
    {
        for(SPM_dynIndex i = 0; i < nElements; i++)
        {
            if( cols[i] < 0 || cols[i] >= ncols_ )
            {
                SPM_PRINTINDEXERRORP(cols[i]);
                return SPM_INDEX_FAULT;
            }
        }
        
        return 0;
    }
    
    
    int checkIndices(const SPM_dynIndex &nElements, const SPM_dynIndex *rows, const SPM_dynIndex *cols)
    {
        int r;
        
        r = checkRowIndices(nElements, rows);
        SPM_IFERRORRETURN(r, r);
        
        r = checkColIndices(nElements, cols);
        SPM_IFERRORRETURN(r, r);
        
        return 0;
    }
    
    
    //this method returns the number os elements
    int getCoefficients(SPM_dynIndex *rows, SPM_dynIndex *cols, SPM_dynValue *values)
    {
        int r;
        SPM_dynIndex nEl = 0;
        
        for( auto it = rows_.begin(); it != rows_.end(); ++it )
        {
            const SPM_dynIndex rowIndex = it->first;
            SPM_dynIndex rowNEl, rowNEl2;
            
            getNumberOfElementsInARow(rowIndex, rowNEl);
            
            
            SPM_setAllArray( rowNEl, &rows[nEl], rowIndex );
            
            r = getCoefficientsInARow(rowIndex, rowNEl2, &cols[nEl], &values[nEl]);
            SPM_IFERRORRETURN(r, r);
            
            #if SPM_DEBUG_MODE
                assert( rowNEl == rowNEl2 );
            #endif
            
            nEl += rowNEl;
        }
        
        return nEl;
    }
    
    
    int setCoefficients(const SPM_dynIndex &nElements, const SPM_dynIndex *rows, const SPM_dynIndex *cols, const SPM_dynValue *values = NULL)
    {
        int r;
        
        r = checkIndices(nElements, rows, cols);
        SPM_IFERRORRETURN(r, r);
        
        for(SPM_dynIndex i = 0; i < nElements; i++)
        {
            const SPM_dynValue value = values ? values[i] : 0.0;
            
            r = setCoefficient( rows[i], cols[i], value);
            SPM_IFERRORRETURN(r, r);
        }
        
        return 0;
    }
    
    
    SPM_dynIndex getNumberOfElements()
    {
        SPM_dynIndex nEl = 0;
        
        for( auto it = rows_.begin(); it != rows_.end(); ++it )
            nEl += (it->second)->getNumberOfElements();
        
        return nEl;
    }
    
    
    int getNumberOfElementsInARow(const SPM_dynIndex &row, SPM_dynIndex &nElements)
    {
        int r = checkRowIndices(1, &row);
        SPM_IFERRORRETURN(r, r);
        
        
        nElements = hasRow(row) ? rows_[row]->getNumberOfElements() : 0;
        
        return 0;
    }
    
    
};


};



#endif
