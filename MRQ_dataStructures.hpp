/*
* MRQ_dataStructure.hpp
*
*  Created on: 27/08/2013
*      Author: yo
*/

#ifndef MRQ_DATASTRUCTURE_HPP_
#define MRQ_DATASTRUCTURE_HPP_



//#include "SPM_SparseMatrix.hpp"
#include "MIP_minlpProblem.hpp"
#include "MRQ_constants.hpp"



namespace muriqui
{
    
    
    typedef minlpproblem::MIP_SparseMatrix MRQ_SparseMatrix;
    
    
    typedef minlpproblem::MIP_VARTYPE MRQ_VARTYPE;
    
    
    #define MRQ_VT_CONTINUOUS minlpproblem::MIP_VT_CONTINUOUS
    #define MRQ_VT_INTEGER minlpproblem::MIP_VT_INTEGER
    
    
    
    /*Class for store the MINLP problem data
    *
    * Min  c'x + 0.5x'Qx + f(x) + d
    * subject to:
    *
    * 	l_c <= a_ix + 0.5x'Q_ix + g_i(x) <= u_c
    *
    *  l_x <= x <= u_x
    *
    *  x_i is integer for i \in I
    *
    */
    
    typedef minlpproblem::MIP_NonLinearEval MRQ_NonLinearEval;
    typedef minlpproblem::MIP_MINLPProb MRQ_MINLPProb;
    
    
    typedef optsolvers::OPT_GeneralSolverParams MRQ_GeneralSolverParams;
    
    
    
    
    
    
    class MRQ_HistorySolution
    {
        
        int allocateSol(const int n);
        
    public:
        int nvars; //total number of variables
        long unsigned int iter; //iteration when that solution was gotten
        double time; //time (wall) when that solution was gotten
        double cputime; //cpu time when that solution was gotten
        double objF;
        double *sol; //solution
        
        MRQ_HistorySolution();
        
        //MRQ_HistorySolution(const int n, const long unsigned int iter, const double time, const double cputime, const double* sol, const double objF);
        
        ~MRQ_HistorySolution();
        
        void freeSolution();
        
        int getnvars();
        
        int getiter();
        
        double getobjvalue();
        
        double gettime();
        
        double getcputime();
        
        int getsolution(double *solution);
        
        int set(const int n, const long unsigned int iter, const double time, const double cputime, const double* sol, const double objF);
    };
    
    
    
    class MRQ_SolutionHistory
    {
        unsigned int nsols;
        MRQ_HistorySolution **hsols; //we use a double pointer to turn easier and faster the reallocation process...
        
    public:
        
        MRQ_SolutionHistory();
        
        ~MRQ_SolutionHistory();
        
        void desallocate();
        
        unsigned int getnsols() const;
        
        int addSolution(const int n, const long unsigned int iter, const double time, const double cputime, const double *sol, const double objF);
        
        //it is only a pointer, not a copy. Do not free. If you want a copy of solution use the method getsol of MRQ_HistorySolution pointer
        MRQ_HistorySolution * getHistSolPointer(const unsigned int index);
        
    };
    
    
    
    
}








#endif /* MRQ_DATASTRUCTURE_HPP_ */
