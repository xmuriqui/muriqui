

#ifndef MIP_GAMS_HPP
#define MIP_GAMS_HPP

#include <string>

#include "MIP_config.hpp"
#include "MIP_minlpProblem.hpp"


#if MIP_HAVE_GAMS
extern "C"{
    #include "gmomcc.h"
    #include "gevmcc.h"
    //#include "gclgms.h"
}
#endif



namespace minlpproblem
{
    
    class MIP_NLEvalGams : public MIP_NonLinearEval
    {
        
        unsigned int nEvalErros;
        int nthreads;
        
    public:
        
        const char *stub;
        double **auxValues_;
        MIP_MINLPProb *prob;
        
        #if MIP_HAVE_GAMS
            struct gmoRec        **gmos_;                // GAMS modeling object
            struct gevRec        **gevs_;                // GAMS environment
        #endif
        
        
        
        MIP_NLEvalGams(const char *stub, MIP_MINLPProb *prob);
        
        
        virtual ~MIP_NLEvalGams();
        
        
        void desallocate();
        
        
        virtual int initialize(const int nthreads, const int n, const int m, const int nzJac, const int nzLagHess);
        
        
        virtual int eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value);
        
        
        virtual int eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *values);
        
        
        virtual int eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *ctrEval, const double *x, double *values);
        
        
        virtual int eval_grad_nl_constrs_part( const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *ctrEval, const double *x, MIP_SparseMatrix& J);
        
        
        virtual int eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double obj_factor, const double *lambda, MIP_SparseMatrix& hessian);
        
        
        virtual void finalize(const int nthreads, const int n, const int mnl, const int nzNLJac, const int nzNLLagHess);
        
    };
    
    
    #define MIP_GAMS_BUFFER_SIZE 1024
    
    
    class MIP_GamsModelReader
    {
    public:
        
        #if MIP_HAVE_GAMS
        struct gmoRec*        gmo;                // GAMS modeling object
        struct gevRec*        gev;                // GAMS environment
        //struct optRec*        opt; // GAMS options object
        #endif
        
        
        MIP_GamsModelReader();
        
        virtual ~MIP_GamsModelReader();
        
        void desallocate();
        
        int getInitialSolution(double *solution);
        
        int getDblOption(const char *optionName, double &value);
        
        //to gen an integer option (parameter)
        int getIntOption(const char *optionName, int &value);
        
        int getStrOption(const char *optionName, std::string &value);
        
        bool isMaximizationProblem();
        
        //if there is no option file, return a error code.
        int getOptionFileName(std::string &fileName);
        
        int readProblem(const char* stub, MIP_MINLPProb& prob, MIP_NonLinearEval* &eval);
        
        //modelStatus should be a enum gmoModelStatus and solverStatus shoul be a enum gmoSolverStatus
        int setSolution( const int modelStatus, const int solverStatus, const double *sol );
    };
    
    
    
}

#endif	
