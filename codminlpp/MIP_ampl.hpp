




#ifndef MIP_AMPL_HPP
#define MIP_AMPL_HPP


#include "MIP_config.hpp"



#include "MIP_minlpProblem.hpp"








#if MIP_DEBUG_MODE
    #define MIP_AMPL_DEBUG_MODE  0
#else
    #define MIP_AMPL_DEBUG_MODE  0
#endif


#define MIP_READ_AMPL_PARAMS 1



extern "C"{
    //we have to include ASL after iostream (it is included in sparse matrix inside of MIP_minlpProblem.hpp). Otherwise, we will got an error.
    #if MIP_HAVE_ASL
        #include "asl.h"
        #include "getstub.h"

        #define asl cur_ASL
    #endif
};



namespace minlpproblem
{
    
    class MIP_NLEvalAmlp : public MIP_NonLinearEval
    {
        
        unsigned int nthreads;
        //int nzJac, nzObjHess, nzCtrHess;
        #if MIP_HAVE_ASL
            ASL **aslv;
        #endif
        //double **ctrAux;
        double **auxVals;
        
        #if MIP_AMPL_DEBUG_MODE
            bool checkJacInd, checkHessInd;
        #endif
        
        unsigned int nEvalErrors;
        
    public:
        bool quadObj, quadConstrs;
        char *stub;
        
        
        MIP_NLEvalAmlp(const bool quadObj, const bool quadConstrs, char* stub);
        
        
        virtual ~MIP_NLEvalAmlp();
        
        
        void desallocate();
        
        
        virtual int initialize(const int nthreads, const int n, const int m, const int nzJac, const int nzLagHess) override;
        
        
        virtual int eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value) override;
        
        
        virtual int eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *values) override;
        
        
        virtual int eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *ctrEval, const double *x, double *values) override;
        
        
        virtual int eval_grad_nl_constrs_part( const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *ctrEval, const double *x, MIP_SparseMatrix& J) override;
        
        
        virtual int eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double obj_factor, const double *lambda, MIP_SparseMatrix& hessian) override;
        
        
        virtual void finalize(const int nthreads, const int n, const int mnl, const int nzNLJac, const int nzNLLagHess) override;
        
    };
    
    
    #if !MIP_HAVE_ASL
        typedef void Option_Info;
    #endif
    
    class MIP_ReadAmplModel
    {
        
    public:
        
        #if MIP_HAVE_ASL
            ASL *asl;
        #endif	
        
        
        MIP_ReadAmplModel();
        
        virtual ~MIP_ReadAmplModel();
        
        void desallocate();
        
        bool isMaximizationProblem();
        
        void putInfoToAmpl(const char *msg, const int return_code, double *primalSol, double *dualSol);
        
        void readParameters(char **argv, Option_Info *opinfo);
        
        int readProblem(char* stub, MIP_MINLPProb& prob, MIP_NonLinearEval* &eval);
        
    };
    
    
    
};




#endif

































