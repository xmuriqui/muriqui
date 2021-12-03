/*
 *
 */


#ifndef __NMC_NUMCOMP_HPP__
#define __NMC_NUMCOMP_HPP__


#include <cmath>


#include "NMC_config.hpp"



namespace numcomp
{


#define NMC_DEBUG_MODE 1



enum NMC_ERROR_CODES
{
	NMC_EC_EVAL_ERROR				= -1,
	NMC_EC_BAD_INPUT				= -2,
	NMC_EC_MAX_ITERATIONS_STOP 		= -3,
	NMC_EC_INTERNAL_ERROR			= -4
};


enum NMC_SUB_ERROR_CODES
{
	NMC_SEC_ROOT_FOUND 			= 1
};



class NMC_Function
{
	
public:
	
	unsigned int nevals;
	
	
	
	NMC_Function();
	
	virtual ~ NMC_Function();
	
	void resetEvalCounter();
	
	//evaluates and increase counter
	int evalIncCounter(double x, double &value);
	
	virtual int eval(double x, double &value) = 0;
};



class NMC_PointR2
{
	
public:
	
	double x;
	double y;
	
	
	
	NMC_PointR2(double x = NAN, double y = NAN);
	
	int calcy(NMC_Function &f, bool incCounter = true);
	
	NMC_PointR2 copy();
	
};




class NMC_ZeroOfFunctionBase
{
protected:
	
	//arguments in order
	int bracket(const NMC_PointR2 &a, const NMC_PointR2 &b, const NMC_PointR2 &c, NMC_PointR2 &bara, NMC_PointR2 &barb, NMC_PointR2 &d);
	
	void newtonQuadratic(const NMC_PointR2 &a, const NMC_PointR2 &b, const NMC_PointR2 &d, int k, double &rk);
	
	//arguments in order
	void ipzero(const NMC_PointR2 &a, const NMC_PointR2 &b, const NMC_PointR2 &c, const NMC_PointR2 &d, double &barx);
	
	
	bool hasSomeEqualPair(const double tolerance, unsigned int nfp, const double **fp);
	
	
public:
	
	unsigned int maxIters;
	double zeroTol;
	NMC_Function *f;
	
	
	NMC_ZeroOfFunctionBase(double zeroTol, NMC_Function *f, unsigned int maxIters);
	
	virtual ~NMC_ZeroOfFunctionBase();
	
	virtual int findZero(double begin, double end, NMC_PointR2 &root, unsigned int *niters = NULL) = 0;
	
};


// Reference: Alefeld, Potra & Shi. Algorithm 748 enclosing zeros of continuous functions. ACM Transactions on Mathematical Software, Vol. 21, No. 3, September 1995, Pages 327-344.
class NMC_ZeroOfFunction748_1 : public NMC_ZeroOfFunctionBase
{
	
	
public:
	
	double mu;
	
	
	NMC_ZeroOfFunction748_1(double zeroTol, NMC_Function *f, unsigned int maxIters);
	
	
	virtual int findZero(double begin, double end, NMC_PointR2 &root, unsigned int *niters = NULL) override;
	
};


// Reference: Alefeld, Potra & Shi. Algorithm 748 enclosing zeros of continuous functions. ACM Transactions on Mathematical Software, Vol. 21, No. 3, September 1995, Pages 327-344
class NMC_ZeroOfFunction748_2 : public NMC_ZeroOfFunctionBase
{
	
public:
	
	double mu;
	
	
	NMC_ZeroOfFunction748_2(double zeroTol, NMC_Function *f, unsigned int maxIters);
	
	
	virtual int findZero(double begin, double end, NMC_PointR2 &root, unsigned int *niters = NULL) override;
	
	
};














}






#endif
