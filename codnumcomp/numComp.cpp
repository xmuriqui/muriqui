


#include <iostream>


#include "NMC_numComp.hpp"
#include "NMC_tools.hpp"



using namespace numcomp;







NMC_Function::NMC_Function()
{
	nevals = 0;
}


NMC_Function::~ NMC_Function()
{
}


int NMC_Function::evalIncCounter(double x, double &value)
{
	nevals++;
	return eval(x, value);
}


void NMC_Function::resetEvalCounter()
{
	nevals = 0;
}





NMC_PointR2::NMC_PointR2(double x, double y)
{
	this->x = x;
	this->y = y;
}


int NMC_PointR2::calcy(NMC_Function &f, const bool incCounter)
{
	if(incCounter)
		return f.evalIncCounter(x, y);
	else
		return f.eval(x, y);
}


NMC_PointR2 NMC_PointR2::copy()
{
	return NMC_PointR2(x, y);
}







int NMC_ZeroOfFunctionBase::bracket(const NMC_PointR2 &a, const NMC_PointR2 &b, const NMC_PointR2 &c, NMC_PointR2 &bara, NMC_PointR2 &barb, NMC_PointR2 &d)
{
	
	if( NMC_abs(c.y) < zeroTol )
	{
		return NMC_SEC_ROOT_FOUND; //c is the zero of function
	}
	
	//std::cout << "\nbracket a.x: " << a.x << " a.y: " << a.y << " b.x: " << b.x << " b.y: " << b.y << " c.x: " << c.x << " c.y: " << c.y ;
	
	if( a.y*c.y < 0.0 )
	{
		bara = a;
		barb = c;
		d = b;
	}
	else 
	{
		if( !(b.y*c.y < 0.0) )
		{ //we have a serious error here
			
			NMC_PRINTERRORMSG("Serious error about interval reduction.");
			
			return NMC_EC_INTERNAL_ERROR;
		}
		
		bara = c;
		barb = b;
		d = a; 
	}
	
	//std::cout << " d.x: " << d.x << " d.y: " << d.y;
	//std::cout << "\n\n";
	
	
	return 0;
}


void NMC_ZeroOfFunctionBase::newtonQuadratic(const NMC_PointR2 &a, const NMC_PointR2 &b, const NMC_PointR2 &d, int k, double &rk)
{
	const double fab = NMC_fbrack2(a, b);
	//const double fbd = _fbrack2(b, d);
	const double fabd = NMC_fbrack3(a, b, d);
	const double &A = fabd, &B = fab;
	
	
	//std::cout << "newton a.x: " << a.x << " a.y: " << a.y << " b.x: " << b.x << " b.y: " << b.y << " d.x: " << d.x << " d.y: " << d.y << "           A: " << A << " B: " << B << " fbd: " << fbd;
	
	
	if( NMC_abs(A) <= zeroTol )
	{
		rk = a.x - a.y/B;
	}
	else
	{
		if( A*a.y > 0 )
			rk = a.x;
		else
			rk = b.x;
		
		
		for(int i = 0; i < k; i++)
		{
			//rk = rk - B*rk/(B + A*(2*rk - a.x - b.x));
			rk = rk - NMC_P(a, b, d, rk)/(B + A*(2*rk - a.x - b.x));
		}
	}
	
	//std::cout << " rk: " << rk << "\n";
	
	//value = rk;
	//return rk;
}


//arguments in order
void NMC_ZeroOfFunctionBase::ipzero(const NMC_PointR2 &a, const NMC_PointR2 &b, const NMC_PointR2 &c, const NMC_PointR2 &d, double &barx)
{
	double q11, q21, q31;
	double d21, d31;
	double q22, q32, d32, q33;
	
	
	q11 = (c.x - d.x)* c.y/(d.y - c.y);
	q21 = (b.x - c.x)* b.y/(c.y - b.y);
	q31 = (a.x - b.x)* a.y/(b.y - a.y);
	
	
	d21 = (b.x - c.x)* c.y/(c.y - b.y);
	d31 = (a.x - b.x)* b.y/(b.y - a.y);
	
	q22 = (d21 - q11)* b.y/(d.y - b.y);
	q32 = (d31 - q21)* a.y/(c.y - a.y);
	
	d32 = (d31 - q21)* c.y/(c.y - a.y);
	q33 = (d32 - q22)* a.y/(d.y - a.y);
	
	barx = a.x + (q31 + q32 + q33);
}


bool NMC_ZeroOfFunctionBase::hasSomeEqualPair(const double tolerance, unsigned int nfp, const double **fp)
{
	for(unsigned int i = 0; i < nfp; i++)
	{
		for(unsigned int j = i+1; j < nfp; j++)
		{
			if( NMC_abs(*fp[i] - *fp[j]) <= tolerance )
			{
				return true;
				//equalf = true;
				//breaking the loops
				//i = nfp;
				//break;
			}
		}
	}
	
	return false;
}



NMC_ZeroOfFunctionBase::NMC_ZeroOfFunctionBase(double zeroTol, NMC_Function *f, unsigned int maxIters)
{
	this->zeroTol = zeroTol;
	this->f = f;
	this->maxIters = maxIters;
}


NMC_ZeroOfFunctionBase::~NMC_ZeroOfFunctionBase()
{
}









NMC_ZeroOfFunction748_1::NMC_ZeroOfFunction748_1(double zeroTol, NMC_Function *f, unsigned int maxIters) : NMC_ZeroOfFunctionBase(zeroTol, f, maxIters)
{
	mu = 0.5;
}


int NMC_ZeroOfFunction748_1::findZero(double begin, double end, NMC_PointR2 &root, unsigned int *niters)
{
	int ret;
	unsigned int myiters;
	
	
	NMC_PointR2 a, b, c, d, e;
	
	const int nfp = 4;
	const double* fp[nfp] = {&a.y, &b.y, &d.y, &e.y};
	
	double fab;
	
	NMC_PointR2 bara, barb, barc, bard;
	NMC_PointR2 hata, hatb, hatc, hatd;
	NMC_PointR2 u;
	
	
	
	if(niters == NULL)
		niters = &myiters;
	
	*niters = 1;
	
	
	a.x = begin;
	b.x = end;
	
	
	ret = a.calcy(*f);
	
	if( ret != 0 )
	{
		#if NMC_DEBUG_MODE
			NMC_PRINTCALLBACKERRORNUMBER(ret);
		#endif
		
		return NMC_EC_EVAL_ERROR;
	}
	
	ret = b.calcy(*f);
	if( ret != 0 )
	{
		#if NMC_DEBUG_MODE
			NMC_PRINTCALLBACKERRORNUMBER(ret);
		#endif
		
		return NMC_EC_EVAL_ERROR;
	}
	
	
	if( !(a.y*b.y < 0) ) //we use ! to treat nan case
	{
		#if NMC_DEBUG_MODE
			std::cout << NMC_PREPRINT "Function must have oposite signals in the bounds of interval. a.x: " << a.x << " a.y: " << a.y << " b.x: " << b.x << " b.y: " << b.y << "\n";
		#endif
		return NMC_EC_BAD_INPUT;
	}
	
	
	fab = NMC_fbrack2(a, b);
	
	c.x = a.x - a.y/fab;
	ret = c.calcy(*f);
	
	if( ret != 0 )
	{
		#if NMC_DEBUG_MODE
			NMC_PRINTCALLBACKERRORNUMBER(ret);
		#endif
		
		return NMC_EC_EVAL_ERROR;
	}
	
	
	//std::cout << "iaia 1\n";
	ret = bracket(a.copy(), b.copy(), c, a, b, d);
	if(ret != 0)
	{
		if(ret == NMC_SEC_ROOT_FOUND)
		{
			root = c;
			return 0;
		}
		
		return NMC_EC_INTERNAL_ERROR;
	}
	
	
	while(*niters < maxIters)
	{
		bool equalf = false;
		double fbarabarb;
		
		
		(*niters)++;
		
		
		equalf = hasSomeEqualPair(zeroTol, nfp, fp);
		
		
		if( *niters == 2 || equalf )
		{
			newtonQuadratic(a, b, d, 2, c.x);
			
			ret = c.calcy(*f);
			if(ret != 0)
			{
				#if NMC_DEBUG_MODE
					NMC_PRINTCALLBACKERRORNUMBER(ret);
				#endif
				
				return NMC_EC_EVAL_ERROR;
			}
		
			//std::cout << "ola 0 c.x: " << c.x << " c.y: " << c.y << "\n";
		}
		else
		{
			//std::cout << "ola 1 a.x: " << a.x << " a.y: " << a.y << " b.x: " << b.x << " b.y: " << b.y <<  " c.x: " << c.x << " c.y: " << c.y << " d.x: " << d.x << " d.y: " << d.y << " e.x: " << e.x << " e.y: " << e.y << "\n";
			
			ipzero(a, b, d, e, c.x);
			
			if(c.calcy(*f) != 0)
			{
				#if NMC_DEBUG_MODE
					NMC_PRINTCALLBACKERRORNUMBER(ret);
				#endif
				
				return NMC_EC_EVAL_ERROR;
			}
			
			//std::cout << "ola 2 c.x: " << c.x << " c.y: " << c.y << "\n";
			
			if( (c.x - a.x)*(c.x - b.x) >= 0.0 )
			{
				newtonQuadratic(a, b, d, 2, c.x);
				
				int ret = c.calcy(*f);
				
				if(ret != 0)
				{
					#if NMC_DEBUG_MODE
						NMC_PRINTCALLBACKERRORNUMBER(ret);
					#endif
					
					return NMC_EC_EVAL_ERROR;
				}
			}
			
			//std::cout << "ola 3 c.x: " << c.x << " c.y: " << c.y << "\n";
		}
		
		
		
		
		//std::cout << "iaia 2\n";
		
		ret = bracket(a, b, c, bara, barb, bard);
		if(ret != 0)
		{
			if(ret == NMC_SEC_ROOT_FOUND)
			{
				root = c;
				return 0;
			}
			
			return NMC_EC_INTERNAL_ERROR;
		}
		
		if( NMC_abs(bara.y) < NMC_abs(barb.y) )
			u = bara;
		else
			u = barb;
		
		
		fbarabarb = NMC_fbrack2(bara, barb);
		
		
		barc.x = u.x - 2.0*u.y/fbarabarb;
		
		ret = barc.calcy(*f);
		if( ret != 0) 
		{
			#if NMC_DEBUG_MODE
				NMC_PRINTCALLBACKERRORNUMBER(ret);
			#endif
			
			return NMC_EC_EVAL_ERROR;
		}
		
		
		
		if( NMC_abs(barc.x - u.x) > 0.5*(barb.x - bara.x) )
		{
			hatc = 0.5*(barb.x + bara.x);
			
			ret = hatc.calcy(*f);
			if( ret != 0)
			{
				#if NMC_DEBUG_MODE
					NMC_PRINTCALLBACKERRORNUMBER(ret);
				#endif
				
				return NMC_EC_EVAL_ERROR;
			}
			
		}
		else
		{
			hatc = barc;
		}
		
		if( NMC_abs(hatc.y) <= zeroTol )
		{
			root = hatc;
			return 0;
		}
		
		
		//std::cout << "iaia 3\n";
		
		ret = bracket(bara, barb, hatc, hata, hatb, hatd);
		if(ret != 0)
		{
			if(ret == NMC_SEC_ROOT_FOUND)
			{
				root = hatc;
				return 0;
			}
			
			return NMC_EC_INTERNAL_ERROR;
		}
		
		
		if( hatb.x - hata.x < mu*(b.x - a.x) )
		{
			a = hata;
			b = hatb;
			d = hatd;
			e = bard;
		}
		else
		{
			NMC_PointR2 t(0.5*(hata.x + hatb.x) );
			
			ret = t.calcy(*f);
			if( ret != 0 )
			{
				#if NMC_DEBUG_MODE
					NMC_PRINTCALLBACKERRORNUMBER(ret);
				#endif
				
				return NMC_EC_EVAL_ERROR;
			}
			
			
			e = hatd;
			
			
			//std::cout << "iaia 4\n";
			
			ret = bracket(hata, hatb, t, a, b, d);
			if(ret != 0)
			{
				if(ret == NMC_SEC_ROOT_FOUND)
				{
					root = t;
					return 0;
				}
				
				return NMC_EC_INTERNAL_ERROR;
			}
			
			
		}
		
	}
	
	
	//if(iter >= maxIters)
	{
		return NMC_EC_MAX_ITERATIONS_STOP;
	}
	
	
	return 0;
}








NMC_ZeroOfFunction748_2::NMC_ZeroOfFunction748_2(double zeroTol, NMC_Function *f, unsigned int maxIters) : NMC_ZeroOfFunctionBase(zeroTol, f, maxIters)
{ 
	mu = 0.5;
}


int NMC_ZeroOfFunction748_2::findZero(double begin, double end, NMC_PointR2 &root, unsigned int *niters)
{
	int ret;
	unsigned int myiters;
	
	
	NMC_PointR2 a, b, c, d, e;
	
	const int nfp = 4;
	const double* fp[nfp] = {&a.y, &b.y, &d.y, &e.y};
	
	NMC_PointR2 bara, barb, barc, bard, bare;
	NMC_PointR2 hata, hatb, hatc, hatd;
	NMC_PointR2 u;
	const int nfpbar = nfp;
	const double* fpbar[nfpbar] = {&bara.y, &barb.y, &bard.y, &bare.y};
	
	
	
	if(niters == NULL)
		niters = &myiters;
	
	*niters = 1;
	
	
	a.x = begin;
	b.x = end;
	
	
	ret = a.calcy(*f);
	
	if( ret != 0 )
	{
		#if NMC_DEBUG_MODE
			NMC_PRINTCALLBACKERRORNUMBER(ret);
		#endif
		
		return NMC_EC_EVAL_ERROR;
	}
	
	
	ret = b.calcy(*f);
	if( ret != 0 )
	{
		#if NMC_DEBUG_MODE
			NMC_PRINTCALLBACKERRORNUMBER(ret);
		#endif
		
		return NMC_EC_EVAL_ERROR;
	}
	
	
	if( !(a.y*b.y < 0) ) //we use ! to treat nan case
	{
		#if NMC_DEBUG_MODE
			std::cout << "Function must have oposite signals in the bounds of interval. a.x: " << a.x << " a.y: " << a.y << " b.x: " << b.x << " b.y: " << b.y << "\n";
		#endif
		
		return NMC_EC_BAD_INPUT;
	}
	
	
	c.x = a.x - a.y/NMC_fbrack2(a, b);
	
	ret = c.calcy(*f);
	if( ret != 0 )
	{
		#if NMC_DEBUG_MODE
			NMC_PRINTCALLBACKERRORNUMBER(ret);
		#endif
		
		return NMC_EC_EVAL_ERROR;
	}
	
	
	//std::cout << "iaia 1\n";
	ret = bracket(a.copy(), b.copy(), c, a, b, d);
	if(ret != 0)
	{
		if(ret == NMC_SEC_ROOT_FOUND)
		{
			root = c;
			return 0;
		}
		
		return NMC_EC_INTERNAL_ERROR;
	}
	
	
	//std::cout << "opa 1   a.x: " << a.x << " a.y: " << a.y << " b.x: " << b.x << " b.y: " << b.y <<  " c.x: " << c.x << " c.y: " << c.y << " d.x: " << d.x << " d.y: " << d.y << " e.x: " << e.x << " e.y: " << e.y << "\n";
	
	
	
	
	while(*niters < maxIters)
	{
		bool equalf;
		double fbarabarb;
		
		
		(*niters)++;
		
		//std::cout << "iter: " << *niters << " interval: [" << a.x << " " << b.x << "]\n";
		
		
		equalf = hasSomeEqualPair(zeroTol, nfp, fp);
		
		
		if( *niters == 2 || equalf )
		{
			newtonQuadratic(a, b, d, 2, c.x);
			
			int ret = c.calcy(*f);
			if(ret != 0)
			{
				#if NMC_DEBUG_MODE
					NMC_PRINTCALLBACKERRORNUMBER(ret);
				#endif
				
				return NMC_EC_EVAL_ERROR;
			}
		
			//std::cout << "ola 0 c.x: " << c.x << " c.y: " << c.y << "\n";
		}
		else
		{
			//std::cout << "ola 1 a.x: " << a.x << " a.y: " << a.y << " b.x: " << b.x << " b.y: " << b.y <<  " c.x: " << c.x << " c.y: " << c.y << " d.x: " << d.x << " d.y: " << d.y << " e.x: " << e.x << " e.y: " << e.y << "\n";
			
			ipzero(a, b, d, e, c.x);
			
			ret = c.calcy(*f);
			if(ret != 0)
			{
				#if NMC_DEBUG_MODE
					NMC_PRINTCALLBACKERRORNUMBER(ret);
				#endif
				
				return NMC_EC_EVAL_ERROR;
			}
			
			//std::cout << "ola 2 c.x: " << c.x << " c.y: " << c.y << "\n";
			
			if( (c.x - a.x)*(c.x - b.x) >= 0.0 )
			{
				newtonQuadratic(a, b, d, 2, c.x);
				
				ret = c.calcy(*f);
				if(ret != 0)
				{
					#if NMC_DEBUG_MODE
						NMC_PRINTCALLBACKERRORNUMBER(ret);
					#endif
					
					return NMC_EC_EVAL_ERROR;
				}
			}
			
			//std::cout << "ola 3 c.x: " << c.x << " c.y: " << c.y << "\n";
		}
		
		bare = d;
		
		//std::cout << "iaia 2\n";
		//std::cout << "a.x: " << a.x << " a.y: " << a.y << " b.x: " << b.x << " b.y: " << b.y <<  " c.x: " << c.x << " c.y: " << c.y << " d.x: " << d.x << " d.y: " << d.y << " e.x: " << e.x << " e.y: " << e.y << "\n";
		ret = bracket(a, b, c, bara, barb, bard);
		if(ret != 0)
		{
			if(ret == NMC_SEC_ROOT_FOUND)
			{
				root = c;
				return 0;
			}
			
			return NMC_EC_INTERNAL_ERROR;
		}
		
		
		equalf = hasSomeEqualPair(zeroTol, nfpbar, fpbar);
		
		if(equalf)
		{
			newtonQuadratic(bara, barb, bard, 3, barc.x);
			
			int ret = barc.calcy(*f);
			
			if(ret != 0)
			{
				#if NMC_DEBUG_MODE
					NMC_PRINTCALLBACKERRORNUMBER(ret);
				#endif
				
				return NMC_EC_EVAL_ERROR;
			}
			
			//std::cout << "ola 0 barc.x: " << barc.x << " barc.y: " << barc.y << "\n";
		}
		else
		{
			ipzero(bara, barb, bard, bare, barc.x);
			
			ret = barc.calcy(*f);
			if(ret != 0)
			{
				#if NMC_DEBUG_MODE
					NMC_PRINTCALLBACKERRORNUMBER(ret);
				#endif
				
				return NMC_EC_EVAL_ERROR;
			}
			
			//std::cout << "ola 2 c.x: " << c.x << " c.y: " << c.y << "\n";
			
			if( (barc.x - bara.x)*(barc.x - barb.x) >= 0.0 )
			{
				newtonQuadratic(bara, barb, bard, 3, barc.x);
				
				ret = barc.calcy(*f);
				if(ret != 0)
				{
					#if NMC_DEBUG_MODE
						NMC_PRINTCALLBACKERRORNUMBER(ret);
					#endif
					
					return NMC_EC_EVAL_ERROR;
				}
			}
			
			//std::cout << "ola 3 barc.x: " << barc.x << " barc.y: " << barc.y << "\n";
		}
		
		//std::cout << "iaia 3\n";
		ret = bracket(bara.copy(), barb.copy(), barc, bara, barb, bard);
		if(ret != 0)
		{
			if(ret == NMC_SEC_ROOT_FOUND)
			{
				root = barc;
				return 0;
			}
			
			return NMC_EC_INTERNAL_ERROR;
		}
		
		
		if( NMC_abs(bara.y) < NMC_abs(barb.y) )
			u = bara;
		else
			u = barb;
		
		
		fbarabarb = NMC_fbrack2(bara, barb);
		
		barc.x = u.x - 2.0*u.y/fbarabarb;
		
		ret = barc.calcy(*f);
		if( ret != 0) 
		{
			#if NMC_DEBUG_MODE
				NMC_PRINTCALLBACKERRORNUMBER(ret);
			#endif
			
			return NMC_EC_EVAL_ERROR;
		}
		
		
		if( NMC_abs(barc.x - u.x) > 0.5*(barb.x - bara.x) )
		{
			hatc = 0.5*(barb.x + bara.x);
			
			ret = hatc.calcy(*f);
			if( ret != 0)
			{
				#if NMC_DEBUG_MODE
					NMC_PRINTCALLBACKERRORNUMBER(ret);
				#endif
				
				return NMC_EC_EVAL_ERROR;
			}
		}
		else
		{
			hatc = barc;
		}
		
		
		//std::cout << "iaia 4\n";
		ret = bracket(bara, barb, hatc, hata, hatb, hatd);
		if(ret != 0)
		{
			if(ret == NMC_SEC_ROOT_FOUND)
			{
				root = hatc;
				return 0;
			}
			
			return NMC_EC_INTERNAL_ERROR;
		}
		
		
		if( hatb.x - hata.x < mu*(b.x - a.x) )
		{
			a = hata;
			b = hatb;
			d = hatd;
			e = bard;
		}
		else
		{
			NMC_PointR2 t(0.5*(hata.x + hatb.x) );
			
			ret = t.calcy(*f);
			if( ret != 0 )
			{
				#if NMC_DEBUG_MODE
					NMC_PRINTCALLBACKERRORNUMBER(ret);
				#endif
				
				return NMC_EC_EVAL_ERROR;
			}
			
			e = hatd;
			
			
			//std::cout << "iaia 5\n";
			ret = bracket(hata, hatb, t, a, b, d);
			if(ret != 0)
			{
				if(ret == NMC_SEC_ROOT_FOUND)
				{
					root = t;
					return 0;
				}
				
				return NMC_EC_INTERNAL_ERROR;
			}
			
			
		}
		
		
	}
	
	
	return NMC_EC_MAX_ITERATIONS_STOP;
}











