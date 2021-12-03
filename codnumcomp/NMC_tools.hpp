

#ifndef __NMC_TOOLS_HPP__
#define __NMC_TOOLS_HPP__



#include "NMC_numComp.hpp"



#ifdef __FILE__
	#ifdef __LINE__
	    #define NMC_DEF_GETFILELINE 1
	#endif
#endif


#ifdef NMC_DEF_GETFILELINE
	     
	#define NMC_GETFILELINE  \
	     " on file: '" __FILE__  "' line: '" << __LINE__ <<  "' function: " << __func__
#else
    #define NMC_GETFILELINE ""
#endif


#define NMC_PREPRINT "numComp: "


#define NMC_PRINTMEMERROR std::cerr << NMC_PREPRINT "Memory error" << NMC_GETFILELINE << "\n"


#define NMC_PRINTERROR std::cerr << NMC_PREPRINT "Error" << NMC_GETFILELINE << "\n"


#define NMC_PRINTERRORNUMBER(number) std::cerr << NMC_PREPRINT "Error " << number << NMC_GETFILELINE << "\n"


#define NMC_PRINTCALLBACKERRORNUMBER(number) std::cerr << NMC_PREPRINT << "Callback function error " << number << NMC_GETFILELINE << "\n"


#define NMC_PRINTERRORMSG(msg) std::cerr << NMC_PREPRINT << msg << NMC_GETFILELINE << "\n"

#define NMC_PRINTERRORMSGP(msg, arg) std::cerr << NMC_PREPRINT << msg << arg << NMC_GETFILELINE << "\n"

#define NMC_PRINTMSG(msg) std::cout << NMC_PREPRINT << msg;


#define NMC_getchar()  NMC_PRINTERRORMSG("stopped in a getchar"), getchar()




namespace numcomp{


template <class myClass>
inline myClass NMC_abs(const myClass &v)
{
	return v >= 0 ? v : -v;
}



//evaluate f[a, b]
inline double NMC_fbrack2(const NMC_PointR2 &a, const NMC_PointR2 &b)
{
	return (b.y - a.y)/(b.x - a.x);
}


//evaluate f[a, b, d]
inline double NMC_fbrack3(const NMC_PointR2 &a, const NMC_PointR2 &b, const NMC_PointR2 &d)
{
	const double fab = NMC_fbrack2(a, b);
	const double fbd = NMC_fbrack2(b, d);
	
	return (fbd - fab)/(d.x - a.x);
}


//evaluate P(a, b, d)(x)
inline double NMC_P(const NMC_PointR2 &a, const NMC_PointR2 &b, const NMC_PointR2 &d, double x)
{
	const double fab = NMC_fbrack2(a, b);
	//const double fbd = _fbrack2(b, d);
	const double fabd = NMC_fbrack3(a, b, d);
	
	
	return a.y + fab*(x-a.x) + fabd*(x - a.x)*(x - b.x);
}



}




























#endif
