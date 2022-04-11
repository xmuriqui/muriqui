
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cassert>
#include <climits>
#include <cstring>

#include <iostream>
#include <new>

#include "muriqui.hpp"
#include "MRQ_tools.hpp"
//#include "MRQ_nlpSolvers.hpp"
#include "MRQ_solvers.hpp"

#include "BBL_tools.hpp"
#include "NMC_numComp.hpp"


using namespace std;
using namespace optsolvers;
using namespace muriqui;
using namespace numcomp;



#if MRQ_BB_SUPER_THREAD_DEBUG_MODE
    std::map< std::thread::id, std::ostream* >  muriqui::MRQ_thsOut;
    MRQ_Mutex SEMAPH_thsOut;
    
    
    void muriqui::MRQ_createFileToThread( const std::thread::id &tid )
    {
        std::stringstream fname;
        fname << "mrq_thread_output_" << tid << ".dbg";
        
        
        SEMAPH_thsOut.lock(2);
        
        if( MRQ_thsOut.count( tid ) == 0 )
        {
            std::ofstream *fout = new std::ofstream;
            
            fout->open( fname.str() );
            
            MRQ_thsOut[ tid ] = fout;
        }
        
        SEMAPH_thsOut.unlock(2);
    }
    
    
    void muriqui::MRQ_closeFileToThread( const std::thread::id &tid )
    {
        
        SEMAPH_thsOut.lock(2);
        
        if( MRQ_thsOut.count( tid ) > 0 )
        {
            auto p = MRQ_thsOut[tid];
            MRQ_thsOut.erase(tid); //removing the file from map. In this, it will be created again in the next execution of some function being debuged
            delete p; //it will call close method (I hope)
        }
        
        SEMAPH_thsOut.unlock(2);
    }
    
    
#endif




int muriqui::MRQ_getNumCores(void)
{
    return branchAndBound::BBL_getNumCores();
}



double muriqui::MRQ_getTime(void)
{
    return branchAndBound::BBL_getTime();
}





inline void MRQ_hello(std::ostream &out = std::cout)
{
    out << "\n----------------------------------------------------------------------------------\n"
    "Muriqui MINLP optimizer, version " MRQ_VERSION " compiled at " __DATE__ " " __TIME__ ".\n"
    "By " MRQ_AUTHOR ", " MRQ_AUTHOR_FILIATION  "\n"
    "collaborators: " MRQ_COLLABORATORS "\n\n"
    "if you use, please cite: \n    " MRQ_BASE_PAPER
    "\n----------------------------------------------------------------------------------\n\n";
}




void muriqui::MRQ_helloOnceTime(std::ostream &out)
{
    static bool printed = false;
    
    if(!printed)
    {
        MRQ_hello();
        printed = true;
    }
}



int muriqui::MRQ_tryReadAlgChoiceFile(const char *algChoiceFileName, int printLevel, MRQ_ALG_CODE &readAlgCode)
{
    FILE *algFile;
    
    
    if( printLevel > 0 )
        std::cout << MRQ_PREPRINT "Trying to read algorithm choice file " << algChoiceFileName << " . ";
    
    
    algFile = fopen( MRQ_MURIQUI_ALG_CHOICE_FILE, "r" );
    if( algFile )
    {
        const int maxSizeOption = 256;
        char option[maxSizeOption], option2[maxSizeOption];
        
        option[0] = MRQ_CHARAC_COMENT_ON_PARAMS_FILE;
        option[1] = '\0';
        
        option2[0] = MRQ_CHARAC_COMENT_ON_PARAMS_FILE;
        option2[1] = '\0';
        
        do
        {
            //fscanf( algFile, "%100[^\n]s", option );
            char *p = fgets( option, maxSizeOption, algFile  );
            if(p == NULL)
                break;
        }while( !feof(algFile) && (option[0] == MRQ_CHARAC_COMENT_ON_PARAMS_FILE || OPT_isEmptyString(option) ) );
        
        
        //removing the \n and whitespaces
        sscanf(option, "%s", option2);
        
        
        fclose(algFile);
        algFile = NULL;
        
        if(printLevel > 0)
            std::cout << "Option read: " << option2 << ". Trying set."; 
        
        const int r = MRQ_strToEnum(option2, readAlgCode);
        if( r == 0 )
        {
            if(printLevel > 0)
                std::cout << " Success!\n";
            return 0;
        }
        else
        {
            if(printLevel > 0)
                std::cout << " Invalid code!\n";
            return MRQ_VALUE_ERROR;
        }
    }
    else
    {
        if(printLevel > 0)
            std::cout << "Not found!\n";
        return MRQ_NAME_ERROR;
    }
    
}



int muriqui::MRQ_writeProblemParameters(const char *probName, MRQ_MINLPProb &prob, FILE *outFile )
{
    int ml, mq, mnl;
    int ret;
    const int ny = prob.getNumberOfIntegerVars();
    const int nbin = prob.getNumberOfBinaryVars();
    
    prob.getConstraintStatistcs(&ml, &mq, &mnl);
    
    
    ret = fprintf(outFile, "\n%s%s%s%s", probName, MRQ_CHAR_SEP, prob.hasNlObj ? "nl": ( prob.Q.getNumberOfElements() == 0 ? "linear" : "quad"), MRQ_CHAR_SEP );
    
    
    //nx, ny, m, m non linear, m quadratic, m linear, nzQ, nzJac, nzLagHess, objective (0 - min, 1 - max)
    ret += fprintf(outFile, "%d%s%d%s%d%s%d%s%d%s%d%s%d%s%d%s%d%s%d%s%d%s", prob.n -ny, MRQ_CHAR_SEP, ny, MRQ_CHAR_SEP, nbin, MRQ_CHAR_SEP, prob.m, MRQ_CHAR_SEP, mnl, MRQ_CHAR_SEP, mq, MRQ_CHAR_SEP,  ml, MRQ_CHAR_SEP,  prob.Q.getNumberOfElements() , MRQ_CHAR_SEP, prob.J.getNumberOfElements() , MRQ_CHAR_SEP, prob.lagH.getNumberOfElements() , MRQ_CHAR_SEP, (int) (prob.objFactor < 0.0) , MRQ_CHAR_SEP);
    
    fflush(outFile);
    
    return ret;
}





MRQ_Random::MRQ_Random()
{
}


MRQ_Random::MRQ_Random(const long int seed)
{
    setSeed(&seed);
}


long int MRQ_Random::getSeed()
{
    return seed;
}


long int MRQ_Random::setSeed(const long int *seed)
{   
    this->seed = seed ? *seed : time(NULL) ;
    
    #if MRQ_HAVE_CPP_2011
        gen.seed(this->seed);
    #else
        srand( this->seed );
    #endif
    
    //this->seed = s;
    
    return this->seed;
}


//generates a random integer in a interval [begin  end]
int MRQ_Random::randInt(const int begin, const int end)
{
    int v;
    
    #if MRQ_HAVE_CPP_2011
        std::uniform_int_distribution<int> dist(begin, end);
        v = dist(gen);
    #else
        v = begin + rand()%(end - begin + 1);
    #endif
    
    return v;
}


//generates a random integer in a interval [begin  end]
unsigned int MRQ_Random::randUInt(const unsigned int begin, const unsigned int end)
{
    int v;
    
    #if MRQ_HAVE_CPP_2011
        std::uniform_int_distribution<unsigned int> dist(begin, end);
        v = dist(gen);
    #else
        v = begin + rand()%(end - begin + 1);
    #endif
    
    return v;
}


//generates true with probability prob
bool MRQ_Random::randBool(const double prob)
{
    return random() <= prob;
}


//generates a random real in the interval [0 1)
double MRQ_Random::random()
{
    double v;
    
    #if MRQ_HAVE_CPP_2011
        static std::uniform_real_distribution<double> distUni(0.0, 1.0);
        v = distUni(gen);
    #else
        v = ((double) rand())/((double) RAND_MAX);
    #endif
    
    return v;
}


//generates a random real in the interval [begin end)
double MRQ_Random::random(const double begin, const double end)
{
    return begin + (end - begin)*random();
}



//generates a normal random real in the interval [0 1)
double MRQ_Random::randomNormal()
{
    double v;
    
    #if MRQ_HAVE_CPP_2011
        static std::normal_distribution<double> distNormal(0.0, 1.0);
        v = distNormal(gen);
    #else
        //it is not ideal, but...
        v = ((double) rand())/((double) RAND_MAX);
    #endif
    
    return v;
}


//generates a random real in the interval [begin end)
double MRQ_Random::randomNormal(const double begin, const double end)
{
    return begin + (end - begin)*randomNormal();
}



/*MRQ_GradientsEvaluation::MRQ_GradientsEvaluation(const int thnumber, MRQ_MINLPProb* prob)
{
    initialize(thnumber, prob);
}*/


void MRQ_GradientsEvaluation::desallocate()
{
    J.desallocateMemory();
    prob = NULL;
}


int MRQ_GradientsEvaluation::initialize(const unsigned int thnumber, MRQ_MINLPProb* prob)
{
    this->thnumber = thnumber;
    this->prob = prob;
    
    const int r = J.copyStructureFrom(prob->J);
    if(r != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
        
        return MRQ_MEMORY_ERROR;
    }
    
    return 0;
}


int MRQ_GradientsEvaluation::evaluateJacobian(const bool newx, const bool* constrEval, const double* x)
{
    const int r = prob->nlJacobianEval(thnumber, newx, constrEval, x, J);
    
    #if MRQ_DEBUG_MODE
    if(r != 0)
    {
        MRQ_PRINTERRORNUMBER(r);
    }
    #endif
    
    return r;
}


void MRQ_GradientsEvaluation::constraintCompleteGradient( const int constr, const double *x, double *grad)
{
    MIP_constrCompleteGrad(*prob, J, constr, x, grad, true);
}







MRQ_LAAPointsStoring::MRQ_LAAPointsStoring(const unsigned int n)
{
    initialize(n);
}


MRQ_LAAPointsStoring::~MRQ_LAAPointsStoring()
{
    desallocate();
}


int MRQ_LAAPointsStoring::allocateMorePoints( const unsigned int nnewpoints )
{
    //bool *auxb;
    int *auxui;
    int r;
    
    r = points.addNewRows(nnewpoints);
    //r += objAp.addNewRows(nnewpoints);
    
    //auxb = (bool *) realloc( used, (nallocated + nnewpoints) * sizeof(bool)  );
    auxui = (int *) realloc( indMaster, (nallocated + nnewpoints) * sizeof(int *) );
    
    if( r != 0 || !auxui )
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTMEMERROR;
        #endif
        return MRQ_MEMORY_ERROR;
    }
    
    //used = auxb;
    indMaster = auxui;
    
    r = nallocated + nnewpoints;
    for(int i = nallocated; i < r; i++)
        indMaster[i] = INT_MIN; //used[i] = false;
    
    
    nallocated += nnewpoints;
    
    return 0;
}


void MRQ_LAAPointsStoring::desallocate()
{
    points.desallocateMemory();
    //objAp.desallocateMemory();
    
    //MRQ_secFree(used);
    MRQ_secFree(indMaster);
    
    npoints = 0;
    nallocated = 0;
}


void MRQ_LAAPointsStoring::initialize(const unsigned int n)
{
    blockSize = 100;
    npoints = 0;
    nallocated = 0;
    nObjLinRem = 0u;
    //used = NULL;
    indMaster = NULL;
    
    points.setNumberOfColumns(n);
    //objAp.setNumberOfColumns(n+1);
}


//n is the dimention of point. Note objApp is the linear approximation of objective fucntion in point, and it should have dimension n+1, and the last position is the RHS of the linear approximation
int MRQ_LAAPointsStoring::addPoint( unsigned int n, const double *point)//, const double *objApp )
{
    int r;
    
    #if MRQ_DEBUG_MODE
        assert( npoints <= (int) nallocated );
    #endif
    
    if( npoints == (int) nallocated )
    {
        r = allocateMorePoints(blockSize);
        if( r != 0 )
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERRORNUMBER(r);
            #endif
            
            return MRQ_MEMORY_ERROR;
        }
    }
    
    
    r = points.setRowStructureAndValues( npoints, point, n );
    //r += objAp.setRowStructureAndValues( npoints, objApp, n+1 );
    
    if( r != 0 )
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(r);
        #endif
        
        return MRQ_MEMORY_ERROR;
    }
    
    
    npoints++;
    
    return 0;
}




//that method should be called when we have an update in zu value. Here, we check if we can remove linearization on some points using the linearization on other points...
int MRQ_LAAPointsStoring::updateObjLinearizationsByNonObjCutPoints2( optsolvers::OPT_LPSolver &master, const double zu, double *auxVals )
{
    const int n = points.getNumberOfColumns();
    int r, code = 0;
    double lb, ub;
    
    
    for(int i = npoints-2; i >= 0; i--)
    {
        if( indMaster[i] < 0)
            continue; //if indMilp is negative, linearization on point i is not in master problem milp more...
        
        //points.getFullRow(i, auxVals);
        
        for(int j = i+1; j < npoints; j++)
        {
            if( indMaster[j] < 0 )
                continue;
            
            int &ind = indMaster[j];
            
            r = master.getFullConstraintLinearPart( ind, auxVals );
            
            r+= master.getConstraintBounds(ind, lb, ub);
            
            if( r != 0 )
            {
                #if MRQ_DEBUG_MODE
                    MRQ_PRINTERRORNUMBER(r);
                #endif
                
                code = MRQ_MILP_SOLVER_ERROR;
                continue;
            }
            
            #if MRQ_DEBUG_MODE
                assert( lb <= -OPT_INFINITY && ub < OPT_INFINITY );
            #endif
            
            
            const double pointValue = points.rowEvaluation(i, auxVals) - ub;//  points[i].evalTimesxt(auxVals) - ub;
            
            
            
            r = _checkByNonObjCut( master, zu, n, indMaster[i], pointValue );
            
            if( r != 0 )
            {
                #if MRQ_DEBUG_MODE
                    MRQ_PRINTERRORNUMBER(r);
                #endif
                
                code = r;
            }
            
            if( indMaster[i] < 0 )
            {
                //std::cout << "removing obj linearization of point " << i << std::endl;
                //getchar();
                
                break;
            }
        }
        
    }
    
    return code;
}



int MRQ_LAAPointsStoring:: checkByNonObjCut(optsolvers::OPT_LPSolver &master, const double zu, const int n, const double *objLin, const double rhs, int &indMaster, const double *point)
{
    double pointValue = 0.0;
    
    
    for(int j = 0; j < n; j++)
        pointValue += objLin[j]*point[j];
    
    pointValue -= rhs; //we perform multiplication separated to avoid numerical problems
    
    
    return _checkByNonObjCut( master, zu, n, indMaster, pointValue);
}


//Point value is the approximation of objective function...
int MRQ_LAAPointsStoring::_checkByNonObjCut(optsolvers::OPT_LPSolver &master, const double zu, const int n, int &indMaster, const double pointValue)
{
    int r, code = 0;
    
    
    
    if( pointValue >= zu )
    {
        //we can descart point i because when we apply this point in the linearisation on point j, it overcome the upper bound.
        
        //If we removed the constraint, we would need upadate indMIlp of other constraints. And the solver would have to mo more opeartions about its internal indices. So, we remove all coeficients and constraint get empty...
        r = master.resetConstraintLinearPart( indMaster, 0, NULL, NULL );
        
        r += master.setConstraintBounds( indMaster, -OPT_INFINITY, OPT_INFINITY );
        
        if( r != 0 )
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERRORNUMBER(r);
            #endif
            
            code = MRQ_MILP_SOLVER_ERROR;
        }
        
        indMaster = -indMaster -1;
        
        nObjLinRem++;
    }
    
    
    return code;
}


int MRQ_LAAPointsStoring:: updateObjLinearizationsByNonObjCutPointsByNewPoint( optsolvers::OPT_LPSolver &master, const double zu, const double *objLin)
{
    
    const unsigned int n = points.getNumberOfColumns();
    const double rhs = objLin[n];
    
    int r, code = 0;
    
    
    
    if( zu >= MRQ_INFINITY )
        return 0;
    
    
    for( int i = 0; i < npoints; i++ )
    {
        if( indMaster[i] < 0 )
            continue; //if indMilp is negative, linearization on point i is not in master problem milp more...
        
        //points.getFullRow( i, auxVals);
        
        const double pointValue = points.rowEvaluation(i, objLin) - rhs; //points[i].evalTimesxt( objLin ) - rhs;
        
        
        
        r = _checkByNonObjCut( master, zu, n, indMaster[i], pointValue );
        
        if( r != 0 )
        {
            #if MRQ_DEBUG_MODE
                MRQ_PRINTERRORNUMBER(r);
            #endif
            
            code = r;
        }
        
        
        /*if( indMaster[i] < 0 )
        {
            std::cout << "removing obj linearization of point " << i << std::endl;
            getchar();
        } */
        
    }
    
    
    
    return code;
}







class MRQ_LineSearchFunction : public NMC_Function
{
    
public:
    
    unsigned int thnumber;
    const bool* cEval;
    const double *startSol;
    const double *endSol;
    double *auxSol;
    double *auxConstr;
    MRQ_MINLPProb *prob;
    
    
    
    MRQ_LineSearchFunction( MRQ_MINLPProb *prob, unsigned int thnumber, const bool* cEval, const double *startSol, const double *endSol, double *auxSol, double *auxConstr );
    
    
    virtual int eval(double lambda, double &value) override;
    
};



MRQ_LineSearchFunction::MRQ_LineSearchFunction( MRQ_MINLPProb *prob, unsigned int thnumber, const bool* cEval, const double *startSol, const double *endSol, double *auxSol, double *auxConstr ) : NMC_Function()
{
    this->prob = prob;
    this->thnumber = thnumber;
    
    this->startSol = startSol;
    this->endSol = endSol;
    
    this->cEval = cEval;
    
    this->auxSol = auxSol;
    this->auxConstr = auxConstr;
    
}



static inline void MRQ_calcSolByLambda( const int n, const double lambda, const double *startSol, const double *endSol, double *sol )
{
    #pragma ivdep
    #pragma GCC ivdep
    for( int i = 0; i < n; i++ )
        sol[i] = lambda * startSol[i] + (1.0 - lambda)*endSol[i];
}


//that evaluation may not work so well if we have equality constraints being considered in ceval
int MRQ_LineSearchFunction::eval(double lambda, double &value)
{
    const int n = prob->n;
    const int m = prob->m;
    
    const double *lc = prob->lc;
    const double *uc = prob->uc;
    
    int ret;
    
    
    MRQ_calcSolByLambda(n, lambda, startSol, endSol, auxSol);
    
    
    ret = prob->constraintsEval(thnumber, true, cEval, auxSol, auxConstr);
    
    if(ret != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTCALLBACKERRORNUMBER(ret);
        #endif
        
        return ret;
    }
    
    
    value = -INFINITY; //note, if solution is interior, we need put a negative number in value to our method to find a zero of function work
    
    for(int i = 0; i < m; i++)
    {
        if(cEval[i] == false)
            continue;
        
        
        if( lc[i] > -MIP_INFINITY )
        {
            const double viol = lc[i] - auxConstr[i];
            
            if( viol > value )
                value = viol;
        }
        
        if( uc[i] < MIP_INFINITY )
        {
            const double viol = auxConstr[i] - uc[i];
            
            if( viol > value )
                value = viol;
        }
    }
    
    
    return 0;
}




MRQ_RETURN_CODE muriqui::MRQ_lineSearch2( const int thnumber,  MRQ_MINLPProb &prob, const double epsToLambda, const double *startSol, const double *endSol, const double absFeasTol, const double relFeasTol, const bool *cEval, double *auxConstr, double *outSol, double *outLambda )
{
    const unsigned int MAXITERS = 1000;
    const double llambda = 0.0, ulambda = 1.0;
    
    const int n = prob.n;
    
    int ret;
    unsigned int iters;
    
    
    MRQ_LineSearchFunction feval(&prob, thnumber, cEval, startSol, endSol, outSol, auxConstr);
    
    NMC_ZeroOfFunction748_2 zeroOfFunction(epsToLambda, &feval, MAXITERS);
    
    NMC_PointR2 root;
    
    
    
    /*{
        int ret = prob.constraintsEval(thnumber, true, cEval, startSol, auxConstr);
        
        assert(ret == 0);
        
        for(int i = 0; i < prob.m; i++)
        {
            std::cout << "cEval["<<i<<"]: " << cEval[i];
            
            if(cEval[i])
                std::cout << " c["<<i<<"]: " << auxConstr[i];
            
            std::cout << "\n";
        }
    } */
    
    
    ret = zeroOfFunction.findZero(llambda, ulambda, root, &iters);
    
    if( outLambda )
        *outLambda = root.x;
    
    MRQ_calcSolByLambda(n, root.x, startSol, endSol, outSol);
    
    //std::cout << "MRQ_lineSearch2: ret: " << ret << " root.x: " << root.x << " root.y: " << root.y << "\n";
    
    if(ret != 0)
    {
        #if MRQ_DEBUG_MODE
            MRQ_PRINTERRORNUMBER(ret);
        #endif
        
        if( ret == NMC_EC_EVAL_ERROR )
            return MRQ_CALLBACK_FUNCTION_ERROR;
        else if( ret == NMC_EC_MAX_ITERATIONS_STOP )
            return MRQ_MAX_ITERATIONS_STOP;
        else
            return MRQ_UNDEFINED_ERROR;
    }
    
    
    
    return MRQ_SUCCESS;
}












































