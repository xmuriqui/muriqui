/* That file contains a implementation of
* a version of diving heuristic for convex MINLP problems.
* It is a heuristic. Its compromise is give us a feasible solution...
*
* References:
*
* 	Bonami & Goncalves, Heuristics for convex mixed integer nonlinear programs.
* Comput Optim Appl, 2010 (online).
*
* Author: Wendel Melo
*
* Date: 19-Feb-2014
*
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <new>

#include "MRQ_solvers.hpp"
#include "MRQ_algClasses.hpp"
#include "MRQ_tools.hpp"


using namespace optsolvers;
using namespace muriqui;

//using namespace std;


inline double MRQ_intGap(const double value)
{
    return MRQ_abs( value - round(value) );
}



MRQ_Diving::MRQ_Diving():MRQ_Heuristic()
{
    resetParameters();
    resetOutput();
    out_algorithm = MRQ_DIVE_HEUR_ALG;
}


MRQ_Diving::~MRQ_Diving(){}



void MRQ_Diving::printParameters(std::ostream &out) const
{
    char strValue[100];
    MRQ_Heuristic::printParameters(out);
    out << "\n"
    
    MRQ_STRFFATT(in_consider_relax_infeas_if_solver_fail) << "\n";
    
    MRQ_enumToStr(in_dive_selec_strategy, strValue);
    out << MRQ_STR(in_dive_selec_strategy) " " << strValue << "\n"
    
    MRQ_STRFFATT(in_percentual_of_add_var_fixing) << "\n";
}



void MRQ_Diving::resetParameters()
{
    MRQ_Heuristic::resetParameters();
    
    in_consider_relax_infeas_if_solver_fail = true;
    in_dive_selec_strategy = MRQ_DIVE_SS_FRACTIONAL;
    in_percentual_of_add_var_fixing = 0.2;
}



int MRQ_Diving::setDoubleParameter(const char *name, const double value)
{
    int ret = MRQ_Heuristic::setDoubleParameter(name, value);
    
    if(ret == 0)
        return 0;
    
    ret = 0;
    
    if( MRQ_setAtt( MRQ_STRATT(in_percentual_of_add_var_fixing), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


int MRQ_Diving::setIntegerParameter(const char *name, const long int value)
{
    int ret = MRQ_Heuristic::setIntegerParameter(name, value);
    
    if(ret == 0)
        return 0;
    
    ret = 0;
    
    
    if( MRQ_setAtt<bool>( MRQ_STRATT(in_consider_relax_infeas_if_solver_fail), name, value ) == 0 );
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


int MRQ_Diving::setStringParameter(const char *name, const char *value)
{
    int r, ret = MRQ_Heuristic::setStringParameter(name, value);
    
    if(ret == 0)
        return 0;
    
    
    
    if( (r = MRQ_setStrAtt( MRQ_STRATT(in_dive_selec_strategy), name, value ) ) >= 0 )
    {
        ret = r == 0 ? 0 : MRQ_VALUE_ERROR;
    }
    else
        ret = MRQ_NAME_ERROR;
    
    
    return ret;
}


MRQ_RETURN_CODE MRQ_Diving::selectFracVariable(MRQ_MINLPProb &prob, const int *intVars, const double *sol, const int *nrowsVars, double *auxVars, bool &intSol, unsigned int &index, bool &roundUp)
{
    //const int n = prob.n;
    const unsigned int nI = prob.getNumberOfIntegerVars();
    unsigned int i, j;
    bool gradCalc = false;
    MRQ_RETURN_CODE code;
    double aux, best = INFINITY;
    const double epsilon = 1.0e-6;
    
    
    if( in_dive_selec_strategy == MRQ_DIVE_SS_FRACTIONAL )
    {
        
        for(i = 0; i < nI; i++)
        {
            aux = MRQ_intGap( sol[ intVars[i] ] ); 
            if( aux > in_integer_tol && aux < best )
            {
                best = aux;
                index = intVars[i];
            }
        }
        
        if(best > 1.0)
        {
            intSol = true;
        }
        else
        {
            intSol = false;
            roundUp = sol[index] - floor(sol[index]) >= 0.5;
        }
        
    }
    else
    {
        
        for(i = 0; i < nI; i++)
        {
            j = intVars[i];
            
            if( MRQ_intGap( sol[j] ) > in_integer_tol )
            {
                
                if( !gradCalc )
                {
                    int r = prob.objGradEval(thnumber, true, sol, auxVars);
        
                    if(r != 0)
                    {
                        code = MRQ_CALLBACK_FUNCTION_ERROR;
                        goto termination;
                    }
                    
                    gradCalc = true;
                }
                
                
                if( auxVars[j] >= 0.0 )
                    aux = ceil(sol[j]);
                else
                    aux = floor(sol[j]);
                
                aux = ( (aux - sol[j])*auxVars[j] + epsilon )/( nrowsVars[j] ); //we already sum 1.0 in nrowsVars on its calculation
                
                if( aux < best )
                {
                    best = aux;
                    index = j;
                }
            }
        }
        
        if( std::isinf(best) )
        {
            intSol = true;
        }
        else
        {
            intSol = false;
            roundUp = auxVars[index] >= 0.0;
        }
        
    }
    
    
    code = MRQ_SUCCESS;
    
termination:
    
    return code;
}




int MRQ_Diving::run(MRQ_MINLPProb &prob, MRQ_GeneralSolverParams* milpSolverParams, MRQ_GeneralSolverParams* nlpSolverParams)
{
    clock_t clockStart = clock();
    double timeStart = MRQ_getTime();
    
    const int n = prob.n;
    const int m = prob.m;
    const int nI= prob.getNumberOfIntegerVars();
    const int K = ceil( nI * in_percentual_of_add_var_fixing );
    const bool preprocess = in_preprocess_lin_constr || in_preprocess_quad_constrs || in_preprocess_quad_constrs;
    
    bool isIntSol, roundUp, updtConstrBounds;
    unsigned int j, jfix;
    int tamIFix;
    int ret;
    unsigned long int iter = 0;
    double rounded, ljfix, ujfix;
    
    
    MRQ_Preprocessor preprocessor(&prob);
    MRQ_NLPSolver *nlp = NULL;
    
    int *intVars = NULL, *revIntVars = NULL, *IFix = NULL;
    int *nrowsVars = NULL; //number of rows where each variable appears
    double *auxVars = NULL;
    double *lintv = NULL, *uintv;
    double *lintb = NULL, *uintb;
    double *lx = run_by_inside ? nlx : prob.lx;
    double *ux = run_by_inside ? nux : prob.ux;
    double *flx = NULL, *fux;
    double *plc = NULL, *puc = NULL;
    
    
    {
        auto ret = algorithmInitialization(1, (in_preprocess_lin_constr || in_preprocess_obj_function || in_preprocess_quad_constrs), milpSolverParams, nlpSolverParams, prob, lx, ux, &preprocessor, &updtConstrBounds, &plc, &puc); //that algorithm is monothread...
        
        if(ret != 0)
        {
            if(in_print_level > 0 )
                MRQ_PRINTERRORNUMBER(ret);
            
            out_return_code = ret;
            goto termination;
        }
    }
    
    
    if(in_print_level > 1)
        std::cout << "\n" MRQ_PREPRINT "Starting Diving heuristic\n\n";
    
    if(in_print_level > 3)
        printSubSolvers(false, true, false);
    
    auxVars = (double *) malloc( n * sizeof(double) );
    revIntVars = (int *) malloc( n * sizeof(int) );
    IFix = (int *) malloc( K * sizeof(int) );
    intVars = (int *) malloc( nI * sizeof(int) );
    lintv = (double *) malloc( 2*nI*sizeof(double) );
    lintb = (double *) malloc( 2*K*sizeof(double) );
    nlp = OPT_newNLPSolver( in_nlp_solver );
    if( !auxVars || !revIntVars || !IFix || !intVars || !lintv || !lintb || !nlp )
    {
        if(in_print_level > 0)
            MRQ_PRINTMEMERROR;
        out_return_code = MRQ_MEMORY_ERROR;
        goto termination;
    }
    
    uintv = &lintv[nI];
    uintb = &lintb[K];
    
    if( preprocess )
    {
        flx = (double *) malloc( 2*n * sizeof(double) );
        if( !flx )
        {
            if(in_print_level > 0)
                MRQ_PRINTMEMERROR;
            out_return_code = MRQ_MEMORY_ERROR;
            goto termination;
        }
        fux = &flx[n];
        
        MRQ_copyArray(n, lx, flx);
        MRQ_copyArray(n, ux, fux);
    }
    
    
    prob.getIntegerIndices(intVars);
    prob.getReverseIntegerIndices(revIntVars);
    for(int i = 0; i < nI; i++)
    {
        const int ind = intVars[i];
        lintv[i] = lx[ind];
        uintv[i] = ux[ind];
    }
    
    
    
    if( in_dive_selec_strategy == MRQ_DIVE_SS_VECTORLENGHT )
    {
        nrowsVars = (int *) malloc( n * sizeof(int) );
        
        if(!nrowsVars)
        {
            if(in_print_level > 0)
                std::cerr << MRQ_PREPRINT << "Memory error " << MRQ_GETFILELINE << std::endl;
        
            out_return_code = MRQ_MEMORY_ERROR;
            goto termination;
        }
        
        MRQ_setAllArray(n, nrowsVars, 1); //we have to sum 1.0 in the denominator to calc the radio. So, we sum already here
        //for(int i = 0; i < n; i++)
            //nrowsVars[i] = 1; //we have to sum 1.0 in the denominator to calc the radio. So, we sum already here
        
        prob.A.countRowsEachColumn(nrowsVars, true);
        prob.J.countRowsEachColumn(nrowsVars, true);
    }
    
    
    ret = MRQ_setNLPRelaxProb(prob, lx, ux, plc, puc, nlp, true, true, true, false, thnumber, in_set_special_nlp_solver_params, nlpSolverParams, in_number_of_threads, in_max_cpu_time, in_max_time, 0, 0);
    if( ret != 0 )
    {
        if(in_print_level > 0)
            MRQ_PRINTERRORNUMBER(ret);
        
        out_return_code = MRQ_NLP_SOLVER_ERROR;
        goto termination;
    }
    
    
    if( in_use_initial_solution )
    {
        nlp->setInitialSolution(xInit, NULL, NULL);
        prob.objEval(thnumber, true, xInit, nlp->objValue); //I am sorry for set nlp->objValue...
        
        MRQ_copyArray(n, xInit, nlp->sol); //I am sorry for set nlp->sol...
    }
    else
    {
        
        ret = nlp->solve(false);
        
        //std::cout << "first nlp. ret: " << ret << " orig ret code: " << nlp->origSolverRetCode << "\n";
        
        if(ret != OPT_OPTIMAL_SOLUTION  &&  !nlp->feasSol)
        {
            if(in_consider_relax_infeas_if_solver_fail)
                out_return_code = MRQ_INFEASIBLE_PROBLEM;
            else
                out_return_code = MRQ_NLP_SOLVER_ERROR;
            
            if(in_print_level > 1)
                std::cerr << MRQ_PREPRINT  "Fail to get initial solution at MRQ_Diving::run!\n";
            
            goto termination;
        }
        
        if( ret == OPT_OPTIMAL_SOLUTION )
            zl = nlp->objValue;
        
        
        nlp->setInitialSolution( nlp->sol, nlp->dualSolC, nlp->dualSolV );
    }
    
    
    
    
    while( true )
    {
        const double *sol = nlp->sol;
        
        iter++;
        
        
        {
            auto ret = selectFracVariable(prob, intVars, sol, nrowsVars, auxVars, isIntSol, jfix, roundUp);
        
            if(ret != 0)
            {
                if( in_print_level > 0 )
                    MRQ_PRINTERRORNUMBER(ret);
                
                out_return_code = ret;
                break;
            }
        }
        
        
        
        if(isIntSol)
        {//integer solution found
            
            if( nlp->objValue < zu )
            {
                zu = out_best_obj = nlp->objValue;
                MRQ_copyArray(n, sol, out_best_sol);
                
                out_return_code = MRQ_HEURISTIC_SUCCESS;
            }
            else
                out_return_code = MRQ_HEURISTIC_FAIL;
            
            break;
        }
        
        
        ljfix = lintv[ revIntVars[jfix] ];
        ujfix = uintv[ revIntVars[jfix] ];
        
        if(roundUp)
        {
            rounded = ceil( sol[jfix] );
            lintv[ revIntVars[jfix] ] = rounded;
            
            //printf("if rounded: %f lintv[ revIntVars[jfix] ]: %f\n", rounded, lintv[ revIntVars[jfix] ]);
        }
        else
        {
            rounded = floor( sol[jfix] );
            uintv[ revIntVars[jfix] ] = rounded;
            
            //printf("else rounded: %f uintv[ revIntVars[jfix] ]: %f\n", rounded, uintv[ revIntVars[jfix] ]);
        }
        
        
        if( preprocess )
        {
            //MRQ_copyArray( n, lx, flx );
            //MRQ_copyArray( n, ux, fux );
            
            //flx[jfix] = lintv[ revIntVars[jfix] ];
            //fux[jfix] = uintv[ revIntVars[jfix] ];
        }
        else
        {
            ret = nlp->setVariableBounds(jfix, lintv[ revIntVars[jfix] ], uintv[ revIntVars[jfix] ] );
            
            #if MRQ_DEBUG_MODE
                assert(ret == 0);
            #endif
        }
        
        
        //fixing possible integer variables (at most K)
        tamIFix = 0;
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            if(lintv[i] != uintv[i])
            {
                if( sol[ind] - lintv[i] <= in_integer_tol ) //var ind is on lower bound
                {
                    lintb[ tamIFix ] = lintv[i]; //backing up lintv[i]
                    uintb[ tamIFix ] = uintv[i]; //backing up uintv[i]
                    
                    uintv[i] = lintv[i];
                    
                    if( preprocess )
                    {
                        //flx[ind] = lintv[i];
                        //fux[ind] = lintv[i];
                    }
                    else
                        nlp->setVariableBounds(ind, lintv[i], lintv[i]);
                    
                    
                    IFix[ tamIFix ] = i;
                    tamIFix++;
                    
                    if( tamIFix >= K)
                        break;
                }
                else if( uintv[i] - sol[ind] <= in_integer_tol ) //var ind is on upper bound
                {
                    lintb[ tamIFix ] = lintv[i]; //backing up lintv[i]
                    uintb[ tamIFix ] = uintv[i]; //backing up uintv[i]
                    
                    lintv[i] = uintv[i];
                    
                    if( preprocess )
                    {
                        //flx[ind] = uintv[i];
                        //fux[ind] = uintv[i];
                    }
                    else
                        nlp->setVariableBounds(ind, uintv[i], uintv[i]);
                    
                    
                    IFix[ tamIFix ] = i;
                    tamIFix++;
                    
                    if( (int) tamIFix >= K)
                        break;
                }
            }
        }
        
        
        if( in_print_level > 3 && tamIFix > 0 )
            std::cout << "Fixing " << tamIFix << " variables at their actual integer values" << std::endl;
        
        
        for(int i = 0; i < 3; i++)
        {
            
            if( preprocess )
            {
                bool updtvb, updtcb;
                
                
                //we do not need copy more bounds of integers vars because we do it before this loop and when change some value of fix vars...
                //MRQ_copyArray(n, lx, flx);
                //MRQ_copyArray(n, ux, fux);
                
                #pragma ivdep
                #pragma GCC ivdep
                for(int j = 0; j < nI; j++)
                {
                    const int ind = intVars[j];
                    flx[ind] = lintv[j];
                    fux[ind] = uintv[j];
                }
                
                ret = preprocessor.preprocess( in_preprocess_quad_constrs, in_preprocess_obj_function && zu < MRQ_INFINITY, zu, flx, fux, updtvb, updtcb, NULL, NULL, plc, puc );
                
                if( ret == minlpproblem::MIP_INFEASIBILITY )
                {
                    ret = OPT_INFEASIBLE_PROBLEM;
                    nlp->retCode = OPT_INFEASIBLE_PROBLEM;
                    nlp->feasSol = false;
                }
                else
                {
                    ret = nlp->setnVariablesBounds(n, flx, fux);
                    for(int i = 0; i < m; i++)
                        ret += nlp->setConstraintBounds(i, plc[i], puc[i]);
                    
                    #if MRQ_DEBUG_MODE
                        if(ret != 0)
                        {
                            if(in_print_level > 0)
                                MRQ_PRINTERROR;
                            out_return_code = MRQ_NLP_SOLVER_ERROR;
                            goto termination;
                        }
                    #endif
                }
                
            }
            else
            {
                ret = 0;
            }
            
            if( ret != OPT_INFEASIBLE_PROBLEM )
            {
                ret = nlp->solve(false);
                
                /*printf("nlp->solve ret: %d obj: %f iter: %d\n", ret, nlp->objValue, (int) iter);
                double *sol = nlp->sol;
                for(int w = 0; w < n; w++)
                    printf("flx[%d]: %f fux[%d]: %f sol[%d]: %f\n", w, flx[w], w, fux[w], w, sol[w]);
                MRQ_getchar();*/
            }
            
            
            if( ret == OPT_OPTIMAL_SOLUTION )
            {
                if( nlp->objValue >= zu ) //our solution cannot be better than zu
                {
                    out_return_code = MRQ_HEURISTIC_FAIL;
                    goto termination;
                }
                
                break;
            }
            else if(nlp->feasSol)
            {
                break;
            }
            else
            {
                if(i == 0 && tamIFix <= 0)
                    i = 1;//we do not have variables to unfix. So, we jump to reverse j branching
                
                if(i == 0 )
                {
                    //unfixing fixed variables...
                    for(int j = 0; j < tamIFix; j++)
                    {
                        const int ind = intVars[ IFix[j] ];
                        
                        lintv[ IFix[j] ] = lintb[j]; //lx[p];
                        uintv[ IFix[j] ] = uintb[j]; //ux[p];
                        
                        if( preprocess )
                        {
                            //flx[ind] = lintb[j];
                            //fux[ind] = uintb[j];
                        }
                        else
                            nlp->setVariableBounds(ind, lintb[j], uintb[j]);
                    }
                }
                else if(i == 1)
                {
                    j = revIntVars[jfix];
                    
                    if( roundUp )
                    {
                        if( preprocess )
                        {
                            //flx[jfix] = ljfix;
                            //fux[jfix] = rounded - 1;
                        }
                        else
                            nlp->setVariableBounds(jfix, ljfix, rounded - 1);
                        
                        lintv[j] = ljfix;
                        uintv[j] = rounded - 1;
                    }
                    else
                    {
                        if( preprocess )
                        {
                            //flx[jfix] = rounded + 1;
                            //fux[jfix] = ujfix;
                        }
                        else
                            nlp->setVariableBounds(jfix, rounded + 1, ujfix);
                        
                        lintv[j] = rounded + 1;
                        uintv[j] = ujfix;
                    }
                }
                else
                {
                    out_return_code = MRQ_HEURISTIC_FAIL;
                    goto termination;
                }
                
                
                //if we are preprocessing, we have to restore original bounds of continuous values because they can have been changed by preprocess due to fix integer variables. Since we change this fixing, we restore the bounds...
                
                //Note, we should restore only bounds of continuous vars, but we retsore all bounds. However there is no problem since before we solve, we update bounds of integers vars
                if( preprocess )
                {
                    MRQ_copyArray(n, lx, flx);
                    MRQ_copyArray(n, ux, fux);
                }
                
            }
            
        }
        
        
        if( checkTerminationCriterions(thnumber, zl, zu, iter, timeStart, clockStart, out_return_code) )
        {
            goto termination;
        }
        
        //getchar();
    }
    
    
termination:
    
    if(in_print_level > 1)
    {
        if( out_return_code == MRQ_HEURISTIC_SUCCESS )
            std::cout << "Diving heuristic found a feasible solution! Obj Function: " << out_best_obj << " ";
        else
            std::cout << "Diving heuristic did not find a feasible solution! ";
    }
    
    
    if(plc)			free(plc);
    if(flx)			free(flx);
    
    if(nlp)			delete nlp;
    if(revIntVars)	free(revIntVars);
    if(intVars) 	free(intVars);
    if(IFix)		free(IFix);
    if(nrowsVars)	free(nrowsVars);
    if(auxVars)		free(auxVars);
    if(lintv)		free(lintv);
    if(lintb)		free(lintb);
    
    
    out_feasible_solution = out_best_obj < MRQ_INFINITY;
    out_number_of_iterations = iter;
    //out_algorithm = muriqui::MRQ_FP_HEUR_ALG;
    out_lower_bound = zl;
    out_upper_bound = zu;
    
    algorithmFinalization(1, prob, lx, ux);
    
    out_cpu_time = ( (double) (clock() - clockStart) )/CLOCKS_PER_SEC;
    out_clock_time = MRQ_getTime() - timeStart;
    
    if(in_print_level > 1)
        std::cout << "cpu time: " << out_cpu_time << std::endl;
    
    
    
    return out_return_code;
}

