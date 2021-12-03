/*
 * Muriqui MINLP solver
 * 
 * By Wendel Melo: autonomous researcher
 * 
 * usage:
 * 
 * 
 * muriqui <input nl file> [options]
 * 
 * 
 * options:
 *  -p : print algoritmh parameter values.
 * 
 * algorithm choice can be made by means of the file muriqui_algorithm.opt in the current directory
 * 
 * algorithm and subsolver parameters can be set by means of the file file muriqui_params.opt
 * 
 * 
 * 
 */

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>

#include <iostream>

#include "muriqui.hpp"
#include "MRQ_bb.hpp"
#include "MRQ_ampl.hpp"
#include "MRQ_tools.hpp"

#include "MIP_gams.hpp"


using namespace muriqui;



int MRQ_gams(char *stub, const bool printAlgParameters, const bool printProblem);



int main(int argc, char **argv)
{
    bool amplModel = false, gamsModel = false; //by default, we will assume we have a ampl model
    
    
    
    MRQ_helloOnceTime();
    
    if( argc > 1 )
    {
        bool printParam = false;
        bool printProblem = false;
        char *inputFile = NULL;
        
        
        for(int i = 1; i < argc; i++)
        {
            if( strcmp("-p", argv[i]) == 0 )
                printParam = true;
            
            if( strcmp("-f", argv[i]) == 0 )
                printProblem = true;
            
            else if( strcmp("-gams", argv[i]) == 0 ) //the flag -gams is passed by our own file gmsmq_ux.out
                gamsModel = true;
            else if( strcmp("-AMPL", argv[i]) == 0 )
                amplModel = true;
            
            if( argv[i][0] != '-' && inputFile == NULL)
                inputFile = argv[i];
            
            std::cout << MRQ_PREPRINT "parameter " << i << ": " << argv[i] << "\n";
        }
        
        if( inputFile )
        {
            std::cout << MRQ_PREPRINT "input file: " << inputFile << "\n";
            
            if( gamsModel )
            {
                #if MRQ_HAVE_GAMS
                    MRQ_gams(inputFile, printParam, printProblem);
                #else
                    MRQ_PRINTERRORMSG("Error! GAMS Joat Library was not available in the compilation time!\n");
                #endif
            }
            else //by default, we assume an AMPL model
            {
                #if MRQ_HAVE_ASL
                    MRQ_ampl(inputFile, printParam, printProblem);
                #else
                    MRQ_PRINTERRORMSG("Error! AMPL Solver Library was not available in the compilation time!\n");
                #endif
            }
        }
    }
    else
    {
        std::cout << 
        "usage:\n\n"
        
        "\t\t" << argv[0] << " <input nl/gams file> [options]\n\n"
        
        "options:\n"
        "\t-p : print algorithm parameter values.\n"
        "\t-AMPL : input file cames from AMPL environment, i.e., is a nl file (default).\n"
        "\t-gams : input file cames from GAMS environment.\n"
        
        "\nAlgorithm choice can be made by means of the file \"" MRQ_MURIQUI_ALG_CHOICE_FILE "\" in the current directory.\n"
        
        "\nAlgorithm parameters can be set by means of the file \"" MRQ_MURIQUI_PARAMS_FILE "\" in the current directory.\n"
        
        "\nMilp solver parameters can be set by means of the file \"" MRQ_MILP_SOLVER_PARAMS_FILE "\" in the current directory\n"
        
        "\nNlp solver parameters can be set by means of the file \"" MRQ_NLP_SOLVER_PARAMS_FILE "\" in the current directory\n"
        
        "\n"
        ;
    }
    
    return 0;
}


