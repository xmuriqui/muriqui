
#include <cstring>
#include <iostream>


#include "muriqui.hpp"
#include "DCT_bbclient.hpp"
#include "MRQ_tools.hpp"


using namespace muriqui;
using namespace dctools;


#define MRQ_COMPONENT_FILE_FLAG "-c"
#define MRQ_ALGORITHM_FLAG "-a"
#define MRQ_AMPL_MODEL_FILE_FLAG "-AMPL"
#define MRQ_GAMS_MODEL_FILE_FLAG "-gams"



int main(int argc, char **argv)
{
	const unsigned int defaultPort = 22022;
	const unsigned int defaultMaxThreads = 0;
	
	
	int r;
	int retCode;
	char *componentFileName = NULL;
	
	DCT_UInt64 nBytesBasicParameters;
	DCT_Int32 basicParameters[10];
	
	
	DCT_Components components;
	DCT_SeveralNumberOfComponents allNofComponents;
	DCT_ComponentsFileReader reader(defaultPort, defaultMaxThreads);
	
	DCT_GeneralParams algorithmParams, milpParams, nlpParams, globalsolverParams;
	
	bool lastParamIsInputFile = true;
	int algOption = MRQ_LP_BB_ECP_BASED_ALG;
	DCT_FileNames inputFiles;
	//DCT_VarBounds varsBounds;
	DCT_AllGeneralParams allGeneralParams;
	
	
	//MRQ_helloOnceTime();
	
	if( argc < 4 )
	{
		std::cout << "usage:\n\n"
		
		<< argv[0] << " [-AMPL | -gams] <input problem file> " MRQ_COMPONENT_FILE_FLAG " <input component file>" " [" MRQ_ALGORITHM_FLAG " <algorithm>]"  "\n\n"
		
		"\talgorithm option (with flag" MRQ_ALGORITHM_FLAG ") can be:\n"
		"\t\tlpbb: (default) lp based branch-and-bound (gecpbb)\n"
		"\t\tgecpbb: the same that lpbb\n"
		"\t\tlpnlpbb: lp\\nlp based based branch-and-bound\n"
		"\t\toabb: the same that lpnlpbb\n"
		"\t\teshbb: extended supportig hiperplane based branch-and-bound\n"
		"\t\tnlpbb: standard nonlinear branch-and-bound\n"
		"\t\tbb: the same that nlpbb\n"
		;
		
		retCode = MRQ_VALUE_ERROR;
		goto termination;
	}
	
	
	for(int i = 1; i < argc; i++)
	{
		if(argv[i][0] == '-' )
		{
			if(i == argc-1)
			{
				std::cout << "missing value for option " << argv[i] << "\n";
				retCode = MRQ_VALUE_ERROR;
				goto termination;
			}
			
			if( strcmp(MRQ_COMPONENT_FILE_FLAG, argv[i]) == 0)
			{
				componentFileName = argv[i+1];
				//we increment i to skip tis parameter 
				
				if( i == argc-2 )
					lastParamIsInputFile = false;
				
				i++;
			}
			else if( (strcmp(MRQ_AMPL_MODEL_FILE_FLAG, argv[i]) == 0) || (strcmp(MRQ_GAMS_MODEL_FILE_FLAG, argv[i]) == 0) ) 
			{
				MRQ_INPUT_FILE_TYPE inFileType = strcmp(MRQ_GAMS_MODEL_FILE_FLAG, argv[i]) == 0 ? MRQ_IFT_GAMS_MODEL_FILE : MRQ_IFT_AMPL_MODEL_FILE;
				
				if( inputFiles.count(inFileType) == 0 )
				{
					try{
						inputFiles.insert( {inFileType, argv[i + 1]} );
					}
					catch( std::bad_alloc &ba )
					{
						retCode = MRQ_MEMORY_ERROR;
						goto termination;
					}
				}
			}
			else if( strcmp(MRQ_ALGORITHM_FLAG, argv[i]) == 0 )
			{
				char * strAlgOption = argv[i+1];
				
				if( strcmp(strAlgOption, "lpbb") == 0 || strcmp(strAlgOption, "gecpbb") == 0 )
				{
					algOption = MRQ_LP_BB_ECP_BASED_ALG;
				}
				else if( strcmp(strAlgOption, "lpnlpbb") == 0 || strcmp(strAlgOption, "oabb") == 0 )
				{
					algOption = MRQ_LP_NLP_BB_OA_BASED_ALG;
				}
				else if( strcmp(strAlgOption, "eshbb") == 0 )
				{
					algOption = MRQ_LP_BB_ESH_BASED_ALG;
				}
				else if( strcmp(strAlgOption, "nlpbb") == 0 || strcmp(strAlgOption, "bb") == 0 )
				{
					algOption = MRQ_BB_ALG;
				}
				else
				{
					std::cout << "invalid value: " << strAlgOption << "for option " << argv[i] << "\n";
					retCode = MRQ_VALUE_ERROR;
					goto termination;
				}
				i++;
			}
			else
			{
				std::cout << "invalid option: " << argv[i] << "\n";
				retCode = MRQ_VALUE_ERROR;
				goto termination;
			}
		}
		else
		{
			//we assume this file is part of problem definition, and the format is AMPL
			try{
				inputFiles.insert( {MRQ_IFT_AMPL_MODEL_FILE, argv[i]} );
			}
			catch( std::bad_alloc &ba )
			{
				retCode = MRQ_MEMORY_ERROR;
				goto termination;
			}
			
		}
	}
	
	
	
	if(lastParamIsInputFile)
	{
		//last argument: we assume this file is part of problem definition
		try{
			inputFiles.insert( {MRQ_IFT_AMPL_MODEL_FILE, argv[argc-1]} );
		}
		catch( std::bad_alloc &ba )
		{
			retCode = MRQ_MEMORY_ERROR;
			goto termination;
		}
		
	}
	
	
	if( inputFiles.count(MRQ_IFT_AMPL_MODEL_FILE) > 0 && inputFiles.count(MRQ_IFT_GAMS_MODEL_FILE) > 0 )
	{
		DCT_PRINTERRORMSG("Warning: both model file types AMPL and GAMS were specified. Running only AMPL file.");
	}
	
	
	
	if( componentFileName == NULL )
	{
		MRQ_PRINTERRORMSG("no component file name");
		retCode = MRQ_VALUE_ERROR;
		goto termination;
	}
	
	
	
	
	
	
	r = reader.read(componentFileName, components, allNofComponents);
	MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
	
	
	
	std::cout << "componentes lidos: \n";
	components.print();
	
	std::cout << "\nconjuntos de numeros de componentes lidos: \n";
	allNofComponents.print();
	
	
	{
		//just to test our server and client we pass a list of integer parameters
		
		basicParameters[0] = algOption;  //algorithm option
		basicParameters[1] = MRQ_CPLEX;		//milp solver
		basicParameters[2] = MRQ_NLP_MOSEK;	//nlp solver
		basicParameters[3] = MRQ_NLP_KNITRO; //global solver
		basicParameters[4] = algOption == MRQ_BB_ALG ? 2000 : 500; //maximum number of nodes to be sent from a server to another
		basicParameters[5] = algOption == MRQ_BB_ALG ? 20000 : 200000; //frequency to server send lower bounds to client
		
		nBytesBasicParameters = 6 * sizeof(basicParameters[0]);
	}
	
	
	{
		//just to test our server, we set some parameters
		algorithmParams.addIntegerParameter("in_print_level", 10);
		//algorithmParams.addIntegerParameter("in_printing_frequency", 100);
		algorithmParams.addIntegerParameter("in_use_outer_app", 0);
		algorithmParams.addIntegerParameter("in_use_outer_app_as_heuristic", 0);
		
		algorithmParams.addIntegerParameter("in_print_parameters_values", 1);
		//algorithmParams.addDoubleParameter("in_max_cpu_time", 4*60*60);
		//algorithmParams.addDoubleParameter("in_max_time", 24*60*60);
		//algorithmParams.addDoubleParameter("in_relative_convergence_tol", 1.0E-2);
		//algorithmParams.addDoubleParameter("in_absolute_convergence_tol", 1.991E-4);
		
		algorithmParams.addStringParameter("in_constr_branching_strategy", "MRQ_BB_CBS_NO_CONSTRAINT_BRANCH");
		//algorithmParams.addStringParameter("in_igma2_strategy", "MRQ_BB_I2S_NO_IGMA2");
		algorithmParams.addStringParameter("in_igma2_gap_min_solver", "MRQ_NLP_IPOPT");
		
		
		milpParams.addIntegerParameter("CPX_PARAM_THREADS", 1);
		milpParams.addDoubleParameter("CPX_PARAM_TILIM", 40*60*60);
		
		
		nlpParams.addIntegerParameter("MSK_IPAR_INTPNT_MAX_ITERATIONS", 1000);
		//nlpParams.addStringParameter("derivative_test", "first-order");
		
		
		//globalsolverParams.addStringParameter("linear_solver", "ma27");
		
		
		allGeneralParams[MRQ_DC_GPT_MURIQUI] = &algorithmParams;
		allGeneralParams[MRQ_DC_GPT_MILP_SOLVER] = &milpParams;
		allGeneralParams[MRQ_DC_GPT_NLP_SOLVER] = &nlpParams;
		allGeneralParams[MRQ_DC_GPT_GLOBAL_SOLVER] = &globalsolverParams;
	}
	
	
	
	/*{
		//just to test our server, we pass some variable bounds
		DCT_Bounds b;
		
		b.lb = -0.9; b.ub = 1.9;
		varsBounds[29] = b;
		
		b.lb = -0.8; b.ub = 1.8;
		varsBounds[28] = b;
		
		b.lb = -0.7; b.ub = 1.7;
		varsBounds[27] = b;
		
		b.lb = -0.6; b.ub = 1.6;
		varsBounds[26] = b;
		
		b.lb = -0.5; b.ub = 1.5;
		varsBounds[25] = b;
	}*/
	
	
	
	
	{
		std::string strOptCode;
		DCT_BBClient bbclient;
		DCT_FinalResults &fr = bbclient.finalResults;
		
		r = bbclient.setSigIntHandler();
		if(r != 0)
		{
			MRQ_PRINTERRORNUMBER(r);
		}
		
		
		r = bbclient.run(&components, &allNofComponents, &inputFiles, nBytesBasicParameters, (DCT_Byte*) basicParameters, &allGeneralParams);
		MRQ_IFERRORGOTOLABEL(r, retCode, r, termination);
		
		
		//getting the solution
		
		r = DCT_enumToStr( (DCT_OPTIMIZATION_RETURN_CODE) fr.optCode, strOptCode);
		if(r != 0)
			MRQ_PRINTERRORNUMBER(r);
		
		std::cout << "muriqui client - algorithm: " << fr.algorithm << " code: " << fr.optCode << " (" << strOptCode << ") lb: " << fr.lowerBound << " obj: " << fr.objBestSol << " cputime: " << fr.cpuTime << " clocktime: " << fr.clockTime << " iters: " << fr.nIters << " server calls: " << fr.nServerCalls << "\n" ;
		
		
		if( fr.feasSolFound )
		{
			double *sol = (double *) bbclient.bestSol;
			const unsigned int n = bbclient.nvars;
			
			std::cout << "obj: " << bbclient.bestObjValue << "\n";
			std::cout << "solution:\n";
			for(unsigned int i = 0; i < n; i++)
				printf("\t%0.16f\n", sol[i]);//std::cout << "x["<<i<<"]: " << sol[i] << std::endl;
		}
		
	}
	
	
	
	retCode = 0;
termination:
	
	
	return retCode;
}

