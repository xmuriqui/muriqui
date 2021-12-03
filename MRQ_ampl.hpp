


#ifndef MRQ_AMPL_HPP_
#define MRQ_AMPL_HPP_



#include "MIP_ampl.hpp"
#include "muriqui.hpp"
#include "MRQ_tools.hpp"



#define MRQ_PRINT_AMPL_EVAL_ERROR 1

#if MRQ_DEBUG_MODE
    #define MRQ_AMPL_DEBUG_MODE  1
#else
    #define MRQ_AMPL_DEBUG_MODE  0
#endif



#define MRQ_SAVE_OUTPUT_FILE 1
#define MRQ_CHAR_SEP "#"		//character separator to output file
#define MRQ_OUT_FILE_NAME "muriqui_output.txt"







namespace muriqui 
{
    
    class MRQ_AMPLParamData
    {
    public:
        
        MRQ_ALG_CODE algCode;
        MRQ_Algorithm *alg;
        MRQ_GeneralSolverParams *milpParams;
        MRQ_GeneralSolverParams *nlpParams;
        MRQ_GeneralSolverParams *minlpParams;
        MRQ_GeneralSolverParams *globalParams;
        
        
        MRQ_AMPLParamData()
        {
            algCode = MRQ_UNDEFINED_ALG;
            /*alg = NULL;
            milpParams = NULL;
            nlpParams = NULL;
            minlpParams = NULL;
            globalParams = NULL;*/
        }
    };
    
#if MRQ_HAVE_ASL
    
    //new functions....
    
    inline char * MRQ_readSolverAMPLParams(Option_Info *oi, keyword *kw, char *value, MRQ_GeneralSolverParams *solverParam)
    {
        const char *type = kw->name;
        char subname[100], strValue[100];
        long int iValue;
        int aux, ret;
        double dValue;
        
        
        
        aux = sscanf(value, "%s %s", subname, strValue);
        
        
        if( aux == 2 )
        {
            
            if(type[0] == 'i') //if( strcmp(type, "int") == 0 )
            {
                ret = sscanf(strValue, "%ld", &iValue);
                
                if( ret == 1 )
                    ret = solverParam->storeIntegerParameter(subname, iValue);
                else
                    ret = -1;
            }
            else if(type[0] == 'd') //if( strcmp(type, "dbl") == 0 )
            {
                ret = sscanf(strValue, "%lf", &dValue);
                
                if( ret == 1 )
                    ret = solverParam->storeDoubleParameter(subname, dValue);
                else
                    ret = -1;
            }
            else if(type[0] == 's') //if( strcmp(type, "str") == 0 ) //could be only else, but ok...
            {
                ret = solverParam->storeStringParameter(subname, strValue);
            }
            else
            {
                ret = -1;
                std::cerr << MRQ_PREPRINT
                "Invalid type of solver parameter specification: " << type << ". It should be [int, dbl, str]\n";
            }
            
            if(ret == 0)
            {
                std::cerr << MRQ_PREPRINT "Operation of setting subsolver parameter " << subname << " to " << strValue << " stored.\n";
            }
            else
            {
                std::cerr << MRQ_PREPRINT "Error at setting " << type << " parameter " << subname << " to " << strValue << ".\n";
            }
            
        }
        else
        {
            std::cerr << MRQ_PREPRINT "Error at reading value of parameter: " << kw->name;
            if(aux >= 1)
                std::cerr << subname;
            
            std::cerr << "\n";
        }
        
        
        
        if(aux >= 1)
        {
            //value should point to rest of string parameter
            value = value + strlen(subname); //pointer Arithmetic
            
            //avoiding white space between subname and pValue
            while( *value != '\0' && (*value == ' ' || *value == '\t') )
                value++; //pointer Arithmetic
        }
        
        
        if(aux >= 2)
        {
            value = value + strlen(strValue); //pointer Arithmetic
        }
        
        
        return value;
    }
    
    
    
    inline char* MRQ_readAlgMuriquiAMPLParams(Option_Info *oi, keyword *kw, char *value, MRQ_ALG_CODE &algCode)
    {
        char strValue[100];
        int aux, ret;
        
        aux = sscanf(value, "%s", strValue);
        
        if(aux ==1)
        {
            //algorithm parameter should be str. Anyway, even user pass wrong type, we set anyway... :)
                
            ret = MRQ_strToEnum( strValue, algCode );
            
            if(ret == 0)
                std::cerr << MRQ_PREPRINT "Muriqui algorithm set to " << strValue << "\n";
            else
                std::cerr << MRQ_PREPRINT "Invalid algorithm type to muriqui: " << strValue << "\n";
        }
        else
        {
            std::cerr << MRQ_PREPRINT "Error at reading value of algorithm parameter: " << kw->name << "\n";
        }
        
        
        
        if(aux >= 1)
        {
            value = value + strlen(strValue); //pointer Arithmetic
        }
        
        
        return value;
    }
    
    
    inline char* MRQ_readMuriquiAMPLParams(Option_Info *oi, keyword *kw, char *value, MRQ_Algorithm *alg)
    {
        const char *type = kw->name;
        char subname[100], strValue[100];
        long int iValue;
        int aux, ret;
        double dValue;
        
        
        
        
        aux = sscanf(value, "%s %s", subname, strValue);
        
        if(aux == 2)
        {
            if(type[0] == 'i') //if( strcmp(type, "int") == 0 )
            {
                ret = sscanf(strValue, "%ld", &iValue);
                
                if(ret == 1)
                    ret = alg->setIntegerParameter(subname, iValue);
                else
                    ret = -1;
            }
            else if(type[0] == 'd') //if( strcmp(type, "dbl") == 0 )
            {
                ret = sscanf(strValue, "%lf", &dValue);
                    
                if(ret == 1)
                    ret = alg->setDoubleParameter(subname, dValue);
                else
                    ret = -1;
                
            }
            else if(type[0] == 's') //if( strcmp(type, "str") == 0 ) //could be only else, but ok...
            {
                ret = alg->setStringParameter(subname, strValue);
            }
            else
            {
                ret = -1;
                std::cerr << MRQ_PREPRINT "Invalid type of solver parameter specification: " << type << ". It should be [int, dbl, str]\n" ;
            }
            
            if(ret == 0)
            {
                std::cerr << MRQ_PREPRINT << type << " Muriqui parameter " << subname << " set to " << strValue << "\n";
            }
            else //if(ret < 0) //if ret == 1, we should not show error messge
            {
                std::cerr << MRQ_PREPRINT << "Error at setting " << type << " Muriqui parameter " << subname << " to " << strValue << ".\n";
            }
        }
        else
        {
            std::cerr << MRQ_PREPRINT "Error at reading value of parameter: " << kw->name;
            if(aux >= 1)
                std::cerr <<  " " << subname;
            
            std::cerr << "\n";
            //getchar();
        }
        
        
        
        
        if(aux >= 1)
        {
            //value should point to rest of string parameter
            value = value + strlen(subname); //pointer Arithmetic
            
            //avoiding white space between subname and pValue
            while( *value != '\0' && (*value == ' ' || *value == '\t') )
                value++; //pointer Arithmetic
        }
        
        
        if(aux >= 2)
        {
            value = value + strlen(strValue); //pointer Arithmetic
        }
        
        
        return value;
    }
    
    
    inline char * MRQ_readAMPLParams(Option_Info *oi, keyword *kw, char *value)
    {
        MRQ_AMPLParamData *data = (MRQ_AMPLParamData *) kw->info;
        char firstEnv = oi->opname[8];
        
        //do not echo parameter name and value...
        oi->option_echo &= ~ASL_OI_echo;
        
        
        //check if the end of the string is "options" 
        if( firstEnv == 'o' ) //( strcmp(&(oi->opname[8]), "options" ) ) 
        {//muriqui_options
            
            return MRQ_readMuriquiAMPLParams(oi, kw, value, data->alg);
        }
        else if(firstEnv == 'n') //( strcmp(&(oi->opname[8])), "nlp_options" )
        {//muriqui_nlp_options
            return MRQ_readSolverAMPLParams(oi, kw, value, data->nlpParams);
        }
        else if(firstEnv == 'm')
        {
            if( oi->opname[10] == 'l' )
                return MRQ_readSolverAMPLParams(oi, kw, value, data->milpParams);
            else
            {
                #if MRQ_DEBUG_MODE
                if(oi->opname[10] == 'n')
                #endif
                    return MRQ_readSolverAMPLParams(oi, kw, value, data->minlpParams);
            }
        }
        else if(firstEnv == 'g')
        {
            return MRQ_readSolverAMPLParams(oi, kw, value, data->globalParams);
        }
        else if(firstEnv == 'a')
        {
            return MRQ_readAlgMuriquiAMPLParams(oi, kw, value, data->algCode);
        }
        
        std::cerr << MRQ_PREPRINT << "Error! Unknow solver environment at MRQ_readAMPLParam!\n";
        MRQ_getchar();
        
        return NULL;
    }
    
    

#endif




    



    class MRQ_ReadAmplModel : public minlpproblem::MIP_ReadAmplModel
    {
        
    public:
        
        
        void readAlgChoice(MRQ_AMPLParamData *paramData);
        
        void readMuriquiParameters(MRQ_AMPLParamData *paramData);
        
    };



    int MRQ_ampl(char* stub, const bool printAlgParameters = false, const bool printProblem = false);
    

}


#endif
