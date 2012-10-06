/*
 *  options.cpp
 *  HMM
 *
 *  Created by Paul C Lott on 8/13/09.
 *  Copyright 2009 University of California, Davis. All rights reserved.
 *
 */

#include "options.h"

namespace StochHMM{


    //----------------------------------------------------------------------------//
    // Description: opt_data initialization                                                              
    //
    //
    // 
    //----------------------------------------------------------------------------//
    opt_data::opt_data(){
        type=UNDEF_OPT;           //what datatype does the option store
        set=false;           //Has the parameter been set
        required=false;      //Is the parameter required, else exit
        restricted=false;    //Are the entries restricted to allowable
    }


    //!Setup options to handle the parameters and options defined by the user
    //!\param param pointer to opt_parameters(user defined options and defaults)
    //!\param size  # of parameters that have been defined
    //!\param usageStatement User defined usage statement for the program
    void options::set_parameters(opt_parameters *param,int size,const char* usageStatement){
        
        usage=usageStatement;
        
        //Process paramaters
        for(int i=0;i<size;i++){
            //Get possible commandline tags and split possible tags;
            std::string opt_tag=param[i].commandline_tag;
            std::vector<std::string> tags;
            size_t found=opt_tag.find_first_of(":");
            if (found==std::string::npos){
                tags.push_back(opt_tag);
            }
            else{
                while (found!=std::string::npos){
                    tags.push_back(opt_tag.substr(0,found));
                    opt_tag.erase(0,found+1);
                    found=opt_tag.find_first_of(":");
                    if (found==std::string::npos){
                        tags.push_back(opt_tag);
                    }
                }
            }
            
            // Check allowable;
            int allowable_defined=0;
            for (int j=0;j<MAX_ALLOWABLE;j++){
                if (param[i].allowable[j]!=""){
                    allowable_defined++;
                }
                else{
                    break;
                }
            }
            
            //Make assignments to opt_data 
            opt_data *tp= new(std::nothrow) opt_data;
            if (tp==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            tp->type=param[i].type;
            tp->required=param[i].required;
            
            
            //Process restricted/allowable secondary options parameters
            if (allowable_defined>0){
                tp->restricted =true;
                for(int k=0;k<allowable_defined;k++){
                    tp->flags[param[i].allowable[k]]=false;
                }
            }
            else {
                tp->restricted =false;
            }
            
            
            //Determine presets and assign
            if (param[i].preset.compare("")!=0){
                if (tp->type==OPT_INT){
                    
                    double tempValue;
                    std::istringstream input(param[i].preset);
                    if(!(input>>tempValue)){
                        std::cerr << "Couldn't convert preset to numerical integer: " << param[i].preset << std::endl;
                        tp->set=false;
                    }
                    else{
                        tp->data.i=tempValue;
                        tp->set=true;
                    }
                }
                else if (tp->type==OPT_DOUBLE){
                    
                    double tempValue;
                    std::istringstream input(param[i].preset);
                    if(!(input >> tempValue)){
                        std::cerr << "Couldn't convert preset to numerical double: " << param[i].preset << std::endl;
                        tp->set=false;
                    }
                    else{
                        tp->data.d=tempValue;
                        tp->set=true;
                    }
                }
                else if (tp->type==OPT_STRING){
                    tp->str=param[i].preset;
                    tp->set=true;
                }
                else if (tp->type==OPT_FLAG){
                    tp->set=true;
                    for(int k=0;k<param[i].preset.size();k++){
                        if (k>=MAX_ALLOWABLE){
                            std::cerr << "More presets than allowed OPT_FLAG" <<std::endl;
                            exit(0);
                        }
                        char c=param[i].preset[k];
                        if (c=='1'){
                            tp->flags[param[i].allowable[k]]=true;
                        }
                        else{
                            tp->flags[param[i].allowable[k]]=false;
                        }
                    }
                }
            }
            
            
            //Assign the possible tags to map. Linked to pointer to opt_data
            for(size_t j=0; j<tags.size();j++){
                //cout << tags[j] << endl;
                opts[tags[j]]=tp;
                alternatives[tags[j]]=tags[0];
            }
        }
        
        return;
    }



    //!Parse the commandline arguments and save them in the options clas
    //!\param argc Number of commandline arguments;
    //!\param argv Commandline argurments
    void options::parse_commandline(int argc, const char * argv[]){
        
        for(int i=1;i<argc;i++){
            std::string tag=argv[i];
            if (tag[0]!='-'){
                std::cout << "Unpaired/Unknown commandline argument:\t" << tag <<std::endl;
                std::cout << usage << std::endl; //Print Usage statement
                exit(1);
            }
            else if (opts.count(tag)){
                std::string secondary;
                opts[tag]->set=true;
                switch (opts[tag]->type) {
                    case OPT_NONE:
                        break;
                    case OPT_INT:
                        if (i+1==argc){
                            std::cerr << "Command:\t" << tag <<" is missing a secondary command line argument\n";
                            std::cout << usage << std::endl; //Print Usage statement
                            exit(2);
                        }
                        //Probably need to check to make sure that the argument is a number;
                        opts[tag]->data.i=atoi(argv[++i]);
                        break;
                    case OPT_DOUBLE:
                        if (i+1==argc){
                            std::cerr << "Command:\t" << tag <<" is missing a argument\n";
                            std::cout << usage << std::endl; //Print Usage statement
                            exit(2);
                        }
                        
                        opts[tag]->data.i=atof(argv[++i]);
                        break;
                    case OPT_STRING:
                        //Check that secondary doesn't start with a dash;
                        if (i+1==argc){
                            secondary="";
                        }
                        else{
                            secondary=argv[i+1];
                        }
                        
                        if (secondary[0]=='-'){
                            opts[tag]->str="";
                        }
                        else{
                            i++;
                            if (opts[tag]->restricted){
                                if (opts[tag]->flags.count(secondary)){
                                    opts[tag]->str=secondary;
                                }
                                else{
                                    std::cerr << "Secondary option for:\t" << tag << "  ->" << secondary << " is not allowed.\n";
                                    std::cout << usage << std::endl; //Print Usage statement
                                    exit(2);
                                }
                            }
                            else{
                                opts[tag]->str=secondary;
                            }
                        }
                        break;
                    case OPT_FLAG:
                        if (i+1==argc){
                            std::cerr << "Command:\t" << tag <<" is missing a secondary command line argument\n";
                            std::cout << usage << std::endl; //Print Usage statement
                            exit(2);
                        }
                        
                        secondary=argv[++i];
                        if (opts[tag]->flags.count(secondary)){
                            opts[tag]->flags[secondary]=true;
                        }
                        else{
                            std::cerr << "Secondary option for:\t" << tag << "  ->" << secondary << " is not allowed.\n";
                            std::cout << usage << std::endl; //Print Usage statement
                            exit(2);
                        }
                        break;
                    default:
                        break;
                }
            }
        }
        
        //check that all required are now set;
        bool set=true;
        std::map<std::string,bool> RequiredNotSet;
        std::map<std::string,opt_data*>::iterator it;
        
        for(it=opts.begin();it!=opts.end();it++){
            if ((*it).second->required && !(*it).second->set){
                set=false;
                RequiredNotSet[alternatives[(*it).first]]=false;
            }
        }
        
        if (!set){
            std::map<std::string,bool>::iterator it;
            for(it=RequiredNotSet.begin();it!=RequiredNotSet.end();it++){
                std::cout << "Required option:\t" << (*it).first <<" not set on command-line\n";
            }
            std::cout << usage <<std::endl;
            exit(1);
        }
        
        return;
    }

    
    //!Get the integer value for option
    //!\param [in] param the option name
    //!\param [out] value integer is to be returned to
    //!\return true if option existed and is returned
    //!\return false if option does not exist
    bool options::getopt(const char* param,int &value){
        if (opts.count(param)){
            if (opts[param]->set && opts[param]->type==OPT_INT){
                value=opts[param]->data.i;
                return true;
            }
            else{
                std::cerr << "Option: " << param << " not set or not of type OPT_INT\n";
                return false;
            }
        }
        else{
            std::cerr << "Option: " << param <<" type doesn't exist\n";
            return false;
        }
    }

    
    //!Get double type option values
    //!\param [in] param the option name
    //!\param [out] value double is to be returned to
    //!\return true if option existed and is returned
    //!\return false if option does not exist
    bool options::getopt(const char* param,double &value){
        if (opts.count(param)){
            if (opts[param]->set && opts[param]->type==OPT_DOUBLE){
                value=opts[param]->data.d;
                return true;
            }
            else{
                std::cerr << "Option: " << param << " not set or not of type OPT_DOUBLE\n";
                return false;
            }
        }
        else{
            std::cerr << "Option: " << param <<" type doesn't exist\n";
            return false;
        }
    }

    
    //!Get string type options values
    //!\param [in] param the option name
    //!\param [out] value std::string is to be returned to
    //!\return true if option existed and is returned
    //!\return false if option does not exist
    bool options::getopt(const char* param,std::string &value){
        if (opts.count(param)){
            if (opts[param]->set && opts[param]->type==OPT_STRING){
                value=opts[param]->str;
                return true;
            }
            else{
                std::cerr << "Option: " << param << " not set or not of type OPT_STRING\n";
                return false;
            }
        }
        else{
            std::cerr << "Option: " << param <<" type doesn't exist\n";
            return false;
        }
    }

    
    //!Check primary and secondary options and returns value of bool flag
    //! \param param primary option name
    //! \param sec secondary option name
    //! \return true if secondary option is set
    //! \return false if secondary option is not set
    bool options::getopt(const char* param,const char* sec){
        std::string secondary=sec;
        if (opts.count(param)){
            if (opts[param]->set && opts[param]->type==OPT_FLAG){
                if (opts[param]->flags.count(secondary)){
                    return opts[param]->flags[secondary];
                }
                else {
                    std::cerr << "Secondary command: " << secondary << " doesn't exist.\n";
                    return false;
                }
            }
            else{
                //cerr << "Option: " << param << " not set or not of type OPT_FLAG\n";
                return false;
            }
        }
        else{
            std::cerr << "Option: " << param <<" type doesn't exist\n";
            exit(2);
        }
    }


    
    //!Returns whether an option is set 
    //!Used also for OPT_NONE option types
    //!\param param parameter name
    //!\return true if option is set or exists in commandline argument
    //!\return false if option is not set or didn't exist in commandline argument
    bool options::optset(const char* param){
        if (opts.count(param)){
            return opts[param]->set;
        }
        else{
            std::cerr << "Option: " << param <<" type doesn't exist\n";
            return false;
        }
    }

    //!Get string value for option
    //!If option is not an OPT_STRING, this will produce an error
    //!\param param parameter name
    //!\return std::string 
    std::string &options::sopt(const char * param){
        if (opts.count(param)){
            if (opts[param]->set && opts[param]->type==OPT_STRING){
                return opts[param]->str;
            }
            else{
                std::cerr << "Option: " << param << " not set or not of type OPT_STRING\n";
                exit(2);
            }
        }
        else{
            std::cerr << "Option: " << param <<" type doesn't exist\n";
            exit(2);
        }
    }


    //!Get integer value for option
    //!If option is not an OPT_INT, this will produce an error
    //!\param param parameter name
    //!\return int
    int &options::iopt(const char * param){
        if (opts.count(param)){
            if (opts[param]->set && opts[param]->type==OPT_INT){
                return opts[param]->data.i;
            }
            else{
                std::cerr << "Option: " << param << " not set or not of type OPT_INT\n";
                exit(2);
            }
        }
        else{
            std::cerr << "Option: " << param <<" type doesn't exist\n";
            exit(2);
        }
    }

    //!Get double value for option
    //!If option is not an OPT_DOUBLE, this will produce an error
    //!\param param parameter name
    //!\return double
    double &options::dopt(const char * param){
        if (opts.count(param)){
            if (opts[param]->set && opts[param]->type==OPT_INT){
                return opts[param]->data.d;
            }
            else{
                std::cerr << "Option: " << param << " not set or not of type OPT_DOUBLE\n";
                exit(2);
            }
        }
        else{
            std::cerr << "Option: " << param <<" type doesn't exist\n";
            exit(2);
        }
    }
    
}
