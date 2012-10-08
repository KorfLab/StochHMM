/*
 *  options.h
 *  HMM
 *
 *  Created by Paul C Lott on 8/13/09.
 *  Copyright 2009 University of California, Davis. All rights reserved.
 *
 */

#ifndef OPTIONS_H
#define OPTIONS_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

namespace StochHMM{


//Options Types
enum optionType {UNDEF_OPT , OPT_NONE , OPT_INT , OPT_DOUBLE, OPT_STRING , OPT_FLAG };
    
#define MAX_ALLOWABLE 10



    
    /*!User defined option parameters
     opt_parameters set up in beginning of your program and defines what options
     are accepted by the program.   
     Used for designing allowed commandline options
     option Types:  
     UNDEF_OPT   - initialized value not defined 
     OPT_NONE    - Functions as a bool (was it set or not)
     OPT_INT     - Integer value supplied as option
     OPT_DOUBLE  - Double value supplied as option
     OPT_STRING  - String value supplied as option
     OPT_FLAG    - Predefined allowed values (each functions as a bool)
     
     Colums of table:  
     Option tag- all possible alternatives separated by colon
     Option type- what type will it be passed
     Option required - If option is required to run the it should be true, else false
     Option default - Default value for option
     Secondary Tag - List of possible secondary tags.  
     
     Usage:
     opt_parameters commandline[] - {
     {"-help:-h" ,OPT_NONE, false, "", {}},
     {"-model:-m",OPT_STRING, true, "",{}},
     {"-play:-p", OPT_FLAG, false, "", {"Ball", "Train"}
     }
     int opts_size=sizeof(commandline)/sizeof(commandline[0]); //required to automatically figure # of options defined
     
     options opt;
     int main (int argc,const char *argv[]){
     opt.set_paramenters(commandline,opt_size,usage); //sets up the options you've defined in commandline variable
     opt.parse_commandline(argc,argv); //Parses the commandline and assigns the variables;
     //To get the string value from model
     string filename=sopt("-model");  
     //To see if scale was set
     if (optset("-scale")) 
     double scale=sopt("-scale"); //Return value of scale option
     return 1;
     };
     */
    struct opt_parameters{
        std::string commandline_tag;  //Tag used on the commandline
        optionType type;  //Option Tag type
        bool required;  //Is the option required to be set in order to run
        std::string preset;  //Preset for the option
        std::string allowable[MAX_ALLOWABLE]; //allowable options secondary tags
    };

        
        
    // TODO:  Need to make opt_data private and allow options to access
    //Actual commandline options data stored in this struct
    
    //!Store the actual command line options data for each option defined by user
    class opt_data{
    public:
        opt_data();
        optionType type;           //what datatype does the option store
        bool set;           //Has the parameter been set (OPT_NONE AND EVERYTHING ELSE)
        bool required;      //Is the parameter required, else exit
        bool restricted;    //Are the entries restricted to allowable
        union {   //If its a double or int for OPT_INT or OPT_DOUBLE;
            int i;
            double d;
        }data;
        std::string str;     //String value for OPT_STRING
        std::map<std::string,bool> flags;  //Allowable flags or entries  for OPT_FLAG
    };



    //!Parse, Store, and access commandline options data
    class options{
    public:
        options();
        
        void set_parameters(opt_parameters*,int,const char*);  //Setup the parameter you've defined in the options variable 
        void parse_commandline(int, const char**);   //Parses the commandline
        
        bool getopt(const char*,int&);          //Get int for option type OPT_INT
        bool getopt(const char*,double&);       //Get double for option OPT_DOUBLE
        bool getopt(const char*,std::string&);       //Get string from option OPT_STRING

    //Accessor functions
        bool getopt(const char*, const char*);  //Get status of secondary option OPT_FLAG
        bool optset(const char*); // check to see if options has been set
        
        
        std::string  &sopt(const char*);  //return string for option
        int     &iopt(const char*);  //return int for option
        double  &dopt(const char*);  //return double for option
        
    private:
        std::map<std::string,opt_data*> opts;  //map each option to string tag
        std::map<std::string,std::string> alternatives; //Store the tag alternatives as keys and main as value
        const char* usage;

    };

}
#endif /*OPTIONS_H*/
