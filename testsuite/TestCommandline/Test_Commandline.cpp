//
//  Test_Commandline.cpp
//  StochHMMme_rc
//
//  Created by Paul C Lott on 7/1/11.
//  Copyright 2011 University of California, Davis. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <time.h>
#include "options.h"
#include "TC_usage.h"

using namespace std;


#pragma mark Options
//Sets the command-line options for the program
opt_parameters commandline[]={
    //Help
    {"-help:-h"     ,OPT_NONE       ,false  ,"",    {}},
    //Required
    {"-model:-m"    ,OPT_STRING     ,true   ,"",    {}},
    {"-seq:-s:-track",OPT_STRING    ,true   ,"",    {}},
    //HMMER Search file
    {"-hmmer"       ,OPT_STRING     ,false  ,"",    {}},
    {"-scale"       ,OPT_DOUBLE     ,false  ,"0.2",{}},
    //Debug
    {"-debug:-d"    ,OPT_FLAG ,false  ,"",    {"model","seq","paths"}},
    {"-pathThresh"  ,OPT_DOUBLE     ,false  ,"",    {}},
    {"-maxTime"     ,OPT_INT        ,false  ,"3600",{}},
    //Non-Stochastic Decoding
    {"-viterbi"     ,OPT_NONE       ,false  ,"",    {}},
    {"-nbest"       ,OPT_INT        ,false  ,"1",   {}},
    //Stochastic Decoding
    {"-stoch"       ,OPT_FLAG ,false  ,"",    {"forward","viterbi"}},
    {"-repetitions:-rep",OPT_INT    ,false  ,"1000",{}},
    {"-report:-rpt" ,OPT_INT        ,false  ,"1000",{}},
    //Output Files and Formats
    {"-gff:-g"      ,OPT_STRING     ,false  ,"",    {}},
    {"-path:-p"     ,OPT_STRING     ,false  ,"",    {}},
    {"-label:-l"    ,OPT_STRING     ,false  ,"",    {}},
    {"-posterior:-post",OPT_STRING  ,false  ,"",    {}},
    {"-heat"        ,OPT_STRING     ,false  ,"",    {}},
    {"-trellis"     ,OPT_STRING     ,false  ,"",    {}},
    {"-traceback"   ,OPT_STRING     ,false  ,"",    {}},
    {"-sidd"        ,OPT_STRING     ,false  ,"",    {}},
    {"-oReport"     ,OPT_STRING     ,false  ,"",    {}},
    //Additional parameters
    {"-threshold"   ,OPT_DOUBLE     ,false  ,"0.001",{}},
    {"-width"       ,OPT_INT        ,false  ,"80",  {}}
};

int opt_size=sizeof(commandline)/sizeof(commandline[0]);

options opt;

int main (int argc, const char * argv[]) {
    
    opt.set_parameters(commandline,opt_size,usage);
	opt.parse_commandline(argc,argv);
	
	return 0;
}
