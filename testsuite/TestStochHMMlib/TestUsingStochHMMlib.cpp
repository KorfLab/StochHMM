//
//  TestUsingStochHMMlib.cpp
//  StochHMMme_rc
//
//  Created by Paul Lott on 8/13/11.
//  Copyright 2011 University of California, Davis. All rights reserved.
//

#include "StochHMMlib.h"
#include <iostream>
#include <math.h>
#include <time.h>

using namespace std;

#pragma mark Options
//Sets the command-line options for the program
opt_parameters commandline[]={
    //Help
    {"-help:-h"     ,OPT_NONE       ,false  ,"",    {}},
    //Required
    {"-model:-m"    ,OPT_STRING     ,true   ,"",    {}},
    {"-models"      ,OPT_STRING     ,false  ,"",    {}},  // TODO: Split out models into models option rather than looking at file first
    {"-seqs"        ,OPT_STRING     ,false  ,"",    {}},  // TODO: create sequences file if using multiple sequences for separate tracks
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
seqTracks seqs;  //Initialize Queue

#ifdef SIDD
sidd_param parameters;
#endif

const char usage[] = "Test Usage"; 

#pragma mark MAIN
int main (int argc, const char * argv[]) {
    multiHMM models;
    opt.set_parameters(commandline,opt_size,usage);
	opt.parse_commandline(argc,argv);
	
    if (!opt.optset("-seq")){
		cerr << "No sequence file defined in command line input\n Please check syntax.\n";
		//cout << usage << endl;
	}
	
	if (!opt.optset("-model")){
		cerr << "No model/models file defined in command line input\n Please check syntax.\n";
		//cout << usage << endl;
	}
    else{
        models.importModels(opt.sopt("-model"));
    }
    
    //for(int i=0;i<models.size();i++){
    //    models[i].print_model();
    //}
    
    seqs.load_fastq(models, opt.sopt("-seq"));
    
    //int count=1;
    while(1){
        seqSet *job = seqs.getJob();
        job->getSidd();
        
        //cout << "Seq:\t" << count <<endl;
        //count++;
        job->print_seq();
        //cout << job->model <endl;
        delete job;
    }
    
    
	return 0;
}
