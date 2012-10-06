//
//  CountSequence.cpp
//  StochHMM
//
//  Created by Paul Lott on 3/5/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#include <iostream>
#include <vector>
#include <string>
#include "options.h"
#include "sequences.h"
#include "sequence.h"
#include "Counter.h"
#include "usage.h"

using namespace std;
//
//opt_parameters commandline[]={
//    //Help
//    {"-help:-h"     ,OPT_NONE       ,false      ,"",    {}},
//    {"-seq:-s"      ,OPT_STRING     ,true       ,"",    {}},
//    {"-order:-o"    ,OPT_INT        ,false      ,"0",   {}},
//    {"-output"      ,OPT_STRING     ,false      ,"",    {}},
//    {"-direction"   ,OPT_FLAG       ,false      ,"forward", {"forward","reverse","both"}},
//    {"-periodic"    ,OPT_INT        ,false      ,"",    {}},
//    {"-pwm"         ,OPT_NONE       ,false      ,"",    {}},
//    {"-upstream"    ,OPT_INT        ,false      ,"",    {}},
//    {"-downstream"  ,OPT_INT        ,false      ,"",    {}},
//    {"-pseudocount" ,OPT_INT        ,false      ,"",    {}},
//    {"-withWords"   ,OPT_NONE       ,false      ,"",    {}},
//    {"-masked"      ,OPT_INT        ,false      ,"",    {}},
//    {"-quiet"       ,OPT_NONE       ,false      ,"",    {}}
//};
//
//int opts_size=sizeof(commandline)/sizeof(commandline[0]);
//
//
//options opt;
//bool quiet;

int main (int argc, const char * argv[]) {
    
//    opt.set_parameters(commandline,opts_size,usage);
//	opt.parse_commandline(argc,argv);
//    
//    //Parse Options
//    countDirection direction= (opt.getopt("-direction", "forward")) ? FORWARD :
//    (opt.getopt("-direction", "reverse")) ? REVERSE :
//    (opt.getopt("-direction", "both")) ? BOTH : FORWARD;
//    
//    
//    string seqFile=opt.sopt("-seq");
//    bool withWords=opt.optset("-withWords");  //Do we want words printed
//    
//    int period=0;
//    int order=0;
//    
//    countType type(PERIODIC);
//    
//    if (!opt.optset("-periodic") && !opt.optset("-pwm") && !opt.optset("-masked")){
//        cout << "Defaulting to Periodic Counting with period =1\n";
//        type=PERIODIC;
//        period=1;
//    }
//    else if(opt.optset("-periodic") && opt.optset("-pwm")){
//        cout << "Both PWM and Periodic set on commandline. Only able to do one at a time.  Defaulting to Periodic Counting\n";
//        period=opt.iopt("-periodic");
//        type=PERIODIC;
//    }
//    else if (opt.optset("-periodic")){
//        period=opt.iopt("-periodic");
//        type=PERIODIC;
//    }
//    else if (opt.optset("-pwm")){
//        type=PWM;
//    }
//    else if (opt.optset("-masked")){
//        //seq.setMasked();
//        period=opt.iopt("-masked");
//        type=PERIODIC;
//    }
//    
//    if (type!=PWM && opt.optset("-order")){
//        order = opt.iopt("-order");
//    }
//    
//    int upstream=(opt.optset("-upstream")) ? opt.iopt("-upstream") : 0;
//    int downstream = (opt.optset("-downstream")) ? opt.iopt("-downstream") : 0;
//    int pseudocount= (opt.optset("-pseudocount")) ? opt.iopt("-pseudocount") : 0;
//    
//    string outputFilename=(opt.optset("-output"))? opt.sopt("-output") : "";
//    
//    if (opt.optset("-quiet")){ 
//        quiet=true;
//    }
//    else{
//        quiet=false;
//    }
    
//    //Initialize Sequence file;
//    sequence seq(opt.sopt("-seq"));
//    
//    if (opt.optset("-masked")){
//        seq.setMasked();
//    }
//    
//    
//    counts cnt(direction,type,period,upstream,downstream,order,pseudocount);
//    
//    while(seq.next()){
//        cnt.count(seq);
//    }
//    
//    cnt.print(withWords);
	
	exit(0);
}