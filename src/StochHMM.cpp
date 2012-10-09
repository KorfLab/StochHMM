#include <iostream>
#include <string>
#include "hmm.h"
#include "sequence.h"
#include "seqTracks.h"
#include "trellis.h"
#include "options.h"
#include <time.h>
#include "stochhmmUsage.h"
using namespace StochHMM;

trellis* perform_decoding(model*, sequences*);

void perform_stochastic_traceback(trellis*);
void perform_nth_traceback(trellis*);
void perform_traceback(trellis*,std::string&);


#pragma mark Options
//Sets the command-line options for the program
opt_parameters commandline[]={
//Help
    {"-help:-h"     ,OPT_NONE       ,false  ,"",    {}},
//Required
    {"-model:-m"    ,OPT_STRING     ,true   ,"",    {}},
    {"-seq:-s:-track",OPT_STRING    ,true   ,"",    {}},
//Debug
    {"-debug:-d"    ,OPT_FLAG ,false  ,"",    {"model","seq","paths"}},
//Non-Stochastic Decoding
    {"-viterbi"     ,OPT_NONE       ,false  ,"",    {}},
    {"-nbest"       ,OPT_INT        ,false  ,"1",   {}},
    {"-posterior"   ,OPT_NONE       ,false  ,"",    {}},
//Stochastic Decoding
    {"-stochastic"      ,OPT_NONE       ,false  ,"",    {}},
    {"-repetitions:-rep",OPT_INT    ,false  ,"1000",{}},
//Output Files and Formats
    {"-gff:-g"      ,OPT_STRING     ,false  ,"",    {}},
    {"-path:-p"     ,OPT_STRING     ,false  ,"",    {}},
    {"-label:-l"    ,OPT_STRING     ,false  ,"",    {}},
    {"-heat"        ,OPT_STRING     ,false  ,"",    {}},
    {"-trellis"     ,OPT_STRING     ,false  ,"",    {}},
};

int opt_size=sizeof(commandline)/sizeof(commandline[0]);

options opt;


int main(int argc, const char * argv[])
{
    
    opt.set_parameters(commandline,opt_size,usage);
    opt.parse_commandline(argc,argv);
    
    model hmm;
    seqTracks jobs;
    
    if (!opt.optset("-model")){
        std::cerr <<"No model file provided.\n" << usage << std::endl;
    }
    else{
        hmm.import(opt.sopt("-model"));
    }
    
    if (opt.getopt("-debug","model")){
        hmm.print();
    }
    
    
    if (!opt.optset("-seq")){
        std::cerr << "No sequence file provided.\n" << usage << std::endl;
    }
    else{
        jobs.loadSeqs(hmm, opt.sopt("-seq"), FASTA);
    }
    
    seqJob *job=jobs.getJob();
    
    
    if (opt.getopt("-debug","seq")){
        job->getSeqs()->print();
        //job->printSeq();
    }
    
    
    //Determine stochastic vs regular
    
    trellis *trell = perform_decoding(job->getModel(), job->getSeqs());
    
    //Outputs
    
    if (opt.optset("-stochastic")){
        std::cout << "Need to finish main implementation of stochastic traceback output\n";
        //perform_stochastic_traceback(trell);
    }
    //else if (opt.optset("-nbest")){
    //    std::cout << "Need to finish main implementation of nth-best viterbi traceback output\n";
    //}
    else{
        perform_traceback(trell,job->getSeqs()->getHeader());
    }

    
    if (opt.optset("-trellis")){
        trell->print();
    }
    return 0;
    
}


trellis* perform_decoding(model* hmm, sequences* seqs){
    trellis *trell = (opt.optset("-stochastic")) ?  new trellis(hmm, seqs, STOCH) :
                                                    new trellis(hmm, seqs, SIMPLE);
    if (opt.optset("-posterior") && opt.optset("-viterbi")){
        trell->decodeAll();
    }
    else if (opt.optset("-posterior")){
        trell->posterior();
    }
    else if (opt.optset("-viterbi")) {
        trell->viterbi();
    }
    else if (opt.optset("-nbest")){
        trell->nthViterbi(opt.iopt("-nbest"));
    }
    else{
        trell->forward();
    }
    
    return trell;
}


void perform_stochastic_traceback(trellis* trell){
    multiTraceback* multi_tb = (opt.optset("-viterbi")) ? trell->stochasticViterbiTraceback(opt.iopt("-rep")) :
                                                          trell->stochasticForwardTraceback(opt.iopt("-rep")) ;
    std::cout << "Need to finish outputs\n";
    return;
}

void perform_traceback(trellis* trell, std::string& name){
    traceback_path* path= trell->traceback();
    
    if (opt.optset("-label") || (!opt.optset("-path") && !(opt.optset("-gff")))){
        path->print_label();
    }
    
    if (opt.optset("-path")){
        path->print_path();
    }
    
    if (opt.optset("-gff")){
        path->print_gff(name);
    }
    
    return; 
}







