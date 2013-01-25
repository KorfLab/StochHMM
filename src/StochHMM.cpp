//Copyright (c) 2007-2012 Paul C Lott
//University of California, Davis
//Genome and Biomedical Sciences Facility
//UC Davis Genome Center
//Ian Korf Lab
//Website: www.korflab.ucdavis.edu
//Email: lottpaul@gmail.com
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of
//this software and associated documentation files (the "Software"), to deal in
//the Software without restriction, including without limitation the rights to
//use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
//the Software, and to permit persons to whom the Software is furnished to do so,
//subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
//FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
//COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
//IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "hmm.h"
#include "sequence.h"
#include "seqTracks.h"
#include "trellis.h"
#include "options.h"
#include <time.h>
#include "stochhmmUsage.h"
using namespace StochHMM;

trellis* perform_decoding(model*, sequences*);
trellis* perform_posterior(model*, sequences*);

void import_model(model&);
void import_sequence(model&);

multiTraceback* perform_stochastic_traceback(trellis*);
void perform_nth_traceback(trellis*, std::vector<traceback_path>&);
traceback_path* perform_traceback(trellis*);


void print_output(multiTraceback*, std::string&);
void print_output(std::vector<traceback_path>&, std::string&);
void print_output(traceback_path*, std::string&);
void print_posterior(std::vector<std::vector<double> >*, model*, std::string&,double);


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
    {"-nbest"       ,OPT_INT        ,false  ,"3",   {}},
    {"-posterior"   ,OPT_NONE       ,false  ,"",    {}},
//Stochastic Decoding
    {"-stochastic"  ,OPT_FLAG       ,false  ,"",    {"viterbi","forward","posterior"}},
    {"-repetitions:-rep",OPT_INT    ,false  ,"1000",{}},
//Output Files and Formats
    {"-gff:-g"      ,OPT_STRING     ,false  ,"",    {}},
    {"-path:-p"     ,OPT_STRING     ,false  ,"",    {}},
    {"-label:-l"    ,OPT_STRING     ,false  ,"",    {}},
    {"-hits"        ,OPT_STRING     ,false  ,"",    {}},
    {"-trellis"     ,OPT_STRING     ,false  ,"",    {}},
};

int opt_size=sizeof(commandline)/sizeof(commandline[0]);

options opt;


seqTracks jobs;


int main(int argc, const char * argv[])
{
    
    //Parse commandline arguments
    opt.set_parameters(commandline,opt_size,usage);
    opt.parse_commandline(argc,argv);
    
    
    model hmm;
    //Check that model argument is defined and import the model
    import_model(hmm);
    
    if (opt.isFlagSet("-debug", "model")){
        hmm.print();
        hmm.writeGraphViz("Temp.viz", true);
    }
   
    
    //Check and import sequence(s)
    import_sequence(hmm);
    
    seqJob *job=jobs.getJob();
    
    
    //Print sequences if -debug seq option defined
    if (opt.isFlagSet("-debug","seq")){
        job->getSeqs()->print();
    }
    
    
	if (opt.isSet("-posterior")){
		trellis* trell = perform_posterior(job->getModel(), job->getSeqs());
		
		//If we need a posterior traceback b/c path,label,or GFF is defined
		if (opt.isSet("-gff") || opt.isSet("-path") || opt.isSet("-label")){
			traceback_path* path = trell->posteriorTraceback();
			print_output(path, job->getSeqs()->getHeader());
		}
		else{
			std::vector<std::vector<double> >* posterior_table = trell->getPosterior();
			print_posterior(posterior_table,job->getModel(),job->getSeqs()->getHeader(),trell->getProbOfSeq());
			
		}
		
		exit(1);
	}
	
    //Perform decoding
    trellis *trell = perform_decoding(job->getModel(), job->getSeqs());
    
    if (trell == NULL){
        std::cerr << "No decoding performed\n";
        exit(0);
    }
    
    if (opt.isSet("-trellis")){
        trell->print();
    }
    
    
    multiTraceback* multi_tb(NULL);
    traceback_path* tb(NULL);
    std::vector<traceback_path> nth;
    
    //Perform Tracebacks and print outputs
    if (opt.isSet("-stochastic")){
        multi_tb = perform_stochastic_traceback(trell);
        print_output(multi_tb, job->getSeqs()->getHeader());
    }
    else if (opt.isSet("-nbest")){
        perform_nth_traceback(trell,nth);
        print_output(nth, job->getSeqs()->getHeader());
    }
    else{
        tb = perform_traceback(trell);  //job->getSeqs()->getHeader()
        print_output(tb, job->getSeqs()->getHeader());
    }

    
    if (opt.isSet("-trellis")){
        trell->print();
    }
    return 0;
    
}

void import_model(model& hmm){
    if (!opt.isSet("-model")){
        std::cerr <<"No model file provided.\n" << usage << std::endl;
    }
    else{
        hmm.import(opt.sopt("-model"));
    }
    
    //If -debug model is defined then print model to stdout;
    if (opt.isFlagSet("-debug","model")){
        hmm.print();
    }
}

void import_sequence(model& hmm){
    
    if (!opt.isSet("-seq")){
        std::cerr << "No sequence file provided.\n" << usage << std::endl;
    }
    else{
        jobs.loadSeqs(hmm, opt.sopt("-seq"), FASTA);
    }
}


trellis* perform_decoding(model* hmm, sequences* seqs){
    trellis *trell;
    
    bool stochastic = (opt.isSet("-stochastic")) ? true : false;
    
    
    bool viterbi = (opt.isFlagSet("-stochastic", "viterbi") || opt.isSet("-viterbi")) ? true : false;
    bool forward = (opt.isFlagSet("-stochastic", "forward")) ? true : false;
    bool posterior = (opt.isFlagSet("-stochastic", "posterior") || opt.isSet("-posterior")) ? true : false;
    bool nbest = (opt.isSet("-nbest")) ? true : false;
   
    trell = (stochastic) ? new trellis(hmm,seqs,STOCH) : new trellis(hmm,seqs,SIMPLE);
    
    if (nbest){
        trell->nthViterbi(opt.iopt("-nbest"));
    }
    else if (viterbi){
        if (forward || posterior){
            trell->decodeAll();
        }
        else{
            trell->viterbi();
        }
    }
    else if (posterior){
        trell->posterior();
    }
    else if (forward){
        trell->forward();
    }
    else{
        std::cerr << usage << "\nNo Decoding option set\n";
        return NULL;
    }
    
    return trell;
}

trellis* perform_posterior(model* hmm, sequences* seqs){
	trellis* trell = new trellis(hmm,seqs,SIMPLE);
	trell->posterior();
	
	return trell;
}


multiTraceback* perform_stochastic_traceback(trellis* trell){
    
    multiTraceback* multi_tb =  (opt.isFlagSet("-stochastic" , "forward"))     ? trell->stochasticForwardTraceback(opt.iopt("-rep")) :
                                (opt.isFlagSet("-stochastic" , "posterior"))   ? trell->stochasticPosteriorTraceback(opt.iopt("-rep")) :
                                                                                 trell->stochasticViterbiTraceback(opt.iopt("-rep"));
    return multi_tb;
}



void perform_nth_traceback(trellis* trell, std::vector<traceback_path>& nth){
    for(size_t number_of_paths = 0; number_of_paths<opt.iopt("-nbest"); number_of_paths++){
        traceback_path* pth = trell->nth_traceback(number_of_paths);
        nth.push_back(*pth);
    }
    return;
}



traceback_path* perform_traceback(trellis* trell){
    traceback_path* path= trell->traceback();
    
    return path;
}


void print_output(multiTraceback* tb, std::string& header){
    
    tb->finalize();
    
    bool previous(true);
    
    if (opt.isSet("-hits")){
        tb->print_hits();
        previous=false;
    }
    
    if (opt.isSet("-gff")){
        tb->print_gff(header);
        previous=false;
    }
    
    if (opt.isSet("-label")){
        tb->print_label();
        previous=false;
    }
    
    //Print path by default if nothing else is set
    if (opt.isSet("-path") || previous){
        tb->print_path();
    }
    
    return;
}


void print_output(traceback_path* tb, std::string& header){
    
    bool previous(true);
    
    if (opt.isSet("-gff")){
		if (!opt.isSet("-posterior")){
			std::cout << "Viterbi Score: " << tb->getScore() << std::endl;
		}
        
        tb->print_gff(header);
        previous=false;
    }
    
    if (opt.isSet("-label")){
        if (!opt.isSet("-posterior")){
			std::cout << "Viterbi Score: " << tb->getScore() << std::endl;
		}
        tb->print_label();
        previous=false;
    }
    
    if (opt.isSet("-path") || previous){
        if (!opt.isSet("-posterior")){
			std::cout << "Viterbi Score: " << tb->getScore() << std::endl;
		}
        tb->print_path();
    }
    
    return;
}


void print_output(std::vector<traceback_path>& tb, std::string& header){
    
    bool previous(true);
    
    if (opt.isSet("-gff")){
        for(size_t i=0;i<tb.size();i++){
            std::cout << "Nth Best: " << i+1 << std::endl;
            std::cout << "Score: " << tb[i].getScore() << std::endl;
            tb[i].print_gff(header);
        }
        previous=false;
    }
    
    if (opt.isSet("-label")){
        for(size_t i=0;i<tb.size();i++){
            std::cout << "Nth Best: " << i+1 << std::endl;
            std::cout << "Score: " << tb[i].getScore() << std::endl;
            tb[i].print_label();
        }
        previous=false;
    }
    
    if (opt.isSet("-path") || previous){
        for(size_t i=0;i<tb.size();i++){
            std::cout << "Nth Best: " << i+1 << std::endl;
            std::cout << "Score: " << tb[i].getScore() << std::endl;
            tb[i].print_path();
        }
    }
    
    return;
}

void print_posterior(std::vector<std::vector<double> >* table , model* hmm, std::string& header,double val){
	std::cout <<"Posterior Probabilities Table" << std::endl;
	std::cout <<"Sequence:\t" << header << std::endl;
	std::cout <<"Probability of Sequence:\t" << exp(val) << std::endl;

	size_t state_size = hmm->state_size();
	
	//Print State names for header
	std::cout << "Position";
	for(size_t i=0; i < state_size; i++){
		std::cout << "\t";
		std::cout << hmm->getStateName(i);
	}
	std::cout << std::endl;
	
	for(size_t position = 0; position < table->size(); position++){
		std::cout << position+1;
		for (size_t st = 0 ; st < state_size ; st++){
			std::cout << "\t";
			std::cout << std::setprecision (3) << (*table)[position][st];
		}
		std::cout << std::endl;
	}
	
	return;
	
}



