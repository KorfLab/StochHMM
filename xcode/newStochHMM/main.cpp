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
#include <string>
#include <iomanip>
#include "hmm.h"
#include "sequence.h"
#include "seqTracks.h"
#include "new_trellis.h"
#include "state.h"
#include "seqTracks.h"
#include "traceback_path.h"
#include "options.h"
#include <time.h>
#include "newStochHMM_Usage.h"
using namespace StochHMM;

trellis* perform_decoding(model*, sequences*);
trellis* perform_posterior(model*, sequences*);

void import_model(model&);
void import_sequence(model&);

traceback_path* perform_traceback(trellis*);
multiTraceback* perform_stochastic_traceback(trellis*);

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
    {"-posterior"   ,OPT_NONE       ,false  ,"",    {}},
	//Output Files and Formats
    {"-gff:-g"      ,OPT_STRING     ,false  ,"",    {}},
    {"-path:-p"     ,OPT_STRING     ,false  ,"",    {}},
    {"-label:-l"    ,OPT_STRING     ,false  ,"",    {}},
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
	
	while (job != NULL){
		//Print sequences if -debug seq option defined
		if (opt.isFlagSet("-debug","seq")){
			job->getSeqs()->print();
		}
		
		
		if (opt.isSet("-posterior")){
			//trellis* trell = perform_posterior(job->getModel(), );
			
			trellis* trell = new trellis(&hmm,job->getSeqs());
			trell->simple_posterior_second();
//			double post = trell -> getProbOfSeq();
//			
//			float_2D* forw  = trell->getForwardTable();
//			float_2D* backw = trell->getBackwardTable();
//			
//			double forward_sum(0);
//			for (size_t i=0;i<forw->size();i++){
//				for(size_t j=0;j< (*forw)[i].size();j++){
//					std::cout << exp((*forw)[i][j]) << std::endl;
//					forward_sum +=exp((*forw)[i][j]);
//				}
//			}
//			std::cout << "Forward Sum:\t" << forward_sum << std::endl;
//			
//			
//			double backward_sum(0);
//			for (size_t i=0;i<backw->size()-1;i++){
//				for(size_t j=0;j< (*backw)[i].size();j++){
//					std::cout << exp((*backw)[i][j]) << std::endl;
//					backward_sum +=exp((*backw)[i][j]);
//				}
//			}
//			
//			double temp = forward_sum + backward_sum;
//			std::cout << "Backward Sum:\t" << backward_sum << std::endl;
//			
//			std::cout << "Combined Sum:\t" << temp << std::endl;
//			std::cout << "Prob:\t" << post << std::endl;
			
			
			if (trell == NULL){
				std::cerr << "No posterior performed\n";
				exit(0);
			}
			
//			std::vector<std::vector<double> >* posterior_table = trell->getPosteriorScores();
//			print_posterior(posterior_table,job->getModel(),job->getSeqs()->getHeader(),trell->getProbOfSeq());
//			delete trell;
		}
		
		if(opt.isSet("-viterbi")){
			//Perform decoding
			trellis *trell = perform_decoding(job->getModel(), job->getSeqs());
			
			if (trell == NULL){
				std::cerr << "No decoding performed\n";
				exit(0);
			}
			
			
			traceback_path* tb(NULL);
			tb = perform_traceback(trell);  //job->getSeqs()->getHeader()
			print_output(tb, job->getSeqs()->getHeader());
			delete trell;
			delete tb;
		}
		
		

		job = jobs.getJob();
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
    
	
    trell =  new(std::nothrow) trellis(hmm,seqs);
	trell->simple_viterbi();
    
    return trell;
}

trellis* perform_posterior(model* hmm, sequences* seqs){
	trellis* trell = new trellis(hmm,seqs);
	trell->simple_posterior();
	
	return trell;
}


void print_output(traceback_path* tb, std::string& header){
    
    bool previous(true);
    
    if (opt.isSet("-gff")){
        std::cout << "Score: " << tb->getScore() << std::endl;
        tb->print_gff(header);
        previous=false;
    }
    
    if (opt.isSet("-label")){
		std::cout << ">" << header ;
        std::cout << "\tScore: " << tb->getScore() << std::endl;
        tb->print_label();
        previous=false;
    }
    
    if (opt.isSet("-path") || previous){
		std::cout << ">" << header ;
        std::cout << "\tScore: " << tb->getScore() << std::endl;
        tb->print_path();
    }
    
    return;
}

traceback_path* perform_traceback(trellis* trell){
	traceback_path* path = new(std::nothrow) traceback_path(trell->get_model());
	trell->traceback(*path);
    return path;
}


void print_posterior(std::vector<std::vector<double> >* table , model* hmm, std::string& header,double val){
	std::cout <<"Posterior Probabilities Table" << std::endl;
	std::cout <<"Sequence:\t" << header << std::endl;
	std::cout <<"Probability of Sequence:Natural Log'd\t" << val << std::endl;
	
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
			std::cout << std::setprecision(3) << exp((*table)[position][st]);
		}
		std::cout << std::endl;
	}
	
	return;
	
}


