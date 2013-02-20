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
#include <time.h>
#include <fstream>
#include "hmm.h"
#include "sequence.h"
#include "seqTracks.h"
#include "trellis.h"
#include "state.h"
#include "seqTracks.h"
#include "traceback_path.h"
#include "options.h"

#include "StochHMM_usage.h"
using namespace StochHMM;

#define STATE_MAX 1024


void import_model(model&);
void import_sequence(model&);

void perform_viterbi_decoding(model* hmm, sequences* seqs);
void perform_nbest_decoding(model* hmm, sequences* seqs);
void perform_posterior(model* hmm, sequences* seqs);
void perform_stochastic_decoding(model* hmm, sequences* seqs);


traceback_path* perform_traceback(trellis&);
multiTraceback* perform_stochastic_traceback(trellis&);

void print_output(multiTraceback*, std::string&);
void print_output(std::vector<traceback_path>&, std::string&);
void print_output(traceback_path*, std::string&);
void print_posterior(trellis&);


#pragma mark Options
//Sets the command-line options for the program
opt_parameters commandline[]={
	//Help
    {"-help:-h"     ,OPT_NONE       ,false  ,"",    {}},
	//Required
    {"-model:-m"    ,OPT_STRING     ,true   ,"",    {}},
    {"-seq:-s:-track",OPT_STRING    ,true   ,"",    {}},
	{"-fastq"		,OPT_NONE		,false	,"",	{}},
	//Debug
    {"-debug:-d"    ,OPT_FLAG		,false  ,"",    {"model","seq","paths"}},
	//Non-Stochastic Decoding
    {"-viterbi"     ,OPT_NONE       ,false  ,"",    {}},
	{"-nbest"       ,OPT_INT        ,false  ,"3",   {}},
    {"-posterior"   ,OPT_STRING		,false  ,"",    {}},
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
    }
	
	if (opt.isFlagSet("-debug", "paths")){
		hmm.writeGraphViz("Model_Path-debug.viz", true);
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
			perform_posterior(&hmm, job->getSeqs());
			return 0;
		}
		else if(opt.isSet("-viterbi")){
			perform_viterbi_decoding(job->getModel(), job->getSeqs());
		}
		else if (opt.isSet("-nbest")){
			perform_nbest_decoding(job->getModel(), job->getSeqs());
		}
		else if (opt.isSet("-stochastic")){
			perform_stochastic_decoding(job->getModel(), job->getSeqs());
		}
		
		job = jobs.getJob();
	}
    
    return 0;
    
}

//Import the model from file
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

//Load sequences from file (Default FASTA)
void import_sequence(model& hmm){
    
    if (!opt.isSet("-seq")){
        std::cerr << "No sequence file provided.\n" << usage << std::endl;
    }
	else if (opt.isSet("-fastq")){
		jobs.loadSeqs(hmm,opt.sopt("-seq"),FASTQ);
	}
    else{
        jobs.loadSeqs(hmm, opt.sopt("-seq"), FASTA);
    }
}

//Perform Viterbi decoding and print the output
void perform_viterbi_decoding(model* hmm, sequences* seqs){
    trellis trell(hmm,seqs);
	trell.viterbi();
	
	traceback_path* tb(NULL);
	tb = perform_traceback(trell);
	print_output(tb, seqs->getHeader());
    return;
}

//Perform nth-best decoding and print the output
void perform_nbest_decoding(model* hmm, sequences* seqs){
	trellis trell(hmm,seqs);
	size_t nth = opt.iopt("-nbest");
	
	trell.naive_nth_viterbi(nth);
	
	for(size_t i=0;i<nth;i++){
		traceback_path path(hmm);
		trell.traceback_nth(path, i);
		print_output(&path, seqs->getHeader());
	}
}


//Perform stochastic decoding
void perform_stochastic_decoding(model* hmm, sequences* seqs){
    
    
    bool viterbi	= (opt.isFlagSet("-stochastic", "viterbi") || opt.isSet("-viterbi")) ? true : false;
    bool forward	= (opt.isFlagSet("-stochastic", "forward")) ? true : false;
    bool posterior	= (opt.isFlagSet("-stochastic", "posterior")) ? true : false;
	
    trellis trell(hmm,seqs);
	
	int repetitions = opt.iopt("-rep");
    
    
    if (viterbi){
		trell.stochastic_viterbi();
		multiTraceback paths;
		trell.stochastic_traceback(paths, repetitions);
		print_output(&paths, seqs->getHeader());
    }
    else if (forward){
		trell.stochastic_forward();
		multiTraceback paths;
		trell.stochastic_traceback(paths, repetitions);
		print_output(&paths, seqs->getHeader());
    }
	else if (posterior){
		trell.posterior();
		multiTraceback paths;
		trell.traceback_stoch_posterior(paths, repetitions);
		print_output(&paths, seqs->getHeader());
	}
    else{
        std::cerr << usage << "\nNo Stochastic decoding option set\n";
        return;
    }
	
    return;
}


//Perform posterior decoding and print the output
void perform_posterior(model* hmm, sequences* seqs){
	trellis trell(hmm,seqs);
	
	//TODO: posterior should check model and choose the appropriate algorithm
	trell.posterior();
	
	//If we need a posterior traceback b/c path,label,or GFF is defined
	if (opt.isSet("-gff") || opt.isSet("-path") || opt.isSet("-label")){
		traceback_path path(hmm);
		trell.traceback_posterior(path);
		print_output(&path, seqs->getHeader());
	}
	else{
		print_posterior(trell);
	}
	
	return;
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


traceback_path* perform_traceback(trellis& trell){
	traceback_path* path = new(std::nothrow) traceback_path(trell.get_model());
	trell.traceback(*path);
    return path;
}

//Print the posterior probabilities for each state at each position
//Each state is in separate column
//Each row is on different row
void print_posterior(trellis& trell){
	model* hmm = trell.getModel();
	float_2D* table = trell.getPosteriorTable();
	
	std::ofstream file;
	std::string filename;
	std::streambuf* oldCoutStream(NULL);
	bool file_open(false);
	opt.getopt("-posterior",filename);
	
	if (!filename.empty()){
		oldCoutStream = std::cout.rdbuf();
		
		file.open(filename.c_str());
		std::cout << std::setprecision(3);
		std::cout << std::fixed;
		if (!file.is_open()){
			std::cerr << "Couldn't open file for posterior" << std::endl;
			exit(1);
		}
		file_open = true;
		std::cout.rdbuf(file.rdbuf());
	}
	
	
	std::cout <<"Posterior Probabilities Table\n";
	std::cout <<"Model:\t" << hmm->getName() << "\n";
	std::cout <<"Sequence:\t" << trell.getSeq()->getHeader() << "\n";
	std::cout <<"Probability of Sequence from Forward: Natural Log'd\t" << trell.getForwardProbability() << "\n";
	std::cout <<"Probability of Sequence from Backward:Natural Log'd\t" << trell.getBackwardProbability() << "\n";
	
	size_t state_size = hmm->state_size();
	
	//Print State names for header
	std::cout << "Position";
	for(size_t i=0; i < state_size; i++){
		std::cout << "\t";
		std::cout << hmm->getStateName(i);
	}
	std::cout << "\n";
	
	for(size_t position = 0; position < table->size(); position++){
		std::cout << position+1;
		for (size_t st = 0 ; st < state_size ; st++){
			std::cout << "\t";
			
			std::cout << exp((*table)[position][st]);
		}
		std::cout << "\n";
	}
	
	if (file_open){
		std::cout.rdbuf(oldCoutStream);
		file.close();
	}
	
	
	return;
	
}


