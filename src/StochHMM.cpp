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
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <time.h>
#include <fstream>
#include "StochHMMlib.h"

#include "StochHMM_usage.h"
using namespace StochHMM;

#define STATE_MAX 1024


void import_model(model&);
void import_sequence(model&);

void perform_viterbi_decoding(model* hmm, sequences* seqs);
void perform_nbest_decoding(model* hmm, sequences* seqs);
void perform_posterior(model* hmm, sequences* seqs);
void perform_stochastic_decoding(model* hmm, sequences* seqs);

void print_output(multiTraceback*, std::string&);
void print_output(std::vector<traceback_path>&, std::string&);
void print_output(traceback_path*, std::string&);
void print_posterior(trellis&);


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

//Stores the number of options in opt
int opt_size=sizeof(commandline)/sizeof(commandline[0]);

//Global options for parsed command-line options
options opt;

//seqTracks stores multiple jobs(model and multiple sequences)
seqTracks jobs;

//Create and initialize StateFuncs
//This will automatically initialize all the Univariate and Multivariate
//PDFs
StateFuncs default_functions;


int main(int argc, const char * argv[])
{
    srand(std::time(NULL));
    //Parse commandline arguments defined
    opt.set_parameters(commandline,opt_size,usage);
	
	//Get and parse command line options supplied by user
    opt.parse_commandline(argc,argv);
    
    
	//Create a model
    model hmm;
	
    //import the model
    import_model(hmm);
    
	
//	if (opt.isFlagSet("-debug", "paths")){
//		hmm.writeGraphViz("Model_Path-debug.viz", true);
//	}
	
    
    //Check and import sequence(s)
	//These will be imported into seqTracks jobs
    import_sequence(hmm);
    
	//Get the job (model and associated sequences)
    seqJob *job=jobs.getJob();
	
	
	//If filename is set for any of the following
	//options we need to re-direct the stdout to the file
	std::string filename;
	if (opt.isSet("-posterior")){
		opt.getopt("-posterior",filename);
	}
	else if (opt.isSet("-gff")){
		opt.getopt("-gff", filename);
	}
	else if (opt.isSet("-path")){
		opt.getopt("-path", filename);
	}
	else if (opt.isSet("-label")){
		opt.getopt("-label", filename);
	}
	
	
	//Create file handle
	std::ofstream file;
	
	//Create ptr to store STDOUT if we redirect 
	std::streambuf* oldCoutStream(NULL);
	bool file_open(false);

	
	//If we defined a filename then redirect stdout
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
	
	
	// Fore each job(sequence) perform the analysis
	while (job != NULL){
		
		//Print sequences if -debug seq option defined
		if (opt.isFlagSet("-debug","seq")){
			job->getSeqs()->print();
		}
		
		//Perform posterior analysis
		if (opt.isSet("-posterior")){
			perform_posterior(&hmm, job->getSeqs());
		}
		
		//Perform viterbi analysis
		else if(opt.isSet("-viterbi")){
			perform_viterbi_decoding(job->getModel(), job->getSeqs());
		}
		
		//Perform nbest Viterbi decoding
		else if (opt.isSet("-nbest")){
			perform_nbest_decoding(job->getModel(), job->getSeqs());
		}
		
		//Perform stochastic decoding
		else if (opt.isSet("-stochastic")){
			perform_stochastic_decoding(job->getModel(), job->getSeqs());
		}
		
		//Get next job
		job = jobs.getJob();
	}
	
	
	//Close file and redirect stdout to original place
	if (file_open){
		std::cout.rdbuf(oldCoutStream);
		file.close();
	}
    
    return 0;
    
}

//Import the model from file
void import_model(model& hmm){
    if (!opt.isSet("-model")){
        std::cerr <<"No model file provided.\n" << usage << std::endl;
    }
    else{

        
		//Import the model using string supplied from commandline
		//Pass StateFuncs.   This will allow the user to define
		//the functions within the model.
		hmm.import(opt.sopt("-model"),&default_functions);
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
		//Fastq import is still preliminary, it will only
		//import the sequence, not the qualities
		jobs.loadSeqs(hmm,opt.sopt("-seq"),FASTQ);
	}
    else{
		//Import the sequence form fasta file
        jobs.loadSeqs(hmm, opt.sopt("-seq"), FASTA);
    }
}

//Perform Viterbi decoding and print the output
void perform_viterbi_decoding(model* hmm, sequences* seqs){
	//Setup the trellis with the model and sequence
    trellis trell(hmm,seqs);
	
	//Perform viterbi decoding
	trell.viterbi();
	
	//Create a traceback path ptr to store traceback from perform_traceback
	//function
	traceback_path path(hmm);
	trell.traceback(path);
		
	//Call print_output (below) to print the traceback in the required format
	print_output(&path, seqs->getHeader());
	
    return;
}

//Perform nth-best decoding and print the output
void perform_nbest_decoding(model* hmm, sequences* seqs){
	//Setup the trellis with the model and sequence
	trellis trell(hmm,seqs);
	
	//Get the number of paths to get
	size_t nth = opt.iopt("-nbest");
	
	//Perform nth viterbi
	trell.naive_nth_viterbi(nth);
	
	//Get the N tracebacks and output them
	for(size_t i=0;i<nth;i++){
		traceback_path path(hmm);
		trell.traceback_nth(path, i); //ith path
		print_output(&path, seqs->getHeader());
	}
}


//Perform stochastic decoding
void perform_stochastic_decoding(model* hmm, sequences* seqs){
    
    //Determine which type of stochastic algorithm to perform
    bool viterbi	= (opt.isFlagSet("-stochastic", "viterbi") || opt.isSet("-viterbi")) ? true : false;
    bool forward	= (opt.isFlagSet("-stochastic", "forward")) ? true : false;
    bool posterior	= (opt.isFlagSet("-stochastic", "posterior")) ? true : false;
	
	//Setup the trellis with the model and sequence
    trellis trell(hmm,seqs);
	
	//Number of times to traceback over path
	int repetitions = opt.iopt("-rep");
    
    if (viterbi){
		clock_t start = clock();
		trell.stochastic_viterbi();
		//create multiple paths object to stor
		multiTraceback paths;
		trell.stochastic_traceback(paths, repetitions);
//		clock_t stop = clock();
//		
//		std::cout << (double) stop-start/ (double) CLOCKS_PER_SEC << std::endl;
		
		//print_output(&paths, seqs->getHeader());
		
		
//		start = clock();
//		trellis alt_trell(hmm,seqs);
//		alt_trell.simple_alt_stochastic_viterbi(hmm,seqs);
//		multiTraceback alt_paths;
//		trell.stochastic_traceback(alt_paths, repetitions);
//		stop = clock();
//		
//		std::cout << (double) stop-start/ (double) CLOCKS_PER_SEC << std::endl;
		
		
//		start = clock();
//		trellis simple_trellis(hmm,seqs);
//		simple_trellis.simple_simple_stochastic_viterbi(hmm,seqs);
//		multiTraceback simple_paths;
//		trell.stochastic_traceback(simple_paths, repetitions);
//		stop = clock();
//		
//		std::cout << (double) stop-start/ (double) CLOCKS_PER_SEC << std::endl;
//		//print_output(&simple_paths, seqs->getHeader());
//		//create multiple paths object to stor
//		multiTraceback paths;
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
        std::cout << "#Score: " << tb->getScore() << std::endl;
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


//Print the posterior probabilities for each state at each position
//Each state is in separate column
//Each row is on different row
void print_posterior(trellis& trell){
	model* hmm = trell.getModel();
	double_2D* table = trell.getPosteriorTable();
	size_t state_size = hmm->state_size();
	char cstr[200];
	
	std::string output;
	output+="Posterior Probabilities Table\n";
	output+="Model:\t" + hmm->getName() + "\n";
	output+="Sequence:\t" + trell.getSeq()->getHeader() + "\n";
	sprintf(cstr, "Probability of Sequence from Forward: Natural Log'd\t%f\n",trell.getForwardProbability());
	output+= cstr;
	sprintf(cstr, "Probability of Sequence from Backward:Natural Log'd\t%f\n",trell.getBackwardProbability());
	output+= cstr;
	output+= "Position";
	for(size_t i=0;i< state_size; ++i){
		output+= "\t" + hmm->getStateName(i);
	}
	output+="\n";
	
	std::cout <<  output;
	

	for(size_t position = 0; position < table->size(); ++position){
		sprintf(cstr, "%ld", position+1);
		output= cstr;
		for (size_t st = 0 ; st < state_size ; st++){
			float val  = exp((*table)[position][st]);
			if (val<= 0.001){
				output+="\t0";
			}
			else if (val == 1.0){
				output+="\t1";
			}
			else{
				sprintf(cstr,"\t%.3f", exp((*table)[position][st]));
				output+= cstr;
			}

		}
		output+="\n";
		std::cout << output;
	}

	std::cout << std::endl;
	
	return;
	
}


