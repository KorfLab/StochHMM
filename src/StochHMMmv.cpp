#include <iostream>
#include <math.h>
#include <time.h>
#include "options.h"
#include "decoding.h"
#include "trellis.h"
#include "stochhmmUsage.h"

// TODO: Move HMMER header and structs into the HMMER_Interface headers.
//HMMER Headers
#include "config.h"		/* compile-time configuration constants */
#include "squidconf.h"
#include "squid.h"		/* general sequence analysis library    */
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
using namespace std;

struct plan7_s *Plan7HMM;
struct dpmatrix_s *Plan7Mat;
map<string,float> Scored;
bool HMMER_init=false;



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


//Need to think about how to deal with these multiple models and tracks
//Currently, it's kind of a hack job
HMM tempModel;
vector<HMM> models;  //Ultimately, I need to create a class to define some stats for each model
//Depending on selection.   //How do we define it 

vector<track> SEQtracks;
multiTrack SEQ;
bool multiModel;



int main (int argc, const char * argv[]) {
	srand(time(NULL));
    
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
        import_model(opt.sopt("-model"),opt.sopt("-seq"));
    }
    
    
    //Choose the model to use 
    //If model is defined use it, else choose what to do.
    HMM& final_model = (multiModel) ? models[0] : tempModel;
    
    
    final_model.print_tracks();
    final_model.print_model();
    
    if (HMMER_init){
        init_hmmsearch(opt.sopt("-hmmer"));
        HMMER_init=false;
    }
    
     
    
    ///Run Debug option
    if (opt.getopt("-debug","model")){
		final_model.print_model();
	}
	
	if (opt.getopt("-debug","seq")){
		final_model.print_tracks();
	}
    
    
    
    //Run non-stochastic options
    if (opt.optset("-viterbi")){
        
    }
    
    if (opt.optset("-nbest")){
        
    }
    
    
    
    //Run stochastic options
    if (opt.optset("-stoch")){
        trellis trell(final_model);
        
    }
    
	
	//for (int i=0;i<5;i++){
	//	cout << final_model.states[0].get_emission(i) <<endl;
	//}
	
	//clock_t start=clock();
	traceback_path x;

	//x=viterbi(final_model);
	//x.print_path();
	//cout << clock()-start << endl;
	
	//vector<traceback_path> list;
	//list=nth_best(final_model, 3);
	//for(int i=0;i<list.size();i++){
	//	cout << "PATH:\t" << i << endl;
	//	list[i].print_path();
	//	cout << endl;
	//}
	
	//start=clock();
	//traceback_path y;
	//y=viterbi_test_one(final_model);
	//y.print_path();
	//cout << clock()-start << endl;
	
	//start=clock();
	//traceback_path z;
	//z=viterbi_test_two(final_model);
	//z.print_path();
	//cout << clock()-start << endl;
	
	//start=clock();
	trellis trell(final_model);

	///
	cout << "begin test: " << endl;
	traceback_path path_test;
	long q = 100, w = 1; short r = 0;   //w=-1 means the init  //r=0 start with state (in position 0 of HMM vector<states>....) and start position(q) = 150
    cout << "START with 0:\t" << final_model.state_name[r] << endl;
	path_test=trell.traceback(q, r, w);
	path_test.print_path();
	//ofstream testout("arguments.txt", ios::out);
	//path_test.fprint_path(testout);
	
	//cout << "beginprinttracks" << " ";
	//HMM tracks_test;
	//tracks_test.print_tracks();
	//tracks_test._print_tracks();
	//cout << "end";
//	track track_test;
//	cout << track_test.alphabet.size();
//	cout << track_test.seq.size();
//	track_test.out_partial_track(testout);
	//testout.close();
	///
	
	//cout << "TRELLIS TIME:\t" << clock()-start << endl;

	x=trell.trace_viterbi();
	x.print_path();
	
	//trell.print(opt);
	

	
	return 0;
}
