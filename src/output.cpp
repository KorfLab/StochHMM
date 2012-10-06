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

///*
// *  output.cpp
// *  StochHMM3
// *
// *  Created by Quidage on 8/27/08.
// *  Copyright 2008 University of California, Davis. All rights reserved.
// *
// */
//
//#include "output.h"
//namespace StochHMM{
//
////Takes trellis once all filled and outputs the tracebacks as label
//
//
//void map_label (trellis &scores, const HMM &model, options &opt){
//	int seq_size=scores.trell[0].size();
//	string modl=opt.sopt("-model");
//    string seq=opt.sopt("-seq");
//    string lbl=opt.sopt("-label");
//
//    //Reassign stdout to file if output filename is defined
//	streambuf* sbuf=cout.rdbuf();
//	ofstream file;
//    
//	if (lbl.compare("")!=0){
//		file.open(lbl.c_str());
//		cout.rdbuf(file.rdbuf());		
//	}
//    
//    
//	string line (50,'#');
//	
//	map<string,string> label;
//	for(int state=0;state<model.states.size();state++){
//		//cout << model.state[state].label <<endl;
//		//cout << model.state[state].name <<endl;
//		label[model.states[state].label]+=model.states[state].name;
//	}
//
//	cout << "#Model: " << modl <<endl;
//	cout << "#Sequence: " << seq << endl;
//	cout << "#Algorithm: Viterbi" <<endl;
//	cout << "#Output: Label" << endl;
//	cout << "#Label Key: ";
//	
//	map<string,string>::iterator itera;	
//	for (itera=label.begin(); itera!=label.end(); itera++){
//		cout << "(" << (*itera).first << "=" << (*itera).second << "),";
//	}
//	
//	cout <<endl;
//	
//	traceback_path (trellis::*funct)();
//	
//	if (opt.getopt("-stoch","viterbi")){
//		funct=&trellis::trace_stoch_viterbi;
//	}
//	else if (opt.getopt("-stoch","forward")){
//		funct=&trellis::trace_stoch_forward;
//	}
//    
//	
//	map<traceback_path,int> paths;
//	int max=0;
//	
//	double posterior=scores.trell[0][seq_size-1].forw;
//	for(unsigned int r=1;r<model.states.size(); r++){
//		posterior+=log(1+exp(scores.trell[r][seq_size-1].forw-posterior));		
//	}
//	
//	for(int i=0;i<opt.iopt("-repetitions");i++){
//		traceback_path theta=(scores.*funct)();
//		paths[theta]++;
//		max = (paths[theta]>max) ? paths[theta] : max;
//	}
//	
//	
//	//print top n paths type vector<int>
//	int sum=0;
//	int i=0;
//	
//	map<traceback_path, int>::iterator it;
//	//int ranking=1;
//	
//	for (int k=max;k>0;k--){
//		for( it=paths.begin(); it!=paths.end();it++){
//			if ((*it).second==k){
//				cout << line << endl;
//
//				cout << "Path occured " << k << " times:" <<endl;
//				traceback_path final=(*it).first;
//				
//				//Determine Path probability
//				//double path_probi=final.path_prob(model,seq);
//				
//				//if (path_probi==-INFINITY){
//				//	continue;
//				//}
//				//cout << "Path probability:\texp(" << path_probi<<")"<<endl;
//				//double post=100*exp(path_probi-posterior);
//				
//				sum+=k;
//		
//				final.print_label(model);
//				
//								
//				i++;
//				if (i>=opt.iopt("-report") || sum==opt.iopt("-repetitions")){
//					goto LAST_PRINT;
//				}
//			}
//		}
//	};
//LAST_PRINT:
//	cout <<endl;
//	
//    //Re-assign the cout to stdout
//	if (lbl.compare("")!=0){
//		cout.rdbuf(sbuf);
//		file.close();
//	}
//}
//
//
////Mapping stochastic traceback
//void map_path (trellis &scores, const HMM &model, options &opt){
//	int seq_size=scores.trell[0].size();
//	string modl=opt.sopt("-model");
//    string seq=opt.sopt("-seq");
//    string path=opt.sopt("-path");
//    
//    
//	streambuf* sbuf=cout.rdbuf();
//	ofstream file;
//	
//	
//	if (path.compare("")!=0){
//		file.open(path.c_str());
//		cout.rdbuf(file.rdbuf());
//	}
//	
//	string line (50,'#');
//	
//	cout << "#Model: " << modl <<endl;
//	cout << "#Sequence: " << seq<< endl;
//	cout << "#Algorithm: Viterbi" <<endl;
//	cout << "#Output: Path" << endl;
//	cout << "#Path Key: ";
//	
//	for (int i=0; i<model.state_name.size()-1; i++){
//		cout << "(" << i << "=" << model.state_name[i] << "),";
//	}
//	
//	cout <<endl;
//	
//	
//	traceback_path (trellis::*funct)();
//	
//	if (opt.getopt("-stoch","viterbi")){
//		funct=&trellis::trace_stoch_viterbi;
//	}
//	else {
//		funct=&trellis::trace_stoch_forward;
//	}	
//	
//	map<traceback_path,int> paths;
//	int max=0;
//	
//	
//
//	double posterior=scores.trell[0][seq_size-1].forw;
//	for(unsigned int r=1;r<model.states.size(); r++){
//		posterior+=log(1+exp(scores.trell[r][seq_size-1].forw-posterior));		
//	}
//	//cout << "Posterior:\texp(" << posterior <<")"<<endl;
//
//	//cout << "Viterbi:\texp("<< scores.trell[scores.vit_end][sequence.size()-1].viti <<")"<< endl<<endl;
//	
//	for(int i=0;i<opt.iopt("-repetitions");i++){
//		traceback_path theta=(scores.*funct)();
//		paths[theta]++;
//		max = (paths[theta]>max) ? paths[theta] : max;
//	}
//	
//	
//	//print top n paths type vector<int>
//	int sum=0;
//	int i=0;
//	
//	map<traceback_path, int>::iterator it;
//	//int ranking=1;
//
//	for (int k=max;k>0;k--){
//		for( it=paths.begin(); it!=paths.end();it++){
//			if ((*it).second==k){
//				cout << line << endl;
//
//				cout << "Path occured " << k << " times:" <<endl;
//				traceback_path final=(*it).first;
//
//				//Determine Path probability
//				//double path_probi=final.path_prob(model,seq);
//				
//				//if (path_probi==-INFINITY){
//				//	continue;
//				//}
//				//cout << "Path probability:\texp(" << path_probi<<")"<<endl;   //Need to rewrite function for path probability
//				//double post=100*exp(path_probi-posterior);
//				
//				sum+=k;
//
//				final.print_path();
//				
//				/*
//				 if (opt.gff_flag==1){
//					final.print_gff(model,opt.sequence,path_probi,ranking, k,posterior);
//					ranking++;
//				}
//				*/
//				
//				i++;
//				if (i>=opt.iopt("-report") || sum==opt.iopt("-repetitions")){
//					goto LAST_PRINT;
//				}
//			}
//		}
//	};
//LAST_PRINT:
//	cout <<endl;
//	if (path.compare("")!=0){
//		cout.rdbuf(sbuf);
//		file.close();
//	}
//}
//
//
////NEED to account for posterior based traceback in map_GFF
////Mapping stochastic viterbi-based traceback
//void map_gff (trellis &scores, const HMM &model, options& opt){
//	int seq_size=scores.trell[0].size();
//	string modl=opt.sopt("-model");
//    string seq=opt.sopt("-seq");
//    string gff=opt.sopt("-gff");
//
//	streambuf* sbuf=cout.rdbuf();
//	ofstream file;
//	
//	if (gff.compare("")!=0){
//		file.open(gff.c_str());
//		cout.rdbuf(file.rdbuf());
//	}
//	
//	string line (50,'#');
//	int repetitions=opt.iopt("-repetitions");
//	
//	map<vector<char>,int> paths;
//	map<char,string> label_key;
//	
//	for(int i=0;i<model.states.size();i++){
//		if (!model.states[i].desc.empty()){
//			label_key[model.states[i].label[0]]=model.states[i].desc;
//			//cout << i << "\t" << model.state[i].label[0] << "\t" << model.state[i].desc <<endl;
//		}
//	}
//	
//	int max=0;
//	
//	//Calculate total posterior probability
//	
//	double posterior=scores.trell[0][seq_size-1].forw;
//	for(unsigned int r=1;r<model.states.size(); r++){
//		posterior+=log(1+exp(scores.trell[r][seq_size-1].forw-posterior));		
//	}
//	
//	//Perform repetitions by tracingback trellis
//	for(int i=0;i<repetitions;i++){
//		vector<char> theta=scores.trace_stoch_viterbi_char();
//		paths[theta]++;
//		max = (paths[theta]>max) ? paths[theta] : max;
//	}
//	
//	//print top n paths type: GFF
//	int sum=0;
//	int i=0;
//	
//	map<vector<char>, int>::iterator it;
//	int ranking=1;
//	
//	for (int k=max;k>0;k--){
//		for( it=paths.begin(); it!=paths.end();it++){
//			if ((*it).second==k){
//
//				//print_gff(model,(*it).first,opt.sequence,strand,path_probi,ranking, k,posterior);
//				convert_label(label_key,(*it).first,opt,ranking,k,posterior);
//				ranking++;
//				i++;
//				if (sum==repetitions){
//					goto LAST_PRINT;
//				}
//			}
//		}
//	};
//LAST_PRINT:
//	cout <<endl;
//	
//	if (gff.compare("")!=0){
//		cout.rdbuf(sbuf);
//		file.close();
//	}
//}
//
//
//
//void convert_label( map<char,string> label_key, vector<char> path, options opt, int ranking,int times,double posterior){
//	string sequence_name=opt.sopt("-seq");
//	char current_label;
//	double start=0;
//	int path_size=path.size();
//	
//	for(int k=path_size-1;k>=0;k--){
//		
//		char new_label=path[k];
//		
//		if (label_key.count(new_label)== 0){
//			if (start>0){
//				//cout << current_label << endl;
//				cout << sequence_name << "\tStochHMM\t" << label_key[current_label] <<"\t"<< start << "\t" << path_size-(k+1) << "\t.\t+\t.\tRank:"<<ranking<<",Counts:" << times<<",Posterior:"<<posterior<<endl;
//				start=0;
//				current_label=new_label;
//			}
//			else {
//				continue;
//			}
//		}
//		else {
//			if(k==0){
//				//cout << current_label << endl;
//				cout << sequence_name << "\tStochHMM\t" << label_key[current_label] <<"\t"<< start << "\t" << path_size << "\t.\t+\t.\tRank:"<<ranking<<",Counts:" << times<<",Posterior:"<<posterior<<endl;
//			}
//			else if (start==0){
//				start=path_size-k;
//				current_label=new_label;
//			}
//			else if (new_label==current_label){
//				continue;
//			}
//			else {
//				//cout << current_label << endl;
//				cout << sequence_name << "\tStochHMM\t" << label_key[current_label] <<"\t"<< start << "\t" << path_size-(k+1) << "\t.\t+\t.\tRank:"<<ranking<<",Counts:" << times<<",Posterior:"<<posterior<<endl;
//				start=path_size-k;
//				current_label=new_label;
//			}
//		}	
//	}
//}
//
//
////Outputs the traceback state hit file for heat map.
////Need to integrate with previous map functions
////void map_heat (trellis &scores, const HMM &model, sequence &seq, int repetitions, int top_N, int flag, options &opt, string strand){
//void map_heat (trellis &scores, const HMM &model, options &opt){
//	
//
//	int states=scores.trell.size();
//	//int seq_length=model.tracks.size();
//	int seq_length=scores.trell[0].size();
//	
//	map<int,int> st;
//	vector<map<int,int> > cell (seq_length,st);
//	vector<vector< map<int,int> > > output (model.states.size(),cell);
//	
//	traceback_path (trellis::*funct)();
//	
//	if (opt.getopt("-stoch","viterbi")){
//		funct=&trellis::trace_stoch_viterbi;
//	}
//	else {
//		funct=&trellis::trace_stoch_forward;
//	}	
//	
//		
//	//cout << "Heatmap" << opt.heatmap <<endl;
//	
//	streambuf* sbuf=cout.rdbuf();
//	ofstream file;
//	file.open(opt.sopt("-heat").c_str());
//	cout.rdbuf(file.rdbuf());
//	
//	for(int i=0;i<opt.iopt("-repetitions");i++){
//		traceback_path theta= (scores.*funct)();
//		for(int j=theta.size()-1;j>=0;j--){
//			if (j==0){
//				output[theta.val(j)][j][-1]++;
//			}
//			else{
//				output[theta.val(j)][j][theta.val(j-1)]++;
//			}
//			//cout << theta[j] << endl;
//		}
//	}
//	
//	cout << "STATE KEY:\t(-1,END),";
//	for(int k=0;k<model.states.size();k++){
//		cout << "(" << k << "," << model.states[k].name << ")";
//		if (k+1<states){
//			cout << ",";
//		}
//	}
//	cout <<endl;
//	
//	/*
//	cout << "SEQUENCE:\t";
//	
//	for(int j=0;j<seq.size();j++){
//		cout << model.alpha[seq.val(j)];
//	}
//	cout <<endl;
//	*/
//	
//	
//	cout << "HEADER KEY:\tSTATE,SEQUENCE,(PREV_STATE,TRACEBACK_COUNTS)" << endl;
//	
//	/*  //Sequence first
//	for (int j=seq_length-1;j>=0;j--){
//		for (int i=0;i<output.size();i++){
//			map<int,int>::iterator it;
//			for(it=output[i][j].begin();it!=output[i][j].end();it++){
//				int position=seq_length-1-j;
//				cout <<position<< "," << i <<","<< it->first <<"," << it->second << endl;
//			}
//		}
//	}
//	*/
//	
//	//State first
//	for (int i=0;i<output.size();i++){
//		for (int j=output[i].size()-1;j>=0;j--){
//			map<int,int>::iterator it;
//			if (output[i][j].size()>0){
//				int position=seq_length-1-j;
//				cout << i << "," << position;
//			}
//			for(it=output[i][j].begin();it!=output[i][j].end();it++){
//				cout  <<",("<< it->first <<"," << it->second << ")";
//			}
//			if (output[i][j].size()>0){
//				cout << endl;
//			}
//			
//		}
//	}
//	cout.rdbuf(sbuf);
//	file.close();
//}
//}
