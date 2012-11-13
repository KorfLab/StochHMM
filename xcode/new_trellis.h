//
//  new_trellis.h
//  StochHMM
//
//  Created by Paul Lott on 11/13/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef __StochHMM__new_trellis__
#define __StochHMM__new_trellis__

#include <iostream>
#include <vector>
#include <stdint.h>
#include <vector>
#include <string>
#include "stochTypes.h"


namespace StochHMM {

	struct tb_score{
		float score;
		signed int16_t tb_ptr;
	};
	
	
	typedef std::vector<std::vector<uint16_t> > int_2D;
	typedef std::vector<std::vector<float> > float_2D;
	typedef std::vector<std::vector<std::vector<uint16_t> > > int_3D;
	typedef std::vector<std::vector<std::vector<float> > > float_3D;

	class trellis{
	public:
		trellis();
		trellis(model* hmm , sequences* sq);
		~trellis();
		
		void viterbi();
		void forward();
		void backward();
		void posterior();
		void stochastic_viterbi();
		void stochastic_forward();
		void nth_viterbi();
		
		void viterbi(model* hmm, sequences *sq);
		void forward(model* hmm, sequences *sq);
		void backward(model* hmm, sequences *sq);
		void posterior(model* hmm, sequences *sq);
		void stochastic_viterbi(model* hmm, sequences *sq);
		void stochastic_forward(model* hmm, sequences *sq);
		void nth_viterbi(model* hmm, sequences *sq);
		
		void traceback(traceback_path&);
        void traceback(traceback_path&,size_t);
		void traceback_stoch_forward(multiTraceback&,size_t);
		void traceback_stoch_viterbi(multiTraceback&,size_t);
		void traceback_nth_viterbi(multiTraceback&);
		
		void baum_welch();
		
		inline bool store(){return store_values;}
		inline void store(bool val){store_values=val; return;}

		void print();
		std::string stringify();
		void export_trellis(ifstream&);
		void export_trellis(std::string& file);
		
	private:
		void init_table();
		double getEndingTransition(size_t);
        double getTransition();
        int traceback_length();
        double exFuncTraceback(transitionFuncParam*);
		
		
		model* hmm;		//HMM model
        sequences* seqs; //Digitized Sequences
		size_t state_size;	//Number of States
		size_t seq_size;	//Length of Sequence
		
		trellisType type;
		
		bool store_values;
		bool exDef_defined;
		
		state* initial;
		state* ending;
		
		
		//Traceback Tables
		int_2D*		traceback;          //Simple traceback table
		float_3D*	stoch_traceback;    //Stochastic traceback table
		int_3D*		nth_traceback;      //Nth-Viterbi traceback table
		
		//Score Tables
		float_2D*	viterbi_score;      //Storing viterbi scores
		float_2D*	forward_score;      //Storing Forward scores
		float_2D*	backward_score;     //Storing Backward scores
		float_2D*	posterior;			//Store posterior scores
		
		//Ending Cells
		double   ending_viterbi_score;
		uint16_t ending_viterbi_tb;
		double   ending_forward_score;

		std::vector<tb_score>*    ending_stoch_tb;
		std::vector<tb_score>*    ending_nth_viterbi;
		
		std::vector<double>* current_calc_cells;
		std::vector<double>* previous_calc_cells;
		std::vector<double>* swap_cells;
	};
	
}

#endif /* defined(__StochHMM__new_trellis__) */
