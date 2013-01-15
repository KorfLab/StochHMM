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
#include <fstream>
#include <bitset>
#include "stochTypes.h"
#include "sequences.h"
#include "hmm.h"
#include "traceback_path.h"
#include "stochMath.h"
#include <stdint.h>


namespace StochHMM {

//	struct tb_score{
//		float score;
//		int16_t tb_ptr;
//	};
	
	
	
	class stochTable{
	public:
		
		struct stoch_value{
			stoch_value(uint16_t id, uint16_t prev, float p): state_id(id), state_prev(prev), prob(p){}
			uint16_t state_id;
			uint16_t state_prev;
			uint16_t prev_cell;
			float prob;
		};
		
		
		stochTable(size_t);
		~stochTable();
		void push(size_t pos, size_t st, size_t st_to, float val);
		void print();
		void finalize();
		uint16_t get_state_position(size_t pos,uint16_t);
		
		void traceback(traceback_path& path);
		
	private:
		size_t last_position;
		std::vector<stoch_value>* state_val;
		std::vector<size_t>* position;
	};
	
	
	typedef std::vector<std::vector<uint16_t> > int_2D;
	typedef std::vector<std::vector<float> > float_2D;
	typedef std::vector<std::vector<std::vector<uint16_t> > > int_3D;
	typedef std::vector<std::vector<std::vector<float> > > float_3D;

	class trellis{
	public:
		trellis();
		trellis(model* h , sequences* sqs);
		~trellis();
		void reset();
		
		
		/*-----------   Decoding Algorithms ------------*/
		
		//TODO: Fix these functions so that they evaluate the model and choose a
		// a default algorithm based on the whether the model defines explicit duration
		// states
		
		void viterbi();
		void viterbi(model* h, sequences* sqs);
		
		void forward();
		void forward(model* h, sequences* sqs);
		
		void forward_viterbi();
		void forward_viterbi(model* h, sequences* sqs);
		
		void backward();
		void backward(model* h, sequences* sqs);
		
		void posterior();
		void posterior(model* h, sequences* sqs);
		
		void stochastic_viterbi();
		void stochastic_viterbi(model* h, sequences* sqs);
		
		void stochastic_forward();
		void stochastic_forward(model* h, sequences* sqs);
		
		void nth_viterbi();
		void nth_viterbi(model* h, sequences *sqs);
		
		void baum_welch();
		
		
		/*-----------   Simple Model Decoding Algorithms ------------*/
		/* These algorithms are for use with models that do not define explicit
			duration states or external functions.   Therefore, probabilities can
			be easily computed or referenced.
		 
		 */
		
		void simple_viterbi();
		void simple_viterbi(model* h, sequences* sqs);
				
		void simple_forward_viterbi();
		void simple_forward_viterbi(model* h, sequences* sqs);
		
		void simple_forward();
		void simple_forward(model* h, sequences* sqs);
		
		void simple_backward();
		void simple_backward(model* h, sequences* sqs);
		
		void simple_posterior();
		void simple_posterior(model* h, sequences* sqs);
		
		void simple_stochastic_viterbi();
		void simple_stochastic_viterbi(model* h, sequences* sqs);
		
		void simple_stochastic_forward();
		void simple_stochastic_forward(model* h, sequences* sqs);
		
		void simple_baum_welch();
		void simple_baum_welch(model* h, sequences* sqs);

		
		/*-----------   Fast Complex Model Decoding Algorithms  ----------*/
		/*	These algorithms are for use with models that define external functions
			or explicit duration states.
		 
			These algorithms store the associated emissions and transitions
			probabilities that were calculated during viterbi algorith.  For later
			use in the forward/backward algorithms.
		 
			The values are stored in dense arrays, so the memory is accessed faster
			at the cost of amount of memory required.
		 */
		
		void fast_complex_viterbi();
		void fast_complex_viterbi(model* h, sequences* sqs);
		
		void fast_complex_forward_viterbi();
		void fast_complex_forward_viterbi(model* h, sequences* sqs);
		
		void fast_complex_backward();
		void fast_complex_backward(model* h, sequences* sqs);
		
		void fast_complex_stochastic_viterbi();
		void fast_complex_stochastic_viterbi(model* h, sequences* sqs);
		
		void fast_complex_stochastic_forward();
		void fast_complex_stochastic_forward(model* h, sequences* sqs);
		
		void fast_complex_baum_welch();
		void fast_complex_baum_welch(model* h, sequences* sqs);
		
		
		/*-----------   Sparse Complex Model Decoding Algorithms  --------*/
		/*	These algorithms store the associated emissions and transitions
			probabilities that were calculated during viterbi algorith.  For later
			use in the forward/backward algorithms.
		
			The values are stored in sparse tables (hash map), so access to values 
			is slower but with lower amount of memory required.
		 */
		
		void sparse_complex_viterbi();
		void sparse_complex_viterbi(model* h, sequences* sqs);
		
		void sparse_complex_foward_viterbi();
		void sparse_complex_foward_viterbi(model* h, sequences* sqs);
		
		void sparse_complex_backward();
		void sparse_complex_backward(model* h, sequences* sqs);
		
		void sparse_complex_stochastic_forward();
		void sparse_complex_stochastic_forward(model* h, sequences* sqs);
		
		void sparse_complex_baum_welch();
		void sparse_complex_baum_welch(model* h, sequences* sqs);
		
		
		/*---------   Standard Traceback Algorithms --------*/
				
		void traceback(traceback_path& path);
        void traceback(traceback_path&,size_t);
		void traceback_stoch_forward(multiTraceback&,size_t);
		void traceback_stoch_viterbi(multiTraceback&,size_t);
		void traceback_nth_viterbi(multiTraceback&);
		
		void stochastic_traceback(traceback_path& path);
		
		
		inline bool store(){return store_values;}
		inline void store(bool val){store_values=val; return;}

		void print();
		std::string stringify();
		void export_trellis(std::ifstream&);
		void export_trellis(std::string& file);
		inline model* get_model(){return hmm;}
		
	private:
		double getEndingTransition(size_t);
        double getTransition(state* st, size_t trans_to_state, size_t sequencePosition);
        size_t get_explicit_duration_length(transition* trans, size_t sequencePosition,size_t state_iter, size_t to_state);
        double exFuncTraceback(transitionFuncParam*);
		
		
		model* hmm;		//HMM model
        sequences* seqs; //Digitized Sequences
		
		size_t state_size;	//Number of States
		size_t seq_size;	//Length of Sequence
		
		trellisType type;
		
		bool store_values;
		bool exDef_defined;
		
		//Traceback Tables
		int_2D*		traceback_table;          //Simple traceback table
		int_3D*		nth_traceback;      //Nth-Viterbi traceback table
		stochTable* stochastic_table;
		
		//Score Tables
		float_2D*	viterbi_score;      //Storing viterbi scores
		float_2D*	forward_score;      //Storing Forward scores
		float_2D*	backward_score;     //Storing Backward scores
		float_2D*	posterior_score;			//Store posterior scores
		
		//Ending Cells
		double   ending_viterbi_score;
		uint16_t ending_viterbi_tb;
		
		double   ending_posterior;

//		std::vector<tb_score>*    ending_stoch_tb;
//		std::vector<tb_score>*    ending_nth_viterbi;
		
		//Cells used for calculating the Viterbi
		std::vector<double>* viterbi_current;
		std::vector<double>* viterbi_previous;
		
		//Array used for calculating the Backward Scores
		std::vector<double>* backward_current;
		std::vector<double>* backward_previous;
		
		std::vector<size_t>* explicit_duration_current;
		std::vector<size_t>* explicit_duration_previous;
		std::vector<size_t>* swap_ptr_duration;
		
		std::vector<double>* swap_ptr;
	};
	
}

#endif /* defined(__StochHMM__new_trellis__) */
